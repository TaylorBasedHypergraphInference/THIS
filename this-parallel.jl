using Distributed

@everywhere include("this.jl")

# ================================================================================
"""
	this_par(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)

Parallel version of THIS. Infers the hypergraph underlying the dynamics of its vertices with knowledge of the states 'X' and of the derivatives 'Y' at each vertex.

Requires to store some data in a directory named `data`.

_INPUT_:\\
`X`: Time series of the system's state. Each row is the time series of the state of one agent.\\
`Y`: Time series of the system's velocity. Each row is the time series of the velocity of on agent.\\
`ooi`: Orders of interest. Vector of integers listing the orders of interactions that we analyze.\\
`dmax`: Maximal degree to be considered in the Taylor expansion. Typically, dmax=maximum(ooi).\\
`λ`: SINDy's threshold deciding whether an hyperedge exists or not.\\
`ρ`: Regularization parameters promoting sparsity.\\
`niter`: Maximal number of iterations for THIS algorithm.

_OUTPUT_:\\
`Ainf`: Dictionary associating the inferred coefficient of the Taylor expansion to a pair (node,hyperedge). Namely, 'Ainf[(i,h)]' is the coefficient corresponding to the hyperedge 'h' in the Taylor expansion of the dynamics of node 'i'. If agents have internal dimensions larger than 1, `Ainf` is just a boolean.\\
`coeff`: Matrix of coefficents obtained by SINDy. The row indices are the agents' indices and the columns indices are the indices of the monomials in the Taylor series.\\
`relerr`: Relative error, i.e., `err` normalized by the magnitude of `Y`.
"""
function this_par(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	# Size of the data matrix. 
	n,T = size(X)

	if size(X) != size(Y)
		@info "Dimensions of states and derivatives do not match."
		return nothing
	end

	# Running THIS ##############################################################
	coeff = Dict{Int64,Vector{Float64}}()
	err = 0.

	args = [(i,X,Y[[i,],:],dmax,λ,ρ,niter) for i in 1:n]
	pmap(mySINDy_par,args)

	for i in 1:n
		coeff[i] = readdlm("data/temp-sindy-Xi-$i.csv",',')
		rm("data/temp-sindy-Xi-$i.csv")
		err += readdlm("data/temp-sindy-err-$i.csv",',')[1]
		rm("data/temp-sindy-err-$i.csv")
	end

	relerr = err/norm(Y,1)

	@info "THIS completed."
	#############################################################################
	
	# Reconstructing the adjacency tensors ######################################
	d = get_d(n,dmax)
	idx_mon = Dict{Int64,Vector{Int64}}()
	for i in 1:size(d)[1]
		mon = d[i,:][d[i,:] .!= 0]
		if length(mon) == length(union(mon))
			idx_mon[i] = sort(mon)
		end
	end

	Ainf = Dict{Int64,Matrix{Float64}}(o => zeros(0,o+1) for o in 1:dmax+1) # For each order o, associates a dictionary associating the pair (agents, hyperedge) to the inferred weight, as seen from the agent.
	for i in 1:n
		c = coeff[i]
		for j in 1:length(c)
			mon = idx_mon[j]
			if !(i in mon)
				Ainf[length(mon)+1] = vcat(Ainf[length(mon)+1],[i mon' c[j]])
			end
		end
	end

	@info "Dictionary of inferred adjacency tensors built."
	#############################################################################

	return Ainf, coeff, relerr
end

# ================================================================================
"""
	this_par_filter(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, par_keep::Union{Int64,Float64}=.9, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)

Parallel verion of `this_filter`. Infers the hypergraph underlying the dynamics of its vertices with knowledge of the states 'X' and of the derivatives 'Y' at each vertex.

Requires to store some data in a directory named `data`.

_INPUT_:\\
`X`: Time series of the system's state. Each row is the time series of the state of one agent.\\
`Y`: Time series of the system's velocity. Each row is the time series of the velocity of on agent.\\
`ooi`: Orders of interest. Vector of integers listing the orders of interactions that we analyze.\\
`dmax`: Maximal degree to be considered in the Taylor expansion. Typically, dmax=maximum(ooi).\\
`par_keep`: Parameter for the pairs to keep. If `Float64`, correlation threshold above which the pairs of time series should be kept. If `Int64`, number of pairs of time series to be kept by filtering (based on correlation).\\
`λ`: SINDy's threshold deciding whether an hyperedge exists or not.\\
`ρ`: Regularization parameters promoting sparsity.\\
`niter`: Maximal number of iterations for THIS algorithm.

_OUTPUT_:\\
`Ainf`: Dictionary associating the inferred coefficient of the Taylor expansion to a pair (node,hyperedge). Namely, 'Ainf[(i,h)]' is the coefficient corresponding to the hyperedge 'h' in the Taylor expansion of the dynamics of node 'i'. If agents have internal dimensions larger than 1, `Ainf` is just a boolean.\\
`coeff`: Matrix of coefficents obtained by SINDy. The row indices are the agents' indices and the columns indices are the indices of the monomials in the Taylor series.\\
`relerr`: Relative error, i.e., `err` normalized by the squared magnitude of `Y`.
"""
function this_par_filter(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, par_keep::Union{Int64,Float64}=.9, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	if size(X) != size(Y)
		@info "Dimensions of states and derivatives do not match."
		return nothing
	end

	# Listing pairs of agents whose correlation coefficient is larger than 'α'.
	keep = keep_correlated(X,par_keep)

	# Running THIS ##############################################################
	n,T = size(X)

	d = get_d(n,dmax) # Lists the indices of the agents involved in each monomial.
	idx_mon = Dict{Int64,Vector{Int64}}()
	for i in 1:size(d)[1]
		mon = d[i,:][d[i,:] .!= 0]
		if length(mon) == length(union(mon))
			idx_mon[i] = sort(mon)
		end
	end
	d2i = Dict{Vector{Int64},Int64}(sort(d[i,:]) => i for i in 1:size(d)[1]) # Returns the index of each element of 'd'.

	i2keep = Dict{Int64,Vector{Int64}}(i => Int64[] for i in 1:n) 
	for k in 1:length(keep)
		i,j = keep[k]
		push!(i2keep[i],k)
	end

	coeff = Dict{Int64,Matrix{Float64}}()
	ids = Dict{Int64,Vector{Int64}}(i => Int64[] for i in 1:n) # For each agent, contains the idinces of the relevant monomials (according to 'keep').
	err = 0.

	args = [(i,X,Y[[i,],:],i2keep[i],keep,d2i,λ,ρ,niter) for i in 1:n]
	pmap(mySINDy_par_filter,args)

	for i in 1:n
		coeff[i] = readdlm("data/temp-sindy-Xi-$i.csv",',')
		rm("data/temp-sindy-Xi-$i.csv")
		ids[i] = Int64.(vec(readdlm("data/temp-sindy-id-$i.csv",',')))
		rm("data/temp-sindy-id-$i.csv")
		err += readdlm("data/temp-sindy-err-$i.csv",',')[1]
		rm("data/temp-sindy-err-$i.csv")
	end

	relerr = err/norm(Y,1)

	@info "THIS completed."
	#############################################################################
	
	# Reconstructing the adjacency tensors ######################################
	Ainf = Dict{Int64,Matrix{Float64}}(o => zeros(0,o+1) for o in 1:dmax+1) # For each order o, associates a dictionary associating the pair (agents, hyperedge) to the inferred weight, as seen from the agent.
	for i in 1:n
		for j in 1:length(ids[i])
			id = ids[i][j]
			jj = idx_mon[id]
			o = length(jj)+1
			Ainf[o] = vcat(Ainf[o],[i jj' coeff[i][j]])
		end
	end

	@info "Dictionary of inferred adjacency tensors built."
	#############################################################################

	return Ainf, coeff, ids, relerr
end


# ================================================================================
"""
	mySINDy_par(arg::Tuple{Int64,Matrix{Float64},Matrix{Float64},Int64,Float64,Float64,Int64})

Parallel version of mySINDy, to be used without filtering!

_INPUT_ (in arg):\\
`i`: Index of the time series analyzed.\\
`X`: Time series of the agents' states.\\
`Y`: Time series of the time derivative of agent `i`.\\
`dmax`: Maximal degree to consider in the monomial library.\\
`λ`: SINDy's threshold deciding whether an hyperedge exists or not.\\
`ρ`: Regularization parameters promoting sparsity.\\
`niter`: Maximal number of iterations for THIS algorithm.

_OUTPUT_ (written in a file):\\
`Ξ`: Matrix of coefficient inferred by SINDy. The index of each row is the index of an agent and the the columns corresponds to elements of the Taylor basis.\\
`err`: 1-norm of the difference between `Y` and the inferred dynamics.
"""
function my_SINDy_par(arg::Tuple{Int64,Matrix{Float64},Matrix{Float64},Int64,Float64,Float64,Int64})
	i,X,Y,dmax,λ,ρ,niter = arg

	n,T = size(X)
	θ,d = get_θd(X,dmax)

	k = 1
	Ξ = Y*θ'*inv(θ*θ' + ρ*Id(m))
	while k < niter
		k += 1
		smallinds = (abs.(Ξ) .< λ)
		Ξ[smallinds] .= 0.
		biginds = (.~smallinds)[1,:]
		Ξ[1,biginds] = Y*θ[biginds,:]'*inv(θ[biginds,:]*θ[biginds,:]' + ρ*Id(sum(biginds)))
	end
	err = norm(Y - Ξ*θ,1)

	writedlm("data/temp-sindy-Xi-$i.csv",Ξ,',')
	writedlm("data/temp-sindy-err-$i.csv",err,',')
end

"""
	mySINDy_par_filter(arg::Tuple{Int64,Matrix{Float64},Matrix{Float64},Vector{Int64},Vector{Vector{Int64}},Dict{Vector{Int64},Int64},Float64,Float64,Int64})

Parallel version of mySINDy, to be used with filtering!

_INPUT_ (in arg):\\
`i`: Index of the time series analyzed.\\
`X`: Time series of the agents' states.\\
`Y`: Time series of the time derivative of agent `i`.\\
`i2keepi`: Filtering parameter.\\
`keep`: List of pairs of agents to be considered.\\
`d2i`: Filtering parameter.\\
`λ`: SINDy's threshold deciding whether an hyperedge exists or not.\\
`ρ`: Regularization parameters promoting sparsity.\\
`niter`: Maximal number of iterations for THIS algorithm.

_OUTPUT_ (written in a file):\\
`Ξ`: Matrix of coefficient inferred by SINDy. The index of each row is the index of an agent and the the columns corresponds to elements of the Taylor basis.\\
`ids`: List of monomials indices that are present in `Ξ`.\\
`err`: 1-norm of the difference between `Y` and the inferred dynamics.
"""
function mySINDy_par_filter(arg::Tuple{Int64,Matrix{Float64},Matrix{Float64},Vector{Int64},Vector{Vector{Int64}},Dict{Vector{Int64},Int64},Float64,Float64,Int64})
	i,X,Y,i2keepi,keep,d2i,λ,ρ,niter = arg

	n,T = size(X)
	θ = ones(1,T)
	id = [1,]
	for k in i2keepi
		j = setdiff(keep[k],[i,])[1]
		θ = vcat(θ,reduce(vcat,[prod(X[v,:],dims=1) for v in [setdiff([j,l],[0,]) for l in 0:n]]))
		append!(id,[d2i[v] for v in [sort([j,l]) for l in 0:n]])
	end

	iii = unique(a -> id[a], eachindex(id))
	id = id[iii]
	θ = θ[iii,:]
	m,T= size(θ)

	@info "Running mySINDy for agent $i"
	if isempty(id)
		coeff = zeros(0,0)
		err = norm(Y,1)
	else
		k = 1
		Ξ = Y*θ'*inv(θ*θ' + ρ*Id(m))
		while k < niter
			k += 1
			smallinds = (abs.(Ξ) .< λ)
			Ξ[smallinds] .= 0.
			biginds = (.~smallinds)[1,:]
			Ξ[1,biginds] = Y*θ[biginds,:]'*inv(θ[biginds,:]*θ[biginds,:]' + ρ*Id(sum(biginds)))
		end
		err = norm(Y - Ξ*θ,1)
	end

	writedlm("data/temp-sindy-Xi-$i.csv",Ξ,',')
	writedlm("data/temp-sindy-id-$i.csv",id,',')
	writedlm("data/temp-sindy-err-$i.csv",err,',')
end




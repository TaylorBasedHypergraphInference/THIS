include("this-tools.jl")

# ================================================================================
"""
	this(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)

Infers the hypergraph underlying the dynamics of its vertices with knowledge of the states 'X' and of the derivatives 'Y' at each vertex.

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
function this(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	# Size of the data matrix. 
	n,T = size(X)

	if size(X) != size(Y)
		@info "Dimensions of states and derivatives do not match."
		return nothing
	end

	# Running THIS ##############################################################
	# Retrieve the values of the monomials at each time step.
	θ,d = get_θd(X,dmax)
	idx_mon = Dict{Int64,Vector{Int64}}()
	for i in 1:size(d)[1]
		mon = d[i,:][d[i,:] .!= 0]
		if length(mon) == length(union(mon))
			idx_mon[i] = sort(mon)
		end
	end

	# Running SINDy
	coeff, err = mySINDy(θ,Y,λ,ρ,niter)
	relerr = err/norm(Y,1)

	@info "THIS completed."
	#############################################################################

	# Reconstructing the adjacency tensors ######################################
	Ainf = Dict{Int64,Matrix{Float64}}(o => zeros(0,o+1) for o in 1:dmax+1) # For each order o, associates a dictionary associating the pair (agents, hyperedge) to the inferred weight, as seen from the agent.
	idx_coeff = Dict{Int64,Vector{Int64}}()
	nz_idx = Int64[]
	for i in keys(idx_mon)
		aaa = setdiff((1:n)[abs.(coeff[:,i]) .> 1e-8],idx_mon[i])
		if !isempty(aaa)
			push!(nz_idx,i)
			idx_coeff[i] = aaa
		end
	end

	for id in nz_idx
		ii = idx_coeff[id]
		jj = idx_mon[id]
		o = length(jj)+1
		Ainf[o] = vcat(Ainf[o],[ii repeat(jj',length(ii),1) coeff[ii,id]])
	end

	@info "Dictionary of inferred adjacency tensors built."
	#############################################################################

	return Ainf, coeff, relerr
end

# ================================================================================
"""
	this_filter(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, par_keep::Union{Int64,Float64}, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)

Infers the hypergraph underlying the dynamics of its vertices with knowledge of the states 'X' and of the derivatives 'Y' at each vertex.

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
function this_filter(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, par_keep::Union{Int64,Float64}, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	if size(X) != size(Y)
		@info "Dimensions of states and derivatives do not match."
		return nothing
	end
	n,T = size(X)

	# Listing pairs of agents whose trajectories have a correlation coefficient above 'α'.
	keep = keep_correlated(X,par_keep)

	# Running THIS ##############################################################
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
	ids = Dict{Int64,Vector{Int64}}(i => Int64[] for i in 1:n) # For each agent, contains the indices of the relevant monomials (according to 'keep').
	err = 0.
	energy = 0.
	for i in 1:n
		θ = ones(1,T)
		id = [1,]
		for k in i2keep[i]
			j = setdiff(keep[k],[i,])[1]
			θ = vcat(θ,reduce(vcat,[prod(X[v,:],dims=1) for v in [setdiff([j,l],[0,]) for l in 0:n]]))
			append!(id,[d2i[v] for v in [sort([j,l]) for l in 0:n]])
		end

		iii = unique(a -> id[a], eachindex(id))
		id = id[iii]
		θ = θ[iii,:]

		@info "Running mySINDy for agent $i"
		if isempty(id)
			coef = zeros(0,0)
			er = norm(Y[i,:],1)
		else
			coef,er = mySINDy(θ,Y[[i,],:],λ,ρ,niter)
			ids[i] = id
		end

		coeff[i] = coef
		err += er
		energy += norm(Y[i,:],1)
	end

	relerr = err/energy

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
	mySINDy(θ::Matrix{Float64}, Y::Matrix{Float64}, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)

Own implementation of SINDy. Adapted from [https://github.com/eurika-kaiser/SINDY-MPC/blob/master/utils/sparsifyDynamics.m], accessed on December 27, 2023.

_INPUT_:\\
`θ`: Values of the basis functions at each time steps. Here the basis functions are the monomials up to degree `dmax`.\\
`Y`: Time series of the system's velocity. Each row is the time series of the velocity of on agent.\\
`λ`: SINDy's threshold deciding whether an hyperedge exists or not.\\
`ρ`: Regularization parameters promoting sparsity.\\
`niter`: Maximal number of iterations for THIS algorithm.

_OUTPUT_:\\
`Ξ`: Matrix of coefficient inferred by SINDy. The index of each row is the index of an agent and the the columns corresponds to elements of the Taylor basis.\\
`err`: 1-norm of the difference between `Y` and the inferred dynamics.
"""
function mySINDy(θ::Matrix{Float64}, Y::Matrix{Float64}, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	# Sizes of the data matrices.
	n,T = size(Y)
	m,T = size(θ)

	# Initialization of the coefficient matrix.
	Ξ = Y*θ'*inv(θ*θ' + ρ*Id(m)) # Least square with Tikhonov regularization
	
	nz = 1e6 # Number of nonzero elements in Ξ.
	k = 1 # Counter of iterations.

	while k < niter && sum(abs.(Ξ) .> 1e-6) != nz
		k += 1
		nz = sum(abs.(Ξ) .> 1e-6)

		# Getting rid of the small (< λ) components and re-optimizing.
		smallinds = (abs.(Ξ) .< λ)
		Ξ[smallinds] .= 0.
		for i in 1:n
			biginds = .~smallinds[i,:]
			Ξ[i,biginds] = Y[[i,],:]*θ[biginds,:]'*inv(θ[biginds,:]*θ[biginds,:]' + ρ*Id(sum(biginds)))
		end
	end

	# Absolute error
	err = norm(Y - Ξ*θ,1)

	return Ξ, err
end



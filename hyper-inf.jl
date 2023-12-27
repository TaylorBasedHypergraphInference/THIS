include("hyper-inf-tools.jl")

# ================================================================================
"""
	hyper_inf(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	hyper_inf(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, d::Int64, λ::Float64=.1, ρ::Float64=1.)

Infers the hypergraph underlying the dynamics of its vertices with knowledge of the states 'X' and of the derivatives 'Y' at each vertex.

_INPUT_:\\
`X`: Time series of the system's state. Each row is the time series of the state of one agent.\\
`Y`: Time series of the system's velocity. Each row is the time series of the velocity of on agent.\\
`ooi`: Orders of interest. Vector of integers listing the orders of interactions that we analyze.\\
`dmax`: Maximal degree to be considered in the Taylor expansion. Typically, dmax=maximum(ooi).\\
[`d`: Internal dimension of the agents. It is assmued that the components of an agents are consecutive in the state vector.]\\
`λ`: SINDy's threshold deciding whether an hyperedge exists or not.\\
`ρ`: Regularization parameters promoting sparsity.\\
`niter`: Maximal number of iterations for THIS algorithm.

_OUTPUT_:\\
`Ainf`: Dictionary associating the inferred coefficient of the Taylor expansion to a pair (node,hyperedge). Namely, 'Ainf[(i,h)]' is the coefficient corresponding to the hyperedge 'h' in the Taylor expansion of the dynamics of node 'i'. If agents have internal dimensions larger than 1, `Ainf` is just a boolean.\\
[`AAinf`: Contains the inferred coefficients for edges in systems where agents have internal dimension larger than 1.]\\
`coeff`: Matrix of coefficents obtained by SINDy. The row indices are the agents' indices and the columns indices are the indices of the monomials in the Taylor series.\\
`err`: Absolute error, i.e., squared difference between `Y` and the inferred corresponding time series.\\
`relerr`: Relative error, i.e., `err` normalized by the squared magnitude of `Y`.
"""
function hyper_inf(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	# Size of the data matrix. 
	n,T = size(X)

	if size(X) != size(Y)
		@info "Dimensions of states and derivatives do not match."
		return nothing
	end

	# Defining the basis of functions to use, i.e., the monomials up to order 'dmax'.
	@variables x[1:n]
	prebasis = polynomial_basis([x[i] for i in 1:n],dmax)
	basis = Basis(prebasis,[x[i] for i in 1:n])

	# Perform THIS
	coeff,err,relerr = this(X,Y,ooi,dmax,λ,ρ,niter)

	# Retrieving the results of SINDy and doing the inference.
	idx_o = Dict{Int64,Vector{Int64}}() # For each order o, contains the indices of the basis elements in x corresponding to hyperedges of order o and involving distinct agents. 
	agents_o = Dict{Int64,Vector{Vector{Int64}}}() # Lists the sets of agents involved in each hyperedge from idx_o. 
	Ainf = Dict{Int64,Any}() # For each order o, associates a dictionary associating the pair (agents, hyperedge) to the inferred weight, as seen from the agent.
	for o in sort(ooi,rev=true)
		Ainf[o] = Dict{Tuple{Int64,Vector{Int64}},Float64}() # Inferred hyperedges of order o
		idx_o[o],agents_o[o] = get_idx_o(o-1,x,prebasis)
		for k in 1:length(idx_o[o])
			id = idx_o[o][k]
			agents = agents_o[o][k]
			cagents = setdiff(1:n,agents)
			y = coeff[cagents,id]
			ynz = y[Int64.(setdiff((1:length(y)),[0.,]))]
			yds = Int64.(setdiff(cagents,[0,]))
			for j in 1:length(yds)
				Ainf[o][(yds[j],sort([yds[j];agents]))] = ynz[j]
			end
		end
	end

	return Ainf, coeff, err, relerr
end

function hyper_inf(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Union{Int64,Vector{Int64}}, dmax::Int64, d::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	Ainf, coeff, err, relerr = hyper_inf(X,Y,ooi,dmax,λ,ρ)

	Ainf,AAinf = one2dim(Ainf,d)

	return Ainf, AAinf, coeff, err, relerr
end

# ================================================================================
"""
	this(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)

Runs THIS algorithm, returning the matrix of inferred coefficients, and the absolute and relative errors.

_INPUT_:\\
`X`: Time series of the system's state. Each row is the time series of the state of one agent.\\
`Y`: Time series of the system's velocity. Each row is the time series of the velocity of on agent.\\
`ooi`: Orders of interest. Vector of integers listing the orders of interactions that we analyze.\\
`dmax`: Maximal degree to be considered in the Taylor expansion. Typically, dmax=maximum(ooi).\\
`λ`: SINDy's threshold deciding whether an hyperedge exists or not.\\
`ρ`: Regularization parameters promoting sparsity.\\
`niter`: Maximal number of iterations for THIS algorithm.

_OUTPUT_:\\
`coeff`: Matrix of coefficient inferred by SINDy. The index of each row is the index of an agent and the the columns corresponds to elements of the Taylor basis.\\
`err`: Squared difference between the inferred time series for `Y` and `Y` itself.\\
`relerr`: `err` normalized by the squared magnitude of `Y`.
"""
function this(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	if size(X) != size(Y)
		@info "Time series' sizes do not match."
	end

	# Retrieve the values of the monomials at each time step.
	θ = get_θ(X,dmax)

	return mySINDy(θ,Y,λ,ρ,niter)
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
`coeff`: Matrix of coefficient inferred by SINDy. The index of each row is the index of an agent and the the columns corresponds to elements of the Taylor basis.\\
`err`: Squared difference between the inferred time series for `Y` and `Y` itself.\\
`relerr`: `err` normalized by the squared magnitude of `Y`.
"""
function mySINDy(θ::Matrix{Float64}, Y::Matrix{Float64}, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)
	# Sizes of the data matrices.
	n,T = size(Y)
	m,T = size(θ)
	# Squared magnitude of Y.
	energy = sum(Y.^2)

	# Initialization of the coefficient matrix.
	Ξ = Y*θ'*pinv(θ*θ' + ρ*Id(m)) # Least square with Tikhonov regularization
	
	nz = 1e6 # Number of nonzero elements in Ξ.
	k = 1 # Counter of iterations.

	while k < niter && sum(abs.(Ξ) .> 1e-6) != nz
		k += 1
		nz = sum(abs.(Ξ) .> 1e-6)
		err = sum((Y - Ξ*θ).^2) # Absolute error
		@info "iter $k: $(nz) nonzero coefficients, error = $(round(err)), rel-err = $(round(err/energy*100,digits=2))%"

		# Getting rid of the small (< λ) components and re-optimizing.
		smallinds = (abs.(Ξ) .< λ)
		Ξ[smallinds] .= 0.
		for i in 1:n
			biginds = .~smallinds[i,:]
			Ξ[i,biginds] = Y[[i,],:]*θ[biginds,:]'*pinv(θ[biginds,:]*θ[biginds,:]' + ρ*Id(sum(biginds)))
		end
	end

	# Absolute error
	err = sum((Y - Ξ*θ).^2)

	return Ξ, err, err/energy
end



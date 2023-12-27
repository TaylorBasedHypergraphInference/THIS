using DataDrivenDiffEq, ModelingToolkit, LinearAlgebra, DataDrivenSparse, LinearAlgebra, PyPlot, Combinatorics, Statistics

# ================================================================================
"""
	get_idx_o(o::Int64, x::Symbolics.Arr{Num,1}, prebasis::Vector{Num})

Retrieves the indices (in 'prebasis') of the monomials of order 'o' in the variables 'x', involving distincts agents.

_INPUT_:\\
`o`: Order of the monomials.\\
`x`: Symbolic variables of the agents' states.\\
`prebasis`: List of equations defining the monomials in terms of the variables in `x`.

_OUTPUT_:\\
`idx`: List of indices of the elements of `prebasis` corresponding to hyperedges of order `o`, involving distinct agents.\\
`agents`: List of the agents involved in each of the above hyperedges.
"""
function get_idx_o(o::Int64, x::Symbolics.Arr{Num,1}, prebasis::Vector{Num})
	comb = combinations(1:length(x),o)
	idx = Int64[]
	agents = Vector{Int64}[]
	for c in comb
		mon = get_monomial(x,c)
		m,i = findmax(isequal.(prebasis,mon))
		if m > .1
			push!(idx,i)
			push!(agents,c)
		end
	end
	return idx, agents
end

# ================================================================================
"""
	get_monomial(x::Symbolics:Arr{Num,1}, c::Vector{Int64})

Constructing the monomial of variables in `x` with indices of `c`.

_INPUT_:\\
`x`: Symbolic variables for the states of agents.\\
`c`: Indices of the variables involved in the monomial.

_OUTPUT_:\\
`mon`: Symbolic monomial.
"""
function get_monomial(x::Symbolics.Arr{Num,1}, c::Vector{Int64})
	mon = 1
	for i in c
		mon *= x[i]
	end

	return mon
end

# ================================================================================
"""
	get_monomials(x::Symbolics:Arr{Num,1}, o::Int64)

Constructing the list of monomial of order `o` with variables in `x`.

_INPUT_:\\
`x`: Symbolic variables for the states of agents.\\
`o`: Order of the monomials.

_OUTPUT_:\\
`mon`: Symbolic monomials.
"""
function get_monomials(x::Symbolics.Arr{Num,1}, o::Int64)
    n = length(x)
    mon = Num[]
    comb = collect(combinations(1:n,o))

    for c in combi
	    push!(mon,get_monomial(x,c))
    end

    return mon
end

# ================================================================================
"""
	get_θ(X::Matrix{Float64}, dmax::Int64)
	get_theta(X::Matrix{Float64}, dmax::Int64)

Returns the matrix of values of the monomials up to degree `dmax` for each time step of `X`.

_INPUT_:\\
`X`: Time series of the agents' states. Rows indices are the agents' indices and the columns indices are the time steps.\\
`dmax`: Maximal monomial degree to be considered. 

_OUTPUT_:\\
`θ`: Matrix of the monomial values. Rows indices are the monomial indices and columns indices are the time steps.
"""
function get_θ(X::Matrix{Float64}, dmax::Int64)
	n,T = size(X)

	@variables x[1:n]
	prebasis = polynomial_basis([x[i] for i in 1:n],dmax)
	basis = Basis(prebasis,[x[i] for i in 1:n])

	θ = zeros(length(basis),0)
	for t in 1:T
		θ = [θ basis(X[:,t])]
	end

	return θ
end

function get_theta(X::Matrix{Float64}, dmax::Int64)
	return get_θ(X,dmax)
end

# ================================================================================
"""
	Id(n::Int64)

Returns the n-dimensional identity matrix.

_INPUT_:\\
`n`: Size of the matrix.

_OUTPUT_:\\
`Id`: Identity matrix of size `n`.
"""
function Id(n::Int64)
	return diagm(0 => ones(n))
end

# ================================================================================
"""
	inferred_adj_2nd(Ainf::Dict{Tuple{Int64,Vector{Int64}},Float64}, n::Int64, thr::Float64=0.)

Returns the inferred 2nd-order adjacency matrix from `Ainf`. 

_INPUT_:\\
`Ainf`: Dictionary of the inferred hyperedges (output of `hyper_inf`).\\
`n`: Number of vertices.\\
`thr`: Threshold below which hyperedges are discarded.

_OUTPUT_:\\
`At_bool`: Boolean adjacency matrix.\\
`At_float`: Weighted adjacency matrix.
"""
function inferred_adj_2nd(Ainf::Dict{Tuple{Int64,Vector{Int64}},Float64}, n::Int64, thr::Float64=0.)
	At_bool = zeros(n,n)
	At_float = zeros(n,n)

	for k in keys(Ainf)
		i = k[1]
		ic = setdiff(k[2],[i,])
		for p in permutations(ic)
			At_bool[i,p[1]] = (abs.(Ainf[k]) > thr)
			At_float[i,p[1]] = Ainf[k]
		end
	end

	return At_bool, At_float
end

"""
	inferred_adj_3rd(Ainf::Dict{Tuple{Int64,Vector{Int64}},Float64}, n::Int64, thr::Float64=0.)

Returns the inferred 3rd-order adjacency tensor from `Ainf`. 

_INPUT_:\\
`Ainf`: Dictionary of the inferred hyperedges (output of `hyper_inf`).\\
`n`: Number of vertices.\\
`thr`: Threshold below which hyperedges are discarded.

_OUTPUT_:\\
`At_bool`: Boolean adjacency tensor.\\
`At_float`: Weighted adjacency tensor.
"""
function inferred_adj_3rd(Ainf::Dict{Tuple{Int64,Vector{Int64}},Float64}, n::Int64, thr::Float64=0.)
	At_bool = zeros(n,n,n)
	At_float = zeros(n,n,n)

	for k in keys(Ainf)
		i = k[1]
		ic = setdiff(k[2],[i,])
		for p in permutations(ic)
			At_bool[i,p[1],p[2]] = (abs.(Ainf[k]) > thr)
			At_float[i,p[1],p[2]] = Ainf[k]
		end
	end

	return At_bool, At_float
end

"""
	inferred_adj_4th(Ainf::Dict{Tuple{Int64,Vector{Int64}},Float64}, n::Int64, thr::Float64=0.)

Returns the inferred 4th-order adjacency tensor from `Ainf`.

_INPUT_:\\
`Ainf`: Dictionary of the inferred hyperedges (output of `hyper_inf`).\\
`n`: Number of vertices.\\
`thr`: Threshold below which hyperedges are discarded.

_OUTPUT_:\\
`At_bool`: Boolean adjacency tensor.\\
`At_float`: Weighted adjacency tensor.
"""
function inferred_adj_4th(Ainf::Dict{Tuple{Int64,Vector{Int64}},Float64}, n::Int64, thr::Float64=0.)
	At_bool = zeros(n,n,n,n)
	At_float = zeros(n,n,n,n)

	for k in keys(Ainf)
		i = k[1]
		ic = setdiff(k[2],[i,])
		for p in permutations(ic)
			At_bool[i,p[1],p[2],p[3]] = (abs.(Ainf[k]) > thr)
			At_float[i,p[1],p[2],p[3]] = Ainf[k]
		end
	end

	return At_bool, At_float
end

# ================================================================================
"""
	one2dim(Ainf::Dict{Int64,Any}, d::Int64=1)

Transforms the hypergraph inferred for `n`*`d` agents with a single internal dimension to a hypergraph of `n` agents with `d` internal dimension.

_INPUT_:\\
`Ainf`: Dictionary of inferred hyperedges (output of `hyper_inf`).\\
`d`: Internal dimension of agents.

_OUTPUT_:\\
`A`: For each hyperedge order `o`, contains a dictionary associating pairs (vertex,hyperedge) to a boolean, indicating the existence of the inferred hyperedge, as viewed from the given vertex.\\
`AA`: For each hyperedge order `o`, contains a dictionary associating the triplet `(j,u,w)` to the corresponding inferred hyperedge weight. The triplet is composed as follows: `j` is the agent' index, `u` is the hyperedge of interaction, `w` indicates which degree of freedom of the agents interact through this hyperedge.
"""
function one2dim(Ainf::Dict{Int64,Any}, d::Int64=1)
	ooi = sort(keys(Ainf))
	A = Dict{Int64,Any}()
	AA = Dict{Int64,Any}()
	for o in ooi
		A[o] = Dict{Tuple{Int64,Vector{Int64}},Float64}()
		AA[o] = Dict{Tuple{Int64,Vector{Int64},Vector{Int64}},Float64}() # Description of the keys: 1. index of the source node of the interaction, 2. list of the interacting agents, 3. list of the components interacting.
		for k in keys(Ainf[o])
			i,v = k
			j = ceil(Int64,i/d)
			u = ceil.(Int64,v./d)
			w = (v .- 1).%d .+ 1
			
			if length(union(u)) == length(u)
				A[o][(j,u)] = 1.
				AA[o][(j,u,w)] = Ainf[o][k]
			end
		end
	end

	return A, AA
end


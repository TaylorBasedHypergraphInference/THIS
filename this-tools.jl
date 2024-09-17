using LinearAlgebra, Combinatorics, Statistics

# ================================================================================
"""
	get_θd(X::Matrix{Float64}, dmax::Int64, i0::Int64=1)
	get_thetad(X::Matrix{Float64}, dmax::Int64, i0::Int64=1)

Returns the matrix of values of the monomials up to degree `dmax` for each time step of `X`. As well as the list of the agents involved in each of the monomials (mostly for indexing purpose).

	get_θ(X::Matrix{Float64}, dmax::Int64, i0::Int64=1)
	get_theta(X::Matrix{Float64}, dmax::Int64, i0::Int64=1)

Returns only the monomials.

	get_d(n::Int64, dmax::Int64, i0::Int64=1)

Returns only the list of agents involved.

_INPUT_:\\
`X`: Time series of the agents' states. Rows indices are the agents' indices and the columns indices are the time steps.\\
`dmax`: Maximal monomial degree to be considered.\\ 
`i0`: Index of the first agent to consider in the list (mostly for recursive calling purpose).

_OUTPUT_:\\
`θ`: Matrix of the monomial values. Rows indices are the monomial indices and columns indices are the time steps.\\
`d`: List of the agents involved in each monomial of `θ`. Each row has `dmax` elements. If one element appears multiple times, it means that its degree is larger that 1 in the monomial. If there are some zeros in the row, it means that the monomial has degree smaller than `dmax`. 
"""
function get_θd(X::Matrix{Float64}, dmax::Int64, i0::Int64=1)
	n,T = size(X)

	θ = ones(1,T)
	d = zeros(Int64,1,dmax)

	if dmax == 0
		return θ,d
	else
		for i in 1:n
			θ0,d0 = get_θd(X[i:n,:],dmax-1,i)
			θ = vcat(θ,repeat(X[[i,],:],size(θ0)[1],1).*θ0)
			d = vcat(d,[d0 i*ones(Int64,size(d0)[1],1)])
		end

		d += (i0-1)*(d .> 0)

		return θ,d
	end
end

function get_thetad(X::Matrix{Float64}, dmax::Int64, i0::Int64=1)
	return get_θd(X,dmax,i0)
end

function get_θ(X::Matrix{Float64}, dmax::Int64, i0::Int64=1)
	n,T = size(X)

	θ = ones(1,T)

	if dmax == 0
		return θ
	else
		for i in 1:n
			θ0 = get_θ(X[i:n,:],dmax-1,i)
			θ = vcat(θ,repeat(X[[i,],:],size(θ0)[1],1).*θ0)
		end

		return θ
	end
end

function get_theta(X::Matrix{Float64}, dmax::Int64, i0::Int64=1)
	return get_θ(X,dmax,i0)
end

function get_d(n::Int64, dmax::Int64, i0::Int64=1)
	d = zeros(Int64,1,dmax)

	if dmax == 0
		return d
	else
		for i in 1:n
			d0 = get_d(n-i+1,dmax-1,i)
			d = vcat(d,[d0 i*ones(Int64,size(d0)[1],1)])
		end

		d += (i0-1)*(d .> 0)

		return d
	end
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
	keep_correlated(X::Matrix{Flaot64}, α::Float64)

Returns the list of pairs of time series (rows of `X`) whose Pearson correlation factor is larger than `α`.

	keep_correlated(X::Matrix{Float64}, k0::Int64)

Returns the list of the `k0` pairs of time series (rows of `X`) that have the largest Pearson correlation. 

_INPUT_:\\
`X`: Time series with rows indexing the degrees of freedom and the columns indexing the time steps.\\
`α`: Correlation threshold above which we keep the pairs of time series.\\
`k0`: Number of pairs of time series that are kept.

_OUTPUT_:\\
`keep`: List of the pairs of time series that are kept.
"""
function keep_correlated(X::Matrix{Float64}, α::Float64)
	C = cor(X')
	n = size(C)[1]

	c = [abs(C[i,j]) for i in 1:n-1 for j in i+1:n]
	ids = [[i,j] for i in 1:n-1 for j in i+1:n]
	keep = ids[c .> α]

	@info "Filtering: $(length(keep))/$(Int(n*(n-1)/2)) pairs kept."

	return keep
end

function keep_correlated(X::Matrix{Float64}, k0::Int64)
	C = cor(X')
	n = size(C)[1]

	k = min(k0,Int64(n*(n-1)/2))

	c = [abs(C[i,j]) for i in 1:n-1 for j in i+1:n]
	ids = [i for i in 1:n-1 for j in i+1:n]
	jds = [j for i in 1:n-1 for j in i+1:n]
	sorted = sortslices([c ids jds],dims=1,rev=true)

	keep = [[Int64(sorted[l,2]),Int64(sorted[l,3])] for l in 1:k]

	return keep
end

# ================================================================================
"""
	one2dim(Ainf::Dict{Int64,Matrix{Float64}}, d::Int64=1)

Transforms the hypergraph inferred for `n`*`d` agents with a single internal dimension to a hypergraph of `n` agents with `d` internal dimension.

_INPUT_:\\
`Ainf`: Dictionary of the inferred adjacency lists (output of `hyper_inf`).\\
`d`: Internal dimension of agents.

_OUTPUT_:\\
`AAinf`: Dictionary of the inferred adjacency lists with aggregted node dimensions.
"""
function one2dim(Ainf::Dict{Int64,Matrix{Float64}}, d::Int64=1)
	ooi = sort(keys(Ainf))
	AAinf = Dict{Int64,Matrix{Float64}}(o => zeros(0,o+1) for o in ooi)
	for o in ooi
		L,x = size(Ainf[o])
		list = Vector{Int64}[]
		for l in 1:L
			ii = ceil.(Int64,Ainf[o][l,1:o]./d)
			a = Ainf[o][l,o+1]
			if ii in list
				j = findall(x -> isequal(x,ii),list)[1]
				AAinf[o][j,end] = max(abs(a),AAinf[o][j,end])
			else
				AAinf[o] = vcat(AAinf[o],[ii' abs(a)])
				push!(list,ii)
			end
		end
	end

	return AAinf
end





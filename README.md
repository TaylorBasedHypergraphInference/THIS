# THIS: Taylor-based Hypergraph Inference using SINDy
Companion codes for the manuscript *Delabays, De Pasquale, D&ouml;rfler, and Zhang (2024)*.

## Summary of the files
- **eeg-data**: Contains the necessary material for hyperedge inference from EEG time series.
- **eeg-run.jl**: Runs THIS on EEG time series. By default, the time series are assumed to be save in **eeg-data** as *.edf* files. As long as all files are correctly provided, this script can be run with *"include("eeg-run.jl")* in a Julia terminal.
- **eeg-tools.jl**: Contains all side functions necessary to run THIS on EEG time series. 
- **hyper-inf.jl**: Contains the core functions to run THIS. 
- **hyper-inf-tools.jl**: Contains all the side functions necessary to run THIS.

## List of functions
Some functions are associated to the inference problem in general, while others are dedicated to the implementation of THIS to the hypergraph inference from EEG data. 

### Functions for inference
- [get\_idx\_o](#get\_idx\_o)
- [get\_monomial](#get\_monomial)
- [get\_monomials](#get\_monomials)
- [get\_θ](#get\_θ)
- [hyper\_inf](#hyper\_inf)
- [Id](#Id)
- [inferred\_adj\_nth](#inferred\_adj\_nth)
- [mySINDy](#mySINDy)
- [this](#this)
 
### Functions for EEG data
- [average\_over\_zones](#average\_over\_zones)
- [denoise\_fourier](#denoise\_fourier)
- [list\_all\_subjects](#list\_all\_subjects)
- [load\_and\_save\_eeg\_avg](#load\_and\_save\_eeg\_avg)
- [read\_eeg](#read\_eeg)
- [restrict\_box\_size](#restrict\_box\_size)

## Detailed documentation for the functions

### average\_over\_zones
*./eeg-tools.jl*

- `average\_over\_zones(s2signal::Dict{String,Vector{Float64}}, s2z::Dict{String,Int64})`

Averages the time series of all the sensors within each zone. Average is taken at each time step independently. 

**INPUT**:\
`s2signal`: Contains the time series (output of `read_eeg`).\
`s2z`: Indicates how sensors are gathered into zones. Each sensor is associated to a zone.

**OUTPUT**:\
`asig`: Matrix of the time series. Row indices are zone indices and columns indices are time steps.

---

### denoise\_fourier
*./eeg-tools.jl*

- `denoise\_fourier(x::Vector{Float64}, nmodes::Int64)`
- `denoise\_fourier(X::Matrix{Float64}, nmodes::Int64)`

Gets rid of the noise in time series by low-pass filtering. Only the `nmodes` lowest modes of the FFT are kept.

**INPUT**:\
`x`, `X`: Time series. Either a single one or a matrix thereof.\
`nmodes`: Number of modes of the FFT to be kept.

**OUTPUT**:\
`y`, `Y`: Denoised time series. 

---

### get\_idx\_o
*./hyper-inf-tool.jl*

- `get\_idx\_o(o::Int64, x::Symbolics.Arr{Num,1}, prebasis::Vector{Num})`

Retrieves the indices (in 'prebasis') of the monomials of order 'o' in the variables 'x', involving distincts agents.

**INPUT**:\
`o`: Order of the monomials.\
`x`: Symbolic variables of the agents' states.\
`prebasis`: List of equations defining the monomials in terms of the variables in `x`.

**OUTPUT**:\
`idx`: List of indices of the elements of `prebasis` corresponding to hyperedges of order `o`, involving distinct agents.\
`agents`: List of the agents involved in each of the above hyperedges.

---

### get\_monomial
*./hyper-inf-tool.jl*
- `get\_monomial(x::Symbolics:Arr{Num,1}, c::Vector{Int64})`

Constructing the monomial of variables in `x` with indices of `c`.

**INPUT**:\
`x`: Symbolic variables for the states of agents.\
`c`: Indices of the variables involved in the monomial.

**OUTPUT**:\
`mon`: Symbolic monomial.

---

### get\_monomials
*./hyper-inf-tool.jl*
- `get\_monomials(x::Symbolics:Arr{Num,1}, o::Int64)`

Constructing the list of monomial of order `o` with variables in `x`.

**INPUT**:\
`x`: Symbolic variables for the states of agents.\
`o`: Order of the monomials.

**OUTPUT**:\
`mon`: Symbolic monomials.

---

### get\_θ
*./hyper-inf-tool.jl*

- `get\_θ(X::Matrix{Float64}, dmax::Int64)`

Returns the matrix of values of the monomials up to degree `dmax` for each time step of `X`.

**INPUT**:\
`X`: Time series of the agents' states. Rows indices are the agents' indices and the columns indices are the time steps.\
`dmax`: Maximal monomial degree to be considered. 

**OUTPUT**:\
`θ`: Matrix of the monomial values. Rows indices are the monomial indices and columns indices are the time steps.

### get\_theta
*./hyper-inf-tool.jl*

- `get\_theta(X::Matrix{Float64}, dmax::Int64)`

Same as `get\_θ`. 

---

### hyper\_inf
*./hyper-inf.jl*

- `hyper\_inf(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)`
- `hyper\_inf(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, d::Int64, λ::Float64=.1, ρ::Float64=1.)`

Infers the hypergraph underlying the dynamics of its vertices with knowledge of the states 'X' and of the derivatives 'Y' at each vertex.

**INPUT**:\
`X`: Time series of the system's state. Each row is the time series of the state of one agent.\
`Y`: Time series of the system's velocity. Each row is the time series of the velocity of on agent.\
`ooi`: Orders of interest. Vector of integers listing the orders of interactions that we analyze.\
`dmax`: Maximal degree to be considered in the Taylor expansion. Typically, dmax=maximum(ooi).\
[`d`: Internal dimension of the agents. It is assmued that the components of an agents are consecutive in the state vector.]\
`λ`: SINDy's threshold deciding whether an hyperedge exists or not.\
`ρ`: Regularization parameters promoting sparsity.\
`niter`: Maximal number of iterations for THIS algorithm.

**OUTPUT**:\
`Ainf`: Dictionary associating the inferred coefficient of the Taylor expansion to a pair (node,hyperedge). Namely, 'Ainf[(i,h)]' is the coefficient corresponding to the hyperedge 'h' in the Taylor expansion of the dynamics of node 'i'. If agents have internal dimensions larger than 1, `Ainf` is just a boolean.\
[`AAinf`: Contains the inferred coefficients for edges in systems where agents have internal dimension larger than 1.]\
`coeff`: Matrix of coefficents obtained by SINDy. The row indices are the agents' indices and the columns indices are the indices of the monomials in the Taylor series.\
`err`: Absolute error, i.e., squared difference between `Y` and the inferred corresponding time series.\
`relerr`: Relative error, i.e., `err` normalized by the squared magnitude of `Y`.

---

### Id
*./hyper-inf-tool.jl*

- `Id(n::Int64)`

Returns the n-dimensional identity matrix.

**INPUT**:\
`n`: Size of the matrix.

**OUTPUT**:\
`Id`: Identity matrix of size `n`.

---

### inferred\_adj\_nth
*./hyper-inf-tool.jl*

- `inferred\_adj\_2nd(Ainf::Dict{Tuple{Int64,Vector{Int64}},Float64}, n::Int64, thr::Float64=0.)`
- `inferred\_adj\_3rd(Ainf::Dict{Tuple{Int64,Vector{Int64}},Float64}, n::Int64, thr::Float64=0.)`
- `inferred\_adj\_4th(Ainf::Dict{Tuple{Int64,Vector{Int64}},Float64}, n::Int64, thr::Float64=0.)`

Returns the inferred nth-order adjacency tensor from `Ainf`.

**INPUT**:\
`Ainf`: Dictionary of the inferred hyperedges (output of `hyper\_inf`).\
`n`: Number of vertices.\
`thr`: Threshold below which hyperedges are discarded.

**OUTPUT**:\
`At_bool`: Boolean adjacency tensor.\
`At_float`: Weighted adjacency tensor.

---

### list\_all\_subjects
*./eeg-tools.jl*

- `list\_all\_subjects(n::Int64)`

Returns a `String` list of all subjects indices with three digits (from "001" up to "999"). 

**INPUT**:\
`n`: Number of subjects.

**OUTPUT**:\
`subjects`: List of subjects.

---

### load\_and\_save\_eeg\_avg
*./eeg-tools.jl*

- `load\_and\_save\_eeg\_avg(subjects::Vector{String}, states::Vector{String}, nz::Int64)`

Loads and saves the time series of `subjects` for tasks `states` in "eeg-data" as .csv files. Each time series is saved under the file name "SxxxRyy-X.csv".

**INPUT**:\
`subjects`: List of subjects whose time series will be saved.\
`states`: List of tasks' indices that will be saved.\
`nz`: Number of zones to gather time series.

---

### mySINDy
*./hyper-inf.jl*

- `mySINDy(θ::Matrix{Float64}, Y::Matrix{Float64}, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)`

Own implementation of SINDy. Adapted from [this link](https://github.com/eurika-kaiser/SINDY-MPC/blob/master/utils/sparsifyDynamics.m), accessed on December 27, 2023.

**INPUT**:\
`θ`: Values of the basis functions at each time steps. Here the basis functions are the monomials up to degree `dmax`.\
`Y`: Time series of the system's velocity. Each row is the time series of the velocity of on agent.\
`λ`: SINDy's threshold deciding whether an hyperedge exists or not.\
`ρ`: Regularization parameters promoting sparsity.\
`niter`: Maximal number of iterations for THIS algorithm.

**OUTPUT**:\
`coeff`: Matrix of coefficient inferred by SINDy. The index of each row is the index of an agent and the the columns corresponds to elements of the Taylor basis.\
`err`: Squared difference between the inferred time series for `Y` and `Y` itself.\
`relerr`: `err` normalized by the squared magnitude of `Y`.

---

### read\_eeg
*./eeg-tools.jl*

- `read\_eeg(file::String)`

Reads the .edf file `file` and returns its content as dictionary. 

**INPUT**:\
`file`: Name of the file to be loaded.

**OUTPUT**:\
`s2signal`: Dictionary associating the time series to each sensor name.

---

### restrict\_box\_size
*./eeg-tools.jl*

- `restrict\_box\_size(X::Matrix{Float64}, nstep::Int64)`

Restrict the "box" size (see manuscript), by keeping only the `nstep` time steps that are the closest to the median value of the time series. 

**INPUT**:\
`X`: Time series. Row indices are the agents time series and columns are the time steps.\
`nstep`: Number of time steps to be kept.

**OUTPUT**:\
`Y`: Relevant time steps.\
`ids`: Indices of the kept time steps.

---

### this
*./hyper-inf.jl*

- `this(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)`

Runs THIS algorithm, returning the matrix of inferred coefficients, and the absolute and relative errors.

**INPUT**:\
`X`: Time series of the system's state. Each row is the time series of the state of one agent.\
`Y`: Time series of the system's velocity. Each row is the time series of the velocity of on agent.\
`ooi`: Orders of interest. Vector of integers listing the orders of interactions that we analyze.\
`dmax`: Maximal degree to be considered in the Taylor expansion. Typically, dmax=maximum(ooi).\
`λ`: SINDy's threshold deciding whether an hyperedge exists or not.\
`ρ`: Regularization parameters promoting sparsity.\
`niter`: Maximal number of iterations for THIS algorithm.

**OUTPUT**:\
`coeff`: Matrix of coefficient inferred by SINDy. The index of each row is the index of an agent and the the columns corresponds to elements of the Taylor basis.\
`err`: Squared difference between the inferred time series for `Y` and `Y` itself.\
`relerr`: `err` normalized by the squared magnitude of `Y`.

---

# THIS: Taylor-based Hypergraph Inference using SINDy
Companion codes for the manuscript *Delabays, De Pasquale, D&ouml;rfler, and Zhang (2024)*.

**Please contact us if you notice any typo or bug!**

[![DOI](https://zenodo.org/badge/736235218.svg)](https://zenodo.org/doi/10.5281/zenodo.10530470)

## Summary of the files
- **eeg-data**: Contains the necessary material for hyperedge inference from EEG time series. Also stores the EDF data files of the EEG measurements. 
- **eeg-run.jl**: Runs THIS on EEG time series. By default, the time series are assumed to be save in **eeg-data** as *.edf* files. As long as all files are correctly provided, this script can be run with *"include("eeg-run.jl")* in a Julia terminal. Running the script for all 218 time series take approximately 60 seconds on a standard laptop.
- **eeg-tools.jl**: Contains all side functions necessary to run THIS on EEG time series. 
- **this.jl**: Contains the core functions to run THIS. 
- **this-parallel.jl**: Contains the core function to run THIS in parallel.
- **this-tools.jl**: Contains all the side functions necessary to run THIS.

## Prerequisite
The code has been developed with Julia 1.9. 
The necessary packages are *Combinatorics*, *DelimitedFiles*, *EDF*, *FFTW*, *LinearAlgebra*, and *Statistics*. 

## Hypergraph inference
The main function for inference is 'this', which performs the whole inference. 

## Hypergraph inference from EEG data
The script `eeg-run.jl` infers hyperedges from EEG data. By default, the time series are assumed to be save in **eeg-data** as *.edf* files. As long as all files are correctly provided, this script can be run with *include("eeg-run.jl")* in a Julia terminal.

The parameters that can be tuned in the script are:
- `nz` [line 6], the number of aggregating zones on the scalp. If the value is changed, the corresponding files `sensors-nz.csv` and `zones-nz.csv` should be provided. 
- `λ` and `ρ` [lines 9 and 10], the sparsity threshold and the regularization parameter in SINDy. 
- `subjects` [line 13], the list of subjects considered. Subjects should be labeled by a 3-digit string index `xxx` appearing in the data file name `SxxxRyy.edf`. 
- `states` [line 14], the list of states considered. States should be labeled by a 2-digit string index `yy` appearing in the data file name `SxxxRyy.ed`. 

## List of functions
Some functions are associated to the inference problem in general, while others are dedicated to the implementation of THIS to the hypergraph inference from EEG data. 

### Functions for inference
- [get\_d](#get\_θd)
- [get\_θ](#get\_θd)
- [get\_θd](#get\_θd)
- [Id](#Id)
- [keep\_correlated](#keep\_correlated)
- [mySINDy](#mySINDy)
- [mySINDy\_par](#mySINDy\_par)
- [mySINDy\_par\_filter](#mySINDy\_par\_filter)
- [one2dim](#one2dim)
- [this](#this)
- [this\_filter](#this\_filter)
- [this\_par](#this\_par)
- [this\_par\_filter](#this\_par\_filter)
 
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

- `average_over_zones(s2signal::Dict{String,Vector{Float64}}, s2z::Dict{String,Int64})`

Averages the time series of all the sensors within each zone. Average is taken at each time step independently. 

**INPUT**:\
`s2signal`: Contains the time series (output of `read_eeg`).\
`s2z`: Indicates how sensors are gathered into zones. Each sensor is associated to a zone.

**OUTPUT**:\
`asig`: Matrix of the time series. Row indices are zone indices and columns indices are time steps.

---

### denoise\_fourier
*./eeg-tools.jl*

- `denoise_fourier(x::Vector{Float64}, nmodes::Int64)`
- `denoise_fourier(X::Matrix{Float64}, nmodes::Int64)`

Gets rid of the noise in time series by low-pass filtering. Only the `nmodes` lowest modes of the FFT are kept.

**INPUT**:\
`x`, `X`: Time series. Either a single one or a matrix thereof.\
`nmodes`: Number of modes of the FFT to be kept.

**OUTPUT**:\
`y`, `Y`: Denoised time series. 

---

### get\_θd
*./this-tools.jl*

- `get\_θd(X::Matrix{Float64}, dmax::Int64, i0::Int64=1)`
- `get\_thetad(X::Matrix{Float64}, dmax::Int64, i0::Int64=1)`

Returns the matrix of values of the monomials up to degree `dmax` for each time step of `X`. As well as the list of the agents involved in each of the monomials (mostly for indexing purpose).

- `get\_θ(X::Matrix{Float64}, dmax::Int64, i0::Int64=1)`
- `get\_theta(X::Matrix{Float64}, dmax::Int64, i0::Int64=1)`

Returns only the monomials.

- `get\_d(n::Int64, dmax::Int64, i0::Int64=1)`

Returns only the list of agents involved.

**INPUT**:\
`X`: Time series of the agents' states. Rows indices are the agents' indices and the columns indices are the time steps.\
`dmax`: Maximal monomial degree to be considered. \
`i0`: Index of the first agent to consider in the list (mostly for recursive calling purpose).

**OUTPUT**:\
`θ`: Matrix of the monomial values. Rows indices are the monomial indices and columns indices are the time steps.\
`d`: List of the agents involved in each monomial of `θ`. Each row has `dmax` elements. If one element appears multiple times, it means that its degree is larger that 1 in the monomial. If there are some zeros in the row, it means that the monomial has degree smaller than `dmax`. 

---

### Id
*./this-tool.jl*

- `Id(n::Int64)`

Returns the n-dimensional identity matrix.

**INPUT**:\
`n`: Size of the matrix.

**OUTPUT**:\
`Id`: Identity matrix of size `n`.

---
### keep\_correlated
*./this-tools.jl*

- `keep\_correlated(X::Matrix{Flaot64}, α::Float64)`

Returns the list of pairs of time series (rows of `X`) whose Pearson correlation factor is larger than `α`.

- `keep\_correlated(X::Matrix{Float64}, k0::Int64)`

Returns the list of the `k0` pairs of time series (rows of `X`) that have the largest Pearson correlation. 

**INPUT**:\
`X`: Time series with rows indexing the degrees of freedom and the columns indexing the time steps.\
`α`: Correlation threshold above which we keep the pairs of time series.\
`k0`: Number of pairs of time series that are kept.

**OUTPUT**:\
`keep`: List of the pairs of time series that are kept.

---

### list\_all\_subjects
*./eeg-tools.jl*

- `list_all_subjects(n::Int64)`

Returns a `String` list of all subjects indices with three digits (from "001" up to "999"). 

**INPUT**:\
`n`: Number of subjects.

**OUTPUT**:\
`subjects`: List of subjects.

---

### load\_and\_save\_eeg\_avg
*./eeg-tools.jl*

- `load_and_save_eeg_avg(subjects::Vector{String}, states::Vector{String}, nz::Int64)`

Loads and saves the time series of `subjects` for tasks `states` in "eeg-data" as .csv files. Each time series is saved under the file name "SxxxRyy-X.csv".

**INPUT**:\
`subjects`: List of subjects whose time series will be saved.\
`states`: List of tasks' indices that will be saved.\
`nz`: Number of zones to gather time series.

---

### mySINDy
*./this.jl*

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
`err`: Squared difference between the inferred time series for `Y` and `Y` itself.

---

### mySINDy\_par
*./this-parallel.jl*

- `mySINDy\_par(arg::Tuple{Int64,Matrix{Float64},Matrix{Float64},Int64,Float64,Float64,Int64})`

Parallel version of mySINDy, to be used without filtering!

**INPUT** (in arg):\
`i`: Index of the time series analyzed.\
`X`: Time series of the agents' states.\
`Y`: Time series of the time derivative of agent `i`.\
`dmax`: Maximal degree to consider in the monomial library.\
`λ`: SINDy's threshold deciding whether an hyperedge exists or not.\
`ρ`: Regularization parameters promoting sparsity.\
`niter`: Maximal number of iterations for THIS algorithm.

**OUTPUT** (written in a file):\
`Ξ`: Matrix of coefficient inferred by SINDy. The index of each row is the index of an agent and the the columns corresponds to elements of the Taylor basis.\
`err`: 1-norm of the difference between `Y` and the inferred dynamics.

---

### mySINDy\_par\_filter
*./this-parallel.jl*

- `mySINDy\_par\_filter(arg::Tuple{Int64,Matrix{Float64},Matrix{Float64},Vector{Int64},Vector{Vector{Int64}},Dict{Vector{Int64},Int64},Float64,Float64,Int64})`

Parallel version of mySINDy, to be used with filtering!

**INPUT** (in arg):\
`i`: Index of the time series analyzed.\
`X`: Time series of the agents' states.\
`Y`: Time series of the time derivative of agent `i`.\
`i2keepi`: Filtering parameter.\
`keep`: List of pairs of agents to be considered.\
`d2i`: Filtering parameter.\
`λ`: SINDy's threshold deciding whether an hyperedge exists or not.\
`ρ`: Regularization parameters promoting sparsity.\
`niter`: Maximal number of iterations for THIS algorithm.

**OUTPUT** (written in a file):\
`Ξ`: Matrix of coefficient inferred by SINDy. The index of each row is the index of an agent and the the columns corresponds to elements of the Taylor basis.\
`ids`: List of monomials indices that are present in `Ξ`.\
`err`: 1-norm of the difference between `Y` and the inferred dynamics.

---

### read\_eeg
*./eeg-tools.jl*

- `read_eeg(file::String)`

Reads the .edf file `file` and returns its content as dictionary. 

**INPUT**:\
`file`: Name of the file to be loaded.

**OUTPUT**:\
`s2signal`: Dictionary associating the time series to each sensor name.

---

### restrict\_box\_size
*./eeg-tools.jl*

- `restrict_box_size(X::Matrix{Float64}, nstep::Int64)`

Restrict the "box" size (see manuscript), by keeping only the `nstep` time steps that are the closest to the median value of the time series. 

**INPUT**:\
`X`: Time series. Row indices are the agents time series and columns are the time steps.\
`nstep`: Number of time steps to be kept.

**OUTPUT**:\
`Y`: Relevant time steps.\
`ids`: Indices of the kept time steps.

---

### one2dim
*./this-tool.jl*

- `one2dim(Ainf::Dict{Int64,Matrix{Float64}}, d::Int64=1)`

Transforms the hypergraph inferred for `n``d` agents with a single internal dimension to a hypergraph of `n` agents with `d` internal dimension.

**INPUT**:\
`Ainf`: Dictionary of the inferred adjacency lists (output of `hyper_inf`).\
`d`: Internal dimension of agents.

**OUTPUT**:\
`AAinf`: Dictionary of the inferred adjacency lists with aggregted node dimensions.

---

### this
*./this.jl*

- `this(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)`

Infers the hypergraph underlying the dynamics of its vertices with knowledge of the states 'X' and of the derivatives 'Y' at each vertex.

**INPUT**:\
`X`: Time series of the system's state. Each row is the time series of the state of one agent.\
`Y`: Time series of the system's velocity. Each row is the time series of the velocity of on agent.\
`ooi`: Orders of interest. Vector of integers listing the orders of interactions that we analyze.\
`dmax`: Maximal degree to be considered in the Taylor expansion. Typically, dmax=maximum(ooi).\
`λ`: SINDy's threshold deciding whether an hyperedge exists or not.\
`ρ`: Regularization parameters promoting sparsity.\
`niter`: Maximal number of iterations for THIS algorithm.

**OUTPUT**:\
`Ainf`: Dictionary associating the inferred coefficient of the Taylor expansion to a pair (node,hyperedge). Namely, 'Ainf[(i,h)]' is the coefficient corresponding to the hyperedge 'h' in the Taylor expansion of the dynamics of node 'i'. If agents have internal dimensions larger than 1, `Ainf` is just a boolean.\
`coeff`: Matrix of coefficents obtained by SINDy. The row indices are the agents' indices and the columns indices are the indices of the monomials in the Taylor series.\
`relerr`: Relative error, i.e., `err` normalized by the magnitude of `Y`.

---

### this\_filter
*./this.jl*

- `this\_filter(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, par_keep::Union{Int64,Float64}, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)`

Infers the hypergraph underlying the dynamics of its vertices with knowledge of the states 'X' and of the derivatives 'Y' at each vertex.

**INPUT**:\
`X`: Time series of the system's state. Each row is the time series of the state of one agent.\
`Y`: Time series of the system's velocity. Each row is the time series of the velocity of on agent.\
`ooi`: Orders of interest. Vector of integers listing the orders of interactions that we analyze.\
`dmax`: Maximal degree to be considered in the Taylor expansion. Typically, dmax=maximum(ooi).\
`par_keep`: Parameter for the pairs to keep. If `Float64`, correlation threshold above which the pairs of time series should be kept. If `Int64`, number of pairs of time series to be kept by filtering (based on correlation).\
`λ`: SINDy's threshold deciding whether an hyperedge exists or not.\
`ρ`: Regularization parameters promoting sparsity.\
`niter`: Maximal number of iterations for THIS algorithm.

**OUTPUT**:\
`Ainf`: Dictionary associating the inferred coefficient of the Taylor expansion to a pair (node,hyperedge). Namely, 'Ainf[(i,h)]' is the coefficient corresponding to the hyperedge 'h' in the Taylor expansion of the dynamics of node 'i'. If agents have internal dimensions larger than 1, `Ainf` is just a boolean.\
`coeff`: Matrix of coefficents obtained by SINDy. The row indices are the agents' indices and the columns indices are the indices of the monomials in the Taylor series.\
`relerr`: Relative error, i.e., `err` normalized by the squared magnitude of `Y`.

---

### this\_par
*./this-parallel.jl*

- `this\_par(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)`

Parallel version of THIS. Infers the hypergraph underlying the dynamics of its vertices with knowledge of the states `X` and of the derivatives `Y` at each vertex.

Requires to store some data in a directory named `data`.

**INPUT**:\
`X`: Time series of the system's state. Each row is the time series of the state of one agent.\
`Y`: Time series of the system's velocity. Each row is the time series of the velocity of on agent.\
`ooi`: Orders of interest. Vector of integers listing the orders of interactions that we analyze.\
`dmax`: Maximal degree to be considered in the Taylor expansion. Typically, dmax=maximum(ooi).\
`λ`: SINDy's threshold deciding whether an hyperedge exists or not.\
`ρ`: Regularization parameters promoting sparsity.\
`niter`: Maximal number of iterations for THIS algorithm.

**OUTPUT**:\
`Ainf`: Dictionary associating the inferred coefficient of the Taylor expansion to a pair (node,hyperedge). Namely, 'Ainf[(i,h)]' is the coefficient corresponding to the hyperedge 'h' in the Taylor expansion of the dynamics of node 'i'. If agents have internal dimensions larger than 1, `Ainf` is just a boolean.\
`coeff`: Matrix of coefficents obtained by SINDy. The row indices are the agents' indices and the columns indices are the indices of the monomials in the Taylor series.\
`relerr`: Relative error, i.e., `err` normalized by the magnitude of `Y`.

---

### this\_par\_filter
*./this-parallel.jl*

- `this\_par\_filter(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, par_keep::Union{Int64,Float64}=.9, λ::Float64=.1, ρ::Float64=1., niter::Int64=10)`

Parallel verion of `this_filter`. Infers the hypergraph underlying the dynamics of its vertices with knowledge of the states `X` and of the derivatives `Y` at each vertex.

Requires to store some data in a directory named `data`.

**INPUT**:\
`X`: Time series of the system's state. Each row is the time series of the state of one agent.\
`Y`: Time series of the system's velocity. Each row is the time series of the velocity of on agent.\
`ooi`: Orders of interest. Vector of integers listing the orders of interactions that we analyze.\
`dmax`: Maximal degree to be considered in the Taylor expansion. Typically, dmax=maximum(ooi).\
`par_keep`: Parameter for the pairs to keep. If `Float64`, correlation threshold above which the pairs of time series should be kept. If `Int64`, number of pairs of time series to be kept by filtering (based on correlation).\
`λ`: SINDy's threshold deciding whether an hyperedge exists or not.\
`ρ`: Regularization parameters promoting sparsity.\
`niter`: Maximal number of iterations for THIS algorithm.

**OUTPUT**:\
`Ainf`: Dictionary associating the inferred coefficient of the Taylor expansion to a pair (node,hyperedge). Namely, 'Ainf[(i,h)]' is the coefficient corresponding to the hyperedge 'h' in the Taylor expansion of the dynamics of node 'i'. If agents have internal dimensions larger than 1, `Ainf` is just a boolean.\
`coeff`: Matrix of coefficents obtained by SINDy. The row indices are the agents' indices and the columns indices are the indices of the monomials in the Taylor series.\
`relerr`: Relative error, i.e., `err` normalized by the squared magnitude of `Y`.

---



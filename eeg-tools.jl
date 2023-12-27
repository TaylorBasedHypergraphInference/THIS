using EDF, DelimitedFiles, FFTW, PyPlot, Statistics

include("hyper-inf.jl")

# ================================================================================
"""
	average_over_zones(s2signal::Dict{String,Vector{Float64}}, s2z::Dict{String,Int64})

Averages the time series of all the sensors within each zone. Average is taken at each time step independently. 

_INPUT_:\\
`s2signal`: Contains the time series (output of `read_eeg`).\\
`s2z`: Indicates how sensors are gathered into zones. Each sensor is associated to a zone.

_OUTPUT_:\\
`asig`: Matrix of the time series. Row indices are zone indices and columns indices are time steps.
"""
function average_over_zones(s2signal::Dict{String,Vector{Float64}}, s2z::Dict{String,Int64})
	l = maximum(values(s2z))
	T = length(s2signal[collect(keys(s2signal))[1]])

	c = zeros(l)
	asig = zeros(l,T)
	for k in keys(s2signal)
		z = s2z[k]
		c[z] += 1
		asig[z,:] = asig[z,:]*(c[z]-1)/c[z] + s2signal[k]/c[z]
	end

	return asig
end

# ================================================================================
"""
	denoise_fourier(x::Vector{Float64}, nmodes::Int64)\\
	denoise_fourier(X::Matrix{Float64}, nmodes::Int64)

Gets rid of the noise in time series by low-pass filtering. Only the `nmodes` lowest modes of the FFT are kept.

_INPUT_:\\
`x`, `X`: Time series. Either a single one or a matrix thereof.\\
`nmodes`: Number of modes of the FFT to be kept.

_OUTPUT_:\\
`y`, `Y`: Denoised time series. 
"""
function denoise_fourier(x::Vector{Float64}, nmodes::Int64)
	T = length(x)
	f = fft(x)
	g = zeros(Complex{Float64},T)
	g[1:nmodes+1] = f[1:nmodes+1]
	g[T-nmodes+1:T] = f[T-nmodes+1:T]
	return real.(ifft(g))
end

function denoise_fourier(X::Matrix{Float64},nmodes::Int64)
	n,T = size(X)
	Y = zeros(n,T)
	for i in 1:n
		Y[i,:] = denoise_fourier(X[i,:],nmodes)
	end
	return Y
end

# ================================================================================
"""
	list_all_subjects(n::Int64)

Returns a `String` list of all subjects indices with three digits (from "001" up to "999"). 

_INPUT_:\\
`n`: Number of subjects.

_OUTPUT_:\\
`subjects`: List of subjects.
"""
function list_all_subjects(n::Int64)
	subjects = String[]
	for i in 1:min(9,n)
		push!(subjects,"00"*string(i))
	end
	for i in 10:min(99,n)
		push!(subjects,"0"*string(i))
	end
	for i in 100:min(999,n)
		push!(subjects,string(i))
	end

	return subjects
end

# ================================================================================
"""
	load_and_save_eeg_avg(subjects::Vector{String}, states::Vector{String}, nz::Int64)

Loads and saves the time series of `subjects` for tasks `states` in "eeg-data" as .csv files. Each time series is saved under the file name "SxxxRyy-X.csv".

_INPUT_:\\
`subjects`: List of subjects whose time series will be saved.\\
`states`: List of tasks' indices that will be saved.\\
`nz`: Number of zones to gather time series.
"""
function load_and_save_eeg_avg(subjects::Vector{String}, states::Vector{String}, nz::Int64)
	s = readdlm("eeg-data/sensors-$nz.csv",',',String)
	z = readdlm("eeg-data/zones-$nz.csv",',',Int64)
	s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))
	for su in subjects
		for st in states
			@info "Working on S"*su*"R"*st
			s2signal = read_eeg("eeg-data/S"*su*"R"*st*".edf")
			asig = average_over_zones(s2signal,s2z)
			writedlm("eeg-data/S"*su*"R"*st*"-X.csv",asig,',')
		end
	end

	return nothing
end

# ================================================================================
"""
	read_eeg(file::String)

Reads the .edf file `file` and returns its content as dictionary. 

_INPUT_:\\
`file`: Name of the file to be loaded.

_OUTPUT_:\\
`s2signal`: Dictionary associating the time series to each sensor name.
"""
function read_eeg(file::String)
	xxx = EDF.read(file)

	s2signal = Dict{String,Vector{Float64}}()

	for sig in xxx.signals
		if typeof(sig) == EDF.Signal{Int16}
			s2signal[sig.header.label] = EDF.decode(sig)
		end
	end

	return s2signal
end

# ================================================================================
"""
	restrict_box_size(X::Matrix{Float64}, nstep::Int64)

Restrict the "box" size (see manuscript), by keeping only the `nstep` time steps that are the closest to the median value of the time series. 

_INPUT_:\\
`X`: Time series. Row indices are the agents time series and columns are the time steps.\\
`nstep`: Number of time steps to be kept.

_OUTPUT_:\\
`Y`: Relevant time steps.\\
`ids`: Indices of the kept time steps.
"""
function restrict_box_size(X::Matrix{Float64}, nstep::Int64)
	n,T = size(X)
	ns = min(T-1,nstep)

	mX = median(X,dims=2)
	dX = X - repeat(mX,1,T)
	nX = vec(sum(dX.^2,dims=1))
	sX = sort(nX)
	th = (sX[ns+1] + sX[ns])/2
	
	ids = (1:T)[nX .< th]

	return X[:,ids], ids
end


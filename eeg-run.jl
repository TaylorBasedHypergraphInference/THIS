include("eeg-tools.jl")

# Data are taken from the PhysioNet data set: https://physionet.org/content/eegmmidb/1.0.0/
# It is assumed that the .edf files with EEG data are stored in "./eeg-data/", under the names "SxxxRyy.edf". 

nz = 7 # Number of zones on the scalp.

# Parameters for SINDy
λ = .1
ρ = .1

# Data to be loaded. We recommend to load a subset of the subjects or of the states.
subjects = list_all_subjects(109) # Runs THIS on all the time series (takes ~60 seconds on a standard laptop).
#subjects = ["001","002"]
states = ["01","02"] # Resting state
#states = ["03","07","11"] # Task 1

# Pairing between sensors to zones
s = readdlm("eeg-data/sensors-$nz.csv",',',String)
z = readdlm("eeg-data/zones-$nz.csv",',',Int64)
s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))

# Storage of inferences outputs
AA2 = Dict{String,Matrix{Float64}}()
AA3 = Dict{String,Matrix{Float64}}()
AA4 = Dict{String,Matrix{Float64}}()
relerr = Dict{String,Float64}()

for su in subjects
	for st in states
		@info "Running S"*su*"R"*st

		# Loading data
		file = "eeg-data/S"*su*"R"*st*".edf"
		s2signal = read_eeg(file)
		
		# Average signal over zones
		asig = average_over_zones(s2signal,s2z)
		
		# Finite differences
		dt = 1/160
		truncat = findmin(vec(maximum(abs.([asig zeros(nz)]),dims=1)) .> 1e-6)[2] - 1
		
		# We get rid of the noise by low-pass filtering the time series, keeping only the 100 lowest components of the FFT.
		X0 = denoise_fourier(asig[:,1:truncat],100)
		Y0 = (X0[:,2:end]-X0[:,1:end-1])./dt
		X0 = X0[:,1:end-1]

		# We normalize all time series so all orders are of comparable magnitudes.
		X = X0./mean(abs.(X0))
		Y = Y0./mean(abs.(Y0))

		# We restrict the inference to a limited "box" size (see manuscript). 
		X,ids = restrict_box_size(X,1000)
		Y = Y[:,ids]

		# Inference
		ooi = [2,3]
		dmax = 4
		
		Ainf,coeff,re = this(X,Y,ooi,dmax,λ,ρ)
		
		if 2 in keys(Ainf)
			global AA2["S"*su*"R"*st] = Ainf[2]
		end
		if 3 in keys(Ainf)
			global AA3["S"*su*"R"*st] = Ainf[3]
		end
		if 4 in keys(Ainf)
			global AA4["S"*su*"R"*st] = Ainf[4]
		end
		global relerr["S"*su*"R"*st] = re
	end
end



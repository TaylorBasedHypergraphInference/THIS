include("eeg-tools.jl")

# Data are taken from the PhysioNet data set: https://physionet.org/content/eegmmidb/1.0.0/
# It is assumed that the .edf files with EEG data are stored in "./eeg-data/", under the names "SxxxRyy.edf". 

nz = 7 # Number of zones on the scalp.

# Parameters for SINDy
λ = .1
ρ = .1

# Data to be loaded. We recommend to load a subset of the subjects or of the states.
#subjects = list_all_subjects(109) # Runs THIS on all the time series (can be long).
subjects = ["001","002"]
states = ["01","02"] # Resting state
#states = ["03","07","11"] # Task 1

# Pairing between sensors to zones
s = readdlm("eeg-data/sensors-$nz.csv",',',String)
z = readdlm("eeg-data/zones-$nz.csv",',',Int64)
s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))

# Storage of inferences outputs
AA2 = zeros(Int64,nz,nz)
AA3 = zeros(Int64,nz,nz,nz)
relerr = zeros(length(states),0)

for subject in subjects
	re = Float64[]
	for state in states
		@info "Running S"*subject*"R"*state

		# Loading data
		file = "eeg-data/S"*subject*"R"*state*".edf"
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
		xxx = hyper_inf(X,Y,ooi,dmax,λ,ρ)
		push!(re,xxx[4])
		
		# Retrieve adjacency tensors	
		A2 = inferred_adj_2nd(xxx[1][2],nz)[2]
		A3 = inferred_adj_3rd(xxx[1][3],nz)[2]
		B2 = (abs.(A2) .> 1e-8)
		B3 = (abs.(A3) .> 1e-8)

		# Collect all boolean adjacency tensors
		global AA2 += B2
		global AA3 += B3

		writedlm("eeg-data/S"*subject*"R"*state*"-A2.csv",A2,',')
		x = zeros(nz,0)
		for i in 1:nz
			x = [x A3[:,:,i]]
		end
		writedlm("eeg-data/S"*subject*"R"*state*"-A3.csv",x,',')

		writedlm("eeg-data/S"*subject*"R"*state*"-coeff.csv",xxx[2],',')
		end
	global relerr = [relerr re]
end

writedlm("eeg-data/relative-error.csv",relerr,',')


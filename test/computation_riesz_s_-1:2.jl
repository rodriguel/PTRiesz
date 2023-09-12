include("../src/PTRiesz.jl")

PATH_OUTPUT = "test/datas/computation_riesz_s_-0.5"


# Number of chains 
n_chains = Threads.nthreads()
# Number of particles 
N = 10
# Number of bins
nbins = div(N, 2) * 10
# Value of s (here, s = 0)
s = 0.
# Temperatures 
Ts = [.1, .2, .5, 1., 5., 10]
# Inverse temperatures 
βs = reverse(1 ./ Ts)
# Creation of the Hamiltonian
logp = build_riesz_hamiltonian(s, N)

# List of histograms for each temperatures
pchs = [RWM.PCHist(N, s, β, nbins) for β in βs]

# Number of samples in Metropolis-Hastings
n_samples = 10_000

# We populate the histogram by Metropolis-Hastings 
for j in eachindex(βs)
    β = βs[j]
    println("Computing for beta = $(β) and s = $s")
    pchs[j] = RWM.compute_pair_correlation(N, s, logp, β, n_samples, n_chains; nbins)
end

# Save the histograms 
RWM.save_histogram(pchs, PATH_OUTPUT, "histograms")
include("../src/PTRiesz.jl")
include("plotting.jl")

PATH_OUTPUT = "git/test/datas/computation_several_N"


# Number of chains 
n_chains = Threads.nthreads()
# Number of particles 
Ns = [40, 60, 80, 100, 120, 150]
# Value of s (here, s = 0)
s = -0.5
# Temperatures 
β = 1/1.2

# List of histograms for each temperatures
pchs = [RWM.PCHist(N, s, β, div(N, 2) * 10) for N in Ns]

# Number of samples in Metropolis-Hastings
n_samples = 10_000_000

# We populate the histogram by Metropolis-Hastings 
for j in eachindex(Ns)
    N = Ns[j]
    nbins = div(N, 2) * 10
    logp = build_riesz_hamiltonian(s, N)
    println("Computing for beta = $(β) and s = $s and N = $N")
    pchs[j] = RWM.compute_pair_correlation(N, s, logp, β, n_samples, n_chains; nbins)
end

# Save the histograms 
RWM.save_histogram(pchs, PATH_OUTPUT, "histograms")

# Reading 
p = my_plot()
for j in eachindex(Ns)
    N = Ns[j]
    rs = create_rs(N, div(N, 2) * 10)
    plot!(p, rs, pchs[j].g)
end 
display(p)
module RWM

    include("utility.jl")
    include("potential.jl")
    include("display.jl")
    
    using Random
    using JLD2          # Saving to .jld2 format
    using ProgressBars  # Show progress bar in verbose mode
    export RWMState
    export PCHist

    """
        RWMState

        Defines a structure to hold the chain. It contains:
            - N : the number of particles (i.e. dimension of the space)
            - x : the current state
            - y : the proposal state
            - u : the proposed jump, i.e. y = x + u
            - q! : the proposal distribution (i.e. centered Gaussian in R^N)
            - σ : the 'size' of the proposal jump (i.e. the variance of the Gaussian)
            - rnd : a random number generator (e.g. MersenneTwister)
    """
struct RWMState
    N::Int;
    x::Vector{Float64};
    y::Vector{Float64}; 
    u::Vector{Float64};
    q!::Function;
    σ::Float64;
    rnd::Random.AbstractRNG;
end

"""
    RWMState(x; σ = 1., rnd = MersenneTwister(), q! = randn!)

    Constructs a structure RWM to hold the chain.
"""
function RWMState(x::Vector{Float64}; σ::Float64=1., rnd::Random.AbstractRNG=MersenneTwister(), 
    q!::Function=randn!)
    N = length(x)
    y = copy(x) 
    u = zeros(N)
    return RWMState(N, copy(x), y, u, q!, σ, rnd)
end


"""
    The inlined functions (draw_u!, update_y!, accept!) are self-explainable.
"""
@inline function draw_u!(s::RWMState)
    s.q!(s.rnd, s.u)
    nothing 
end 

@inline function update_y!(s::RWMState)
    @inbounds s.y .= s.x .+ s.σ*s.u
    wrap!(s.y)
    nothing 
end 

@inline function accept!(s::RWMState)
    @inbounds s.x .= s.y
    nothing
end 

"""
    compute_α(energy_y, energy_x, β)

    Computes the acceptance ratio from the proposal 'energy_y'
    and the current 'energy_x' at temperature 1/β.
"""
@inline compute_α(energy_y, energy_x, β) = exp(β*(energy_y - energy_x))


@inline accepted(s, α) = rand(s.rnd) < α
"""
    sample!(s, log_p, β)

    Samples a single (!) state from the chain 's' with target distribution
    'log_p' (i.e. the total energy) at temperature '1/β' (default β = 1).

"""
function sample!(s::RWMState,log_p::Function, β::Float64 = 1)
    energy_x = - log_p(s.x); energy_y = zero(energy_x)
    draw_u!(s); update_y!(s)
    energy_y = -log_p(s.y)
    α = compute_α(energy_y, energy_x, β)
    acptd = false
    if accepted(s, α)
        acptd = true
        accept!(s)
        energy_x = energy_y 
    end
    return acptd, energy_x
end
"""
    sample!(s, log_p, nsamples, β)

    Samples 'nsamples' state from the chain 's' with target distribution
    'log_p' (i.e. the total energy) at temperature '1/β' (default β = 1).
    This is just a recursive call to the previous function. 

"""
function sample!(s::RWMState,log_p::Function, nsamples::Int64, β::Float64 = 1)
    for _ = 1:nsamples
        sample!(s, log_p, β)
    end
end 


struct PCHist
    cur_hist::Vector{Int64};
    hist::Vector{Int64};
    g :: Vector{Float64};
    nbins::Int64;
    N::Int64;
    s::Float64;
    β::Float64;
end

function PCHist(N, s, β, nbins)
    cur_hist = zeros(nbins)
    hist = zeros(nbins)
    g = zeros(nbins)
    return PCHist(cur_hist, hist, g, nbins, N, s, β)
end 

function save_histogram(pch, folder, file_name)
    path = "$(folder)/$(file_name).jld2"
    save_object(path, pch)
end


function load_histogram(path)
    pch = load_object("$(path).jld2")
    return pch
end

""" 
    normalize_histogram(hist, N)

    Returns the 'normalized' histogram (i.e. the pair correlation),
    that this which tends to 1 at infinity and/or which has integral (N-1)/2.
    The version 'normalize_historgram(pch :: PCHist)' can directly be applied
    to the structure.

"""
function normalize_histogram(hist, N)
    nbins = length(hist)
    return hist/sum(hist) * (N-1)/N * nbins
end

function normalize_histogram(pch::PCHist)
    pch.g .= normalize_histogram(pch.hist, pch.N)
    nothing
end 

""" 
compute_fourier_histogram(hist, N)

    Returns the Fourier transform of the histogram (i.e. the structure factor),
    or to be more precise of the truncated pair correlation ('histogram - 1').
    The version 'compute_fourier_histogram(pch :: PCHist)' can directly be applied
    to the structure.

"""

function compute_fourier_histogram(hist, N)
    nbins = length(hist)
    xs = [N/(2*nbins) * (i-1/2) for i=1:nbins]
    S(q) = 1 + 2 * N/(2*nbins) * sum(cos.(2*π*q*xs) .* (hist .- 1))
    return S    
end 

compute_fourier_histogram(pch::PCHist) = compute_fourier_histogram(pch.g, pch.N)

"""
    update_hist!(pch)

    Updates the histogram by adding a new sample. 

"""
function update_hist!(pch::PCHist)
    @inbounds pch.hist .+= pch.cur_hist  
    nothing
end


"""
    compute_hist!(h, nbins, x, lx)

    Compute the histogram from a configuration of N points.
    Since one needs to iterate over all the pair of particles
    in the configuration, the cost is of order O(N^2). 
    This is a bottleneck function.

"""
function compute_hist!(h, nbins, x, lx)
    for i=1:(lx-1)
        for j=(i+1):lx
            d = distance_circle(lx, x[i], x[j])
            pos = ceil(Int64, 2 * nbins * d/lx)
            pos = (pos == 0) ? 1 : pos
            h[pos] += 1
        end
    end
    nothing
end


"""
    sample_hist!(s, pch, log_p, nsamples, β)

    This is the main function. It simulates the chain 's' with 
    target distribution 'log_p' (i.e. the total energy) at temperature 1/β,
    during 'nsamples' times. It updates along the histogram (e.g. pair correlation)
    contained in the structure 'pch'.

"""
function sample_hist!(s::RWMState, pch::PCHist, log_p::Function, nsamples, β = 1.; verbose = true)
    N = s.N
    cnt = 0
    pb = (verbose) ? ProgressBar(1:nsamples) : 1:nsamples
    for _ in pb 
        accepted, _ = sample!(s, log_p, β)
        cnt = accepted ? cnt + 1 : cnt
        compute_hist!(pch.cur_hist, pch.nbins, s.x, N)
        update_hist!(pch)
    end
    return cnt/nsamples
end 

function sample_energy!(s::RWMState, log_p::Function, nsamples, β = 1.; verbose = true)
    N = s.N
    cnt = 0
    energies = zeros(nsamples)
    pb = (verbose) ? ProgressBar(1:nsamples) : 1:nsamples
    for k in pb
        accepted, e = sample_single!(s, log_p, β)
        cnt = accepted ? cnt + 1 : cnt
        energies[k] = e
    end
    return energies, cnt/nsamples
end

function sample_energy_parallel!(chains::Vector{RWMState}, log_p, nsamples_by_chain, β = 1.)
    energies = zeros(nsamples_by_chain, length(chains))
    Threads.@threads for c in eachindex(chains)
        e = @view energies[:, c]
        e .= sample_energy!(chains[c], log_p, nsamples_by_chain, β; verbose = false)[1]
    end
    return energies
end 

"""
    sample_hist_parallel!(chains, pchs, log_p, nsamples_by_chain , β)

    Parallel (i.e. multithreading) version of the 'sample_hist!', as 
    defined above. 

"""
function sample_hist_parallel!(chains::Vector{RWMState}, pchs::Vector{PCHist}, log_p, nsamples_by_chain, β = 1.)
    acrs = zeros(length(chains)) 
    Threads.@threads for c in eachindex(chains)
        acrs[c] = sample_hist!(chains[c], pchs[c], log_p, nsamples_by_chain, β; verbose = false)
    end
    return sum(acrs)/length(chains)
end


function concatenate_histograms(pchs::Vector{PCHist})
    N, nbins, npch = pchs[1].N, pchs[1].nbins, length(pchs)
    hist = zeros(nbins)
    g = zeros(nbins)
    for i in eachindex(pchs)
        hist .+= pchs[i].hist
        g .+= normalize_histogram(pchs[i].hist, N)
    end
    return hist, g./npch
end 

function concatenate_histograms(pch::PCHist, pchs::Vector{PCHist})
    hist, g = concatenate_histograms(pchs)
    pch.hist .= hist; pch.g .= g 
    nothing
end

function creating_chains_histograms(nchains, N, s, β, nbins; crystal = true)
    x = (crystal) ? init_crystal(N) : init_random(N); σ = best_sigma(β, N)
    chains = [RWM.RWMState(x; σ) for c = 1:nchains]
    pchs = [RWM.PCHist(N, s, β, nbins) for c = 1:nchains]
    return chains, pchs
end 


function compute_pair_correlation(N, s, logp, β, nsamples, nchains; nbins = div(N,2) * 10)
    nsamples_by_chain = div(nsamples, nchains)
    pch = RWM.PCHist(N, s, β, nbins)
    chains, pchs = creating_chains_histograms(nchains, N, s, β, nbins)
    acc = RWM.sample_hist_parallel!(chains, pchs, logp, nsamples_by_chain, β)
    RWM.concatenate_histograms(pch, pchs)
    return pch
end

function compute_pair_correlation(N, s, β, nsamples, nchains)
    log_p = build_riesz_hamiltonian(s, N)
    return compute_pair_correlation(N, s, log_p, β, nsamples, nchains)
end







#########################
# Some utility function #


"""
    Bunch of small functions not needed. Kept 'in case'.
"""
@inline function update_mean(μ::Float64, x, k)
    μ = (k-1)/k*μ + 1/k*x
    return μ
end 

@inline function update_mean!(μ::Vector{Float64}, x, k)
    @inbounds μ .= (k-1)/k*μ + 1/k*x
    nothing
end

@inline function update_var(σ::Float64, old_μ, μ, k)
    σ = (k-1)/k*σ + (k+1)*(μ .- old_μ)^2
    return σ
end

@inline function update_var(σ::Vector{Float64}, old_μ, μ, k)
    @inbounds σ .= (k-1)/k*σ + (k+1)*(μ - old_μ).^2
    nothing
end 

end
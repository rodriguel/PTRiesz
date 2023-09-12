using LinearRegression
#using Optim
#using CurveFit
#using LinearAlgebra
"""
    mod_circle(N, r)

    Computes the projection of the point r onto the circle of length N.
    Here, the fundamental domain is taken to be the segment [-N/2, N/2].
"""
@inline @fastmath mod_circle(N, r) = @. mod(r + N/2, N) - N/2

"""
    wrap!(x)

    Projects the vector x (e.g. inplace) onto its image into circle of length N.
"""
function wrap!(x)
    N = length(x)
    for i = 1:N
        x[i] = mod_circle(N, x[i])
    end
end

"""
    distance_circle(N, r, s)

    Compute the distance of two points on the circle of size N. Beware that
    the points need be projected onto the circle before using this function !
"""
@inline @fastmath distance_circle(N, r, s) = min(abs(r-s), N - abs(r-s))

"""
    init_random(N)

    Returns a random configuration of N points onto the circle of length.
    This function is used for troubleshooting. Indeed, from the perspective
    of MCMC, it is better to start directly onto the Wigner crystal 
    (See 'init_crystal' function).  

"""
init_random(N) = N*rand(N) .- N/2 

"""
    init_crystal(N)

    Returns the Wigner crystal configuration onto the circle of length N,
    which is known to be the minimizer of the energy (hence, the mode of
    the canonical ensemble at any temperature.)

"""
init_crystal(N) = Array(LinRange(-N/2 + 1/2, N/2 - 1/2, N))


"""
        best_sigma(β, N)

        Returns a 'good' proposal jump for the MCMC, that is
        to yield a good acceptance rate. Although this 
        has been determined phenomenologically, the expression is
        rather logical (i.e. Gaussian in R^N with variance sqrt(N*T)...)
"""
best_sigma(β, N) = 0.4*sqrt(2/β)*sqrt(1/N)

"""
    regress_S(S, N; n_min = N/8, n_max = N/4)

    Regresses the structure factor of the interval
    n_min/N to n_max/N with step size 1/N.
"""


"""
function regress_S_(S, N; n_min, n_max)
    ks = n_min:1/N:n_max
    ys = S.(ks)
    lr = linregress(log.(ks), log.(ys))
    a, b = LinearRegression.coef(lr)
    return a, exp(b)
end 

"""

function regress_S_2(S, N; n_min, n_max)
    ks = n_min/N:1/N:n_max/N
    ys = S.(ks)
    b, a = power_fit(ks, ys)
    return b, a
end 



function regress_S(S, N, β; n_min, n_max)
    ks = n_min/N:1/N:n_max/N
    n = length(ks); ys = S.(ks)
    w = [1/k^2 for k=1:n]
    function L(θ) 
        local ℓ(k) = θ[1] * abs(k)^θ[2]
        tmp = @. w * abs(ℓ(ks) - ys)
        return sum(tmp)/n
    end
    res = Optim.optimize(L, ones(Float64, 2))
    C, a = Optim.minimizer(res)
    min = Optim.minimum(res)
    return C, a, min
end

function regress_S_(S, N, β, lmax; n_min)
    etas = zeros(length(lmax))
    for i in eachindex(lmax)
        n_max = lmax[i]
        _, etas[i], _ = regress_S(S, N, β; n_min, n_max)
    end
    return etas
end 
"""
function regress_S_(S, N, lmax; n_min)
    etas = zeros(length(lmax))
    for i in eachindex(lmax)
        n_max = lmax[i]
        etas[i], _ = regress_S_(S, N; n_min, n_max)
    end
    return etas
end 
"""

function regress_S_bkt(S, N, n_min, n_max)
    function L(θ) 
        local ℓ(k) = θ[1] * (abs(k - 1)^θ[2])
        ks = n_min/N:1/N:n_max/N
        ys = S.(ks)
        tmp = @. abs2(ℓ(ks) - ys)
        return sum(tmp)/length(ks)
    end
    res = Optim.optimize(L, rand(2))
    C, a = Optim.minimizer(res)
    return C, a
end



function regress_S_peak(S, N; x_min)
    ks = 1/N:1/N:(1-x_min)
    s(x) = S(1 - x); ys = s.(ks)
    lr = linregress(log.(ks), log.(ys))
    a, b = LinearRegression.coef(lr)
    return a, exp(b)
end

function slope_S(S, N; n_max = 5)
    ks = 1/N:1/N:n_max/N
    ys = S.(ks)
    lr = linregress(ks, ys)
    a, b = LinearRegression.coef(lr)
    return a
end 
"""
    regress_S(SS, N; n_min, n_max)

    Same as previous function, but takes an array 
    of histogram and returns an array of function.
    Not used, be kept 'in case'.

"""
function regress_S(Ss::Array{Function, 1}, N;
    n_min, n_max)
    coefs = Array{Array{Float64}, 1}(undef, length(Ss))
    ks = n_min:1/N:n_max
    for s in eachindex(Ss)
        coefs[s] = regress_S(S, N; n_min, n_max)
    end
    return coefs
end 
"""
function regress_peak(S, N; n_min, x_max = 1 - 1/N)
    ks = n_min/N:1/N:x_max; 
    n = length(ks)
    ys = S.(ks)
    function L(θ) 
        local ℓ(k) = 1/abs(k - 1)^θ[1] + θ[2]
        tmp = @. abs2(ℓ(ks) - ys)
        return sum(tmp)/n
    end
    res = Optim.optimize(L, zeros(Float64, 2))
    return Optim.minimizer(res)
end
"""

"""
    Self-explanory functions. 
    Not used, kept 'in case'.
"""
@inline positive_part(x) =  x[x .>= 0]; 
@inline negative_part(x) = x[x .<= 0];

"""
    amplitude_of_g(g, N, n_period_average = 3)

    Return the mean amplitude of the oscillation (if any)
    of the pair correlation 'g' over 'n_period_average'
    period far from the origin.
"""

function amplitude_of_g(g, N, n_period_average = 3)
    nbins = length(g)
    nbins_per_period = 2 * div(nbins, N)
    amplitudes = zeros(n_period_average)
    for i = 1:n_period_average
        a = nbins - i*nbins_per_period
        b = nbins - (i - 1)*nbins_per_period
        h_ = g[a:b] .- 1
        amplitude_ = maximum(positive_part(h_), init = 0) - minimum(negative_part(h_), init = 0)
        amplitude = (amplitude_ < 0) ? 0 : amplitude_
        amplitudes[i] = amplitude
    end
    res = sum(amplitudes)/n_period_average
    return res
end

function create_rs(N, nbins)
    return [N/(2*nbins)*(i-1/2) for i=1:nbins]
end 



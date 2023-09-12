using SpecialFunctions

"""
    polylogarithm(u, z)

    Compute the polylogarithm function Li_u(z) using the Hurwitz zeta function,
    as implemented in the SpecialFunctions.jl package. 

function polylogarithm(u, z)
    a_u = gamma(1-u)/(2π)^(1-u)
    b_u = im^(1-u); c_u = im^(u-1)
    z_u1 = zeta(1-u, 1/2 + log(Complex(-z))/(2π*im))
    z_u2 = zeta(1-u, 1/2 - log(Complex(-z))/(2π*im))
    return real(a_u * (b_u * z_u1 + c_u * z_u2))
end

"""

function polylogarithm(u, z)
    a = 1/(1 + im^(2*u)) * (2*π*im)^u/gamma(u)
    z_u1 = zeta(1-u, 1/2 + log(Complex(-z))/(2π*im))
    z_u2 = zeta(1-u, 1/2 + log(Complex(-1/z))/(2π*im))
    return real(a * (z_u1 + z_u2))/2
end


"""
    periodic_riesz_(s, N, x)

    Compute the periodic s-Riesz potential for N particles using the polylogarithm.
    This function is rather time consuming, so that we instead use a interpolated version.
    (See 'build_periodic_riesz' function)
"""
function periodic_riesz_(s::Float64, N::Int64, x)
    #c_s = sign(s) * 2^(1.5 - s) * π * gamma((1-s)/2)/gamma(s/2)
    #c_s = N^(-s) * c_s * (2π)^(s-1-1/2)/abs(s)
    return  N^(-s) * real(polylogarithm(1-s, exp(-im*2*π*x/N)) + polylogarithm(1-s, exp(im*2*π*x/N)))/2
end


"""
    log_potential(s, N, x)

    Returns the periodic logarithmic potential (s=0). 
"""
function log_potential(N, x)
    return -log(abs(sin(π*x/N)))
end 

"""
    madelung(s, N)

    Compute the Madelung constant. In MCMC, we only care about the difference
    of energies, therefore we do need to incorporate the Madelung constant inside
    the Hamiltonian. We leave it here for troubleshooting matter.
"""
function madelung(s::Float64, N::Int)
    res = (s < 0) ? -2*zeta(s) : 2*zeta(s)
    return res/N^s 
end

"""
    riesz_potential(s, x)

    Compute the non-periodic Riesz potential for any (non-zero) s.
    We do not use this function, but we leave it here for troubleshooting matter.
"""
function riesz_potential(s::Float64, x)
    return sign(s)*1/abs(x)^s
end

"""
    build_periodic_riesz(s)

    Returns a function (i.e. closure), namely an (linear) interpolated version
    of the periodic Riesz potential of unit period — that is a tabulated version.
    The interpolation is simply taken to be a linear one. Beware that the grid size
    needs to be adjusted when the number of particles N is big for accuracy. Taking
    N^2 ensures that the error should goes to zero away from the boundaries.
"""

function build_periodic_riesz(s::Float64, N::Int64, gridSize = N*10_000)
    xx = LinRange(0, 1, gridSize)
    yy = [periodic_riesz_(s, 1, xx[step]) for step=1:gridSize]
    @fastmath function tab_periodic_riesz(x) 
        x = mod(abs(x), 1)
        x = (x < 1/gridSize) ? 1/gridSize : x # PQ ?
        i_after = ceil(Int64, gridSize * x)
        λ = gridSize*(xx[i_after] - x)
        res = (i_after > 1) ? λ*yy[i_after] + (1-λ)*yy[i_after - 1] : yy[i_after]
        res
    end
    return tab_periodic_riesz
end 

"""
    build_periodic_riesz(s, N)

    Returns a function (i.e. closure), namely the full Hamiltonian with the tabulated
    version the periodic Riesz potential for N particles. We do not appeal to some sort
    of FMM, so that the cost is of order O(N^2). This is the main bottleneck. 
"""
function build_riesz_hamiltonian(s::Float64, N::Int, gridSize = N*10_000)
    if s == 0 return x-> log_hamiltonian(N, x) end
    periodic_riesz = build_periodic_riesz(s, N, gridSize)
    @inbounds @fastmath function riesz_hamiltonian(x)
        res = 0;
        for i = 1:(N-1) 
            for j = (i+1):N
            res += periodic_riesz((x[i] - x[j])/N)
            end 
        end
        return 1/N^s * res
        end
    return riesz_hamiltonian
end

@inbounds @fastmath function log_hamiltonian(N::Int, x)
    res = 0;
    for i = 1:(N-1) 
        for j = (i+1):N
        res += log_potential(N, x[i] - x[j])
        end 
    end
    return res
end

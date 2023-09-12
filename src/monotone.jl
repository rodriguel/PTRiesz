
include("RWM.jl")
using DSP



"""
    is_increasing(y, threshold)

    Return a boolean to whether or not the structure
    factor is monotone (i.e. increasing) near k = 1,
    using the preceding functions. 

"""
function is_increasing(y, thr = 0.5)
    res = true
    for i in eachindex(y)
        if (y[i] < thr)
            res = false
            break
        end
    end
    return res
end 


"""
    smooth_signal(S, N, sw)
    
    Returns a smooth version of the structure factor
    by convoluting with a rectangle filter of size 'sw'.
    The convolution is made on [1/N, 1.5] (arbitrary)

"""
function smooth_signal(S, N, sw)
    ks = 1/N:1/N:1.5
    ys = S.(ks)
    rect_filter = rect(sw)/sw
    ys_smooth = conv(ys, rect_filter)
    return ys_smooth
end


"""
    deriv_signal(xs)
    
    Returns the discrete derivative of the signal 'xs'

"""
function deriv_signal(xs)
    d_filter = [1, -1]
    return conv(xs, d_filter)
end


function is_monotone(pch::RWM.PCHist; dn = 10, sw = 2, thr = 0.1, sz = 300)
    N = pch.N; 
    S = RWM.compute_fourier_histogram(pch.g, N)
    p = plot_S(S, N, pch.s, pch.β, :black, low=N-dn, up=N+dn)
    ys_smooth = smooth_signal(S, N, sw)[1:end-sw+1]
    y_deriv = deriv_signal(ys_smooth)[1:end-1]
    inc = is_increasing(y_deriv, thr)
    c = (inc) ? :blue : :red
    xx = 1/N:1/N:length(ys_smooth)/N
    plot!(p, xx, ys_smooth, c = c, ls = :dash, lw = 2,
    xlims = (1-(dn+1)/N, 1 + (dn + 1)/N), ylims = (0.8, 1.1),
    label = "Smoothed", size = (sz, sz))
    return p
end

function is_monotone(pchs::Vector{RWM.PCHist}; dn = 10, sw = 2, thr = 0.1, sz = 100)
    p = Any[]; s = pchs[1].s; npch = length(pchs)
    for j = 1:npch
        p_ = is_monotone(pchs[j]; dn, sw, thr, sz)
        plot!(p_, xlabel = "", ylabel = "", legend=false, title = "\$\\beta = $(pchs[j].β), s = $s\$")
        push!(p, p_)
    end
    p_final = plot(map(plot, p)..., size = (1000, 1000), dpi = 1000)
    return p_final
end



# Old version — kept in case

"""
    smooth_signal_near_1(S, N; number_neighbors = 10)
    
    Returns the smoothed structure factor near k = 1, using
    the 'number_neighbors' to the left and right, by appling
    a Gaussian convolution

"""
function smooth_signal_near_1(S, N; number_neighbors = 10)
    ks = 1-number_neighbors/N:1/N:1+number_neighbors/N
    l = length(ks);
    ys = S.(ks)
    gaussian_filter = gaussian(l, 2/N)
    ys_smooth = conv(ys, gaussian_filter)
    return ys_smooth
end


"""
    smooth_deriv_signal_near_1(S, N; number_neighbors = 10)

    Returns the discrete derivative of the structure factor
    near k = 1, which has been smoothed before (see above).
    The derivate filter is applied by considerind the 
    'number_neighbors' to the left and right. 

"""
function smooth_deriv_signal_near_1(S, N; number_neighbors = 10)
    ks = 1-number_neighbors/N:1/N:1+number_neighbors/N
    l = length(ks);
    ys = S.(ks)
    ys_smooth = smooth_signal_near_1(S, N; number_neighbors)
    ys_deriv =  conv(ys_smooth, N*[1,0,-1])
    return ys_deriv
end 
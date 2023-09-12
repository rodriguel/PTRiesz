using Plots
using LaTeXStrings
using ColorSchemes
using Colors

red = colorant"#ff0000"
red_dark = colorant"#990000"
blue = colorant"#0000ff"
green = colorant"#66ff66"

function plot_g(g, N, s, β, c)
    p = plot(size = (500, 500), tickfontsize = 6, legend_position = (0.15,0.20),
    grid = false, xlabel = "\$r\$", ylabel = "\$g(r)\$")
    title = latexstring("Pair correlation (\$s = $(s)\$)")
    title!(p, titlefontsize = 9, subplot = 1, title)
    
    nbins = length(g)
    xs = [N/(2*nbins)*(i-1/2) for i=1:nbins]
    plot!(p, xs, g, color = c, linewidth = 1, label = "\$\\beta = $β \$")
    return p
end 


function plot_g_multiple(g, N, s, temperatures; n1 = N/4, n2 = 3*N/8)
    p = plot(size = (500, 500), tickfontsize = 6, legend_position = (0.15,0.20),
    grid = false, xlabel = "\$r\$", ylabel = "\$g(r)\$")
    title = latexstring("Pair correlation (\$s = $(s)\$)")
    title!(p, titlefontsize = 9, subplot = 1, title)
    nbins, ntemperatures = size(g)
    xs = [N/(2*nbins)*(i-1/2) for i=1:nbins]
    for i = 1:ntemperatures
        c = get(ColorSchemes.algae,i./ntemperatures)
        plot!(p, xs, g[:, i], color = c, linewidth = 1, label = "\$T = $(temperatures[i]) \$")
    end

    # Adding inset 
    l = lens!(p, subplot = 2, [n1, n2], [0.75, 1.25], 
    inset = (1, bbox(.5, .5, 0.4, 0.4)), 
    ticks=false, xlabel = "", ylabel = "",
    framestyle=:box)
    return p
end

function plot_S(S, N, s, β, c; ymax = 4, low = 1, up = 3*N)
    p = plot(size = (500, 500), tickfontsize = 6, legend_position = (0.85,0.90),
    grid = false, xlabel = "\$k\$", ylabel = "\$S(k)\$", ylims = (0, ymax))
    title = latexstring("Structure factor (\$s = $(s)\$)")
    title!(p, titlefontsize = 9, subplot = 1, title)
    ks = low/N:1/N:up/N
    plot!(p, ks, S.(ks), markershape = :circle,
        markersize=2, markeralpha = 0.5, color = c, linewidth = 1,
        label="\$\\beta  = $β \$")
    return p
end


function plot_S(Ss::Vector{Function}, N, s, temperatures :: Vector{Float64}, lens = false; low = 1, up = 3*N, n1 = N/8, n2 = N/4)
    p = plot(size = (500, 500), tickfontsize = 6, legend_position = (0.85,0.90),
    grid = false, xlabel = "\$k\$", ylabel = "\$S(k)\$")
    title = latexstring("Structure factor (\$s = $(s)\$)")
    title!(p, titlefontsize = 9, subplot = 1, title)
    ks = low/N:1/N:up/N
    ntemperatures = length(temperatures)
    for t in 1:ntemperatures
        c = get(ColorSchemes.algae,t./ntemperatures)
        plot!(p, ks, (Ss[t]).(ks), markershape = :circle,
        markersize=2, markeralpha = 0.5, color = c, linewidth = 1,
        label="\$T = $(temperatures[t])\$")
    end
    if lens
        l = lens!(p, subplot = 2, [n1/N, n2/N], [0, 0.4], 
        inset = (1, bbox(.5, .3, 0.4, 0.4)), 
        ticks=false, xlabel = "", ylabel = "",
        framestyle=:box)
    end

    return p
end



function plot_phase_diagram(βs, ss)
    ticks = round.(1 ./ βs, digits = 2)
    p = plot(size = (500, 500), tickfontsize = 5, legend =:false,
    grid = false, xlabel = "\$s\$", ylabel = "\$T_{s}\$", aspect_ratio = :auto,
    ymirror = true, widen =:false, yticks = ticks, xticks = ss,
    tick_direction =:out, tickfonthalign = :right, title="Transition BKT")
    title!(p, titlefontsize = 9, subplot = 1,"BKT phase transition")
    plot!(p, ss, 1 ./βs, c = :black, markershape = :star, markersize = 6,fill = true,
    fillopacity = 0.2, fillcolor = :red)
    for j in eachindex(ss)
        plot!(p, [0, ss[j]], [1/βs[j], 1/βs[j]], c =:black,
        ls = :dot)
        plot!(p, [ss[j], ss[j]], [0, 1/βs[j]], c =:black,
        ls = :dot)
    end
    return p
end 

using Plots
# default(aspect_ratio = :equal)

function plot_generate(size = 500; legend = false, ylim = false)
    # function plot_generate(size = 500; legend = false, ylim = (-ymax-0.1,ymax+0.1))
    if ylim == false
        return plot(size = (size,size), aspect_ratio = :equal , legend = legend)
    else
        return plot(size = (size,size), aspect_ratio = :equal , legend = legend, ylim = ylim)
    end
end

function plot_edges(p, g, linecolor = "pink", linewidth = 2)
    x = real(g.pts)
    y = imag(g.pts)
    xall = []
    yall = []
    for edges in g.edges
        i, j = edges
        push!(xall, [x[i], x[j]])
        push!(yall, [y[i], y[j]])
    end
    plot!( xall, yall, linecolor=linecolor, linewidth=linewidth, label="")

    # plot!( [x[i], x[j]],[y[i] , y[j]], linecolor=linecolor, linewidth=linewidth, label="")
    return p
end

function plot_sites(p, g, mc=:blue, markerstrokewidth=0; annotate = false, markersize=  6)
    if annotate
        scatter!(p, real(g.pts), imag(g.pts), mc=mc, markerstrokewidth=markerstrokewidth, series_annotations = text.(1:g.n), label = "",  markersize=markersize)
    else
        scatter!(p, real(g.pts), imag(g.pts), mc=mc, markerstrokewidth=markerstrokewidth, label = "",  markersize=markersize)
    end
    return p
end

function plot_single_point(p, point, mc=:magenta, markersize=6)
    scatter!(p, [real(point)], [imag(point)], mc=mc, markerstrokewidth=0, label="", markersize=markersize)
    return p
end

function plot_graph(p,g, edge_color="pink", linewidth = 2, point_color=:blue)
        p = plot_edges(p,g, edge_color, linewidth)
        p = plot_sites(p, g, point_color)
    return p
end

function plot_membrane(p, membrane_edges) 
    #membrane_edges = get_membrane_edges_theory()
    for l in membrane_edges
        i, j = l
        plot!(p, [pts_x[i], pts_x[j]],[pts_y[i] , pts_y[j]], linecolor="brown", linewidth=5, label="")
    end
    return p
end

function plot_config(p, config,spin_to_G, so_lat, annotate = false)
    for i in 1:length(config)
        if config[i] == true
            a, b = spin_to_G[i]
            p1, p2 = so_lat.pts[a], so_lat.pts[b] 
            x, y = real.([p1,p2]) ,  imag.([p1, p2])
            plot!(p,x, y , linecolor="teal", linewidth=4)
        end
    end
    if annotate
        scatter!(p, real(pts), imag(pts), mc=:blue, markerstrokewidth=0, series_annotations = text.(1:npts))
    end
    return p
end

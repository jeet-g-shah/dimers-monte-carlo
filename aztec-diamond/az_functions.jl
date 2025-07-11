
function aztec_pts(N::Int64)
    pts = Complex{Float64}[]
    for i in -N:N
        for j in -N:N
            pt = 1/2 + i + 1/2*1im + 1im *j
            if abs(real(pt)) + abs(imag(pt)) <= N-2
                push!(pts, pt)
            end
        end
    end
    return pts
end

function aztec_adj_edges(pts::Vector{Complex{Float64}})
    adj, edges = get_adj(1, pts)
    return adj, edges
end

function aztec_plaq(pts, N)
    asq1 = 1.0
    asq2 = 1.0im
    square = [-0.5 - 0.5im, 0.5 - 0.5im,  0.5 + 0.5im , -0.5 + 0.5im]
    plaq = Vector{Int64}[]
    for i in -N:N
        for j in -N:N
            curr_square = square .+ asq1 * i .+ asq2 * j
            indices = [index_pts(pts, curr_square[i]) for i in 1:length(square)]
            if !(0 ∈ indices)
                push!(plaq, indices)
            end
        end
    end
    return vecvec_to_mat(plaq)
end

function generate_aztec(N::Int64)
    pts = aztec_pts(N)
    adj, edges = aztec_adj_edges(pts)
    nbr = get_nbr(adj)
    plaq = aztec_plaq(pts, N)
    npts = length(pts)
    az_lat = Graph(pts, adj, nbr, npts, plaq, fill(0,1,1), edges)
    return az_lat
end

function aztec_spin_plaq(spin_pts, N)
    asq1 = 1.0
    asq2 = 1.0im
    square = [0.5 , 0.5im,  -0.5 , -0.5im]
    plaq = Vector{Int64}[]
    for i in -N:N
        for j in -N:N
            curr_square = square .+ asq1 * i .+ asq2 * j
            indices = [index_pts(spin_pts, curr_square[i]) for i in 1:length(square)]
            if !(0 ∈ indices)
                push!(plaq, indices)
            end
        end
    end
    return vecvec_to_mat(plaq)
end

function generate_spin_lat(az_lat::Graph, N::Int64)
    npts = az_lat.n
    spin_pts = Complex{Float64}[]
    spin_to_G = Vector{Int64}[]
    for edge in az_lat.edges
        site1, site2 = edge
        push!(spin_pts, (az_lat.pts[site1] + az_lat.pts[site2])/2)
        push!(spin_to_G, [site1, site2])
    end
    
    nspin = length(spin_pts) 
    spin_nbr = [Int64[] for i in 1:nspin]
    spin_adj = fill(0,nspin, nspin)
    spin_edges = Vector{Int64}[]

    for i in 1:nspin-1
        for j in i+1:nspin
            if abs(spin_pts[i] - spin_pts[j]) < 1/sqrt(2) + 1e-10
                push!(spin_nbr[i], j)
                push!(spin_nbr[j], i)
                spin_adj[i,j] = 1
                spin_adj[j,i] = 1
                push!(spin_edges, [i,j])
            end
        end
    end
    
    plaq_sq_spin = aztec_spin_plaq(spin_pts, N)
    

    spins = Graph(spin_pts, spin_adj, spin_nbr, nspin, plaq_sq_spin, fill(0,1,1), spin_edges)
    return spins, spin_to_G
end

function generate_config_vert(az_lat::Graph, az_lat_spin::Graph, spin_to_G::Vector{Vector{Int64}})
    config = fill(false, az_lat_spin.n)
    for (i,pos) in enumerate(az_lat_spin.pts)
        site1, site2 = spin_to_G[i]
        site1_pos, site2_pos = az_lat.pts[site1], az_lat.pts[site2]
        x, y = real(pos), imag(pos)
        xmax = maximum(real.(pos))
        if (Int(round(abs(x) + abs(y) - 0.5)) % 2 == 1) && (real(site1_pos) == real(site2_pos))
            config[i] = true
        end  
    end
    return config
end

function generate_config_hori(az_lat::Graph, az_lat_spin::Graph, spin_to_G::Vector{Vector{Int64}})
    config = fill(false, az_lat_spin.n)
    for (i,pos) in enumerate(az_lat_spin.pts)
        site1, site2 = spin_to_G[i]
        site1_pos, site2_pos = az_lat.pts[site1], az_lat.pts[site2]
        x, y = real(pos), imag(pos)
        xmax = maximum(real.(pos))
        if (Int(round(abs(x) + abs(y) - 0.5)) % 2 == 1) && (imag(site1_pos) == imag(site2_pos))
            config[i] = true
        end  
    end
    return config
end

function az_randomize_config(nstep::Int64, config::Vector{Bool}, plaq_sq_spin::Matrix{Int64})
    nsq = size(plaq_sq_spin)[2]
    for i in 1:nstep
        while plaq_sq_move(config,plaq_sq_spin, nsq ) 
        end
    end
    return config
end

function plot_arctic_circle(p, N::Int64)
    r = (N-2)/sqrt(2)             # Radius of the circle
    theta = range(0, 2π, length=100)  # Angle range

    # Parametric equations
    x = r * cos.(theta)
    y = r * sin.(theta)

    # Plot the circle
    p = plot!(p, x, y, aspect_ratio=1, label="Circle", lw=2, color=:black)
    return p
end
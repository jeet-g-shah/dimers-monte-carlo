

function get_pts(N::Int64)
    a1 = 1.0 + 0.0im
    a2 = 0.0 + 1.0im
    pts = ComplexF64[]
    for i in 1:N
        for j in 1:N
            push!(pts, a1 * i + a2 * j )
        end
    end
    return pts, length(pts)
end

function get_adj(N::Int64, pts::Vector{ComplexF64}; pbc = false)
    npts = length(pts)
    adj = fill(0,npts, npts)
    edges = Vector{Int64}[]
    for i in 1:npts-1
        for j in i+1: npts
            if abs(abs(pts[i] - pts[j]) - 1) < 1e-5
                adj[i,j] = 1
                adj[j,i] = 1
                push!(edges, [i,j])
            end
        end
    end
    if pbc
        for i in 1:npts-1
            for j in i+1:npts
                if (abs( abs(pts[i] - pts[j]) - N + 1 ) < 1e-5)
                    adj[i,j] = 1
                    adj[j,i] = 1
                    push!(edges, [i,j])
                end
            end
        end
    end
    return adj, edges
end
function get_nbr(adj::Matrix{Int64})
    nbr = Vector{Int64}[]
    npts = size(adj)[1]
    for i in 1:npts
        currnbr = Int64[]
        for j in 1:npts
            if adj[i,j] == 1
                push!(currnbr, j)
            end
        end
        push!(nbr, currnbr)
    end
    return nbr
end

function get_plaq(pts, N;pbc=false)
    asq1 = 1.0
    asq2 = 1.0im
    square = [0.0, 1.0,  1.0 + 1.0im, 1.0im]
    plaq = Vector{Int64}[]
    for i in 1:N
        for j in 1:N
            curr_square = square .+ asq1 * i .+ asq2 * j
            indices = [index_pts(pts, curr_square[i]) for i in 1:length(square)]
            if !(0 ∈ indices)
                push!(plaq, indices)
            end
        end
    end
    if pbc
        hor_periodic_sq = [0.0, N-1.0, N-1.0 + 1im, 1im] 
        for j in 1:N
            curr_square = hor_periodic_sq .+ asq1 * 1 .+ asq2 * j
            indices = [index_pts(pts, curr_square[i]) for i in 1:length(square)]
            if !(0 ∈ indices)
                push!(plaq, indices)
            end
        end
        ver_periodic_sq = [0.0, 1.0, 1.0 + (N-1.0)*1im, (N-1.0)*1im]
        for i in 1:N
            curr_square = ver_periodic_sq .+ asq1 * i .+ asq2 * 1
            indices = [index_pts(pts, curr_square[i]) for i in 1:length(square)]
            if !(0 ∈ indices)
                push!(plaq, indices)
            end
        end
        periodic_sq = [0.0, (N-1), (N-1) * (1+1im), (N-1)*1im] .+ 1.0 .+ 1im
        indices = [index_pts(pts, periodic_sq[i]) for i in 1:length(square)]
        if !(0 ∈ indices)
            push!(plaq, indices)
        end
    end
    # return plaq
    return vecvec_to_mat(plaq)
end

function get_sq_lat(N::Int64;pbc=false)
    pts, npts = get_pts(N)
    adj, edges = get_adj(N, pts, pbc=pbc)
    nbr = get_nbr(adj)
    plaq = get_plaq(pts, N, pbc=pbc)
    sq_lat = Graph(pts, adj, nbr, npts, plaq, fill(0,1,1), edges)
    return sq_lat
end

function get_vert_col_state(sq_lat_spin::Graph, N::Int64)
    config = fill(false,sq_lat_spin.n)
    for i in 1:sq_lat_spin.n
        spin_pt = sq_lat_spin.pts[i]
        x, y = real(spin_pt), imag(spin_pt)
        if (abs(modf(x)[1]) < 1e-5) && (Int(modf(y)[2]) % 2 == 1)
            config[i] = true
        end
    end
    return config
end 

function get_hori_col_state(sq_lat_spin::Graph, N::Int64)
    config = fill(false,sq_lat_spin.n)
    for i in 1:sq_lat_spin.n
        spin_pt = sq_lat_spin.pts[i]
        x, y = real(spin_pt), imag(spin_pt)
        if (abs(modf(y)[1]) < 1e-5) && (Int(modf(x)[2]) % 2 == 1)
            config[i] = true
        end
    end
    return config
end 

function get_sq_lat_spin(sq_lat::Graph)
    latA = Int64[1]
    latB = Int64[]
    npts = sq_lat.n
    lat = fill(0, npts)
    lat[1] = 1
    cont = 0
    currpt = 1
    g = sq_lat
    nbr = g.nbr
    pts = g.pts
    stack = [1]
    while (stack != []) && (cont<10000) 
        currpt = pop!(stack)
        currnbr = nbr[currpt]
        for n in currnbr
            if lat[n] == 0
                lat[n] = -lat[currpt]
                if lat[n] > 0
                    push!(latA, n)
                elseif lat[n] < 0
                    push!(latB, n)
                end
                push!(stack, n)
            end
        end
        cont = cont + 1
    end
    sub_lat = fill(0,npts)
    for i in 1:length(latA)
        sub_lat[latA[i]] = i
    end
    for i in 1:length(latB)
        sub_lat[latB[i]] = i
    end
    spin_pts = ComplexF64[]
    spin_to_G = Vector{Int64}[]
    G_to_spin = fill(0, npts, npts)
    cnt = 1
    for i in 1:npts
        currnbr = nbr[i]
        for n in currnbr
            if G_to_spin[i,n] == 0
                # push!(spin_pts, (pts[i] + pts[n])/2)
                p1 = pts[i]
                p2 = pts[n]
                if abs( abs(p1 - p2) - 1.0) < 1e-5
                    # i and n are neighbors on open BC
                    push!(spin_pts, (p1+ p2)/2)
                else
                    # i and n are neighbors on open BC
                    dir = (abs(imag(p1 - p2)) < 1e-5) ? "x" : "y"
                    if dir == "x"
                        push!(spin_pts, maximum( real.([p1, p2]) ) + 1/2 + imag(p1)*1im )
                    elseif dir == "y"
                        push!(spin_pts, maximum( imag.([p1, p2]) )*1im  + 1im/2 + real(p1))
                    end
                end
                push!(spin_to_G,[i,n])
                G_to_spin[i,n] = cnt
                G_to_spin[n,i] = cnt
                cnt =cnt + 1
            end
        end
    end
    nspin = length(spin_pts)

    aa, bb = size(g.sq_plaq)
    plaq_sq_spin = fill(0,aa, bb)
    spin_nbr = [Int64[] for i in 1:nspin]
    spin_adj = fill(0,nspin, nspin)
    spin_edges = Vector{Int64}[]

    for i in 1:bb
        p = g.sq_plaq[:,i]
        push!(p, g.sq_plaq[1,i])
        s = [G_to_spin[p[j], p[j+1]] for j in 1:4]  
        for j in 1:4
            plaq_sq_spin[j,i] = s[j]
        end
        push!(s, s[1])
        
        for j in 1:4
            spin_adj[ s[j], s[j+1] ] = 1
            spin_adj[ s[j+1], s[j] ] = 1
            push!(spin_edges, [s[j], s[j+1]])
            push!(spin_nbr[s[j]], s[j+1])
            push!(spin_nbr[s[j+1]], s[j])
        end
        
    end

    sq_lat_spin = Graph(spin_pts, spin_adj, spin_nbr, nspin,plaq_sq_spin, fill(0,1,1), spin_edges)
    return sq_lat_spin, latA, latB, sub_lat, G_to_spin, spin_to_G
end

function torus_distance(r, N)
    x = min(abs(real(r)) ,  N - abs(real(r)))
    y = min(abs(imag(r)) ,  N - abs(imag(r)))
    return sqrt(x^2 + y^2)
end

function randomize_config_sq(nstep::Int64, config::Vector{Bool}, plaq_sq_spin::Matrix{Int64})
    nsq = size(plaq_sq_spin)[2]
    for i in 1:nstep
        while plaq_sq_move(config,plaq_sq_spin, nsq )
        end
    end
    return config
end
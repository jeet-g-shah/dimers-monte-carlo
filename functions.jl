
using BipartiteMatching
using ProgressMeter
using SparseArrays

function mean(vec)
    return sum(vec)/ length(vec)
end


cell = [1/2 + 1im/2, -1/2 + 1im/2,  -1/2 - 1im/2, 1/2 - 1im/2]

oct = [
    1/2 + 1im/2,
    1/2 + 1/sqrt(2) + 1im/2 + 1im/sqrt(2),
    3/2 + 1/sqrt(2) + 1im/2 + 1im/sqrt(2),
    3/2 + sqrt(2) + 1im/2 ,
    3/2 + sqrt(2) - 1im/2 ,
    3/2 + 1/sqrt(2) - 1im/2 - 1im/sqrt(2),
    1/2 + 1/sqrt(2) - 1im/2 - 1im/sqrt(2),
    1/2 - 1im/2,
]


a1 = (1 + sqrt(2))*(1/sqrt(2)+ 1im/sqrt(2))
a2 = (1 + sqrt(2))*(1/sqrt(2) -1im/sqrt(2))

struct Graph
    pts::Vector{ComplexF64}
    adj::Matrix{Int64}
    nbr::Vector{Vector{Int64}}
    n::Int64
    sq_plaq::Matrix{Int64}
    oct_plaq::Matrix{Int64}
    edges::Vector{Vector{Int64}}
end

# Base.copy(g::Graph) = Graph(g.pts, g.adj, g.nbr, g.n, g.sq_plaq, g.oct_plaq)
Base.copy(g::Graph) = Graph(g.pts, g.adj, g.nbr, g.n, g.sq_plaq, g.oct_plaq, g.edges)

function reflect_y(vec)
    return conj.(vec)
end
function reflect_x(vec)
    return reflect_y( vec .* 1im ) ./ 1im
end
function invert(vec)
    return reflect_x(reflect_y(vec))
end

function clean_up_lattice(pts::Vector{ComplexF64}, N::Int64)
    new_pts = ComplexF64[]
    xmax = maximum(real.(pts))
    xmin = minimum(real.(pts))
    ymax = maximum(imag.(pts))
    ymin = minimum(imag.(pts))
    extrm = ((1+sqrt(2)) * (N-1) / 2 + 1/sqrt(2))*(1/sqrt(2) + 1im/sqrt(2))
    extrms = [extrm + 2*a2*i for i in (-((N-1)÷2)):((N-1)÷2)]
    exclude = vcat(extrms, reflect_x(extrms), reflect_y(extrms), invert(extrms))
    # println(exclude)
    for p in pts
        skip = false
        for ex in exclude
            if abs(ex - p) < 1e-5
                skip = true
                continue
            end
        end
        if !skip
            push!(new_pts, p)
        end
    end
    return new_pts, exclude
end

function clean_up_lattice_rect(pts::Vector{ComplexF64}, N::Int64)
    new_pts = ComplexF64[]
    xmax = maximum(real.(pts)) - 0.25
    xmin = minimum(real.(pts)) + 0.25
    ymax = maximum(imag.(pts)) - 0.25 
    ymin = minimum(imag.(pts)) + 0.25
    diag_extrm1 = 1/2 + 1im/2
    diag_extrm2 = -1/2- 1im/2 + (2N - 1) * a1
    for p in pts
        if ((xmin < real(p) < xmax) && (ymin < imag(p) < ymax)) && !( abs(p-diag_extrm1) < 1e-5 ) && !( abs(p-diag_extrm2) < 1e-5 )
            push!(new_pts, p)
        end
    end

    return new_pts
end

function index_pts(pts::Vector{ComplexF64}, point::ComplexF64)
    npts = length(pts)
    for i in 1:npts
        if abs(pts[i] - point) < 1e-5
            return i
        end
    end
    return 0
end

function index_pts(pts::Vector{ComplexF64}, point_vec::Vector{ComplexF64})
    npts = length(pts)
    nvec = length(point_vec)
    nlocated = 0
    loc_in_pts = fill(0, nvec)
    status = [i for i in 1:nvec]
    for i in 1:npts
        for si in status
            if abs(pts[i] - point_vec[si]) < 1e-5
                loc_in_pts[si] = i
                status = filter(x -> x != si, status)
                nlocated = nlocated + 1
                if nlocated == nvec
                    return loc_in_pts
                end
            end
        end
    end
    return loc_in_pts
end

function vecvec_to_mat(vecvec::Vector{Vector{Int64}})
    mat = vecvec[1]
    for i in 2:length(vecvec)
        mat =hcat(mat, vecvec[i])
    end
    return mat
end



function get_pts_adj(N::Int64, cell::Vector{ComplexF64})
    pts = ComplexF64[]
    for i in 0:N-1
        for j in 0:N-1
            r = a1*(i) + a2*(j)
            for v in cell
                push!(pts, r + v)
            end
        end
    end
    cm = sum(pts) / length(pts)
    pts = pts .- [cm for i in 1:length(pts)]
    pts, exclude = clean_up_lattice(pts, N)
    npts = length(pts)
    adj = spzeros(Int64, npts, npts)
    edges = Vector{Int64}[]
    for i in 1:npts-1
        for j in i+1:npts
            if abs(abs(pts[i] - pts[j]) -1.0) < 1e-5
                adj[i,j] = 1
                adj[j,i] = 1
                push!(edges, [i,j])
            end
        end
    end
    return npts, pts, adj, edges
end





function clean_up_lattice_rect_paper(pts::Vector{ComplexF64}, N::Int64)
    new_pts = ComplexF64[]
    cutoff = N - 1/2 + N/sqrt(2) + 1e-5
    for p in pts
        if abs(real(p)) < cutoff && abs(imag(p)) < cutoff 
            push!(new_pts, p)
        end
    end

    ret_pts = ComplexF64[]
    xmax = maximum(real.(new_pts)) 

    q1 = xmax * (1 + 1im)
    q2 = xmax * (-1 + 1im)
    q3 = -xmax * (1 + 1im)
    q4 = -xmax * (-1 + 1im)
    for p in new_pts
        if !( abs(p-q1) < 1e-5 ) && !( abs(p-q2) < 1e-5 ) && !( abs(p-q3) < 1e-5 ) && !( abs(p-q4) < 1e-5 )
            push!(ret_pts, p)
        end
    end
    return ret_pts
end



function get_pts_adj_rect_paper(N::Int64, cell::Vector{ComplexF64})
    b1 = 1 + sqrt(2) 
    b2 = 1im * b1

    pts1 = ComplexF64[]
    for i in -N:N
        for j in -N:N
            r = b1*(i) + b2*(j)
            for v in cell
                push!(pts1, r + v)
            end
        end
    end
    pts = pts1 .* (1/sqrt(2) + 1im/sqrt(2))
    cm = sum(pts) / length(pts)
    pts = pts .- [cm for i in 1:length(pts)]
    pts = clean_up_lattice_rect_paper(pts, N)
    npts = length(pts)

    adj = spzeros(Int64, npts, npts)

    edges = Vector{Int64}[]
    for i in 1:npts-1
        for j in i+1:npts
            if abs(abs(pts[i] - pts[j]) -1.0) < 1e-5
                adj[i,j] = 1
                adj[j,i] = 1
                push!(edges, [i,j])
            end
        end
    end

    return npts, pts, adj, edges
end

function get_pts_adj_rect(N::Int64, cell::Vector{ComplexF64})
    b1 = 2 + sqrt(2)
    b2 = 1im * b1
    
    pts1 = ComplexF64[]
    for i in 0:N-1
        for j in 0:N-1
            r = b1*(i) + b2*(j)
            for v in cell
                push!(pts1, r + v)
            end
        end
    end
    pts2 = pts1 .+ (1 + sqrt(2))*(1/sqrt(2)+ 1im/sqrt(2))
    pts = vcat(pts1, pts2)
    
    pts = clean_up_lattice_rect(pts, N)
    cm = sum(pts) / length(pts)
    pts = pts .- [cm for i in 1:length(pts)]
    npts = length(pts)

    adj = spzeros(Int64, npts, npts)

    edges = Vector{Int64}[]
    for i in 1:npts-1
        for j in i+1:npts
            if abs(abs(pts[i] - pts[j]) -1.0) < 1e-5
                adj[i,j] = 1
                adj[j,i] = 1
                push!(edges, [i,j])
            end
        end
    end
    
    return npts, pts, adj, edges
end


function get_sq_oct_plaq(pts::Vector{ComplexF64}, cell::Vector{ComplexF64}, oct::Vector{ComplexF64}, N::Int64, shift= 0.0 + 0.0im)
    # -------------------------------------------
    # Sort points into bin to speed up finding index of a point.
    min_re_pts = minimum(real.(pts))
    max_re_pts = maximum(real.(pts))
    ϵ = 1e-5
    nbins = Int64(floor(max_re_pts - min_re_pts + ϵ)) + 1

    ord_pts = [ComplexF64[] for i in 1:nbins]
    ord_ind = [Int64[] for i in 1:nbins]
    for i in 1:length(pts)
        pt = pts[i]
        bin_no = Int64( floor(real(pt) - min_re_pts + ϵ) ) + 1
        push!(ord_pts[bin_no], pt)
        push!(ord_ind[bin_no], i)
    end
    # -------------------------------------------

    
    pts_sq = []
    for i1 in -(N-1) : (N-1) 
        for i2 in -(N-1) : (N-1) 
            push!(pts_sq, cell .+ a1*i1 .+ a2*i2 .+ shift)
        end
    end

    plaq_sq = Vector{Int64}[]
    for p in pts_sq
        valid_plaq = true
        plaq_ind = Int64[]
        for pt in p 
            bin_no = Int64( floor(real(pt) - min_re_pts + ϵ) ) + 1
            if bin_no < 1 || bin_no > nbins
                valid_plaq = false
                break
            end
            bin_ind = index_pts(ord_pts[bin_no], pt)
            if bin_ind == 0
                valid_plaq = false
                break
            end
            ind = ord_ind[bin_no][bin_ind]
            push!(plaq_ind, ind)
        end
        if valid_plaq
            push!(plaq_sq, plaq_ind)
        end
    end

    pts_oct = []
    for i1 in -(N) : (N) 
        for i2 in -(N) : (N) 
            push!(pts_oct, oct .+ a1*i1 .+ a2*i2 .+ shift)
        end
    end

    plaq_oct = Vector{Int64}[]

    for p in pts_oct
        valid_plaq = true
        plaq_ind = Int64[]
        for pt in p 
            bin_no = Int64( floor(real(pt) - min_re_pts + ϵ) ) + 1
            if bin_no < 1 || bin_no > nbins
                valid_plaq = false
                break
            end
            bin_ind = index_pts(ord_pts[bin_no], pt)
            if bin_ind == 0
                valid_plaq = false
                break
            end
            ind = ord_ind[bin_no][bin_ind]
            push!(plaq_ind, ind)
        end
        if valid_plaq
            push!(plaq_oct, plaq_ind)
        end
    end

    vec_plaq_sq, vec_plaq_oct = vecvec_to_mat(plaq_sq), vecvec_to_mat(plaq_oct)
    return vec_plaq_sq, vec_plaq_oct 
end


function get_nbr(adj::Union{Matrix{Int64}, SparseMatrixCSC{Int64, Int64}})
    npts = size(adj)[1]
    nbr = [Int64[] for i in 1:npts]
    for i in 1:npts
        for j in 1:npts
            if adj[i,j] == 1
                push!(nbr[i], j)
            end
        end
    end
    return nbr
end



function get_so_lat_rect(N::Int64, oct::Vector{ComplexF64}; paper=false)

    cell = [1/2 + 1im/2, -1/2 + 1im/2,  -1/2 - 1im/2, 1/2 - 1im/2]
    cell_paper = [1.0 + 0im, 0 + 1im, -1 + 0im, 0 - 1im] / sqrt(2)

    # if mod(N,2) != 0
    #     println("Change N")
    # end

    if paper
        npts, pts, adj, edges = get_pts_adj_rect_paper(N, cell_paper); # Fast step
    else
        npts, pts, adj, edges = get_pts_adj_rect(N, cell); # Fast step
    end

    plaq_sq, plaq_oct = get_sq_oct_plaq(pts, cell, oct, N,  0.0+0.0im ) # Slowest step
    nbr = get_nbr(adj) # Fast step
    so_lat_rect = Graph(pts, adj, nbr, npts,plaq_sq, plaq_oct, edges)

    # Check if the lattice is agrees with expected numbers.
    nvertices_expected = 4*(N^2 +(N-1)^2) + 8*(N-1) 
    nsq_expected = (N)^2 + (N-1)^2 
    noct_expected = nsq_expected - 1

    nsq = size(so_lat_rect.sq_plaq)[2]
    noct = size(so_lat_rect.oct_plaq)[2]

    if (so_lat_rect.n != nvertices_expected)
        printstyled("Vertices not as expected.\n";color=:red)
    end
    if (nsq != nsq_expected)
        printstyled("Number of squares not as expected.\n";color=:red)
    end
    if (noct != noct_expected)
        printstyled("Number of octagons not as expected.\n";color=:red)
    end

    return so_lat_rect
end



function get_so_lat(N::Int64, cell::Vector{ComplexF64}, oct::Vector{ComplexF64})
    @time npts, pts, adj, edges = get_pts_adj(N, cell);
    @time plaq_sq, plaq_oct = get_sq_oct_plaq(pts, cell, oct, N)
    @time nbr = get_nbr(adj)
    so_lat = Graph(pts, adj, nbr, npts,plaq_sq, plaq_oct, edges)
    return so_lat
end

function get_perfect_covering(latA::Vector{Int64}, latB::Vector{Int64},sub_lat::Vector{Int64}, nbr::Vector{Vector{Int64}} )
    
    bi_adj = falses(length(latA), length(latB))
    for i in 1:length(latA)
        A = latA[i]
        for B in nbr[A]
            j = sub_lat[B]
            bi_adj[i,j] = true
        end
    end
    matching, matched = findmaxcardinalitybipartitematching(bi_adj)
    covering = []
    for A in 1:length(latA)
        B = matching[A]
        edg = [latA[A], latB[B]]
        push!(covering, edg)
    end
    return covering
end

function convering_to_config(covering, so_lat_spin, G_to_spin)
    config = fill(0,so_lat_spin.n)
    for e in covering
        p1, p2 = e
        s = G_to_spin[p1, p2]
        config[s] = 1
    end
    bin_config = convert(Vector{Bool}, config);
    return bin_config
end

function reset_config(covering, so_lat_spin, G_to_spin)
    return convering_to_config(covering, so_lat_spin, G_to_spin)
end

function get_so_lat_spin(so_lat::Graph)
    latA = Int64[1]
    latB = Int64[]
    npts = so_lat.n
    lat = fill(0, npts)
    lat[1] = 1
    cont = 0
    currpt = 1
    g = so_lat
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
    G_to_spin = spzeros(Int64,npts, npts)
    cnt = 1
    for i in 1:npts
        currnbr = nbr[i]
        for n in currnbr
            if G_to_spin[i,n] == 0
                push!(spin_pts, (pts[i] + pts[n])/2)
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
    spin_adj = spzeros(Int64,nspin, nspin)
    spin_edges = Vector{Int64}[]

    for i in 1:size(g.sq_plaq)[2]
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

    aa, bb = size(g.oct_plaq)
    plaq_oct_spin = fill(0,aa, bb)

    for i in 1:size(g.oct_plaq)[2]
        p = g.oct_plaq[:,i]
        push!(p, g.oct_plaq[1,i])
        s = [G_to_spin[p[j], p[j+1]] for j in 1:8]    
        for j in 1:8
            plaq_oct_spin[j,i] = s[j]
        end
        push!(s, s[1])
        for j in 1:8
            spin_adj[ s[j], s[j+1] ] = 1
            spin_adj[ s[j+1], s[j] ] = 1
            push!(spin_edges, [s[j], s[j+1]])
            push!(spin_nbr[s[j]], s[j+1])
            push!(spin_nbr[s[j+1]], s[j])
        end
    end
    so_lat_spin = Graph(spin_pts, spin_adj, spin_nbr, nspin,plaq_sq_spin, plaq_oct_spin, spin_edges)
    return so_lat_spin, latA, latB, sub_lat, G_to_spin, spin_to_G
end

function randomize_config(nstep::Int64, config::Vector{Bool}, plaq_sq_spin::Matrix{Int64}, plaq_oct_spin::Matrix{Int64})
    nsq = size(plaq_sq_spin)[2]
    noct = size(plaq_oct_spin)[2]
    # printstyled("Using either a plaq move or nothing as a step. This is the correct method. \n";color=:green)
    # printstyled("Randomized algorithm need to be fixed\n";color=:red)

    for i in 1:nstep
        # while !( (plaq_oct_move(config,plaq_oct_spin, noct ) ==false) || (plaq_sq_move(config,plaq_sq_spin, nsq ) == false) )
        # end
        if rand(Bool)
            plaq_oct_move(config,plaq_oct_spin, noct )
        else
            plaq_sq_move(config,plaq_sq_spin, nsq )
        end
    end
    return config
end

function plaq_sq_move(config::Vector{Bool}, plaq_sq_spin::Matrix{Int64}, nsq::Int64) #rng::UnitRange{Int64})
    # Returns true if the plaquette is not flippable
    # Returns false if the plaquette if flippable. Also flips the plaquette in this case
    i = rand(1:nsq)
    @inbounds s1, s2, s3, s4 = plaq_sq_spin[1,i], plaq_sq_spin[2,i], plaq_sq_spin[3,i], plaq_sq_spin[4,i]
    @inbounds cs1, cs2, cs3, cs4 = config[s1], config[s2], config[s3], config[s4]
    if cs1 == cs2
        return true
    elseif cs2 == cs3
        return true
    elseif cs3 == cs4
        return true
    elseif cs4 == cs1
        return true
    end
    @inbounds config[s1] = cs2
    @inbounds config[s2] = cs1
    @inbounds config[s3] = cs4
    @inbounds config[s4] = cs3
    return false
end

function plaq_oct_move(config::Vector{Bool}, plaq_oct_spin::Matrix{Int64}, noct::Int64) #rng::UnitRange{Int64})
    # Returns true if the plaquette is not flippable
    # Returns false if the plaquette if flippable. Also flips the plaquette in this case
    i = rand(1:noct)
    @inbounds s1, s2, s3, s4, s5, s6, s7, s8 = plaq_oct_spin[1,i], plaq_oct_spin[2,i], plaq_oct_spin[3,i], plaq_oct_spin[4,i], plaq_oct_spin[5,i], plaq_oct_spin[6,i], plaq_oct_spin[7,i], plaq_oct_spin[8,i]
    @inbounds cs1, cs2, cs3, cs4, cs5, cs6, cs7, cs8 = config[s1], config[s2], config[s3], config[s4], config[s5], config[s6], config[s7], config[s8]
    if cs1 == cs2 || cs2 == cs3 ||  cs3 == cs4 || cs4 == cs5 || cs5 == cs6 || cs6 == cs7 || cs7 == cs8 || cs8 == cs1 
        # Plaquette is not flippable
        return true
    end    
    # Plaquette is flippable
    @inbounds config[s1] = cs2
    @inbounds config[s2] = cs1
    @inbounds config[s3] = cs4
    @inbounds config[s4] = cs3
    @inbounds config[s5] = cs6
    @inbounds config[s6] = cs5
    @inbounds config[s7] = cs8
    @inbounds config[s8] = cs7
    return false
end

function conv_to_arr(b::String, nline::Int64)
    arr = Bool[]
    for i in 1:length(b)
        if b[i] == '\0'
            push!(arr,false )
        else
            push!(arr, true)
        end
    end
    arr = reshape(arr, (nline, length(b) ÷ nline))
    return arr
end




function vis_frm_arr(arr::Matrix{Bool})
    l, t = size(arr)
    vis_list = Float64[]
    for i in 1:l
        B = @views arr[1:i,:]
        vis = abs( sum((-1) .^ (sum(B ; dims = 1)))) / t
        push!(vis_list, vis)
    end
    return vis_list
end

function vis_from_arr_ij(arr::Matrix{Bool}, i::Int64, j::Int64)
    i, j = min(i, j), max(i, j)
    l, t = size(arr)
    B = @views arr[i:j,:]
    vis = abs( sum((-1) .^ (sum(B ; dims = 1)))) / t
    return vis
end

function szsz_frm_arr(arr::Matrix{Bool}, i, j)
    si = @views arr[i, :] 
    sj = @views arr[j, :]
    szsz = sum( si .& sj )/size(arr)[2]
    return szsz
end 
function sz_frm_arr(arr::Matrix{Bool}, i)
    si = @views arr[i, :] 
    sz = sum( si )/size(arr)[2]
    return sz
end 

function linear_fit(x, y)
    Lx = length(x)
    Ly = length(y)
    L = Lx
    if Lx != Ly
        println("X and Y have different dimensions")
    end
    x_bar = sum(x)/L
    y_bar = sum(y)/L
    b = sum( (x .- x_bar) .* (y .- y_bar) ) / sum( (x .- x_bar) .^ 2 )
    a = y_bar - b * x_bar
    return [a, b]
end 

function get_vert_col_state(so_lat_rect_spin, so_lat_rect)
    config = fill(false, so_lat_rect_spin.n)
    ymax = maximum(imag.(so_lat_rect.pts))
    ymin = minimum(imag.(so_lat_rect.pts))
    for i in 1:so_lat_rect_spin.n
        s1 = so_lat_rect_spin.pts[i]
        p1, p2 = spin_to_G[i]
        vec = so_lat_rect.pts[p1] - so_lat_rect.pts[p2]
        if (abs( real.(vec) ) < 1e-5) || (abs( imag.(s1) - ymax) < 1e-5) || (abs( imag.(s1) - ymin) < 1e-5) 
            config[i] = true
        end
    end
    return config
end

function get_hori_col_state(so_lat_rect_spin, so_lat_rect)
    config = fill(false, so_lat_rect_spin.n)
    xmax = maximum(real.(so_lat_rect.pts))
    xmin = minimum(real.(so_lat_rect.pts))
    for i in 1:so_lat_rect_spin.n
        s1 = so_lat_rect_spin.pts[i]
        p1, p2 = spin_to_G[i]
        vec = so_lat_rect.pts[p1] - so_lat_rect.pts[p2]
        if (abs( imag.(vec) ) < 1e-5) || (abs( real.(s1) - xmax) < 1e-5) || (abs( real.(s1) - xmin) < 1e-5) 
            config[i] = true
        end
    end
    return config
end

function ordered_state(so_lat_spin::Graph, N::Int64, order::String)
    # This function generate the two extremela ordered states
    # The two ordering could be horizontal or vertical
    green_plus = [ (1/2 + 1/2/sqrt(2)) * i for i in  1 : 4 : 2*(N-1)+1]
    green_minus = [ (-1/2 - 1/2/sqrt(2)) * i for i in  1 : 4 : 2*(N-1)+1]
    green_list = vcat(green_plus, green_minus)
    
    if order == "vert"
        pts_relevant_coord = real.(so_lat_spin.pts)
    elseif order == "hori"
        pts_relevant_coord = imag.(so_lat_spin.pts)
    else
        error("order must be either 'vert' or 'hori'")
    end
    green_spins = [findall(p -> abs(p - green_list[i]) < 1e-5, pts_relevant_coord) for i in 1:length(green_list)]

    green_spins = vcat(green_spins...)

    orange_plus1 = [(1/2 + 1/2/sqrt(2)) * i - 1/2/sqrt(2) for i in  3 : 4 : 2*N+2]
    orange_plus2 = [(1/2 + 1/2/sqrt(2)) * i + 1/2/sqrt(2) for i in  3 : 4 : 2*N+2]
    orange_minus1 = [ (-1/2 - 1/2/sqrt(2)) * i + 1/2/sqrt(2) for i in  3 : 4 : 2*N+2]
    orange_minus2 = [ (-1/2 - 1/2/sqrt(2)) * i - 1/2/sqrt(2) for i in  3 : 4 : 2*N+2]
    orange_list = vcat(orange_plus1, orange_plus2, orange_minus1, orange_minus2)
    orange_spins = [findall(p -> abs(p - orange_list[i]) < 1e-5, pts_relevant_coord) for i in 1:length(orange_list)]
    orange_spins = vcat(orange_spins...)

    config = fill(false, so_lat_spin.n)
    config[green_spins] .= true
    config[orange_spins] .= true
    return config
end

    

function sort_path_real(path, so_lat_rect_spin)
    path_pts = so_lat_rect_spin.pts[path]
    path_pts_re = real.(path_pts)
    perm = sortperm(path_pts_re)
    path_sorted = path[perm]
    return path_sorted
end

function avg(x...)
    mn = zeros(length(x[1]))
    for i in 1:length(x)
        mn = mn .+ x[i]
    end
    return mn ./ length(x)
end
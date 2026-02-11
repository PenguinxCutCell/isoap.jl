"""
    isoap(poly::Polyhedron, phi::AbstractVector{Float64}, phiiso::Float64) -> IsoResult

Extract an isosurface from a general polyhedron.

# Arguments
- `poly`: the polyhedron
- `phi`: scalar values at each vertex of the polyhedron (length == nvertices(poly))
- `phiiso`: the iso-value

# Returns
An `IsoResult` containing the iso-polygons, iso-vertex coordinates, and the
face indices over which each iso-edge is constructed.
"""
function isoap(poly::Polyhedron, phi::AbstractVector{<:Real}, phiiso::Real)
    ntp = nvertices(poly)
    nts = nfaces(poly)

    # Tag vertices: 1 if phi > phiiso, 0 otherwise
    ia = fill(-1, ntp)
    icontp = 0
    icontn = 0

    for is in 1:nts
        for ip in poly.ipv[is]
            if ia[ip] == -1
                if phi[ip] > phiiso
                    ia[ip] = 1
                    icontp += 1
                else
                    ia[ip] = 0
                    icontn += 1
                end
            end
        end
    end

    # If all vertices are on the same side, no isosurface
    if icontp == 0 || icontn == 0
        return IsoResult(Vector{Int}[], Int[], NTuple{3,Float64}[])
    end

    # Insert and arrange the iso-vertices
    ipviso, isoeface_vec, ipia0, ipia1, nipnew = isopol(ia, poly.ipv, nts)

    if isempty(ipviso)
        return IsoResult(Vector{Int}[], Int[], NTuple{3,Float64}[])
    end

    # Compute iso-vertex coordinates by linear interpolation
    vertiso = Vector{NTuple{3,Float64}}(undef, nipnew)
    for iso_poly in ipviso
        for ip in iso_poly
            ip0 = ipia0[ip]
            ip1 = ipia1[ip]
            t = (phiiso - phi[ip0]) / (phi[ip1] - phi[ip0])
            x0, y0, z0 = poly.vertp[ip0]
            x1, y1, z1 = poly.vertp[ip1]
            vertiso[ip] = (
                x0 + t * (x1 - x0),
                y0 + t * (y1 - y0),
                z0 + t * (z1 - z0)
            )
        end
    end

    return IsoResult(ipviso, isoeface_vec, vertiso)
end

"""
    isopol(ia, ipv, nts) -> (ipviso, isoeface, ipia0, ipia1, nipnew)

Insert and arrange iso-vertices on the faces of the polyhedron cut by the isosurface.

Internal function used by `isoap`.
"""
function isopol(ia::Vector{Int}, ipv::Vector{Vector{Int}}, nts::Int)
    # Determine which faces are cut and count intersected edges per face
    nedge = zeros(Int, nts)
    iscut = zeros(Int, nts)
    niscut = 0

    for is in 1:nts
        nipv_is = length(ipv[is])
        if nipv_is > 0
            for iv in 1:nipv_is
                ip = ipv[is][iv]
                iv1 = iv == nipv_is ? 1 : iv + 1
                ip1 = ipv[is][iv1]
                if ia[ip] != ia[ip1]
                    iscut[is] = 1
                    nedge[is] += 1
                end
            end
            if iscut[is] == 1
                niscut += 1
            end
        end
    end

    if niscut == 0
        return (Vector{Int}[], Int[], Int[], Int[], 0)
    end

    # Iso-vertex insertion
    # For each cut face, find intersection points on edges and record them
    nipnew = 0
    ipia0 = Dict{Int,Int}()  # iso-vertex -> polyhedron vertex with ia=0
    ipia1 = Dict{Int,Int}()  # iso-vertex -> polyhedron vertex with ia=1
    
    # ISE[is][edge_count] = iso-vertex index for the edge_count-th cut edge of face is
    ise = Dict{Int,Vector{Int}}()
    # IPVINT[is] = list of iso-vertices for face is (in face-traversal order)
    ipvint = Dict{Int,Vector{Int}}()
    # IVISE[is][ipnew] = position within ipvint[is] of iso-vertex ipnew
    ivise = Dict{Int,Dict{Int,Int}}()
    # IPISE[ipnew] = (face_for_type1, face_for_type2) 
    # type 1: edge goes from ia=0 to ia=1, type 2: edge goes from ia=1 to ia=0
    ipise = Dict{Int,Vector{Int}}()

    for is in 1:nts
        if iscut[is] == 1
            niv = 0
            nint = 0
            ise[is] = Int[]
            ipvint[is] = Int[]
            ivise[is] = Dict{Int,Int}()

            nipv_is = length(ipv[is])
            for iv in 1:nipv_is
                ip = ipv[is][iv]
                iv1 = iv == nipv_is ? 1 : iv + 1
                ip1_edge = ipv[is][iv1]

                if ia[ip] != ia[ip1_edge]
                    nint += 1
                    niv += 1
                    if ia[ip] == 1
                        ip1i = ip
                        ip0i = ip1_edge
                        itype = 2
                    else
                        ip1i = ip1_edge
                        ip0i = ip
                        itype = 1
                    end

                    # Check if this edge was already seen in a previous face
                    found = false
                    for is1 in 1:is-1
                        if haskey(ise, is1)
                            for ie in 1:length(ise[is1])
                                ipnew_candidate = ise[is1][ie]
                                ip0n = ipia0[ipnew_candidate]
                                ip1n = ipia1[ipnew_candidate]
                                if ip0n == ip0i && ip1n == ip1i
                                    push!(ise[is], ipnew_candidate)
                                    push!(ipvint[is], ipnew_candidate)
                                    ivise[is][ipnew_candidate] = niv
                                    if !haskey(ipise, ipnew_candidate)
                                        ipise[ipnew_candidate] = [0, 0]
                                    end
                                    ipise[ipnew_candidate][itype] = is
                                    found = true
                                    break
                                end
                            end
                            if found
                                break
                            end
                        end
                    end

                    if !found
                        nipnew += 1
                        ipia0[nipnew] = ip0i
                        ipia1[nipnew] = ip1i
                        push!(ipvint[is], nipnew)
                        push!(ise[is], nipnew)
                        ivise[is][nipnew] = niv
                        if !haskey(ipise, nipnew)
                            ipise[nipnew] = [0, 0]
                        end
                        ipise[nipnew][itype] = is
                    end
                end
            end
        end
    end

    # Iso-vertex arrangement: trace connected iso-polygons
    ipmark = zeros(Int, nipnew)
    ipviso_result = Vector{Vector{Int}}()
    isoeface_result = zeros(Int, nipnew)

    ivnewt = 0
    # Start from the first unmarked iso-vertex
    ipnew_start = 1
    while ipnew_start <= nipnew
        if ipmark[ipnew_start] != 0
            ipnew_start += 1
            continue
        end

        # Start a new iso-polygon
        current_polygon = Int[]
        ipini = ipnew_start
        ipnew = ipnew_start

        push!(current_polygon, ipnew)
        isoeface_result[ipnew] = ipise[ipnew][1]
        ipmark[ipnew] = 1
        ivnewt += 1

        while true
            is = ipise[ipnew][1]
            iv = ivise[is][ipnew]
            iv1 = iv - 1
            if iv1 == 0
                iv1 = length(ipvint[is])
            end
            ipnew = ipvint[is][iv1]

            if ipnew == ipini
                break
            end

            push!(current_polygon, ipnew)
            isoeface_result[ipnew] = ipise[ipnew][1]
            ipmark[ipnew] = 1
            ivnewt += 1

            if ivnewt == nipnew
                break
            end
        end

        push!(ipviso_result, current_polygon)
        ipnew_start += 1
    end

    # Convert ipia0, ipia1 dicts to vectors
    ipia0_vec = zeros(Int, nipnew)
    ipia1_vec = zeros(Int, nipnew)
    for i in 1:nipnew
        ipia0_vec[i] = ipia0[i]
        ipia1_vec[i] = ipia1[i]
    end

    return (ipviso_result, isoeface_result, ipia0_vec, ipia1_vec, nipnew)
end

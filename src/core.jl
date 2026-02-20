"""
    IsoapWorkspace(poly::Polyhedron)

Reusable scratch buffers for allocation-free `isopol!` and `isoap!`.
The workspace size is tied to the provided polyhedron topology.
"""
mutable struct IsoapWorkspace
    ntp::Int
    nts::Int
    max_face_vertices::Int
    max_edges::Int
    ia::Vector{Int}
    nedge::Vector{Int}
    iscut::Vector{Int}
    edge_map::Matrix{Int}
    face_iso::Matrix{Int}
    face_iso_count::Vector{Int}
    ivise::Matrix{Int}
    ipise::Matrix{Int}
    ipia0::Vector{Int}
    ipia1::Vector{Int}
    ipmark::Vector{Int}
    isoeface::Vector{Int}
    vertiso::Vector{NTuple{3,Float64}}
    polys::Matrix{Int}
    poly_lens::Vector{Int}
    nipnew::Int
    npoly::Int
end

function IsoapWorkspace(poly::Polyhedron)
    ntp = nvertices(poly)
    nts = nfaces(poly)
    max_face_vertices = maximum(length(face) for face in poly.ipv)
    max_edges = sum(length(face) for face in poly.ipv)
    return IsoapWorkspace(
        ntp,
        nts,
        max_face_vertices,
        max_edges,
        fill(-1, ntp),
        zeros(Int, nts),
        zeros(Int, nts),
        zeros(Int, ntp, ntp),
        zeros(Int, nts, max_face_vertices),
        zeros(Int, nts),
        zeros(Int, nts, max_edges),
        zeros(Int, max_edges, 2),
        zeros(Int, max_edges),
        zeros(Int, max_edges),
        zeros(Int, max_edges),
        zeros(Int, max_edges),
        fill((0.0, 0.0, 0.0), max_edges),
        zeros(Int, max_edges, max_edges),
        zeros(Int, max_edges),
        0,
        0,
    )
end

"""
    isopol!(ws, ia, ipv, nts) -> nothing

Allocation-free isopol arrangement using a preallocated workspace.
Results are stored in `ws`:
- `ws.npoly`, `ws.poly_lens`, `ws.polys`
- `ws.nipnew`, `ws.ipia0`, `ws.ipia1`, `ws.isoeface`
"""
function isopol!(ws::IsoapWorkspace, ia::Vector{Int}, ipv::Vector{Vector{Int}}, nts::Int)
    fill!(ws.nedge, 0)
    fill!(ws.iscut, 0)
    fill!(ws.edge_map, 0)
    fill!(ws.face_iso_count, 0)
    fill!(ws.ivise, 0)
    fill!(ws.ipise, 0)

    niscut = 0
    for is in 1:nts
        nipv_is = length(ipv[is])
        if nipv_is == 0
            continue
        end
        nedge_is = 0
        for iv in 1:nipv_is
            ip = ipv[is][iv]
            iv1 = iv == nipv_is ? 1 : iv + 1
            ip1 = ipv[is][iv1]
            if ia[ip] != ia[ip1]
                nedge_is += 1
            end
        end
        if nedge_is > 0
            ws.iscut[is] = 1
            ws.nedge[is] = nedge_is
            niscut += 1
        end
    end

    if niscut == 0
        ws.nipnew = 0
        ws.npoly = 0
        return nothing
    end

    nipnew = 0
    for is in 1:nts
        if ws.iscut[is] == 0
            continue
        end
        nint = 0
        nipv_is = length(ipv[is])
        for iv in 1:nipv_is
            ip = ipv[is][iv]
            iv1 = iv == nipv_is ? 1 : iv + 1
            ip1_edge = ipv[is][iv1]
            if ia[ip] == ia[ip1_edge]
                continue
            end

            nint += 1
            if ia[ip] == 1
                ip1i = ip
                ip0i = ip1_edge
                itype = 2
            else
                ip1i = ip1_edge
                ip0i = ip
                itype = 1
            end

            ipnew = ws.edge_map[ip0i, ip1i]
            if ipnew == 0
                nipnew += 1
                ipnew = nipnew
                ws.edge_map[ip0i, ip1i] = ipnew
                ws.ipia0[ipnew] = ip0i
                ws.ipia1[ipnew] = ip1i
            end

            ws.face_iso[is, nint] = ipnew
            ws.ivise[is, ipnew] = nint
            ws.ipise[ipnew, itype] = is
        end
        ws.face_iso_count[is] = nint
    end

    fill!(ws.ipmark, 0)
    fill!(ws.poly_lens, 0)
    fill!(ws.isoeface, 0)

    npoly = 0
    ivnewt = 0
    for ipnew_start in 1:nipnew
        if ws.ipmark[ipnew_start] != 0
            continue
        end

        npoly += 1
        ipini = ipnew_start
        ipnew = ipnew_start
        len = 1

        ws.polys[npoly, len] = ipnew
        ws.isoeface[ipnew] = ws.ipise[ipnew, 1]
        ws.ipmark[ipnew] = 1
        ivnewt += 1

        while true
            is = ws.ipise[ipnew, 1]
            iv = ws.ivise[is, ipnew]
            iv1 = iv - 1
            if iv1 == 0
                iv1 = ws.face_iso_count[is]
            end
            ipnew = ws.face_iso[is, iv1]

            if ipnew == ipini
                break
            end

            len += 1
            ws.polys[npoly, len] = ipnew
            ws.isoeface[ipnew] = ws.ipise[ipnew, 1]
            ws.ipmark[ipnew] = 1
            ivnewt += 1
            if ivnewt == nipnew
                break
            end
        end

        ws.poly_lens[npoly] = len
    end

    ws.nipnew = nipnew
    ws.npoly = npoly
    return nothing
end

"""
    isoap!(ws, poly, phi, phiiso) -> nothing

Allocation-free isosurface extraction. Results are written into `ws`.
"""
function isoap!(ws::IsoapWorkspace, poly::Polyhedron, phi::AbstractVector{<:Real}, phiiso::Real)
    ntp = nvertices(poly)
    nts = nfaces(poly)

    icontp = 0
    icontn = 0
    fill!(ws.ia, -1)

    for is in 1:nts
        for ip in poly.ipv[is]
            if ws.ia[ip] == -1
                if phi[ip] > phiiso
                    ws.ia[ip] = 1
                    icontp += 1
                else
                    ws.ia[ip] = 0
                    icontn += 1
                end
            end
        end
    end

    if icontp == 0 || icontn == 0
        ws.nipnew = 0
        ws.npoly = 0
        return nothing
    end

    isopol!(ws, ws.ia, poly.ipv, nts)
    if ws.npoly == 0
        return nothing
    end

    for ip in 1:ws.nipnew
        ip0 = ws.ipia0[ip]
        ip1 = ws.ipia1[ip]
        t = (phiiso - phi[ip0]) / (phi[ip1] - phi[ip0])
        x0, y0, z0 = poly.vertp[ip0]
        x1, y1, z1 = poly.vertp[ip1]
        ws.vertiso[ip] = (
            x0 + t * (x1 - x0),
            y0 + t * (y1 - y0),
            z0 + t * (z1 - z0),
        )
    end

    return nothing
end

"""
    isoap(poly::Polyhedron, phi::AbstractVector{Float64}, phiiso::Float64) -> IsoResult

Extract an isosurface from a general polyhedron.
"""
function isoap(poly::Polyhedron, phi::AbstractVector{<:Real}, phiiso::Real)
    ws = IsoapWorkspace(poly)
    isoap!(ws, poly, phi, phiiso)

    if ws.npoly == 0
        return IsoResult(Vector{Int}[], Int[], NTuple{3,Float64}[])
    end

    ipviso = Vector{Vector{Int}}(undef, ws.npoly)
    for i in 1:ws.npoly
        len = ws.poly_lens[i]
        poly_i = Vector{Int}(undef, len)
        for j in 1:len
            poly_i[j] = ws.polys[i, j]
        end
        ipviso[i] = poly_i
    end

    isoeface = Vector{Int}(undef, ws.nipnew)
    vertiso = Vector{NTuple{3,Float64}}(undef, ws.nipnew)
    for i in 1:ws.nipnew
        isoeface[i] = ws.isoeface[i]
        vertiso[i] = ws.vertiso[i]
    end

    return IsoResult(ipviso, isoeface, vertiso)
end

"""
    isopol(ia, ipv, nts) -> (ipviso, isoeface, ipia0, ipia1, nipnew)

Insert and arrange iso-vertices on the faces of the polyhedron cut by the isosurface.
"""
function isopol(ia::Vector{Int}, ipv::Vector{Vector{Int}}, nts::Int)
    ntp = 0
    for face in ipv
        for ip in face
            if ip > ntp
                ntp = ip
            end
        end
    end

    poly_dummy = Polyhedron(ipv, fill((0.0, 0.0, 0.0), ntp))
    ws = IsoapWorkspace(poly_dummy)
    isopol!(ws, ia, ipv, nts)

    if ws.npoly == 0
        return (Vector{Int}[], Int[], Int[], Int[], 0)
    end

    ipviso = Vector{Vector{Int}}(undef, ws.npoly)
    for i in 1:ws.npoly
        len = ws.poly_lens[i]
        poly_i = Vector{Int}(undef, len)
        for j in 1:len
            poly_i[j] = ws.polys[i, j]
        end
        ipviso[i] = poly_i
    end

    isoeface = Vector{Int}(undef, ws.nipnew)
    ipia0 = Vector{Int}(undef, ws.nipnew)
    ipia1 = Vector{Int}(undef, ws.nipnew)
    for i in 1:ws.nipnew
        isoeface[i] = ws.isoeface[i]
        ipia0[i] = ws.ipia0[i]
        ipia1[i] = ws.ipia1[i]
    end

    return (ipviso, isoeface, ipia0, ipia1, ws.nipnew)
end

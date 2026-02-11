# -----------------------------------------------------------------------
# Grid infrastructure for isosurface extraction on computational grids
# -----------------------------------------------------------------------

using Printf

"""
    Grid

Representation of a computational grid for isosurface extraction.

# Fields
- `ncell::Int`: number of grid cells
- `nface::Int`: number of grid faces
- `npoint::Int`: number of grid nodes
- `icface::Matrix{Int}`: indices of the two cells sharing each face (`nface × 2`).
  `icface[f,2] == 0` for boundary faces.
- `ipcell::Vector{Vector{Int}}`: node indices for each cell
- `ipface::Vector{Vector{Int}}`: ordered node indices for each face
- `iscell::Vector{Vector{Int}}`: face indices for each cell
- `vnode::Matrix{Float64}`: node coordinates (`npoint × 3`)
"""
struct Grid
    ncell::Int
    nface::Int
    npoint::Int
    icface::Matrix{Int}         # nface × 2
    ipcell::Vector{Vector{Int}} # ipcell[cell] = [node indices...]
    ipface::Vector{Vector{Int}} # ipface[face] = [node indices...]
    iscell::Vector{Vector{Int}} # iscell[cell] = [face indices...]
    vnode::Matrix{Float64}      # npoint × 3
end

"""
    constgrid(igrid::Int, nx::Int, ny::Int, nz::Int, xlast::Float64, ylast::Float64, zlast::Float64) -> Grid

Construct a uniform structured cubic grid (`igrid == 1`).

# Arguments
- `igrid`: grid type (must be 1 for this method)
- `nx, ny, nz`: number of cells along each coordinate axis
- `xlast, ylast, zlast`: domain extents

# Returns
A `Grid` struct.

Vertex numbering within each cell follows the isoap convention:
```
       7/----------/3
       /|         /|
      / |        / |
    8/__|______4/  |
     |  |       |  |
     |  /6------|--/2
     | /        | /
     |/_________|/
     5           1
```
"""
function constgrid(::Val{1}, nx::Int, ny::Int, nz::Int,
                   xlast::Float64, ylast::Float64, zlast::Float64)
    ncell = nx * ny * nz
    npoint = (nx + 1) * (ny + 1) * (nz + 1)
    nface = 3 * ncell + nx * ny + nx * nz + ny * nz

    dx = xlast / nx
    dy = ylast / ny
    dz = zlast / nz

    nxp = nx + 1
    nyp = ny + 1
    nzp = nz + 1

    # Node coordinates
    vnode = zeros(npoint, 3)
    for ixp in 1:nxp
        for iyp in 1:nyp
            for izp in 1:nzp
                ip = (izp - 1) * nxp * nyp + (iyp - 1) * nxp + ixp
                vnode[ip, 1] = (ixp - 1.0) * dx
                vnode[ip, 2] = (iyp - 1.0) * dy
                vnode[ip, 3] = (izp - 1.0) * dz
            end
        end
    end

    # Initialize cell data
    nipcell = fill(8, ncell)
    niscell_count = zeros(Int, ncell)
    nipface = fill(4, nface)

    icface = zeros(Int, nface, 2)
    ipcell_raw = zeros(Int, ncell, 50)   # temporary max NP
    ipface_raw = zeros(Int, nface, 50)   # temporary max NP
    iscell_raw = zeros(Int, ncell, 30)   # temporary max NF

    is = 0

    for ix in 1:nx
        for iy in 1:ny
            for iz in 1:nz
                ic = (iz - 1) * nx * ny + (iy - 1) * nx + ix

                # Cell vertices (isoap convention)
                ip1 = (iz) * nxp * nyp + (iy - 1) * nxp + (ix + 1)    # IXP=IX+1, IYP=IY, IZP=IZ+1
                ip2 = (iz - 1) * nxp * nyp + (iy - 1) * nxp + (ix + 1) # IXP=IX+1, IYP=IY, IZP=IZ
                ip3 = (iz - 1) * nxp * nyp + (iy) * nxp + (ix + 1)     # IXP=IX+1, IYP=IY+1, IZP=IZ
                ip4 = (iz) * nxp * nyp + (iy) * nxp + (ix + 1)         # IXP=IX+1, IYP=IY+1, IZP=IZ+1
                ip5 = (iz) * nxp * nyp + (iy - 1) * nxp + ix           # IXP=IX, IYP=IY, IZP=IZ+1
                ip6 = (iz - 1) * nxp * nyp + (iy - 1) * nxp + ix       # IXP=IX, IYP=IY, IZP=IZ
                ip7 = (iz - 1) * nxp * nyp + (iy) * nxp + ix           # IXP=IX, IYP=IY+1, IZP=IZ
                ip8 = (iz) * nxp * nyp + (iy) * nxp + ix               # IXP=IX, IYP=IY+1, IZP=IZ+1

                ipcell_raw[ic, 1] = ip1
                ipcell_raw[ic, 2] = ip2
                ipcell_raw[ic, 3] = ip3
                ipcell_raw[ic, 4] = ip4
                ipcell_raw[ic, 5] = ip5
                ipcell_raw[ic, 6] = ip6
                ipcell_raw[ic, 7] = ip7
                ipcell_raw[ic, 8] = ip8

                # Internal faces
                if iz != 1
                    is += 1
                    ipface_raw[is, 1:4] .= (ip2, ip6, ip7, ip3)
                    icface[is, 1] = ic
                    niscell_count[ic] += 1
                    iscell_raw[ic, niscell_count[ic]] = is
                    ic_nb = ic - nx * ny
                    icface[is, 2] = ic_nb
                    niscell_count[ic_nb] += 1
                    iscell_raw[ic_nb, niscell_count[ic_nb]] = is
                end

                if iy != 1
                    is += 1
                    ipface_raw[is, 1:4] .= (ip1, ip5, ip6, ip2)
                    icface[is, 1] = ic
                    niscell_count[ic] += 1
                    iscell_raw[ic, niscell_count[ic]] = is
                    ic_nb = ic - nx
                    icface[is, 2] = ic_nb
                    niscell_count[ic_nb] += 1
                    iscell_raw[ic_nb, niscell_count[ic_nb]] = is
                end

                if ix != 1
                    is += 1
                    ipface_raw[is, 1:4] .= (ip5, ip8, ip7, ip6)
                    icface[is, 1] = ic
                    niscell_count[ic] += 1
                    iscell_raw[ic, niscell_count[ic]] = is
                    ic_nb = ic - 1
                    icface[is, 2] = ic_nb
                    niscell_count[ic_nb] += 1
                    iscell_raw[ic_nb, niscell_count[ic_nb]] = is
                end
            end
        end
    end

    # Boundary faces

    # IX=1 boundary
    ix = 1
    for iy in 1:ny
        for iz in 1:nz
            ip5 = (iz) * nxp * nyp + (iy - 1) * nxp + ix
            ip6 = (iz - 1) * nxp * nyp + (iy - 1) * nxp + ix
            ip7 = (iz - 1) * nxp * nyp + (iy) * nxp + ix
            ip8 = (iz) * nxp * nyp + (iy) * nxp + ix
            is += 1
            ipface_raw[is, 1:4] .= (ip5, ip8, ip7, ip6)
            ic = (iz - 1) * nx * ny + (iy - 1) * nx + ix
            icface[is, 1] = ic
            niscell_count[ic] += 1
            iscell_raw[ic, niscell_count[ic]] = is
            icface[is, 2] = 0
        end
    end

    # IX=NX boundary
    ix = nx
    for iy in 1:ny
        for iz in 1:nz
            ip1 = (iz) * nxp * nyp + (iy - 1) * nxp + (ix + 1)
            ip2 = (iz - 1) * nxp * nyp + (iy - 1) * nxp + (ix + 1)
            ip3 = (iz - 1) * nxp * nyp + (iy) * nxp + (ix + 1)
            ip4 = (iz) * nxp * nyp + (iy) * nxp + (ix + 1)
            is += 1
            ipface_raw[is, 1:4] .= (ip1, ip2, ip3, ip4)
            ic = (iz - 1) * nx * ny + (iy - 1) * nx + ix
            icface[is, 1] = ic
            niscell_count[ic] += 1
            iscell_raw[ic, niscell_count[ic]] = is
            icface[is, 2] = 0
        end
    end

    # IY=1 boundary
    iy = 1
    for ix in 1:nx
        for iz in 1:nz
            ip1 = (iz) * nxp * nyp + (iy - 1) * nxp + (ix + 1)
            ip2 = (iz - 1) * nxp * nyp + (iy - 1) * nxp + (ix + 1)
            ip5 = (iz) * nxp * nyp + (iy - 1) * nxp + ix
            ip6 = (iz - 1) * nxp * nyp + (iy - 1) * nxp + ix
            is += 1
            ipface_raw[is, 1:4] .= (ip1, ip5, ip6, ip2)
            ic = (iz - 1) * nx * ny + (iy - 1) * nx + ix
            icface[is, 1] = ic
            niscell_count[ic] += 1
            iscell_raw[ic, niscell_count[ic]] = is
            icface[is, 2] = 0
        end
    end

    # IY=NY boundary
    iy = ny
    for ix in 1:nx
        for iz in 1:nz
            ip3 = (iz - 1) * nxp * nyp + (iy) * nxp + (ix + 1)
            ip4 = (iz) * nxp * nyp + (iy) * nxp + (ix + 1)
            ip7 = (iz - 1) * nxp * nyp + (iy) * nxp + ix
            ip8 = (iz) * nxp * nyp + (iy) * nxp + ix
            is += 1
            ipface_raw[is, 1:4] .= (ip3, ip7, ip8, ip4)
            ic = (iz - 1) * nx * ny + (iy - 1) * nx + ix
            icface[is, 1] = ic
            niscell_count[ic] += 1
            iscell_raw[ic, niscell_count[ic]] = is
            icface[is, 2] = 0
        end
    end

    # IZ=1 boundary
    iz = 1
    for ix in 1:nx
        for iy in 1:ny
            ip2 = (iz - 1) * nxp * nyp + (iy - 1) * nxp + (ix + 1)
            ip3 = (iz - 1) * nxp * nyp + (iy) * nxp + (ix + 1)
            ip6 = (iz - 1) * nxp * nyp + (iy - 1) * nxp + ix
            ip7 = (iz - 1) * nxp * nyp + (iy) * nxp + ix
            is += 1
            ipface_raw[is, 1:4] .= (ip2, ip6, ip7, ip3)
            ic = (iz - 1) * nx * ny + (iy - 1) * nx + ix
            icface[is, 1] = ic
            niscell_count[ic] += 1
            iscell_raw[ic, niscell_count[ic]] = is
            icface[is, 2] = 0
        end
    end

    # IZ=NZ boundary
    iz = nz
    for ix in 1:nx
        for iy in 1:ny
            ip1 = (iz) * nxp * nyp + (iy - 1) * nxp + (ix + 1)
            ip4 = (iz) * nxp * nyp + (iy) * nxp + (ix + 1)
            ip5 = (iz) * nxp * nyp + (iy - 1) * nxp + ix
            ip8 = (iz) * nxp * nyp + (iy) * nxp + ix
            is += 1
            ipface_raw[is, 1:4] .= (ip1, ip4, ip8, ip5)
            ic = (iz - 1) * nx * ny + (iy - 1) * nx + ix
            icface[is, 1] = ic
            niscell_count[ic] += 1
            iscell_raw[ic, niscell_count[ic]] = is
            icface[is, 2] = 0
        end
    end

    @assert is == nface "Expected $nface faces, got $is"

    # Convert to dynamic arrays
    ipcell = [ipcell_raw[ic, 1:nipcell[ic]] for ic in 1:ncell]
    ipface = [ipface_raw[f, 1:nipface[f]] for f in 1:nface]
    iscell = [iscell_raw[ic, 1:niscell_count[ic]] for ic in 1:ncell]

    return Grid(ncell, nface, npoint, icface, ipcell, ipface, iscell, vnode)
end

"""
    constgrid(igrid::Int, nx::Int, ny::Int, nz::Int,
              xlast::Float64, ylast::Float64, zlast::Float64) -> Grid

Construct a computational grid.

For `igrid == 1`: constructs a uniform structured cubic grid.
For `igrid == 2`: reads an OpenFOAM-format grid from files `faces`, `owner`, `neighbour`, `points`.
"""
function constgrid(igrid::Int, nx::Int, ny::Int, nz::Int,
                   xlast::Float64, ylast::Float64, zlast::Float64)
    return constgrid(Val(igrid), nx, ny, nz, xlast, ylast, zlast)
end

# Convenience constructor for uniform grids
constgrid(nx::Int, ny::Int, nz::Int,
          xlast::Float64, ylast::Float64, zlast::Float64) =
    constgrid(1, nx, ny, nz, xlast, ylast, zlast)

"""
    constgrid(::Val{2}, ncell::Int, nface_hint::Int, npoint_hint::Int;
              basedir::AbstractString=".") -> Grid

Construct a grid from OpenFOAM-format files (`faces`, `owner`, `neighbour`, `points`).

# Arguments
- `ncell`: number of cells
- `nface_hint`: expected number of faces
- `npoint_hint`: expected number of grid nodes
- `basedir`: directory containing the OpenFOAM mesh files (default: current directory)
"""
function constgrid(::Val{2}, ncell::Int, nface_hint::Int, npoint_hint::Int;
                   basedir::AbstractString=".")
    # Read faces
    nface, ipface = _read_openfoam_faces(joinpath(basedir, "faces"))

    # Read owner
    nface2, owner = _read_openfoam_intarray(joinpath(basedir, "owner"))
    @assert nface2 == nface "Number of faces in 'faces' ($nface) and 'owner' ($nface2) differ"

    # Initialize cell arrays
    niscell_count = zeros(Int, ncell)
    nipcell_count = zeros(Int, ncell)
    iscell_raw = [Int[] for _ in 1:ncell]
    ipcell_raw = [Int[] for _ in 1:ncell]
    icface = zeros(Int, nface, 2)

    # Process owner
    for iface in 1:nface
        ic = owner[iface] + 1  # convert 0-based to 1-based
        icface[iface, 1] = ic
        icface[iface, 2] = 0
        push!(iscell_raw[ic], iface)
        niscell_count[ic] += 1
        # Add unique face nodes to cell
        for ip in ipface[iface]
            if !(ip in ipcell_raw[ic])
                push!(ipcell_raw[ic], ip)
            end
        end
    end

    # Read neighbour
    nfaceint, neighbour = _read_openfoam_intarray(joinpath(basedir, "neighbour"))

    for iface in 1:nfaceint
        ic = neighbour[iface] + 1  # convert 0-based to 1-based
        icface[iface, 2] = ic
        push!(iscell_raw[ic], iface)
        niscell_count[ic] += 1
        for ip in ipface[iface]
            if !(ip in ipcell_raw[ic])
                push!(ipcell_raw[ic], ip)
            end
        end
    end

    # Read points
    npoint, vnode = _read_openfoam_points(joinpath(basedir, "points"))

    return Grid(ncell, nface, npoint, icface, ipcell_raw, ipface, iscell_raw, vnode)
end

"""
    constgrid_openfoam(ncell::Int, nface::Int, npoint::Int;
                       basedir::AbstractString=".") -> Grid

Convenience wrapper for constructing a grid from OpenFOAM files.
"""
function constgrid_openfoam(ncell::Int, nface::Int, npoint::Int;
                            basedir::AbstractString=".")
    return constgrid(Val(2), ncell, nface, npoint; basedir=basedir)
end

# -----------------------------------------------------------------------
# OpenFOAM file readers
# -----------------------------------------------------------------------

"""Parse an OpenFOAM faces file, returns (nface, ipface)."""
function _read_openfoam_faces(filename::AbstractString)
    lines = readlines(filename)
    idx = 1
    nface = parse(Int, strip(lines[idx])); idx += 1
    # Skip header line (e.g. "(")
    idx += 1

    ipface = Vector{Vector{Int}}(undef, nface)
    for iface in 1:nface
        line = strip(lines[idx]); idx += 1
        # Try to parse OpenFOAM format: N(v1 v2 ... vN)
        m = match(r"^(\d+)\s*\((.+)\)", line)
        if m !== nothing
            n = parse(Int, m.captures[1])
            verts = parse.(Int, split(strip(m.captures[2])))
            ipface[iface] = verts .+ 1  # 0-based → 1-based
        else
            # Alternative: multi-line format
            n = parse(Int, line)
            # Next line(s) may contain the vertex list
            vline = strip(lines[idx]); idx += 1
            verts = parse.(Int, split(vline))
            ipface[iface] = verts .+ 1
        end
    end

    return nface, ipface
end

"""Parse an OpenFOAM owner/neighbour file, returns (count, values)."""
function _read_openfoam_intarray(filename::AbstractString)
    lines = readlines(filename)
    idx = 1
    count = parse(Int, strip(lines[idx])); idx += 1
    # Skip header line
    idx += 1
    values = parse.(Int, split(strip(lines[idx])))
    return count, values
end

"""Parse an OpenFOAM points file, returns (npoint, vnode)."""
function _read_openfoam_points(filename::AbstractString)
    lines = readlines(filename)
    idx = 1
    npoint = parse(Int, strip(lines[idx])); idx += 1
    # Skip header line
    idx += 1

    vnode = zeros(npoint, 3)
    for ip in 1:npoint
        line = strip(lines[idx]); idx += 1
        m = match(r"\((.+)\)", line)
        if m !== nothing
            coords = parse.(Float64, split(strip(m.captures[1])))
            vnode[ip, :] .= coords
        end
    end

    return npoint, vnode
end

# -----------------------------------------------------------------------
# CELLGRID: Extract a polyhedral cell from the grid
# -----------------------------------------------------------------------

"""
    cellgrid(grid::Grid, ic::Int) -> Polyhedron

Extract the polyhedral cell `ic` from the grid as a `Polyhedron`.

This translates the Fortran `CELLGRID` subroutine. The face vertex ordering
is reversed for faces where the cell is the neighbour (not the owner) to 
ensure consistent outward normals.
"""
function cellgrid(grid::Grid, ic::Int)
    ntp = length(grid.ipcell[ic])

    # Copy node coordinates to local polyhedron vertices
    vertp = Vector{NTuple{3,Float64}}(undef, ntp)
    ipg = grid.ipcell[ic]  # global node indices for this cell
    for ip in 1:ntp
        g = ipg[ip]
        vertp[ip] = (grid.vnode[g, 1], grid.vnode[g, 2], grid.vnode[g, 3])
    end

    nts = length(grid.iscell[ic])
    ipv = Vector{Vector{Int}}(undef, nts)

    for is_idx in 1:nts
        iface = grid.iscell[ic][is_idx]
        nipv_face = length(grid.ipface[iface])

        # Determine traversal direction based on owner/neighbour
        if ic == grid.icface[iface, 1]
            # Cell is owner → normal order
            face_verts = Int[]
            for iv in 1:nipv_face
                global_ip = grid.ipface[iface][iv]
                # Find local index
                iplocal = findfirst(==(global_ip), ipg)
                push!(face_verts, iplocal)
            end
        else
            # Cell is neighbour → reverse order
            face_verts = Int[]
            for iv in nipv_face:-1:1
                global_ip = grid.ipface[iface][iv]
                iplocal = findfirst(==(global_ip), ipg)
                push!(face_verts, iplocal)
            end
        end

        ipv[is_idx] = face_verts
    end

    return Polyhedron(ipv, vertp)
end

# -----------------------------------------------------------------------
# GETFPHI: Assign scalar field values to grid nodes
# -----------------------------------------------------------------------

"""
    getfphi(grid::Grid, ifunc::Int) -> Vector{Float64}
    getfphi(grid::Grid, func::Function) -> Vector{Float64}
    getfphi(grid::Grid, filename::AbstractString) -> Vector{Float64}

Assign scalar field values to the grid nodes.

# Methods
- `ifunc == 0`: read from file `"phigrid"`
- `ifunc == 1`: sphere function
- `ifunc == 2`: torus function
- `ifunc == 3`: orthocircle function
- `func::Function`: custom function `f(x, y, z) -> Float64`
- `filename::AbstractString`: read from specified file
"""
function getfphi(grid::Grid, ifunc::Int; phifile::AbstractString="phigrid")
    fphi = zeros(grid.npoint)
    if ifunc == 0
        fphi = _read_phi_file(phifile, grid.npoint)
    elseif ifunc == 1
        for ip in 1:grid.npoint
            fphi[ip] = func_sphere(grid.vnode[ip, 1], grid.vnode[ip, 2], grid.vnode[ip, 3])
        end
    elseif ifunc == 2
        for ip in 1:grid.npoint
            fphi[ip] = func_torus(grid.vnode[ip, 1], grid.vnode[ip, 2], grid.vnode[ip, 3])
        end
    elseif ifunc == 3
        for ip in 1:grid.npoint
            fphi[ip] = func_orthocircle(grid.vnode[ip, 1], grid.vnode[ip, 2], grid.vnode[ip, 3])
        end
    else
        error("Invalid ifunc=$ifunc. Use 0 (file), 1 (sphere), 2 (torus), or 3 (orthocircle).")
    end
    return fphi
end

function getfphi(grid::Grid, func::Function)
    fphi = Vector{Float64}(undef, grid.npoint)
    for ip in 1:grid.npoint
        fphi[ip] = func(grid.vnode[ip, 1], grid.vnode[ip, 2], grid.vnode[ip, 3])
    end
    return fphi
end

function getfphi(grid::Grid, filename::AbstractString)
    return _read_phi_file(filename, grid.npoint)
end

function _read_phi_file(filename::AbstractString, npoint::Int)
    fphi = zeros(npoint)
    open(filename, "r") do io
        for ip in 1:npoint
            line = readline(io)
            try
                fphi[ip] = parse(Float64, strip(line))
            catch
                fphi[ip] = 0.0
            end
        end
    end
    return fphi
end

# -----------------------------------------------------------------------
# VTKGRID: Print the computational grid in VTK format
# -----------------------------------------------------------------------

"""
    vtkgrid(filename::AbstractString, grid::Grid)
    vtkgrid(grid::Grid)

Write the computational grid faces in VTK legacy polydata format.

If no filename is given, writes to `"grid.vtk"`.
"""
function vtkgrid(filename::AbstractString, grid::Grid)
    open(filename, "w") do io
        println(io, "# vtk DataFile Version 2.0")
        println(io, "Grid")
        println(io, "ASCII")
        println(io)
        println(io, "DATASET POLYDATA")
        @printf(io, "POINTS%9d float\n", grid.npoint)
        for ip in 1:grid.npoint
            @printf(io, "%12.6f%12.6f%12.6f\n",
                    grid.vnode[ip, 1], grid.vnode[ip, 2], grid.vnode[ip, 3])
        end
        ndata = 0
        for iface in 1:grid.nface
            ndata += 1 + length(grid.ipface[iface])
        end
        @printf(io, "POLYGONS%9d%9d\n", grid.nface, ndata)
        for iface in 1:grid.nface
            nipf = length(grid.ipface[iface])
            @printf(io, "%7d\n", nipf)
            for ip in 1:nipf
                @printf(io, "%7d\n", grid.ipface[iface][ip] - 1)  # 0-based for VTK
            end
        end
    end
    return nothing
end

vtkgrid(grid::Grid) = vtkgrid("grid.vtk", grid)

# -----------------------------------------------------------------------
# ISOVTKGRID: Extract isosurface on grid and write VTK output
# -----------------------------------------------------------------------

"""
    isovtkgrid(grid::Grid, fphi::Vector{Float64}, phiiso::Float64;
               ifile::Int=0, filename::Union{Nothing,AbstractString}=nothing) -> Nothing

Extract the isosurface from the entire grid and write the result as a VTK file.

This is the grid-level equivalent of running `isoap` on every cell that is
cut by the isosurface, then combining the results into a single VTK polydata file.

# Arguments
- `grid`: the computational grid
- `fphi`: scalar field values at each grid node
- `phiiso`: iso-value
- `ifile`: file number for output naming (default: 0 → `isogrid-0000.vtk`)
- `filename`: explicit output filename (overrides `ifile` if given)
"""
function isovtkgrid(grid::Grid, fphi::Vector{Float64}, phiiso::Float64;
                    ifile::Int=0, filename::Union{Nothing,AbstractString}=nothing)
    # Mark cells that are cut by the isosurface
    ncelliso = 0
    iciso = Int[]
    for ic in 1:grid.ncell
        phimin = Inf
        phimax = -Inf
        for iv in 1:length(grid.ipcell[ic])
            ip = grid.ipcell[ic][iv]
            phimin = min(phimin, fphi[ip])
            phimax = max(phimax, fphi[ip])
        end
        if phimin < phiiso && phimax > phiiso
            ncelliso += 1
            push!(iciso, ic)
        end
    end

    # Process each cut cell
    all_polygons = Vector{Vector{NTuple{3,Float64}}}()  # polygon vertices
    for iiso in 1:ncelliso
        ic = iciso[iiso]

        # Get phi values for cell vertices
        phi_cell = Float64[]
        for iv in 1:length(grid.ipcell[ic])
            ip = grid.ipcell[ic][iv]
            push!(phi_cell, fphi[ip])
        end

        # Extract polyhedral cell from grid
        poly = cellgrid(grid, ic)

        # Run isoap
        result = isoap(poly, phi_cell, phiiso)

        # Collect iso-polygons with absolute vertex coordinates
        for iso_idx in 1:niso(result)
            verts = NTuple{3,Float64}[]
            for ip in result.ipviso[iso_idx]
                push!(verts, result.vertiso[ip])
            end
            push!(all_polygons, verts)
        end
    end

    # Write combined VTK file
    if filename === nothing
        filename = @sprintf("isogrid-%04d.vtk", ifile)
    end

    # Count total vertices and build output arrays
    nvlc = 0
    for poly in all_polygons
        nvlc += length(poly)
    end

    npol = length(all_polygons)
    ndata = 0
    for poly in all_polygons
        ndata += 1 + length(poly)
    end

    open(filename, "w") do io
        println(io, "# vtk DataFile Version 2.0")
        println(io, "LC")
        println(io, "ASCII")
        println(io)
        println(io, "DATASET POLYDATA")
        @printf(io, "POINTS%10d float\n", nvlc)

        for poly in all_polygons
            for v in poly
                @printf(io, "%12.6f%12.6f%12.6f\n", v[1], v[2], v[3])
            end
        end

        @printf(io, "POLYGONS%10d%10d\n", npol, ndata)
        offset = 0
        for poly in all_polygons
            @printf(io, "%7d\n", length(poly))
            for iv in 1:length(poly)
                @printf(io, "%7d\n", offset + iv - 1)  # 0-based
            end
            offset += length(poly)
        end
    end

    return nothing
end

# -----------------------------------------------------------------------
# ISOAPVARDEF: Read single-cell test case parameters
# -----------------------------------------------------------------------

"""
    isoapvardef(configfile::AbstractString="isoapvardef") -> (icellgeom, ifunc, phiiso)

Read the configuration file for single-cell test cases.

Returns:
- `icellgeom`: cell geometry index
- `ifunc`: implicit function index (0 = read from file)
- `phiiso`: iso-value
"""
function isoapvardef(configfile::AbstractString="isoapvardef")
    lines = readlines(configfile)
    icellgeom = parse(Int, strip(lines[2]))
    ifunc = parse(Int, strip(lines[4]))
    phiiso = parse(Float64, strip(lines[6]))
    return (icellgeom, ifunc, phiiso)
end

"""
    get_cell_geometry(icellgeom::Int) -> Polyhedron

Return the polyhedron corresponding to the given `icellgeom` code.

Geometry codes:
- 1: cube, 2: tetrahedron, 3: dodecahedron, 4: icosahedron, 5: complex cell
- 101: distorted cube, 102: pentagonal pyramid, 103: cut cube, 104: stellated cube
- 105: non-convex hexahedron, 106: stellated dodecahedron, 107: stellated icosahedron
- 108: hollowed cube, 109: drilled cube, 110: zig-zag prism
"""
function get_cell_geometry(icellgeom::Int)
    geom_map = Dict(
        1 => cube, 2 => tetrahedron, 3 => dodecahedron, 4 => icosahedron,
        5 => complexcell, 101 => distortedcube, 102 => pentapyramid,
        103 => cutcube, 104 => stellatedcube, 105 => nchexahedron,
        106 => stellateddodecahedron, 107 => stellatedicosahedron,
        108 => hollowedcube, 109 => drilledcube, 110 => zigzagcell,
    )
    if !haskey(geom_map, icellgeom)
        error("Invalid ICELLGEOM=$icellgeom. Valid values: 1-5, 101-110.")
    end
    return geom_map[icellgeom]()
end

"""
    get_phi_function(ifunc::Int) -> Function

Return the implicit function corresponding to the given `ifunc` code.

- 1: sphere, 2: torus, 3: orthocircle
"""
function get_phi_function(ifunc::Int)
    func_map = Dict(1 => func_sphere, 2 => func_torus, 3 => func_orthocircle)
    if !haskey(func_map, ifunc)
        error("Invalid IFUNC=$ifunc. Valid values: 1, 2, 3.")
    end
    return func_map[ifunc]
end

"""
    assign_phi(poly::Polyhedron, ifunc::Int, phiiso::Float64;
               phifile::AbstractString="phi") -> Vector{Float64}

Assign scalar phi values to each vertex of the polyhedron.

- `ifunc == 0`: read from file
- `ifunc == 1,2,3`: evaluate implicit function
"""
function assign_phi(poly::Polyhedron, ifunc::Int;
                    phifile::AbstractString="phi")
    ntp = nvertices(poly)
    phi = zeros(ntp)
    if ifunc == 0
        open(phifile, "r") do io
            for ip in 1:ntp
                line = readline(io)
                try
                    phi[ip] = parse(Float64, strip(line))
                catch
                    phi[ip] = 0.0
                end
            end
        end
    else
        func = get_phi_function(ifunc)
        for ip in 1:ntp
            x, y, z = poly.vertp[ip]
            phi[ip] = func(x, y, z)
        end
    end
    return phi
end

# -----------------------------------------------------------------------
# ISOAPVARDEFGRID: Read grid test case parameters
# -----------------------------------------------------------------------

"""
    isoapvardefgrid(configfile::AbstractString="isoapvardefgrid")

Read the configuration file for grid test cases.

Returns a named tuple with fields:
- `ifunc`: test case number
- `igrid`: grid type (1=uniform, 2=OpenFOAM)
- `nx, ny, nz`: cells per axis (for igrid==1)
- `phiiso`: iso-value
- `xlast, ylast, zlast`: domain extents (for igrid==1)
- `ncell, nface, npoint`: mesh sizes
"""
function isoapvardefgrid(configfile::AbstractString="isoapvardefgrid")
    lines = readlines(configfile)
    ifunc = parse(Int, strip(lines[2]))
    igrid = parse(Int, strip(lines[4]))
    nx = parse(Int, strip(lines[6]))
    ny = parse(Int, strip(lines[8]))
    nz = parse(Int, strip(lines[10]))
    phiiso = parse(Float64, strip(lines[12]))

    if igrid == 1
        xlast = parse(Float64, strip(lines[14]))
        ylast = parse(Float64, strip(lines[16]))
        zlast = parse(Float64, strip(lines[18]))
        ncell = nx * ny * nz
        npoint = (nx + 1) * (ny + 1) * (nz + 1)
        nface = 3 * ncell + nx * ny + nx * nz + ny * nz
    else
        # Read meshdata for OpenFOAM
        xlast = 0.0
        ylast = 0.0
        zlast = 0.0
        # Would need to read meshdata file
        ncell = 0
        nface = 0
        npoint = 0
    end

    return (ifunc=ifunc, igrid=igrid, nx=nx, ny=ny, nz=nz, phiiso=phiiso,
            xlast=xlast, ylast=ylast, zlast=zlast, ncell=ncell, nface=nface,
            npoint=npoint)
end

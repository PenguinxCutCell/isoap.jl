"""
    Polyhedron

Representation of a general polyhedron for isosurface extraction.

# Fields
- `ipv::Vector{Vector{Int}}`: vertex indices for each face (1-based)
- `vertp::Vector{NTuple{3,Float64}}`: coordinates of the polyhedron vertices
"""
struct Polyhedron
    ipv::Vector{Vector{Int}}       # ipv[face] = [vertex indices...]
    vertp::Vector{NTuple{3,Float64}}  # vertp[vertex] = (x, y, z)
end

"""
    nfaces(p::Polyhedron)

Number of faces of the polyhedron.
"""
nfaces(p::Polyhedron) = length(p.ipv)

"""
    nvertices(p::Polyhedron)

Number of vertices of the polyhedron.
"""
nvertices(p::Polyhedron) = length(p.vertp)

"""
    IsoResult

Result of isosurface extraction from a polyhedron.

# Fields
- `ipviso::Vector{Vector{Int}}`: vertex indices for each iso-polygon (1-based)
- `isoeface::Vector{Int}`: face index of the polyhedron over which each iso-edge is constructed
- `vertiso::Vector{NTuple{3,Float64}}`: coordinates of the iso-vertices
"""
struct IsoResult
    ipviso::Vector{Vector{Int}}       # ipviso[polygon] = [iso-vertex indices...]
    isoeface::Vector{Int}             # isoeface[iso-vertex] = face index
    vertiso::Vector{NTuple{3,Float64}}  # vertiso[iso-vertex] = (x, y, z)
end

"""
    niso(r::IsoResult)

Number of iso-polygons in the result.
"""
niso(r::IsoResult) = length(r.ipviso)

"""
    ISOAP

A Julia package for isosurface extraction on arbitrary polyhedra.

Translated from the Fortran code by J. Lopez and J. Hernandez (2020).

Reference:
    J. Lopez, J. Hernandez, P. Gomez, R. Zamora, C. Zanzi, F. Faura,
    "A new isosurface extraction method for VOF interface reconstruction
    on arbitrary grids", Computer Methods for Applied Mechanics and Engineering.
"""
module ISOAP

# Types
export Polyhedron, IsoResult, Grid

# Core algorithm
export isoap, isopol

# Geometries
export cube, tetrahedron, dodecahedron, icosahedron, complexcell
export distortedcube, pentapyramid, cutcube, stellatedcube
export nchexahedron, stellateddodecahedron, stellatedicosahedron
export hollowedcube, drilledcube, zigzagcell

# Implicit functions
export func_sphere, func_torus, func_orthocircle

# VTK output (single cell)
export polvtk, isovtk

# Grid infrastructure
export constgrid, constgrid_openfoam, cellgrid, getfphi, vtkgrid, isovtkgrid

# Configuration readers
export isoapvardef, isoapvardefgrid, get_cell_geometry, get_phi_function, assign_phi

include("types.jl")
include("core.jl")
include("geometries.jl")
include("vtk.jl")
include("grid.jl")

end # module ISOAP

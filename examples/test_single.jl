# -----------------------------------------------------------------------
# test_single.jl
# -----------------------------------------------------------------------
# Julia equivalent of PROGRAM TEST (test.f / test.c)
#
# Test program for iso-surface extraction on single cells using ISOAP.
#
# Usage:
#   julia --project=. examples/test_single.jl [icellgeom] [ifunc] [phiiso]
#
# Default: cube (1) with sphere function (1) at phiiso=0.0
# -----------------------------------------------------------------------

using ISOAP

function main()
    println("-----------------------------------------------------")
    println("|---------------------------------------------------|")
    println("|      TEST PROGRAM FOR SINGLE CELLS OF ISOAP       |")
    println("|                                                   |")
    println("|             Julia version                         |")
    println("|                                                   |")
    println("|           Original Fortran code by                |")
    println("|            J. Lopez and J. Hernandez              |")
    println("|---------------------------------------------------|")
    println("-----------------------------------------------------")
    println()

    # Parse command-line arguments or use defaults
    icellgeom = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
    ifunc     = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
    phiiso    = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.0

    println("ICELLGEOM: ", icellgeom)
    println("IFUNC: ", ifunc)
    println("PHIISO: ", phiiso)
    println()

    # Cell geometry selection
    # ICELLGEOM:
    #   1 = cube, 2 = tetrahedron, 3 = dodecahedron, 4 = icosahedron,
    #   5 = complex cell (18 faces, 32 vertices)
    # 101 = distorted cube, 102 = pentagonal pyramid, 103 = cut cube,
    # 104 = stellated cube, 105 = non-convex hexahedron,
    # 106 = stellated dodecahedron, 107 = stellated icosahedron,
    # 108 = hollowed cube, 109 = drilled cube, 110 = zig-zag prism
    poly = get_cell_geometry(icellgeom)
    println("Cell: $(ISOAP.nfaces(poly)) faces, $(ISOAP.nvertices(poly)) vertices")

    # Assign scalar phi using implicit function
    # IFUNC:
    #   0 = read from file "phi"
    #   1 = sphere with radius 0.325 centered at (0.5,0.5,0.5)
    #   2 = torus with major radius 0.2, minor radius 0.1, centered at (0.5,0.5,0.5)
    #   3 = orthocircle surface centered at (1.25,1.25,1.25)
    phi = assign_phi(poly, ifunc)

    # Print the polyhedral cell to 'geo00000.vtk'
    ifile = 0
    polvtk(ifile, poly)
    println("Wrote geo00000.vtk")

    # Iso-surface extraction
    result = isoap(poly, phi, phiiso)
    println("Iso-surface: $(ISOAP.niso(result)) polygon(s), $(length(result.vertiso)) vertices")

    # Print the iso-polygons to 'iso00000.vtk'
    isovtk(ifile, result)
    println("Wrote iso00000.vtk")
end

main()

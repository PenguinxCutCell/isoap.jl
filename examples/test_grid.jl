# -----------------------------------------------------------------------
# test_grid.jl
# -----------------------------------------------------------------------
# Julia equivalent of PROGRAM TESTGRID (testgrid.f / testgrid.c)
#
# Test program for iso-surface extraction on computational grids
# using ISOAP.
#
# Usage:
#   julia --project=. examples/test_grid.jl [ifunc] [nx] [ny] [nz] [phiiso]
#
# Default: sphere function (1) on 10×10×10 grid with phiiso=0.0
# -----------------------------------------------------------------------

using ISOAP

function main()
    println("-----------------------------------------------------")
    println("|---------------------------------------------------|")
    println("|         TEST PROGRAM FOR GRIDS OF ISOAP           |")
    println("|                                                   |")
    println("|             Julia version                         |")
    println("|                                                   |")
    println("|           Original Fortran code by                |")
    println("|            J. Lopez and J. Hernandez              |")
    println("|---------------------------------------------------|")
    println("-----------------------------------------------------")
    println()

    # Parse command-line arguments or use defaults
    ifunc  = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
    nx     = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 10
    ny     = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 10
    nz     = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 10
    phiiso = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 0.0
    xlast  = length(ARGS) >= 6 ? parse(Float64, ARGS[6]) : 1.0
    ylast  = length(ARGS) >= 7 ? parse(Float64, ARGS[7]) : 1.0
    zlast  = length(ARGS) >= 8 ? parse(Float64, ARGS[8]) : 1.0

    # Iso-surface extraction test cases:
    # IFUNC = 1: sphere with radius 0.325 centered at (0.5,0.5,0.5)
    #       = 2: torus with major radius 0.2, minor radius 0.1,
    #            centered at (0.5,0.5,0.5)
    #       = 3: orthocircle surface centered at (1.25,1.25,1.25)
    println("IFUNC: ", ifunc)
    println("Grid: $(nx) × $(ny) × $(nz)")
    println("Domain: [0, $xlast] × [0, $ylast] × [0, $zlast]")
    println("PHIISO: ", phiiso)
    println()

    # Construct the grid
    tstart = time()
    grid = constgrid(1, nx, ny, nz, xlast, ylast, zlast)
    tgrid = time() - tstart
    println("Grid constructed: $(grid.ncell) cells, $(grid.nface) faces, $(grid.npoint) nodes")

    # Uncomment to print the grid in VTK format:
    # vtkgrid(grid)
    # println("Wrote grid.vtk")

    # Assign scalar field to the grid points
    tstart = time()
    fphi = getfphi(grid, ifunc)
    tphi = time() - tstart

    # Extract the iso-surface
    tstart = time()
    isovtkgrid(grid, fphi, phiiso; ifile=0)
    tiso = time() - tstart
    println("Wrote isogrid-0000.vtk")

    println()
    println("-------------------------------------------------")
    println("|              EXECUTION TIMES:                 |")
    println("-------------------------------------------------")
    println("GRID CONSTRUCTION       = $(round(tgrid; digits=6)) secs")
    println("SCALAR PHI ASSIGNMENT   = $(round(tphi; digits=6)) secs")
    println("ISO-SURFACE EXTRACTION  = $(round(tiso; digits=6)) secs")
    println("-------------------------------------------------")
end

main()

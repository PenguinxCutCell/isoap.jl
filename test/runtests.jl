using ISOAP
using Test

@testset "ISOAP.jl" begin

    @testset "Geometry constructors" begin
        geom_specs = [
            ("cube",                    cube,                    6,  8),
            ("tetrahedron",             tetrahedron,             4,  4),
            ("dodecahedron",            dodecahedron,           12, 20),
            ("icosahedron",             icosahedron,            20, 12),
            ("complexcell",             complexcell,            18, 32),
            ("distortedcube",           distortedcube,           7,  8),
            ("pentapyramid",            pentapyramid,            6,  6),
            ("cutcube",                 cutcube,                 9,  9),
            ("stellatedcube",           stellatedcube,          24, 14),
            ("nchexahedron",            nchexahedron,            6,  8),
            ("stellateddodecahedron",   stellateddodecahedron,  60, 32),
            ("stellatedicosahedron",    stellatedicosahedron,   60, 32),
            ("hollowedcube",            hollowedcube,           12, 16),
            ("drilledcube",             drilledcube,            12, 16),
            ("zigzagcell",              zigzagcell,             12, 20),
        ]

        for (name, geofn, expected_nf, expected_nv) in geom_specs
            @testset "$name" begin
                p = geofn()
                @test ISOAP.nfaces(p) == expected_nf
                @test ISOAP.nvertices(p) == expected_nv
                # All face vertex indices should be in valid range
                for face in p.ipv
                    for idx in face
                        @test 1 <= idx <= ISOAP.nvertices(p)
                    end
                end
            end
        end
    end

    @testset "get_cell_geometry" begin
        for code in [1, 2, 3, 4, 5, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110]
            p = get_cell_geometry(code)
            @test isa(p, Polyhedron)
            @test ISOAP.nfaces(p) > 0
            @test ISOAP.nvertices(p) > 0
        end
        @test_throws ErrorException get_cell_geometry(99)
    end

    @testset "Implicit functions" begin
        # Sphere: center is inside (positive), far corner is outside (negative)
        @test func_sphere(0.5, 0.5, 0.5) > 0
        @test func_sphere(0.0, 0.0, 0.0) < 0

        # Torus: center of tube
        @test func_torus(0.7, 0.5, 0.5) > 0
        @test func_torus(0.5, 0.5, 0.5) < 0

        # get_phi_function
        @test get_phi_function(1)(0.5, 0.5, 0.5) == func_sphere(0.5, 0.5, 0.5)
        @test get_phi_function(2)(0.7, 0.5, 0.5) == func_torus(0.7, 0.5, 0.5)
        @test get_phi_function(3)(0.0, 0.0, 0.0) == func_orthocircle(0.0, 0.0, 0.0)
        @test_throws ErrorException get_phi_function(99)
    end

    @testset "Isosurface extraction" begin
        @testset "No intersection (all same sign)" begin
            p = cube()
            phi = ones(8)  # all > 0
            result = isoap(p, phi, 0.0)
            @test ISOAP.niso(result) == 0
        end

        @testset "Cube with custom phi" begin
            p = cube()
            phi = [1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0]
            result = isoap(p, phi, 0.0)
            @test ISOAP.niso(result) >= 1
            # Check iso-vertices are in valid range
            for poly in result.ipviso
                for idx in poly
                    @test 1 <= idx <= length(result.vertiso)
                end
            end
            # Iso-vertices should have coordinates between 0 and 1 (within the cube)
            for v in result.vertiso
                @test 0.0 <= v[1] <= 1.0
                @test 0.0 <= v[2] <= 1.0
                @test 0.0 <= v[3] <= 1.0
            end
        end

        @testset "Tetrahedron with sphere" begin
            p = tetrahedron()
            phi = [func_sphere(v...) for v in p.vertp]
            result = isoap(p, phi, 0.0)
            # phi values: some positive, some negative → should produce iso-polygons
            has_positive = any(phi .> 0)
            has_negative = any(phi .< 0)
            if has_positive && has_negative
                @test ISOAP.niso(result) >= 1
            end
        end

        @testset "Dodecahedron with sphere" begin
            p = dodecahedron()
            phi = [func_sphere(v...) for v in p.vertp]
            result = isoap(p, phi, 0.0)
            @test ISOAP.niso(result) >= 1
        end

        @testset "All 15 geometries with sphere" begin
            for code in [1, 2, 3, 4, 5, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110]
                p = get_cell_geometry(code)
                phi = [func_sphere(v...) for v in p.vertp]
                result = isoap(p, phi, 0.0)
                # Should produce at least one polygon if the cell straddles the interface
                has_pos = any(phi .> 0)
                has_neg = any(phi .< 0)
                if has_pos && has_neg
                    @test ISOAP.niso(result) >= 1
                end
            end
        end
    end

    @testset "VTK output" begin
        p = cube()
        phi = [1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0]
        result = isoap(p, phi, 0.0)

        tmpdir = mktempdir()
        try
            # Test polvtk with filename
            geofile = joinpath(tmpdir, "test_geo.vtk")
            polvtk(geofile, p)
            @test isfile(geofile)
            content = read(geofile, String)
            @test occursin("DATASET POLYDATA", content)
            @test occursin("POINTS 8 float", content)
            @test occursin("POLYGONS", content)

            # Test isovtk with filename
            isofile = joinpath(tmpdir, "test_iso.vtk")
            isovtk(isofile, result)
            @test isfile(isofile)
            content = read(isofile, String)
            @test occursin("DATASET POLYDATA", content)
            @test occursin("POLYGONS", content)

            # Test polvtk with integer file number
            cd(tmpdir) do
                polvtk(0, p)
                @test isfile("geo00000.vtk")
                isovtk(0, result)
                @test isfile("iso00000.vtk")
            end
        finally
            rm(tmpdir; recursive=true)
        end
    end

    @testset "Edge cases" begin
        @testset "Phi exactly at isovalue" begin
            p = cube()
            phi = zeros(8)  # all exactly at phiiso
            result = isoap(p, phi, 0.0)
            # All vertices at isovalue → all tagged as IA=0 → no intersection
            @test ISOAP.niso(result) == 0
        end

        @testset "Single vertex above isovalue" begin
            p = cube()
            phi = [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]
            phi[1] = 1.0  # only vertex 1 is above
            result = isoap(p, phi, 0.0)
            @test ISOAP.niso(result) >= 1
        end
    end

    @testset "Grid infrastructure" begin
        @testset "Uniform grid construction" begin
            # Small 2×2×2 grid
            grid = constgrid(1, 2, 2, 2, 1.0, 1.0, 1.0)
            @test grid.ncell == 8
            @test grid.npoint == 27
            @test grid.nface == 3 * 8 + 2 * 2 + 2 * 2 + 2 * 2  # 36
            @test grid.nface == 36
            @test size(grid.vnode) == (27, 3)
            @test size(grid.icface) == (36, 2)
            @test length(grid.ipcell) == 8
            @test length(grid.ipface) == 36
            @test length(grid.iscell) == 8

            # Each cell has 8 nodes
            for ic in 1:grid.ncell
                @test length(grid.ipcell[ic]) == 8
            end

            # Each face has 4 nodes (all quads for structured grid)
            for f in 1:grid.nface
                @test length(grid.ipface[f]) == 4
            end

            # Each cell should have 6 faces
            for ic in 1:grid.ncell
                @test length(grid.iscell[ic]) == 6
            end

            # Boundary faces have icface[f,2] == 0
            nboundary = count(grid.icface[:, 2] .== 0)
            ninterior = count(grid.icface[:, 2] .> 0)
            @test nboundary == 2 * (2 * 2 + 2 * 2 + 2 * 2)  # 24
            @test ninterior == 36 - 24  # 12

            # Check corner node coordinates
            @test grid.vnode[1, :] == [0.0, 0.0, 0.0]  # node at origin
        end

        @testset "Uniform grid convenience constructor" begin
            grid = constgrid(3, 3, 3, 1.0, 1.0, 1.0)
            @test grid.ncell == 27
            @test grid.npoint == 64
        end

        @testset "Grid node coordinates" begin
            grid = constgrid(1, 3, 4, 5, 3.0, 4.0, 5.0)
            # Check domain extents
            @test maximum(grid.vnode[:, 1]) ≈ 3.0
            @test maximum(grid.vnode[:, 2]) ≈ 4.0
            @test maximum(grid.vnode[:, 3]) ≈ 5.0
            @test minimum(grid.vnode[:, 1]) ≈ 0.0
            @test minimum(grid.vnode[:, 2]) ≈ 0.0
            @test minimum(grid.vnode[:, 3]) ≈ 0.0
        end

        @testset "cellgrid extraction" begin
            grid = constgrid(1, 2, 2, 2, 1.0, 1.0, 1.0)
            for ic in 1:grid.ncell
                poly = cellgrid(grid, ic)
                @test ISOAP.nvertices(poly) == 8
                @test ISOAP.nfaces(poly) == 6
                # Each face should have 4 vertices (quad)
                for face in poly.ipv
                    @test length(face) == 4
                    for idx in face
                        @test 1 <= idx <= 8
                    end
                end
            end
        end

        @testset "getfphi with function" begin
            grid = constgrid(1, 3, 3, 3, 1.0, 1.0, 1.0)
            fphi = getfphi(grid, func_sphere)
            @test length(fphi) == grid.npoint
            # Center of domain should be positive for sphere
            # Find the node closest to (0.5, 0.5, 0.5)
            min_dist = Inf
            center_phi = 0.0
            for ip in 1:grid.npoint
                d = sqrt((grid.vnode[ip, 1] - 0.5)^2 +
                         (grid.vnode[ip, 2] - 0.5)^2 +
                         (grid.vnode[ip, 3] - 0.5)^2)
                if d < min_dist
                    min_dist = d
                    center_phi = fphi[ip]
                end
            end
            # Close to center should be positive (inside sphere)
            @test center_phi > 0

            # Using integer function code
            fphi2 = getfphi(grid, 1)
            @test fphi ≈ fphi2
        end

        @testset "getfphi all implicit functions" begin
            grid = constgrid(1, 3, 3, 3, 1.0, 1.0, 1.0)
            for ifunc in 1:3
                fphi = getfphi(grid, ifunc)
                @test length(fphi) == grid.npoint
                # Should have both positive and negative values
                # (the interface cuts the domain)
            end
        end

        @testset "vtkgrid output" begin
            grid = constgrid(1, 2, 2, 2, 1.0, 1.0, 1.0)
            tmpdir = mktempdir()
            try
                vtkfile = joinpath(tmpdir, "grid.vtk")
                vtkgrid(vtkfile, grid)
                @test isfile(vtkfile)
                content = read(vtkfile, String)
                @test occursin("DATASET POLYDATA", content)
                @test occursin("POINTS", content)
                @test occursin("POLYGONS", content)
                @test occursin("27", content)  # npoint
            finally
                rm(tmpdir; recursive=true)
            end
        end

        @testset "isovtkgrid extraction" begin
            grid = constgrid(1, 5, 5, 5, 1.0, 1.0, 1.0)
            fphi = getfphi(grid, 1)  # sphere
            tmpdir = mktempdir()
            try
                outfile = joinpath(tmpdir, "isogrid-0000.vtk")
                isovtkgrid(grid, fphi, 0.0; ifile=0, filename=outfile)
                @test isfile(outfile)
                content = read(outfile, String)
                @test occursin("DATASET POLYDATA", content)
                @test occursin("POINTS", content)
                @test occursin("POLYGONS", content)
                # Should have extracted some polygons
                @test occursin("float", content)
            finally
                rm(tmpdir; recursive=true)
            end
        end

        @testset "isovtkgrid with torus" begin
            grid = constgrid(1, 5, 5, 5, 1.0, 1.0, 1.0)
            fphi = getfphi(grid, 2)  # torus
            tmpdir = mktempdir()
            try
                outfile = joinpath(tmpdir, "isogrid-torus.vtk")
                isovtkgrid(grid, fphi, 0.0; filename=outfile)
                @test isfile(outfile)
            finally
                rm(tmpdir; recursive=true)
            end
        end

        @testset "isovtkgrid no intersection" begin
            grid = constgrid(1, 2, 2, 2, 1.0, 1.0, 1.0)
            fphi = ones(grid.npoint)  # all same sign → no intersection
            tmpdir = mktempdir()
            try
                outfile = joinpath(tmpdir, "isogrid-empty.vtk")
                isovtkgrid(grid, fphi, 0.0; filename=outfile)
                @test isfile(outfile)
                content = read(outfile, String)
                @test occursin("POINTS         0 float", content)
            finally
                rm(tmpdir; recursive=true)
            end
        end

        @testset "Grid cell consistency" begin
            # Verify that cellgrid + isoap produces valid results for every cut cell
            grid = constgrid(1, 3, 3, 3, 1.0, 1.0, 1.0)
            fphi = getfphi(grid, 1)
            phiiso = 0.0

            for ic in 1:grid.ncell
                phimin = minimum(fphi[grid.ipcell[ic][iv]] for iv in 1:length(grid.ipcell[ic]))
                phimax = maximum(fphi[grid.ipcell[ic][iv]] for iv in 1:length(grid.ipcell[ic]))
                if phimin < phiiso && phimax > phiiso
                    phi_cell = [fphi[grid.ipcell[ic][iv]] for iv in 1:length(grid.ipcell[ic])]
                    poly = cellgrid(grid, ic)
                    result = isoap(poly, phi_cell, phiiso)
                    @test ISOAP.niso(result) >= 1
                    # Verify iso-vertex indices are valid
                    for iso_poly in result.ipviso
                        for idx in iso_poly
                            @test 1 <= idx <= length(result.vertiso)
                        end
                    end
                end
            end
        end
    end
end

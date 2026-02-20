#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using ISOAP
using Profile

const AVAILABLE_CASES = (
    :isoap_cube,
    :isoap_tetrahedron,
    :isoap_dodecahedron,
    :isoap_icosahedron,
    :isoap_complexcell,
    :isoap_distortedcube,
    :isoap_pentapyramid,
    :isoap_cutcube,
    :isoap_stellatedcube,
    :isoap_nchexahedron,
    :isoap_stellateddodecahedron,
    :isoap_stellatedicosahedron,
    :isoap_hollowedcube,
    :isoap_drilledcube,
    :isoap_zigzagcell,
    :isopol_cube,
)

function build_case(case::Symbol)
    if case == :isoap_cube
        poly = cube()
        phi = assign_phi(poly, 1)
        phiiso = 0.0
        ws = IsoapWorkspace(poly)
        return () -> isoap!(ws, poly, phi, phiiso)
    elseif case == :isoap_tetrahedron
        poly = tetrahedron()
        phi = assign_phi(poly, 1)
        phiiso = 0.0
        ws = IsoapWorkspace(poly)
        return () -> isoap!(ws, poly, phi, phiiso)
    elseif case == :isoap_dodecahedron
        poly = dodecahedron()
        phi = assign_phi(poly, 1)
        phiiso = 0.0
        ws = IsoapWorkspace(poly)
        return () -> isoap!(ws, poly, phi, phiiso)
    elseif case == :isoap_icosahedron
        poly = icosahedron()
        phi = assign_phi(poly, 1)
        phiiso = 0.0
        ws = IsoapWorkspace(poly)
        return () -> isoap!(ws, poly, phi, phiiso)
    elseif case == :isoap_complexcell
        poly = complexcell()
        phi = assign_phi(poly, 1)
        phiiso = 0.0
        ws = IsoapWorkspace(poly)
        return () -> isoap!(ws, poly, phi, phiiso)
    elseif case == :isoap_distortedcube
        poly = distortedcube()
        phi = assign_phi(poly, 1)
        phiiso = 0.0
        ws = IsoapWorkspace(poly)
        return () -> isoap!(ws, poly, phi, phiiso)
    elseif case == :isoap_pentapyramid
        poly = pentapyramid()
        phi = assign_phi(poly, 1)
        phiiso = 0.0
        ws = IsoapWorkspace(poly)
        return () -> isoap!(ws, poly, phi, phiiso)
    elseif case == :isoap_cutcube
        poly = cutcube()
        phi = assign_phi(poly, 1)
        phiiso = 0.0
        ws = IsoapWorkspace(poly)
        return () -> isoap!(ws, poly, phi, phiiso)
    elseif case == :isoap_stellatedcube
        poly = stellatedcube()
        phi = assign_phi(poly, 1)
        phiiso = 0.0
        ws = IsoapWorkspace(poly)
        return () -> isoap!(ws, poly, phi, phiiso)
    elseif case == :isoap_nchexahedron
        poly = nchexahedron()
        phi = assign_phi(poly, 1)
        phiiso = 0.0
        ws = IsoapWorkspace(poly)
        return () -> isoap!(ws, poly, phi, phiiso)
    elseif case == :isoap_stellateddodecahedron
        poly = stellateddodecahedron()
        phi = assign_phi(poly, 1)
        phiiso = 0.0
        ws = IsoapWorkspace(poly)
        return () -> isoap!(ws, poly, phi, phiiso)
    elseif case == :isoap_stellatedicosahedron
        poly = stellatedicosahedron()
        phi = assign_phi(poly, 1)
        phiiso = 0.0
        ws = IsoapWorkspace(poly)
        return () -> isoap!(ws, poly, phi, phiiso)
    elseif case == :isoap_hollowedcube
        poly = hollowedcube()
        phi = assign_phi(poly, 1)
        phiiso = 0.0
        ws = IsoapWorkspace(poly)
        return () -> isoap!(ws, poly, phi, phiiso)
    elseif case == :isoap_drilledcube
        poly = drilledcube()
        phi = assign_phi(poly, 1)
        phiiso = 0.0
        ws = IsoapWorkspace(poly)
        return () -> isoap!(ws, poly, phi, phiiso)
    elseif case == :isoap_zigzagcell
        poly = zigzagcell()
        phi = assign_phi(poly, 1)
        phiiso = 0.0
        ws = IsoapWorkspace(poly)
        return () -> isoap!(ws, poly, phi, phiiso)
    elseif case == :isopol_cube
        poly = cube()
        phi = assign_phi(poly, 1)
        ia = [v > 0.0 ? 1 : 0 for v in phi]
        nts = ISOAP.nfaces(poly)
        ws = IsoapWorkspace(poly)
        return () -> isopol!(ws, ia, poly.ipv, nts)
    else
        error("Unknown case: $(case). Use one of $(join(string.(AVAILABLE_CASES), ", ")) or :all.")
    end
end

function _benchmark_case(f, current_case::Symbol, n_eval::Int, profile_n::Int)
    # Warm-up compilation before measuring runtime/allocations.
    f()

    elapsed = @elapsed begin
        for _ in 1:n_eval
            f()
        end
    end

    total_alloc_bytes = @allocated begin
        for _ in 1:n_eval
            f()
        end
    end

    println("Case: $(current_case)")
    println("Runs: $(n_eval)")
    println("Total time (s): ", elapsed)
    println("Avg time (us/run): ", elapsed * 1e6 / n_eval)
    println("Total allocated (bytes): ", total_alloc_bytes)
    println("Avg allocated (bytes/run): ", total_alloc_bytes / n_eval)

    Profile.clear()
    @profile begin
        for _ in 1:profile_n
            f()
        end
    end

    println("\nCPU profile (flat, by sample count):")
    Profile.print(format=:flat, sortedby=:count, mincount=5)
    println("\n", "-"^70, "\n")
    return nothing
end

function benchmark_profile(; case::Symbol=:all, n_eval::Int=5000, profile_n::Int=1000)
    cases = case == :all ? collect(AVAILABLE_CASES) : [case]

    for current_case in cases
        f = build_case(current_case)
        _benchmark_case(f, current_case, n_eval, profile_n)
    end

    return nothing
end

case = length(ARGS) >= 1 ? Symbol(ARGS[1]) : :all
n_eval = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 5000
profile_n = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1000

benchmark_profile(; case=case, n_eval=n_eval, profile_n=profile_n)

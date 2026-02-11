# -----------------------------------------------------------------------
# VTK output routines
# -----------------------------------------------------------------------

using Printf

"""
    polvtk(filename::AbstractString, poly::Polyhedron)
    polvtk(ifile::Int, poly::Polyhedron)

Write the geometry of a polyhedron to a VTK-format file.

If an integer `ifile` is given, the file is named `geo<ifile>.vtk` (zero-padded to 5 digits).
"""
function polvtk(filename::AbstractString, poly::Polyhedron)
    ntp = nvertices(poly)
    nts = nfaces(poly)

    open(filename, "w") do io
        println(io, "# vtk DataFile Version 2.0")
        println(io, "File: ", filename)
        println(io, "ASCII")
        println(io)
        println(io, "DATASET POLYDATA")
        println(io, "POINTS ", ntp, " float")
        for ip in 1:ntp
            x, y, z = poly.vertp[ip]
            @printf(io, "%12.6f%12.6f%12.6f\n", x, y, z)
        end

        npoly = 0
        ndata = 0
        for is in 1:nts
            n = length(poly.ipv[is])
            if n > 0
                npoly += 1
                ndata += n + 1
            end
        end

        println(io, "POLYGONS ", npoly, " ", ndata)
        for is in 1:nts
            n = length(poly.ipv[is])
            if n > 0
                println(io, n)
                for iv in 1:n
                    println(io, poly.ipv[is][iv] - 1)  # 0-based for VTK
                end
            end
        end
    end
    return nothing
end

function polvtk(ifile::Int, poly::Polyhedron)
    filename = @sprintf("geo%05d.vtk", ifile)
    polvtk(filename, poly)
end

"""
    isovtk(filename::AbstractString, result::IsoResult)
    isovtk(ifile::Int, result::IsoResult)

Write the geometry of extracted iso-polygons to a VTK-format file.

If an integer `ifile` is given, the file is named `iso<ifile>.vtk` (zero-padded to 5 digits).
"""
function isovtk(filename::AbstractString, result::IsoResult)
    niso_count = niso(result)
    if niso_count == 0
        return nothing
    end

    ntp = sum(length(p) for p in result.ipviso)

    open(filename, "w") do io
        println(io, "# vtk DataFile Version 2.0")
        println(io, "File: ", filename)
        println(io, "ASCII")
        println(io)
        println(io, "DATASET POLYDATA")
        println(io, "POINTS ", ntp, " float")
        for ip in 1:ntp
            x, y, z = result.vertiso[ip]
            @printf(io, "%12.6f%12.6f%12.6f\n", x, y, z)
        end

        npoly = 0
        ndata = 0
        for is in 1:niso_count
            n = length(result.ipviso[is])
            if n > 0
                npoly += 1
                ndata += n + 1
            end
        end

        println(io, "POLYGONS ", npoly, " ", ndata)
        for is in 1:niso_count
            n = length(result.ipviso[is])
            if n > 0
                println(io, n)
                for iv in 1:n
                    println(io, result.ipviso[is][iv] - 1)  # 0-based for VTK
                end
            end
        end
    end
    return nothing
end

function isovtk(ifile::Int, result::IsoResult)
    filename = @sprintf("iso%05d.vtk", ifile)
    isovtk(filename, result)
end

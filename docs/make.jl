using Documenter
using ISOAP

makedocs(
    modules = [ISOAP],
    authors = "Fastaxx and contributors",
    repo = "https://github.com/PenguinxCutCell/isoap.jl/blob/{commit}{path}#{line}",
    sitename = "isoap.jl",
    format = Documenter.HTML(
        canonical = "https://PenguinxCutCell.github.io/isoap.jl",
        repolink = "https://github.com/PenguinxCutCell/isoap.jl",
        collapselevel = 2,
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "Geometry Types" => "geometry-types.md",
        "Core Algorithms" => "core-algorithms.md",
        "Grid Handling" => "grid-handling.md",
        "VTK Export" => "vtk-export.md",
        "Reference" => "95-reference.md",
    ],
    pagesonly = true,
    warnonly = true,
)

deploydocs(
    repo = "github.com/PenguinxCutCell/isoap.jl",
    push_preview = true,
)

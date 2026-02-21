```@meta
CurrentModule = ISOAP
```

# isoap.jl

## Installation

```julia
using Pkg
Pkg.add("ISOAP")
```

## Overview

`isoap.jl` is a Julia package for computing the Isoap (Isosurface Approximation) algorithm, which is used for extracting isosurfaces from volumetric data.

The package includes :
- Core algorithms for isosurface extraction
- Support for various geometry types
- Grid handling utilities
- VTK export functionality for visualization

## Quick example

```julia
using ISOAP
# Create a sample grid
grid = constgrid((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (10, 10, 10))
# Define a scalar field on the grid
f = (x, y, z) -> x^2 + y^2 + z^2
# Compute the isosurface for a given isovalue
isovalue = 0.5
iso_surface = isoap(grid, f, isovalue)
# Export the isosurface to VTK format for visualization
polvtk(iso_surface, "isosurface.vtk")
```

## Main Sections

- [Getting Started](@ref getting-started)
- [Geometry Types](@ref geometry-types)
- [Core Algorithms](@ref core-algorithms)
- [Grid Handling](@ref grid-handling)
- [VTK Export](@ref vtk-export)
- [Reference](@ref reference)

## Performance

Allocations free and fast implementations.

## Reference

López, Joaquín; Hernández, Julio (2021), “isoap: A software for isosurface extraction on arbitrary polyhedra”, Mendeley Data, V1, doi: 10.17632/4rcf98s74c.1
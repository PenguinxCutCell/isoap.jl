# isoap.jl

![CI](https://github.com/PenguinxCutCell/isoap.jl/workflows/CI/badge.svg)


Julia translation of ISOAP — isosurface extraction on arbitrary polyhedra and grids.

This repository contains:
- `src/` — package source (core algorithms, mesh generators, utilities)
- `test/` — unit tests that match the original Fortran test suite

## Performance

Allocations free and fast `isoap!`and `isopol!`functions.

## License and credits

All credits go to the original ISOAP authors (J. Lopez, J. Hernandez, and collaborators).
Original Fortran code and authorship details remain with the original authors (see source headers).

## References

López, Joaquín; Hernández, Julio (2021), “isoap: A software for isosurface extraction on arbitrary polyhedra”, Mendeley Data, V1, doi: 10.17632/4rcf98s74c.1
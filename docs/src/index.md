# ThreeBodyTB.jl Documentation


**A three-body tight binding program written in Julia**

*Primary Author: Kevin F. Garrity, [NIST](https://www.nist.gov/people/kevin-garrity)*

ThreeBodyTB.jl is a package for tight-binding, written in [Julia](https://julialang.org/).

## Package Features

- Run tight-binding calculations with near DFT level accuracy (PBEsol).
- Get results in seconds based on pre-fit parameters from across periodic table.
- Calculate band structures and total energies.
- Get forces, stresses, and relax structures.
- Parameters based on two- and **three-body** interactions.
- Includes self-consistent treatment of long-range Coulomb interaction.
- Plotting based on interface of [Plots.jl](http://docs.juliaplots.org/latest/)

### New features

- Larger and better tested set of coefficients
- Coefficients for three-atom interactions for thousands of common combinations
- Magnetism (set `nspin=2` in `scf_energy` and similar functions)
- Charged unit cells
- Sparse matrix implementation of many key functions
- Improved performance

## User's guide

```@contents
Pages = ["ug_run.md", "ug_fit.md"]
Depth = 2
```

## Index

```@index
```



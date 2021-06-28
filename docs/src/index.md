# ThreeBodyTB.jl Documentation


**A three-body tight binding program written in Julia**

*Primary Author: Kevin F. Garrity, [NIST](https://www.nist.gov/people/kevin-garrity)*

!!! note

    This package is currently under development, but basic functionality should work. A manuscript is under preparation.

ThreeBodyTB.jl is a package for tight-binding, written in [Julia](https://julialang.org/).

## Package Features

- Run tight-binding calculations with near DFT level accuracy (PBEsol).
- Get results in seconds based on pre-fit parameters from across periodic table.
- Calculate band structures and total energies.
- Get forces, stresses, and relax structures.
- Parameters based on two- and **three-body** interactions.
- Includes self-consistent treatment of long-range Coulomb interaction.
- Plotting based on interface of [Plots.jl](http://docs.juliaplots.org/latest/)

## User's guide

```@contents
Pages = ["ug_run.md", "ug_fit.md"]
Depth = 2
```

## Index

```@index
```



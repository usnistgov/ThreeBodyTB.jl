# ThreeBodyTB.jl Documentation


**A three-body tight binding program written in Julia**

*Primary Maintainer: Kevin F. Garrity, [NIST](https://www.nist.gov/people/kevin-garrity)*

Other Contributors: Kamal Choudhury, NIST

!!! note

    This package is currently under development, but basic functionality should work. A manuscript has been submitted, a preprint link will be posted very shortly.

ThreeBodyTB.jl is a package for tight-binding, written in [Julia](https://julialang.org/).

## Package Features

- Run tight-binding calculations with near DFT level accuracy (PBEsol functional).
- Get results in seconds based on pre-fit parameters from across periodic table.
- Coefficients **pre-computed** for 65 elements and **any** binary combination of those elements
- Calculate band structures and total energies.
- Get forces, stresses, and relax structures.
- Parameters based on two- and **three-body** interactions.
- Includes self-consistent treatment of long-range Coulomb interaction.
- See also our python interface to this Julia code: [TB3py](https://github.com/usnistgov/tb3py)
- Based on over **800,000 DFT** calculations available at [jarvis-qetb](https://jarvis.nist.gov/jarvisqetb/)
- Plotting based on interface of [Plots.jl](http://docs.juliaplots.org/latest/)

## User's guide

```@contents
Pages = ["ug_run.md", "ug_fit.md"]
Depth = 2
```

## Index

```@index
```



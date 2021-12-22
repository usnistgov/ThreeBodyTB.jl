![asdf](https://github.com/usnistgov/ThreeBodyTB.jl/workflows/CI/badge.svg)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://pages.nist.gov/ThreeBodyTB.jl/)

<!--
[![codecov](https://codecov.io/gh/kfgarrity/ThreeBodyTB.jl/branch/main/graph/badge.svg?token=U8COIKIWG6)](https://codecov.io/gh/kfgarrity/ThreeBodyTB.jl)
-->


<!--  
[![Coverage Status](https://coveralls.io/repos/github/kfgarrity/ThreeBodyTB.jl/badge.svg?branch=main)](https://coveralls.io/github/kfgarrity/ThreeBodyTB.jl?branch=main)
[![Build Status](https://travis-ci.com/kfgarrity/ThreeBodyTB.jl.svg?branch=main)](https://travis-ci.com/kfgarrity/ThreeBodyTB.jl)
-->

# ThreeBodyTB.jl

<img align="right" src="https://github.com/kfgarrity/ThreeBodyTB.jl/blob/main/docs/src/assets/logo.svg" alt="logo" width="200" >

Code for accurate and efficient electronic structure calculations
using tight-binding (TB), including three-body interactions. Run TB
calculations with nearly DFT (pbesol) accuracy in seconds, using a
pre-fit and tested model for most elemental or binary systems.

Or fit your own coefficients!

Calculate energies, forces, and band structures.

ThreeBodyTB is written in the Julia programming language, for easy and efficient computations.

## Easy Installation

After installing julia (see [here](https://julialang.org/downloads/)), start a REPL and install using the package manager:

```
using Pkg
Pkg.add("ThreeBodyTB")
```

Alternatively, you can install from github directly

```
using Pkg
Pkg.add(url="https://github.com/usnistgov/ThreeBodyTB.jl")
```
## Documentation

See [documentation](https://pages.nist.gov/ThreeBodyTB.jl/) for more details.

## Example code

Example of creating a crystal structure of AlP (zincblende) and calculating the energy and band structure.

```
using ThreeBodyTB
c = makecrys([2.6 2.6 0;2.6 0 2.6;0 2.6 2.6], [0 0 0; 0.25 0.25 0.25], ["Al", "P"])
energy, tbc, flag = scf_energy(c)
plot_bandstr(tbc)
```

See also the [user's guide](https://pages.nist.gov/ThreeBodyTB.jl/ug_run/) for more detail.


*Performance Note:* julia compiles code the first time you run a function, but
subsequent runs of this code will take under 1 second. Also,
ThreeBodyTB.jl can take advantage of multiple processors if you define
the environment variable ```JULIA_NUM_THREADS=8;export JULIA_NUM_THREADS``` 
or start julia with ```julia --threads 8``` where
the number of threads is set as appropriate for your system.


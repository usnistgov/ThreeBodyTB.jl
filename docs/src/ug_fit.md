# Fitting tight-binding coefficients from Quantum Espresso.

ThreeBodyTB has a set of pre-fit coefficients that are sufficient for
running elemental or binary systems without doing your own fitting. If
you still want to do the fitting yourself, read on...

Currently, ThreeBodyTB is set up to fit coefficients from [Quantum
Espresso](https://www.quantum-espresso.org/) (QE) DFT calculations from `pw.x` using
the `projwfc.x` code to get atomic-wavefunction-projected band structures.

It may be easiest to consider the fitting example in the `examples/`
folder while reading this.

A brief overview of the steps:


1. Tell `ThreeBody.jl` where you QE code is, and optionally pseudopotentials, etc.

2. Run self-consistent-field (SCF) non-spin-polarized DFT total energy calculations.

3. Run [`ThreeBodyTB.AtomicProj.projwfc_workf`](@ref) to get k-space tight-binding models for each SCF calculation. This runs a non-SCF DFT calculation with `pw.x` and atomic projections with `projwfc.x`, and then calculates the TB model.

4. Run the fitting code [`ThreeBodyTB.FitTB.do_fitting_recursive`](@ref)


## 0. Preliminaries

You need to install Quantum Espresso (make pw; make pp) and let the
code know where the bin/ directory is. For example, [`set_bin_dirs(qe="/home/kfg/codes/q-e-qe-6.5/bin/")`](@ref).  You also
need to let it know about any mpi commands you may need to run in
parallel unless you use mpirun. `set_bin_dirs(qe="/home/kfg/codes/q-e-qe-6.5/bin/", mpi="mpirun -np ")`, where you must insert your appropriate command.

### Using different pseudopotential or convergence settings.

By default, the code is set up for
the slightly modified [GBRV](https://www.physics.rutgers.edu/gbrv/)
pseudopotential set available in this distribution under `pseudo/gbrv_pbesol/`.
The default template QE input files are in `template_inputs/`.
You can change those directories as well with `set_bin_dirs(pseudodir="mypsp/", templatedir="mytemp/")`

If you do change those, you will also have to edit
`ThreeBodyTB.Atomdata:atoms`. In particular, you need to perform some non-spin-polarized isolated atom
calculations to get atomic data. From those, you need the total energy, the number of
semicore states, the types of valence orbitals, the energies of
valence orbitals, and the U values. You can get U from calculations
with variable numbers of electrons.

## 1. Do DFT calculations.

You need SCF DFT data to fit to. Create crystal structures with [`makecrys`](@ref) and then run dft calculations. For example:

```
c = makecrys([10 0 0; 0 10 0; 0 0 10], [0 0 0], ["H"]);
dft = ThreeBodyTB.DFT.runSCF(c, directory="dirname", tmpdir="dirname",nprocs=8 )
```

You can look in `code_for_dataset_gen/Prototypes.jl`
and `reference_structures/` for examples of how to create a database of example structures.

You can load `dft` variables from the
"prefix.save/data-file-schema.xml" output file as 
[`dft = ThreeBodyTB.QE.loadXML("prefix.save")`]. You will need the charge
density in the next step to run the non-SCF calculation, but after
that you don't need the wavefunctions or charge density from either
the DFT or NSCF run.

## 2. Create the tight-binding Hamiltonian for each DFT calculation

The function [`ThreeBodyTB.AtomicProj.projwfc_workf`](@ref) will perform all of the
calculations to create a projected tight-binding Hamiltonian from a
DFT calculation (workf stands for workflow). The basic steps it will perform are:


1. Run a non-SCF calculation with more bands, so that every atomic orbital can project onto bands.

2. Run projwfc.x to get the atomic wavefunction projected band structure (the file `prefix.save/atomic_proj.xml`).

3. Get the projected Hamiltonian in k-space.

4. *(Optional)* Fourier transform to get the Hamiltonian in r-space

The optional Fourier transform only works if the k-grid from the NSCF is a regular gamma centered MP
grid without symmetry. However, this step isn't necessary for the main
fitting procedure. You can use the Hamiltonian in k-space directly, saving computing time.


A typical call is
```
tbc, tbck, projection_warning = ThreeBodyTB.AtomicProj.projwfc_workf(
     dft,
     directory="dirname",
     nprocs=8,
     only_kspace=true)
```

`tbck` is the important return, which is the tight binding Hamiltonian
in k-space [`tb_crys_kspace`](@ref). `projection_warning` will warn you if the quality of the
projection is detected to be low. This can be caused by atoms that are
too close together or not including enough empty bands to guarantee that
all of the atomic wavefunctions can project onto some states.

In order to
fit a SCF model, it is necessary to subtract the self-consistent part
from the TB Hamiltonian with [`ThreeBodyTB.SCF.remove_scf_from_tbc`](@ref):

``` tbck_scf = ThreeBodyTB.SCF.remove_scf_from_tbc(tbck); ```

You can read/write either `tbck` or `tbck_scf` with [`ThreeBodyTB.TB.read_tb_crys_kspace`](@ref) or [`ThreeBodyTB.TB.write_tb_crys_kspace`](@ref)

If `only_kspace=false` you will also get the real-space `tbc` output,
which is similar to the "prefix\_hr.dat" from Wannier90, and you can also use
`remove_scf_from_tbc` on that struct and use it in the fitting. The
read/write functions are [`ThreeBodyTB.read_tb_crys`](@ref) and
[`ThreeBodyTB.write_tb_crys`](@ref). These files can also be used to
interpolate band structures or DOS to arbitrary k-points like a
Wannier Hamiltonian, while the k-space versions are limited to fixed
k-point grids.


## 3. Do actual fitting

An example command to actually do the fitting is

```
database = ThreeBodyTB.FitTB.do_fitting_recursive(
	 tbck_scf_list,
	 dft_list = dft_list,
	 weights_list=weights_list,
	 starting_dict=starting_dict,
	 NLIM=50)
```

Here, `tbck_scf_list` is an array of `tbck_scf` or `tbc_scf` variables
(type `TB.tb_crys_kspace`), and `dft_list` is an array of `dft`
variables (type `DFToutMod.dftout`). `weights_list` is an optional
array of the real numbers with relative weights of the structures for
the fitting. `NLIM` is the maximum number of k-points to use for each
structure. Higher values will require more time to do the fitting.

`starting_dict` is a dictionary of previously fit coefficients that
are not included in the fitting and are kept frozen. For example, when
fitting binaries, typically the elemental fitting coefficients are
kept frozen to their elemental values.

There are many other options to this function. `ks_weight` and
`rs_weight` are the k-space and real-space weights for
fitting of the tight-binding matrix elements. Default is to set both to zero and only fit
eigenvalues/energies. `energy_weight` is the weight of the total energy,
relative to the eigenvalues. Default is `20.0`.  `niters` is the number
of recursive iterations. `lambda` is a regularization parameter,
default is `0.0` (no regularization). `RW_PARAM` is the weight of the fully empty eigenvalues
relative to the filled ones. Default is `0.0` (no weight to high energy empty states, although there is some weight to states near the Fermi level).

## 4. Run calculations with new `database`

The dictionary returned by the fitting function has coefficient
objects (type [`ThreeBodyTB.CalcTB.coefs`](@ref)) in it, as well as a variable on
whether the coefficients require self-consistency or not. For example,
if we fit Al coefficients, there will be keys with `(:Al, :Al)`, with
the two-body interaction coefs, and `(:Al, :Al, :Al)` with the three
body interactions, and a `bool` called `"SCF"`.

You can save the coefficients in xml files for use later with
`ThreeBodyTB.CalcTB.write_coefs("Al_2bdy.xml", database[(:Al,:Al)])` and read them with
`ThreeBodyTB.CalcTB.read_coefs("Al_2bdy.xml")`

You can use the coefficients to run a calculation by using
[`scf_energy(c; database = database)`](@ref), which will override the default
database and use your database. The force/stress and relaxations
similarly take `database` variables.

There is also a way to manage a directory with coefficients in it. You
have to setup your directory structure the same way as
`dats/pbesol/v1.2/`. Then, you can also load coefficients as needed
for crystal `c` stored in a directory using
[`ThreeBodyTB.ManageDatabase.prepare_database(c; directory="dirname")`](@ref) into the
internal database cache, instead of loading from the default directory
with pre-fit coefficients. Note, this will only work if the coefficients follow my naming
conventions.





# Running tight-binding calculations

How to run tight-binding calculations using the pre-fit tight-binding
coefficients. Note, only elemental and binary systems are currently
supported.

!!! note

    Running a julia function for the first time will compile the function. Future runs will be *much* faster.

## **Create a crystal object**

Use [`makecrys`](@ref) to create a [`crystal`](@ref) from lattice vectors, atomic positions, and atom types:

```@example 1
using ThreeBodyTB
A = [2.1 2.1 0.0;2.1 0.0 2.1;0.0 2.1 2.1];
pos = [0.0 0.0 0.0];
types =        ["Al"];
fcc_al = makecrys(A, pos, types)
```
Current default units are Angstrom and eV. You can change the global units to atomic units with [`set_units(both="atomic")`](@ref) if you prefer.

Alternatively, you can read the positions from a simple POSCAR or Quantum Espresso input file.

```@example 1
rbcl = makecrys("../src/POSCAR_rbcl")
```

## **Do a self-consistent calculation.**

Gets the energy and charge density with [`scf_energy`](@ref):

```@example 1
alp = makecrys("../src/POSCAR_alp")
energy, tbc_alp = scf_energy(alp); 
println("The energy is $energy eV")
```
This returns the (non-magnetic) atomization energy, and a tight-binding object with the TB matrix elements and SCF electron density calculated for post-processing.

## **Plot the band structure.**

Using the [`tb_crys`](@ref) tight-binding object `tbc_alp` from above. Note: SCF must be done first.

```@example 1
using Plots #hide
gr() #hide
ENV["GKSwstype"] = "100" #hide
plot_bandstr(tbc_alp, do_display=false); 
savefig("alp.png"); #hide
```

![AlP plot](alp.png)

Use `do_display=true` (the default) to produce an interactive plot. Here `do_display` is set to `false` because we are saving a static figure with `savefig` for the docs.

The default [`plot_bandstr`](@ref) just picks some random kpoints, but you can add your own kpath. We can also project onto the *s* orbital of Al.

```@example 1
kpath=[0.0 0.0 0.0; 0.5 0.5 0.5; 0.0 0.5 0.5];
knames=["Γ", "X", "V"];
plot_bandstr(tbc_alp, kpath=kpath, names=knames, npts=100, proj_orbs=[:s], proj_types=["Al"], do_display=false);
savefig("alp2.png"); #hide
```

![AlP plot 2](alp2.png)

You can also plot the DOS using [`dos`](@ref)

```@example 1
dos(tbc_alp, do_display=false);
savefig("alp_dos.png"); #hide
```

![AlP DOS](alp_dos.png)

Project onto orbitals instead with `proj_type=:orbs`

### Symmetry-enabled plotting

Plot bandstructure path with symmetry:
```
plot_bandstr_sym(tbc_alp, kpath=kpath, names=knames, npts=100, proj_orbs=[:s], proj_types=["Al"], do_display=false);
savefig("alp3.png"); #hide
```

![AlP plot 3](alp3.png)

Plot bandstructure with DOS:
```
plot_bandstr_dos(tbc_alp, kpath=kpath, names=knames, npts=100, proj_orbs=[:s], proj_types=["Al"], do_display=false);
savefig("alp4.png"); #hide
```

![AlP plot 4](alp4.png)



## **Calculate force / stress**

With [`scf_energy_force_stress`](@ref)

```@example 1
energy, force, stress, tbc = scf_energy_force_stress(tbc_alp);

println("energy $energy")
println()
println("Forces")
show(stdout, "text/plain", force)
println()
println("Stress")
show(stdout, "text/plain", stress)
nothing #hide
```
Can also be called directly on a new crystal structure instead of a `tb_crys` object.

## **Relax structure**

Using [`relax_structure`](@ref)

```@example 1
crys_new, tbc_updated, energy, force, stress = relax_structure(alp);

println("Energy new $energy")
println()
println("Force")
show(stdout, "text/plain", force)
println()
println("Stress")
show(stdout, "text/plain", stress)
nothing #hide
```
Energy is lower, stress is near zero, forces are zero by symmetry in Zinc Blende structure.

Force/Stress defaults are eV/Ang and eV/Ang^3.

## **Magnetic calculations**

Set `nspin=2`, default is ferromagnetic.

```@example 1
c_hatom = makecrys([10.0 0 0; 0 10.0 0; 0  0 10.0], [0 0 0], [:H])
en, tbc, flag = scf_energy(c_hatom, nspin=2)
plot_bandstr(tbc, do_display=false)
savefig("h_bands.png"); #hide
```
![H atom bands](h_bands.png)


using ThreeBodyTB
using Plots

#setup chosen Plots backend, i recommend pyplot
#pyplot()
#gr()

#make the crystal object
#here we choose fcc Al

types=["Al"];

#positions, crystal units
pos=zeros((1,3));

#lattice vectors, in Angstrom by default
A=[ [4.0 4.0 0]; [4.0 0 4.0 ]; [ 0 4.0 4.0]] * 0.529177;

#makes the crystal
c=makecrys(A, pos, types)

println("Initial Crystal:")
println()
println(c)
println()

energy, forces, stress, tbc = scf_energy_force_stress(c)


println("Initial Forces--")
println(ThreeBodyTB.Utility.arr2str(forces))
println()
println("Initial Stress--")
println(ThreeBodyTB.Utility.arr2str(stress))
println()


cfinal, tbc, energy, forces, stress = relax_structure(c)

println("Final crystal:")
println(cfinal)
println()
println("Final Forces--")
println(ThreeBodyTB.Utility.arr2str(forces))
println()
println("Finial Stress--")
println(ThreeBodyTB.Utility.arr2str(stress))
println()

#plot band structure
#kpath=[0.0 0.0 0.0;0.5 0 0; 0.5 0.25 0.75; 0.0 0.0 0.0]
#knames=["Γ", "X", "W", "Γ"]
#plot_bandstr(tbc, kpath=kpath, names=knames)




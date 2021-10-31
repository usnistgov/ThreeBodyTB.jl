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

#lattice vectors, in Angstrom units units currently
A=[ [3.8 3.8 0]; [3.8 0 3.8 ]; [ 0 3.8 3.8]] * 0.529177;

#makes the crystal
c=makecrys(A, pos, types)

println("The crystal is")
println(c)

# Do scf calculation for the energy
# returns energy, tight-binding-crystal object, and a flag (true is good)
energy, tbc, flag = scf_energy(c)

println()
println("The energy is $energy")
println()
println(tbc)


#plot band structure
kpath=[0.0 0.0 0.0;0.5 0 0; 0.5 0.25 0.75; 0.0 0.0 0.0]
knames=["Γ", "X", "W", "Γ"]

plot_bandstr(tbc, kpath=kpath, names=knames)

savefig("al_fcc.pdf")


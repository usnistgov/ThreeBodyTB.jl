using ThreeBodyTB
using Plots

#setup chosen Plots backend, i recommend pyplot
#pyplot()
#gr()


#make the crystal object
#here we choose Sc P rocksalt

types=["Sc", "P"];

#positions, crystal units
pos = [0 0 0 ; 0.5000000000  0.5000000000  0.5000000000]
#lattice vectors, in Angstrom units currently
A=[ [4.7 4.7 0]; [4.7 0 4.7 ]; [ 0 4.7 4.7]] * 1.05 * 0.529177;

#makes the crystal
c=makecrys(A, pos, types)

println("Starting crystal is")
println(c)

println("SCF ENERGY ---------------------------------------------------------------")

energy, tbc, flag = scf_energy(c)

println("BAND STRUCTURE ---------------------------------------------------------------")

#plot band structure
kpath=[0.0 0.0 0.0;0.5 0 0; 0.5 0.25 0.75; 0.0 0.0 0.0]
knames=["Γ", "X", "W", "Γ"]

plot_bandstr(tbc, kpath=kpath, names=knames, proj_types=["P"])

savefig("ScP_rocksalt.pdf")

println()
println("DOS ---------------------------------------------------------------")
#Density of states. Default is atomic projection.
dos(tbc)

savefig("ScP_rocksalt_DOS.pdf")

#dos(tbc, proj_type=:orbs) #orbital projection. :none also allowed.
#savefig("ScP_rocksalt_DOS_orbital.pdf")


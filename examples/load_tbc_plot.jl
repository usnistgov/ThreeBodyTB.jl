using ThreeBodyTB
using Plots

#setup chosen backend
#pyplot()
#gr()

# load tb_crys from file

TESTDIR = ThreeBodyTB.TESTDIR

filname = "$TESTDIR/data_forces/znse.in_vnscf_vol_2/projham.xml.gz"
tbc = read_tb_crys(filname)

println(tbc)

#these are all legal plotting commands

#plot_bandstr(tbc) #default plotting

#plot_bandstr(tbc, efermi=-0.3, align="vbm", proj_types=["Mg"]) #projection align
#plot_bandstr(tbc, efermi=-0.3, align="min") #align
#plot_bandstr(tbc, efermi=-0.3, align="fermi") #align

kpath=[0.0 0.0 0.0;0.5 0 0; 0.5 0.25 0.75; 0.0 0.0 0.0]
knames=["Γ", "X", "W", "Γ"]

plot_bandstr(tbc, kpath=kpath, names=knames, efermi=-0.3, align="vbm", proj_types=["Mg"], proj_orbs=[:s])

savefig("mgs_proj.pdf")

#compare two tb_crys
plot(reuse=false)
plot_compare_tb(tbc,tbc, efermi=-0.3, align="fermi") # compare to self

savefig("mgs_compare_self.pdf")


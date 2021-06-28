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

#load dft calculation
filname_dft = "$TESTDIR/data_forces/znse.in_vnscf_vol_2/qe.save/"
dft = ThreeBodyTB.QE.loadXML(filname_dft)

println(dft)

#compare to dft calculation
plot_compare_dft(tbc, dft)

savefig("mgs_compare_dft.pdf")




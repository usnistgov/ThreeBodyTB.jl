#directory info

SRCDIR = dirname(pathof(ThreeBodyTB))
EXAMPLESDIR = joinpath(dirname(pathof(ThreeBodyTB)), "..", "examples")
TESTDIR = joinpath(dirname(pathof(ThreeBodyTB)), "..", "test")
TEMPLATEDIR = joinpath(dirname(pathof(ThreeBodyTB)), "..", "template_inputs")
STRUCTDIR = joinpath(dirname(pathof(ThreeBodyTB)), "..", "reference_structures")

PSEUDODIR = joinpath(dirname(pathof(ThreeBodyTB)), "..", "pseudo", "gbrv_pbesol")


DATSDIR1 = joinpath(dirname(pathof(ThreeBodyTB)), "..", "dats", "pbesol", "v1.3")

#DATSDIR1 = joinpath(dirname(pathof(ThreeBodyTB)), "..", "dats", "pbesol", "v0.9")

#DATSDIR1 = joinpath(dirname(pathof(ThreeBodyTB)), "..", "dats", "pbesol", "v0.4")

DATSDIR2 = joinpath(dirname(pathof(ThreeBodyTB)), "..", "dats", "pbesol", "v1.0")

DOCSDIR = joinpath(dirname(pathof(ThreeBodyTB)), "..", "docs")


#DATSDIR1 = "/home/kfg/codes/TB_fit/fit_dimer/datab_2/"
#DATSDIR2 = ""

#DATSDIR2 = "/home/kfg/codes/TB_fit/binary_v10/datab/"


global MPI_STRING="mpirun -np"
#global QE_BIN_DIR_STRING="/home/kfg/codes/q-e-qe-6.5/bin/"
global QE_BIN_DIR_STRING="/usr/local/almalinux9/qe/6.8/openmpi-4.1.5-gcc-9/bin/"
global PSEUDOS=PSEUDODIR
global TEMPLATES=TEMPLATEDIR

global WANNIER_BIN_DIR_STRING="/users/kfg/codes/wannier90-2.1.0/" #not necessary


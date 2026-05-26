
#default values

using ..ThreeBodyTB:QE_BIN_DIR_STRING
using ..ThreeBodyTB:WANNIER_BIN_DIR_STRING
using ..ThreeBodyTB:MPI_STRING

"""
    function make_commands(nprocs=1)

Returns a dictionary with command lines to call external programs on the command line
nprocs is the number of processors for parallel execution
Needs to be changed for you specfic program locations 
and mpi commands (if any)

This needs to be edited to actually run QE yourself.
Running wannier90 is optional, not part of current code.
"""
function make_commands(nprocs=1)
"""
Returns a dictionary with command lines to call external programs on the command line
nprocs is the number of processors for parallel execution
Needs to be changed for you specfic program locations 
and mpi commands (if any)
"""



    #    mpi=`mpirun -n `
#    mpi=`$MPI_STRING`    
    mpi=split(MPI_STRING)
#    qebin = "/users/kfg/codes/q-e-qe-6.3.rc1/bin/"
#    w90bin= "/users/kfg/codes/wannier90-2.1.0/"

    qebin = QE_BIN_DIR_STRING
    w90bin= WANNIER_BIN_DIR_STRING


    juiladir = "../"
    

    
    #main QE SCF driver
    if length(QE_BIN_DIR_STRING) == 0
        pw = "pw.x"
    else
        pw = "$qebin/pw.x"
    end


    #OMP="export OMP_NUM_THREADS=1;"
    OMP=split("OMP_NUM_THREADS=1 ")
#    pwscf_command_serial=`OMP_NUM_THREADS=1\; $mpi 1  $pw -input `
#    pwscf_command_parallel=` OMP_NUM_THREADS=1\;  $mpi $nprocs $pw -npool 2 -input `#

#    pwscf_command_parallel_backup=`OMP_NUM_THREADS=1\; $mpi $nprocs   $pw -ndiag 1 -npool 1 -input `


#    pwscf_command_serial= (fil, procs) -> `bash -c 'export OMP_NUM_THREADS=1; $mpi 1  $pw -input $fil '`
#    pwscf_command_parallel= (fil, procs) -> `bash -c 'export OMP_NUM_THREADS=1; $mpi $procs  $pw -npool 2  -input $fil '`
#    pwscf_command_parallel_backup= (fil, procs) -> `bash -c 'export OMP_NUM_THREADS=1; $mpi $procs  $pw -ndiag 1 -input $fil' `

    pwscf_command_serial= (fil, procs) -> `bash -c 'export OMP_NUM_THREADS=1; mpirun -np 1  pw.x -input $fil '`

    function f(fil, procs); println("in f $fil $procs"); st = "export OMP_NUM_THREADS=1; mpirun -np $procs  pw.x -ndiag 1 -input $fil"; println("st"); println(st); println(); return `bash -c $st `; end
    pwscf_command_parallel=f
    #pwscf_command_parallel= (fil, procs) -> `bash -c 'export OMP_NUM_THREADS=1; mpirun $procs  pw.x -npool 2  -input $fil '`
    pwscf_command_parallel_backup= (fil, procs) -> `bash -c 'export OMP_NUM_THREADS=1; mpirun -np $procs  pw.x -ndiag 1 -input $fil' `
    
    #pwscf_command_serial=`bash -c 'export OMP_NUM_THREADS=1; $mpi 1  $pw -input `
    #pwscf_command_parallel=` OMP_NUM_THREADS=1\;  $mpi $nprocs $pw -npool 2 -input `#
    #pwscf_command_parallel_backup=`OMP_NUM_THREADS=1\; $mpi $nprocs   $pw -ndiag 1 -npool 1 -input `
    
    
    #qe-to-wannier90 code
    pw2wan_command_serial=`$qebin/pw2wannier90.x -input `
    pw2wan_command_parallel=`$mpi $nprocs $qebin/pw2wannier90.x -input `

    #uses symmetry to transform QE scf calculation into full k-point grid, so that w90 can understand (not currently in use)
    og_command_serial=`$qebin/open_grid.x -input `
    og_command_parallel=`$qebin/open_grid.x -input `
    #    og_command_parallel=`$mpi $nprocs $qebin/open_grid.x -input `

    #qe-project-wavefunction code

    if length(QE_BIN_DIR_STRING) == 0
        proj = "projwfc.x"
    else
        proj = "$qebin/projwfc.x"
    end

    

    function fproj(fil, procs); procs = 1; println("in fproj $fil $procs"); st = "export OMP_NUM_THREADS=1; mpirun -np $procs  projwfc.x -nd 1 -input $fil"; println("st"); println(st); println(); return `bash -c $st `; end
    proj_command_serial=fproj
    proj_command_parallel=fproj
    
#    proj_command_serial=` $mpi 1 OMP_NUM_THREADS=1  $proj -nd 1 -input `
#    proj_command_parallel=` $mpi $nprocs OMP_NUM_THREADS=1  $proj -nd 1 -input `

    proj_command_serial_backup=` $proj  -input `
    proj_command_parallel_backup=` $mpi 1 $proj  -input `
    
    
    #w90 (serial)
    wannier90_command=`$w90bin/wannier90.x `


    ############################################# Change above here for specific computer system

    command_dict = Dict()

#    command_dict["juliadir"] = juiladir

    if nprocs == 1

        command_dict["qe"] = pwscf_command_serial
        command_dict["qe_backup"] = pwscf_command_serial
        command_dict["pw2wan"] = pw2wan_command_serial        
        command_dict["og"] = og_command_serial
        command_dict["proj"] = proj_command_serial
        command_dict["proj_serial"] = proj_command_parallel
        command_dict["proj_backup"] = proj_command_serial
        command_dict["proj_serial_backup"] = proj_command_parallel

    else

        command_dict["qe"] = pwscf_command_parallel
        command_dict["qe_backup"] = pwscf_command_parallel_backup
        command_dict["pw2wan"] = pw2wan_command_parallel        
        command_dict["og"] = og_command_parallel
        command_dict["proj"] = proj_command_parallel        
        command_dict["proj_serial"] = proj_command_serial
        command_dict["proj_backup"] = proj_command_serial
        command_dict["proj_serial_backup"] = proj_command_parallel

    end

    command_dict["wannier90"] = wannier90_command


    return command_dict

end    


        

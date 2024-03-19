
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
    pwscf_command_serial=`$qebin/pw.x -input `
    pwscf_command_parallel=`$mpi $nprocs $qebin/pw.x -npool 2 -input `

    pwscf_command_parallel_backup=`$mpi $nprocs $qebin/pw.x -ndiag 1 -npool 1 -input `

    
    #qe-to-wannier90 code
    pw2wan_command_serial=`$qebin/pw2wannier90.x -input `
    pw2wan_command_parallel=`$mpi $nprocs $qebin/pw2wannier90.x -input `

    #uses symmetry to transform QE scf calculation into full k-point grid, so that w90 can understand (not currently in use)
    og_command_serial=`$qebin/open_grid.x -input `
    og_command_parallel=`$qebin/open_grid.x -input `
    #    og_command_parallel=`$mpi $nprocs $qebin/open_grid.x -input `

    #qe-project-wavefunction code
    proj_command_serial=`$qebin/projwfc.x -nd 1 -input `
    proj_command_parallel=`$mpi $nprocs $qebin/projwfc.x -nd 1 -input `

    
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
        #command_dict["proj"] = proj_command_serial
        command_dict["proj"] = proj_command_parallel

    else

        command_dict["qe"] = pwscf_command_parallel
        command_dict["qe_backup"] = pwscf_command_parallel_backup
        command_dict["pw2wan"] = pw2wan_command_parallel        
        command_dict["og"] = og_command_parallel
        command_dict["proj"] = proj_command_parallel        
        #command_dict["proj"] = proj_command_serial

    end

    command_dict["wannier90"] = wannier90_command


    return command_dict

end    


        

"""
    module ThreeBodyTB

Main Module
"""
module ThreeBodyTB

#println("develop3 branch")

function get_ver()
    println("Main branch. v1.0")
end

#get_ver()


include("GlobalUnits.jl")
include("SetDir.jl")




include("Utility.jl")
include("BandTools.jl")
include("Atomic.jl")
include("Atomdata.jl")
include("AtomicMag.jl")
include("Crystal.jl")
include("Symmetry.jl")


"""
    function set_units(;energy=missing, length=missing, both=missing)

Set global units for energy/length. Run with no arguments to check/return current units.

- Default units are `"eV"` and `"Angstrom"` (or `"Å"` or `"Ang"` ).
- Choose atomic units by setting `energy="Ryd"` and `length="Bohr"`.
- Set both at the same time with `both="atomic"` or `both="eVAng"`
- Internally, all units are atomic. Only main public facing functions actually change units.
"""
function set_units(a=missing;energy=missing, length=missing, both=missing)

    if !ismissing(a)
        both = a
    end
    
    if !ismissing(both)
        both = String(both)
        if both == "atomic" || both == "Atomic" || both == "au" || both == "a.u."
            energy="Ryd."
            length="Bohr"
        elseif both == "eVang" || both == "AngeV" || both == "eVAng" || both == "evang"
            energy="eV"
            length="Å"
        else
            println("I don't understand both : ", both)
        end
            
    end
    
    if !ismissing(length)
        length = String(length)
        if length == "A" || length == "Ang" || length == "Ang." || length == "Angstrom" || length == "a" || length == "ang" || length == "ang." || length == "angstrom" || length == "Å"

            global global_length_units="Å"
        elseif length == "Bohr" || length == "bohr" || length == "au" || length == "a.u." || length == "atomic"
            global global_length_units="Bohr"
        else
            println("I don't understand length : ", length)
        end
    end

    if !ismissing(energy)
        energy = String(energy)
        if energy == "eV" || energy == "ev" || energy == "EV" || energy == "electronvolts"
            global global_energy_units="eV"
        elseif energy == "Rydberg" || energy == "Ryd." || energy == "Ryd" || energy == "atomic" || energy == "au" || energy == "a.u." || energy == "rydberg" || energy == "ryd" || energy == "ryd."
            global global_energy_units="Ryd."
        else
            println("I don't understand energy : ", energy)
        end
            
    end

    println("Units are now $global_energy_units and $global_length_units ")

    return deepcopy(global_energy_units), deepcopy(global_length_units)
    
end


using Suppressor

using Printf
using .CrystalMod:crystal
using .CrystalMod:makecrys

export makecrys
export crystal

include("Ewald.jl")
include("DFToutMod.jl")
include("RunDFT.jl")

using .DFToutMod:dftout
using .DFToutMod:makedftout
export dftout
export makedftout

include("TB.jl")
###include("Magnetic.jl") #now inside TB.jl

using .TB:tb_crys
using .TB:tb_crys_dense
using .TB:tb_crys_sparse
using .TB:tb_crys_kspace
using .TB:tb

export tb_crys
export tb_crys_kspace
export tb

using .TB:Hk
using .TB:calc_bands
using .TB:read_tb_crys
using .TB:write_tb_crys
using .TB:get_formation_energy

include("AtomicProj.jl")

include("CalcTB_laguerre.jl")

#using .CalcTB:calc_tb_fast
#export calc_tb_fast
using .CalcTB:calc_twobody

include("ManageDatabase.jl")


include("Classical.jl")

using .Classical:energy_force_stress_cl

include("SCF.jl")
include("Force_Stress.jl")

include("DOS.jl")
using .DOS:dos
using .DOS:dos_realspace

include("BandStruct.jl")

using .BandStruct:plot_compare_tb
using .BandStruct:plot_bandstr
using .BandStruct:plot_compare_dft
#using .BandStruct:set_no_display
using .BandStruct:band_summary
using .BandStruct:plot_bandstr_dos
using .BandStruct:plot_bandstr_sym



export Hk
export calc_bands
export plot_compare_tb
export plot_bandstr
export plot_compare_dft
export plot_bandstr_dos
export plot_bandstr_sym
export read_tb_crys
export dos
export dos_realspace
export write_tb_crys
export band_summary

using .TB:write_hr_dat
export write_hr_dat


export set_units
export set_bin_dirs

using .CrystalMod:plot
export plot

#include("RunWannier90.jl")


include("FitTB_laguerre.jl")





include("MyOptim.jl")

include("Relax.jl")
using .CrystalMod:print_with_force_stress

include("compile.jl")

export scf_energy
export scf_energy_force_stress
export relax_structure



#function amiworking()
#    println("yes")
#    pwd()
#    println(pwd())
#end



"""
    function set_bin_dirs(;qe=missing, mpi=missing, pseudodir=missing, templatedir=missing, wannier=missing)

Set directories where things like quantum espresso bins are located. If run with everything missing, instead print current dirs.

- `qe` - set bin directory of quantum espresso. Needs pw and pp installed. No useful default.
- `mpi` - set mpi command. something like "mpirun -np "
- `pseudo` - set directory for pseudopotentials. Default is to use ../pseudo/gbrv_pbesol/ as distributed with this code
- `template` - set directory for quantum espresso template files. Uses ../template_inputs/ by default.

"""
function set_bin_dirs(;qe=missing, mpi=missing, pseudodir=missing, templatedir=missing, wannier=missing)

    if !ismissing(qe)
        global QE_BIN_DIR_STRING = qe
        println("set qe bin dir to $QE_BIN_DIR_STRING")
    end

    if !ismissing(mpi)
        global MPI_STRING = mpi
        println("set mpi string to $MPI_STRING")        
    end

    if !ismissing(pseudodir)
        global PSEUDOS = pseudodir
        println("set pseudopotential dir to $PSEUDOS")
    end

    if !ismissing(templatedir)
        global TEMPLATEDIR = templatedir
        println("set directory of QE templates to $TEMPLATEDIR")
    end

    if !ismissing(wannier)
        global WANNIER_BIN_DIR_STRING = wannier
        println("defunct: set wannier bin dir to $WANNIER_BIN_DIR_STRING")
    end

    if ismissing(qe) && ismissing(mpi) && ismissing(pseudodir) && ismissing(templatedir) && ismissing(wannier)
        println()
        println("Report:")
        println()
        println(QE_BIN_DIR_STRING)
        println(MPI_STRING)
        println(PSEUDOS)
        println(TEMPLATEDIR)
#        println(WANNIER_BIN_DIR_STRING)
        println()
    end
    
end



"""
    relax_structure(c::crystal; mode="vc-relax")

Find the lowest energy atomic configuration of crystal `c`.

# Arguments
- `c::crystal`: the structure to relax, only required argument
- `mode="vc-relax"`: Default (variable-cell relax) will relax structure and cell, anything else will relax structure only.
- `database=missing`: coefficent database, default is to use the pre-fit pbesol database
- `smearing=0.01`: smearing temperature (ryd), default = 0.01
- `grid=missing`: k-point grid, e.g. [10,10,10], default chosen automatically
- `nsteps=100`: maximum iterations
- `update_grid=true`: update automatic k-point grid during relaxation
- `conv_thr = 2e-3 `: Convergence threshold for gradient
- `energy_conv_thr = 2e-4 `: Convergence threshold for energy in Ryd
- `sparse = :auto`: Default is to use dense matricies for `nat < 100`. Can be `true` or `false` to force choice.
"""
function relax_structure(c::crystal; database=missing, smearing = 0.01, grid = missing, mode="vc-relax", nsteps=50, update_grid=true, conv_thr = 2e-3, energy_conv_thr = 2e-4, nspin=1, repel=true, do_tb=true, database_classical=missing, do_classical=true, tot_charge=0.0, sparse=:auto)

    if ismissing(database_classical)
        do_classical=false
    end

    if sparse == :auto
        if c.nat >= 100
            sparse = true
            println("auto use sparse matricies, set sparse=false to avoid")
        else
            sparse = false
        end
    end
    
    
    if ismissing(database)
        ManageDatabase.prepare_database(c)
        database = ManageDatabase.database_cached
    end
    if !ismissing(grid)
        update_grid = false
    end

    cfinal, tbc, energy, force, stress = Relax.relax_structure(c, database, smearing=smearing, grid=grid, mode=mode, nsteps=nsteps, update_grid=update_grid, conv_thr=conv_thr, energy_conv_thr = energy_conv_thr, nspin=nspin, repel=repel, do_tb=do_tb, database_classical=database_classical, do_classical=do_classical, tot_charge=tot_charge, sparse=sparse)

   
    println("Relax done")

    println("energy $energy")
    println("force ", force)
    println("stress ", stress)
    
#    println("Calculate final energy")
#
#    energy_tot, tbc, conv_flag = scf_energy(cfinal; database=database, smearing=0.01, grid = missing)
#    energy_tot, f_cart, stress = scf_energy_force_stress(tbc, database=database, smearing=smearing, grid=grid)
#
#    println("done with all relax")
#    println()


#    energy = convert_energy(energy)
    println()
    println("Formation energy: " , round(convert_energy(get_formation_energy(energy, cfinal)), digits=3) , " $global_energy_units/atom" )
    println()


    println()
    println("---------------------------------")
    println("Final Energy $energy ")
    
    print_with_force_stress(cfinal, force, stress)
    println()

    energy = convert_energy(energy)
    force = convert_force(force)
    stress = convert_stress(stress)



    

    

    println()
    
    
    return cfinal, tbc, energy, force, stress

end

"""
    scf_energy_force_stress(c::crystal; database = missing, smearing = 0.01, grid = missing)

Calculate energy, force, and stress for a crystal.

`return energy, force, stress, tight_binding_crystal_object`


# Arguments
- `c::crystal`: the structure to calculate on. Only required argument.
- `database=missing`: Source of coeficients. Will be loaded from pre-fit coefficients if missing.
- `smearing=0.01`: Gaussian smearing temperature, in Ryd. Usually can leave as default.
- `grid=missing`: k-point grid, e.g. [10,10,10], default chosen automatically
- `sparse = :auto`: Default is to use dense matricies for `nat < 100`. Can be `true` or `false` to force choice.
"""
function scf_energy_force_stress(c::crystal; database = missing, smearing = 0.01, grid = missing, nspin=1, repel=true , use_sym=true, verbose=false, do_classical=true, database_classical=missing, do_tb=true, tot_charge=0.0, sparse = :auto)

    #nothing case
    if !do_tb && !do_classical
        println("WARNING Nothing to do; do_tb=$do_tb do_classical=$do_classical")
        return 0.0, zeros(c.nat, 3), zeros(3,3), missing
    end

    if ismissing(database_classical)
        do_classical=false
    end

    if sparse == :auto
        if c.nat >= 100
            sparse = true
            println("auto use sparse matricies, set sparse=false to avoid")
        else
            sparse = false
        end
    end
    
    if sparse
        use_sym = true
    end
    
    if do_classical
        energy_cl, force_cl, stress_cl = energy_force_stress_cl(c, database=database_classical)
    else
        energy_cl=0.0
        force_cl=zeros(c.nat, 3)
        stress_cl=zeros(3,3)
    end
    
    #classical only case
    if !do_tb
        if ismissing(database_classical)
            println("WARNING ismissing database_classical")
        end
        return energy_cl, force_cl, stress_cl, missing
    end
#    println("tot_charge before 00000 ", tot_charge)
    energy_tot, tbc, conv_flag = scf_energy(c; database=database, smearing=smearing, grid = grid, nspin=nspin, conv_thr=1e-6, verbose=verbose, repel=repel, use_sym=use_sym, tot_charge=tot_charge, sparse=sparse )
#    println("tot charge 11111111111 ", tbc.tot_charge, " " , tbc.nelec)
    
    if ismissing(database)
        database = ManageDatabase.database_cached
    end

    #this is a compatibility issue, haven't programmed symmetry for NON-scf case.
    if database["SCF"] == false
        use_sym=false
    end
    
    println()
    println("Calculate Force, Stress")
    if use_sym
#        println("use sym")
        energy_tot, f_cart, stress = Force_Stress.get_energy_force_stress_fft_LV_sym_SINGLE(tbc, database, do_scf=false, smearing=smearing, grid=grid, nspin=nspin, repel=repel)
    else
        energy_tot, f_cart, stress = Force_Stress.get_energy_force_stress_fft_LV(tbc, database, do_scf=false, smearing=smearing, grid=grid, nspin=nspin, repel=repel)        
    end
#    println("tot charge 22222222 ", tbc.tot_charge, " " , tbc.nelec)

    if do_classical
        energy_tot += energy_cl
        f_cart += force_cl
        stress += stress_cl
    end
    
    println("done")
    println("----")
    println()
    println("Formation energy: " , round(convert_energy(get_formation_energy(energy_tot, c)), digits=3) , " $global_energy_units/atom" )
    println()


    
    energy_tot = convert_energy(energy_tot)
    f_cart = convert_force(f_cart)
    stress = convert_stress(stress)

    print_with_force_stress(c, f_cart, stress)

#    println("tot charge 333333 ", tbc.tot_charge, " " , tbc.nelec)
    
    return energy_tot, f_cart, stress, tbc

end

"""
    scf_energy_force_stress(tbc::tb_crys; database = missing, smearing = 0.01, grid = missing)

Calculate energy, force, and stress for a tight binding crystal
object. This allows the calculation to run without re-doing the SCF
calculation. Assumes SCF already done!

returns energy, force, stress, tight_binding_crystal_object

"""
function scf_energy_force_stress(tbc::tb_crys; database = missing, smearing = 0.01, grid = missing, do_scf=false, repel=true, use_sym=true,do_classical=true, database_classical=missing)

    
    if ismissing(database_classical)
        do_classical=false
    end

    if typeof(tbc) <: tb_crys_sparse
        use_sym=true
    end
    
    if ismissing(database)
        ManageDatabase.prepare_database(tbc.crys)
        database = ManageDatabase.database_cached
    end

    println()
    println("Calculate Force, Stress (no scf)")

    if use_sym
        energy_tot, f_cart, stress = Force_Stress.get_energy_force_stress_fft_LV_sym_SINGLE(tbc, database, do_scf=true, smearing=smearing, grid=grid, nspin=size(tbc.eden)[1], repel=repel)
    else
        energy_tot, f_cart, stress = Force_Stress.get_energy_force_stress_fft_LV(tbc, database, do_scf=true, smearing=smearing, grid=grid, nspin=size(tbc.eden)[1], repel=repel)
    end

    if do_classical && !ismissing(database_classical)
        energy_cl, force_cl, stress_cl = energy_force_stress_cl(tbc.crys, database=database_classical)
        energy_tot += energy_cl
        f_cart += force_cl
        stress += stress_cl
        
    end

    
    println("done")
    println("----")
    println()
    println("Formation energy: " , round(convert_energy(get_formation_energy(energy_tot, tbc.crys)),digits=3) , " $global_energy_units/atom" )
    println()

    energy_tot = convert_energy(energy_tot)
    f_cart = convert_force(f_cart)
    stress = convert_stress(stress)

    print_with_force_stress(tbc.crys, f_cart, stress)


    return energy_tot, f_cart, stress, tbc

end




"""
    scf_energy(c::crystal)

Calculate energy, force, and stress for a crystal. Does self-consistent-field (SCF) calculation if using self-consistent electrostatics.

returns energy, tight-binding-crystal-object, error-flag

# Arguments
- `c::crystal`: the structure to calculate on. Only required argument.
- `database=missing`: Source of coeficients. Will be loaded from pre-fit coefficients if missing.
- `smearing=0.01`: Gaussian smearing temperature, in Ryd. Usually can leave as default.
- `grid=missing`: k-point grid, e.g. [10,10,10], default chosen automatically
- `conv_thr = 1e-5`: SCF convergence threshold (Ryd).
- `sparse = :auto`: Default is to use dense matricies for `nat < 100`. Can be `true` or `false` to force choice.
- `iter = 75`: number of iterations before switch to more conservative settings.
- `mix = -1.0`: initial mixing. -1.0 means use default mixing. Will automagically adjust mixing if SCF is failing to converge. Starting default is smaller for larger unit cells.
- `mixing_mode = :simple`: default is simple. Other options are `:simple` and `:DIIS` / `:pulay` (direct inversion of iterative subspace). Will automatically switch to simple if Pulay fails. 
"""
function scf_energy(c::crystal; database = missing, smearing=0.01, grid = missing, conv_thr = 2e-5, iters = 100, mix = -1.0, mixing_mode=:simple, nspin=1, eden=missing, verbose=false, repel=true, tot_charge=0.0, use_sym=true, do_classical=true, do_tb=true, database_classical=missing, sparse=:auto)
    println()
#    println("Begin scf_energy-------------")
#    println()
    if ismissing(database)
        if verbose; println("Load TB parameters from file"); end
        ManageDatabase.prepare_database(c)
        database = ManageDatabase.database_cached
        println()
    end
    
    if ismissing(database_classical)
        do_classical=false
    end
        
    if sparse == :auto
        if c.nat >= 100
            sparse = true
            println("auto use sparse matricies, set sparse=false to avoid")
        else
            sparse = false
        end
    end
    
    energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbc = SCF.scf_energy(c, database, smearing=smearing, grid = grid, conv_thr = conv_thr, iters = iters, mix = mix,  mixing_mode=mixing_mode, nspin=nspin, e_den0=eden, verbose=verbose, repel=repel, tot_charge=tot_charge, use_sym=use_sym, do_classical=do_classical, database_classical=database_classical, do_tb=do_tb, sparse=sparse)

    conv_flag = !error_flag
    if do_tb
        if tbc.within_fit == false
            conv_flag = false
        end
    end
    if conv_flag == true
        println("scf_energy success, done")
    else
        println("scf_energy error detected somewhere!!!!!!!!!!, done")
    end
    println()
    println("Formation energy: " , round(convert_energy(get_formation_energy(energy_tot, c)), digits=3) , " $global_energy_units/atom" )
    println()



    energy_tot = convert_energy(energy_tot)

    return energy_tot, tbc, conv_flag

end


"""
    scf_energy(d::dftout)

    SCF energy using crystal structure from DFT object.
"""
function scf_energy(d::dftout; database = Dict(), smearing=0.01, grid = missing, conv_thr = 2e-5, iters = 75, mix = -1.0, mixing_mode=:DIIS, nspin=1, verbose=true, repel=true, use_sym=true, do_classical=true, database_classical=missing, do_tb=true, sparse=:auto)

    if ismissing(database_classical)
        do_classical=false
    end
        
    
    return scf_energy(d.crys, smearing=smearing, grid = grid, conv_thr = conv_thr, iters = iters, mix = mix, mixing_mode=mixing_mode, nspin=nspin, verbose=verbose, repel=repel, tot_charge=dft.tot_charge, use_sym=use_sym, do_classical=do_classical, database_classical=database_classical,do_tb=do_tb, sparse=sparse)

end


"""
    scf_energy(tbc::tb_crys)

SCF energy using crystal structure from TBC object.
"""
function scf_energy(tbc::tb_crys; smearing=0.01, grid = missing, e_den0 = missing, conv_thr = 2e-5, iters = 75, mix = -1.0, mixing_mode=:DIIS, nspin=1, verbose=true, tot_charge=missing, use_sym=true, do_classical=true, database_classical=missing, repel=true)

    if ismissing(database_classical)
        do_classical=false
    end
        
    
    energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbc = SCF.scf_energy(tbc; smearing=smearing, grid = grid, e_den0 = e_den0, conv_thr = conv_thr, iters = iters, mix = mix, mixing_mode=mixing_mode, nspin=nspin, verbose=verbose, tot_charge=missing, use_sym=use_sym, database_classical=database_classical)
    println()
    println("Formation energy: " , round(convert_energy(get_formation_energy(energy_tot, tbc.crys)), digits=3) , " $global_energy_units/atom" )
    println()

    conv_flag = !error_flag
    if tbc.within_fit == false
        conv_flag = false
    end
    if conv_flag == true
        println("success, done")
    else
        println("error detected!!!!!!!!!!, done")
    end


    energy_tot = convert_energy(energy_tot)
    
    return energy_tot, tbc, conv_flag

end

"""
    function precompile()

For precompiling in order to make a sysimage
"""
function my_precompile()

    set_no_display(true)
    for t in ["calculate_energy_plot_fcc_al.jl", "load_dft_plot.jl", "load_tbc_plot.jl", "plot_Sc_P.jl", "relax_fcc_al.jl"]
        println("precomp example $t")
        @suppress include(joinpath( ThreeBodyTB.EXAMPLESDIR ,t))
    end
    set_no_display(false)

    Nothing
end

"""
    function calc_twobody(t1::String,t2::String,orb1::String,orb2::String,dist,lmn; database=missing)

Get twobody terms. Returns hamiltonian and overlap (H,S) between atom t1 orb1 and atom t2 orb2, that are dist appart with direction cosines lmn
"""
function get_twobody(t1::String,t2::String,orb1::String,orb2::String,dist::Number,lmn; database=missing)

    if ismissing(database)
        ManageDatabase.prepare_database([t1,t2])
        database = ManageDatabase.database_cached
    end

    if global_length_units == "Å"
        dist = dist / Ang
    end
    h,s = calc_twobody(Symbol(t1),Symbol(t2),Symbol(orb1), Symbol(orb2), dist, lmn, database)
    return convert_energy(h), s
end

"""
    function calc_twobody(t1::String,t2::String,orb1::String,orb2::String,R; database=missing)

Get twobody terms. Returns hamiltonian and overlap (H,S) between atom t1 orb1 and atom t2 orb2, that are positioned at R = R2 - R1 (R is an array with 3 reals). Does not apply to onsite matrix elements
"""

function get_twobody(t1::String,t2::String,orb1::String,orb2::String, R; database=missing)

    dist = sum(R.^2)^0.5
    if dist > 1e-5
        lmn = R / dist
    else
        lmn = zero(3)
        println("warning, doesn't apply to onsite $dist")
    end
    return get_twobody(t1,t2,orb1, orb2, dist, lmn, database=database)
    
end


include("ClassicalFit.jl")


end #end module


"""
    module ThreeBodyTB

Main Module
"""
module ThreeBodyTB


include("GlobalUnits.jl")
include("SetDir.jl")

include("Utility.jl")
include("BandTools.jl")
include("Atomic.jl")
include("Atomdata.jl")
include("AtomicMag.jl")
include("Crystal.jl")

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
using .TB:tb_crys_kspace
using .TB:tb

export tb_crys
export tb_crys_kspace
export tb

using .TB:Hk
using .TB:calc_bands
using .TB:read_tb_crys
using .TB:write_tb_crys

include("AtomicProj.jl")
include("BandStruct.jl")

using .BandStruct:plot_compare_tb
using .BandStruct:plot_bandstr
using .BandStruct:plot_compare_dft
using .BandStruct:set_no_display
using .BandStruct:band_summary

include("DOS.jl")
using .DOS:dos

export Hk
export calc_bands
export plot_compare_tb
export plot_bandstr
export plot_compare_dft
export read_tb_crys
export dos
export write_tb_crys
export band_summary


export set_units
export set_bin_dirs


#include("RunWannier90.jl")
include("CalcTB_laguerre.jl")
using .CalcTB:calc_tb_fast
export calc_tb_fast
using .CalcTB:calc_twobody
include("SCF.jl")


include("FitTB_laguerre.jl")
include("Force_Stress.jl")




include("ManageDatabase.jl")

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
"""
function relax_structure(c::crystal; database=missing, smearing = 0.01, grid = missing, mode="vc-relax", nsteps=50, update_grid=true, conv_thr = 2e-3, energy_conv_thr = 2e-4, nspin=1)

    if ismissing(database)
        ManageDatabase.prepare_database(c)
        database = ManageDatabase.database_cached
    end
    if !ismissing(grid)
        update_grid = false
    end

    cfinal, tbc, energy, force, stress = Relax.relax_structure(c, database, smearing=smearing, grid=grid, mode=mode, nsteps=nsteps, update_grid=update_grid, conv_thr=conv_thr, energy_conv_thr = energy_conv_thr, nspin=nspin)

   
    println("Relax done")

#    println("Calculate final energy")
#
#    energy_tot, tbc, conv_flag = scf_energy(cfinal; database=database, smearing=0.01, grid = missing)
#    energy_tot, f_cart, stress = scf_energy_force_stress(tbc, database=database, smearing=smearing, grid=grid)
#
#    println("done with all relax")
#    println()


#    energy = convert_energy(energy)

    energy = convert_energy(energy)
    force = convert_force(force)
    stress = convert_stress(stress)

    println()
    println("---------------------------------")
    println("Final Energy $energy ")
    
    print_with_force_stress(cfinal, force, stress)
    println()


    

    

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
"""
function scf_energy_force_stress(c::crystal; database = missing, smearing = 0.01, grid = missing, nspin=1)
    
    energy_tot, tbc, conv_flag = scf_energy(c; database=database, smearing=smearing, grid = grid, nspin=nspin, conv_thr=1e-7, verbose=false)

    if ismissing(database)
        database = ManageDatabase.database_cached
    end

    println()
    println("Calculate Force, Stress")
    
    energy_tot, f_cart, stress = Force_Stress.get_energy_force_stress_fft(tbc, database, do_scf=false, smearing=smearing, grid=grid, nspin=nspin)

    println("done")
    println("----")

    
    energy_tot = convert_energy(energy_tot)
    f_cart = convert_force(f_cart)
    stress = convert_stress(stress)

    print_with_force_stress(c, f_cart, stress)
    
    return energy_tot, f_cart, stress, tbc

end

"""
    scf_energy_force_stress(tbc::tb_crys; database = missing, smearing = 0.01, grid = missing)

Calculate energy, force, and stress for a tight binding crystal
object. This allows the calculation to run without re-doing the SCF
calculation. Assumes SCF already done!

returns energy, force, stress, tight_binding_crystal_object

"""
function scf_energy_force_stress(tbc::tb_crys; database = missing, smearing = 0.01, grid = missing, do_scf=false)
    
    if ismissing(database)
        ManageDatabase.prepare_database(tbc.crys)
        database = ManageDatabase.database_cached
    end

    println()
    println("Calculate Force, Stress (no scf)")
    
    energy_tot, f_cart, stress = Force_Stress.get_energy_force_stress_fft(tbc, database, do_scf=false, smearing=smearing, grid=grid, nspin=size(tbc.eden)[1])

    println("done")
    println("----")

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
- `iter = 75`: number of iterations before switch to more conservative settings.
- `mix = -1.0`: initial mixing. -1.0 means use default mixing. Will automagically adjust mixing if SCF is failing to converge.
- `mixing_mode =:pulay`: default is Pulay mixing (DIIS). Other option is :simple, for simple linear mixing of old and new electron-density. Will automatically switch to simple if Pulay fails.
"""
function scf_energy(c::crystal; database = missing, smearing=0.01, grid = missing, conv_thr = 1e-5, iters = 75, mix = -1.0, mixing_mode=:pulay, nspin=1, eden=missing, verbose=false)
    println()
#    println("Begin scf_energy-------------")
#    println()
    if ismissing(database)
        if verbose; println("Load TB parameters from file"); end
        ManageDatabase.prepare_database(c)
        database = ManageDatabase.database_cached
        println()
    end

    energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbc = SCF.scf_energy(c, database, smearing=smearing, grid = grid, conv_thr = conv_thr, iters = iters, mix = mix,  mixing_mode=mixing_mode, nspin=nspin, e_den0=eden, verbose=verbose)

    conv_flag = !error_flag
    if tbc.within_fit == false
        conv_flag = false
    end
    if conv_flag == true
        println("scf_energy success, done")
    else
        println("scf_energy error detected somewhere!!!!!!!!!!, done")
    end
    println()


    energy_tot = convert_energy(energy_tot)

    return energy_tot, tbc, conv_flag

end


"""
    scf_energy(d::dftout)

    SCF energy using crystal structure from DFT object.
"""
function scf_energy(d::dftout; database = Dict(), smearing=0.01, grid = missing, conv_thr = 1e-5, iters = 75, mix = -1.0, mixing_mode=:pulay, nspin=1, verbose=true)

    return scf_energy(d.crys, smearing=smearing, grid = grid, conv_thr = conv_thr, iters = iters, mix = mix, mixing_mode=mixing_mode, nspin=nspin, verbose=verbose)

end


"""
    scf_energy(tbc::tbc_crys)

    SCF energy using crystal structure from TBC object.
"""
function scf_energy(tbc::tb_crys; smearing=0.01, grid = missing, e_den0 = missing, conv_thr = 1e-5, iters = 75, mix = -1.0, mixing_mode=:pulay, nspin=1, verbose=true)

    energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbc = SCF.scf_energy(tbc; smearing=smearing, grid = grid, e_den0 = e_den0, conv_thr = conv_thr, iters = iters, mix = mix, mixing_mode=mixing_mode, nspin=nspin, verbose=verbose)

    conv_flag = !error_flag
    if tbc.within_fit == false
        conv_flag = false
    end
    if conv_flag == true
        println("success, done")
    else
        println("error detected!!!!!!!!!!, done")
    end
    println()

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


end #end module


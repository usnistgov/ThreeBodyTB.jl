module Relax

using ..SCF:scf_energy
#using ..Force_Stress:get_energy_force_stress_fft
using ..Force_Stress:get_energy_force_stress_fft_LV_sym_SINGLE
using ..Force_Stress:safe_mode_energy
using ..CrystalMod:get_grid
using ..CrystalMod:crystal
using ..CrystalMod:print_with_force_stress
using ..CrystalMod:makecrys
using ..CrystalMod:write_axsf
using ..CalcTB:calc_tb_lowmem2
using ..CalcTB:calc_tb_LV
using ..CalcTB:calc_tb_LV_sparse
using ..Atomdata:atom_radius

using LinearAlgebra
using ..Utility:inv_reshape_vec
using ..Utility:reshape_vec

##using Optim
using LineSearches
using ..ManageDatabase:prepare_database
using ..ManageDatabase:database_cached
using ..MyOptim:conjgrad

using ..Classical:calc_energy_cl
using ..Classical:energy_force_stress_cl

export relax_structure

#using ..ThreeBodyTB:global_length_units
#using ..ThreeBodyTB:global_energy_units
using ..ThreeBodyTB:set_units


"""
    function relax_structure(crys::crystal, database; smearing = 0.01, grid = missing, mode="vc-relax", nsteps=50, update_grid=true, conv_thr=2e-4)

Relax structure. Primary user function is relax_structure in ThreeBodyTB.jl, which calls this one.
"""
function relax_structure(crys::crystal, database; smearing = 0.01, grid = missing, mode="vc-relax", nsteps=50, update_grid=true, conv_thr = 1e-2, energy_conv_thr = 2e-4, filename="t.axsf", nspin=1, repel=true, do_tb=true, database_classical=missing, do_classical=true, tot_charge=0.0, sparse=false)

#    println("relax_structure conv_thr $conv_thr energy_conv_thr (Ryd) $energy_conv_thr ")

    
#    old_units = [deepcopy(global_energy_units), deepcopy(global_length_units)]
    
    old_units = set_units()
    set_units(both="atomic")
#    set_units(both="atomic")

    
    if update_grid==false && ismissing(grid)
        grid = get_grid(crys)
    end

    if (!do_classical && !do_tb) || (!do_tb && ismissing(database_classical))
        println("WARNING, no tb and no classical relax makes no sense")
        return c, missing, 0.0, zeros(c.nat, 3), zeros(3,3)
    end

    

    eden = missing #starts off missing

    error_flag = false
    
    #do this ahead of first iteration, to get memory in correct place
    if do_tb
        if sparse == false
            tbc = calc_tb_LV(deepcopy(crys), database, verbose=false, repel=repel)
        else
            tbc = calc_tb_LV_sparse(deepcopy(crys), database, verbose=false, repel=repel)
        end
        energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=eden, conv_thr=1e-7,nspin=nspin, verbose=false,database_classical=database_classical, do_classical=do_classical, tot_charge=tot_charge )
        eden = deepcopy(tbc.eden)
    else
        tbc = missing
        eden = missing
    end
    
    if error_flag
        println("warning error computing scf in relax_structure, zeroth iteration")
    end



    A0 = deepcopy(crys.A)
    strain = zeros(3,3)
    
    x0 = inv_reshape_vec(crys.coords, strain, crys.nat, strain_mode=true)
    crys_working = deepcopy(crys)

    fcall = 0
    firstiter = true

    nat = crys.nat


    function fix_strain(s)
        for i = 1:3
            for j = 1:3
                s[i,j] = min(s[i,j], 2.0)
                s[i,j] = max(s[i,j], -2.0)
            end
        end
        return 0.5*(s'+s)
    end

    energy_global = -99.0

    CRYSTAL = []
    FORCES = []
    
    function fn(x)
#        println("CALL FN FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF")
#        println(x)
#        println()
        coords, strain = reshape_vec(x, nat, strain_mode=true)

        strain=fix_strain(strain)
        coords = coords .% 1.0

#        println("strain")
#        println(strain)
        
        A = A0 * (I(3) + strain)
        
        crys_working.coords = coords

        
        if  mode == "vc-relax"
            crys_working.A = A
        end

#        println("crys_working")
#        println(crys_working)
#        println()

        if do_tb
            if crys_working == tbc.crys
                #println("SKIP")
                #we already have energy from calling grad, we don't need to call again.
                return energy_global, true
            end
        end
        
        #        println("FN crys")
#        println(crys_working)
#        println("before short ")
        
        tooshort, energy_short = safe_mode_energy(crys_working, database)
        if tooshort
            println("too short $tooshort   $energy_short")
            return energy_short, false
        end



        if do_tb
            if crys_working != tbc.crys
                if sparse == false
                    tbc = calc_tb_LV(deepcopy(crys_working), database, verbose=false, repel=repel)
                else
                    tbc = calc_tb_LV_sparse(deepcopy(crys_working), database, verbose=false, repel=repel)
                end
                tbc.tot_charge=tot_charge
            end
        end

        if do_tb
            energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=eden, verbose=false, conv_thr=5e-7,nspin=nspin,database_classical=database_classical, do_classical=do_classical   )
            eden = deepcopy(tbcx.eden)
        else
            energy_tot, error_flag = calc_energy_cl(crys_working, database=database_classical)

        end           

#        println("fn $energy_tot fnffnfnfnffnfnfnfnfnfnfnfnfnfnfffffff")

#        println("ENERGY $energy_tot $energy_global")
        
        energy_global=energy_tot
#        println("FN end")

        return energy_tot, !error_flag

    end

    f_cart_global = []
    stress_global = []

#    function grad(x)
    function grad(storage, x)

#        println("CALL GRAD", x)

#        println("typeof x ", typeof(x))
        fcall += 1

        coords, strain = reshape_vec(x, nat, strain_mode=true)

        strain=fix_strain(strain)
        coords = coords .% 1.0

        A = A0 * (I(3) + strain)

        crys_working.coords = coords

        if  mode == "vc-relax"
            crys_working.A = A
        end

#        println("GRAD crys")
#        println(crys_working)

        
#        println("crys_working $fcall")
#        println(crys_working)

        tooshort, energy_short = safe_mode_energy(crys_working, database)

#        println("too short ", tooshort)

        if do_tb
            if sparse == false
                tbc = calc_tb_LV(deepcopy(crys_working), database, verbose=false, repel=repel)
            else
                tbc = calc_tb_LV_sparse(deepcopy(crys_working), database, verbose=false, repel=repel)                
            end
            tbc.tot_charge=tot_charge
        end

        #        if crys_working != tbc.crys && !tooshort
        if !tooshort        
#            println("yes calc forces -----------------------------------------------------------------------------------------")
#            println(crys_working)

            if do_tb
                energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=eden, verbose=false, conv_thr=1e-6,nspin=nspin,database_classical=database_classical, do_classical=do_classical    )
                eden = deepcopy(tbcx.eden)

            else
                energy_tot, error_flag  = calc_energy_cl(crys_working, database=database_classical)
            end


            #        else
#            println("no calc forces ------------------------------------------------------------------------------------------")
#            println(crys_working)
#            println(tbc.crys)
#            println("__")
        elseif tooshort
            energy_tot=energy_short
        end
        
        energy_global=energy_tot

        f_cart = zeros(crys_working.nat,3)
        stress = zeros(3,3)
        if do_tb
            if database["SCF"] == true
                energy_tmp,  f_cart_tb, stress_tb =  get_energy_force_stress_fft_LV_sym_SINGLE(tbc, database; do_scf = false, smearing = smearing, grid = grid, vv=[VECTS, VALS, efermi] ,nspin=nspin, repel=repel)
            else
                energy_tmp,  f_cart_tb, stress_tb =  get_energy_force_stress_fft_LV(tbc, database; do_scf = false, smearing = smearing, grid = grid, vv=[VECTS, VALS, efermi] ,nspin=nspin, repel=repel)
            end
            
            f_cart += f_cart_tb
            stress += stress_tb
        end
        if do_classical && !ismissing(database_classical)
            energy_tmp, f_cart_cl, stress_cl = energy_force_stress_cl(crys_working, database=database_classical)
            f_cart += f_cart_cl
            stress += stress_cl
        end
            
        

        #energy_tmp,  f_cart, stress =  get_energy_force_stress_fft(tbc, database; do_scf = true, smearing = smearing, grid = grid )


        push!(CRYSTAL, deepcopy(crys_working)) #correct?
        push!(FORCES, f_cart)
        
        f_cart_global = f_cart
        stress_global = stress
        
        fsum = sum((f_cart).^2)^0.5
        ssum= sum((stress).^2)^0.5

#        println("f_cart")
#        println(f_cart)

        println()
        println("FCALL $fcall en:  $energy_tot (Ryd)  fsum:  $fsum  ssum:  $ssum    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        println()

        print_with_force_stress(crys_working, f_cart_global, stress_global)
        println()
#        println("grad stress: ", stress)
        
#        stress_units = cell_force( crys_working.A, stress)
        #        stress_units = (stress * abs(det(crys_working.A)) *  inv(crys_working.A))'
        #stress_units = (stress * abs(det(crys_working.A)) )'        

#        println("grad stress units: ", -stress_units)

        #        f_crys = f_cart * inv(crys_working.A)

        f_crys = f_cart * crys_working.A'

        g = inv_reshape_vec(-f_crys, -stress* abs(det(crys_working.A))  , crys_working.nat, strain_mode=true)

#        println("g ", g)
        
        if mode != "vc-relax" #don't need stress
            g[3*nat+1:end] = zeros(6)
        end

        storage[:] = g

        neaten(storage)
        
        return energy_global, storage
#        println("ret storage $fcall")
#        println(storage)
    end


    println("starting vec ")
    println(x0)
    println()


#    opts = Optim.Options(g_tol = 7e-4,f_tol = 7e-4, x_tol = 7e-4,
#                         iterations = nsteps,
#                         store_trace = true,
#                         show_trace = false)

    
#    opts = Optim.Options(g_tol = conv_thr,f_tol = conv_thr, x_tol = conv_thr,
#                             iterations = nsteps,
#                             store_trace = true,
#                             show_trace = false)


    #res = optimize(fn,grad, x0, ConjugateGradient(), opts)



    #=
    #preconditioner, guess for inv Hess is based on the  metric. see qe bfgs_module.f90
    function init(x)
        num = 3*nat + 6

#        factor=0.1
        factor = 1.0
        
        P = zeros(eltype(x), num, num)

        A = crys.A
        g = A'*A
        ginv = inv(g)
        
        for a = 1:nat
            aa = (a-1)*3
            for i = 1:3
                for j = 1:3
                    P[aa+i ,aa+j] = g[i,j] * factor
                end
            end
        end

        vol = abs(det(A))
        
#        if mode == "vc-relax"
        for b = 1:2
            bb = 3 * nat + (b - 1)*3
            for i = 1:3
                for j = 1:3
                    P[bb+i,bb+j] = 0.04 * vol * ginv[i,j] * factor
                end
            end
        end
#    end
        
        return inv(P)
    end
=#
    
    res = missing
    #    if false

        #=
        res = optimize(fn,grad, x0, ConjugateGradient( linesearch=LineSearches.BackTracking( maxstep=0.05  ) ) , opts)

        minvec = Optim.minimizer(res)    

        coords, strain = reshape_vec(minvec, nat, strain_mode=true)
        coords= coords .% 1.0

        A = A0 *( I(3) + strain)

        crys = makecrys(A, coords, crys.types, units="Bohr")

        tbc = calc_tb_fast(deepcopy(crys), database)
        energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=eden, conv_thr=1e-7)

        eden = deepcopy(tbc.eden)

        A0 = deepcopy(crys.A)
        strain = zeros(3,3)
        
        x0 = inv_reshape_vec(crys.coords, strain, crys.nat, strain_mode=true)
        crys_working = deepcopy(crys)

        fcall = 1000
        firstiter = true

        minvec = Optim.minimizer(res)    
        =#
        ##    else

#    println("fn ", fn(x0))
    

    minvec, energy, grad = conjgrad(fn, grad, x0; maxstep=1.0, niters=nsteps, conv_thr = conv_thr, fn_conv = energy_conv_thr, verbosity="low")

##    end


#    println("minvec")
#    println(minvec)

    coords, strain = reshape_vec(minvec, nat, strain_mode=true)
    coords = coords .% 1.0
    
    strain = fix_strain(strain)


    A = A0 *( I(3) + strain)


    cfinal = makecrys(A, coords, crys.types, units="Bohr")


    if !ismissing(filename)
        write_axsf(CRYSTAL, filename=filename, FORCES=FORCES)
    end

    

    set_units(energy=old_units[1], length=old_units[2])
    #global_energy_units = old_units[1]
    #global_length_units = old_units[2]
    
    return cfinal, tbc, energy, f_cart_global, stress_global

#    return res

end

"""
    function neaten(storage)

Detect/Enforce symmetries in forces/stress, tight tolerance. Deal with minor numerical issues causing symmetry breaking
"""
function neaten(storage)
#    println("n before ", storage)
    n = length(storage)
    for i in 1:n-6
        if abs(storage[i]) < 1e-8
            storage[i] = 0.0
        end
    end
        
    for i in n-2:n
        if abs(storage[i]) < 1e-7
            storage[i] = 0.0
        end
    end

    for i in 1:n
        for j in i+1:n
            if abs(storage[j] - storage[i]) < 1e-8
                storage[j] = storage[i]
            elseif abs(storage[j] + storage[i]) < 1e-8
                storage[j] = -storage[i]
            elseif abs(storage[j] + 2.0*storage[i]) < 1e-8
                storage[j] = -2.0*storage[i]
            elseif abs(storage[j] + 0.5*storage[i]) < 1e-8
                storage[j] = -0.5*storage[i]
            elseif abs(storage[j] - 2.0*storage[i]) < 1e-8
                storage[j] = 2.0*storage[i]
            elseif abs(storage[j] - 0.5*storage[i]) < 1e-8
                storage[j] = 0.5*storage[i]
            end
        end
    end
#    println("n after ", storage)
end


#=
function finite_diff_grad(fn, x0, step)
    x = deepcopy(x0)
    g = zeros(size(x))
    for i in 1:length(x)
        x[:] = x0[:]
        x[i] += step
        a = fn(x)

        x[:] = x0[:]
        x[i] -= step
        b = fn(x)

        g[i] = (a-b)/(2*step)
    end
    
    return g

end
=#


function make_random_crys(types; database=missing)
    return make_random_crystal(types; database=database)
end

"
    function make_random_crystal(types)

Make a random crystal with types t
"
function make_random_crystal(types; database=missing)

    nat = length(types)

    db = database
    if ismissing(database)
        prepare_database(types)
        db = database_cached
    end

    total_vol = 0.0
    for t in types
        total_vol += 4 * pi / 3 * (atom_radius[t]/100 / 0.529177)^3 * 1.0
    end
    println("total_vol $total_vol")
    
    for i = 1:200

        A = zeros(3,3)
        c = rand(nat, 3)
        c[1,:] = [0.0 0.0 0.0]
        
        A[1,1] = 1.0
        
        l2 = 1.0 + (rand(1)[1] - 0.5)
        th2 =  pi*(1.0 + 0.5*(rand(1)[1] - 0.5))
        
        l3 = 1.0 + (rand(1)[1] - 0.5)
        th3 =  pi*(1.0 + 0.5*(rand(1)[1] - 0.5))

        phi = pi*(0.25*(rand(1)[1] - 0.5))
        
        A[2,:] = [l2 * sin(th2)  l2 * cos(th2) 0.0]
        A[3,:] = [sin(phi)  l3 * sin(th3) *cos(phi)  l3 * cos(th3)* cos(phi)]
        
#        println("A before ", det(A), "  ", total_vol^(1/3) * 1.5 )
        A = A / (abs(det(A)))^(1/3) * total_vol^(1/3) * 1.2
#        println("A after ", det(A))
        
        A = A * det(A) / abs(det(A))
        
        crys = makecrys(A, c, types, units="Bohr")
        
        within_fit, rsum = calc_tb_LV(crys*0.97, db, check_only=true)
        if within_fit
            return crys * 1.03
        end

        total_vol = total_vol * 1.04
        
    end
    println("failed to generate crystal")
    return missing
        
end


end #end module

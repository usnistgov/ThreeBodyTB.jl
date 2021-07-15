###include("Crystal.jl")
###include("DFToutMod.jl")
#using XMLDict



"""
    module Force_Stress

Module for calculating force and stress
"""
module Force_Stress
"""
Scripts to calculate force and stress
"""

#calc for testing non-autodiff forces only 
#using Calculus

using Base.Threads
using FFTW
using LinearAlgebra
using ForwardDiff
#using ReverseDiff
###using Optim
using LineSearches
using ..CrystalMod:crystal
using ..CrystalMod:makecrys


using ..CalcTB:calc_tb_fast
using ..CalcTB:distances_etc_3bdy_parallel
using ..TB:calc_energy_charge_fft
using ..TB:tb_crys
using ..TB:types_energy
using ..TB:make_kgrid
using ..TB:get_dq
using ..TB:get_h1
using ..TB:ewald_energy
using ..CalcTB:calc_frontier

using ..Ewald:electrostatics_getgamma
using ..Ewald:estimate_best_kappa
using ..SCF:scf_energy

using ..BandTools:gaussian
using ..CrystalMod:get_grid

export get_energy_force_stress
#export relax_structure


"""
    function get_energy_force_stress(crys::crystal, database; smearing = 0.01, grid = missing)

Get force and stress, non-fft algorithm. Generally use the fft algorithm.

`return energy_tot,  f_cart, stress`

Returns Ryd units. Generally users should use the `scf_energy_force_stress` function

Uses automatic differentation for gradient.
"""
function get_energy_force_stress(crys::crystal, database; smearing = 0.01, grid = missing)

    println("crys")
#    tbc = []
    tbc = calc_tb_fast(crys, database)
    return get_energy_force_stress_fft(tbc, database, do_scf=tbc.scf, grid = grid, smearing=smearing)
end

"""
    function get_energy_force_stress(crys::crystal, database; smearing = 0.01, grid = missing)

Get force and stress, non-fft algorithm
"""
function get_energy_force_stress(tbc::tb_crys, database; do_scf=false, smearing = 0.01, grid = missing, e_den0=missing, vv = missing, cs = 4)

    if ismissing(grid)
        grid = get_grid(tbc.crys)
    end
    
    kgrid, kweights = make_kgrid(grid)
    nk = size(kgrid)[1]

    tooshort, energy_tot = safe_mode_energy(tbc.crys, database)

    if !(tooshort)
        if !ismissing(vv)
            VECTS, VALS, efermi = vv
            energy_tot = 0.0
        else
        
            #prepare eigenvectors / values
            error_flag = false
            if do_scf
                energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=e_den0, conv_thr = 1e-9)
            else
                energy_tot, efermi, e_den, VECTS, VALS, error_flag =  calc_energy_charge_fft(tbc, grid=grid, smearing=smearing)
            end
            if error_flag
                println("warning, trouble with eigenvectors/vals in initial step get_energy_force_stress")
            end
            
        end

        h1, dq = get_h1(tbc)
        
#        println("energy_tot $energy_tot")

        OCCS = gaussian.(VALS.-efermi, smearing)
    
    end


    ct = deepcopy(tbc.crys)
#=
    #reshape function
    function reshape_vec(x, nat)
        T=typeof(x[1])

        print("size x ", size(x))
        x_r = zeros(T, ct.nat, 3)
        for n = 1:ct.nat
            for j = 1:3
                x_r[n,j] = x[3*(n-1) + j]
            end
        end

        x_r_strain = zeros(T, 3,3)


        x_r_strain[1,1] = x[3*nat+1]
        x_r_strain[2,2] = x[3*nat+2]
        x_r_strain[3,3] = x[3*nat+3]

        x_r_strain[2,3] = x[3*nat+4]
        x_r_strain[3,2] = x[3*nat+4]

        x_r_strain[1,3] = x[3*nat+5]
        x_r_strain[3,1] = x[3*nat+5]

        x_r_strain[1,2] = x[3*nat+6]
        x_r_strain[2,1] = x[3*nat+6]

        return x_r, x_r_strain
    end

    #reshape function
    function inv_reshape_vec(x, strain, nat)
        T=typeof(x[1])
        x_r = zeros(T, nat*3 + 6)
        for n = 1:nat
            for j = 1:3
                x_r[3*(n-1) + j] = x[n,j] 
            end
        end

        x_r[3*ct.nat + 1] = strain[1,1]
        x_r[3*ct.nat + 2] = strain[2,2]
        x_r[3*ct.nat + 3] = strain[3,3]
        x_r[3*ct.nat + 4] = strain[2,3]+strain[3,2]
        x_r[3*ct.nat + 5] = strain[1,3]+strain[3,1]
        x_r[3*ct.nat + 6] = strain[1,2]+strain[2,1]

        return x_r
    end
=#

    function f(x::Vector)
        T=typeof(x[1])

        x_r, x_r_strain = reshape_vec(x, ct.nat, strain_mode=true)
        
        A = ct.A * (I(3) + x_r_strain)
        #A = deepcopy(ct.A)

        crys_dual = makecrys( A , ct.coords + x_r, ct.types, units="Bohr")


        #this deals with cases where the distances between atoms become very short, which can happen during relaxations
        #currently we just have an artificial repulsive force in this case
        if tooshort
            tooshort, energy_short = safe_mode_energy(crys_dual, database, var_type=T)
            return energy_short
        end


        
        if database["scf"] == true
            scf = true
            kappa = estimate_best_kappa(ct.A)
            gamma_dual = electrostatics_getgamma(crys_dual, kappa=kappa)
        else
            scf = false
            gamma_dual=zeros(T, ct.nat,ct.nat)
        end
        

        tbc_dual = calc_tb_fast(crys_dual, database; verbose=false, var_type=T, use_threebody=true, use_threebody_onsite=true, gamma=gamma_dual, check_frontier=false)

        nwan = tbc.tb.nwan

        hk = zeros(Complex{T}, nwan, nwan)
        hk0 = zeros(Complex{T}, nwan, nwan)
        hka = zeros(Complex{T}, nwan, nwan)
        sk = zeros(Complex{T}, nwan, nwan)

        twopi_i = -1.0im*2.0*pi

        VALS0 = zeros(T, nk, nwan)


        #analytic fourier transform
        for k = 1:nk
            #vect, vals, hk, sk, vals0 = Hk(hktemp, sktemp, tbc.tb, grid[k,:])


            hk0[:,:] .= 0.0
            sk[:,:] .= 0.0
            
            kmat = repeat(kgrid[k,:]', tbc.tb.nr,1)
            exp_ikr = exp.(twopi_i * sum(kmat .* tbc.tb.ind_arr,dims=2))

            for m in 1:tbc.tb.nwan
                 for n in 1:tbc.tb.nwan
                     hk0[m,n] = tbc_dual.tb.H[m,n,:]'*exp_ikr[:]
                     sk[m,n]  = tbc_dual.tb.S[m,n,:]'*exp_ikr[:]
                 end
            end
            hk0 = 0.5*(hk0 + hk0')
            sk = 0.5*(sk + sk')


            
            if scf
                hk0 = hk0 + h1 .* sk
            end

            for a = 1:nwan
                
                hka[:,:] = hk0 - ( VALS[k,a] )  * sk
                #hka[:,:] = hk0 - ( VALS[k,a] )  * sk                

                VALS0[k,a] += real.(VECTS[k,:,a]'*hka*VECTS[k,:,a])
            end
        end            


        energy0 = sum(OCCS .* VALS0) / nk * 2.0

        if scf
            eewald, pot = ewald_energy(crys_dual, gamma_dual, dq)
        else
            eewald = 0.0
        end


        etypes = types_energy(tbc.crys)
        

        return energy0 + etypes + eewald


    end

#    x0 = inv_reshape_vec(ct.coords, ct.nat)

    chunksize=min(cs, 3*ct.nat + 6)
    cfg = ForwardDiff.GradientConfig(f, zeros(3*ct.nat + 6), ForwardDiff.Chunk{chunksize}())
    g = ForwardDiff.gradient(f, zeros(3*ct.nat + 6)  )

    x, stress = reshape_vec(g, ct.nat)

    f_cart = -1.0 * x
    f_cart = f_cart * inv(ct.A)'

    stress = -stress / abs(det(ct.A))

    for i = 1:3
        for j = 1:3
            if abs(stress[i,j]) < 1e-12
                stress[i,j] = 0.0
            end
        end
    end

    return energy_tot,  f_cart, stress


end


#primarily for testing
"""
    function finite_diff(crys::crystal, database, ind1, ind2; stress_mode=false, step = 0.0002, smearing = 0.01, grid = missing)

Finite differences force/stress, for testing.
# Arguments
- `crys::crystal` Crystal structure
- `database` Database of fitting coefficents.
- `ind1`  atom index for first stress index
- `ind2` cartesian index or second stress index
- `stress_mode=false` true for stress, otherwise force.
- `step = 0.0002` step_size for finite steps.
- `smearing = 0.01` smearing energy
- `grid = missing` kpoint grid
"""
function finite_diff(crys::crystal, database, ind1, ind2; stress_mode=false, step = 0.0002, smearing = 0.01, grid = missing)
    if ismissing(grid)
        grid = get_grid(crys)
    end


    tbc0 = calc_tb_fast(crys, database, verbose=false, check_frontier=false)

    energy_tot0, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc0, smearing=smearing, grid=grid, conv_thr = 1e-9)


    if stress_mode == false
        
        println("force mode")

        crys1 = deepcopy(crys)
        cart1 = crys1.coords * crys1.A
        cart1[ind1, ind2] += step

        crys1.coords = cart1 * inv(crys1.A)

        tbc1 = calc_tb_fast(crys1, database, verbose=false, check_frontier=false)

        energy_tot1, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc1, smearing=smearing, grid=grid, conv_thr = 1e-9)


        crys2 = deepcopy(crys)
        cart2 = crys2.coords * crys2.A
        cart2[ind1, ind2] -= step

        crys2.coords = cart2 * inv(crys2.A)

        tbc2 = calc_tb_fast(crys2, database, verbose=false, check_frontier=false)

        energy_tot2, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc2, smearing=smearing, grid=grid, conv_thr = 1e-9)
        
        force = - (energy_tot1 - energy_tot2) / (2 * step)

        return energy_tot0, force


    else
        println("stress mode")

        crys1 = deepcopy(crys)

        strain = zeros(3,3)
        strain[ind1,ind2] = step
        strain[ind2,ind1] = step
       
        crys1.A = crys1.A *(I(3) + strain)

        tbc1 = calc_tb_fast(crys1, database, verbose=false, check_frontier=false)

        energy_tot1, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc1, smearing=smearing, grid=grid, conv_thr = 1e-7)


        crys2 = deepcopy(crys)
        strain = zeros(3,3)
        strain[ind1,ind2] = -step
        strain[ind2,ind1] = -step
       
        crys2.A = crys2.A *(I(3) + strain)

        tbc2 = calc_tb_fast(crys2, database, verbose=false, check_frontier=false)

        energy_tot2, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc2, smearing=smearing, grid=grid, conv_thr = 1e-7)

        stress = -1.0* (energy_tot1 - energy_tot2) / (2 * step) / abs(det(crys.A))

        if ind1 != ind2
            stress = stress / 2.0
        end

        return energy_tot0, stress
    end
end


function cell_force(A, stress)

    stress = (stress + stress')/2.0

#    println("cf stress ", stress)
#    println("cf A      ", A)

    Ainv = inv(A)
    cf = zeros(3,3)
    for j= 1:3
        for i = 1:3
            cf[i,j] = sum(Ainv[j,:].*stress[i,:])
        end
    end

    return cf * abs(det(A))

end



######################################################
#=
"""
    function relax_structure(crys::crystal, database; smearing = 0.01, grid = missing, mode="vc-relax", nsteps=100, update_grid=true, conv_thr=2e-4)

Relax structure. Primary user function is relax_structure in ThreeBodyTB.jl, which calls this one.
"""
function relax_structure(crys::crystal, database; smearing = 0.01, grid = missing, mode="vc-relax", nsteps=100, update_grid=true, conv_thr = 2e-4)

    if update_grid==false
        grid = get_grid(crys)
    end

    eden = missing #starts off missing

    #do this ahead of first iteration, to get memory in correct place
    tbc = calc_tb_fast(deepcopy(crys), database)
    energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=eden, conv_thr=1e-7)

    if error_flag
        println("warning error computing scf in relax_structure, zeroth iteration")
    end

    eden = deepcopy(tbc.eden)

    x0 = inv_reshape_vec(crys.coords, crys.A, crys.nat, strain_mode=false)
    crys_working = deepcopy(crys)

    fcall = 0
    firstiter = true

    nat = crys.nat

    function fn(x)
        println("CALL FN", x)
        coords, A = reshape_vec(x, nat)

#        println("coords ", coords)
#        println("A ", A)

        crys_working.coords = coords
        if  mode == "vc-relax"
            crys_working.A = A
        end
        println("FN crys")
        println(crys_working)
        
        tooshort, energy_short = safe_mode_energy(crys_working, database)
#        println("fn too short ", tooshort)
        if tooshort
            return energy_short
        end


        if firstiter
            firstiter = false
        else
            if crys_working != tbc.crys
#                println("yes calc xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
#                println(crys_working)
                tbc = calc_tb_fast(deepcopy(crys_working), database, verbose=false)
            else
#                println("nocalc xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
#                println(crys_working)
#                println(tbc.crys)
                
            end
#            println(crys_working)
            energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=eden, verbose=false, conv_thr=1e-6)
            eden = deepcopy(tbcx.eden)
        end
#        println("fn $energy_tot fnffnfnfnffnfnfnfnfnfnfnfnfnfnfffffff")
        return energy_tot

    end

    f_cart_global = []
    stress_global = []
    
    function grad(storage, x)
        println("CALL GRAD", x)

#        println("typeof x ", typeof(x))
        fcall += 1

        coords, A = reshape_vec(x, nat)
        crys_working.coords = coords

        if  mode == "vc-relax"
            crys_working.A = A
        end

        println("GRAD crys")
        println(crys_working)

        
#        println("crys_working $fcall")
#        println(crys_working)

        tooshort, energy_short = safe_mode_energy(crys_working, database)
        
#        println("too short ", tooshort)

        if crys_working != tbc.crys && !tooshort
#            println("yes calc forces -----------------------------------------------------------------------------------------")
#            println(crys_working)
            tbc = calc_tb_fast(deepcopy(crys_working), database, verbose=false)
            energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=eden, verbose=false, conv_thr=1e-6)

            eden = deepcopy(tbcx.eden)

            #        else
#            println("no calc forces ------------------------------------------------------------------------------------------")
#            println(crys_working)
#            println(tbc.crys)
#            println("__")
        end


        energy_tmp,  f_cart, stress =  get_energy_force_stress_fft(tbc, database; do_scf = false, smearing = smearing, grid = grid, vv=[VECTS, VALS, efermi] )
        #energy_tmp,  f_cart, stress =  get_energy_force_stress_fft(tbc, database; do_scf = true, smearing = smearing, grid = grid )

        f_cart_global = f_cart
        stress_global = stress
        
        fsum = sum((f_cart).^2)^0.5
        ssum= sum((stress).^2)^0.5

#        println("f_cart")
#        println(f_cart)
        
        println("FCALL $fcall en:  $energy_tot (Ryd)  fsum:  $fsum  ssum:  $ssum    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")

#        println("grad stress: ", stress)
        
        stress_units = cell_force( crys_working.A, stress)
        #        stress_units = (stress * abs(det(crys_working.A)) *  inv(crys_working.A))'
        #stress_units = (stress * abs(det(crys_working.A)) )'        

#        println("grad stress units: ", -stress_units)

        #        f_crys = f_cart * inv(crys_working.A)
        f_crys = f_cart * crys_working.A

        g = inv_reshape_vec(-f_crys, -stress_units, crys_working.nat, strain_mode=false)

#        println("g ", g)
        
        if mode != "vc-relax" #don't need stress
            g[3*nat+1:end] = zeros(9)
        end
        storage[:] = g
#        println("ret storage $fcall")
#        println(storage)
    end

#=
    function fn_fake(x)
        println("fn_fake ", x[1]*2)
        return sum(x.^2)
    end

    function grad_fake(storage, x)
        storage[:] =  2.0 * x

        println("grad_fake  ", x[1]*2)
        println("grad_fakex ", x)
        println("grad_fakes ", storage)
        return storage

    end
=#


#    opts = Optim.Options(g_tol = 1e-3,f_tol = 1e-3, x_tol = 1e-3,
#                         iterations = 10,
#                         store_trace = true,
#                         show_trace = false)


#    println("X000000000000000000, ", x0)
#    println()
#    println(fn(x0))
#    println("X000000000000000001, ", x0)
#    storage = zeros(size(x0))
#    println(grad(x0, storage))

#    println("storage")
#    println(storage)


    

    println("starting vec ")
    println(x0)
    println()


#    opts = Optim.Options(g_tol = 7e-4,f_tol = 7e-4, x_tol = 7e-4,
#                         iterations = nsteps,
#                         store_trace = true,
#                         show_trace = false)

    
    opts = Optim.Options(g_tol = conv_thr,f_tol = conv_thr, x_tol = conv_thr,
                             iterations = nsteps,
                             store_trace = true,
                             show_trace = false)


    #res = optimize(fn,grad, x0, ConjugateGradient(), opts)


    #preconditioner, guess for inv Hess is based on the  metric. see qe bfgs_module.f90
    function init(x)
        num = 3*nat + 9

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
        for b = 1:3
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

#    P = init(x0)
#    for i = 1:size(P)[1]
#        for j = 1:size(P)[2]
#            println("P $i $j ", P[i,j])
#        end
#    end
#    println("P")
#    println(P)

    #x -> Matrix{eltype(x)}(I, length(x), length(x))

#    res = optimize(fn,x0, NelderMead(), opts)

    
    
    #    res = optimize(fn,grad, x0, ConjugateGradient())
    
#    res = optimize(fn,grad, x0, BFGS( ), opts)

    #res = optimize(fn,grad, x0, BFGS(  initial_invH = init ), opts)

    res = missing
    
    try    
        res = optimize(fn,grad, x0, BFGS(  initial_invH = init, linesearch=LineSearches.MoreThuente() ), opts)

        #res = optimize(fn,grad, x0, BFGS(initial_invH = init ), opts)
        
    catch e
        if e isa InterruptException
            println("user interrupt")

            return tbc.crys
        else
            println("unknown error")
            return tbc.crys
            
        end
    end
    
    #res = optimize(fn,grad, x0, BFGS(  initial_invH = init, linesearch=LineSearches.HagerZhang() ), opts)    

    energy = res.minimum
    
    # bad : res = optimize(fn,grad, x0, ConjugateGradient(linesearch=LineSearches.BackTracking(order=3) ), opts)    
    #linesearch=LineSearches.BackTracking(order=2),
    #LineSearches.MoreThuente()
    #LineSearches.HagerZhang()
    #linesearch=LineSearches.BackTracking(order=2)
    #linesearch=LineSearches.StrongWolfe()
    
#for testing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#    storage = zeros(size(x0))
#    g  = grad(storage, x0)

#    println()

#    res = Calculus.gradient(fn, x0)

#    println("grad x0")
#    println(storage)
#    println("calc")
#    println(res)

#    println([storage res])
#
#    return [storage res]
#end for testing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#    res = optimize(fn, x0, LBFGS())
#
#    res = optimize(fn,grad, x0, LBFGS())

#    res = optimize(fn, x0, ConjugateGradient())

#    res = optimize(fn,grad, x0, LBFGS(), opts)

#    res = optimize(fn, x0, LBFGS(), opts)

#    f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
#    function g!(storage, x)
 #       println("g $x");
 #       storage[1] = -2.0 * (1.0 - x[1]) - 400.0 * (x[2] - x[1]^2) * x[1]
 #       storage[2] = 200.0 * (x[2] - x[1]^2)
 #       println("g storage", storage)
 #   end

#    res = optimize(f, g!,  [0.0, 0.0], LBFGS())

    println()
    println("res")
    println(res)

    minvec = Optim.minimizer(res)    
#    println("minvec")
#    println(minvec)

    coords, A = reshape_vec(minvec, nat)
    cfinal = makecrys(A, coords, crys.types, units="Bohr")
    return cfinal, tbc, energy, f_cart_global, stress_global

#    return res

end
=#

"""
    function reshape_vec(x, nat; strain_mode=false)

The force and relax algorithms from outside codes take in vectors, not `crystal` 's. So we have to reshape
vectors to and from `crystals`
"""    
function reshape_vec(x, nat; strain_mode=false)
#    println("RESHAPEVEC RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRr ", strain_mode)

    T=typeof(x[1])
    
#    println("size x ", size(x))
    x_r = zeros(T, nat, 3)
    for n = 1:nat
        for j = 1:3
            x_r[n,j] = x[3*(n-1) + j]
        end
    end
    
    x_r_strain = zeros(T, 3,3)

    if length(x) == 3*nat+6 && strain_mode
        x_r_strain[1,1] = x[3*nat+1]
        x_r_strain[2,2] = x[3*nat+2]
        x_r_strain[3,3] = x[3*nat+3]

        x_r_strain[2,3] = 0.5*x[3*nat+4]
        x_r_strain[3,2] = 0.5*x[3*nat+4]

        x_r_strain[1,3] = 0.5*x[3*nat+5]
        x_r_strain[3,1] = 0.5*x[3*nat+5]

        x_r_strain[1,2] = 0.5*x[3*nat+6]
        x_r_strain[2,1] = 0.5*x[3*nat+6]
    elseif length(x) == 3*nat+6 
        x_r_strain[1,1] = x[3*nat+1]
        x_r_strain[2,2] = x[3*nat+2]
        x_r_strain[3,3] = x[3*nat+3]

        x_r_strain[2,3] = x[3*nat+4]
        x_r_strain[3,2] = x[3*nat+4]

        x_r_strain[1,3] = x[3*nat+5]
        x_r_strain[3,1] = x[3*nat+5]

        x_r_strain[1,2] = x[3*nat+6]
        x_r_strain[2,1] = x[3*nat+6]

    elseif length(x) == 3*nat+9
        x_r_strain[1,1] = x[3*nat+1]
        x_r_strain[1,2] = x[3*nat+2]
        x_r_strain[1,3] = x[3*nat+3]
        x_r_strain[2,1] = x[3*nat+4]
        x_r_strain[2,2] = x[3*nat+5]
        x_r_strain[2,3] = x[3*nat+6]
        x_r_strain[3,1] = x[3*nat+7]
        x_r_strain[3,2] = x[3*nat+8]
        x_r_strain[3,3] = x[3*nat+9]
    else
        println("I'm confusing about the length reshape_vec $nat ", length(x) )
    end

    return x_r, x_r_strain
end

"""
    function inv_reshape_vec(x, strain, nat; strain_mode=true)

The force and relax algorithms from outside codes take in vectors, not `crystal` 's. So we have to reshape
vectors to and from `crystals`
"""
function inv_reshape_vec(x, strain, nat; strain_mode=true)
    T=typeof(x[1])
    if strain_mode
        x_r = zeros(T, nat*3 + 6)
    else
        x_r = zeros(T, nat*3 + 9)
    end

    for n = 1:nat
        for j = 1:3
            x_r[3*(n-1) + j] = x[n,j] 
        end
    end
    if strain_mode
        x_r[3*nat + 1] = strain[1,1]
        x_r[3*nat + 2] = strain[2,2]
        x_r[3*nat + 3] = strain[3,3]
        x_r[3*nat + 4] = (strain[2,3]+strain[3,2])
        x_r[3*nat + 5] = (strain[1,3]+strain[3,1])
        x_r[3*nat + 6] = (strain[1,2]+strain[2,1])
    else
        x_r[3*nat+1] = strain[1,1]  
        x_r[3*nat+2] = strain[1,2]  
        x_r[3*nat+3] = strain[1,3]  
        x_r[3*nat+4] = strain[2,1]  
        x_r[3*nat+5] = strain[2,2]  
        x_r[3*nat+6] = strain[2,3]  
        x_r[3*nat+7] = strain[3,1]  
        x_r[3*nat+8] = strain[3,2]  
        x_r[3*nat+9] = strain[3,3]  
    end

    return x_r
end

"""
    function safe_mode_energy(crys::crystal, database; var_type=Float64)

Relaxation can accidently lead to very small atom-atom distances during the relaxation precedure if too large of step is taken. This function is a repulsive energy
function at short range to make sure the relaxtion doesn't get stuck at very short distances
where the fitting doesn't apply.
"""
function safe_mode_energy(crys::crystal, database; var_type=Float64, check=true)

    diststuff = distances_etc_3bdy_parallel(crys,10.0, 0.0, var_type=var_type)
    R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = diststuff
    nkeep = size(R_keep_ab)[1]
    
    energy = 1.0
    tooshort = false

    warned = zeros(Bool, crys.nat)
    
    for a1 = 1:crys.nat
        t1 = crys.stypes[a1]
        for a2 = 1:crys.nat
            t2 = crys.stypes[a2]
            dmin = database[(t1,t2)].min_dist * 0.979
            for c = 1:nkeep
                cind = R_keep_ab[c,1]
                if dist_arr[a1,a2,cind,1] < dmin*2 && dist_arr[a1,a2,cind,1] > 1e-7
                    energy += 1e-3/dist_arr[a1,a2,cind,1]
                end

                if dist_arr[a1,a2,cind,1] < dmin*1.01999 && dist_arr[a1,a2,cind,1] > 1e-7
                    tooshort = true
                    energy += 0.02 * (dist_arr[a1,a2,cind,1] - dmin)^2 + 0.3 * abs(dist_arr[a1,a2,cind,1] - dmin)
                    if var_type == Float64 && warned[a1] == false
                        println("WARNING, SAFE MODE $a1 $t1 $a2 $t2 $c ", dist_arr[a1,a2,cind,1])
                        warned[a1] = true
                    end
                end
            end
        end
    end
    if tooshort
        println("WARNING, safe mode activated, minimum distances < fitting data * 0.98")
        return tooshort, energy
    elseif check==true
        #violation_list, vio_bool = calc_frontier(crys, database, test_frontier=true, diststuff=diststuff, verbose=false)
        violation_list, vio_bool = calc_frontier(crys, database, test_frontier=true, verbose=false)
#        println("vio_vool $vio_bool")
        if vio_bool == false
            println("ACTIVATE SAFE MODE")
            return true, energy
        else
            return false, energy
        end
    else
        return tooshort, energy
    end

        

    
end

##############################################################################################################
"""
    function get_energy_force_stress_fft(tbc::tb_crys, database; do_scf=false, smearing = 0.01, grid = missing, e_den0=missing, vv = missing)

Calculate energy/force/stress using fft algorithm. Users should use `scf_energy_force_stress`, which calls this. Uses automatic differentation for jacobian.
"""
function get_energy_force_stress_fft(tbc::tb_crys, database; do_scf=false, smearing = 0.01, grid = missing, e_den0=missing, vv = missing)

    println("get_energy_force_stress_fft")

    if ismissing(grid)
        grid = get_grid(tbc.crys)
    end
    
    kgrid, kweights = make_kgrid(grid)
    nk = size(kgrid)[1]

    ct = deepcopy(tbc.crys)
    
#    println("test safe get_energy_force_stress_fft")
    tooshort, energy_tot = safe_mode_energy(tbc.crys, database)
    
    if tooshort ##########################
        println("safemode")
        function f(x :: Vector)
            T=typeof(x[1])
            x_r, x_r_strain = reshape_vec(x, ct.nat, strain_mode=true)
            A = ct.A * (I(3) + x_r_strain)
            crys_dual = makecrys( A , ct.coords + x_r, ct.types, units="Bohr")
            
            tooshort, energy_short = safe_mode_energy(crys_dual, database, var_type=T, check=false)
            return energy_short
        end
        garr = ForwardDiff.gradient(f, zeros(3*ct.nat + 6) )


    else #not too short ##################
        
            
        
        

        
        scf = database["scf"]
        #    println("scf ", scf)
        if !(tooshort)
            if !ismissing(vv)
                VECTS, VALS, efermi = vv
                energy_tot = 0.0
            else
                #prepare eigenvectors / values
                error_flag = false
                if do_scf
                    println("doing scf aaaaaaaaaaa")
                    energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=e_den0, conv_thr = 1e-8)
                else
                    energy_tot, efermi, e_den, VECTS, VALS, error_flag =  calc_energy_charge_fft(tbc, grid=grid, smearing=smearing)
                end
                if error_flag
                    println("warning, trouble with eigenvectors/vals in initial step get_energy_force_stress")
                end
            end
            h1, dq = get_h1(tbc)
            OCCS = gaussian.(VALS.-efermi, smearing)
        end




        size_ret = tbc.tb.nwan * tbc.tb.nwan * tbc.tb.nr 
        
        function ham(x :: Vector)
            T=typeof(x[1])
            

            x_r, x_r_strain = reshape_vec(x, ct.nat, strain_mode=true)

            #x_r, x_r_strain = reshape_vec(x, 0, strain_mode=true)

            A = ct.A * (I(3) + x_r_strain)
            crys_dual = makecrys( A , ct.coords + x_r, ct.types, units="Bohr")

            #crys_dual = makecrys( A , ct.coords , ct.types)

            #        gamma_dual=zeros(T, ct.nat,ct.nat)
            #println("gamma")
            if database["scf"] == true
                scf = true
                kappa = estimate_best_kappa(ct.A)
                gamma_dual = electrostatics_getgamma(crys_dual, kappa=kappa)
            else
                scf = false
                gamma_dual=zeros(T, ct.nat,ct.nat)
            end

            tbc_dual = calc_tb_fast(crys_dual, database; verbose=false, var_type=T, use_threebody=true, use_threebody_onsite=true, gamma=gamma_dual , check_frontier=false)
            ret = zeros(T, size_ret * 2 + 1)

            
            ret[1:size_ret] = real.(tbc_dual.tb.H[:])
            ret[size_ret+1:size_ret*2] = real.(tbc_dual.tb.S[:])


            #we stick eewald in the last entry of ret
            if scf
                eewald, pot = ewald_energy(crys_dual, gamma_dual, dq)
                ret[end] = eewald
            end

            return ret
        end

#        println("jac")
        @time begin

            ret = ham(zeros(3*ct.nat + 6))

            #ret = ham(zeros( 6))

            chunksize=min(9, 3*ct.nat + 6)
            
            cfg = ForwardDiff.JacobianConfig(ham, zeros(3*ct.nat + 6), ForwardDiff.Chunk{chunksize}())
            #cfg = ForwardDiff.JacobianConfig(ham, zeros( 6), ForwardDiff.Chunk{chunksize}())
            g = ForwardDiff.jacobian(ham, zeros(3*ct.nat + 6) , cfg ) ::  Array{Float64,2}
#            g = ForwardDiff.jacobian(ham, zeros( 6) , cfg ) ::  Array{Float64,2}

        end
#        println("end jac")
        #    function f_es(x::Vector)
        #        x_r, x_r_strain = reshape_vec(x, ct.nat, strain_mode=true)
        #        A = ct.A * (I(3) + x_r_strain)
        #
        #        crys_dual = makecrys( A , ct.coords + x_r, ct.types)
        #
        #
        #        
        #        if scf
        #            kappa = estimate_best_kappa(ct.A)
        #            gamma_dual = electrostatics_getgamma(crys_dual, kappa=kappa)
        #            eewald, pot = ewald_energy(crys_dual, gamma_dual, dq)
        #        else
        #            eewald = 0.0
        #        end
        #        return eewald
        #    end

        #    cfg_grad = ForwardDiff.GradientConfig(f_es, zeros(3*ct.nat + 6), ForwardDiff.Chunk{3*ct.nat + 6}())
        #    g_es = ForwardDiff.gradient(f_es, zeros(3*ct.nat + 6) , cfg_grad )

        
        
        #    @time gr = ReverseDiff.jacobian(ham, zeros(3*ct.nat + 6)  )

        #    return 0,0,0
        
#        ham_orig = ham(zeros(3*ct.nat + 6))

        #    println("size of g ", size(g))
        #    println("eltype ", eltype(g))

        hr_g = zeros(Complex{eltype(g)},  tbc.tb.nwan, tbc.tb.nwan,  grid[1], grid[2], grid[3],3*ct.nat + 6 )
        sr_g = zeros(Complex{eltype(g)},  tbc.tb.nwan, tbc.tb.nwan,  grid[1], grid[2], grid[3],3*ct.nat + 6 )

        #    hr = zeros(Complex{eltype(g)}, tbc.tb.nwan, tbc.tb.nwan,  grid[1], grid[2], grid[3])
        
        ind = zeros(Int64, 3)
        new_ind = zeros(Int64, 3)

        for c in 1:size(tbc.tb.ind_arr)[1]

            ind[:] = tbc.tb.ind_arr[c,:]
            
            new_ind[1] = mod(ind[1], grid[1])+1
            new_ind[2] = mod(ind[2], grid[2])+1
            new_ind[3] = mod(ind[3], grid[3])+1

            for na = 1:tbc.tb.nwan
                for nb = 1:tbc.tb.nwan
                    hr_g[na,nb,new_ind[1], new_ind[2], new_ind[3], :] += g[na + (nb-1) * tbc.tb.nwan + tbc.tb.nwan^2 * (c-1) , :]
                    

                    #                hr[na,nb,new_ind[1], new_ind[2], new_ind[3]] += ham_orig[na + (nb-1) * tbc.tb.nwan + tbc.tb.nwan^2 * (c-1) ]

                    sr_g[na,nb,new_ind[1], new_ind[2], new_ind[3], :] += g[size_ret + na + (nb-1) * tbc.tb.nwan + tbc.tb.nwan^2 * (c-1),: ]
                    
                end
            end
        end

        begin 
            hk_g = fft(hr_g, [3,4,5])
            sk_g = fft(sr_g, [3,4,5])
            
            VALS0 = zeros(Float64,  prod(grid), tbc.tb.nwan,3*ct.nat+6)

            #        hka = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)
            #        hka = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan, 3*ct.nat+6)

            # hk = fft(hr, [3,4,5])
            #    VALS00 = zeros(Float64, prod(grid), tbc.tb.nwan)        
        end

        #    if false
        #
        #        pVECTS = permutedims(VECTS, [2,3,1])
        #
        #
        #        hk_g2 = zeros(Complex{Float64}, prod(grid), tbc.tb.nwan, tbc.tb.nwan,  3*ct.nat+6)
        #        sk_g2 = zeros(Complex{Float64}, prod(grid), tbc.tb.nwan, tbc.tb.nwan,  3*ct.nat+6)
        #    
        #
        #        for gc = 1:3*ct.nat+6
        #            c=0
        #            for k1 = 1:grid[1]
        #                for k2 = 1:grid[2]
        #                    for k3 = 1:grid[3]
        #                        c += 1
        #                        hk_g2[c, :,:,gc] = hk_g[:,:,k1,k2,k3,gc]
        #                        sk_g2[c, :,:,gc] = sk_g[:,:,k1,k2,k3,gc]
        #                    end
        #                end
        #            end
        #        end
        #
        #        hka = zeros(Complex{Float64}, prod(grid), tbc.tb.nwan, tbc.tb.nwan,  tbc.tb.nwan,3*ct.nat+6)
        #        for gc = 1:3*ct.nat+6
        #            for a = 1:tbc.tb.nwan
        #                for a2 = 1:tbc.tb.nwan
        #                    for a1 = 1:tbc.tb.nwan
        #                        hka[:,a1,a2,a, gc] = hk_g2[:,a1,a2,gc] - VALS[:,a] .*  sk_g2[:,a1,a2,gc]
        #                    end
        #                end
        #            end
        #        end
        #    end        
        #    hka[:,:] =  hk_g[:,:,k1,k2,k3,gc] - ( VALS[c,a] )  *   sk_g[:,:,k1,k2,k3, gc]

        #    if scf
        #        hk0 = hk0 + h1 .* sk
        #    end

#        println("old vals")
        if false
            pVECTS = permutedims(VECTS, [2,3,1])

            if scf
                @threads for a1 = 1:tbc.tb.nwan
                    for a2 = 1:tbc.tb.nwan
                        hk_g[a2,a1,:,:,:,:] += h1[a2,a1]*sk_g[a2,a1,:,:,:,:]
                    end
                end
            end
            @threads for gc = 1:3*ct.nat+6
                c=0
                for k1 = 1:grid[1]
                    for k2 = 1:grid[2]
                        for k3 = 1:grid[3]
                            c += 1
                            for a =1:tbc.tb.nwan
                                #                        hka[:,:] =  hk_g[:,:,k1,k2,k3,gc] - ( VALS[c,a] )  *   sk_g[:,:,k1,k2,k3, gc]


                                VALS0[c, a,gc] = real(pVECTS[:,a,c]' * ( hk_g[:,:,k1,k2,k3,gc] - ( VALS[c,a] )  *   sk_g[:,:,k1,k2,k3, gc]   )  * pVECTS[:,a,c])                  #+ h1 .*  sk_g[:,:,k1,k2,k3, gc]
                            end
                        end
                    end
                end
            end
        end

#        println("new vals")

        psi_gradH_psi(VALS0, VECTS, hk_g, sk_g, h1, VALS, scf, tbc.tb.nwan, ct.nat, grid, OCCS)

        #        @time psi_gradH_psi(VALS0, VECTS, hk_g, sk_g, h1, VALS, scf, tbc.tb.nwan, ct.nat, grid)        
        
        #    @time for gc = 1:3*ct.nat+6
        #        c=0
        #        for k1 = 1:grid[1]
        #            for k2 = 1:grid[2]
        #                for k3 = 1:grid[3]
        #                    c += 1
        #                    for a =1:tbc.tb.nwan
        #                        hka[:,:] =  hk_g[:,:,k1,k2,k3,gc] - ( VALS[c,a] )  *   sk_g[:,:,k1,k2,k3, gc]
        #                        VALS0[c, a,gc] = real(pVECTS[:,a,c]' * hka[c,:,:,a, gc]  * pVECTS[:,a,c])
        #                    end
        #                end
        #            end 
        #        end
        #    end


        #end
        #    end

        garr =zeros(Float64, 3*ct.nat+6)
        for gc = 1:3*ct.nat+6
            garr[gc] = sum(OCCS .* VALS0[:,:,gc]) / nk * 2.0
        end

        if scf
            garr += g[end,:][:]
        end

    end
        
    x, stress = reshape_vec(garr, ct.nat)
    f_cart = -1.0 * x
    f_cart = f_cart * inv(ct.A)'
    stress = -stress / abs(det(ct.A))

    for i = 1:3
        for j = 1:3
            if abs(stress[i,j]) < 1e-12
                stress[i,j] = 0.0
            end
        end
    end
    for i = ct.nat
        for j = 1:3
            if abs(f_cart[i,j]) < 1e-12
                f_cart[i,j] = 0.0
            end
        end
    end


    
    return energy_tot,  f_cart, stress


end


#############################################################################################################


"""
    function psi_gradH_psi(VALS0, VECTS, hk_g, sk_g, h1, VALS, scf, nwan, nat, grid, OCCS)

Helper function for <psi | grad_Ham | psi >
"""
function psi_gradH_psi(VALS0, VECTS, hk_g, sk_g, h1, VALS, scf, nwan, nat, grid, OCCS)

    pVECTS = permutedims(VECTS, [2,3,1])

    if scf
        @threads for a1 = 1:nwan
            for a2 = 1:nwan
                @inbounds hk_g[a2,a1,:,:,:,:] +=  h1[a2,a1]*(@view sk_g[a2,a1,:,:,:,:])
            end
        end
    end
    occs_max = ones(Int64, prod(grid))*nwan

    @threads for c = 1:prod(grid)
        for a = 1:nwan
            if OCCS[c,a] < 1e-4
                occs_max[c] = max(a-1, 1)
                break
            end
        end
    end
        
    @threads for gc = 1:3*nat+6
        c=0
        for k1 = 1:grid[1]
            for k2 = 1:grid[2]
                for k3 = 1:grid[3]
                    c += 1
                    for a =1:occs_max[c]
                        #                        hka[:,:] =  hk_g[:,:,k1,k2,k3,gc] - ( VALS[c,a] )  *   sk_g[:,:,k1,k2,k3, gc]
                        
                        
                        #VALS0[c, a,gc] = real(pVECTS[:,a,c]' * ( hk_g[:,:,k1,k2,k3,gc] - ( VALS[c,a] )  *   sk_g[:,:,k1,k2,k3, gc]   )  * pVECTS[:,a,c])                  #+ h1 .*  sk_g[:,:,k1,k2,k3, gc]
                        @inbounds VALS0[c, a,gc] = real( (@view pVECTS[:,a,c])' * ( (@view hk_g[:,:,k1,k2,k3,gc]) - ( VALS[c,a] )  *   (@view sk_g[:,:,k1,k2,k3, gc])   )  * (@view pVECTS[:,a,c]))                  #+ h1 .*  sk_g[:,:,k1,k2,k3, gc]
                    end
                end
            end
        end
    end

end


"""
    function psi_gradH_psi(VALS0, VECTS, hk_g, sk_g, h1, VALS, scf, nwan, nat, grid, OCCS)

Helper function for <psi | grad_Ham | psi >

This version isn't used.
"""
function psi_gradH_psi2(VALS0, VECTS, hk_g, sk_g, h1, VALS, scf, nwan, nat, grid, OCCS)

    pVECTS = permutedims(VECTS, [2,3,1])

    if scf
        @threads for a1 = 1:nwan
            for a2 = 1:nwan
                hk_g[a2,a1,:,:,:,:] += h1[a2,a1]*(@view sk_g[a2,a1,:,:,:,:])
            end
        end
    end
    occs_max = ones(Int64, prod(grid))*nwan

    @threads for c = 1:prod(grid)
        for a = 1:nwan
            if OCCS[c,a] < 1e-4
                occs_max[c] = max(a-1, 1)
                break
            end
        end
    end

    
    @threads for gc = 1:3*nat+6
        c=0
        for k1 = 1:grid[1]
            for k2 = 1:grid[2]
                for k3 = 1:grid[3]
                    c += 1
                    for a =1:occs_max[c]
                        #                        hka[:,:] =  hk_g[:,:,k1,k2,k3,gc] - ( VALS[c,a] )  *   sk_g[:,:,k1,k2,k3, gc]
                        
                        
                        #                        VALS0[c, a,gc] = real(pVECTS[:,a,c]' * ( hk_g[:,:,k1,k2,k3,gc] - ( VALS[c,a] )  *   sk_g[:,:,k1,k2,k3, gc]   )  * pVECTS[:,a,c])                  #+ h1 .*  sk_g[:,:,k1,k2,k3, gc]

                        VALS0[c, a,gc] = real( (@view pVECTS[:,a,c])' * ( (@view hk_g[:,:,k1,k2,k3,gc]) - ( VALS[c,a] )  *   (@view sk_g[:,:,k1,k2,k3, gc])   )  * (@view pVECTS[:,a,c]))                  #+ h1 .*  sk_g[:,:,k1,k2,k3, gc]                        
                    end
                end
            end
        end
    end

end



end #end module

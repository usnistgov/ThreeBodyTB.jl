###include("Crystal.jl")
###include("DFToutMod.jl")
#using XMLDict



"""
    module CalcTB

Create TB matrix sets from coefficients, or prepare to fit the coefficients.
"""
module ClassicalFit

using ..Utility:arr2str
using ..Utility:str_w_spaces
using ..Utility:parse_str_ARR_float
using ..Utility:dict2str
using ..Utility:str2tuplesdict

using Random
using LinearAlgebra
using SpecialFunctions
using GZip
using EzXML
using XMLDict
using Printf
using SparseArrays
using LoopVectorization
using ..CrystalMod:distances_etc_3bdy_parallel_LV
using ..CalcTB:laguerre_fast!
using ..CalcTB:calc_frontier_list
using ..CalcTB:calc_frontier


using Suppressor

#using Statistics

using Base.Threads
#using Base.lock

using ..AtomicMod:atom
using ..CrystalMod:crystal
#using ..TB:orbital_index
using ..CrystalMod:orbital_index

using ..CrystalMod:get_dist
using ..CrystalMod:get_dist_only

using ..Utility:cutoff_fn
using ..Utility:cutoff_fn_fast
using ForwardDiff
using ..CrystalMod:makecrys

using ..ThreeBodyTB:scf_energy_force_stress
using ..ThreeBodyTB:scf_energy

using Random


using ..CrystalMod:cutoff2X
using ..CrystalMod:cutoff3bX
using ..CrystalMod:cutoff_onX
using ..CrystalMod:cutoff_length
using ..CrystalMod:cutoff4X

using ..Force_Stress:reshape_vec
using ..Atomdata:get_cutoff

using ..DFToutMod:dftout


using ..Classical:make_coefs_cl
using ..Classical:calc_energy_cl
using ..Classical:energy_force_stress_cl


using ..Classical:n_2body_cl
using ..Classical:n_2body_em
using ..Classical:n_3body_cl_same
using ..Classical:n_3body_cl_pair
using ..Classical:n_3body_cl_diff

function kfold(L, n)

    r = randperm(L)
    l = L ÷ n
    IND_TEST = []
    IND_TRAIN = []
    for nn = 1:n
        ind_test = r[(1:l) .+ (nn-1)*l ]
        ind_train = r[setdiff(1:L, (1:l) .+ (nn-1)*l)]
        push!(IND_TEST, ind_test)
        push!(IND_TRAIN, ind_train)
    end

    return IND_TRAIN, IND_TEST
    

end

function atom_perm_to_dist_perm(perm)

    d12 = sort([perm[1], perm[2]])
    d13 = sort([perm[1], perm[3]])
    d23 = sort([perm[2], perm[3]])

    normal_order = [[1,2], [1,3], [2,3]]

    new_order = [findfirst(x->x == normal_order[1], [d12,d13,d23]),
                 findfirst(x->x == normal_order[2], [d12,d13,d23]),
                 findfirst(x->x == normal_order[3], [d12,d13,d23])]

    

    return new_order
    
    
end


function do_fit_cl_RECURSIVE(DFT_start::Array{Any,1}; GEN_FN=missing, procs=procs, thresh=thresh, use_threebody=true, energy_weight=1.0, use_energy=true, use_force=true, use_stress = true , database_start=missing, lambda = -1.0, use_fourbody=false, use_em = true, subtract_scf=false, niters=20 )

    vtot = missing
    rtot = missing
    DFT = deepcopy(DFT_start)
    
    CRYS_tot = []

    database =Dict()
    passtest = false
    for iter = 1:niters

        println("ITER $iter -----------------------------")
        database, vtot, rtot = do_fit_cl(DFT, Vtot_start = vtot, Rtot_start=rtot, use_threebody=use_threebody, energy_weight=energy_weight, use_energy=use_energy, use_force=use_force, use_stress = use_stress , database_start=missing, lambda = lambda, use_fourbody=use_fourbody, use_em = use_em, subtract_scf=subtract_scf, return_mats = true, CRYS_tot = CRYS_tot)

        DFT, NAMES, passtest = GEN_FN(database, do_tb=(subtract_scf), procs=4, ITER=iter )


        println("passtest $passtest")
        if passtest
            println("passtest done")
            break
        end
        
    end

    
    
    return database, passtest
    
end


function do_fit_cl(DFT; Vtot_start = missing, Rtot_start=missing, use_threebody=true, energy_weight=1.0, use_energy=true, use_force=true, use_stress = true , database_start=missing, lambda = -1.0, use_fourbody=false, use_em = true, subtract_scf=false, return_mats = false, CRYS_tot=missing)

    println("DFT version")
    
    if !use_energy && !use_force && !use_stress
        println("everything is false")
        return Dict()
    end
    
    CRYS = crystal[]
    if ismissing(CRYS_tot)
        CRYS_tot = deepcopy(CRYS)
    end

    if use_energy
        ENERGIES = Float64[] #per atom
    else
        ENERGIES=missing
    end
    if use_force
        FORCES = Float64[]
    else
        FORCES = missing
    end
    if use_stress
        STRESSES = Float64[]
    else
        STRESSES = missing
    end

    DFT_energy = Float64[]
    
    for dft in DFT

        energyT = 0.0
        forceT = zeros(dft.crys.nat, 3)
        stressT= zeros(3,3)

        
        #subtract old forces
        if !ismissing(database_start)

            energyT, flag = calc_energy_cl(dft.crys, database=database_start, use_threebody=use_threebody, verbose=false, use_fourbody=use_fourbody, use_em = use_em, turn_off_warn=true)
            if flag == true
                continue
            end
            energyTT, forceTT, stressTT = energy_force_stress_cl(dft.crys, database=database_start, use_threebody=use_threebody, verbose=false, use_fourbody=use_fourbody, use_em = use_em, turn_off_warn=true)
            energyT += energyTT
            forceT += forceTT
            stressT += stressTT
            
        end
            
        
        if subtract_scf
            energyTT, tbc, flag = scf_energy(dft.crys, do_classical=false)
            if flag == false
                println("skip -----------------------")
                println(dft.crys)
                println()
                continue
            end
            if use_stress || use_force
                energyTT, forceTT, stressTT = scf_energy_force_stress(tbc, do_classical=false)   #true, database_classical=database_start)
            else
                forceTT = zeros(dft.crys.nat, 3)
                stressTT = zeros(3,3)
            end
            energyT += energyTT
            forceT += forceTT
            stressT += stressTT
        end

        push!(CRYS, dft.crys)
        push!(CRYS_tot, dft.crys)
        push!(DFT_energy, dft.atomize_energy/dft.crys.nat) 
        if use_energy
            ENERGIES = vcat(ENERGIES, [dft.atomize_energy/dft.crys.nat - energyT/dft.crys.nat])
        end
        if use_force
            FORCES = vcat(FORCES, dft.forces[:] - forceT[:])
        end
        if use_stress
            stress = dft.stress - stressT
            STRESSES = vcat(STRESSES, [stress[1,1], stress[1,2],stress[1,3],stress[2,2],stress[2,3],stress[3,3]])
        end
    end

        
        
    min_en = abs(minimum(DFT_energy))^1.5
    min_w = min_en * 0.4
    weights_train = (abs.(DFT_energy).^1.5 .+ min_w)/(min_en .+ min_w)

#    weights_train = ones(length(CRYS))
    
    #per atom
    ENERGIES = weights_train .* ENERGIES * energy_weight

#    println("ENERGIES ", ENERGIES[1:5])
    
    return do_fit_cl(CRYS, Vtot_start=Vtot_start, Rtot_start=Rtot_start, use_threebody=use_threebody, energy_weight = energy_weight, weights_train = weights_train, ENERGIES=ENERGIES, FORCES=FORCES, STRESSES=STRESSES, database_start=database_start, lambda=lambda, use_fourbody=use_fourbody, use_em = use_em, return_mats=return_mats, CRYS_tot=CRYS_tot)
    
end

function do_fit_cl(CRYS::Array{crystal,1}; Vtot_start = missing, Rtot_start = missing, use_threebody=true,ENERGIES=missing, FORCES=missing, STRESSES=missing, energy_weight = 1.0, database_start=missing, lambda = -1.0, use_fourbody=false, use_em = true, return_mats=false, weights_train = missing, CRYS_tot = missing)

    #ENERGIES PER ATOM. 
    
    println("CRYS version")

    if ismissing(CRYS_tot)
        CRYS_tot = deepcopy(CRYS)
    end
    
    if ismissing(weights_train)
        weights_train = ones(length(CRYS))
    end
    
    if !ismissing(FORCES) || !ismissing(STRESSES)
        get_force=true
    else
        get_force = false
    end

    Ve,Vf,Vs, vars, at_types, ind_set = prepare_fit_cl(CRYS; use_threebody=use_threebody, get_force=get_force, database=database_start, use_fourbody=use_fourbody, use_em=use_em)

    ntot = max(size(Ve)[2], size(Vf)[2], size(Vs)[2])
    if ismissing(Vtot_start)
        Vtot_start = zeros(0,ntot)
    end
    if ismissing(Rtot_start)
        Rtot_start = zeros(0)
    end
    
    if !ismissing(ENERGIES) && !ismissing(FORCES) && !ismissing(STRESSES)
        println("Use energy, force and stress; energy weight $energy_weight")
        Vtot = vcat(Vtot_start, energy_weight*Ve.*weights_train,Vf,Vs)
        Rtot = vcat(Rtot_start, ENERGIES, FORCES, STRESSES)
    elseif !ismissing(ENERGIES) && !ismissing(FORCES)
        println("Use energy and force; energy weight $energy_weight")
        Vtot = vcat(Vtot_start,energy_weight*Ve.*weights_train,Vf)
        Rtot = vcat(Rtot_start,ENERGIES, FORCES)
    elseif !ismissing(FORCES)
        println("Use forces only")
        Vtot = vcat(Vtot_start, Vf)
        Rtot = vcat(Rtot_start, FORCES)
    elseif !ismissing(ENERGIES)
        println("Use energies only")
        Vtot = vcat(Vtot_start, energy_weight*Ve.*weights_train)
        Rtot = vcat(Rtot_start, ENERGIES)
    else
        println("something wrong do_fit_cl ", [ismissing(ENERGIES), ismissing(FORCES), ismissing(STRESSES)])
        return Dict()
    end

    if lambda > 1e-10
        LAM = I(ntot) * lambda
        Vtot_lam = vcat(Vtot, LAM)
        Rtot_lam = vcat(Rtot, zeros(ntot))
        
    elseif lambda < -1e-10

        L = size(Vtot)[1]
        if L < 5
            lambda = 1e-5
            println("WARNING, too few items to auto determine lambda")
        else
            l = L ÷ 5

            #            println("L $L l $l")
            IND_TRAIN, IND_TEST = kfold(L, 5)

            
            LAM = [2.5e-8, 1e-7,2.5e-7,1e-6, 2.5e-6, 1e-5, 2.5e-5, 1e-4, 2.5e-4, 1e-3,2.5e-3, 1e-2, 2.5e-2]
            ERR = zeros(length(LAM))
            
            for n = 1:5  #very simple cross-validation determination of lambda

                ind_test = IND_TEST[n]
                ind_train = IND_TRAIN[n]
                #                r = randperm(L)
#                ind_test = r[1:l]
#                ind_train = r[l+1:end]
                for (cind, lambda) = enumerate(LAM)
                    V = vcat(Vtot[ind_train,:], I(ntot) * lambda)
                    R = vcat(Rtot[ind_train], zeros(ntot))
                    x = V \ R
                    ERR[cind] += sum( (Vtot[ind_test,:]*x - Rtot[ind_test,:]).^2)
#                    println("n $n $cind $cind $lambda ", sum( (Vtot[ind_test,:]*x - Rtot[ind_test,:]).^2))
                end
            end
            for (l,e) in zip(LAM, ERR)
                println("lam   $l  err   $e")
            end
            
#            println("LAM $LAM")
#            println("ERR $ERR")
            
            lambda = LAM[argmin(ERR)] * 1.1
            println("auto lambda $lambda")
        end
        
        LAM = I(ntot) * lambda
        Vtot_lam = vcat(Vtot, LAM)
        Rtot_lam = vcat(Rtot, zeros(ntot))
        
    else
        Vtot_lam = deepcopy(Vtot)
        Rtot_lam = deepcopy(Rtot)

        
    end


    println("size Vtot_lam ", size(Vtot_lam))
    
    x = Vtot_lam \ Rtot_lam

    Vx = Vtot_lam*x

    for (v,r) in zip(Vx,Rtot_lam)
        println("v $v   r $r                $(v-r)")
    end

    
    if ismissing(database_start)
        database = Dict()
    else
        database = deepcopy(database_start)
    end

    frontier = calc_frontier_list(CRYS_tot)
    
    println("put stuff in database")
    for ind in keys(ind_set)
        println("ind $ind")
        if length(ind) == 2
            t1,t2 = ind
#            println("t1 t2 $t1 $t2 ind ", ind_set[ind])
#            println("t1 t2 $t1 $t2 x   ", x[ind_set[ind]])
            if (t1,t2) in keys(frontier)
                c = make_coefs_cl([t1, t2], 2, datH = x[ind_set[ind]], dist_frontier = frontier, version = 1, min_dist = frontier[(t1,t2)])
            else
                c = make_coefs_cl([t1, t2], 2, datH = x[ind_set[ind]], version = 1)
            end                
            database[ind] = c
        elseif length(ind) == 3 && ind[3] == :em
            t1,t2,em_var = ind

            if (t1,t2) in keys(frontier)
                c = make_coefs_cl([t1, t2], 2, datH = x[ind_set[ind]], dist_frontier = frontier, version = 1, min_dist = frontier[(t1,t2)],  em=true)
            end
            database[ind] = c
            
            #database[ind] = make_coefs_cl([t1,t2], x[ind_set[ind]]
            

        elseif length(ind) == 3

            t1,t2,t3 = ind
            if (t1,t2,t3) in keys(frontier)
                c = make_coefs_cl([t1, t2, t3], 3, datH = x[ind_set[ind]], dist_frontier = frontier, version = 1)
            else
                c = make_coefs_cl([t1, t2, t3], 3, datH = x[ind_set[ind]], version = 1)
            end
            database[ind] = c

        elseif length(ind) == 4
            t1,t2,t3,t4 = ind
            if (t1,t2,t3,t4) in keys(frontier)
                c = make_coefs_cl([t1, t2, t3,t4], 4, datH = x[ind_set[ind]], dist_frontier = frontier, version = 1)
            else
                c = make_coefs_cl([t1, t2, t3,t4], 4, datH = x[ind_set[ind]], version = 1)
            end
            database[ind] = c
            
        else
            println("something wrong do_fit_cl $ind")            
        end
            
    end

#    if use_em
#        database[:em] = x[1:5]
#    end

    if return_mats == false
        return database
    else
        return database, Vtot, Rtot
    end
end    

function efs(crys, dat_vars, at_types, ind_set,vars_list, use_threebody, use_fourbody, use_em, DIST)
    var_type = eltype(dat_vars)
#    println("EFS ----------------------------------------------------------------------------------------")
    #energy = calc_energy_cl(crys, dat_vars=dat_vars, at_types=at_types, ind_set=ind_set,vars_list= vars_list, use_threebody=use_threebody)
    energy, force, stress = energy_force_stress_cl(crys, dat_vars=dat_vars, at_types=at_types, ind_set=ind_set,vars_list= vars_list, use_threebody=use_threebody, var_type=var_type, verbose=false, use_fourbody=use_fourbody, use_em = use_em, DIST=DIST)
    #return energy 
    v = vcat(force[:], [stress[1,1], stress[1,2],stress[1,3],stress[2,2],stress[2,3],stress[3,3]])
    #    return v
    #return energy 
end


function prepare_fit_cl(CRYS; use_threebody=true, get_force=true, database=missing, use_fourbody=false, use_em = true)
    println("prepare_fit_cl $use_threebody $use_em ")
    types_dict = Dict()
    types_dict_reverse = Dict()
    types_counter = 0


    vars = Set()
    at_types_set = Set()
    for crys in CRYS
        st = Set(crys.stypes)
        for s1 in st
            push!(at_types_set, s1)
            for s2 in st
                if !ismissing(database)
                    if ! ( (s1,s2) in keys(database))
                        push!(vars,[Set([s1,s2]),2])
                    end 
                else
                    push!(vars,[Set([s1,s2]),2])
                end
                if use_threebody
                    for s3 in st
                        if !ismissing(database)
                            if ! ( (s1,s2,s3) in keys(database))
                                 push!(vars,[Set([s1,s2,s3]), 3])
                            end
                        else
                            push!(vars,[Set([s1,s2,s3]), 3])
                        end
                        if use_fourbody
                            for s4 in st
                                if !ismissing(database)
                                    if ! ( (s1,s2,s3,s4) in keys(database))
                                        push!(vars,[Set([s1,s2,s3,s4]), 4])
                                    end
                                else
                                    push!(vars,[Set([s1,s2,s3,s4]), 4])
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    println("vars ", vars)
    
    vars_list = collect(vars)
    at_types = collect(at_types_set)
    ntot = 0
#    if use_em
#        ntot = 5
#    end
    println("ntot start $ntot")
    at_ind = []

    ind_set = Dict()

    pattern = [ [1,1,1], #these are distances 12, 13, 23

                [2,1,1],
                [1,2,1],
                [1,1,2],

                [3,1,1],
                [1,3,1],
                [1,1,3],

                [2,2,2],

                [1,2,2],
                [2,1,2],
                [2,2,1], 

                [0,1,1], 
                [1,0,1], 
                [1,1,0], 

                [0,1,2], 
                [0,2,1], 
                [2,0,1], 
                [1,0,2], 
                [1,2,0], 
                [2,1,0],

                [0,2,2],
                [2,0,2],
                [2,2,0] ]


    
    for v in vars_list
        println("v $v")
        println()
        if v[2] == 2
            if length(v[1]) == 1
                t1 = collect(v[1])[1]
                t2 = t1
            elseif length(v[1]) == 2
                t1,t2 = collect(v[1])[1:2]
            else
                println("something wrong v length(v[1]) ", length(v[1]) )
            end
            ind_set[(t1,t2)] = collect(1:n_2body_cl) .+ ntot
            ind_set[(t2,t1)] = collect(1:n_2body_cl) .+ ntot
            ntot += n_2body_cl
            if use_em
                ind_set[(t1,t2,:em)] = collect(1:n_2body_em) .+ ntot
                ntot += n_2body_em
                if t1 != t2
                    ind_set[(t2,t1,:em)] = collect(1:n_2body_em) .+ ntot
                    ntot += n_2body_em
                end
            end
                
        elseif v[2] == 3
            println("three ", v[1])            
            if length(v[1]) == 1
                println("elemental")
                t1 = collect(v[1])[1]
                #ind_set[(t1,t1,t1)] = ntot .+ [1,2,2,2,3,3,3,4,5,5,5,6,6,6,7,7,7,7,7,7, 8, 8, 8, 9,9,9,9,9,9]
                ind_set[(t1,t1,t1)] = ntot .+ [1,2,2,2,3,3,3,4,5,5,5,6,6,6,7,7,7,7,7,7, 8, 8, 8]
                #ind_set[(t1,t1,t1)] = ntot .+ [1,2,2,2,3,4,4,4]
                ntot += n_3body_cl_same
                
            elseif length(v[1]) == 2
                t1,t2 = collect(v[1])[1:2]

                #1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20
                ind_set[(t1,t2,t2)  ] = ntot .+ [1, 2, 2, 3, 4, 4, 5, 6, 7, 7, 8, 9, 9,10,13,12,12,13,11,11,14,14,15] #3 special
                ind_set[(t2,t1,t2)  ] = ntot .+ [1, 2, 3, 2, 4, 5, 4, 6, 7, 8, 7, 9,10, 9,12,13,11,11,13,12,14,15,14] #2 special
                ind_set[(t2,t2,t1)  ] = ntot .+ [1, 3, 2, 2, 5, 4, 4, 6, 8, 7, 7,10, 9, 9,11,11,13,12,12,13,15,14,14] #1 special

                ind_set[(t2,t1,t1)  ] = ntot .+ n_3body_cl_pair .+ [1, 2, 2, 3, 4, 4, 5, 6, 7, 7, 8, 9, 9,10,11,12,12,11,13,13,14,14,15       ]  #3 special
                ind_set[(t1,t2,t1)  ] = ntot .+ n_3body_cl_pair .+ [1, 2, 3, 2, 4, 5, 4, 6, 7, 8, 7, 9,10, 9,12,13,11,11,13,12,14,15,14       ] #2 special
                ind_set[(t1,t1,t2)  ] = ntot .+ n_3body_cl_pair .+ [1, 3, 2, 2, 5, 4, 4, 6, 8, 7, 7,10, 9, 9,11,11,13,12,12,13,15,14,14       ] #1 special
                
                            
                ntot += n_3body_cl_pair * 2 #two possible pairs
                            
                
                
            elseif length(v[1]) == 3
                ntot += n_3body_cl_diff
                
                t1,t2,t3 = collect(v[1])[1:3]
                t123 = [t1,t2,t3]
                ind_set[(t1,t2,t3)] = ntot .+ collect(1:n_3body_cl_diff)
                
                PERM = [ [1,2,3],[2,1,3],[3,1,2],[3,2,1],[1,3,2],[2,3,1]] #atom permutations
                for perm in PERM
                    #
                    dist_perm = atom_perm_to_dist_perm(perm)
                    
                    ind_set[ (t123[perm[1]], t123[perm[2]],t123[perm[3]])  ] = Int64[]
                    for p in pattern
                        pnew = p[dist_perm]
                        ind = findfirst(x->x==pnew,pattern)
                        push!(ind_set[ (t123[perm[1]], t123[perm[2]],t123[perm[3]])  ], ind)
                    end
                end
                    
                ntot += n_3body_cl_diff
                    
                    #                for p in pattern
                    #                    
                    #                ind_set[(t2,t1,t3)] = 
                    
            else
                println("something wrong get_fit_cl v[2] == 3  $v")
            end                
        elseif v[2] == 4
            println("four  ", v[1])            
            t1 = collect(v[1])[1]
            if length(v[1]) == 1
                t1 = collect(v[1])[1]
                println("elemental $t1")
                ind_set[(t1,t1,t1,t1)] = ntot .+ [1,2]
                ntot += 2
            end            
            
        else
                println("something wrong get_fit_cl $v")
        end
        println("ntot $ntot ")
        
    end

    println("at_types $at_types")
    println("ind_set")
    println(ind_set)
    println()
#    return ntot, ind_set, vars_list
    
        #    for (ind, crys) in enumerate(CRYS)
        #        en = calc_energy_cl(crys, dat_vars=ones(ntot), at_types=at_types, vars_list= vars_list, use_threebody=use_threebody)
        #        println("ind $ind en $en")
        #    end
        
    V = zeros(length(CRYS), ntot)

    crysX = CRYS[1]
    DIST = []
    function go(x)
        en, _ = calc_energy_cl(crysX, verbose=false, dat_vars=x, at_types=at_types, ind_set=ind_set,vars_list= vars_list, use_threebody=use_threebody, use_fourbody=use_fourbody, use_em = use_em, DIST=DIST, check_frontier=false)
        return en
    end

    chunksize=min(30, ntot)
    cfg = ForwardDiff.GradientConfig(go, ones(ntot), ForwardDiff.Chunk{chunksize}())

    
    for (ind, crys) in enumerate(CRYS)
#        println("ind crys")
        #        println(crys)
        println("indE $ind")
        crysX = crys
#        println(crys)
        begin
            #println("dist")
            if (use_fourbody)
                R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3,array_ind4 = distances_etc_3bdy_parallel_LV(crysX,cutoff2X, cutoff3bX,var_type=Float64, return_floats=false, cutoff4 = cutoff4X)
                DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3, array_ind4
            else
                
                R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel_LV(crysX,cutoff2X, cutoff3bX,var_type=Float64, return_floats=false)
                DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3
            end
            @time V[ind, :] = ForwardDiff.gradient(go, ones(ntot), cfg ) / crys.nat
            #V[ind, :] = ForwardDiff.gradient(x->efs(crys, x, at_types, ind_set,vars_list, use_threebody), ones(ntot))
        end
    end


    if get_force 
        NAT = 0
        for c in CRYS
            NAT += c.nat
        end
        Vf = zeros(NAT*3, ntot)
        Vs = zeros(length(CRYS) * 6, ntot)
        
        counter_f = 0
        counter_s = 0

        function go2(x)
            return efs(crysX, x, at_types, ind_set,vars_list, use_threebody, use_fourbody, use_em, DIST)
        end
        cfg2 = ForwardDiff.JacobianConfig(go2, ones(ntot), ForwardDiff.Chunk{chunksize}())
        
        for (ind, crys) in enumerate(CRYS)
#            println("ind crys")
            #            println(crys)
            println("indEFS  $ind")
            crysX = crys
            begin

                #println("dist")
                if (use_fourbody)
                    R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3,array_ind4 = distances_etc_3bdy_parallel_LV(crys,cutoff2X, cutoff3bX,var_type=Float64, return_floats=false, cutoff4 = cutoff4X)
                    DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3, array_ind4
                else
                    
                    R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel_LV(crys,cutoff2X, cutoff3bX,var_type=Float64, return_floats=false)
                    DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3
                end

                @time ret = ForwardDiff.jacobian(go2, ones(ntot), cfg2 )
                Vf[counter_f .+ (1:crys.nat*3),:] = ret[1:crys.nat*3,:]
                Vs[counter_s .+ (1:6),:] = ret[(crys.nat*3+1):(crys.nat*3+6),:]
                counter_f += 3*crys.nat
                counter_s += 6
            end
        end
    else
        Vf = zeros(0, ntot)
        Vs = zeros(0, ntot)
    end        
        
    
    return V, Vf, Vs, vars, at_types, ind_set
    
end


#function ham(x :: Vector, ct, database, donecheck, DIST, FloatX, use_threebody, dat_vars, at_types, vars_list, ind_set, turn_off_warn, verbose, use_fourbody, use_em)
#    T=eltype(x)
#
#    x_r, x_r_strain = reshape_vec(x, ct.nat, strain_mode=true)
#    A = FloatX.(ct.A) * (I(3) + x_r_strain)
#    crys_dual = makecrys( A , ct.coords + x_r, ct.types, units="Bohr")

 #   energy = calc_energy_cl(crys_dual;  database=database, DIST=DIST, verbose=verbose, use_threebody=use_threebody, var_type=T, dat_vars=dat_vars, at_types=at_types, vars_list=vars_list, ind_set=ind_set, turn_off_warn=turn_off_warn, use_fourbody=use_fourbody , use_em=use_em, check_frontier=false)
##    

#end



end #end module




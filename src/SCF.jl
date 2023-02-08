##include("Crystal.jl")
###include("DFToutMod.jl")
#using XMLDict



####################### Wannier90 specific 
"""
    module SCF

Module for self-consistent field calculations for TB objects.
"""
module SCF
"""
self-consistent field
"""

#include("Atomdata.jl")
using ..Atomdata:atoms
using Printf
using LinearAlgebra
using ..CrystalMod:crystal
using ..TB:make_tb_crys
using ..TB:make_tb
using ..TB:tb
using ..TB:tb_crys
using ..TB:tb_crys_kspace
using ..TB:types_energy

using ..TB:myfft
using ..TB:calc_energy_charge_fft_band
using ..TB:calc_energy_charge_fft_band2
using ..TB:calc_energy_charge_fft
using ..TB:myfft_R_to_K
using ..TB:ewald_energy
using ..TB:magnetic_energy

using ..CrystalMod:orbital_index
using ..TB:get_neutral_eden

using ..TB:get_h1
using ..TB:get_h1_dq
using ..TB:get_spin_h1
using ..TB:get_dq
using ..TB:get_energy_electron_density_kspace
using ..TB:smearing_energy

#using ..CalcTB:calc_tb_lowmem2
using ..CalcTB:calc_tb_LV
#using ..CalcTB:calc_tb_lowmem
using ..TB:get_magmom
using ..CrystalMod:get_grid

using SpecialFunctions
using ..ThreeBodyTB:global_energy_units
using ..ThreeBodyTB:eV

export scf_energy

"""
    function scf_energy(c::crystal, database::Dict; smearing=0.01, grid = missing, conv_thr = 1e-5, iters = 75, mix = -1.0, mixing_mode=:pulay, verbose=true)

Run scf calculation of `c::crystal`, using `database` of `coefs`. The main user version is `scf_energy` in ThreeBodyTB, which calls this one.

- `smearing` is smearing energy in Ryd.
- `grid` is k-point grid (gamma centered MP), will use default.
- `conv_thr` convergence threshold in Ryd.
- `iters` maximum iterations for first attempt
- `mix=-1.0` default is choose mixing for you. Otherwise, set between `0.0` and `1.0`
- `mixing_mode=:pulay` default using Pulay mixing (DIIS). Any other input uses simple mixing.
- `verbose=true` verbosity level.

"""
function scf_energy(c::crystal, database::Dict; smearing=0.01, grid = missing, conv_thr = 1e-5, iters = 100, mix = -1.0, mixing_mode=:pulay, nspin=1, e_den0=missing, verbose=false, repel=true, tot_charge=0.0)

    #println("calc tb")
    tbc = calc_tb_LV(c, database, verbose=verbose, repel=repel, tot_charge=tot_charge);
    #println("asdf ", tbc.eden, ", tc ", tot_charge)
    #println("lowmem")
    #@time tbc = calc_tb_lowmem(c, database, verbose=verbose, repel=repel);
    t = scf_energy(tbc, smearing = smearing, grid=grid, conv_thr = conv_thr, iters=iters, mix=mix,mixing_mode=mixing_mode, nspin=nspin, e_den0=e_den0, verbose=verbose)
    return t
    
end

"""
    function scf_energy(tbc::tb_crys; smearing=0.01, grid = missing, e_den0 = missing, conv_thr = 1e-5, iters = 100, mix = -1.0, mixing_mode=:pulay, verbose=true)
"""
function scf_energy(tbc::tb_crys; smearing=0.01, grid = missing, e_den0 = missing, conv_thr = 0.5e-4, iters = 200, mix = -1.0, mixing_mode=:simple, verbose=true, nspin=1, tot_charge=missing)
"""
Solve for scf energy, also stores the updated electron density and h1 inside the tbc object.
"""

    #println("before before ", tbc.eden)
    
    if !ismissing(tot_charge)
        ind2orb, orb2ind, etotal, nval = orbital_index(tbc.crys)
        tbc.nelec = nval - tot_charge
    end
    
    #println("e_den0 $e_den0")
    if tbc.nspin == 2
        nspin = 2
    end

    if nspin == 2
        magnetic = true
    else
        magnetic = false
    end

#    if magnetic
#        mixing_mode=:simple
#    end
    
    if !ismissing(e_den0)
        if size(e_den0)[1] == 2 && nspin == 1
            println("switch to nspin == 2 because of initial charge density")
            nspin = 2
            ud = sum(e_den0, dims=2)
            if abs(ud[1] - ud[2]) <1e-2
                magnetic = false
            else
                magnetic= true
            end
        end
    end

    a1 = sqrt(sum(tbc.crys.A[1,:].^2))
    a2 = sqrt(sum(tbc.crys.A[2,:].^2))
    a3 = sqrt(sum(tbc.crys.A[3,:].^2))

    extend = false
    if a1 / a2 > 3.0 || a1 / a3 > 3.0 || a2 / a3 > 3.0 || a2 / a1 > 3.0 || a3 / a1 > 3.0 || a3 / a2 > 3.0 #likely a surface, 
        extend=true
#        mixing_mode = :simple
    end


    if mixing_mode != :pulay && mixing_mode != :kfg
        mixing_mode = :simple
        if mix < 0

            if extend
                mix = 0.05
            else
                if tbc.crys.nat <= 10 
                    mix = 0.5
                else
                    mix = 0.3
                end
            end
        end
    else
        if mix < 0

            if extend
                mix = 0.05
            else
                if tbc.crys.nat <= 10 
                    mix= 0.5
                else
                    mix= 0.3
                end
            end
                
        end
    end


    if global_energy_units == "eV"  #for display only, keep units Rydberg internally
        energy_units = eV
    else
        energy_units = 1.0
    end

    
    error_flag = false
    if tbc.scf == false 
        println("doesn't require scf")

        energy, efermi, eden, VECTS, VALS, error_flag =  calc_energy_charge_fft(tbc; grid=grid, smearing=smearing)
#        energy, efermi, eden, VECTS, VALS =  calc_energy_charge_fft(tbc; grid=grid, smearing=smearing)
#        error_flag = false
        dq = get_dq(tbc.crys, eden)

        return energy, efermi, eden, dq, VECTS, VALS, error_flag, tbc
    end

    if verbose
        println("Do SCF")
    end

    tot_charge = 0.0
    for (i, t) in enumerate(tbc.crys.types)
        at = atoms[t]
        z_ion = at.nval
        tot_charge += z_ion / 2.0
    end
#    println("tot_charge calc ", tot_charge, " tbc.nelec ", tbc.nelec)
    
    

    if ismissing(e_den0)
        e_den0 = deepcopy(tbc.eden)
    end
    if nspin == 2 && size(e_den0,1) == 1
        e_den0 = [e_den0;e_den0]
        e_den0[1,:] = e_den0[1,:] + e_den0[1,:]*0.1
        e_den0[2,:] = e_den0[2,:] - e_den0[2,:]*0.1
    end    

    dq = get_dq(tbc.crys, e_den0)

#    println("before ", e_den0)
    if abs(sum(e_den0) - nspin* tbc.nelec/2.0 ) > 1e-5
        if verbose println("bad guess, instead get atoms initial guess") end
        e_den0 = get_neutral_eden(tbc, nspin=nspin, magnetic=magnetic)
        bv = e_den0 .> 1e-5
        e_den0[bv] = e_den0[bv] .- (sum(e_den0) - nspin* tbc.nelec/2.0)/tbc.crys.nat
        
        #        e_den0 = e_den0 .- (sum(e_den0) - nspin* tbc.nelec/2.0)
    end

#    println("e_den0 ", e_den0)
#    println("sum ", sum(e_den0))
    
    if verbose
        if magnetic 
            println("Initial ΔQ: ", (round.(dq; digits=3)))
            println("Initial μB: ", round.(get_magmom(tbc.crys, e_den0), digits=3))
            
        else
            println("Initial ΔQ: ", (round.(dq; digits=3)))
        end
    end

#else
#        if verbose println("Get initial charge density from input") end
#        if abs(sum(e_den0) - nspin*tot_charge) > 1e-5
#            if verbose println("bad guess, instead get neutral atoms initial guess") end
#            e_den0 = get_neutral_eden(tbc, nspin=nspin, magnetic=magnetic)
#        end
#    end
    
    if ismissing(grid)
        grid = get_grid(tbc.crys)
    end

    nk = prod(grid)
    
    if verbose
        println()
        println("Parameters:")
        println("smearing = $smearing conv_thr = $conv_thr, iters = $iters, mix = $mix $mixing_mode , grid = $grid, nspin=$nspin")
        println("nk $nk: $grid")        
        println()
    end

#    if !ismissing(grid)
#        println("grid = $grid")
#    end


    println()
    println("START SCF ----------------")
    e_den = deepcopy(e_den0)

    etypes = types_energy(tbc.crys)
#    println("fft time")
    hk3, sk3 = myfft_R_to_K(tbc, grid)   #we only have to do this once

    thetype=typeof(real(sk3[1,1,1,1,1]))

    

    
    VECTS = zeros(Complex{Float64}, nspin, nk, tbc.tb.nwan, tbc.tb.nwan)
    VALS = zeros(Float64, nspin, nk, tbc.tb.nwan)
                  
    energy_old = 1e12
    delta_energy_old = 1e14

    conv = false
    
    energy_tot = 0.0
    dq = zeros(tbc.crys.nat)
    dq_old = zeros(tbc.crys.nat)
    dq_old2 = zeros(tbc.crys.nat)

    efermi = 0.0


#    println("dq start", round.(dq; digits=2))

    use_kfg = false

    h1, dq = get_h1(tbc, e_den)

    Qpropose = deepcopy(dq)

#=
    nwan = tbc.tb.nwan
    VECTS_w = zeros(Complex{thetype}, nwan, nwan, nk, nspin)
    SK_w = zeros(Complex{thetype}, nwan, nwan, nk)
    DEN_w = zeros(Complex{thetype}, nwan, nwan, nk)
    function test(e_denA)

        h1, dq = get_h1(tbc, e_denA)
        if magnetic 
            h1spin = get_spin_h1(tbc, e_denA)
        else
            h1spin = missing
        end
        energy_band , efermi, e_den_NEW, VECTS, VALS, error_flag = calc_energy_charge_fft_band2(hk3, sk3, tbc.nelec, smearing=smearing, h1=h1, h1spin = h1spin, DEN=DEN_w, VECTS=VECTS_w, SK = SK_w)

        return e_den_NEW
    end

    e_den1 = test(e_den)
    err0 = sum(abs.(e_den1 - e_den))
    println("err 0 $err0")
    
    e_den1_1 = test(e_den1)
    err1 = sum(abs.(e_den1_1 - e_den1))
    println("err 1 $err1")

    e_den1_2 = test(e_den1_1)
    err1_1 = sum(abs.(e_den1_2 - e_den1_1))
    println("err 11 $err1_1")
    
    
    for x = 0.1:0.1:1.5
        e_denN = test((1-x)*e_den +e_den1*x )
        errN = sum(abs.(e_denN -  ((1-x)*e_den +e_den1*x)))
        println("err $x $errN")
    end

    return
=#
    
    function innnerloop(mixA, smearingA, e_denA, conv_thrA, ITERS)
#        println("innnerloop $mixA $ITERS")
        #main SCF loop
        convA = false
        
        if mixing_mode == :pulay
            n1in = copy(e_denA[:])
            n2in = copy(e_denA[:])
            n3in = copy(e_denA[:])
            n1out = copy(e_denA[:])
            n2out = copy(e_denA[:])
            n3out = copy(e_denA[:])

            n1 = copy(e_denA[:])
            n2 = copy(e_denA[:])
            n3 = copy(e_denA[:])

            #B = zeros(3,3)
            #B[1,3] = 1
            #B[2,3] = 1
            #B[3,3] = 0.0
            #B[3,1] = 1
            #B[3,2] = 1

            #B = zeros(4,4)
            B = zeros(7,4)
            B[1,4] = 1
            B[2,4] = 1
            B[3,4] = 1
            B[4,4] = 0.0
            B[4,1] = 1
            B[4,2] = 1            
            B[4,3] = 1
            B[5,1] = 0.00001
            B[6,2] = 0.00001
            B[7,3] = 0.00001

#            R = zeros(3,1)
#            R[3,1] = 1.0

            R = zeros(7,1)
            R[4,1] = 1.0

            R1 = zeros(size(e_denA[:]))
            R2  = zeros(size(e_denA[:]))
            R3  = zeros(size(e_denA[:]))
            n_pulay = zeros(size(e_denA[:]))
        end

        if mixing_mode == :kfg
            Qin = zeros(3,tbc.crys.nat)
            Qout = zeros(3,tbc.crys.nat)
        end
        
        delta_eden = 0.0

#        println("S ", e_denA)
        nreduce = 0
        nkeep = 0
        e_den_NEW = deepcopy(e_denA)

        energy_band = 0.0
        efermi = 0.0

        error_flag = false

        magmom = 0.0
        magmom_old = 0.0
        
        nwan = tbc.tb.nwan
        

        VECTS_w = zeros(Complex{thetype}, nwan, nwan, nk, nspin)
        SK_w = zeros(Complex{thetype}, nwan, nwan, nk)
        DEN_w = zeros(Complex{thetype}, nwan, nwan, nk)

        delta_dq = ones(tbc.crys.nat)*1000.0
        delta_dq2 = ones(tbc.crys.nat)*1000.0
        
        for iter = 1:ITERS

            if false
                println("ΔQ: ", (round.(dq; digits=3)))
                if magnetic 
                    println("μB: ", round.(get_magmom(tbc.crys, e_denA), digits=3))
                end
            end
            
            #            println("e_denA ", e_denA)
            
            delta_dq2[:] = delta_dq[:]
            
            
            dq_old2[:] = dq_old[:]
            dq_old[:] = dq[:]

#            println("e_denA", e_denA)

            if mixing_mode == :kfg && use_kfg
                h1 = get_h1_dq(tbc, Qpropose)
                dq = Qpropose
            else
                h1, dq = get_h1(tbc, e_denA)
            end


#            println("h1 ")
#            println(h1)
            
                
            if mixing_mode == :kfg
                for i = 1:(size(Qin,1)-1)
                    Qin[i,:] = Qin[i+1,:]
                end
#                println("size dq ", size(dq))
#                println("size Qin ", size(Qin))
                Qin[end,:] = dq
            end
            
            if magnetic 
#                h1up, h1dn = get_spin_h1(tbc, e_denA)
#                h1spin = [h1up, h1dn]
                h1spin = get_spin_h1(tbc, e_denA)
            else
                h1spin = missing
            end
#            println("O ", e_denA)
                
#            if iter == 1
#                println("DQ ", round.(dq, digits=2))
#            end

#            try
            #            println("calc_energy_charge_fft_band")

#            println("en1")
            #energy_band , efermi, e_den_NEW, VECTS, VALS, error_flag = calc_energy_charge_fft_band(hk3, sk3, tbc.nelec, smearing=smearingA, h1=h1, h1spin = h1spin)
#            println("en2")

#            println("mix $mixA dq_in   ", round.(dq; digits=3))
            
            energy_band , efermi, e_den_NEW, VECTS, VALS, error_flag = calc_energy_charge_fft_band2(hk3, sk3, tbc.nelec, smearing=smearingA, h1=h1, h1spin = h1spin, DEN=DEN_w, VECTS=VECTS_w, SK = SK_w)

            println("e_den_NEW ", e_den_NEW)
            
#            println("energy band ", energy_band)
            
#            println("check ", energy_band2 - energy_band , " , " , sum(abs.(e_den_NEW  - e_den_NEW2)))
#            println(" e_den_NEW ", e_den_NEW)
            h1NEW, dqNEW = get_h1(tbc, e_den_NEW)
#            println("mix $mixA dq_out  ", round.(dqNEW; digits=3))

            delta_dq[:] = dqNEW - dq

            if mixing_mode == :kfg
                for i = 1:(size(Qout,1)-1)
                    Qout[i,:] = Qout[i+1,:]
                end
                Qout[end,:] = dqNEW
            end
            

#            println("e_den_NEW")
#            println(round.(e_den_NEW[1,:], digits=3))
#            if magnetic
#                println(round.(e_den_NEW[2,:], digits=3))                
#            end
#            catch BadOverlap
#                println("caught bad overlap")
#                energy_tot = -999.0
#                error_flag=true
#                break
#            end

            if error_flag == true
                println("error_flag $error_flag")
                break
            end

            #if iter == 1
            #    e_denA[:] = 0.95 * e_den_NEW[:] .+ 0.05 * e_denA[:]
            #end
            
#            println("N ", e_den_NEW)

#            println("size e_den_NEW $(size(e_den_NEW)) e_denA $(size(e_denA))")

            delta_eden_old = deepcopy(delta_eden)

#            println("size e_den_NEW $(size(e_den_NEW)) e_denA $(size(e_denA))")
            delta_eden = sum(abs.(e_den_NEW - e_denA))
            

#            println("eden_NEW ", round.(e_den_NEW , digits=4))
#            println("edenA    ", round.(e_denA , digits=4))

            
            energy_charge, pot = ewald_energy(tbc, dq)
            if magnetic
                energy_magnetic = magnetic_energy(tbc, e_denA)
                magmom_old = deepcopy(magmom)
                magmom = sum(get_magmom(tbc, e_denA))
            else
                energy_magnetic = 0.0
            end


#            println("energy_magnetic $energy_magnetic")
#            println("energy_charge $energy_charge")
#            println("etypes $etypes")
#            println("energy_band $energy_band")

            energy_tot = etypes + energy_band + energy_charge + energy_magnetic

#            println("energy $energy_tot types $etypes band $energy_band charge $energy_charge magnetic $energy_magnetic")

            #            if iter > 4 && (delta_eden >= delta_eden_old*0.99999 )  #|| delta_energy_old < abs(energy_old - energy_tot)
            if iter == 2 && maximum(delta_dq) > 1.0
                mixA = min(0.1, mixA * 0.5)
            elseif iter > 3 && sum(abs.(delta_dq)) > sum(abs.(delta_dq2)) 
                mixA = max(mixA * 0.7, 0.0001)
                nreduce += 1
                if nreduce > 5 
                    #if nreduce > 3 && mixing_mode == :pulay
                    if mixing_mode != :simple
                        println("switch to :simple")
                        mixing_mode=:simple
                    end
                    #@break
                    #                    mixA = 0.02
#                    nreduce = -5
                end
                @printf("                               reduce mixing: % 6.4f   newerr:  % 10.8f  olderr: % 10.8f \n" , mixA ,sum(abs.(delta_dq)) , sum(abs.(delta_dq2)) )
#                println("delta_energy_old $delta_energy_old new ",  abs(energy_old - energy_tot), " " ,  delta_energy_old < abs(energy_old - energy_tot))

            else
                nkeep +=1
            end
            
            if iter == 100
                mixA = max(mixA * 0.6, 0.001)
                nreduce += 1
                @printf("                               reduce mixing: % 6.4f   olderr:  % 10.8f  newerr: % 10.8f \n" , mixA ,delta_eden_old, delta_eden)
            end

            if nkeep > 5
                mixA = min(0.9, mixA*1.2)
                #                println("keep $mixA")
                nkeep = 0
            end
            

            if mixing_mode == :pulay

                n1in[:] = n2in[:]
                n1out[:] = n2out[:]

                n2in[:] = n3in[:]
                n2out[:] = n3out[:]

                n3in[:] = e_denA[:]
                n3out[:] = e_den_NEW[:]

#                n2in[:] = e_denA[:]
#                n2out[:] = e_den_NEW[:]

                n1[:] = n2[:]
#                n2[:] = e_denA[:]
                n2[:] = n3[:]
                n3[:] = e_denA[:]

            end
                
            
            if mixing_mode == :simple 
                mixA_temp = mixA
                if iter == 1
                    if extend
                        mixA_temp = 0.1
                    else
                        mixA_temp = 0.5
                    end
                end 
                e_denA = e_denA * (1 - mixA_temp ) + e_den_NEW * (mixA_temp )  
                use_kfg = false
            elseif iter <= 3
                if iter == 1
                    if extend
                        mixA_temp = 0.02
                    else
                        if tbc.crys.nat <= 10
                            mixA_temp = 0.05
                        else
                            mixA_temp = 0.05
                        end
                    end
                else
                    if extend
                        mixA_temp = 0.02
                    else
                        mixA_temp = 0.04
                    end

                end

                e_denA = e_denA * (1 - mixA_temp ) + e_den_NEW * (mixA_temp ) 
                use_kfg = false
 
            elseif mixing_mode == :pulay

#                R1 = n1out - n1in
#                R2 = n2out - n2in

#                R1[:] = n2 - n1
 #               R2[:] = e_den_NEW - n2

#                println("dR1 ", sum(abs.(R1t - R1)) , " dR2 ", sum(abs.(R2t - R2)), "n2-n1out ", sum(abs.(n2 - n1out)) , " n1 n1in ", sum(abs.(n1 - n1in)))

#                B[1,1] = R1' * R1
#                B[2,2] = R2' * R2
#                B[1,2] = R1' * R2
#                B[2,1] = R2' * R1


                R1 = (n1out - n1in)
                R2 = (n2out - n2in)
                R3 = (n3out - n3in)

                B[1,1] = R1' * R1
                B[2,2] = R2' * R2
                B[3,3] = R3' * R3
                B[1,2] = R1' * R2
                B[2,1] = R2' * R1
                B[1,3] = R1' * R3
                B[3,1] = R3' * R1
                B[2,3] = R2' * R3
                B[3,2] = R3' * R2

#                c = zeros(2)
                c = zeros(3)
                try
                    c = B \ R
                catch
                    c = [0.333333333,0.3333333333,0.3333333333]
                    c[3] = 1.0 - sum(c[1:2])
                end

                #n_pulay[:] = B \ R
                
#                n_pulay = n2 * c[1,1] + e_den_NEW * c[2,1]
#                n_pulay = e_denA * c[1,1] + e_den_NEW * c[2,1]

#                n_pulay = n1in * c[1,1] + n2in * c[2,1]
#                n_pulay = n1out * c[1,1] + n2out * c[2,1]

#                n_pulay = n1out * c[1,1] + n2out * c[2,1]
                n_pulay = n1out * c[1,1] + n2out * c[2,1] + n3out * c[3,1]

#                println("c $c")
                
                if nspin == 1
                    e_denA[1,:] = (1 - mixA) * e_denA[1,:] + mixA * n_pulay[:]
                elseif nspin == 2
                    e_denA[:,:]  = (1 - mixA) * e_denA[:,:] + mixA * reshape(n_pulay, 2, tbc.tb.nwan)
#                    e_denA[1,:]  = (1 - mixA) * e_denA[1,:] + mixA * n_pulay[1:tbc.tb.nwan]
#                    e_denA[2,:]  = (1 - mixA) * e_denA[2,:] + mixA * n_pulay[(tbc.tb.nwan+1):end]
                end                    
                use_kfg = false

#               e_denA = (1 - mixA) * 

#               
#                R1[:] = n1out - n1in
#                R2[:] = n2out - n2in

#
#                B[1,1] = R1' * R1
#                B[2,2] = R2' * R2
#                B[1,2] = R1' * R2
#                B[2,1] = R2' * R1
#
#                c = B \ R
                
#                println("R1R2 ", sum(abs.(R1)), "   ", sum(abs.(R2)))

#                n_pulay[:] = e_denA * c[1,1] + e_den_NEW * c[2,1]
#                e_denA = (1 - mixA) * e_denA + mixA * n_pulay 
            elseif mixing_mode == :kfg

                #println("mixing_mode == kfg")
                use_kfg = true
                
                dQ = (Qin - Qout).^2
                Qpropose = zeros(tbc.crys.nat)

#                println("dQ")
#                println(dQ)
                
                for i = 1:tbc.crys.nat
                    m = [ones(size(Qin[:,i])) Qin[:,i] Qin[:,i].^2]
                    try
                        c = m \ dQ[:,i]
                        if c[3] > 0.0
                            Qpropose[i] = -c[2] / (2*c[3])
                        else
                            Qpropose[i] = Qout[end,i]
                        end
                    catch
                        Qpropose[i] = Qout[end,i]
                    end
                end

#                println("Qin")
#                println(Qin)
#                println("Qout")
#                println(Qout)
#                println("Qpropose unnorm")
#                println(Qpropose)
#                println("norm er ", sum(Qpropose)/tbc.crys.nat)
                Qpropose = Qpropose .- sum(Qpropose)/tbc.crys.nat
#                println("Qpropose norm")
#                println(Qpropose)
              
                Qpropose = (1 - mixA) * Qin[end,:] + mixA*Qpropose*0.25 + mixA*Qout[end,:]*0.75
                #println("Qpropose ", Qpropose, " mix ", mixA, " dq ", Qin[end,:] - Qpropose, " sum ", sum(Qpropose))
                
            end
            
            

#            if true || iter <= 3
#            else#
#
#                println("pulay/anderson mixing  $mixA")#
#
#            end
            

#            println("sum before e_denA ", sum(e_denA), " ", sum(e_denA, dims=2))
            e_denA = e_denA .- (sum(e_denA) -  nspin* tbc.nelec / 2.0) / (tbc.tb.nwan)  #tot charge is correct            
#            println("sum after  e_denA ", sum(e_denA), " ", sum(e_denA, dims=2))

#            dcharge = sum(abs.(dq - dq_old ))            

            if iter == 1
#                println("SCF CALC $iter energy   $energy_tot                                  $dq ")
                @printf("SCF CALC %04i energy  % 10.8f    \n", iter, energy_tot*energy_units )
#                println("dq   ", round.(dq; digits=3))
                
            else
#                println("SCF CALC $iter energy   $energy_tot   en_diff ", abs(energy_tot - energy_old), "   dq_diff:   $delta_eden    ")
                #                @printf("SCF CALC %04i energy  % 10.8f  en_diff:   %08E  dq_diff:   %08E \n", iter, energy_tot*energy_units, abs(energy_tot - energy_old)*energy_units, delta_eden )
                @printf("SCF CALC %04i energy  % 10.8f  en_diff:   %08E  dq_diff:   %08E    \n", iter, energy_tot*energy_units, abs(energy_tot - energy_old)*energy_units, sum(abs.( delta_dq )) )

                
                
#                println("mix $mixA dq   ", round.(dq; digits=3))

#                println(round.(e_denA[1,:], digits=3))
#                if magnetic
#                    println(round.(e_denA[2,:], digits=3))
#                end
            end
            
            #println("conv ", abs(energy_old - energy_tot), " " ,  sum(abs.( delta_dq )), " " , abs(magmom-magmom_old))
            if abs(energy_old - energy_tot) < conv_thrA * tbc.crys.nat && iter >= 2

                if sum(abs.( delta_dq )) < 0.02 * tbc.crys.nat && abs(magmom-magmom_old) < 0.02 * tbc.crys.nat

                    #                    println("conv thr $(conv_thrA * tbc.crys.nat)   $(0.02 * tbc.crys.nat)   $(0.02 * tbc.crys.nat)")
                    #                if delta_eden < 0.05 * tbc.crys.nat
                #if sum(abs.( delta_dq )) < conv_thrA * tbc.crys.nat * 5000 && abs(magmom-magmom_old) < 0.05 * tbc.crys.nat 
                    convA = true
                    println()
                    eu = energy_tot*energy_units
                    if abs( energy_units - 1)< 1e-5
                        println("YES convergence in $iter iters, energy $eu Ryd. ")
                    else
                        println("YES convergence in $iter iters, energy $eu eV  ")
                    end
                    
                    println("END SCF ------------------")
                    println()
                    return convA, e_denA
                    break
                end
            end

            delta_energy_old = deepcopy(abs(energy_old - energy_tot))
            energy_old = deepcopy(energy_tot)
        end
        return convA, e_denA
    end

#    println("STEP1")
#    conv, e_den = innnerloop(0.2, 0.1, e_den, 1e-4)
#    println("e_den  $e_den")


    if mixing_mode == :pulay
#        println("time pulay")
        conv, e_den = innnerloop(mix, smearing, e_den, conv_thr, iters)  #first step

        if true && conv == false 
            println("pulay mix no convergence, switch to simple mixing")
            mixing_mode = :simple
            e_den = get_neutral_eden(tbc, nspin=nspin, magnetic=magnetic)
#            e_den = 0.5*(e_den + e_denS)
            mix = 0.005

        end
    end

    if mixing_mode == :kfg
#        println("try kfg")
        conv, e_den = innnerloop(mix, smearing, e_den, conv_thr, iters)  #first step
        
#        if conv == false
#            println("kfg NO  converge")
#        else
#            println("kfg YES converge")
#        end            
        
        if true && conv == false 
            println("kfg mix no convergence, switch to simple mixing")
            mixing_mode = :simple
            e_den = get_neutral_eden(tbc, nspin=nspin, magnetic=magnetic)
            mix = 0.01
        end
    end
    

    
    if mixing_mode == :simple
#        println("eden start")
#        println(e_den[1,:])
#        if magnetic
#            println(e_den[2,:])
#        end
        e_den_OLD = deepcopy(e_den)

        conv, e_den = innnerloop(mix, smearing, e_den, conv_thr, iters*5)

    end

#    conv, e_den = innnerloop(mix, smearing, e_den, conv_thr, 200)

    
    if conv == false
        println("WARNING !!!! NO convergence in $iters iters")
        error_flag = true
    end

    tbc.eden = e_den #save charge density!!!!!  side effect!!!!!
    if use_kfg
        h1= get_h1_dq(tbc, Qpropose)
        dq = Qpropose
    else
        h1, dq = get_h1(tbc, e_den)
    end
    tbc.tb.h1 = h1     #moar side effect
    tbc.tb.scf = true  #just double checking
    
    if magnetic
        h1spin = get_spin_h1(tbc, e_den)
        tbc.tb.h1spin = h1spin
        tbc.tb.scfspin = true
    end
    println("ΔQ = ", round.(dq, digits=2))
    if nspin == 2
        println("μB = ", round.(get_magmom(tbc), digits=2))
    end
    println()

    tbc.efermi=efermi
    tbc.energy=energy_tot

    V = permutedims(VECTS[:,:,:,:], [3,4,1,2])

    
    return energy_tot, efermi, e_den, dq, V, VALS, error_flag, tbc

end

#=
#descreening
function remove_scf_from_tbc(tbc::tb_crys; smearing=0.01, grid = missing, e_den = missing)
   """
        this function removes the scf effects from the tight-binding elements
        if you then do an scf solution to the new tbc, you should get the old tbc back.
    """

    if tbc.scf == true
        println("do not need to adjust, already scf")
        return tbc
    end

    hk3, sk3 = myfft_R_to_K(tbc, grid)
    hk3a = deepcopy(hk3)
    
    if ismissing(e_den)
        energy_band_orig , efermi, e_den = calc_energy_charge_fft_band(hk3, sk3, tbc.nelec, smearing=smearing)
    end

    h1, dq = get_h1(tbc, e_den)
    h1 = (h1 + h1')/2.0
    println("dq $dq")

#    println("h1")
#    println(h1[3,3])
#    println(h1[7,7])
#    println(h1[12,12])
#test
    #    dq = zeros(size(dq))
#    h1 = zeros(size(h1))
    
    energy_charge, pot = ewald_energy(tbc, dq)

    println("energy_charge ", energy_charge)
    
    ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(tbc.crys)


    #we need to shift eigenvalues this much
    new_target_energy = energy_band_orig - energy_charge

    
    grid = size(hk3)[3:5]

    sk = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)
    hk = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)


    for k1 = 1:grid[1]
        for k2 = 1:grid[2]
            for k3 = 1:grid[3]

                hk[:,:] = 0.5*(hk3[:,:,k1,k2,k3] + hk3[:,:,k1,k2,k3]')
                sk[:,:] = 0.5*(sk3[:,:,k1,k2,k3] + sk3[:,:,k1,k2,k3]')#

                hk -= sk .* h1 #removing scf term
                
                hk3[:,:,k1,k2,k3] = hk[:,:]#

            end
        end
    end

    energy_band_new , efermi, e_den_NEW = calc_energy_charge_fft_band(hk3, sk3, tbc.nelec, smearing=smearing, h1=h1)
    
    shift = (energy_band_orig - energy_charge - energy_band_new)/nval
#    shift = ( -1.0*energy_charge )/nval
    println("shift $shift")
    for k1 = 1:grid[1]
        for k2 = 1:grid[2]
            for k3 = 1:grid[3]

                hk[:,:] = 0.5*(hk3[:,:,k1,k2,k3] + hk3[:,:,k1,k2,k3]')
                sk[:,:] = 0.5*(sk3[:,:,k1,k2,k3] + sk3[:,:,k1,k2,k3]')

                vals, vects = eigen(hk, sk)
                hk =  sk*vects*diagm(vals .+ shift)*inv(vects)
                hk = (hk + hk')/2.0

                vals_new, vects_new = eigen(hk, sk)

                if sum(abs.(vals_new - (vals .+ shift)  )) > 1e-5
                    println( "ERROR ")
                    println(vals .+ shift)
                    println(vals_new)
                end

                if k1 == 1 && k2 == 1 && k3 == 1
                    println("vals_old ", vals )
                    println("vals_new ", vals_new .- shift)
#                    println("vects_old ")
#                    println(vects)
#                    println("vects_new ")
#                    println(vects_new)
                    
                end
#                hk3[:,:,k1,k2,k3] = hk[:,:]  - sk .* h1


            end
        end
    end
#    energy_band_origX , efermiX, e_denX = calc_energy_charge_fft_band(hk3, sk3, tbc.nelec, smearing=smearing)
#    h1, dq = get_h1(tbc, e_denX)
#    println("dq new", dq)
#    h1 = (h1 + h1')/2.0
#    for k1 = 1:grid[1]
#        for k2 = 1:grid[2]
#            for k3 = 1:grid[3]
#                hk[:,:] = 0.5*(hk3[:,:,k1,k2,k3] + hk3[:,:,k1,k2,k3]')
#                sk[:,:] = 0.5*(sk3[:,:,k1,k2,k3] + sk3[:,:,k1,k2,k3]')#
#                hk -= sk .* h1 #removing scf term
#                hk3[:,:,k1,k2,k3] = hk[:,:]#
#
#            end
#        end
#    end



    energy_band_newX , efermiX, e_den_NEWX = calc_energy_charge_fft_band(hk3, sk3, tbc.nelec, smearing=smearing, h1=h1)
    t = energy_charge + energy_band_newX
    println("energy newold $t $energy_band_orig")

    println("e_den")
    println(e_den)
    println("e_den_NEWX")
    println(e_den_NEWX)



    
#                vals, vects = eigen(hk, sk)
#                hk =  vects*diagm(vals .+ shift)*vects'
                
                
                
#    println("max err ", maximum(abs.(hk3a - hk3)))

    ham_r, S_r, r_dict, ind_arr = myfft(tbc.crys, true, grid, [],hk3, sk3)

    for i = 1:size(S_r)[3]
        S_r[:,:,i] = S_r[:,:,i]'
        ham_r[:,:,i] = ham_r[:,:,i]'
    end       

        
    tb_new = make_tb(ham_r, ind_arr, r_dict, S_r, h1=h1 )


    #with charge density
    tbc_new = make_tb_crys(tb_new, deepcopy(tbc.crys), tbc.nelec, tbc.dftenergy, scf=true, eden=e_den_NEWX, gamma = tbc.gamma)

    h1NEW, dqNEW = get_h1(tbc_new, e_den_NEWX)

#    println("dq")
#    println(dqNEW)
#    println(dq)

#    println("e_den")
#    println(e_den_NEWX)
#    println(e_den_NEW)
#    println(e_den)
#    println(tbc_new.eden)

    return tbc_new
end
=#

"""
    function remove_scf_from_tbc(tbcK::tb_crys_kspace; smearing=0.01, e_den = missing)

This function takes a `tbc_crys_kspace` object that does not require scf and adjusts it 
so that it does require scf, but gives the same energy and band structure.
"""
function remove_scf_from_tbc(tbcK::tb_crys_kspace; smearing=0.01, e_den = missing)
    
    hk3 = deepcopy(tbcK.tb.Hk)
    sk3 = deepcopy(tbcK.tb.Sk)

    
#    energy_orig, e_den_new, VECTS, VALS, error_flag = get_energy_electron_density_kspace(tbcK, smearing=smearing)

    energy_orig, e_den_new, VECTS, VALS, error_flag = get_energy_electron_density_kspace(tbcK, smearing=smearing)
    println("ooooooooo $energy_orig")


    if ismissing(e_den)
        e_den = e_den_new
    end

    h1, dq = get_h1(tbcK, e_den)
    if tbcK.tb.nspin == 2
        h1spin = get_spin_h1(tbcK, e_den)
    else
        h1spin = zeros(2,tbcK.tb.nwan,tbcK.tb.nwan)
    end
    
    energy_charge, pot = ewald_energy(tbcK, dq)
    println("tbcK.tb.nspin ", tbcK.tb.nspin)
    if tbcK.tb.nspin ==2
        energy_magnetic = magnetic_energy(tbcK, e_den)
    else
        energy_magnetic = 0.0
    end

    println("energy_charge ", energy_charge)
    println("energy_magnetic ", energy_magnetic)


    ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(tbcK.crys)
    
    #we need to shift eigenvalues this much
    new_target_energy = energy_orig - energy_charge - energy_magnetic
    
    Hk_new = zeros(Complex{Float64}, size(tbcK.tb.Hk))

    sk = zeros(Complex{Float64}, tbcK.tb.nwan, tbcK.tb.nwan)
                               hk = zeros(Complex{Float64}, tbcK.tb.nwan, tbcK.tb.nwan)

    for spin = 1:tbcK.nspin
        for k = 1:tbcK.tb.nk
            sk = tbcK.tb.Sk[:,:,k]
            #        println(size(sk), " " , size(h1))
            t = sk .* (h1 + h1spin[spin,:,:])
                       #        println(size(t), " t  " , size(Hk_new))
                       
            Hk_new[:,:,k, spin] = tbcK.tb.Hk[:,:,k,spin] - 0.5 * (t + t')
        end
    end

    tbcK_new = deepcopy(tbcK)
    tbcK_new.tb.Hk = Hk_new

    if true
        tbcK_new.scf = true
        tbcK_new.tb.scf = true
        tbcK_new.eden = e_den
        tbcK_new.tb.h1 = h1
#        println("SCF?" , tbcK_new.scf, " " ,tbcK_new.tb.scf)
    end

    tbcK_new.scf = true
    tbcK_new.tb.scf = true
    if tbcK.nspin == 2
        tbcK_new.tb.scfspin = true
    else
        tbcK_new.tb.scfspin = false
    end

    tbcK_new.tb.h1spin = h1spin
    tbcK_new.eden = e_den
    tbcK_new.tb.Hk = Hk_new
    tbcK_new.tb.h1 = h1
    tbcK_new.tb.h1spin = h1spin

#    energy_new, e_den_new, VECTS_new, VALS_new, efermi_new, error_flag = get_energy_electron_density_kspace(tbcK_new.tb, nval, smearing=smearing)
    energy_new, e_den_new, VECTS_new, VALS_new, error_flag = get_energy_electron_density_kspace(tbcK_new, smearing=smearing)
    println("xxxxxxxxxxxxxxx $energy_new")
#    energy_smear = smearing_energy(VALS_new, tbcK.tb.kweights, efermi_new, smearing)
#    println("energy smear " , energy_smear)

#    shift = (energy_orig - energy_charge - energy_magnetic - energy_new)/nval
    shift = (energy_orig - energy_new) / nval

    println("shift $shift   $energy_orig $energy_charge $energy_magnetic $energy_new $nval")
 #    println("shift $shift")

    for spin = 1:tbcK.nspin
        for k = 1:tbcK.tb.nk
            hk[:,:] = Hk_new[:,:,k, spin]
            sk[:,:] = tbcK_new.tb.Sk[:,:,k]
            
            vals, vects = eigen(0.5*(hk+hk'), 0.5*(sk+sk'))
            hk =  sk*vects*diagm(vals .+ shift)*inv(vects)
            
            Hk_new[:,:,k,spin] = 0.5*(hk[:,:]  + hk[:,:]')
            
        end
    end

    tbcK_new.tb.Hk = Hk_new
    tbcK_new.tb.h1 = h1
    tbcK_new.tb.h1spin = h1spin
    tbcK_new.eden = e_den
    tbcK_new.scf = true
    tbcK_new.tb.scf = true
    if tbcK.nspin == 2
        tbcK_new.tb.scfspin = true
    else
        tbcK_new.tb.scfspin = false
    end

#    energy_new, e_den_new, VECTS_new, VALS_new, efermi_new, error_flag = get_energy_electron_density_kspace(tbcK_new.tb, nval, smearing=smearing)
    energy_new, e_den_new, VECTS_new, VALS_new, error_flag = get_energy_electron_density_kspace(tbcK_new, smearing=smearing)
    println("yyyyyyyy $energy_new")

#    println("tbcK_new")
#    println(tbcK_new)
    return tbcK_new

end

"""
    function remove_scf_from_tbc(tbc::tb_crys; smearing=0.01, grid = missing, e_den = missing)

This function takes a `tbc_crys` object that does not require scf and adjusts it 
so that it does require scf, but gives the same energy and band structure.
"""
function remove_scf_from_tbc(tbc::tb_crys; smearing=0.01, grid = missing, e_den = missing)
    if tbc.scf == true
        println("do not need to adjust, already scf")
        return tbc
    end
 


    hk3, sk3 = myfft_R_to_K(tbc, grid)
    #println("size hk3 ", size(hk3))

    hk3, sk3, e_den_NEW, h1, h1spin, efermi = remove_scf_from_tbc(hk3, sk3, tbc; smearing=0.01, e_den = tbc.eden)
#    println("size hk3 v2 ", size(hk3))

    grid = size(hk3)[4:6]

    ctemp = tbc.crys

    hk3a = permutedims(hk3, [3,1,2,4,5,6])
#    println("size hk3a ", size(hk3a))
#    println("before myfft")
    ham_r, S_r, r_dict, ind_arr = myfft(ctemp, true, grid, [],hk3a, sk3)


#    println("after myff")

    if tbc.nspin == 1
        h1spin=missing
    end
        
#    println("size(ham_r) ", size(ham_r))
    tb_new = make_tb(ham_r, ind_arr, r_dict, S_r, h1=h1, h1spin=h1spin )

    #with charge density
    tbc_new = make_tb_crys(tb_new, deepcopy(tbc.crys), tbc.nelec, tbc.dftenergy, scf=true, eden=tbc.eden, gamma = tbc.gamma, background_charge_correction = tbc.background_charge_correction)
    tbc_new.efermi = efermi
    
#    if false
#        println("test")
#        hk3_A, sk3_A = myfft_R_to_K(tbc_new, grid)#
#
#        println("hk3 diff ", sum(abs.(hk3_A - hk3)))
#        println("sk3 diff ", sum(abs.(sk3_A - sk3)))
#    end

    return tbc_new

end


"""
    function remove_scf_from_tbc(hk3, sk3, tbc; smearing=0.01, e_den = missing)

This function takes a hk3, sk3 set of hamiltonian / overlap that does not require scf and adjusts it 
so that it does require scf, but gives the same energy and band structure.
"""
function remove_scf_from_tbc(hk3, sk3, tbc; smearing=0.01, e_den = missing)
   """
        this function removes the scf effects from the tight-binding elements
        if you then do an scf solution to the new tbc, you should get the old tbc back.
    """

    println("1qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq size hk3 ", size(hk3))

#    hk3a = deepcopy(hk3)

    println("aaa ", hk3[1,1,1,1,1,1])
    
    energy_band_orig , efermi, e_den2, VECTS, error_flag = calc_energy_charge_fft_band(hk3, sk3, tbc.nelec, smearing=smearing)


    println("energy_band_orig $energy_band_orig")
    
    println("e_den in ", e_den)
    println("e_den2   ", e_den2)
    
    if ismissing(e_den)
        e_den = e_den2
    end

    h1, dq = get_h1(tbc, e_den)
    if tbc.nspin == 2
        h1spin = get_spin_h1(tbc, e_den)
    else
        h1spin = zeros(2,tbc.tb.nwan, tbc.tb.nwan)
    end
    println("typeof(h1spin) ", typeof(h1spin))

    println("dq $dq")
#test
#    dq = zeros(size(dq))
#    h1 = zeros(size(h1))
    
    energy_charge, pot = ewald_energy(tbc, dq)
    if tbc.nspin == 2
        energy_magnetic = magnetic_energy(tbc, e_den)
        println("en mag $energy_magnetic")

    else
        energy_magnetic = 0.0
    end
    println("energy_charge ", energy_charge)
    println("energy_magnetic ", energy_magnetic)
    
    ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(tbc.crys)


    #we need to shift eigenvalues this much
    new_target_energy = energy_band_orig - energy_charge - energy_magnetic

    
    grid = size(hk3)[4:6]

    sk = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)
    hk = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)

#    println("aaa ", hk3[1,1,1,1,1,1])
#    println("tbc.nspin ", tbc.nspin)
#    println("2qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq size hk3 ", size(hk3))
    
    for spin = 1:tbc.nspin
        for k1 = 1:grid[1]
            for k2 = 1:grid[2]
                for k3 = 1:grid[3]
#                    if k1 == 1 && k2 == 1 && k3 == 1
#                        println("aaa1 ", hk3[1,1,1,k1,k2,k3], " - " , h1[1,1].*sk3[1,1,k1,k2,k3] , " " , h1[1,1], " " , sk3[1,1,k1,k2,k3])
#                    end


                    
                    hk[:,:] = 0.5*(hk3[:,:,spin, k1,k2,k3] + hk3[:,:,spin, k1,k2,k3]')

#                    if k1 == 5 && k2 == 2 && k3 == 3
#                        println("hk SCF xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#xxxxxxxxxx")
#                        println(hk[1:4, 1:4])
 #                       println(sk[1:4, 1:4])
#                        println(h1[1:4, 1:4])
#
 #                   end

                    sk[:,:] = 0.5*(sk3[:,:,k1,k2,k3] + sk3[:,:,k1,k2,k3]')
                    t=sk .* (h1  + h1spin[spin,:,:])
                    hk -= 0.5*(t + t')      #removing scf term
                    
                    hk3[:,:,spin, k1,k2,k3] = hk[:,:]


                    #if k1 == 1 && k2 == 1 && k3 == 1
                    #    println("aaa2 ", hk3[1,1,1,k1,k2,k3])
                    #end
                    
                    
                end
            end
        end
    end

#    println("aaa ", hk3[1,1,1,1,1,1] + h1[1,1] * sk3[1,1,1,1,1], " x ", hk3[1,1,1,1,1,1], " " ,  h1[1,1] * sk3[1,1,1,1,1], " h1 ", h1[1,1], " " , sk3[1,1,1,1,1])
#    println("h1 ", sum(abs.(h1 - h1')))
#    println("sk ", sum(abs.(sk3[:,:,1,1,1] - sk3[:,:,1,1,1]')))
    
    if tbc.nspin == 1
        h1spin = missing
    end
    
    energy_band_new , efermi, e_den_NEW, error_flag = calc_energy_charge_fft_band(hk3, sk3, tbc.nelec, smearing=smearing, h1=h1, h1spin=h1spin)

    println("e_den_NEW ", e_den_NEW)
    println("energy_band_new $energy_band_new")

    
    #    shift = (energy_band_orig - energy_charge - energy_magnetic - energy_band_new)/nval
    sk = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)
    hk = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)

#    shift = (energy_band_orig - energy_charge - energy_magnetic - energy_band_new)/nval
    shift = (energy_band_orig - energy_charge - energy_magnetic - energy_band_new)/nval

    println("shift $shift   $energy_band_orig $energy_charge $energy_band_new $nval")
    for spin = 1:tbc.nspin
        for k1 = 1:grid[1]
            for k2 = 1:grid[2]
                for k3 = 1:grid[3]
                    
                    hk[:,:] = 0.5*(hk3[:,:,spin, k1,k2,k3] + hk3[:,:,spin, k1,k2,k3]')
                    sk[:,:] = 0.5*(sk3[:,:,k1,k2,k3] + sk3[:,:,k1,k2,k3]')

#                    if k1 == 1 && k2 == 1 && k3 == 1
#                        println("aaa3 ", hk3[1,1,1,k1,k2,k3])
#                    end

                    
#                    t=sk .* (h1  + h1spin[spin,:,:])
#                    vals, vects = eigen(hk+t  , sk)

                    
                    vals, vects = eigen(hk  , sk)
                    hk =  sk*vects*diagm(vals .+ shift)*inv(vects)
#                    hk =  sk*vects*diagm(vals .+ shift)*vects'
                    
                    hk3[:,:,spin, k1,k2,k3] = hk[:,:] 
                    
                    
                end
            end
        end
    end

    energy_band_new , efermi, e_den_NEW, error_flag = calc_energy_charge_fft_band(hk3, sk3, tbc.nelec, smearing=smearing, h1=h1, h1spin=h1spin)
    println("energy check ", energy_band_new + energy_charge)

    println("e_den_NEW ", e_den_NEW)
    
#    energy_band_new , efermi, e_den_NEW, error_flag = calc_energy_charge_fft_band2(hk3, sk3, tbc.nelec, smearing=smearing, h1=h1, h1spin=h1spin)
#    println("energy check ", energy_band_new)
    
    return hk3, sk3, e_den_NEW, h1, h1spin, efermi


end

#=
function asdf(tbc::tb_crys; smearing=0.01, grid = missing)
   """
        this function removes the scf effects from the tight-binding elements
        if you then do an scf solution to the new tbc, you should get the old tbc back.
    """


    hk3, sk3 = myfft_R_to_K(tbc, grid)


    grid = size(hk3)[3:5]

    energy_band_orig , efermi, e_den, error_flag = calc_energy_charge_fft_band(hk3, sk3, tbc.nelec, smearing=smearing)

    h1, dq = get_h1(tbc, e_den)

    dq = zeros(size(dq))
    h1 = zeros(size(h1))


    ctemp = tbc.crys
#    ctemp.coords *= -1.0

    ham_r, S_r, r_dict, ind_arr = myfft(ctemp, true, grid, [],hk3, sk3)
    

#
#    for i = 1:size(S_r)[3]
#        S_r[:,:,i] = S_r[:,:,i]'
#        ham_r[:,:,i] = ham_r[:,:,i]'
 #   end       

        
    tb_new = make_tb(ham_r, ind_arr, r_dict, S_r, h1=missing )


    #with charge density
    tbc_new = make_tb_crys(tb_new, deepcopy(tbc.crys), tbc.nelec, tbc.dftenergy, scf=false, gamma = tbc.gamma)


    if true
        println("test")
        hk3_A, sk3_A = myfft_R_to_K(tbc_new, grid)

        println("hk3 diff ", sum(abs.(hk3_A - hk3)))
        println("sk3 diff ", sum(abs.(sk3_A - sk3)))
    end


    return tbc_new
end
=#

end #end module


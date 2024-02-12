
using ..TB:tb_crys_sparse
using ..TB:tb_sparse
using ..TB:myfft_R_to_K_sparse
using ..TB:calc_energy_charge_fft_band2_sym_sparse

using SparseArrays


function scf_energy(tbc::tb_crys_sparse; smearing=0.01, grid = missing, e_den0 = missing, conv_thr = 0.5e-4, iters = 200, mix = -1.0, mixing_mode=:simple, verbose=true, nspin=1, tot_charge=missing, use_sym=true, database_classical=missing, do_classical=true)
"""
Solve for scf energy, also stores the updated electron density and h1 inside the tbc object.
"""
#    println("SCF_ENERGY TOT ", tbc.tot_charge)

    @time begin 
    use_sym = true
    
    if do_classical
        println("DO CLASSICAL $do_classical")
        
        if !ismissing(database_classical)
            energy_classical, _ = calc_energy_cl(tbc.crys, database=database_classical)
            println("energy classical   $energy_classical ")
        else
            energy_classical = 0.0
        end
    else
        energy_classical = 0.0
    end

    
    if mixing_mode == :diis
        mixing_mode = :DIIS
    end


    #println("mixing mode $mixing_mode")
    #println("before before ", tbc.eden)

    if ismissing(tot_charge)
        tot_charge=tbc.tot_charge
    end
    
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

    if mixing_mode == :DIIS && mix < 0.0
        mix = 1.0
    end
    

    if mixing_mode != :pulay  &&  mixing_mode != :DIIS
        mixing_mode = :simple
        if mix < 0

            if extend
                mix = 0.05
            else
                if tbc.crys.nat <= 10 
                    mix = 0.4
                else
                    mix = 0.1
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


    if verbose
        println("Do SCF")
        println()
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
        e_den0[1,:] = e_den0[1,:] + e_den0[1,:]*0.2
        e_den0[2,:] = e_den0[2,:] - e_den0[2,:]*0.2
    end    

    dq = get_dq(tbc.crys, e_den0)
    if tbc.scf == false
        dq .= 0.0
    end
    
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
        println()
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
    
    if use_sym
        nk_red, grid_ind, kpts, kweights = get_kgrid_sym(tbc.crys, grid=grid)
        ind2orb, orb2ind, etotal, nval = orbital_index(tbc.crys)
        sgn, dat, SS, TT, atom_trans = get_symmetry(tbc.crys, verbose=verbose);
    end
                                              

    
    if verbose
        println()
        println("Parameters:")
        println("smearing = $smearing conv_thr = $conv_thr, iters = $iters, mix = $mix $mixing_mode , grid = $grid, nspin=$nspin")
        if use_sym
            println("nk $nk: $grid ; nk_red: $nk_red")
        else
            println("nk $nk: $grid")
        end            
        println()
    end

#    if !ismissing(grid)
#        println("grid = $grid")
#    end


    println()
    println("START SCF ----------------")
    e_den = deepcopy(e_den0)

    etypes = types_energy(tbc.crys)


    thetype=typeof(real(tbc.tb.H[1][1,1]))

    VECTS = zeros(Complex{Float64}, nspin, nk_red, tbc.tb.nwan, tbc.tb.nwan)
    VALS = zeros(Float64, nspin, nk_red, tbc.tb.nwan)
    
                  
    energy_old = 1e12
    delta_energy_old = 1e14

    conv = false
    
    energy_tot = 0.0
    dq = zeros(tbc.crys.nat)
    dq_old = zeros(tbc.crys.nat)
    dq_old2 = zeros(tbc.crys.nat)

    efermi = 0.0


#    println("dq start", round.(dq; digits=2))


    h1, dq = get_h1(tbc, e_den)
    if tbc.scf == false
        dq .= 0.0
    end

    
    Qpropose = deepcopy(dq)

        
    end
    println("done begin")
    #@time hk3, sk3 = myfft_R_to_K(tbc, grid)   #we only have to do this once
    println("fft time")
    @time hk3, sk3 = myfft_R_to_K_sparse(tbc, nk_red, kpts)
    

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

#    DenMat = []
#    HK = []    

    function innnerloop(mixA, smearingA, e_denA, conv_thrA, ITERS)
#        println("innnerloop $mixA $ITERS")
        #main SCF loop
        convA = false

#        println("INNER")
        
        rho_in = []
        rho_out = []        
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
        

        VECTS_w = zeros(Complex{thetype}, nwan, nwan, nk_red, nspin)
        SK_w = zeros(Complex{thetype}, nwan, nwan, nk_red)
        DEN_w = zeros(Complex{thetype}, nwan, nwan, nk_red)
        
        delta_dq = ones(tbc.crys.nat)*1000.0
        delta_dq2 = ones(tbc.crys.nat)*1000.0
        
#        println("error_flag $error_flag")
        for iter = 1:ITERS

#            println("iter $iter $error_flag")
            
            if false
                println("ΔQ: ", (round.(dq; digits=3)), " " , sum(dq))
                if magnetic 
                    println("μB: ", round.(get_magmom(tbc.crys, e_denA), digits=3))
                end
            end
            
            #            println("e_denA ", e_denA)
            
            delta_dq2[:] = delta_dq[:]
            
            
            dq_old2[:] = dq_old[:]
            dq_old[:] = dq[:]

#            println(sum(e_denA) , " e_denA  before ", round.(e_denA, digits=4))

            h1, dq = get_h1(tbc, e_denA)
            if tbc.scf == false
                dq .= 0.0
            end
#            println("new dq ", dq, " --------------------------------------------")


#            println("h1 ")
#            println(h1)
            
                
            
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

#            println("iter $iter add rho_in")
            push!(rho_in, e_denA)
#            println("iter $iter add rho_in $(length(rho_in))")            
            println("calc_energy")


            energy_band , efermi, e_den_NEW, VECTS, VALS, error_flag = calc_energy_charge_fft_band2_sym_sparse(hk3, sk3, tbc.nelec, smearing=smearingA, h1=h1, h1spin = h1spin, DEN=DEN_w, VECTS=VECTS_w, nk_red=nk_red, grid_ind=grid_ind, kweights=kweights)

#            println("e_den_NEW ", e_den_NEW, " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
            #                push!(DenMat, denmat)
#                push!(HK, hk)

#                DIIS(min(iter, 5), HK, SK_w, DenMat, nk_red, grid_ind, kweights)
                

#            println(sum(e_den_NEW), " e_den_NEW  presym ", round.(e_den_NEW, digits=8))
            
            if use_sym
                for spin =1:nspin
                    e_den_NEW[spin,:] = symmetrize_charge_den(tbc.crys, e_den_NEW[spin,:] , SS, atom_trans, orb2ind)
                end
            end

#            println(sum(e_den_NEW), " e_den_NEW  after ", round.(e_den_NEW, digits=8))
            
            push!(rho_out, e_den_NEW)

            h1NEW, dqNEW = get_h1(tbc, e_den_NEW)

            if tbc.scf == false
                dqNEW .= 0.0
            end

            
            delta_dq[:] = dqNEW - dq


            if error_flag == true
                println("error_flag $error_flag")
                break
            end


            delta_eden_old = deepcopy(delta_eden)

            delta_eden = sum(abs.(e_den_NEW - e_denA))
            
#            println("ewald dq $dq")
            energy_charge, pot = ewald_energy(tbc, dq)
#            println("energy_charge, $energy_charge")

            if magnetic
                energy_magnetic = magnetic_energy(tbc, e_denA)
                magmom_old = deepcopy(magmom)
                magmom = sum(get_magmom(tbc, e_denA))
            else
                energy_magnetic = 0.0
            end

#            println("energy_charge $energy_charge energy_band $energy_band etypes $etypes energy_magnetic $energy_magnetic ec $energy_classical")            
            energy_tot = etypes + energy_band + energy_charge + energy_magnetic + energy_classical

            #            if iter > 4 && (delta_eden >= delta_eden_old*0.99999 )  #|| delta_energy_old < abs(energy_old - energy_tot)
            if iter == 2 && maximum(delta_dq) > 1.0
                mixA = min(0.1, mixA * 0.5)
            elseif iter > 3 && sum(abs.(delta_dq)) > sum(abs.(delta_dq2)) 

                if mixing_mode != :DIIS || nreduce > 10
                    mixA = max(mixA * 0.8, 0.0001)
                end

                nreduce += 1
                #println("nreduce $nreduce")
                if nreduce > 20
                    #if nreduce > 3 && mixing_mode == :pulay
                    if mixing_mode != :simple
                        println("switch to :simple")
                        mixing_mode=:simple
                    end
                    #@break
                    #                    mixA = 0.02
#                    nreduce = -5
                end
                #@printf("                               reduce mixing: % 6.4f   newerr:  % 10.8f  olderr: % 10.8f \n" , mixA ,sum(abs.(delta_dq)) , sum(abs.(delta_dq2)) )
#                println("delta_energy_old $delta_energy_old new ",  abs(energy_old - energy_tot), " " ,  delta_energy_old < abs(energy_old - energy_tot))

            else
                nkeep +=1
            end
            
            if iter == 100
                mixA = max(mixA * 0.6, 0.001)
                nreduce += 1
#                @printf("                               reduce mixing: % 6.4f   olderr:  % 10.8f  newerr: % 10.8f \n" , mixA ,delta_eden_old, delta_eden)
            end

            if nkeep > 5
                mixA = min(1.0, mixA*1.2)
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
            elseif iter <= 2
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


            elseif mixing_mode == :DIIS
#                if nreduce > 8
#                    mixA = 0.5
#                else
#                    mixA = 1.0
#                end
                e_denA = DIIS(min(iter, n_diis), nwan, nspin, rho_in, rho_out, mixA)
                
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
#                println("dq   ", round.(dq; digits=5))
                
            else
#                println("SCF CALC $iter energy   $energy_tot   en_diff ", abs(energy_tot - energy_old), "   dq_diff:   $delta_eden    ")
                #                @printf("SCF CALC %04i energy  % 10.8f  en_diff:   %08E  dq_diff:   %08E \n", iter, energy_tot*energy_units, abs(energy_tot - energy_old)*energy_units, delta_eden )
                if nspin == 2
                    @printf("SCF CALC %04i energy  % 10.8f  en_diff:   %08E  dq_diff:   %08E   dmag: %08E   mix: %06E \n", iter, energy_tot*energy_units, abs(energy_tot - energy_old)*energy_units, sum(abs.( delta_dq )),  abs(magmom-magmom_old), mixA )
                else
                    @printf("SCF CALC %04i energy  % 10.8f  en_diff:   %08E  dq_diff:   %08E   mix: %06E \n", iter, energy_tot*energy_units, abs(energy_tot - energy_old)*energy_units, sum(abs.( delta_dq )), mixA )
                end

                
                
#                println("mix $mixA dq   ", round.(dq; digits=3))

#                println( sum(abs.( delta_dq )) ,"  x " , round.(e_denA[1,:], digits=3))
#                if magnetic
#                    println(round.(e_denA[2,:], digits=3))
#                end
            end
            
#            println("conv ", abs(energy_old - energy_tot), " " ,  sum(abs.( delta_dq )), " " , abs(magmom-magmom_old))
#            println([abs(energy_old - energy_tot) < conv_thrA * tbc.crys.nat,  sum(abs.( delta_dq )) < 0.02 * tbc.crys.nat, abs(magmom-magmom_old) < 0.02 * tbc.crys.nat])
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


    println("inner loop")
    @time begin
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

    
    if mixing_mode == :DIIS

        n_diis = 15
        
        conv, e_den = innnerloop(mix, smearing, e_den, conv_thr, iters)  #first step
#        if conv == false
#            println("kfg NO  converge")
#        else
#            println("kfg YES converge")
#        end
        if conv == false
            mixing_mode = :simple
            e_den = get_neutral_eden(tbc, nspin=nspin, magnetic=magnetic)
            mix = 0.05
        end

    end        

    
    if mixing_mode == :simple
#        println("eden start")
#        println(e_den[1,:])
#        if magnetic
#            println(e_den[2,:])
#        end
        e_den_OLD = deepcopy(e_den)

        conv, e_den = innnerloop(mix, smearing, e_den, conv_thr, iters)

    end

#    conv, e_den = innnerloop(mix, smearing, e_den, conv_thr, 200)


    
    
    if conv == false
        println("WARNING !!!! NO convergence in $iters iters")
        error_flag = true
    end
    end
    println("done inn")
    @time begin 
        if size(e_den) == size(tbc.eden)
            tbc.eden[:] = e_den #save charge density!!!!!  side effect!!!!!
        else #i don't think i need this option anymore
            tbc= make_tb_crys(tbc.tb, tbc.crys, tbc.nelec, tbc.dftenergy, scf=true, eden=e_den, gamma=tbc.gamma, background_charge_correction=tbc.background_charge_correction, within_fit=tbc.within_fit, tb_energy=energy_tot, fermi_energy=efermi)
        end
        
        h1, dq = get_h1(tbc, e_den)
        if tbc.scf == false
            dq .= 0.0
        end

        tbc.tb.h1[:,:] = h1     #moar side effect
        tbc.tb.scf = true  #just double checking
        tbc.dq = dq #for convenience
        tbc.tot_charge=-sum(dq) #for convenience
        
        if magnetic
            h1spin = get_spin_h1(tbc, e_den)
            tbc.tb.h1spin[:,:,:] = h1spin
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

        #    println("FINAL error_flag ", error_flag)
    end
    println("done wrap")
        
    return energy_tot, efermi, e_den, dq, V, VALS, error_flag, tbc

end


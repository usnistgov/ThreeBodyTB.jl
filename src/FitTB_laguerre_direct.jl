using Optim
using LineSearches

using ..TB:get_h1_dq

using ..Symmetry:get_symmetry
using ..Symmetry:symmetrize_charge_den

function topstuff_direct(list_of_tbcs, prepare_data; EDEN_input=missing, weights_list=missing, dft_list=missing, kpoints = [0 0 0; 0 0 0.5; 0 0.5 0.5; 0.5 0.5 0.5], starting_database = missing,  update_all = false, fit_threebody=true, fit_threebody_onsite=true, do_plot = false, energy_weight = missing, rs_weight=missing, ks_weight = missing, niters=50, lambda=[0.0, 0.0], leave_one_out=false, RW_PARAM=0.0001, KPOINTS=missing, KWEIGHTS=missing, nk_max = 0)

#    database_linear, ch_lin, cs_lin, X_Hnew_BIG, Y_Hnew_BIG,               X_H,               X_Snew_BIG, Y_H, h_on,              ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3,keepind, keepdata = prepare_data
    
    database_linear, ch_lin, cs_lin, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, Ys_new, cs, ch_refit, SPIN, threebody_inds  = prepare_data

    println("AAAAAAAA ch_lin ", ch_lin)

    if false
    list_of_tbcs_new = []
    for (c, tbck) in enumerate(list_of_tbcs)
        tbck_new = missing
        @suppress begin
            tbck_new = remove_scf_from_tbc(deepcopy(tbck); smearing=0.01, e_den = EDEN_input[c,[1],1:tbck.tb.nwan])
        end
        push!(list_of_tbcs_new, tbck_new)
    end
    else
        list_of_tbcs_new = list_of_tbcs        
    end
    
    ##########    
    
    (ch_keep, keep_inds, toupdate_inds, cs_keep, keep_inds_S, toupdate_inds_S) = keepdata



    list_of_tbcs_new = list_of_tbcs_new[keepind]
    if !ismissing(dft_list)
        dft_list = dft_list[keepind]
    end
    KPOINTS = KPOINTS[keepind]
    KWEIGHTS = KWEIGHTS[keepind]

    if ismissing(dft_list)
        dft_list = []
        for i = 1:length(KPOINTS)
            push!(dft_list, missing)
        end
    end
    
    println("niters $niters")
    
    #SCF MODE
    scf = false
    for m in list_of_tbcs_new
        if !ismissing(m)
            scf = m.scf
            break
        end
    end
    println("scf $scf -----------------------------------------------------------------------------------------------------------------")
#    return
    
#    if list_of_tbcs_new[1].scf == true
#        scf = true
#    end
    println("SCF is $scf")
    
    if ismissing(energy_weight)
        energy_weight = 1.0
    end
    if ismissing(rs_weight)
        rs_weight = 0.0
    end
    if ismissing(ks_weight)
        ks_weight = 0.0
    end


    if ismissing(weights_list)
        weights_list = ones(Float64, length(list_of_tbcs_new))
    else
        weights_list = weights_list[keepind]
    end

    println("update_all $update_all")
    
    println("TOUPDATE_INDS ", length(toupdate_inds))

    #    cs = cs_lin
    #    ch = ch_lin
    
    
    NWAN_MAX = maximum(ind_BIG[:,3])
    SPIN_MAX= maximum(SPIN)
    NAT_MAX = 0
    for tbc in list_of_tbcs_new
        NAT_MAX = max(NAT_MAX, tbc.crys.nat)
    end

#    nk = size(kpoints)[1]
    NCALC = length(list_of_tbcs_new)

    VALS     = zeros(NCALC, nk_max, NWAN_MAX, SPIN_MAX)
    #VECTS     = zeros(Complex{Float64}, NCALC, nk_max, NWAN_MAX, NWAN_MAX, SPIN_MAX)
    VECTS_ref = Dict{Int64, Array{Complex{Float64},4} }()
    S_ref = Dict{Int64, Array{Complex{Float64},3} }()
    
    VALS0     = zeros(NCALC, nk_max, NWAN_MAX, SPIN_MAX)

    E_DEN     = zeros(NCALC, SPIN_MAX, NWAN_MAX)
    H1     = zeros(NCALC, NWAN_MAX, NWAN_MAX)
    H1spin     = zeros(NCALC, 2, NWAN_MAX, NWAN_MAX)
    DQ     = zeros(NCALC, NAT_MAX)

    ENERGY_SMEAR = zeros(NCALC)

    OCCS     = zeros(NCALC, nk_max, NWAN_MAX, SPIN_MAX)
    WEIGHTS     = zeros(NCALC,  nk_max, NWAN_MAX, SPIN_MAX)
    ENERGIES = zeros(NCALC)

    Ys = Ys_new + Xc_Snew_BIG
    X_Snew_BIG = nothing
    Xc_Snew_BIG = nothing

    NCOLS_orig = size(X_Hnew_BIG)[2]
    NCOLS = size(X_Hnew_BIG)[2]

#    println("redo lsq")
#    @time ch = X_H \ Y_H 
 
    ch=deepcopy(ch_lin)

    println("ch start ", ch)
    
    if !ismissing(ch_refit)
        println("using refit")
        for i = 1:length(ch)
            if abs(ch_refit[i]) > 1e-5
                ch[i] = ch_refit[i]
            end
        end
    end



#    if start_small
#        ch = ch / 10.0
#    end
   
    keep_bool = true


    println("NCOLS $NCOLS")
    



    #PREPARE REFERENCE ENERGIES / EIGENVALUES
    println("prepare reference eigs")
    #println([length(list_of_tbcs_new), length(KPOINTS), length(KWEIGHTS), length(dft_list), length(SPIN)])
    c=0
    NVAL = zeros(Float64, length(list_of_tbcs_new))
    NAT = zeros(Int64, length(list_of_tbcs_new))
    @time for (tbc, kpoints, kweights, d, spin ) in zip(list_of_tbcs_new, KPOINTS, KWEIGHTS, dft_list, SPIN)
        c+=1

#        println("c $c")
#        println("kpoints $kpoints")
#        println("kweights ", sum(kweights))
        
        NAT[c] = tbc.crys.nat

        #        nw = ind_BIG[c, 3]
        nw = tbc.tb.nwan
#        println("NW $nw")
        nk = size(kpoints)[1]
        nspin = tbc.nspin
        row1, rowN, nw = ind_BIG[c, 1:3]

        vmat     = zeros(Complex{Float64}, nspin, nk, nw, nw)
        smat     = zeros(Complex{Float64}, nk, nw, nw)


        wan, semicore, nwan, nsemi, wan_atom, atom_wan = tb_indexes(tbc.crys)
        ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(tbc.crys)

        if !ismissing(d)
            nval = d.bandstruct.nelec - nsemi * 2
            println("nval v1 d $nval")
        elseif sum(tbc.eden) > 1e-10
            nval = sum(tbc.eden) / nspin * 2.0
            println("nval v2 tbc $nval")            
        else
            nval = nval
            println("nval v3 default $nval")                        
        end

        println("$c nval $nval tbc")
        println(tbc)
        
        
        NVAL[c] = nval

        if !ismissing(d)
            band_en = band_energy(d.bandstruct.eigs[:,nsemi+1:end, :], d.bandstruct.kweights, nval)
            etypes = types_energy(d.crys)
            etot_dft = d.energy
            e_smear = d.energy_smear
            atomization_energy = etot_dft - etotal_atoms - etypes  - e_smear
            band_en = band_en 
            shift = (atomization_energy - band_en  )/nval
        end
        #println("c atomization $atomization_energy $etot_dft $etotal_atoms $etypes $e_smear $fit_to_dft_eigs")
#        println("nk $nk")
#        println(kpoints)

        println("xx")

        vects_arr = zeros(Complex{Float64}, nk, nw, nw, spin)
        s_arr = zeros(Complex{Float64}, nk, nw, nw)
        for spin in 1:tbc.nspin
            for k in 1:nk
                println("ref k $k spin $spin c $c k $k")
                #                println("nw $nw")
                if !ismissing(tbc) 
                    vects, vals, hk, sk, vals0 = Hk(tbc, kpoints[k,:], spin=spin)  #reference

#                    println("size(vals) ",size(vals))
#                    println("size VALS ", size(VALS))
#                    println("$c $k $nw $spin")

                    VALS[c,k,1:nw, spin] = vals                           #reference
                    VALS0[c,k,1:nw,spin] = vals0                          #reference
                    #VECTS[c,k,1:nw,1:nw,spin] = vects
                    vects_arr[k,:,:,spin] = vects
                    s_arr[k,:,:] = sk
                    if c == 1 && k == 1
                        println("ref vals")
                        println(vals)
                        println("ref vals0")
                        println(vals0)
                    end
                    vmat[spin, k, :,:] = vects
                    smat[k, :,:] = sk
                    #if fit_to_dft_eigs
                    #    n = min(d.bandstruct.nbnd - nsemi, nw)    
                    #    nval_int = Int64(round(nval))
                    #    #                    println("mix fit_to_dft_eigs $fit_to_dft_eigs old $c $k ", VALS[c,k,1:n])
                    #    for k2 in 1:d.bandstruct.nks
                    #        if sum(abs.(d.bandstruct.kpts[k2,:] - kpoints[k,:]) ) < 1e-5
                    #            n = min(d.bandstruct.nbnd - nsemi, nw)    
                    #            VALS[c, k,1+nval_int:n, spin] = 0.5 * VALS[c, k,1+nval_int:n,spin] + 0.5*(d.bandstruct.eigs[k2,nsemi+1+nval_int:nsemi+n, spin] .+ shift)
                    #            #                            println("change")
                    #            break
                    #        end
                    #    end
                    #    #                    println("mix fit_to_dft_eigs $fit_to_dft_eigs new $c $k ", VALS[c,k,1:n])
                    #end
                else
                    
                    #               println("missing tbc $c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!!!!!!x" )
                    for spin = 1:d.nspin
                        for k2 in 1:d.bandstruct.nks
                            if sum(abs.(d.bandstruct.kpts[k2,:] - kpoints[k,:]) ) < 1e-5
                                VALS[c, k,:, spin] .= 100.0
                                VALS0[c, k,:,spin] .= 100.0
                                
                                n = min(d.bandstruct.nbnd - nsemi, nw)
                                VALS[c, k,1:n,spin] = d.bandstruct.eigs[k2,nsemi+1:nsemi+n, spin] .+ shift
                                VALS0[c, k,1:n,spin] = d.bandstruct.eigs[k2,nsemi+1:nsemi+n, spin] .+ shift
                                
                                
                                #                        println("$c DFT_val ",  d.bandstruct.eigs[k2,nsemi+1:end] .+ shift)
                            end
                        end
                    end
                end
                #            VECTS_MIX[c,k, 1:nw, 1:nw] = vects

                
            end

        end
        VECTS_ref[c] = deepcopy(vects_arr)
        S_ref[c] = deepcopy(s_arr)
#        println("VALS ", VALS[c, 1, :,:])
#        println("size ", size(VALS[c, 1:nk,1:nw, 1:d.nspin]), " " , size(kweights), " d.nspin ", d.nspin, " tbc.nspin ", tbc.nspin, " nw $nw nk $nk sum kweights ", sum(kweights))
#        println("kweights[1:6] of ", size(kweights), "  " , kweights[1:6])

        energy_tmp,  efermi = band_energy(VALS[c, 1:nk,1:nw,1:tbc.nspin], kweights, nval, 0.01, returnef=true) 
        println("nval for band energy $nval efermi $efermi")
#        println("energy_tmp $energy_tmp $efermi $efermi nval $nval")
        
        
        occs = gaussian.(VALS[c,1:nk,1:nw,1:tbc.nspin].-efermi, 0.01)
        
#        println("sum occs ", sum(sum(occs[:,:,:], dims=[2,3]).* kweights))

#        println("occs early $efermi $nval ",occs)
#        println("VALS ", VALS)
        
        energy_smear = smearing_energy(VALS[c, 1:nk,1:nw,1:tbc.nspin], kweights, efermi, 0.01)
        

        if !ismissing(tbc) && typeof(tbc) <: tb_crys
            println("tbc.eden ", tbc.eden)
            eden, h1, dq, h1spin = get_electron_density(tbc, kpoints, kweights, vmat, occs, smat)        
            E_DEN[c,1:tbc.nspin, 1:nw] = eden
            println("E_DEN new ", eden)
            H1[c,1:nw, 1:nw] = tbc.tb.h1
            DQ[c,1:tbc.crys.nat] = get_dq(tbc)
            H1spin[c,:,1:nw, 1:nw] = tbc.tb.h1spin

        end
        if !ismissing(tbc) && typeof(tbc) == tb_crys_kspace{Float64}
            eden, h1, dq, h1spin = get_electron_density(tbc, kpoints, kweights, vmat, occs, smat)        
            E_DEN[c,1:tbc.nspin, 1:nw] = eden
#            println("x ", tbc.tb.h1)
            
            H1[c,1:nw, 1:nw] = tbc.tb.h1
            DQ[c,1:tbc.crys.nat] = get_dq(tbc)
#            println("SIZE tbc.tb.h1spin ", size(tbc.tb.h1spin))
            H1spin[c,:,1:nw, 1:nw] = tbc.tb.h1spin
        end

#        println("dq $dq")
        
        if !ismissing(tbc)
#            println("tbc")
#            println(tbc)
            s1 = sum(occs .* VALS0[c,1:nk,1:nw,1:tbc.nspin], dims=[2,3])[:]
            energy_band = sum(s1 .* kweights) #/ tbc.nspin
#            println("ENERGY_BAND ", energy_band, " " , tbc.nspin)
#            println("before ", typeof(tbc), " " , typeof(dq))
            if scf
                energy_charge, pot = ewald_energy(tbc, dq)
            else
                energy_charge = 0.0
            end

#            println("energy magnetic ", tbc.tb.scfspin)
            if tbc.tb.scfspin
                energy_magnetic = magnetic_energy(tbc, eden)
#                println("eden ", eden)
            else
                energy_magnetic = 0.0
            end
#            println("energy_magnetic $energy_magnetic")
            etypes = types_energy(tbc)
            println(min(Int64(ceil(nval/2.0-1e-5)) + 1, nw)," $c CALC ENERGIES $etypes $energy_charge $energy_band $energy_smear $energy_magnetic = ", etypes + energy_charge + energy_band + energy_smear + energy_magnetic)
#            println("VALS ", VALS[c,1:nk,1:nw,1:tbc.nspin])
            
            ENERGIES[c] = etypes + energy_charge + energy_band + energy_smear + energy_magnetic

            s1 = sum(occs .* VALS[c,1:nk,1:nw,1:tbc.nspin], dims=[2,3])[:]
            energy_band_vals = sum(s1 .* kweights) #/ tbc.nspin
            
            shift_value = (energy_charge + energy_band - energy_band_vals)/nval
            println("shift_value c $c = $shift_value")
            
        else
            ENERGIES[c] = etot_dft - etotal_atoms
        end


        #ENERGIES[c] = tbc.dftenergy

        OCCS[c, 1:nk,1:nw, 1:nspin] = occs

        if true

            nweight = min(Int64(ceil(nval/2.0-1e-5)) + 1, nw)
            occs2 = zeros(size(occs))
            occs2[ :,1:nweight,1:nspin   ] .= 0.5
            occs2[ :,nweight  ,1:nspin   ] .= 0.5
            occs2[ :,1        ,1:nspin   ] .= 0.5
#            occs2[1,1] = 5.0

            if !ismissing(tbc)
                etmp, occs3 = band_energy(VALS[c,1:nk,1:nw,1:nspin], kweights, nval+0.25, 0.1, returnocc=true)
            else
                etmp, occs3 = band_energy(VALS[c,1:nk,1:nw,1:nspin], kweights, nval+0.25, 0.1, returnocc=true)
            end
           
            WEIGHTS[c,1:nk,1:nw,1:nspin] = (occs + occs2*0.5 + occs3*0.5)/(1 + 0.5 + 0.5)


#        else
#
#            nweight = min(Int64(ceil(nval/2.0-1e-5)) + 1, nw)
#            occs2 = zeros(size(occs))
#            occs2[:,1:nweight] .= 0.3
#            occs2[:,nweight] .= 0.1
#            occs2[:,1] .= 0.5
#            if !ismissing(tbc)
#                etmp, occs3 = band_energy(VALS[c,1:nk,1:nw], kweights, nval+0.95, 0.1, returnocc=true)
#            else
#                etmp, occs3 = band_energy(VALS[c,1:nk,1:nw], kweights, nval+0.25, 0.1, returnocc=true)
#            end
#           
#            WEIGHTS[c,1:nk,1:nw] = (occs + occs2 + occs3)/3.0

        end


        if !ismissing(tbc)
            WEIGHTS[c,:,:,:] .+= RW_PARAM
        end

#double check
        for s = 1:nspin
            for k = 1:nk
                for i = 1:nw
                    if VALS[c,k,i,s] > 99.0 
                        WEIGHTS[c,k,i,s] = 0.0
                    end
                end
            end
        end

        if ismissing(tbc)
            tbc_fake = calc_tb_fast(d.crys, repel=false)
            tbc = tbc_fake
        end


#        WEIGHTS[c,1,1:nw] = WEIGHTS[c,1,1:nw] * 100

    end


    println("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
    println("REFERENCE ENERGIES ")

    for c = 1:NCALC
        println(ENERGIES[c], "  $c   ", DQ[c,:]  ,  "    EDEN ", E_DEN[c,1, :][:], " sum ", sum(E_DEN[c,1, :][:]), " NVAL ", NVAL[c] )
    end
    println("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE")

#    return 

    return list_of_tbcs_new,KPOINTS, KWEIGHTS, dft_list, scf, energy_weight, rs_weight, ks_weight, weights_list, NWAN_MAX, NCALC, VALS, VALS0, E_DEN, H1, H1spin, DQ, ENERGY_SMEAR, OCCS, WEIGHTS, ENERGIES, NCOLS_orig, NCOLS, ch, NVAL, NAT , SPIN_MAX, Ys, keep_bool, keep_inds, toupdate_inds, ch_keep, keep_inds_S, toupdate_inds_S, cs_keep, VECTS_ref, S_ref

end


function do_fitting_direct(list_of_tbcs_nonscf ; weights_list = missing, dft_list=missing, kpoints = missing, starting_database = missing,  update_all = false, fit_threebody=true, fit_threebody_onsite=true, do_plot = false, energy_weight = missing, rs_weight=missing,ks_weight=missing, niters=50, lambda=[0.0,0.0], leave_one_out=false, prepare_data = missing, RW_PARAM=0.0, NLIM = 100, refit_database = missing, start_small = false, fit_to_dft_eigs=false, fit_eam=false, ch_startX = missing, energy_diff_calc = false, gen_add_ham=false, fitting_version = fitting_version_default, opt_S = true, conjgrad=false, cs_startX = missing, use_sym=true)

    println("do_fitting_direct   niters $niters update_all $update_all fit_threebody $fit_threebody fit_threebody_onsite $fit_threebody_onsite  energy_weight $energy_weight  rs_weight $rs_weight ks_weight $ks_weight lambda $lambda RW_PARAM $RW_PARAM NLIM $NLIM fit_eam $fit_eam energy_diff_calc $energy_diff_calc opt_S $opt_S ")

    if true
    list_of_tbcs = []
    for tbck in list_of_tbcs_nonscf
        tbck_new = missing
        @suppress begin
            tbck_new = remove_scf_from_tbc(tbck; smearing=0.01)
        end
        push!(list_of_tbcs, tbck_new)
    end
    else
list_of_tbcs  = deepcopy(list_of_tbcs_nonscf)        
    end
    
    #list_of_tbcs  = deepcopy(list_of_tbcs_nonscf)
    
    #if ismissing(kpoints)
    #    kpoints, kweights = make_kgrid([2,2,2])
    #end
#    if !ismissing(dft_list)
#        println("top")
    KPOINTS, KWEIGHTS, nk_max = get_k(dft_list, list_of_tbcs, NLIM=NLIM)
#    else
#        println("bot")
#        KPOINTS, KWEIGHTS, nk_max = get_k_simple(kpoints, list_of_tbcs)
##        println("KPOINTS ", KPOINTS)
#    end



    
#    println("KWEIGHTS 3 ", size(KWEIGHTS[3]), " " , KWEIGHTS[3][1:6])
    
    if ismissing(prepare_data)
        println("DO LINEAR FITTING")

        if update_all == true
            starting_database_t = missing #keep all
        else
            starting_database_t = starting_database
        end

        pd = do_fitting_linear(list_of_tbcs; kpoints = KPOINTS, mode=:kspace, dft_list = dft_list,  fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = false, starting_database=starting_database_t, return_database=false, NLIM=NLIM, refit_database=refit_database, fit_eam=fit_eam, ch_startX = ch_startX, fitting_version=fitting_version)
    else
        println("SKIP LINEAR MISSING")
        pd = prepare_data
#        database_linear, ch_lin, cs_lin, X_Hnew_BIG, Y_Hnew_BIG, X_H, X_Snew_BIG, Y_H, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3 = prepare_data
    end

    return do_fitting_direct_main(list_of_tbcs_nonscf,list_of_tbcs, pd; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = ks_weight, niters=niters, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam, ch_startX = ch_startX, energy_diff_calc = energy_diff_calc, gen_add_ham=gen_add_ham, fitting_version=fitting_version, opt_S = opt_S, cg = conjgrad, cs_startX = cs_startX, use_sym=use_sym)
end

function do_fitting_direct_main(list_of_tbcs_nonscf, list_of_tbcs, prepare_data; weights_list=missing, dft_list=missing, kpoints = missing, starting_database = missing,  update_all = false, fit_threebody=true, fit_threebody_onsite=true, do_plot = false, energy_weight = missing, rs_weight=missing, ks_weight = missing, niters=50, lambda=[0.0, 0.0], leave_one_out=false, RW_PARAM=0.0001, KPOINTS=missing, KWEIGHTS=missing, nk_max=0, start_small=false, fit_to_dft_eigs=false, fit_eam=false, optimS = false, top_vars = missing, ch_startX = missing, energy_diff_calc = false, gen_add_ham=false, fitting_version=fitting_version_default, opt_S=true, cg = false, cs_startX = missing, use_sym=true)

    leave_out = -1

    solve_scf_mode = false
    
    println("start do_fitting_recursive_direct")
    if typeof(lambda) <: Float64
        lambda = [lambda, lambda]
    end

    println("test scf ", list_of_tbcs_nonscf[1].scf, " ", list_of_tbcs[1].scf)

    database_linear, ch_lin, cs_lin, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, Ys_new, cs, ch_refit, SPIN, threebody_inds  = prepare_data

#    println("ch_lin ", ch_lin)
#    println("cs_lin ", cs_lin)
#    return missing
    
    list_of_tbcs_nonscf,KPOINTS, KWEIGHTS, dft_list, scf, energy_weight, rs_weight, ks_weight, weights_list, NWAN_MAX, NCALC, VALS, VALS0, E_DEN, H1, H1spin, DQ, ENERGY_SMEAR, OCCS, WEIGHTS, ENERGIES, NCOLS_orig, NCOLS, ch, NVAL, NAT, SPIN_MAX, Ys, keep_bool, keep_inds, toupdate_inds, ch_keep, keep_inds_S, toupdate_inds_S, cs_keep, VECTS_ref, S_ref =
        topstuff_direct(list_of_tbcs_nonscf, prepare_data; weights_list=weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight, ks_weight = ks_weight, niters=niters, lambda=lambda,  leave_one_out=false, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max = nk_max)

    SSS = []
    ORB2IND = []
    ATOMTRANS = []
    if use_sym
        for tbc in list_of_tbcs_nonscf
            sgn, dat, SS, TT, atom_trans = get_symmetry(tbc.crys);
            push!(SSS, SS)
            push!(ATOMTRANS, atom_trans)
            ind2orb, orb2ind, etotal, nval = orbital_index(tbc.crys)
            push!(ORB2IND, orb2ind)
            
        end
    end
    

    
    Ys = Ys_new + Xc_Snew_BIG

    NCOLS_S = size(X_Snew_BIG)[2]

#    println("start ncols ", [NCOLS, NCOLS_S]  , " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    
    weights_list = weights_list ./ maximum(weights_list)

    if !ismissing(ch_startX)
        ch[1:length(ch_startX)] = ch_startX
    end
    if !ismissing(cs_startX)
        cs[1:length(cs_startX)] = cs_startX
    end
    #    ch = ones(size(ch))
    
    println("NOW, DO DIRECT FITTING")
    
    VALS_working = zeros(size(VALS))
    ENERGIES_working = zeros(size(ENERGIES))
    OCCS_working = zeros(size(OCCS))

    energy_diff_weight = 1.0
    verbose=0

    ERROR = zeros(Int64, NCALC)

    println("scf $scf")
    
    EDEN_FITTED = zeros(NCALC, SPIN_MAX, NWAN_MAX)

    
    
    function construct_fitted(ch, cs, solve_self_consistently = false)

#        println("construct_fitted ch $ch cs $cs")
        
        Xc = (X_Hnew_BIG * ch) + Xc_Hnew_BIG
        Ys = X_Snew_BIG * cs + Xc_Snew_BIG
        
        #        VECTS_FITTED     = zeros(Complex{Float64}, NCALC, nk_max, NWAN_MAX, NWAN_MAX)
        VECTS_FITTED = Dict{Int64, Array{Complex{Float64},4} }()
        VALS_FITTED     = zeros(NCALC, nk_max, NWAN_MAX, SPIN_MAX)
        VALS0_FITTED     = zeros(NCALC, nk_max, NWAN_MAX, SPIN_MAX)
        OCCS_FITTED     = zeros(NCALC, nk_max, NWAN_MAX, SPIN_MAX)
        ENERGIES_FITTED = zeros(NCALC)
        c=0


        if solve_self_consistently == true
            #            println("SCF SOLUTIONS")
            #niter_scf = 150
            niter_scf = 500
        else
            niter_scf = 1
        end
        
        for (tbc, kpoints, kweights, dft) in zip(list_of_tbcs, KPOINTS, KWEIGHTS, dft_list)
            c+=1
            #            println("construct_fitted c $c")

            nval = NVAL[c]

            etypes = types_energy(tbc.crys)

            nw = ind_BIG[c, 3]
            nk = size(kpoints)[1]

            row1, rowN, nw = ind_BIG[c, 1:3]

            H0 = zeros(Complex{Float64}, nw, nw)
            H = zeros(Complex{Float64}, nw, nw)
            S = zeros(Complex{Float64}, nw, nw)

            Smat = zeros(Complex{Float64}, nk, nw, nw)
            H0mat = zeros(Complex{Float64}, nk, nw, nw)

            N = hermetian_indpt(nw)

            SS = zeros(Complex{Float64}, nw, nw)

            energy_band = 0.0

            #println("get h")
            for k in 1:nk 
                
                #fitted eigenvectors / vals
                for i = 1:nw
                    for j = 1:nw
                        if i <= j
                            ind = hermetian_index(i,j, nw)
                            H0[i,j] = Xc[row1-1 + (k-1)*N + ind] + im*Xc[row1-1 + (k-1)*N + ind + nk * N]  + h_on[c][i, j]
                            S[i,j] = Ys[row1-1 + (k-1)*N + ind] + im*Ys[row1-1 + (k-1)*N + ind + nk * N]
                            
                            H0[j,i] = Xc[row1-1 + (k-1)*N + ind] - im*Xc[row1-1 + (k-1)*N + ind + nk * N]  + h_on[c][j, i]
                            S[j,i] = Ys[row1-1 + (k-1)*N + ind] - im*Ys[row1-1 + (k-1)*N + ind + nk * N]
                            
                            #                            SS[i,j] = Y_Snew_BIG[row1-1 + (k-1)*N + ind] + im*Y_Snew_BIG[row1-1 + (k-1)*N + ind + nk * N]
                            #                            SS[j,i] = Y_Snew_BIG[row1-1 + (k-1)*N + ind] - im*Y_Snew_BIG[row1-1 + (k-1)*N + ind + nk * N]
                            
                            #                            H0[i,j] = Xc[row1-1 + (k-1)*nw*nw + (i-1)*nw + j] + im*Xc[row1-1 + (k-1)*nw*nw + (i-1)*nw + j + nk*nw*nw]  + h_on[c][i, j]
                            #                            S[ i,j] = Ys[row1-1 + (k-1)*nw*nw + (i-1)*nw + j] + im*Ys[row1-1 + (k-1)*nw*nw + (i-1)*nw + j + nk*nw*nw]
                        end
                    end
                    S[ i,i] += 1.0 #onsite
                end
                H0 = (H0+H0')/2.0
                S = (S+S')/2.0
                
                
                Smat[k,:,:] = S
                H0mat[k,:,:] = H0
                
                
            end
            
            if scf
#                h1 = deepcopy(H1[c,1:nw,1:nw])
                dq = deepcopy(DQ[c,1:tbc.crys.nat])
                h1 = get_h1_dq(tbc, dq)

            else
                h1 = zeros(Float64, nw, nw)
                dq = zeros(tbc.crys.nat)
            end
            if tbc.tb.scfspin
                h1spin = deepcopy(H1spin[c,:,1:nw,1:nw])
            else
                h1spin = zeros(Float64,2, nw, nw)
            end

            
            #            println("STARTING DQ")
            #            println(dq)
            #            println("STARTING H1")
            #            println(h1)
            #            println()
            function get_eigen(h1_in, h1spin_in, nspin, iter)
                vals = zeros(Float64, nw)
                vects = zeros(Complex{Float64}, nw, nw)
                VECTS = zeros(Complex{Float64}, nspin, nk, nw, nw)
                for spin = 1:nspin
                    for k in 1:nk 
                        
                        H0 = H0mat[k,:,:]
                        S  = Smat[k,:,:]
                        if scf
                            H = H0 + h1_in .* S
                        else
                            H[:,:] = H0[:,:]
                        end
#                                                if c == 1 && k == 1
#                                                    println("k $k ------------------------------------------------")
#                                                    println("H = ", H)
#                                                    println("H0 = ", H0)
#                                                    println("S = ", S)                            
#                                                end
                        if tbc.tb.scfspin
                            H += h1spin_in[spin,:,:] .* S
                        end
                        if iter == 1
                            try
                                valsS = eigvals( 0.5*(S+S'))
                                if minimum(real.(valsS)) < 0.001
                                    println("WARNING, S has negative eigenvalues $c $k valsS $valsS , set error to true")
                                    ERROR[c] = 1
                                    continue
                                end
                                
                            catch
                                valsS = zeros(1)
                            end
                        end                             
                        try
                            vals, vects = eigen(0.5*(H+H'), 0.5*(S + S'))
#                            if c == 1
#                                println("H ")
#                                println(H)
#                                println("S")
#                                println(S)
#                                println()
#                            end
                            if iter == 1
                                if maximum(abs.(imag.(vals))) > 1e-10
                                    println("WARNING, S has negative eigenvalues $c $k valsS $valsS , set error to true")
                                    ERROR[c] = 1
                                    continue
                                end
                            end 

                            
#                            if k == 1 && c == 5
#                                println("get eigen h1[1,1] = $(h1_in[1,1]) ")
#                                println("get eigen vals[1] = $(vals[1]) ")
#                                println("get eigen H0[1,1] = $(H0[1,1]) ")
#                                println("get eigen S[1,1] = $(S[1,1]) ")
#                                println("get eigen dq = $(dq) ")
#                                println()
#                            end
#                            
                        catch
                            
                            ERROR[c] = 1
                            vals = zeros(Float64, nw)
                            vects = zeros(Complex{Float64}, nw, nw)
                            println("WARNING, S has negative eigenvalues $c $k")
                            
                        end

                        #                    VECTS_FITTED[c,k,1:nw,1:nw] = vects

                        

                        VECTS[spin, k,1:nw,1:nw] = vects
                        VALS_FITTED[c,k,1:nw, spin] = real(vals)
                        #if k == 1 && c == 18
                        #    println("vals diff $(sum(abs.(VALS_FITTED[c,k,1:nw, spin] - VALS[c,k,1:nw, spin]))) 1 5 $(real(vals[1:4])) $(VALS[c,k,1:4, spin])")
                        #end
                        VALS0_FITTED[c,k,1:nw, spin] =  real.(diag(vects'*H0*vects))

                        

                        #                        println("VALS_FITT $spin $k ", real.(vals))
                        #                    VALS_FITTED[c,k,1:nw] = real.(diag(vects'*H*vects))


                        #                        println("VALS0_FITTED c $c spin $spin k $k ", VALS0_FITTED[c,k,1:nw, spin])
#                                            if k == 1 && c == 1
#                                                println("VALS ", VALS_FITTED[c,k,1:nw])
#                                                println("VALS0", VALS0_FITTED[c,k,1:nw])
#                                            end
                    end #kpoints loop
                end

                
                VECTS_FITTED[c] = deepcopy(VECTS)
                
                
            end
            
            
            energy_old = 1000.0
            energy_diff_old = 0.0
            conv = false

            mix = 0.5
            
            for c_scf = 1:niter_scf
                #                println("$c_scf aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")                

                #println("c_scf $c_scf get_eigen")
#                if c_scf == 1 && c == 1
#                    println("sum h1 ", sum(abs.(h1)) , "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
#                end
                get_eigen(h1, h1spin, tbc.tb.nspin, c_scf)

#                if c_scf == 1 && c == 1
#                    println("tbc.crys ")
#                    println(tbc.crys)
#                    println("eigs ", VALS_FITTED[c, 1,1:nw,1:tbc.tb.nspin])
#                end
                
                
                #                if c_scf == 1
                #                    println("$c VALS_FITTED k0 ", VALS_FITTED[c,1,1:nw])
                #                    println("$c VALS0_FITTED k0 ", VALS0_FITTED[c,1,1:nw])
                #                end

                #if solve_self_consistently == true
                #println("other stuff")
                if true
                    
                    #mix = 0.05 * 0.95^c_scf  #start aggressive, reduce mixing slowly if we need most iterations
                    
                    
                    #                    energy_tmp,efermi   = band_energy(VALS_FITTED[c,1:nk,1:nw], kweights, tbc.nelec, 0.01, returnef=true) 

                    #efermi = calc_fermi(VALS_FITTED[c,1:nk,1:nw, 1:tbc.tb.nspin], kweights, nval, 0.01)
                    energy_tmp,  efermi = band_energy(VALS_FITTED[c, 1:nk,1:nw,1:tbc.tb.nspin], kweights, nval, 0.01, returnef=true) 

                    occs = gaussian.(VALS_FITTED[c,1:nk,1:nw, 1:tbc.tb.nspin].-efermi, 0.01)

                    #                    println("occs late $efermi $nval ",occs)
                    #                    println("valsf ", size(VALS_FITTED), " " , size(VALS_FITTED[c,1:nk,1:nw, 1:tbc.tb.nspin]))
                    
                    energy_smear = smearing_energy(VALS_FITTED[c,1:nk,1:nw, 1:tbc.tb.nspin], kweights, efermi, 0.01)

                    #                    eden, h1_new, dq_new = get_electron_density(tbc, kpoints, kweights, VECTS_FITTED[c,:,1:nw,1:nw], occs, Smat)         #updated h1

                    eden, h1_new, dq_new, h1spin_new = get_electron_density(tbc, kpoints, kweights, VECTS_FITTED[c], occs, Smat)         #updated h1

                    if use_sym
                        for spin =1:tbc.tb.nspin
                            eden[spin,:] = symmetrize_charge_den(tbc.crys, eden[spin,:] , SSS[c], ATOMTRANS[c], ORB2IND[c])
                        end
                        h1_new, dq_new = get_h1(tbc, eden)
                        if tbc.tb.nspin > 1
                            h1spin_new = get_spin_h1(tbc, eden)
                        end
                    end
                    
                    #                    println("EDEN DDDDDDDDDD ", eden)

                    #if solve_self_consistently
                    EDEN_FITTED[c, 1:tbc.nspin,1:nw] = eden
                    #end
                    
                    s1 = sum(occs .* VALS0_FITTED[c,1:nk,1:nw, 1:tbc.tb.nspin], dims=[2,3])
                    energy_band = sum(s1 .* kweights) #/ tbc.tb.nspin
                    #                    println("energy_band ", occs[1] , " " ,  VALS0_FITTED[1], " ", kweights[1])
                    
#                    if maximum(abs.(dq - dq_new)) > 0.1
#                        mix = mix * 0.5
#                    end

#                    h1 = h1*(1-mix) + h1_new * mix
#                    dq = dq*(1-mix) + dq_new * mix
#                    h1spin = h1spin*(1-mix) + h1spin_new * mix

                    if solve_self_consistently == true
                    #if false

                        dq = dq*(1-mix) + dq_new * mix
                        
                        #h1 = h1*(1-mix) + h1_new * mix
#                        h1spin = h1spin*(1-mix) + h1spin_new * mix


#                    else
#                        h1 .= 0.0
#                        h1spin .= 0.0
#                        energy_charge = 0.0
#                        energy_magnetic = 0.0
                    end

                    h1 = get_h1_dq(tbc, dq)
                    
                    if scf
                        energy_charge, pot = ewald_energy(tbc, dq)
                        if tbc.tb.scfspin
                            energy_magnetic = magnetic_energy(tbc, eden)
                        else
                            energy_magnetic = 0.0
                        end
                    else
                        energy_charge = 0.0
                        energy_magnetic = 0.0
                    end
                    
                    energy_new = energy_charge + energy_band + energy_smear + energy_magnetic + etypes
                    #                    if c == 34
                    #                        println( " scf $c_scf $c ", energy_new+etypes, "    $dq   $energy_charge $energy_band $energy_smear $energy_magnetic")
                    #                    end
                    
                    #                    println( " scf $c_scf $c ", energy_new+etypes, "    $dq   $energy_charge $energy_band $energy_smear $energy_magnetic")

#                    if c == 5
#                        println("c $c   $c_scf mix $mix energy_new $energy_new energy_old $energy_old diff $(abs(energy_new  - energy_old))  dq $dq")
#                    end
                    if (abs(energy_new  - energy_old) < 1e-4 && c_scf >= 2) ||  (abs(energy_new  - energy_old) < 1e-3 && c_scf >= 20) 
                        #                        println("scf converged  $c ", energy_new+etypes, "    $dq " )
                        conv = true
                        break
                    end

                    if abs(energy_new - energy_old) > energy_diff_old || maximum(abs.(dq - dq_new)) > 0.1
                        mix = mix * 0.85
                    end
                    energy_diff_old = abs(energy_new - energy_old)
                    energy_old = energy_new
                    
                end
                
            end #scf loop

            if scf && conv && solve_self_consistently
                #                h1 = (h1 + H1[c,1:nw,1:nw]) / 2.0   #mixing
                #                dq = (dq + DQ[c,1:tbc.crys.nat]) / 2.0

                DQ[c,1:tbc.crys.nat] = dq
                H1spin[c,:,1:nw, 1:nw] = h1spin
#                println("update dq c $c ", DQ[c,1:tbc.crys.nat])
            #elseif scf && solve_self_consistently

                #    h1 = (h1 + H1[c,1:nw,1:nw]) / 2.0   #mixing
                #    dq = (dq + DQ[c,1:tbc.crys.nat]) / 2.0
                #h1spin = (h1spin + H1spin[c,:,1:nw,1:nw]) / 2.0
#                H1[c,1:nw,1:nw] =  h1 
#                DQ[c,1:tbc.crys.nat] = dq
#                H1spin[c,:,1:nw, 1:nw] = h1spin


            end

            H1[c,1:nw,1:nw] = get_h1_dq(tbc, dq)

            
#            println("final get eigen")
#            if c == 1
#                println("sum h1 v2 ", sum(abs.(h1)) , "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
#            end
            
            get_eigen(h1, h1spin, tbc.tb.nspin, 2)

            energy_tmp,  efermi = band_energy(VALS_FITTED[c, 1:nk,1:nw,1:tbc.tb.nspin], kweights, nval, 0.01, returnef=true)             
            #            efermi = calc_fermi(VALS_FITTED[c,1:nk,1:nw,1:tbc.tb.nspin], kweights, nval, 0.01)
            occs = gaussian.(VALS_FITTED[c,1:nk,1:nw,1:tbc.tb.nspin].-efermi, 0.01)
            OCCS_FITTED[c,1:nk,1:nw,1:tbc.tb.nspin] = occs

            #            energy_tmp, occs_fitted = band_energy(VALS_FITTED[c,1:nk,1:nw], kweights, tbc.nelec, 0.01, returnocc=true)
            #            OCCS_FITTED[c,1:nk,1:nw] = occs_fitted
            
            energy_smear = smearing_energy(VALS_FITTED[c,1:nk,1:nw,1:tbc.tb.nspin], kweights, efermi, 0.01)
            ENERGY_SMEAR[c] = energy_smear
            
            #            println("$c smear $energy_smear")

            
            #            if scf
            #                energy_charge, pot = ewald_energy(tbc, dq)
            #            else
            #                energy_charge = 0.0
            #            end
            #            
            #            etypes = types_energy(tbc)

            #            s1 = sum(occs .* VALS0_FITTED[c,1:nk,1:nw,1:tbc.tb.nspin], dims=[2,3])
            #            energy_band = sum(s1 .* kweights)
            #            if tbc.tb.nspin == 2
            #                energy_band = energy_band * 2
            #            end
            
            
            if scf
                energy_charge, pot = ewald_energy(tbc, dq)
            else
                energy_charge = 0.0
            end
            if tbc.tb.scfspin
                energy_magnetic = magnetic_energy(tbc.crys, EDEN_FITTED[c, 1:tbc.nspin,1:nw])
                #                if c == 4
                #                    println("magnetic_energy $energy_magnetic")
                #                    println(tbc.crys)
                #                    println(EDEN_FITTED[c, 1:tbc.nspin,1:nw])
                #                    println()
                #                    println("VALS  1", VALS_FITTED[c,1,1:nw, 1])
                #                    println("VALS  2", VALS_FITTED[c,1,1:nw, 2])
                #                    println()
                #                    println("VALS0 1", VALS0_FITTED[c,1,1:nw,1])
                #                    println("VALS0 2", VALS0_FITTED[c,1,1:nw,2])
                #                    println()
                #                end
            else
                energy_magnetic = 0.0
            end



            
            ENERGIES_FITTED[c] = etypes + energy_charge + energy_band + energy_smear + energy_magnetic
            #            println("contr ", round.([etypes , energy_charge , energy_band , energy_smear , energy_magnetic, -99, efermi], digits=8))
            #            println("EDEN 1 ", EDEN_FITTED[c, 1,1:nw])
            #            println("EDEN 2 ", EDEN_FITTED[c, 2,1:nw])
            println("scf $conv $c ", ENERGIES_FITTED[c], " d  ", ENERGIES_FITTED[c] - ENERGIES[c], "    $dq " )
            

        end

        return ENERGIES_FITTED, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED

    end


    if lambda[1] > 1e-10
        NLAM = NCOLS
        #        println("alt lambda ", NLAM, " " , NCOLS)
    else
        NLAM  = 0
    end
    NLAM = Int64(NLAM)


    
    
    function construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED::Array{Float64,4}, ncalc::Int64, ncols::Int64, ncols_S::Int64, nlam::Int64, ERROR::Array{Int64,1}, EDEN_FITTED::Array{Float64,3}; leave_out=-1, cref=missing)
        
#        println("in construct_newXY")
        #        nlam = 0

        #        NEWX = zeros(ncalc*nk_max*NWAN_MAX + ncalc + nlam, ncols)
        #        NEWY = zeros(ncalc*nk_max*NWAN_MAX + ncalc + nlam)

        counter = 0

        temp = 0.0+0.0*im
        
        nonzero_ind = Int64[]
        #count nonzero ahead of time
        for calc = 1:ncalc
            if calc == leave_out #skip this one
                continue
            end
            if ERROR[calc] == 1
                continue
            end
            row1, rowN, nw = ind_BIG[calc, 1:3]
            N = hermetian_indpt(nw)
            vals = ones(nw) * 100.0
            nk = size(KPOINTS[calc])[1]
            for spin = 1:SPIN[calc]
                for k = 1:nk
                    for i = 1:nw
                        counter += 1
                        push!(nonzero_ind, counter)
                    end
                end
            end
            counter += 1
            push!(nonzero_ind, counter)
            if energy_diff_calc == true
                counter += 1
                push!(nonzero_ind, counter)
            end
        end

        NEWX = zeros(counter + nlam, ncols)
        NEWX_S = zeros(counter + nlam, ncols_S)

        NEWY = zeros(counter + nlam)
        
        counter = 0
        energy_counter = []

        CALC_IND = Int64[]
        SPECIAL = Int64[]
        
        (emin, indmin) = findmin( (ENERGIES ./ NAT)[:])
#        println("emin energy_diff_calc $energy_diff_calc $emin $indmin")
        EDIFF_IND = Int64[]
        energy_min_arr_ind = 0
        
        for calc = 1:ncalc
#            println("calc $calc")
            if calc == leave_out #skip this one
                continue
            end
            if ERROR[calc] == 1
                println("ERROR $calc skip")
                continue
            end
#            println("calc2 $calc")

            row1, rowN, nw = ind_BIG[calc, 1:3]
            N = hermetian_indpt(nw)

            vals = ones(nw) * 100.0
            #            H = zeros(Complex{Float64}, nw, nw, ncols)
            #            H = zeros(Complex{Float64}, nw, nw)
            H_cols = zeros(Complex{Float64}, nw, nw, ncols)
            S_cols = zeros(Complex{Float64}, nw, nw, ncols_S)
            
            H_fixed = zeros(Complex{Float64}, nw, nw)
            
            
            #            H_cols = H_COLS[calc]
            
            VECTS = zeros(Complex{Float64}, nw, nw)
            #            S = zeros(Complex{Float64}, nw, nw)
            nk = size(KPOINTS[calc])[1]
            
            X_TOTEN = zeros(ncols)
            SX_TOTEN = zeros(ncols_S)
            Y_TOTEN = ENERGIES[calc]
            
            if scf
                nat = list_of_tbcs[calc].crys.nat
                energy_charge, pot = ewald_energy(list_of_tbcs[calc], DQ[calc,1:nat])
            else
                energy_charge = 0.0
            end
            if list_of_tbcs[calc].tb.scfspin
                energy_magnetic = magnetic_energy(list_of_tbcs[calc], EDEN_FITTED[calc, :, 1:nw])
            else
                energy_magnetic = 0.0
            end
            
            etypes = types_energy(list_of_tbcs[calc].crys)
            
            energy_smear = ENERGY_SMEAR[calc]

            Y_TOTEN -= energy_charge + etypes + energy_smear + energy_magnetic
            
            VECTS = zeros(Complex{Float64}, nw, nw)
            VECTS_p = zeros(Complex{Float64}, nw, nw)

            vals_test_other = zeros(nw, ncols)
            Svals_test_other = zeros(nw, ncols_S)
            vals_test_on = zeros(nw)

            for spin = 1:list_of_tbcs[calc].nspin
                for k = 1:nk
                    for i = 1:nw
                        for j = 1:nw
                            if i <= j
                                #                            H_cols[i,j,:] = X_Hnew_BIG[row1-1 + (k-1)*nw*nw + (i-1)*nw + j,:] + im*X_Hnew_BIG[row1-1 + (k-1)*nw*nw + (i-1)*nw + j + nk*nw*nw,:]
                                ind = hermetian_index(i,j, nw)
                                H_cols[i,j,:] = X_Hnew_BIG[row1-1 + (k-1)*N + ind,:] + im*X_Hnew_BIG[row1-1 + (k-1)*N + ind +  nk*N,:]
                                H_cols[j,i,:] = X_Hnew_BIG[row1-1 + (k-1)*N + ind,:] - im*X_Hnew_BIG[row1-1 + (k-1)*N + ind +  nk*N,:]                            
                                S_cols[i,j,:] = X_Snew_BIG[row1-1 + (k-1)*N + ind,:] + im*X_Snew_BIG[row1-1 + (k-1)*N + ind +  nk*N,:]
                                S_cols[j,i,:] = X_Snew_BIG[row1-1 + (k-1)*N + ind,:] - im*X_Snew_BIG[row1-1 + (k-1)*N + ind +  nk*N,:]                            
                            end
                            
                        end
                    end
                    for i = 1:nw
                        for j = 1:nw
                            if keep_bool
                                if i <= j
                                    ind = hermetian_index(i,j, nw)
                                    H_fixed[i,j] = Xc_Hnew_BIG[row1-1 + (k-1)*N + ind] + im*Xc_Hnew_BIG[row1-1 + (k-1)*N + ind +  nk*N]
                                    H_fixed[j,i] = Xc_Hnew_BIG[row1-1 + (k-1)*N + ind] - im*Xc_Hnew_BIG[row1-1 + (k-1)*N + ind +  nk*N]
                                    
                                end
                                #                            H_fixed[i,j] = Xc_keep[row1-1 + (k-1)*nw*nw + (i-1)*nw + j] + im*Xc_keep[row1-1 + (k-1)*nw*nw + (i-1)*nw + j + nk*nw*nw]
                            end
                        end
                    end


                    #                VECTS[:,:] = VECTS_FITTED[calc,k,1:nw,1:nw]
                    VECTS[:,:] = VECTS_FITTED[calc][spin,k,1:nw,1:nw]
                    VECTS_p[:,:] = VECTS'

                    if keep_bool
                        vals_test_on[:] = real(diag(VECTS_p * (h_on[calc] + H_fixed) * VECTS)) 
                    else
                        vals_test_on[:] = real(diag(VECTS_p * (h_on[calc] ) * VECTS)) 
                    end
                    
                    vals_test_other[:,:] .= 0.0
                    Svals_test_other[:,:] .= 0.0

                    #                @time for i = 1:nw
                    #                    for j = 1:nw
                    #                        for kk = 1:nw
                    #                            for ii = 1:ncols
                    #                                vals_test_other[i,ii] += real(VECTS_p[i,j] .* H_cols[ j,kk,ii]  .* VECTS[kk,i])
                    #                            end
                    #                        end
                    #                    end
                    #                end


                    
                    for ii = 1:ncols
                        for i = 1:nw
                            temp = 0.0+0.0im
                            for kk = 1:nw
                                for j = 1:nw
                                    temp += VECTS_p[i,j] * H_cols[ j,kk,ii]  * VECTS[kk,i]
                                end
                            end
                            vals_test_other[i,ii] += real(temp)
                        end
                    end

                    h1val = zeros(nw)
                    if scf
                        for i = 1:nw
                            for j = 1:nw
                                for k = 1:nw
                                    h1val[i] = real(VECTS_p[i,j] * H1[calc,j,k]  * VECTS[k,i])
                                end
                            end
                        end
                    end                                                
                                

                    
                    for ii = 1:ncols_S
                        for i = 1:nw

                            tempS = 0.0+0.0im
                            for kk = 1:nw
                                for j = 1:nw

                                    tempS += VECTS_p[i,j] * S_cols[ j,kk,ii]  * VECTS[kk,i]
                                end
                            end

                            Svals_test_other[i,ii] += real(tempS)
                        end
                    end
                    
                    for i = 1:nw
                        counter += 1



                        #magic to get lowest energy band correct especially near gamma
                        if k == 1 && nk > 1 && i == 1
                            w_special = 10.0
#                        elseif i == 1
#                            w_special = 2.0
                        else
                            w_special = 1.0
                        end
                        
                        NEWX[counter, 1:ncols] = vals_test_other[i,:] .* WEIGHTS[calc, k, i, spin] * w_special 

                        X_TOTEN[:] +=            vals_test_other[i,:] .* (KWEIGHTS[calc][k] * OCCS_FITTED[calc,k,i, spin])  #* list_of_tbcs[calc].nspin
                        NEWX_S[counter, 1:ncols_S] = VALS_FITTED[calc, k,i,spin] *  WEIGHTS[calc, k, i, spin] * w_special * Svals_test_other[i,:] 

                        #NEWX_S[counter, 1:ncols_S] = VALS_FITTED[calc, k,i,spin] *  WEIGHTS[calc, k, i, spin] * w_special * Svals_test_other[i,:]

                        SX_TOTEN[:] +=   Svals_test_other[i,:] .* (KWEIGHTS[calc][k] * OCCS_FITTED[calc,k,i, spin]) * VALS_FITTED[calc, k,i,spin]
                        
                        #NEWY[counter] =  (VALS[calc,k,i, spin] - vals_test_on[i] - (VALS_FITTED[calc, k,i,spin] - VALS0_FITTED[calc, k,i,spin]) ) .* WEIGHTS[calc, k, i, spin] * w_special
                        #NEWY[counter] =  (VALS0[calc,k,i, spin] - vals_test_on[i] ) .* WEIGHTS[calc, k, i, spin] * w_special

                        vref =@view  VECTS_ref[calc][k, :,:, spin]
                        sref =@view  S_ref[calc][k, :,:]
                        jjj_min = i
                        dist_min = 10^10
#                        println()
                        for jjj = 1:nw
                            dist = (1.0 - abs(VECTS[:,i]'*sref*vref[:,jjj])) + abs(VALS0[calc,k,jjj, spin] - VALS0_FITTED[calc,k,i,spin])

                            if dist < dist_min
                                dist_min = dist
                                jjj_min = jjj
                            end
                            if k == 1 && calc == 1
                                println("i $i jjj $jjj dist $dist dist_min $dist_min jjj_min $jjj_min  $((1.0 - abs(VECTS[:,i]'*sref*vref[:,jjj])))  $(abs(VALS0[calc,k,jjj, spin] - VALS0_FITTED[calc,k,i,spin]))")
                            end
                        end
#                        println()
                        NEWY[counter] =  (VALS0[calc,k,jjj_min, spin] - vals_test_on[i] - h1val[i] ) .* WEIGHTS[calc, k, i, spin] * w_special                         

                        if k == 1 && calc == 1
                            if !ismissing(cref)
                                println("test i $i vals_test_other*c ", vals_test_other[[i],:]*cref, " VALS0 ", VALS0[calc,k,i, spin] - vals_test_on[i] - h1val[i])
                            end
                        end
                        
                        Y_TOTEN += -1.0 * vals_test_on[i] * (KWEIGHTS[calc][k] * OCCS_FITTED[calc,k,i, spin]) #*  list_of_tbcs[calc].nspin

                        #                        if k == 1
                        #                            if i == 1
                        #                                println()
                        #                            end
                        #                            println("typeof $(typeof(NEWX[counter, :])) $typeof(chX[:])")
                        #                            println("test nwan $i calc $calc k $k NEWY $(VALS0[calc,k,i, spin])  NEWX  $(sum(vals_test_other[i,:] .* chX[:]) + vals_test_on[i])       diff $(VALS0[calc,k,i, spin] - (sum(vals_test_other[i,:] .* chX[:]) + vals_test_on[i]))")
                        #                        end
                        #                        println("test test i $i k $k ",  abs(sum(NEWX[counter, :] ) - NEWY[counter]) , " [ $(sum(abs.(NEWX[counter, :] ))), $(NEWY[counter]) ] VALS0 $(VALS0[calc,k,i, spin])")
                        
                        push!(CALC_IND, calc)
                        if calc == 1 && k == 1
                            push!(SPECIAL, counter)
                        end
                        push!(nonzero_ind, counter)

                        ###println("$calc $k $i : ",  vals_test[i], "\t" , VALS_FITTED[calc,k,i],"\t", vals_test_on[i] + vals_test_other[i,:]'*ch)
                    end
                    
                    
                end
#                println("calc3 $calc")
            end
            counter += 1
            NEWX[counter, :] = X_TOTEN[:] * energy_weight * weights_list[calc] / NAT[calc]
            NEWX_S[counter, :] = SX_TOTEN[:] * energy_weight * weights_list[calc] / NAT[calc]

            NEWY[counter] = Y_TOTEN * energy_weight * weights_list[calc]  / NAT[calc]
            push!(CALC_IND, calc)
            if calc == 1 
                push!(SPECIAL, counter)
            end

            
            push!(nonzero_ind, counter)
            push!(energy_counter, counter)

            if energy_diff_calc == true
                counter += 1
                NEWX[counter, :] = X_TOTEN[:] * energy_weight / NAT[calc]
                NEWX_S[counter, :] = SX_TOTEN[:] * energy_weight / NAT[calc]
                NEWY[counter] = Y_TOTEN * energy_weight / NAT[calc]
                push!(CALC_IND, calc)
                push!(nonzero_ind, counter)
                push!(EDIFF_IND, counter)
                if calc == indmin
                    energy_min_arr_ind = counter
                end
            end
            
            
        end


        if energy_diff_calc
            if indmin != 0 && energy_min_arr_ind != 0
                for (calc,ind) in enumerate(EDIFF_IND)
                    NEWX[ind, :] = ( NEWX[ind, :]  - NEWX[energy_min_arr_ind, :] ) * energy_diff_weight  * weights_list[calc]
                    NEWX_S[ind, :] = ( NEWX_S[ind, :] - NEWX_S[energy_min_arr_ind, :] ) * energy_diff_weight  * weights_list[calc]
                    NEWY[ind] = ( NEWY[ind]  - NEWY[energy_min_arr_ind] ) * energy_diff_weight * weights_list[calc]
                end
            end
        end

#=        
        if lambda[1] > 1e-10 || lambda[2] > 1e-10
            for ind3 = 1:nlam
                counter += 1
                if !(ind3 in threebody_inds)
        NEWX[counter,ind3] = lambda[1]
                else
                    NEWX[counter,ind3] = lambda[2]
                end
                push!(nonzero_ind, counter)
            end
            for ind3 = 1:size(NEWX_S,2)
                counter += 1
                NEWX_S[counter,ind3] = lambda[1]
            end
        end
=#
        #        println("len nonzero_ind ", length(nonzero_ind))
        #        println("size NEWX old ", size(NEWX))
        #        NEWX = NEWX[nonzero_ind,:]
        #        NEWY = NEWY[nonzero_ind]
        #        println("size NEWX new ", size(NEWX))
        #        println("size NEWY new ", size(NEWY))

        #        println("c newX neWy ", NEWX \ NEWY)
        #        println("test ", sum(abs.(NEWX * ones(size(NEWX)[2]) - NEWY)))

        #        println("size NEWX end ", size(NEWX))
        return NEWX, NEWY, NEWX_S, energy_counter, CALC_IND, SPECIAL

    end


    
    function error_fn(NEWX, NEWY, ch, cs)
        error = sum( (NEWX *ch - NEWY).^2)
        not_threebody_inds = setdiff(threebody_inds, 1:size(NEWX)[2])
#        error += sum(ch[threebody_inds].^2 * lambda[2]) + sum(ch[not_threebody_inds].^2 * lambda[1])
#        error += sum(cs.^2 * lambda[1])
        return error
        #return cs[1]^2
    end

    optim_S = true
    function grad_error_fn(NEWX, NEWY, NEWX_S, ch, cs)
        nh = length(ch)
        ns = length(cs)
        grad = zeros(nh + ns)
        temp = (NEWX *ch - NEWY)
        for i = 1:size(NEWX)[1]
            grad[1:nh] += (2 * temp[i] * NEWX[i,:])[:]
        end

        if optim_S
            for i = 1:size(NEWX)[1]
                grad[nh .+ (1:ns)] += -(2 * temp[i] * NEWX_S[i,:])[:]
            end
        end
        println("grad $grad")

#        not_threebody_inds = setdiff(threebody_inds, 1:nh)
#        grad[threebody_inds] += 2 * ch[threebody_inds] * lambda[2]
#        grad[not_threebody_inds] += 2 * ch[not_threebody_inds] * lambda[1]
#        grad[nh .+ (1:ns)] += -2 * cs * lambda[1]
        #        return grad
#        grad[5] = 2*cs[1]
        return grad
        
    end
    #solve_scf_mode=scf
    
    
#    ch[1] = ch[1] * 0.95
#    cs[1] = cs[1] * 0.95


    nh = length(ch)
    ns = length(cs)
    
    ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
    NEWX, NEWY, NEWX_S, energy_counter, CALC_IND, SPECIAL = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
    buffer = [ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED, NEWX, NEWY, NEWX_S, energy_counter, CALC_IND]
    
    function calculate_common_part!(chcs, last_chcs, buffer)
        if sum(abs.(chcs - last_chcs)) > 1e-12
            @suppress begin
                ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
                NEWX, NEWY, NEWX_S, energy_counter, CALC_IND, SPECIAL = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
                copy!(buffer, [ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED,NEWX, NEWY, NEWX_S, energy_counter, CALC_IND])
            end
        end
    end

    function f(x,buffer, last_x)
        calculate_common_part!(x, last_x, buffer)
        ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED, NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = buffer
        ch = x[1:nh]
        cs = x[nh .+ (1:ns)]
        err = error_fn(NEWX, NEWY, ch, cs)
        println("err $err  x $x")
        return err
    end
    function g!(x,stor,buffer,last_x)
        calculate_common_part!(x, last_x, buffer)
        ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED, NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = buffer
        ch = x[1:nh]
        cs = x[nh .+ (1:ns)]
        stor[:] = grad_error_fn(NEWX, NEWY, NEWX_S, ch, cs)
    end


    function iterate(ch, cs, updateH, updateS, niters, mix_iter, mix_iterS; adjust_mix=false)
        println("ch 1 ", ch[1] , " sum DQ ", sum(abs.(DQ)))
        scf = true
        solve_self_consistently = false
        solve_scf_mode = false
        println("iterate solve_scf_mode $solve_scf_mode scf $scf")
        err1 = -99.0
        err2 = -99.0
#        mix_iterS = mix_iter
        cs_old = deepcopy(cs)
        for iter = 1:niters

            err_old = error_fn(NEWX, NEWY, ch, cs)

            
            if updateS
                cs_old[:] = cs[:]
                if iter==1
                    iter_s = 2
                else
                    iter_s = 1
                end
                for i = 1:iter_s
                    @suppress begin
                        ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
                        NEWX, NEWY, NEWX_S, energy_counter, CALC_IND, SPECIAL = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
                    end
                    
                    
                    cst = NEWX_S \ ( (NEWX*ch - NEWY ))
                    #cs = (cs + cst)
                    cs_test = (cs + cst) * mix_iterS + cs*(1 - mix_iterS)

                    cs = cs_test
                    
                end
                    
            end
                
            err1 = error_fn(NEWX, NEWY, ch, cs)
            if err1 > err_old  && adjust_mix
                mix_iterS = mix_iterS * 0.8
                mix_iterS = max(mix_iterS, 0.05)
            else
                mix_iterS = mix_iterS * 1.05
            end
            mix_iterS = min(0.5, mix_iterS)

            if updateH
                if iter==1
                    iter_h = 1
                else
                    iter_h = 1
                end
                SPECIAL = []
                for i in 1:iter_h
                    @suppress begin
                        ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
                    end
                        NEWX, NEWY, NEWX_S, energy_counter, CALC_IND, SPECIAL = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out, cref = ch)
                        #println("size NEWX ", size(NEWX))


                    println()
                    println("vals calc ", VALS_FITTED[1,1,:,1])
                    println()
                    
                    ch_new = NEWX \ NEWY

                    #println("ch_new size ", size(ch_new))
                    ch_test = ch_new*mix_iter + ch*(1 - mix_iter)
                    if length(SPECIAL) > 0
                        y = NEWX*ch_test
                        println("sp     y    NEWY")
                        for s in SPECIAL
                            println("special $s  $(y[s])   $(NEWY[s])")
                        end
                        println()
                    end

                    #err_test = error_fn(NEWX, NEWY, ch_test, cs)
                    #println("err_test $err_test")
                    ch = ch_test
                    #                    err_test = error_fn(NEWX, NEWY, ch, cs)
                    #println("err_test2 $err_test")
                    
                end
            end
            err2 = error_fn(NEWX, NEWY, ch, cs)
            if err2 > err1  && adjust_mix
                mix_iter = mix_iter * 0.5
            else
                mix_iter = mix_iter * 1.1
            end
            mix_iter = min(mix_iter, 0.5)
            
#            println("size NEWX ", size(NEWX))
#            println("size NEWY ", size(NEWY))
#            println("size ch ", size(ch))
#            println("size cs ", size(cs))
            err2 = error_fn(NEWX, NEWY, ch, cs)
            println("iter $iter   $err1   $err2  err_old $err_old   ch1 $(ch[1]) cs1 $(cs[1]) mix_iter $mix_iter mix_iterS $mix_iterS")
            flush(stdout)
            if updateS == false
                if abs(err_old - err2) < 1e-4 && iter > 2
                    println("break, no progress")
                    break
                end
                if abs(err_old - err2) < 1e-3 && iter > 5
                    println("break, no progress")
                    break
                end
            else
                if abs(err_old - err2) < 1e-3 && iter > 5 && sum(abs.(cs_old - cs)) < 0.1
                    println("break, no progress")
                    break
                end
            end
            
            if updateS
                if err2 > err_old && iter > 4
                    println("break, update S, negative progress")                    
                    break
                end
            else
                if err2 > err_old && iter > 2
                    println("break, update H only, negative progress")                                        
                    break
                end
            end
            err_old = err2
            
        end
        println("iterate_fn_end $err1 $err2")
        return ch, cs, err2, mix_iter, mix_iterS
    end

    if true
#        println("if true")
#        println("NCOLS, NCOLS_S ", [NCOLS, NCOLS_S])
        println("opt_S $opt_S")
        #if opt_S
            #ch, cs = iterate(ch, cs, true, false, 200, 0.05, adjust_mix=true)
            #ch, cs = iterate(ch, cs, true, true, 550, 0.2, adjust_mix=false)
            #ch, cs = iterate(ch, cs, true, true, 50, 0.02, adjust_mix=false)
        #else

        if true
            solve_scf_mode=true
            scf = true
            begin 
                for init_scf = 1:5
                    DQ_old = deepcopy(DQ)
                    println("DQ")
                    println(DQ)
                    ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
                    if maximum(abs.(DQ_old - DQ)) < 0.01
                        println("break scf early $init_scf $(maximum(abs.(DQ_old - DQ)))")
                        break
                    end
                end
                #NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
                solve_scf_mode=false
                println("solve_scf_mode $solve_scf_mode scf $scf")
                ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
                NEWX, NEWY, NEWX_S, energy_counter, CALC_IND, SPECIAL = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
                println()
                println("final DQ")
                println(DQ)
                
            end
            solve_scf_mode=false

            err = error_fn(NEWX, NEWY, ch, cs)
            println("start err")
            
            
            
#            solve_scf_mode=false
#            ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
#            NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)

            solve_scf_mode=false


            println("DQ BEFORE")
            println(DQ)
            println()


            scf = true
            mix_iterS = 0.005
            mix_iter = 0.004
            err_old_bigiter = 10.0^10
            err = 10.0^9.0
            for big_iter = 1:40
                println("BIG ITER $big_iter solve_scf_mode $solve_scf_mode scf $scf sum DQ $(sum(abs.(DQ))) ---------------------------------------------------------------------------------------------------------------------------- $err_old_bigiter")
                println()

                if false
                @suppress begin
                list_of_tbcs,KPOINTS, KWEIGHTS, dft_list, scf, energy_weight, rs_weight, ks_weight, weights_list, NWAN_MAX, NCALC, VALS, VALS0, E_DEN, H1, H1spin, DQ, ENERGY_SMEAR, OCCS, WEIGHTS, ENERGIES, NCOLS_orig, NCOLS, CH, NVAL, NAT, SPIN_MAX, Ys, keep_bool, keep_inds, toupdate_inds, ch_keep, keep_inds_S, toupdate_inds_S, cs_keep =
                    topstuff_direct(deepcopy(list_of_tbcs_nonscf), prepare_data; EDEN_input = E_DEN, weights_list=weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight, ks_weight = ks_weight, niters=niters, lambda=lambda,  leave_one_out=false, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max = nk_max)
                end
                end

               
#                println("aaaa solve_scf_mode $solve_scf_mode scf $scf")
                

                solve_scf_mode=false
                if opt_S
                    ch, cs, err, mix_iter, mix_iterS = iterate(ch, cs,true , opt_S, 30, max(0.02, mix_iter * 0.8), mix_iterS , adjust_mix=true)
                else
                    ch, cs, err, mix_iter, mix_iterS = iterate(ch, cs,true , opt_S, 30, max(0.02, mix_iter * 0.8), max(0.02, mix_iterS * 0.8) , adjust_mix=true)
                end

 #               println("DQ AFTER ITERATE")
 #               println(DQ)
 #               println()


                H1_bigiter = deepcopy(H1)
                DQ_bigiter = deepcopy(DQ)
#                DQ_EDEN_bigiter = deepcopy(DQ_EDEN)
                H1spin_bigiter = deepcopy(H1spin)
                EDEN_bigiter = deepcopy(E_DEN)
#                println("a solve_scf_mode $solve_scf_mode scf $scf")

                solve_scf_mode=true
                for init_scf = 1:5
                    DQ_old = deepcopy(DQ)
                    #println("DQ")
                    #println(DQ)
                    @suppress begin
                        ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
                    end
                    if maximum(abs.(DQ_old - DQ)) < 0.01
                        println("break scf early $init_scf $(maximum(abs.(DQ_old - DQ)))")
                        break
                    end
                end
                ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)

                NEWX, NEWY, NEWX_S, energy_counter, CALC_IND, SPECIAL = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
                #err = error_fn(NEWX, NEWY, ch, cs)
                #println("err X ", err)
                
                
                mix = 0.5
                #H1 = (1-mix)*H1 + mix*H1_bigiter 
                DQ = (mix)*DQ + (1-mix)*DQ_bigiter
                println("sum DQ after mix  $(sum(abs.(DQ))) ")
#                println(DQ)
#                println()
                E_DEN = (mix)*EDEN_FITTED + (1-mix)*EDEN_bigiter

                #                DQ_EDEN = (1-mix)*DQ_EDEN + mix*DQ_EDEN_bigiter
                #H1spin = (1-mix)*H1spin + mix*H1spin_bigiter

                #solve_scf_mode = false
#                println()
#                println("now false  solve_scf_mode $solve_scf_mode scf $scf")
                
                #ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
#                println("e solve_scf_mode $solve_scf_mode scf $scf")
                #NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
#                println("f solve_scf_mode $solve_scf_mode scf $scf")                
                #err = error_fn(NEWX, NEWY, ch, cs)
                #println("err Y ", err)
                #solve_scf_mode=false
                if abs(err_old_bigiter - err) < 5e-3 && big_iter >= 3
                    println("done, bigiter break err_old $err_old_bigiter err $err diff $(abs(err_old_bigiter - err))")
                    break
                end
                err_old_bigiter = err
            end
            
            #            ch, cs = iterate(ch, cs, true, false, 50, 0.05, adjust_mix=true)
#            ch, cs = iterate(ch, cs, true, false, 50, 0.02, adjust_mix=false)


            #ch, cs = iterate(ch, cs, true, false, 150, 0.2, adjust_mix=false)
            #ch, cs = iterate(ch, cs, true, false, 50, 0.02, adjust_mix=false)
        end
        println("DQ AFTER")
        println(DQ)
        println()
        println("solve_scf_mode $solve_scf_mode scf $scf")
        
        ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
        NEWX, NEWY, NEWX_S, energy_counter, CALC_IND, SPECIAL = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)

        good =  (abs.(ENERGIES - ENERGIES_working) ./ NAT) .< 0.05
        good = good .&& .!(Bool.(ERROR))
        

        ####
        ch_big = zeros(length(keep_inds) + length(toupdate_inds))
        #println(length(keep_inds), " " , length(toupdate_inds))

        ch_big[toupdate_inds] = ch[:]
        ch_big[keep_inds] = ch_keep[:]
        chX2 = ch_big
        
        cs_big = zeros(length(keep_inds_S) + length(toupdate_inds_S))
        cs_big[toupdate_inds_S] = cs[:]
        cs_big[keep_inds_S] = cs_keep[:]
        csX2 = cs_big
        ####

        if cg == false
            database = make_database(chX2, csX2,  KEYS, HIND, SIND,DMIN_TYPES,DMIN_TYPES3, scf=scf, starting_database=starting_database, tbc_list = list_of_tbcs[good], fit_eam=fit_eam, fitting_version=fitting_version)
            return database
        end
    end
#    T = []
#    TG = []
#    for x = -0.1:0.02:0.5
#        cs[1] = x
#        ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
#        NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
#        err = error_fn(NEWX, NEWY, ch, cs)
#        ggg = grad_error_fn(NEWX, NEWY, NEWX_S, ch, cs)
#        println("x $x cs $cs   err $err")
#        push!(T, err)
#        push!(TG, ggg[5])
#    end
#    plot(-0.1:0.02:0.5, T)
#    plot!(-0.1:0.02:0.5, TG)
#    x = -0.1:0.02:0.5
    
#    plot!( 0.5*(x[1:end-1] + x[2:end]), (T[2:end] .- T[1:end-1])/(0.02) , linestyle=:dash)
    

    if cg
        opts = Optim.Options(g_tol = 1e-5, f_abstol = 1e-5, x_abstol = 1e-5, f_calls_limit = 200,
                             iterations = 200,                                                                                                                                                                
                             store_trace = true,                                                                                                                                                                 
                             show_trace = false)                                                                                                                                                                 


        optim_S = true
        initial_x = vcat(chX2,csX2)
        stor = zeros(size(initial_x))
        last_x = similar(initial_x)             

        f0 = f(initial_x, buffer, initial_x)
        g0 = g!(initial_x, stor, buffer, initial_x)

        println("f0 $f0")
        println("g0 $g0")
        
        df = OnceDifferentiable(x -> f(x, buffer, initial_x),(stor, x) -> g!(x, stor, buffer, last_x), deepcopy(initial_x))
        #ret = optimize(df, initial_x, LBFGS(linesearch = BackTracking(order=3)))
        #ret = optimize(df, initial_x, LBFGS(linesearch = BackTracking()))
        ret = optimize(df, initial_x, LBFGS(linesearch = MoreThuente()) , opts)
        #ret = optimize(df, initial_x, LBFGS(linesearch = HagerZhang()) , opts)

        println("RESTART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
#        sleep(0.5)
        optim_S = true
        initial_x = ret.minimizer
        stor = zeros(size(initial_x))
        last_x = similar(initial_x)             

        f0 = f(initial_x, buffer, initial_x)
        g0 = g!(initial_x, stor, buffer, initial_x)

        println("f0 $f0")
        println("g0 $g0")
        
        df = OnceDifferentiable(x -> f(x, buffer, initial_x),(stor, x) -> g!(x, stor, buffer, last_x), deepcopy(initial_x))
        #ret = optimize(df, initial_x, LBFGS(linesearch = BackTracking()))
        #ret = optimize(df, initial_x, LBFGS(linesearch = BackTracking(order=3)))
        ret = optimize(df, initial_x, LBFGS(linesearch = MoreThuente()) , opts) 
        #ret = optimize(df, initial_x, LBFGS(linesearch = HagerZhang()) , opts)
       

        
        #ret = optimize(df, initial_x, ConjugateGradient())
        #ret = optimize(df, initial_x, LBFGS(), opts)

        println("RESTART 222222222222222222 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
#        sleep(0.5)
        optim_S = true
        initial_x = ret.minimizer
        stor = zeros(size(initial_x))
        last_x = similar(initial_x)             

        f0 = f(initial_x, buffer, initial_x)
        g0 = g!(initial_x, stor, buffer, initial_x)

        println("f0 $f0")
        println("g0 $g0")
        
        df = OnceDifferentiable(x -> f(x, buffer, initial_x),(stor, x) -> g!(x, stor, buffer, last_x), deepcopy(initial_x))
        #ret = optimize(df, initial_x, LBFGS(linesearch = BackTracking()))
        #ret = optimize(df, initial_x, LBFGS(linesearch = BackTracking(order=3)))
        ret = optimize(df, initial_x, LBFGS(linesearch = MoreThuente()), opts ) 
        #ret = optimize(df, initial_x, LBFGS(linesearch = HagerZhang()) , opts)
       
        ch = ret.minimizer[1:length(ch)]
        cs = ret.minimizer[length(ch) .+ (1:length(cs))]

        ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
        NEWX, NEWY, NEWX_S, energy_counter, CALC_IND, SPECIAL = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)

        good =  (abs.(ENERGIES - ENERGIES_working) ./ NAT) .< 0.05
        good = good .&& .!(Bool.(ERROR))
        
        database = make_database(ch, cs,  KEYS, HIND, SIND,DMIN_TYPES,DMIN_TYPES3, scf=scf, starting_database=starting_database, tbc_list = list_of_tbcs[good], fit_eam=fit_eam, fitting_version=fitting_version)
        return database

    end


    
#    return ret
    
    ####test stuff    

    if false
        ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
        NEWX, NEWY, NEWX_S, energy_counter, CALC_IND, SPECIAL = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)

        println("ch $ch")
        println("cs $cs")
        println()
        err0 = error_fn(NEWX, NEWY, ch, cs)
        grad0 = grad_error_fn(NEWX, NEWY, NEWX_S, ch, cs)
        println("err0 $err0")

        println()
        err_vec = []
        dx = 1e-7
        for i = 1:length(ch)
            println("ch i $i")
            cht = deepcopy(ch)
            cht[i] += dx
            @suppress begin 
                @time ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(cht, cs, solve_scf_mode)
                @time  NEWX, NEWY, NEWX_S, energy_counter, CALC_IND, SPECIAL = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
            end
            err = error_fn(NEWX, NEWY, cht, cs)
            push!(err_vec, err)
            println("err $err")
            grad = grad_error_fn(NEWX, NEWY, NEWX_S, cht, cs)
            println("grad $grad")
        end



        for i = 1:length(cs)
            println("cs i $i")
            cst = deepcopy(cs)
            cst[i] += dx
            @suppress begin 
                @time ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cst, solve_scf_mode)
                @time  NEWX, NEWY, NEWX_S, energy_counter, CALC_IND, SPECIAL = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
            end
            err = error_fn(NEWX, NEWY, ch, cst)
            push!(err_vec, err)
            println("err $err")
            grad = grad_error_fn(NEWX, NEWY, NEWX_S, ch, cst)
            println("grad $grad")
        end
        println()
        println("test ", (err_vec .- err0) / dx)
        println("grad0 ", grad0)
    end    
    
    #database = make_database(ch, cs,  KEYS, HIND, SIND,DMIN_TYPES,DMIN_TYPES3, scf=scf, starting_database=starting_database, tbc_list = list_of_tbcs[good], fit_eam=fit_eam)

#    return database
    
end


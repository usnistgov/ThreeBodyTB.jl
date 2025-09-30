using Optim
using LineSearches

function do_fitting_direct(list_of_tbcs ; weights_list = missing, dft_list=missing, kpoints = missing, starting_database = missing,  update_all = false, fit_threebody=true, fit_threebody_onsite=true, do_plot = false, energy_weight = missing, rs_weight=missing,ks_weight=missing, niters=50, lambda=[0.0,0.0], leave_one_out=false, prepare_data = missing, RW_PARAM=0.0, NLIM = 100, refit_database = missing, start_small = false, fit_to_dft_eigs=false, fit_eam=false, ch_startX = missing, energy_diff_calc = false, gen_add_ham=false, fitting_version = fitting_version_default, opt_S = true, conjgrad=false)

    println("do_fitting_direct   niters $niters update_all $update_all fit_threebody $fit_threebody fit_threebody_onsite $fit_threebody_onsite  energy_weight $energy_weight  rs_weight $rs_weight ks_weight $ks_weight lambda $lambda RW_PARAM $RW_PARAM NLIM $NLIM fit_eam $fit_eam energy_diff_calc $energy_diff_calc opt_S $opt_S ")
    
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

    return do_fitting_direct_main(list_of_tbcs, pd; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = ks_weight, niters=niters, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam, ch_startX = ch_startX, energy_diff_calc = energy_diff_calc, gen_add_ham=gen_add_ham, fitting_version=fitting_version, opt_S = opt_S, cg = conjgrad)
end

function do_fitting_direct_main(list_of_tbcs, prepare_data; weights_list=missing, dft_list=missing, kpoints = missing, starting_database = missing,  update_all = false, fit_threebody=true, fit_threebody_onsite=true, do_plot = false, energy_weight = missing, rs_weight=missing, ks_weight = missing, niters=50, lambda=[0.0, 0.0], leave_one_out=false, RW_PARAM=0.0001, KPOINTS=missing, KWEIGHTS=missing, nk_max=0, start_small=false, fit_to_dft_eigs=false, fit_eam=false, optimS = false, top_vars = missing, ch_startX = missing, energy_diff_calc = false, gen_add_ham=false, fitting_version=fitting_version_default, opt_S=true, cg = false)

    leave_out = -1
    
    println("start do_fitting_recursive_direct")
    if typeof(lambda) <: Float64
        lambda = [lambda, lambda]
    end
    
    database_linear, ch_lin, cs_lin, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, Ys_new, cs, ch_refit, SPIN, threebody_inds  = prepare_data

    println("ch_lin ", ch_lin)
    println("cs_lin ", cs_lin)
#    return missing
    
    if ismissing(top_vars)
        list_of_tbcs,KPOINTS, KWEIGHTS, dft_list, scf, energy_weight, rs_weight, ks_weight, weights_list, NWAN_MAX, NCALC, VALS, VALS0, E_DEN, H1, H1spin, DQ, ENERGY_SMEAR, OCCS, WEIGHTS, ENERGIES, NCOLS_orig, NCOLS, ch, NVAL, NAT, SPIN_MAX, Ys, keep_bool, keep_inds, toupdate_inds, ch_keep, keep_inds_S, toupdate_inds_S, cs_keep =
            topstuff(list_of_tbcs, prepare_data; weights_list=weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight, ks_weight = ks_weight, niters=niters, lambda=lambda,  leave_one_out=false, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max = nk_max)
    else
        list_of_tbcs,KPOINTS, KWEIGHTS, dft_list, scf, energy_weight, rs_weight, ks_weight, weights_list, NWAN_MAX, NCALC, VALS, VALS0, E_DEN, H1, H1spin, DQ, ENERGY_SMEAR, OCCS, WEIGHTS, ENERGIES, NCOLS_orig, NCOLS, ch, NVAL, NAT, SPIN_MAX, Ys, keep_bool, keep_inds, toupdate_inds, ch_keep, keep_inds_S, toupdate_inds_S, cs_keep = top_vars
    end
    Ys = Ys_new + Xc_Snew_BIG

    NCOLS_S = size(X_Snew_BIG)[2]

    println("start ncols ", [NCOLS, NCOLS_S]  , " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    
    weights_list = weights_list ./ maximum(weights_list)

    if !ismissing(ch_startX)
        ch[1:length(ch_startX)] = ch_startX
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
    
    function construct_fitted(ch, cs, solve_self_consistently = false)

        println("construct_fitted ch $ch cs $cs")
        
        Xc = (X_Hnew_BIG * ch) + Xc_Hnew_BIG
        Ys = X_Snew_BIG * cs
        
        #        VECTS_FITTED     = zeros(Complex{Float64}, NCALC, nk_max, NWAN_MAX, NWAN_MAX)
        VECTS_FITTED = Dict{Int64, Array{Complex{Float64},4} }()
        VALS_FITTED     = zeros(NCALC, nk_max, NWAN_MAX, SPIN_MAX)
        VALS0_FITTED     = zeros(NCALC, nk_max, NWAN_MAX, SPIN_MAX)
        OCCS_FITTED     = zeros(NCALC, nk_max, NWAN_MAX, SPIN_MAX)
        ENERGIES_FITTED = zeros(NCALC)
        EDEN_FITTED = zeros(NCALC, SPIN_MAX, NWAN_MAX)
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
                h1 = deepcopy(H1[c,1:nw,1:nw])
                dq = deepcopy(DQ[c,1:tbc.crys.nat])

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
                            if iter == 1
                                if maximum(abs.(imag.(vals))) > 1e-10
                                    println("WARNING, S has negative eigenvalues $c $k valsS $valsS , set error to true")
                                    ERROR[c] = 1
                                    continue
                                end
                            end 
                            
                            #                            println("get eigen vals $vals")
#                                                        if c == 1 && k == 1
#                                                            println("vals k $k   $vals")
#                                                            println()
#                                                        end
                        catch
                            
                            ERROR[c] = 1
                            vals = zeros(Float64, nw)
                            vects = zeros(Complex{Float64}, nw, nw)
                            println("WARNING, S has negative eigenvalues $c $k")
                            
                        end

                        #                    VECTS_FITTED[c,k,1:nw,1:nw] = vects

                        

                        VECTS[spin, k,1:nw,1:nw] = vects
                        VALS_FITTED[c,k,1:nw, spin] = real(vals)
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
            conv = false

            
            
            for c_scf = 1:niter_scf
                #                println("$c_scf aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")                

                #println("c_scf $c_scf get_eigen")
                get_eigen(h1, h1spin, tbc.tb.nspin, c_scf)
                
                
                #                if c_scf == 1
                #                    println("$c VALS_FITTED k0 ", VALS_FITTED[c,1,1:nw])
                #                    println("$c VALS0_FITTED k0 ", VALS0_FITTED[c,1,1:nw])
                #                end

                #if solve_self_consistently == true
                #println("other stuff")
                if true
                    
                    mix = 0.05 * 0.95^c_scf  #start aggressive, reduce mixing slowly if we need most iterations

                    
                    #                    energy_tmp,efermi   = band_energy(VALS_FITTED[c,1:nk,1:nw], kweights, tbc.nelec, 0.01, returnef=true) 

                    #efermi = calc_fermi(VALS_FITTED[c,1:nk,1:nw, 1:tbc.tb.nspin], kweights, nval, 0.01)
                    energy_tmp,  efermi = band_energy(VALS_FITTED[c, 1:nk,1:nw,1:tbc.tb.nspin], kweights, nval, 0.01, returnef=true) 

                    occs = gaussian.(VALS_FITTED[c,1:nk,1:nw, 1:tbc.tb.nspin].-efermi, 0.01)

                    #                    println("occs late $efermi $nval ",occs)
                    #                    println("valsf ", size(VALS_FITTED), " " , size(VALS_FITTED[c,1:nk,1:nw, 1:tbc.tb.nspin]))
                    
                    energy_smear = smearing_energy(VALS_FITTED[c,1:nk,1:nw, 1:tbc.tb.nspin], kweights, efermi, 0.01)

                    #                    eden, h1_new, dq_new = get_electron_density(tbc, kpoints, kweights, VECTS_FITTED[c,:,1:nw,1:nw], occs, Smat)         #updated h1

                    eden, h1_new, dq_new, h1spin_new = get_electron_density(tbc, kpoints, kweights, VECTS_FITTED[c], occs, Smat)         #updated h1

                    #                    println("EDEN DDDDDDDDDD ", eden)
                    
                    EDEN_FITTED[c, 1:tbc.nspin,1:nw] = eden
                    
                    s1 = sum(occs .* VALS0_FITTED[c,1:nk,1:nw, 1:tbc.tb.nspin], dims=[2,3])
                    energy_band = sum(s1 .* kweights) #/ tbc.tb.nspin
                    #                    println("energy_band ", occs[1] , " " ,  VALS0_FITTED[1], " ", kweights[1])
                    
                    if maximum(abs.(dq - dq_new)) > 0.1
                        mix = mix * 0.5
                    end

#                    h1 = h1*(1-mix) + h1_new * mix
#                    dq = dq*(1-mix) + dq_new * mix
#                    h1spin = h1spin*(1-mix) + h1spin_new * mix

                    if solve_self_consistently == true
                    #if false

                        dq = dq*(1-mix) + dq_new * mix
                        h1 = h1*(1-mix) + h1_new * mix
#                        h1spin = h1spin*(1-mix) + h1spin_new * mix

                        energy_charge, pot = ewald_energy(tbc, dq)
                        if tbc.tb.scfspin
                            energy_magnetic = magnetic_energy(tbc, eden)
                        else
                            energy_magnetic = 0.0
                        end

                    else
                        h1 .= 0.0
                        h1spin .= 0.0
                        energy_charge = 0.0
                        energy_magnetic = 0.0
                    end
                    
                    
                    energy_new = energy_charge + energy_band + energy_smear + energy_magnetic + etypes
                    #                    if c == 34
                    #                        println( " scf $c_scf $c ", energy_new+etypes, "    $dq   $energy_charge $energy_band $energy_smear $energy_magnetic")
                    #                    end
                    
                    #                    println( " scf $c_scf $c ", energy_new+etypes, "    $dq   $energy_charge $energy_band $energy_smear $energy_magnetic")

#                    println("c $c   $c_scf energy_new $energy_new energy_old $energy_old diff $(abs(energy_new  - energy_old))  dq $dq")
                    if (abs(energy_new  - energy_old) < 1e-6 && c_scf >= 2) ||  (abs(energy_new  - energy_old) < 1e-5 && c_scf >= 20) ||  (abs(energy_new  - energy_old) < 1e-4 && c_scf >= 50)
                        #                        println("scf converged  $c ", energy_new+etypes, "    $dq " )
                        conv = true
                        break
                    end
                    energy_old = energy_new
                    
                end
                
            end #scf loop

            if scf && conv
                #                h1 = (h1 + H1[c,1:nw,1:nw]) / 2.0   #mixing
                #                dq = (dq + DQ[c,1:tbc.crys.nat]) / 2.0

                H1[c,1:nw,1:nw] =  h1 
                DQ[c,1:tbc.crys.nat] = dq
                H1spin[c,:,1:nw, 1:nw] = h1spin
                
            elseif scf

                h1 = (h1 + H1[c,1:nw,1:nw]) / 2.0   #mixing
                dq = (dq + DQ[c,1:tbc.crys.nat]) / 2.0
                h1spin = (h1spin + H1spin[c,:,1:nw,1:nw]) / 2.0
                H1[c,1:nw,1:nw] =  h1 
                DQ[c,1:tbc.crys.nat] = dq
                H1spin[c,:,1:nw, 1:nw] = h1spin


            end


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


    
    
    function construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED::Array{Float64,4}, ncalc::Int64, ncols::Int64, ncols_S::Int64, nlam::Int64, ERROR::Array{Int64,1}, EDEN_FITTED::Array{Float64,3}; leave_out=-1)

        println("in construct_newXY")
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

        (emin, indmin) = findmin( (ENERGIES ./ NAT)[:])
        println("emin energy_diff_calc $energy_diff_calc $emin $indmin")
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
                        if k == 1 && nk > 1 && i <= 2
                            w_special = 10.0
                        elseif i == 1
                            w_special = 2.0
                        else
                            w_special = 1.0
                        end
                        
                        NEWX[counter, 1:ncols] = vals_test_other[i,:] .* WEIGHTS[calc, k, i, spin] * w_special
                        X_TOTEN[:] +=            vals_test_other[i,:] .* (KWEIGHTS[calc][k] * OCCS_FITTED[calc,k,i, spin])  #* list_of_tbcs[calc].nspin
                        NEWX_S[counter, 1:ncols_S] = VALS_FITTED[calc, k,i,spin] *  WEIGHTS[calc, k, i, spin] * w_special * Svals_test_other[i,:]
                        #NEWX_S[counter, 1:ncols_S] = VALS_FITTED[calc, k,i,spin] *  WEIGHTS[calc, k, i, spin] * w_special * Svals_test_other[i,:]
                        SX_TOTEN[:] +=   Svals_test_other[i,:] .* (KWEIGHTS[calc][k] * OCCS_FITTED[calc,k,i, spin]) * VALS_FITTED[calc, k,i,spin]
                        
                        NEWY[counter] =  (VALS[calc,k,i, spin] - vals_test_on[i] - (VALS_FITTED[calc, k,i,spin] - VALS0_FITTED[calc, k,i,spin]) ) .* WEIGHTS[calc, k, i, spin] * w_special

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
        return NEWX, NEWY, NEWX_S, energy_counter, CALC_IND

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
    solve_scf_mode=scf
    
#    ch[1] = 0.9
#    cs[1] = 9.0
#    cs[1] = 0.9

    nh = length(ch)
    ns = length(cs)
    
    ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
    NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
    buffer = [ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED, NEWX, NEWY, NEWX_S, energy_counter, CALC_IND]
    
    function calculate_common_part!(chcs, last_chcs, buffer)
        if sum(abs.(chcs - last_chcs)) > 1e-12
            @suppress begin
                ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
                NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
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


    function iterate(ch, cs, updateH, updateS, niters, mix_iter; adjust_mix=false)
        println("iterate_fn $([updateH, updateS])   niters $niters mix_iter $mix_iter")
        err1 = -99.0
        err2 = -99.0
        for iter = 1:niters

            err_old = error_fn(NEWX, NEWY, ch, cs)

            if updateH
                if iter==1
                    iter_h = 1
                else
                    iter_h = 1
                end
                for i in 1:iter_h
                    @suppress begin
                        ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
                        println("size(VALS_FITTED ", size(VALS_FITTED))
                        NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
                        println("size NEWX ", size(NEWX))
                    end
                    
                    ch_new = NEWX \ NEWY
                    #println("ch_new size ", size(ch_new))
                    ch_test = ch_new*mix_iter + ch*(1 - mix_iter)

                    if mod(iter,10) == 1 && adjust_mix
                        for i = 1:20
                            println("current mix $mix_iter")
                            @suppress begin 
                                ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch_test, cs, solve_scf_mode)
                                NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
                            end
                            err_new = error_fn(NEWX, NEWY, ch_test, cs)
                            println("i $i err_old $err_old err_new $err_new")
                            if err_new > err_old
                                mix_iter = mix_iter * 0.6
                                ch_test = ch_new*mix_iter + ch*(1 - mix_iter)
                                if mix_iter < 0.01
                                    mix_iter = 0.01
                                    break
                                end
                            else
                                mix_iter = mix_iter*1.05
                                break
                            end
                        end
                    end
                    ch = ch_test
                    
                end
            end
                
            err1 = error_fn(NEWX, NEWY, ch, cs)

            if updateS
                
                if iter==1
                    iter_s = 2
                else
                    iter_s = 1
                end
                for i = 1:iter_s
                    @suppress begin
                        ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
                        NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
                    end
                    
                    
                    cst = NEWX_S \ ( (NEWX*ch - NEWY ))
                    #cs = (cs + cst)
                    cs_test = (cs + cst) * mix_iter + cs*(1 - mix_iter)

                    
                    if mod(iter,10) == 1 && adjust_mix 
                        for i = 1:20
                            println("current mix $mix_iter")
                            @suppress begin 
                                ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs_test, solve_scf_mode)
                                NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
                            end
                            err_new = error_fn(NEWX, NEWY, ch, cs_test)
                            println("i $i err_old $err_old err_new $err_new")
                            if err_new > err1*0.99
                                mix_iter = mix_iter * 0.6
                                cs_test = (cs + cst) * mix_iter + cs*(1 - mix_iter)
                                if mix_iter < 0.01
                                    mix_iter = 0.01
                                    break
                                end

                            else
                                mix_iter = mix_iter*1.05
                                break
                            end
                        end
                    end
                    cs = cs_test
                    
                end
                    
            end
#            println("size NEWX ", size(NEWX))
#            println("size NEWY ", size(NEWY))
#            println("size ch ", size(ch))
#            println("size cs ", size(cs))
            err2 = error_fn(NEWX, NEWY, ch, cs)
            println("iter $iter   $err1   $err2     ch1 $(ch[1]) cs1 $(cs[1]) ")
            if abs(err_old - err2) < 1e-5 && iter > 11
                break
            end
            err_old = err2
            
        end
        println("iterate_fn_end $err1 $err2")
        return ch, cs
    end

    if true
#        println("if true")
#        println("NCOLS, NCOLS_S ", [NCOLS, NCOLS_S])
        println("opt_S $opt_S")
        if opt_S
            ch, cs = iterate(ch, cs, true, false, 200, 0.05, adjust_mix=true)
            ch, cs = iterate(ch, cs, true, true, 550, 0.2, adjust_mix=false)
            ch, cs = iterate(ch, cs, true, true, 50, 0.02, adjust_mix=false)
        else
            ch, cs = iterate(ch, cs, true, false, 200, 0.05, adjust_mix=true)
            ch, cs = iterate(ch, cs, true, false, 150, 0.2, adjust_mix=false)
            ch, cs = iterate(ch, cs, true, false, 50, 0.02, adjust_mix=false)
        end
        ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
        NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)

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
        NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)

        good =  (abs.(ENERGIES - ENERGIES_working) ./ NAT) .< 0.05
        good = good .&& .!(Bool.(ERROR))
        
        database = make_database(ch, cs,  KEYS, HIND, SIND,DMIN_TYPES,DMIN_TYPES3, scf=scf, starting_database=starting_database, tbc_list = list_of_tbcs[good], fit_eam=fit_eam, fitting_version=fitting_version)
        return database

    end


    
#    return ret
    
    ####test stuff    

    if false
        ENERGIES_working, VECTS_FITTED, VALS_FITTED, OCCS_FITTED, VALS0_FITTED, ERROR, EDEN_FITTED = construct_fitted(ch, cs, solve_scf_mode)
        NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)

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
                @time  NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
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
                @time  NEWX, NEWY, NEWX_S, energy_counter, CALC_IND = construct_newXY(VECTS_FITTED, VALS_FITTED, OCCS_FITTED, NCALC, NCOLS, NCOLS_S, NLAM, ERROR, EDEN_FITTED, leave_out=leave_out)
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


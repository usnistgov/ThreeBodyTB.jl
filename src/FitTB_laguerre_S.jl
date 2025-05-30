


function do_fitting_recursive_S(list_of_tbcs ; weights_list = missing, dft_list=missing, kpoints = [0 0 0; 0 0 0.5; 0 0.5 0.5; 0.5 0.5 0.5; 0 0 0.25; 0 0.25 0; 0.25 0 0 ; 0.25 0.25 0.25; 0.25 0 0.25], starting_database = missing,  update_all = false, fit_threebody=true, fit_threebody_onsite=true, do_plot = false, energy_weight = missing, rs_weight=missing,ks_weight=missing, niters=50, lambda=0.0, leave_one_out=false, prepare_data = missing, RW_PARAM=0.0, NLIM = 100, refit_database = missing, start_small = false, fit_to_dft_eigs=false, fit_eam=false)

    if !ismissing(dft_list)
        println("top")
        KPOINTS, KWEIGHTS, nk_max = get_k(dft_list, length(dft_list), list_of_tbcs, NLIM=NLIM)
    else
        println("bot")
        KPOINTS, KWEIGHTS, nk_max = get_k_simple(kpoints, list_of_tbcs)
    end


    
#    println("KWEIGHTS 3 ", size(KWEIGHTS[3]), " " , KWEIGHTS[3][1:6])
    
    if ismissing(prepare_data)
        println("DO LINEAR FITTING")

        if update_all == true
            starting_database_t = missing #keep all
        else
            starting_database_t = starting_database
        end

        pd = do_fitting_linear(list_of_tbcs; kpoints = KPOINTS, mode=:kspace, dft_list = dft_list,  fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = false, starting_database=starting_database_t, return_database=false, NLIM=NLIM, refit_database=refit_database, fit_eam=fit_eam)
    else
        println("SKIP LINEAR MISSING")
        pd = prepare_data
#        database_linear, ch_lin, cs_lin, X_Hnew_BIG, Y_Hnew_BIG, X_H, X_Snew_BIG, Y_H, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3 = prepare_data
    end

    do_fitting_recursive_main(list_of_tbcs, pd; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = ks_weight, niters=niters, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam)
    
    database_linear, ch_lin, cs_lin, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, Ys_new, cs, ch_refit, SPIN  = pd

    ch_rec = do_fitting_recursive_main(list_of_tbcs, pd; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = ks_weight, niters=niters, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam, optimS=true)

    println("ch_rec ", ch_rec)
    println("reoptimize      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    pd_new = database_linear, ch_rec, cs_lin, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, Ys_new, cs, ch_refit, SPIN  

    database_final = do_fitting_recursive_main(list_of_tbcs, pd_new; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = ks_weight, niters=niters, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam, optimS=false)

    
    return database_final

end

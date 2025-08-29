
using ..CalcTB:n_2body_S_default

function do_fitting_recursive_S(list_of_tbcs ; weights_list = missing, dft_list=missing, kpoints = missing, starting_database = missing,  update_all = false, fit_threebody=true, fit_threebody_onsite=true, do_plot = false, energy_weight = missing, rs_weight=missing,ks_weight=missing, niters=50, lambda=0.0, leave_one_out=false, prepare_data = missing, RW_PARAM=0.0, NLIM = 100, refit_database = missing, start_small = false, fit_to_dft_eigs=false, fit_eam=false, opt=true)
    println("do_fitting_recursive_S")
#    sleep(3)
    
    if !ismissing(dft_list)
        println("top")
        KPOINTS, KWEIGHTS, nk_max = get_k(dft_list, length(dft_list), list_of_tbcs, NLIM=NLIM)
    else
        println("bot")
        KPOINTS, KWEIGHTS, nk_max = get_k_simple(kpoints, list_of_tbcs)
    end

    flush(stdout)
    
#    println("KWEIGHTS 3 ", size(KWEIGHTS[3]), " " , KWEIGHTS[3][1:6])
    pd = missing
    
    if ismissing(prepare_data)
        println("DO LINEAR FITTING")

        if update_all == true
            starting_database_t = missing #keep all
        else
            starting_database_t = starting_database
        end

        @suppress  pd = do_fitting_linear(list_of_tbcs; kpoints = KPOINTS, mode=:kspace, dft_list = dft_list,  fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, starting_database=starting_database_t, return_database=false, NLIM=NLIM, refit_database=refit_database, fit_eam=fit_eam)
    else
        println("SKIP LINEAR MISSING")
        pd = prepare_data
#        database_linear, ch_lin, cs_lin, X_Hnew_BIG, Y_Hnew_BIG, X_H, X_Snew_BIG, Y_H, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3 = prepare_data
    end

#    do_fitting_recursive_main(list_of_tbcs, pd; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = ks_weight, niters=niters, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam)

    flush(stdout)    
    database_linear, ch_lin, cs_lin, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, Ys_new, cs, ch_refit, SPIN,threebody_inds  = pd
    flush(stdout)
    top_vars = missing
    @suppress  top_vars = topstuff(list_of_tbcs,  pd; weights_list=weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight, ks_weight = ks_weight, niters=niters, lambda=lambda,  leave_one_out=false, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max = nk_max)

    ch_rec = missing
    error = missing
    println("DO RECURSIVE FITTING")

    ch_rec, error = do_fitting_recursive_main(list_of_tbcs, pd; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = ks_weight, niters=niters, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam, optimS=true, top_vars=top_vars)

    function update(cs)

        Ys_new = X_Snew_BIG * cs
        
        pd_new = database_linear, ch_rec, cs, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, Ys_new, cs, ch_refit, SPIN  , threebody_inds
        @suppress ch_rec, error = do_fitting_recursive_main(list_of_tbcs, pd_new; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = ks_weight, niters=50, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam, optimS=true, top_vars=top_vars)
        println("error  $error  ! cs $(round.(cs, digits=5))")
        flush(stdout)
        return error
    end


    
    cs_best = deepcopy(cs_lin)
    ch_best = deepcopy(ch_rec)
    error_best = 10.0^10;
    println("PRELIMINARY ")
    #for x = [1.05, 1.0, 0.95]
    for x in 1.1:-0.2:0.5
        for y in 1.2:-0.2:0.8
            println("x $x y $y")
            if length(cs_lin) == n_2body_S_default
                error = update(cs_lin[1:n_2body_S_default] * x*y)
            else
                error = update( vcat(cs_lin[1:n_2body_S_default]*x*y, cs_lin[ (1+n_2body_S_default) : end ] * y))
            end
            if error < error_best
                error_best = error
                cs_best = vcat(cs_lin[1:n_2body_S_default]*x*y, cs_lin[ (1+n_2body_S_default) : end ] * y)
                ch_best = deepcopy(ch_rec)
            end
        end
    end
    cs_final = deepcopy(cs_best)
    println("cs final prelim ", cs_final)
    Ys_new = X_Snew_BIG * cs_final
    println("reoptimize final      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    pd_new = database_linear, ch_rec, cs_final, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, Ys_new, cs_final, ch_refit, SPIN  , threebody_inds
    database_final = missing
    println()
    @suppress database_final = do_fitting_recursive_main(list_of_tbcs, pd_new; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = 0.0, niters=niters, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam, optimS=false, top_vars=top_vars)

    if opt
        println("Optim min ")
        ch_rec = deepcopy(ch_best) #start
        cs_lin = deepcopy(cs_best)

        inds = []
        for x = 1 : Int64(length(cs_lin) / n_2body_S_default)
            for y = 1:n_2body_S_default
                push!(inds, y + (x-1) * n_2body_S_default)
            end

        end
        println("inds $inds")

        error_min = 10.0^20
        cs_min = deepcopy(cs_lin)
        function update2(cs)
            cst = deepcopy(cs_lin)
            cst[inds] = cs
            Ys_new = X_Snew_BIG * cst
            
            pd_new = database_linear, ch_rec, cst, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, Ys_new, cst, ch_refit, SPIN  ,threebody_inds
            @suppress ch_rec, error = do_fitting_recursive_main(list_of_tbcs, pd_new; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = ks_weight, niters=50, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam, optimS=true, top_vars=top_vars)
            println("error  $error  ! cst $(round.(cst, digits=5))")
            flush(stdout)
            if error < error_min
                cs_min = deepcopy(cst)
                error_min = error
            end
            return error
        end

        
        opts = Optim.Options(f_abstol = 1e-2, x_abstol = 1e-2, f_calls_limit = 800,
                             iterations = 800,                                                                                                                                                                
                             store_trace = true,                                                                                                                                                                 
                             show_trace = false)                                                                                                                                                                 

        if maximum(abs.(cs_lin)) > 100.0
            vv = 0.002
        elseif maximum(abs.(cs_lin)) > 10.0
            vv = 0.005
        else
            vv = 0.02
        end
        
        res = optimize(update2, cs_lin[inds],NelderMead(initial_simplex = Optim.AffineSimplexer(vv, vv)),  opts)
        cs_final = deepcopy(cs_final)
        cs_final[inds] = Optim.minimizer(res)
        println("cs_final $cs_final")

        cs_final = deepcopy(cs_min)
        println("cs_final final $cs_final")
        
        Ys_new = X_Snew_BIG * cs_final
        println("reoptimize final      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        pd_new = database_linear, ch_rec, cs_final, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, Ys_new, cs_final, ch_refit, SPIN  , threebody_inds
        database_final_opt = missing
        @suppress database_final_opt = do_fitting_recursive_main(list_of_tbcs, pd_new; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = 0.0, niters=niters, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam, optimS=false, top_vars=top_vars)

    else
        database_final_opt = Dict()
    end
        
    
    return database_final, database_final_opt

end



#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


function do_fitting_recursive_S2(list_of_tbcs ; weights_list = missing, dft_list=missing, kpoints = [0 0 0; 0 0 0.5; 0 0.5 0.5; 0.5 0.5 0.5; 0 0 0.25; 0 0.25 0; 0.25 0 0 ; 0.25 0.25 0.25; 0.25 0 0.25], starting_database = missing,  update_all = false, fit_threebody=true, fit_threebody_onsite=true, do_plot = false, energy_weight = missing, rs_weight=missing,ks_weight=missing, niters=50, lambda=0.0, leave_one_out=false, prepare_data = missing, RW_PARAM=0.0, NLIM = 100, refit_database = missing, start_small = false, fit_to_dft_eigs=false, fit_eam=false, opt=true)
    println("do_fitting_recursive_S")
    sleep(3)
    
    if !ismissing(dft_list)
        println("top")
        KPOINTS, KWEIGHTS, nk_max = get_k(dft_list, length(dft_list), list_of_tbcs, NLIM=NLIM)
    else
        println("bot")
        KPOINTS, KWEIGHTS, nk_max = get_k_simple(kpoints, list_of_tbcs)
    end

    flush(stdout)
    
#    println("KWEIGHTS 3 ", size(KWEIGHTS[3]), " " , KWEIGHTS[3][1:6])
    pd = missing
    
    if ismissing(prepare_data)
        println("DO LINEAR FITTING")

        if update_all == true
            starting_database_t = missing #keep all
        else
            starting_database_t = starting_database
        end

        @suppress pd = do_fitting_linear(list_of_tbcs; kpoints = KPOINTS, mode=:kspace, dft_list = dft_list,  fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = false, starting_database=starting_database_t, return_database=false, NLIM=NLIM, refit_database=refit_database, fit_eam=fit_eam)
    else
        println("SKIP LINEAR MISSING")
        pd = prepare_data
#        database_linear, ch_lin, cs_lin, X_Hnew_BIG, Y_Hnew_BIG, X_H, X_Snew_BIG, Y_H, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3 = prepare_data
    end

#    do_fitting_recursive_main(list_of_tbcs, pd; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = ks_weight, niters=niters, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam)

    flush(stdout)    
    database_linear, ch_lin, cs_lin, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, Ys_new, cs, ch_refit, SPIN, threebody_inds  = pd
    flush(stdout)
    top_vars = missing
    @suppress  top_vars = topstuff(list_of_tbcs,  pd; weights_list=weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight, ks_weight = ks_weight, niters=niters, lambda=lambda,  leave_one_out=false, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max = nk_max)

    ch_rec = missing
    error = missing
    println("DO RECURSIVE FITTING")

    @suppress ch_rec, error = do_fitting_recursive_main(list_of_tbcs, pd; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = ks_weight, niters=niters, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam, optimS=true, top_vars=top_vars)

    function update(cs)

        Ys_new = X_Snew_BIG * cs
        
        pd_new = database_linear, ch_rec, cs, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, Ys_new, cs, ch_refit, SPIN  , threebody_inds
        @suppress ch_rec, error = do_fitting_recursive_main(list_of_tbcs, pd_new; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = ks_weight, niters=12, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam, optimS=true, top_vars=top_vars)
        println("error  $error  ! cs $(round.(cs, digits=5))")
        flush(stdout)
        return error
    end


    
    cs_best = deepcopy(cs_lin)
    ch_best = deepcopy(ch_rec)
    error_best = 10.0^10;
    println("PRELIMINARY ")
    #for x = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for x in 1.05:-0.05:0.0
        println("x $x ")
        error = update(cs_lin * x)
        if error < error_best
            error_best = error
            cs_best = cs_lin * x
            ch_best = deepcopy(ch_rec)
        end
    end
    cs_final = deepcopy(cs_best)

    Ys_new = X_Snew_BIG * cs_final
    println("reoptimize final      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    pd_new = database_linear, ch_rec, cs_final, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, Ys_new, cs_final, ch_refit, SPIN  , threebody_inds
    database_final = missing
    @suppress database_final = do_fitting_recursive_main(list_of_tbcs, pd_new; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = 0.0, niters=niters, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam, optimS=false, top_vars=top_vars)

    if opt
        println("Optim min ")
        ch_rec = deepcopy(ch_best) #start
        cs_lin = deepcopy(cs_best)

        ncopy = Int64(length(cs_lin) / n_2body_S_default)
        
        function update2(X)
            x = X[2]

            t = [1, x^1, x^2, x^3, x^4, x^5, x^6, x^7][1:n_2body_S_default]
            tt = repeat(t, ncopy)
            cst = deepcopy(cs_lin) .* X[1] .* tt
        #    cst[inds] = cs
            Ys_new = X_Snew_BIG * cst
            
            pd_new = database_linear, ch_rec, cst, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, Ys_new, cst, ch_refit, SPIN  , threebody_inds
            @suppress ch_rec, error = do_fitting_recursive_main(list_of_tbcs, pd_new; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = ks_weight, niters=12, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam, optimS=true, top_vars=top_vars)
            println("error  $error !  $X   ! cst $(round.(cst, digits=5))")
            flush(stdout)
            return error
        end

        
        opts = Optim.Options(f_abstol = 1e-2, x_abstol = 1e-2, f_calls_limit = 500,
                             iterations = 500,                                                                                                                                                                
                             store_trace = true,                                                                                                                                                                 
                             show_trace = false)                                                                                                                                                                 

        vv= 0.002
        
        res = optimize(update2, [1.0, 1.0],NelderMead(initial_simplex = Optim.AffineSimplexer(vv, vv)),  opts)
        X = Optim.minimizer(res)
        println("X $X")
        x = X[2]
        t = [1, x^1, x^2, x^3, x^4, x^5, x^6, x^7][1:n_2body_S_default ]
        tt = repeat(t, ncopy)
        cs_final = deepcopy(cs_lin) .* X[1] .* tt

        Ys_new = X_Snew_BIG * cs_final
        println("reoptimize final      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        pd_new = database_linear, ch_rec, cs_final, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, h_on, ind_BIG, KEYS, HIND, SIND, DMIN_TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, Ys_new, cs_final, ch_refit, SPIN  , threebody_inds
        database_final_opt = missing
        @suppress database_final_opt = do_fitting_recursive_main(list_of_tbcs, pd_new; weights_list = weights_list, dft_list=dft_list, kpoints = kpoints, starting_database = starting_database,  update_all = update_all, fit_threebody=fit_threebody, fit_threebody_onsite=fit_threebody_onsite, do_plot = do_plot, energy_weight = energy_weight, rs_weight=rs_weight,ks_weight = 0.0, niters=niters, lambda=lambda, leave_one_out=leave_one_out, RW_PARAM=RW_PARAM, KPOINTS=KPOINTS, KWEIGHTS=KWEIGHTS, nk_max=nk_max,  start_small = start_small , fit_to_dft_eigs=fit_to_dft_eigs, fit_eam=fit_eam, optimS=false, top_vars=top_vars)

    else
        database_final_opt = Dict()
    end
        
    
    return database_final, database_final_opt

end

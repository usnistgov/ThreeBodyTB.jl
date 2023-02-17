using ..Symmetry:get_kgrid_sym

##############################################################################################################
"""
    function get_energy_force_stress_fft(tbc::tb_crys, database; do_scf=false, smearing = 0.01, grid = missing, e_den0=missing, vv = missing)

Calculate energy/force/stress using fft algorithm. Users should use `scf_energy_force_stress`, which calls this. Uses automatic differentation for jacobian.
"""
function get_energy_force_stress_fft_LV_sym(tbc::tb_crys, database; do_scf=false, smearing = 0.01, grid = missing, e_den0=missing, vv = missing, nspin = 1, repel=true)

#    println("get_energy_force_stress_fft")

    do_scf = true
    
    #FloatX = Float32
    FloatX = Float64
    ct = deepcopy(tbc.crys)


    
    #println("dist")
    begin
        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel_LV(ct,cutoff2X,cutoff3bX,var_type=FloatX, return_floats=false)
        DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3
    end
    
    if ismissing(grid)
        grid = get_grid(tbc.crys)
    end
    
    kgrid, kweights = make_kgrid(grid)
    nk = size(kgrid)[1]

    nk_red, grid_ind, kpts, kweights = get_kgrid_sym(tbc.crys, grid=grid)
    println("nk $nk nk_red $nk_red")
    
#    println("get_energy_force_stress_fft2")
    
    #println("test safe get_energy_force_stress_fft")
    tooshort, energy_tot = safe_mode_energy(tbc.crys, database, DIST=DIST)
    
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
        
            
        
        
#        println("not too ")
        
        scf = database["scf"]
#        println("scf ", scf)
        if !(tooshort)
            if !ismissing(vv)
                VECTS, VALS, efermi = vv
                energy_tot = 0.0
                if !ismissing(e_den0)
                    e_den = e_den0
                else
                    e_den = tbc.eden
                end
            else
                #prepare eigenvectors / values
                error_flag = false
                if do_scf
#                    println("do_scf")
                    energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=e_den0, conv_thr = 1e-8, nspin=nspin, verbose=false, use_sym = true)
                else
#                    println("calc_energy_charge_fft")
                    energy_tot, efermi, e_den, VECTS, VALS, error_flag =  calc_energy_charge_fft(tbc, grid=grid, smearing=smearing)
                end
                if error_flag
                    println("warning, trouble with eigenvectors/vals in initial step get_energy_force_stress")
                end
            end
            h1, dq = get_h1(tbc, e_den)
            if nspin == 2
                h1spin = get_spin_h1(tbc, e_den)
            else
                h1spin = zeros(2,tbc.tb.nwan, tbc.tb.nwan)
            end
            OCCS = gaussian.(VALS.-efermi, smearing)
        end
#        println("sum OCCS get_energy_force_stress_fft ", sum(OCCS))
        
    #println("shrink")
        begin
            ht = (abs.(tbc.tb.H[1,:,:,:]) .> 1e-7) .|| (abs.(tbc.tb.S) .> 1e-7)
            num_nonzero = sum(ht)

            nz_arr = zeros(Int64, num_nonzero)

            counter = 0
            for i = 1:tbc.tb.nwan * tbc.tb.nwan * tbc.tb.nr
                if ht[i] == 1
                    counter += 1
                    nz_arr[counter] = i
                end
            end
            size_ret = num_nonzero
#            println("counter $counter size_ret $size_ret normal ", prod(size(tbc.tb.H)))
        end        

        #size_ret = tbc.tb.nwan * tbc.tb.nwan * tbc.tb.nr 

        if database["scf"] == true
            scf = true
        else
            scf = false
        end
        
        dontcheck = calc_tb_LV(ct, database; check_only=true, use_threebody=true, use_threebody_onsite=true, DIST=DIST, verbose=false)
#        println("dontcheck $dontcheck")
        

        chunksize=min(15, 3*ct.nat + 6)        

        g_ew = zeros(FloatX, 1, 3*ct.nat + 6)
        size_ret_full = tbc.tb.nwan * tbc.tb.nwan * tbc.tb.nr 
        g = zeros(FloatX, size_ret_full*2 + 1, 3*ct.nat + 6)

        println("Calculate Jacobian of TB object")
        @time begin

            #println("ew")
            EW = begin
                FN_ew = x->ew(x,ct,FloatX, dq)
                if scf
                    cfg = ForwardDiff.JacobianConfig(FN_ew, zeros(FloatX, 3*ct.nat + 6), ForwardDiff.Chunk{chunksize}())
                    g_ew = ForwardDiff.jacobian(FN_ew, zeros(FloatX, 3*ct.nat + 6) , cfg ) ::  Array{FloatX,2}
                end
            end
#            println("done ew")
#            println("ham")
            HAM = begin

                FN_ham = x->ham(x,ct,database,dontcheck, repel, DIST, nz_arr, FloatX, size_ret)
                
                chunksize=min(15, 3*ct.nat + 6)
                cfg = ForwardDiff.JacobianConfig(FN_ham, zeros(FloatX, 3*ct.nat + 6), ForwardDiff.Chunk{chunksize}())
                #println("jacham")
                g_nz = ForwardDiff.jacobian(FN_ham, zeros(FloatX, 3*ct.nat + 6) , cfg ) ::  Array{FloatX,2}


                counter = 0
                for i = 1:size_ret
                    g[nz_arr[i],:] = g_nz[i,:]
                    g[size_ret_full + nz_arr[i],:] = g_nz[size_ret + i,:]
                end
                g[end,:] = g_nz[end,:]
                g_nz = missing
            end
#            println("done ham")
            #println("wait")
#            wait(HAM)
#            wait(EW)
        end

        #return g
        
        if scf
            g[end,:] = g_ew[1,:]
        end
        size_ret = tbc.tb.nwan * tbc.tb.nwan * tbc.tb.nr 
        println("Calculate Force / Stress")
        @time begin
        
            hr_g = zeros(Complex{eltype(g)},  tbc.tb.nwan, tbc.tb.nwan,  grid[1], grid[2], grid[3] )
            sr_g = zeros(Complex{eltype(g)},  tbc.tb.nwan, tbc.tb.nwan,  grid[1], grid[2], grid[3] )
            hk_g = similar(hr_g)
            sk_g = similar(sr_g)
            
            
            garr =zeros(Float64, 3*ct.nat+6)

            #        pVECTS = permutedims(Complex{FloatX}.(VECTS), [2,3,1])
            pVECTS = permutedims(Complex{FloatX}.(VECTS), [3,4,1,2])
            pVECTS_conj = conj.(pVECTS)


            #        hk_g_temp = deepcopy(hk_g)
            
            #VALS0 = zeros(FloatX,  prod(grid), tbc.tb.nwan, nspin)
            VALS0 = zeros(FloatX,  nk_red, nspin)
            
            OCCS = Float32.(OCCS)

            #println("denmat")
            begin
                
                DENMAT = zeros(Complex{FloatX}, tbc.tb.nwan, tbc.tb.nwan, nk_red, nspin)
                DENMAT_V = zeros(Complex{FloatX}, tbc.tb.nwan, tbc.tb.nwan, nk_red, nspin)

                go_denmat_sym!(DENMAT, DENMAT_V, grid, nspin, OCCS, VALS, pVECTS, pVECTS_conj, tbc, nk_red, grid_ind, kweights )
                println("size DENMAT $(size(DENMAT)) VALS $(size(VALS)) VECTS $(size(VECTS)) OCCS $(size(OCCS))")
            end            

            #println("FIND")
            for FIND in 1:3*ct.nat + 6
                
                VALS0 .= 0.0
                
                hr_g .= 0.0
                sr_g .= 0.0


                forloops!(tbc, hr_g, sr_g, size_ret, FIND, grid, g)
                
                
                begin 

                    h = @spawn begin
                        hk_g .= fft(hr_g, [3,4,5])
                    end
                    s = @spawn begin 
                        sk_g .= fft(sr_g, [3,4,5])
                    end

                    wait(h)
                    wait(s)
                end

                psi_gradH_psi3_sym(VALS0, hk_g, sk_g, Complex{FloatX}.(h1), FloatX.(h1spin), scf, tbc.tb.nwan, ct.nat, grid, DENMAT, DENMAT_V, nk_red, grid_ind, kweights)            

                
                garr[FIND] = sum(VALS0) 

            end #end FIND
            
            
            if scf
                garr += g[end,:][:]
            end
            
        end

    end
    println()
        
    x, stress = reshape_vec(garr, ct.nat)
    f_cart = -1.0 * x
    f_cart = f_cart * inv(ct.A)' / nspin
    stress = -stress / abs(det(ct.A)) /nspin

    println(f_cart)
    
    #neaten
    for i = 1:3
        for j = 1:3
            if abs(stress[i,j]) < 1e-9
                stress[i,j] = 0.0
            end
        end
    end
    for i = 1:ct.nat
        for j = 1:3
            if abs(f_cart[i,j]) < 1e-7
                f_cart[i,j] = 0.0
            end
        end
    end


    
    return energy_tot,  f_cart, stress


end

function go_denmat_sym!(DENMAT, DENMAT_V, grid, nspin, OCCS, VALS, pVECTS, pVECTS_conj, tbc, nk_red, grid_ind, kweights)
    for spin = 1:nspin
        for a =1:tbc.tb.nwan
        @inbounds @fastmath for c1 = 1:tbc.tb.nwan
            for c2 = 1:tbc.tb.nwan
                    for counter = 1:nk_red
                        
                        DENMAT[c1,c2,counter, spin]   += kweights[counter]*OCCS[counter,a, spin] * ( pVECTS[c1,a,counter, spin]) * ( pVECTS_conj[c2,a,counter, spin])                 
                        DENMAT_V[c1,c2,counter, spin]   += kweights[counter]*OCCS[counter,a, spin] * ( pVECTS[c1,a,counter, spin]) * ( pVECTS_conj[c2,a,counter, spin]) * VALS[counter, a, spin]
                
                    end
                end
            end
        end
    end
    DENMAT=DENMAT/2.0
    DENMAT_V=DENMAT_V/2.0
end


#=function ew(x :: Vector, ct, FloatX, dq)
    #            println("ew top")
    T=typeof(x[1])
    ret = zeros(T, 1)


    x_r, x_r_strain = reshape_vec(x, ct.nat, strain_mode=true)
    A = FloatX.(ct.A) * (I(3) + x_r_strain)
    crys_dual = makecrys( A , ct.coords + x_r, ct.types, units="Bohr")
    
    kappa = estimate_best_kappa(FloatX.(ct.A))
    gamma_dual, background_charge_correction = electrostatics_getgamma(crys_dual, kappa=kappa)
    eewald, pot = ewald_energy(crys_dual, gamma_dual,background_charge_correction, dq)
    ret[end] = eewald
    return ret
end
=#

#=function forloops!(tbc, hr_g, sr_g, size_ret, FIND, grid, g)

            
    for c in 1:size(tbc.tb.ind_arr)[1] #@threads 
        
        ind = tbc.tb.ind_arr[c,:]
        new_ind = [mod(ind[1], grid[1])+1, mod(ind[2], grid[2])+1, mod(ind[3], grid[3])+1]
        @inbounds @fastmath @threads for nb = 1:tbc.tb.nwan

            #hr_g[:,nb,new_ind[1], new_ind[2], new_ind[3]] += @view g[(1:tbc.tb.nwan) .+ ((nb-1) * tbc.tb.nwan + tbc.tb.nwan^2 * (c-1)) , FIND]
            #sr_g[:,nb,new_ind[1], new_ind[2], new_ind[3]] += @view g[(1:tbc.tb.nwan) .+  (size_ret + (nb-1) * tbc.tb.nwan + tbc.tb.nwan^2 * (c-1)), FIND ]
            
            for na = 1:tbc.tb.nwan
                hr_g[na,nb,new_ind[1], new_ind[2], new_ind[3]] += g[na + (nb-1) * tbc.tb.nwan + tbc.tb.nwan^2 * (c-1) , FIND]
                sr_g[na,nb,new_ind[1], new_ind[2], new_ind[3]] += g[size_ret + na + (nb-1) * tbc.tb.nwan + tbc.tb.nwan^2 * (c-1), FIND ]
            end
        end
    end

    
end
=#



function psi_gradH_psi3_sym(VALS0, hk_g, sk_g, h1, h1spin, scf, nwan, nat, grid, DENMAT, DENMAT_V, nk_red, grid_ind, kweights)

    nspin = size(VALS0)[end]
    nwan1, nwan2, G1, G2, G3 = size(hk_g)
    #println("typeof(VALS0) ", typeof(VALS0), " size ", size(VALS0))
    VALS0_t = zeros(nk_red)
    #println("typeof(VALS0_t) ", typeof(VALS0_t), " size ", size(VALS0_t))
    for spin = 1:nspin

        hk_g_temp = deepcopy(hk_g)
       # hk_g_temp = zeros(eltype(hk_g), nwan, nwan, nk_red)

        if scf
            @inbounds @fastmath for a1 = 1:nwan
                for a2 = 1:nwan
                    for k = 1:nk_red
                        g1,g2,g3 = grid_ind[k,:]
                        hk_g_temp[a2,a1,g1,g2,g3] +=  (h1[a2,a1] + h1spin[spin,a2,a1]) .* ( sk_g[a2,a1,g1,g2,g3])
                    end
                end
            end
        end
        
        @inbounds @fastmath for a1 = 1:nwan
            for a2 = 1:nwan
                #c=0
                for c = 1:nk_red
                    k1,k2,k3=grid_ind[c,:]
                    #                @threads for k1 = 1:grid[1]
#                    for k2 = 1:grid[2]
#                        for k3 = 1:grid[3]
#                            c = (k1-1)*grid[2]*grid[3] + (k2-1)*grid[3] + k3
                    VALS0_t[c] += real( (DENMAT[a1,a2,c, spin]) * (hk_g_temp[a2,a1,k1,k2,k3]) - (DENMAT_V[a1,a2,c, spin]) * (sk_g[a2,a1,k1,k2,k3]))
                end
            end
        end
        VALS0[:,spin] = VALS0_t[:]
    end
end
using ..TB:tb_crys_sparse
using ..Symmetry:get_kgrid_sym
using ..Symmetry:symmetrize_vector_tensor
using ..CalcTB:calc_tb_LV_sparse

using SparseArrays

function ham_SINGLE_sparse(x :: Vector, ct, database, dontcheck, repel, DIST, FloatX, nwan, nr, atom)
#    println("HS atom $atom")
    begin
#        println("begin")
        begin
            T=typeof(x[1])
            println("atom $atom x $x ")
            println(ct)
            
            if atom >= 1
                println([typeof(x), typeof(atom), typeof(ct.nat)])
                x_r, x_r_strain = reshape_vec_SINGLE(x, atom, ct.nat, strain_mode=false)
                A = FloatX.(ct.A) * (I(3))
            else
                x_r, x_r_strain = reshape_vec_SINGLE(x, -1, ct.nat, strain_mode=true)
                A = FloatX.(ct.A) * (I(3) + x_r_strain)            
            end
        end        
        #println("dual")
        crys_dual = makecrys( A , ct.coords + x_r, ct.types, units="Bohr", type=eltype(x))
        #println("calc_tb_LV")
        H, S, INDarr = calc_tb_LV_sparse(crys_dual, database; verbose=true, var_type=T, use_threebody=true, use_threebody_onsite=true, gamma=zeros(T, ct.nat,ct.nat) , check_frontier= !dontcheck, repel=repel, DIST=DIST, retmat=true, atom = atom)
        
        nnz = 0
        for n = 1:length(H)
            nnz += length(H[n])
        end

        ret = zeros(T, nnz*2)
        for n = 1:length(H)
            ret[1:length(H[n])] = H[n]
            ret[ nnz .+ (1:length(H[n]))] = S[n]
        end

    end

    return ret
end



function get_energy_force_stress_fft_LV_sym_SINGLE(tbc::tb_crys_sparse, database; do_scf=false, smearing = 0.01, grid = missing, e_den0=missing, vv = missing, nspin = 1, repel=true)

#    println("SINGLE tot", tbc.tot_charge)
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

    #println("grid sym")
    nk_red, grid_ind, kpts, kweights = get_kgrid_sym(tbc.crys, grid=grid)
#    println("nk $nk nk_red $nk_red")
    
#    println("get_energy_force_stress_fft2")
    
    #println("test safe get_energy_force_stress_fft")
    tooshort, energy_tot = safe_mode_energy(tbc.crys, database, DIST=DIST)
    
    #println("main if")
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
        println("SINGLE scf ", scf)
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
                    energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=e_den0, conv_thr = 1e-6, nspin=nspin, verbose=false, use_sym = true)
                else
                    println("calc_energy_charge_fft")
                    energy_tot, efermi, e_den, VECTS, VALS, error_flag =  calc_energy_charge_fft(tbc, grid=grid, smearing=smearing, repel=repel)
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

        #size_ret = tbc.tb.nwan * tbc.tb.nwan * tbc.tb.nr 

        if database["scf"] == true
            scf = true
        else
            scf = false
        end
        
        dontcheck, sum_repel = calc_tb_LV(ct, database; check_only=true, use_threebody=true, use_threebody_onsite=true, DIST=DIST, verbose=false)
#        println("dc $dontcheck sum_repel $sum_repel")
        dontcheck = dontcheck && sum_repel

#        println("dontcheck $dontcheck")
        

        chunksize=min(15, 3*ct.nat + 6)        

        g_ew = zeros(FloatX, 1, 3*ct.nat + 6)
        size_ret_full = tbc.tb.nwan * tbc.tb.nwan * tbc.tb.nr 
#        g = zeros(FloatX, size_ret_full*2 + 1, 3*ct.nat + 6)

        println("(Pre-)Calculate Jacobian of TB object")
        begin

            #println("ew")
            EW = begin
                FN_ew = x->ew(x,ct,FloatX, dq)
                if scf
                    cfg = ForwardDiff.JacobianConfig(FN_ew, zeros(FloatX, 3*ct.nat + 6), ForwardDiff.Chunk{chunksize}())
                    g_ew = ForwardDiff.jacobian(FN_ew, zeros(FloatX, 3*ct.nat + 6) , cfg ) ::  Array{FloatX,2}
                end
            end
            
            #println("done ew")
            #println("ham")
            HAM = begin

                
#                Hdual = zeros(ForwardDiff.Dual{FloatX}, tbc.tb.nwan, tbc.tb.nwan, tbc.tb.nr)
#                Sdual = zeros(ForwardDiff.Dual{FloatX}, tbc.tb.nwan, tbc.tb.nwan, tbc.tb.nr)

                memory = Dict()
                ATOM = 1
                FN_ham = x->ham_SINGLE_sparse(x,ct,database,dontcheck, repel, DIST, FloatX, tbc.tb.nwan, tbc.tb.nr, ATOM)
                
                #chunksize=min(15, 3*ct.nat + 6)
                chunksize=3
                cfg3 = ForwardDiff.JacobianConfig(FN_ham, zeros(FloatX, 3), ForwardDiff.Chunk{chunksize}())
                chunksize=6
                cfg6 = ForwardDiff.JacobianConfig(FN_ham, zeros(FloatX, 6), ForwardDiff.Chunk{chunksize}())
                #println("jacham")


                #g_nz = ForwardDiff.jacobian(FN_ham, zeros(FloatX, 3*ct.nat + 6) , cfg ) ::  Array{FloatX,2}


                #counter = 0
                #for i = 1:size_ret
                #    g[nz_arr[i],:] = g_nz[i,:]
                #    g[size_ret_full + nz_arr[i],:] = g_nz[size_ret + i,:]
                #end
                #g[end,:] = g_nz[end,:]
                #g_nz = missing
            end
#            println("done ham")
            #println("wait")
#            wait(HAM)
#            wait(EW)
        end

        ATOM = 1
        println("repel $repel")
        jjj = FN_ham(zeros(FloatX, 3))
        jac = ForwardDiff.jacobian(FN_ham, zeros(FloatX, 3) , cfg3 ) ::  Array{FloatX,2}

        return jjj, jac
        
        #        return g
        
#        if scf
#            g[end,:] = g_ew[1,:]
#        end
        size_ret = tbc.tb.nwan * tbc.tb.nwan * tbc.tb.nr 
        println("Calculate Force / Stress")
        begin
        
            #println("mem")
            begin
#                hr_g = zeros(Complex{eltype(g)},  tbc.tb.nwan, tbc.tb.nwan,  grid[1], grid[2], grid[3] )
#                sr_g = zeros(Complex{eltype(g)},  tbc.tb.nwan, tbc.tb.nwan,  grid[1], grid[2], grid[3] )

                hr_g_r = zeros(FloatX,    grid[1], grid[2], grid[3], tbc.tb.nwan, tbc.tb.nwan )
                sr_g_r = zeros(FloatX,    grid[1], grid[2], grid[3], tbc.tb.nwan, tbc.tb.nwan )
                
                hr_g = zeros(Complex{FloatX},    grid[1], grid[2], grid[3], tbc.tb.nwan, tbc.tb.nwan )
                sr_g = zeros(Complex{FloatX},    grid[1], grid[2], grid[3], tbc.tb.nwan, tbc.tb.nwan )

#                hk_g = similar(hr_g)
#                sk_g = similar(sr_g)
                

                sk_g_r = zeros(FloatX, grid[1], grid[2], grid[3], tbc.tb.nwan, tbc.tb.nwan)
                sk_g_i = zeros(FloatX, grid[1], grid[2], grid[3], tbc.tb.nwan, tbc.tb.nwan)
                hk_g_r = zeros(FloatX, grid[1], grid[2], grid[3], tbc.tb.nwan, tbc.tb.nwan)
                hk_g_i = zeros(FloatX, grid[1], grid[2], grid[3], tbc.tb.nwan, tbc.tb.nwan)

                
                garr =zeros(Float64, 3*ct.nat+6)

                #        pVECTS = permutedims(Complex{FloatX}.(VECTS), [2,3,1])
                pVECTS = permutedims(Complex{FloatX}.(VECTS), [3,4,1,2])
                pVECTS_conj = conj.(pVECTS)


                #        hk_g_temp = deepcopy(hk_g)
                
                #VALS0 = zeros(FloatX,  prod(grid), tbc.tb.nwan, nspin)
                VALS0 = zeros(FloatX,  nk_red, nspin)
                
                OCCS = Float32.(OCCS)
            end
            #println("denmat")
            begin
                
                DENMAT = zeros(Complex{FloatX}, tbc.tb.nwan, tbc.tb.nwan, nk_red, nspin)
                DENMAT_V = zeros(Complex{FloatX}, tbc.tb.nwan, tbc.tb.nwan, nk_red, nspin)
#                println("SUM OCCS ", sum(abs.(OCCS)))
                go_denmat_sym!(DENMAT, DENMAT_V, grid, nspin, OCCS, VALS, pVECTS, pVECTS_conj, tbc, nk_red, grid_ind, kweights )
#                println("size DENMAT $(size(DENMAT)) VALS $(size(VALS)) VECTS $(size(VECTS)) OCCS $(size(OCCS))")

                DENMAT_r = real(DENMAT)
                DENMAT_i = imag(DENMAT)
                DENMAT_V_r = real(DENMAT_V)
                DENMAT_V_i = imag(DENMAT_V)



            end            

            #println("declare mem")
            begin
                gatom = zeros(FloatX,2*size_ret+1 ,3)
                gstress = zeros(FloatX,2*size_ret+1 ,6)
            end
#            print("atom ")
            for atom in 0:ct.nat
#                print("$atom ")
                ATOM=atom
                if atom >= 1
                    ForwardDiff.jacobian!(gatom, FN_ham, zeros(FloatX, 3) , cfg3 ) ::  Array{FloatX,2}
                    g = gatom
                    find= 1:3
#                    println("XXXXXXXXXXXX atom $atom ", sum(abs.(gatom)))
                else
                    ForwardDiff.jacobian!(gstress, FN_ham, zeros(FloatX, 6) , cfg6 ) ::  Array{FloatX,2}
                    g = gstress
                    find = 1:6
                end                    
                #println("sum abs gz ", sum(abs.(g)))
#                println("done atom; find")
                for FIND in find

                    #println("FIND $FIND")
                    VALS0 .= 0.0
                    
                    hr_g .= 0.0
                    sr_g .= 0.0
#                    hr_g_r .= 0.0
#                    sr_g_r .= 0.0


                #    println("forl")
                    #                @time forloops!(tbc, hr_g, sr_g, size_ret, FIND, grid, g)
                    begin
                    #    println("normal")
                        #@time forloops2_SINGLE!(tbc, hr_g_r, sr_g_r, size_ret, FIND, grid, g)
                    #    println("check")
                        forloops2_SINGLE_check!(tbc, hr_g, sr_g, size_ret, FIND, grid, g)                        
#                        hr_g[:] .= hr_g_r[:]
#                        sr_g[:] .= sr_g_r[:]

                    end
                    
                    #println("fft ", size(hr_g))
                    #println("fft")

                    begin 
                        
                        h = @spawn begin
                            #   hk_g .= fft(hr_g, [1,2,3])
                            fft!(hr_g, [1,2,3])
                        end
                        s = @spawn begin 
                            #fft!(sr_g, [3,4,5])
                            fft!(sr_g, [1,2,3])
                        end

                        wait(h)
                        wait(s)
                    end

#                    println("psi")
#                    if atom >= 1
#                        println("ATOM $atom $FIND ", sum(abs.(hr_g)), " ", sum(abs.(sr_g)), " " , sum(abs.(DENMAT_r)), " " , sum(abs.(DENMAT_V_r)))
#                    end
                    #@time psi_gradH_psi3_sym2(VALS0, hr_g, sr_g, FloatX.(h1), FloatX.(h1spin), scf, tbc.tb.nwan, ct.nat, grid, DENMAT_r, DENMAT_i, DENMAT_V_r, DENMAT_V_i, nk_red, grid_ind, kweights, sk_g_r, sk_g_i, hk_g_r, hk_g_i)

                    psi_gradH_psi3_sym2_sparse(VALS0, hr_g, sr_g, FloatX.(h1), FloatX.(h1spin), scf, tbc.tb.nwan, ct.nat, grid, DENMAT_r, DENMAT_i, DENMAT_V_r, DENMAT_V_i, nk_red, grid_ind, kweights, sk_g_r, sk_g_i, hk_g_r, hk_g_i, DENMAT, DENMAT_V)                                

                    #println("atom $atom sum vals ", sum(VALS0))
                    #println("sum")
                    if atom >= 1
                        garr[3*(atom-1) + FIND] += sum(VALS0)
#                        println("ATOM END $atom $FIND ", sum(VALS0))
                    else
                        garr[3*(ct.nat) + FIND] += sum(VALS0) 
                    end
                    
                end #end FIND
            end
            println()            
            if scf
                #println("IF SCF ", scf, " " , sum(abs.(g_ew[1,:][:])))
                #println("ONLY EW")
                garr += g_ew[1,:][:]
                #garr = g_ew[1,:][:]
            end
                
        end #end begin

    end #ends else
    #println("end else")
    #    println()
        
    x, stress = reshape_vec(garr, ct.nat)
    f_cart = -1.0 * x
    f_cart = f_cart * inv(ct.A)' / nspin
    stress = -stress / abs(det(ct.A)) /nspin

#    println("before ")
#    println(f_cart)
#    println()
#    println(stress)


    #important if there is symmetry
    f_cart,stress = symmetrize_vector_tensor(f_cart,stress, ct)
    
    
#    println(f_cart)
    
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


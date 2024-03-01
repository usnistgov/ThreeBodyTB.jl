using ..TB:tb_crys_sparse
using ..Symmetry:get_kgrid_sym
using ..Symmetry:symmetrize_vector_tensor
using ..CalcTB:calc_tb_LV_sparse

using SparseArrays

#contains sparse matrix versions of Force/Stress routines.
#actually uses sparse matrix multiplcation, so there should be a speedup in the limit that the Hamiltonian is actually sparse.
#also uses less memory

function ham_SINGLE_sparse(x :: Vector, ct, database, dontcheck, repel, DIST, FloatX, nwan, nr, atom, INDH, INDS, ret)

    begin

        begin
            T=typeof(x[1])

            if atom >= 1

                x_r, x_r_strain = reshape_vec_SINGLE(x, atom, ct.nat, strain_mode=false)
                A = FloatX.(ct.A) * (I(3))
            else
                x_r, x_r_strain = reshape_vec_SINGLE(x, -1, ct.nat, strain_mode=true)
                A = FloatX.(ct.A) * (I(3) + x_r_strain)            
            end
        end        
        
        crys_dual = makecrys( A , ct.coords + x_r, ct.types, units="Bohr", type=eltype(x))
        H, S = calc_tb_LV_sparse(crys_dual, database; verbose=false, var_type=T, use_threebody=true, use_threebody_onsite=true, gamma=zeros(T, ct.nat,ct.nat) , check_frontier= !dontcheck, repel=repel, DIST=DIST, retmat=true, atom = atom)
        
        
        #we turn sparse matricies back into dense array + index , for the purpose of running ForwardDiff.jacobian on it.
        begin

            c=0            
            for n = 1:length(H)
            
                I,J,data = findnz(H[n])
                push!(INDH, [I,J])
                ret[(1+c):(c+length(data))] = data
                c += length(data)
            end
            for n = 1:length(S)
                I,J,data = findnz(S[n])
                push!(INDS, [I,J])
                ret[(1+c):(c+length(data))] = data
                c += length(data)                
            end
            
            
        end
        
    end


    return ret
end

"""
    function sparsify_jac(jac, INDH, INDS, nwan)

A pain of using ForwardDiff.jacobian on a function that returns a sparse matrix Hamiltonian is that ForwardDiff wants a dense return. So we have to 
do a lot of work returning a dense matrix version of the sparse matrix and the index of the non-zero entries. We then sparsify it here.
"""
function sparsify_jac(jac, INDH, INDS, nwan)


    Harr = []
    Sarr = []
    for find = 1:size(jac,2)
        H = []
        c=0
        for n = 1:length(INDH)
            I,J = INDH[n]
            push!(H, droptol!(sparse(I,J,jac[(1+c):(c+length(I)),find], nwan, nwan), 1e-8))
            c+=length(I)
        end
        push!(Harr, H)
        S = []
        for n = 1:length(INDS)
            I,J = INDS[n]

            push!(S, droptol!(sparse(I,J,jac[(1+c):(c+length(I)),find], nwan, nwan), 1e-8))

            c+=length(I)
            
        end
        push!(Sarr, S)
    end

    return Harr, Sarr
end
    
"""
    function get_energy_force_stress_fft_LV_sym_SINGLE(tbc::tb_crys_sparse)

Main sparse matrix force/stress calculator. Should be faster for large
systems, due to faster matrix/vector and matrix/matrix multiply
(assuming sparscity).
Also requires much less memory.
"""
function get_energy_force_stress_fft_LV_sym_SINGLE(tbc::tb_crys_sparse, database; do_scf=false, smearing = 0.01, grid = missing, e_den0=missing, vv = missing, nspin = 1, repel=true)




    do_scf = true
    
    
    FloatX = Float64
    ct = deepcopy(tbc.crys)

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

    
    tooshort, energy_tot = safe_mode_energy(tbc.crys, database, DIST=DIST)
    
    
    if tooshort 
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


    else 

        
        scf = database["scf"]

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
                
                error_flag = false
                if do_scf
                    
                    energy_tot, efermi, e_den, dq, VECTS, VALS, error_flag, tbcx  = scf_energy(tbc, smearing=smearing, grid=grid, e_den0=e_den0, conv_thr = 1e-6, nspin=nspin, verbose=false, use_sym = true)
                    
                else

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
            h1_sparse = sparse(h1)
            if nspin == 2
                h1spin_sparse = [sparse(h1spin[1,:,:][:,:]), sparse(h1spin[2,:,:][:,:])]
            else
                h1spin_sparse = [h1spin[1,:,:][:,:]]
            end

            
            OCCS = gaussian.(VALS.-efermi, smearing)
        end


        if database["scf"] == true
            scf = true
        else
            scf = false
        end
        
        dontcheck, sum_repel = calc_tb_LV_sparse(ct, database; check_only=true, use_threebody=true, use_threebody_onsite=true, DIST=DIST, verbose=false)

        dontcheck = dontcheck && sum_repel

        chunksize=min(15, 3*ct.nat + 6)        

        g_ew = zeros(FloatX, 1, 3*ct.nat + 6)
        size_ret_full = tbc.tb.nwan * tbc.tb.nwan * tbc.tb.nr 

        begin

            #ewald term, TODO - make sparse and better scaling algorithm
            EW = begin
                FN_ew = x->ew(x,ct,FloatX, dq)
                if scf
                    cfg = ForwardDiff.JacobianConfig(FN_ew, zeros(FloatX, 3*ct.nat + 6), ForwardDiff.Chunk{chunksize}())
                    g_ew = ForwardDiff.jacobian(FN_ew, zeros(FloatX, 3*ct.nat + 6) , cfg ) ::  Array{FloatX,2}
                end
            end

            HAM = begin


                memory = Dict()
                ATOM = 1
                INDH = []
                INDS = []

                nnz = 0
                for n = 1:tbc.tb.nr
                    I,J,nz = findnz(tbc.tb.H[n])
                    nnz += length(nz)
                end

                for n = 1:tbc.tb.nr
                    I,J,nz = findnz(tbc.tb.S[n])
                    nnz += length(nz)
                end

                ret = zeros(nnz*3)
                
                FN_ham = (ret,x)->ham_SINGLE_sparse(x,ct,database,dontcheck, repel, DIST, FloatX, tbc.tb.nwan, tbc.tb.nr, ATOM, INDH, INDS, ret)
                
                
                chunksize=3
                cfg3 = ForwardDiff.JacobianConfig(FN_ham, ret, zeros(FloatX, 3), ForwardDiff.Chunk{chunksize}())
                chunksize=6
                cfg6 = ForwardDiff.JacobianConfig(FN_ham, ret, zeros(FloatX, 6), ForwardDiff.Chunk{chunksize}())
                
            end
        end
        INDH = []
        INDS = []
        ATOM=1
        size_ret = tbc.tb.nwan * tbc.tb.nwan * tbc.tb.nr 
        println("Calculate Force / Stress")
        begin
            
            begin
                
                garr =zeros(Float64, 3*ct.nat+6)
                VALS0 = zeros(FloatX,  nk_red, nspin)
                max_occ = 1
                for a = 2:tbc.tb.nwan
                    if maximum(OCCS[:,a,:]) > 1e-6
                        max_occ = a
                    end
                end
                
            end
            


        end            

        
        begin
            gatom = zeros(FloatX,2*size_ret+1 ,3)
            gstress = zeros(FloatX,2*size_ret+1 ,6)

            twopi_i = 2*pi*im

            pVECTS = permutedims(Complex{FloatX}.(VECTS), [3,4,1,2])
            
        end
        
        for atom in 0:ct.nat
        
            if atom == 0
                print("...doing stress ")
            else
                if atom == 1
                    print("...doing atom $atom ")
                elseif atom == ct.nat
                    println("$atom.")
                else
                    print("$atom ")
                end
            end
            ATOM=atom
            
            INDH = []
            INDS = []

            ret .= 0.0
            if atom >= 1
                jac =ForwardDiff.jacobian(FN_ham,ret, zeros(FloatX, 3) , cfg3 ) ::  Array{FloatX,2}

                
                Harr, Sarr = sparsify_jac(jac, INDH, INDS, tbc.tb.nwan)
                find= 1:3
                
            else
                jac = ForwardDiff.jacobian(FN_ham, ret, zeros(FloatX, 6) , cfg6 ) ::  Array{FloatX,2}
                
#                println("done jac part stress")
                Harr, Sarr = sparsify_jac(jac, INDH, INDS, tbc.tb.nwan)
                
                find = 1:6
            end                    

            f_temp = zeros(Complex{Float64}, nthreads())
            for FIND in find
                f_temp .= 0.0
                
                
                for k = 1:nk_red
                    
                    begin
                        h_k = spzeros(Complex{FloatX},  tbc.tb.nwan, tbc.tb.nwan)
                        s_k = spzeros(Complex{FloatX},  tbc.tb.nwan, tbc.tb.nwan)
                    end


                    begin
                        kpoint = kpts[k,:]
                        kmat = repeat(kpoint', tbc.tb.nr,1)
                        exp_ikr = exp.(twopi_i * sum(kmat .* tbc.tb.ind_arr,dims=2))



                        for nr = 1:tbc.tb.nr
                            h_k += Harr[FIND][nr] * exp_ikr[nr] 
                            s_k += Sarr[FIND][nr] * exp_ikr[nr]
                        end
                    end

                    for spin = 1:nspin

                        h_k_t = h_k + (h1_sparse + h1spin_sparse[spin]) .* s_k
                        @fastmath @inbounds @threads for a = 1:max_occ #sparse matrix multiplies with dense vector.
                            id = threadid()
                            f_temp[id] += (kweights[k]*OCCS[k,a,spin])*((@view pVECTS[:,a,k,spin])'*(h_k_t - VALS[k, a,spin] * s_k)*(@view pVECTS[:,a,k,spin]))
                        end
                    end
                end
                
                if atom >= 1
                    
                    garr[3*(atom-1) + FIND] += real(sum(f_temp) )
                    
                else
                    
                    garr[3*(ct.nat) + FIND] += real(sum(f_temp) )
                end

            end 
#            println("end find")
        end

        if scf
            
            
            garr += g_ew[1,:][:]
            
        end
        
    end 
    #println("end else")
    
        
    x, stress = reshape_vec(garr, ct.nat)
    f_cart = -1.0 * x
    f_cart = f_cart * inv(ct.A)' / nspin
    stress = -stress / abs(det(ct.A)) /nspin

    f_cart,stress = symmetrize_vector_tensor(f_cart,stress, ct)
    
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


module DivideAndConquer

using LinearAlgebra
using SparseArrays
using ..TB:tb_crys_sparse
using ..TB:tb_sparse
using ..CrystalMod:orbital_index
using ..CrystalMod:makecrys
using ..Atomdata:atoms
using ..CrystalMod:distances_etc_3bdy_parallel_LV
using ..BandTools:gaussian
using ..Utility:cutoff_fn_fast

function get_hams(tbc::tb_crys_sparse, efermi; cut_dist=missing, smearing=0.01)

    if ismissing(cut_dist)
        cut_dist = tbc.cut_dist
    end

    #if precalculated cutoff distance is too short, recalc
    @time if cut_dist > (tbc.cut_dist + 1e-5)
        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel_LV(tbc.crys,cut_dist, 0.0,var_type=Float64, return_floats=true, keep_extra=true)
        tbc.cut_dist = cut_dist
        tbc.dist_arr = dist_arr[:,1][:]
        tbc.R_keep_ab = R_keep_ab
    end        

    @time H = spzeros(tbc.tb.nwan, tbc.tb.nwan)
    @time S = spzeros(tbc.tb.nwan, tbc.tb.nwan)
    @time for r = 1:tbc.tb.nr #flatten
        H += tbc.tb.H[r]
        S += tbc.tb.S[r]
    end

    
    @time ind2orb, orb2ind, etotal, nval = orbital_index(tbc.crys)
    ind2orb_reverse = Dict()
    @time for k in keys(ind2orb)
        ind2orb_reverse[ind2orb[k][1], ind2orb[k][3]] = k
    end

    cols, rows, vals = findnz(H)
    colsS, rowsS, valsS = findnz(S)


    
#    at_nz = Dict()
#    for a = 1:tbc.crys.nat
#        at_nz[a] = Set([a])
#    end
           
#    for (c,r) in zip(cols, rows)
#        at1, at_type1, orb_name1 = ind2orb[c]
#        at2, at_type2, orb_name2 = ind2orb[r]
#        push!(at_nz[at1], at2)
#        push!(at_nz[at2], at1)
    #    end

#    at_nz_sort = Dict()
    at_oldnew = Dict()
    at_newold = Dict()
    OTHER_ATOMS = []
    DIST_ARR = Dict()
    @time for a = 1:tbc.crys.nat
        #        at_nz_sort[a] = sort(collect(at_nz[a]))

        #push!(OTHER_ATOMS, Set())
        s = Set()
        for c = 1:size(tbc.R_keep_ab)[1]
            a1 = tbc.R_keep_ab[c,2]
            a2 = tbc.R_keep_ab[c,3]
            dist = tbc.dist_arr[c]
            if a1 == a || a2 == a
                if dist < cut_dist
                    push!(s, a1)
                    push!(s, a2)
                    if (a1,a2) in keys(DIST_ARR)

                        DIST_ARR[a1,a2] = min(dist, DIST_ARR[a1,a2])
                        DIST_ARR[a2,a1] = min(dist, DIST_ARR[a2,a1])
                    else
                        DIST_ARR[a1,a2] = dist
                        DIST_ARR[a2,a1] = dist
                    end
                    
                end
            end
        end

        #println("other atoms $a $(sort(collect(s)))")
        push!(OTHER_ATOMS, sort(collect(s)))
        
        at_oldnew[a] = Dict()
        at_newold[a] = Dict()
        for (anew, aa) in enumerate(OTHER_ATOMS[a])
            at_oldnew[a][aa] = anew
            at_newold[a][anew] = aa
        end
        
    end


    #row cols stuff
    println("rc1")
    ROWSCOLS = []
    ROWSCOLS_S = []
    @time for a in 1:tbc.crys.nat
        push!(ROWSCOLS, Set())
        push!(ROWSCOLS_S, Set())
    end
    println("rc2")
    @time  for (counter, (c,r)) in enumerate(zip(cols, rows))
        at1, at_type1, orb_name1 = ind2orb[c]
        at2, at_type2, orb_name2 = ind2orb[r]
        for a in 1:tbc.crys.nat
            if at1 in OTHER_ATOMS[a] && at2 in OTHER_ATOMS[a] 
                push!(ROWSCOLS[a], counter)
            end
        end
    end
    @time  for (counter, (c,r)) in enumerate(zip(colsS, rowsS))
        at1, at_type1, orb_name1 = ind2orb[c]
        at2, at_type2, orb_name2 = ind2orb[r]
        for a in 1:tbc.crys.nat
            if at1 in OTHER_ATOMS[a] && at2 in OTHER_ATOMS[a] 
                push!(ROWSCOLS_S[a], counter)
            end
        end
    end

    
    HARR = []
    SARR = []

    convert_small_to_big = []

    IND2ORB_SMALL = []

    println("a")
    @time for a = 1:tbc.crys.nat
        
        push!(convert_small_to_big, Dict())
        norb = 0
        for s in OTHER_ATOMS[a]
            norb += Int64(atoms[tbc.crys.stypes[s]].nwan / 2)
        end


        HH = zeros(norb, norb)
        SS = zeros(norb, norb)
        c_small = makecrys(tbc.crys.A, tbc.crys.coords[OTHER_ATOMS[a], :], tbc.crys.stypes[OTHER_ATOMS[a]])
#        println("c_small ", c_small.nat)
        ind2orb_small, orb2ind_small, etotal_small, nval_small = orbital_index(c_small)
        push!(IND2ORB_SMALL, ind2orb_small)
        
        ind2orb_small_reverse = Dict()
        for k in keys(ind2orb_small)
            ind2orb_small_reverse[ind2orb_small[k][1], ind2orb_small[k][3]] = k
        end
        #       println("keys ", keys(orb2ind_small))
        

        #for (c,r, h) in zip(cols, rows, vals)
        for counter in ROWSCOLS[a]
            at1, at_type1, orb_name1 = ind2orb[cols[counter]]
            at2, at_type2, orb_name2 = ind2orb[rows[counter]]

            anew1 = at_oldnew[a][at1]
            anew2 = at_oldnew[a][at2]
            cnew = ind2orb_small_reverse[anew1, orb_name1]
            rnew = ind2orb_small_reverse[anew2, orb_name2]
                
            HH[rnew, cnew] = vals[counter] #* cutoff_fn_fast(DIST_ARR[at1,at2], cut_dist-1.0, cut_dist) #smooth the hamiltonian cutoff
        end
        for counter in ROWSCOLS_S[a]
            at1, at_type1, orb_name1 = ind2orb[colsS[counter]]
            at2, at_type2, orb_name2 = ind2orb[rowsS[counter]]

            anew1 = at_oldnew[a][at1]
            anew2 = at_oldnew[a][at2]
            cnew = ind2orb_small_reverse[anew1, orb_name1]
            rnew = ind2orb_small_reverse[anew2, orb_name2]
                
            SS[rnew, cnew] = valsS[counter] # * cutoff_fn_fast(DIST_ARR[at1,at2], cut_dist-1.0, cut_dist)
        end
        push!(HARR, HH)
        push!(SARR, SS)


        
##        for (c,r, h) in zip(cols, rows, vals)
##            at1, at_type1, orb_name1 = ind2orb[c]
##            at2, at_type2, orb_name2 = ind2orb[r]
##
##            if at1 in OTHER_ATOMS[a] && at2 in OTHER_ATOMS[a] 
##                anew1 = at_oldnew[a][at1]
##                anew2 = at_oldnew[a][at2]
##                cnew = ind2orb_small_reverse[anew1, orb_name1]
##                rnew = ind2orb_small_reverse[anew2, orb_name2]
##                
##                HH[rnew, cnew] = h
##            end
##
##        end
        

        #=


        #S
        for (c,r, h) in zip(colsS, rowsS, valsS)
            at1, at_type1, orb_name1 = ind2orb[c]
            at2, at_type2, orb_name2 = ind2orb[r]
            #            println("c $c r $c h $h at1 $at1 at2 $at2")
            #            if at1 in at_nz[a] && at2 in at_nz[a]
            if at1 in OTHER_ATOMS[a] && at2 in OTHER_ATOMS[a] 
                
                anew1 = at_oldnew[a][at1]
                anew2 = at_oldnew[a][at2]
                cnew = ind2orb_small_reverse[anew1, orb_name1]
                rnew = ind2orb_small_reverse[anew2, orb_name2]
                
                SS[rnew, cnew] = h
            end
        end
=#
        for k in keys(ind2orb_small)
            #            println("k $k  $(typeof(convert_small_to_big[a]))  $(ind2orb_reverse[ind2orb_small[k][1], ind2orb_small[k][3]])")
#            println("k $k")
#            println("k $k ind2orb_small[k][1] $(ind2orb_small[k][1]) ")
#            println("at_newold $(at_newold[a][ind2orb_small[k][1]])")

            convert_small_to_big[a][k] = ind2orb_reverse[at_newold[a][ind2orb_small[k][1]] , ind2orb_small[k][3]]
        end
        
        

    end
#    return
    
    VALS = []
    VECTS = []
    println("eig")
    norb_max = 0
    @time for (h,s) in zip(HARR, SARR)
#        println("h $(sum(abs.(h - h'))) s $(sum(abs.(s - s')))")
        norb_max = max(norb_max, size(h)[1])
        vals, vects = eigen(0.5(h+h'),0.5*(s+s'));
        push!(VALS, vals)
        push!(VECTS, vects)
    end
    println("norb_max $norb_max")

    @time DENMAT = spzeros(tbc.tb.nwan, tbc.tb.nwan)
#    println("size DENMAT ", size(DENMAT))

    println("b")
    @time for (a, (vals, vects, s)) in enumerate(zip(VALS, VECTS, SARR))
#        println("a $a aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
#        denmat = zeros(Complex{Float64}, size(vects))
#        dind = (vals .<= efermi)
        denmat = vects * diagm(gaussian(vals .- efermi, smearing))  * vects'
#        println("trace denmat ", tr(denmat * s))
        c=0
        I = zeros(Int64, prod(size(vects)))
        J = zeros(Int64, prod(size(vects)))
        denmat_arr = zeros(Complex{Float64}, prod(size(vects)))

        for i = 1:size(vects)[1]
            at1, at_type1, orb_name1 = IND2ORB_SMALL[a][i]
            ii = convert_small_to_big[a][i]
            for j = 1:size(vects)[1]
                at2, at_type2, orb_name2 = IND2ORB_SMALL[a][j]

                c += 1
                jj = convert_small_to_big[a][j]
                I[c] = ii
                J[c] = jj
                
                denmat_arr[c] = denmat[i,j] * 0.5 * ( (at_oldnew[a][a] == at1) + (at_oldnew[a][a] == at2 ))

                #                if at_oldnew[a][a] == at1 && at_oldnew[a][a] == at2
#                    denmat_arr[c] = round(denmat[i,j], digits=15)
#                elseif at1 == at_oldnew[a][a] || at2 == at_oldnew[a][a]
#                    denmat_arr[c] = round(denmat[i,j] * 0.5, digits=15)
#                end

            end
        end

        DENMAT[1:maximum(I), 1:maximum(J)] += sparse(I,J,round.(denmat_arr, digits=12))
        
        
    end

    H = 0.5*(H + H')
    S = 0.5*(S + S')
    DENMAT = 0.5*(DENMAT + DENMAT')
    return HARR, SARR, H, S, VALS, VECTS, DENMAT

end



end #end module

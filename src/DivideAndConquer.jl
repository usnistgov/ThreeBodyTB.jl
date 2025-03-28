module DivideAndConquer

using LinearAlgebra
using SparseArrays
using ..TB:tb_crys_sparse
using ..TB:tb_sparse
using ..CrystalMod:orbital_index
using ..CrystalMod:makecrys
using ..Atomdata:atoms

function get_hams(tbc::tb_crys_sparse, efermi)

    H = spzeros(tbc.tb.nwan, tbc.tb.nwan)
    S = spzeros(tbc.tb.nwan, tbc.tb.nwan)
    for r = 1:tbc.tb.nr #flatten
        H += tbc.tb.H[r]
        S += tbc.tb.S[r]
    end

    
    ind2orb, orb2ind, etotal, nval = orbital_index(tbc.crys)
    ind2orb_reverse = Dict()
    for k in keys(ind2orb)
        ind2orb_reverse[ind2orb[k][1], ind2orb[k][3]] = k
    end

    cols, rows, vals = findnz(H)
    colsS, rowsS, valsS = findnz(S)


    
    at_nz = Dict()
    for a = 1:tbc.crys.nat
        at_nz[a] = Set([a])
    end
           
    for (c,r) in zip(cols, rows)
        at1, at_type1, orb_name1 = ind2orb[c]
        at2, at_type2, orb_name2 = ind2orb[r]
        push!(at_nz[at1], at2)
        push!(at_nz[at2], at1)
    end

    at_nz_sort = Dict()
    at_oldnew = Dict()
    at_newold = Dict()
    for a = 1:tbc.crys.nat
        at_nz_sort[a] = sort(collect(at_nz[a]))
        
        at_oldnew[a] = Dict()
        at_newold[a] = Dict()
        for (anew, aa) in enumerate(at_nz_sort[a])
            at_oldnew[a][aa] = anew
            at_newold[a][anew] = aa
        end
        
    end

    HARR = []
    SARR = []

    convert_small_to_big = []

    IND2ORB_SMALL = []
    for a = 1:tbc.crys.nat
        push!(convert_small_to_big, Dict())
        norb = 0
        for aa in at_nz_sort[a]
            norb += Int64(atoms[tbc.crys.stypes[aa]].nwan / 2)
        end
        HH = zeros(norb, norb)
        SS = zeros(norb, norb)
        c_small = makecrys(tbc.crys.A, tbc.crys.coords[at_nz_sort[a], :], tbc.crys.stypes[at_nz_sort[a]])

        ind2orb_small, orb2ind_small, etotal_small, nval_small = orbital_index(c_small)
        push!(IND2ORB_SMALL, ind2orb_small)
        
        ind2orb_small_reverse = Dict()
        for k in keys(ind2orb_small)
            ind2orb_small_reverse[ind2orb_small[k][1], ind2orb_small[k][3]] = k
        end
        #       println("keys ", keys(orb2ind_small))
        
        for (c,r, h) in zip(cols, rows, vals)
            at1, at_type1, orb_name1 = ind2orb[c]
            at2, at_type2, orb_name2 = ind2orb[r]
            #          println("c $c r $c h $h at1 $at1 at2 $at2")
            if at1 in at_nz[a] && at2 in at_nz[a]
                anew1 = at_oldnew[a][at1]
                anew2 = at_oldnew[a][at2]
                cnew = ind2orb_small_reverse[anew1, orb_name1]
                rnew = ind2orb_small_reverse[anew2, orb_name2]
                
                HH[rnew, cnew] = h
            end
        end
        push!(HARR, HH)

        #S
        for (c,r, h) in zip(colsS, rowsS, valsS)
            at1, at_type1, orb_name1 = ind2orb[c]
            at2, at_type2, orb_name2 = ind2orb[r]
            #            println("c $c r $c h $h at1 $at1 at2 $at2")
            if at1 in at_nz[a] && at2 in at_nz[a]
                anew1 = at_oldnew[a][at1]
                anew2 = at_oldnew[a][at2]
                cnew = ind2orb_small_reverse[anew1, orb_name1]
                rnew = ind2orb_small_reverse[anew2, orb_name2]
                
                SS[rnew, cnew] = h
            end
        end
        push!(SARR, SS)

        for k in keys(ind2orb_small)
            #            println("k $k  $(typeof(convert_small_to_big[a]))  $(ind2orb_reverse[ind2orb_small[k][1], ind2orb_small[k][3]])")
            println("k $k")
            println("k $k ind2orb_small[k][1] $(ind2orb_small[k][1]) ")
            println("at_newold $(at_newold[a][ind2orb_small[k][1]])")

            convert_small_to_big[a][k] = ind2orb_reverse[at_newold[a][ind2orb_small[k][1]] , ind2orb_small[k][3]]
        end
        
        

    end

    VALS = []
    VECTS = []
    for (h,s) in zip(HARR, SARR)
        vals, vects = eigen(h,s);
        push!(VALS, vals)
        push!(VECTS, vects)
    end


    DENMAT = spzeros(tbc.tb.nwan, tbc.tb.nwan)
    println("size DENMAT ", size(DENMAT))

    for (a, (vals, vects, s)) in enumerate(zip(VALS, VECTS, SARR))
        println("a $a aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
        denmat = zeros(Complex{Float64}, size(vects))
        dind = (vals .<= efermi)
        denmat = vects[:,dind] * vects[:,dind]'
        println("trace denmat ", tr(denmat * s))
        c=0
        I = zeros(Int64, prod(size(vects)))
        J = zeros(Int64, prod(size(vects)))
        denmat_arr = zeros(Complex{Float64}, prod(size(vects)))
        for i = 1:size(vects)[1]
            for j = 1:size(vects)[1]
                at1, at_type1, orb_name1 = IND2ORB_SMALL[a][i]
                at2, at_type2, orb_name2 = IND2ORB_SMALL[a][j]

                c += 1
                ii = convert_small_to_big[a][i]
                jj = convert_small_to_big[a][j]
                I[c] = ii
                J[c] = jj
#                println("add c $c i $i j $j ii $ii jj $jj")
                if at_oldnew[a][a] == at1 && at_oldnew[a][a] == at2
                    denmat_arr[c] = round(denmat[i,j], digits=15)
                elseif at1 == at_oldnew[a][a] || at2 == at_oldnew[a][a]
                    denmat_arr[c] = round(denmat[i,j] * 0.5, digits=15)
                end
#                println("add denmat_arr c $c ", denmat_arr[c])
            end
        end
#        println("sizes ", [size(I), size(J), size(denmat_arr)])
#        println("size2 ", size(sparse(I,J,denmat_arr)))
#        println(sparse(I,J,denmat_arr) )
#        println()
        DENMAT[1:maximum(I), 1:maximum(J)] += sparse(I,J,denmat_arr)
        
        
    end


    return HARR, SARR, H, S, VALS, VECTS, DENMAT

end



end #end module

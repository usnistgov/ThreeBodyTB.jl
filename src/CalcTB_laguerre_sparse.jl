using SparseArrays

using ..TB:make_tb_crys_sparse
using ..TB:make_tb_sparse

#contains sparse matrix implementations.

function core3b_sparse!(cind1,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_3, memory, DAT_ARR_3, cut_h, H,IND, counter, sym_dat1, sym_dat2, lmn31, lmn32)

#            println("cut_h $cut_h")
    
    @inbounds @simd for o2x = 1:norb[a2]
        o2 = orbs_arr[a2,o2x,1]
        s2 = orbs_arr[a2,o2x,2]
        sum2 = orbs_arr[a2,o2x,3]

        sym32 = symmetry_factor_int(s2,1,lmn32, one ) 

        for o1x = 1:norb[a1]

            o1 = orbs_arr[a1,o1x,1]
            s1 = orbs_arr[a1,o1x,2]
            sum1 = orbs_arr[a1,o1x,3]

            
            sym31 = symmetry_factor_int(s1,1,lmn31, one )    

                      
            temp = 0.0
            for i = 1:DAT_IND_ARR_3[t1,t2,t3, sum1, sum2,1]
                temp += memory[i] * DAT_ARR_3[t1,t2,t3, DAT_IND_ARR_3[t1,t2,t3, sum1, sum2,1+i]]
            end
            #if abs(temp * sym31*sym32*10^3 * cut_h) > 1e-100
#                println("add $o1 $o2 $cind1 ", temp * sym31*sym32*10^3 * cut_h)

            #H_thread[o1,o2,cind1] += temp * sym31*sym32*10^3 * cut_h
#            H_thread[cind1] += temp * sym31*sym32*10^3 * cut_h

            H[counter[cind1]] = temp * sym31*sym32*10^3 * cut_h
            IND[counter[cind1],1] = o1
            IND[counter[cind1],2] = o2
            counter[cind1] += 1
            
            # end
            #println("add $o1 $o2 $cind1   ", temp * sym31*sym32*10^3 * cut_h)
            
            #temp = temp * sym31*sym32*10^3 * cut_h
            
        end
    end
            
end

function core_onsite3b_sparse!(c_zero,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_onsite_3, memory, DAT_ARR_3, cut_o, H, IND, counter)

    #@inbounds @fastmath @simd     
    @inbounds @fastmath @simd for o1x = 1:norb[a1]

        o1 = orbs_arr[a1,o1x,1]
        s1 = orbs_arr[a1,o1x,2]
        sum1 = orbs_arr[a1,o1x,3]
        
        
        temp = 0.0
        for i = 1:DAT_IND_ARR_onsite_3[t1,t2,t3, sum1,1 ]
            temp += memory[i] * DAT_ARR_3[t1,t2,t3, DAT_IND_ARR_onsite_3[t1,t2,t3, sum1, 1+i]]
        end
#            println(DAT_IND_ARR_3[t1,t2,t3, sum1, sum2,1], ", $a1 $a2 $a3 | $o1x $o2x | temp $temp | sum(mem) ", sum(memory), " sum DAT_ARR " , sum(DAT_ARR_3[t1,t2,t3,sum1, sum2, :]))
#        if abs(temp *10^3 * cut_o) > 1e-7
#        println("$a1 $o1 o  $(temp *10^3 * cut_o)   $cut_o  memory $(memory[1])   $(DAT_ARR_3[t1,t2,t3, DAT_IND_ARR_onsite_3[t1,t2,t3, sum1, 1+1]])")
        #            H_thread[o1,o1,c_zero] += temp *10^3 * cut_o

        H[counter[c_zero]] += temp *10^3 * cut_o        
        IND[counter[c_zero],1] = o1
        IND[counter[c_zero],2] = o1
        counter[c_zero] += 1
        #        end
    end
    
end

function core_sparse!(cham, a1, a2, t1, t2, norb, orbs_arr, DAT_IND_ARR, lag_arr, DAT_ARR, cut_a, H, S, IND, counter, lmn, sym_dat, sym_datS)
    
    @inbounds @simd for o2x = 1:norb[a2]
        o2 = orbs_arr[a2,o2x,1]
        s2 = orbs_arr[a2,o2x,2]
        sum2 = orbs_arr[a2,o2x,3]
        #        t2 = orbs_arr[a2,o2x,3]
        for o1x = 1:norb[a1]
            o1 = orbs_arr[a1,o1x,1]
            s1 = orbs_arr[a1,o1x,2]
            sum1 = orbs_arr[a1,o1x,3]
            #            t1 = orbs_arr[a1,o1x,3]

            
            nind = DAT_IND_ARR[t1,t2,1,sum1,sum2,1]
            #DAT_IND_ARR[t1,t2,1,orbs_arr[a1,o1x,2],orbs_arr[a2,o2x,2],1]
            
            sym_dat[1] = 0.0
            sym_datS[1] = 0.0
            
            for n = 1:6 #DAT_IND_ARR[t1,t2,1,orbs_arr[a1,o1x,2],orbs_arr[a2,o2x,2],1]
                sym_dat[1] +=  lag_arr[n]*DAT_ARR[t1,t2,1,DAT_IND_ARR[t1,t2,1,sum1,sum2,n+1]  ]
                sym_datS[1] +=  lag_arr[n]*DAT_ARR[t1,t2,2,DAT_IND_ARR[t1,t2,2,sum1,sum2,n+1]  ]
            end

            if nind >= 12
                sym_dat[2] = 0.0
                sym_datS[2] = 0.0

                for n = 7:12 #DAT_IND_ARR[t1,t2,1,orbs_arr[a1,o1x,2],orbs_arr[a2,o2x,2],1]
                    sym_dat[2] +=  lag_arr[n-6]*DAT_ARR[t1,t2,1,DAT_IND_ARR[t1,t2,1,sum1,sum2,n+1]]
                    sym_datS[2] +=  lag_arr[n-6]*DAT_ARR[t1,t2,2,DAT_IND_ARR[t1,t2,2,sum1,sum2,n+1]]
                end
                if nind >= 18
                    sym_dat[3] = 0.0
                    sym_datS[3] = 0.0
                    for n = 13:18 #DAT_IND_ARR[t1,t2,1,orbs_arr[a1,o1x,2],orbs_arr[a2,o2x,2],1]
                        sym_dat[3] +=  lag_arr[n-12]*DAT_ARR[t1,t2,1,DAT_IND_ARR[t1,t2,1,sum1,sum2,n+1]]
                        sym_datS[3] +=  lag_arr[n-12]*DAT_ARR[t1,t2,2,DAT_IND_ARR[t1,t2,2,sum1,sum2,n+1]]
                    end
                end
            end

            H[cham][counter[cham]] = symmetry_factor_int(s1, s2, lmn, sym_dat)  * cut_a
            S[cham][counter[cham]] = symmetry_factor_int(s1, s2, lmn, sym_datS)  * cut_a
            IND[cham][counter[cham],1] = o1
            IND[cham][counter[cham],2] = o2
            counter[cham] += 1
            #            H[ o1, o2] = symmetry_factor_int(s1, s2, lmn, sym_dat)  * cut_a
            #            S[ o1, o2] = symmetry_factor_int(s1, s2, lmn, sym_datS)  * cut_a
            
        end
    end
end

function core_onsite_sparse!(c_zero, a1, a2, t1, t2, norb, orbs_arr, DAT_IND_ARR_O, lag_arr, DAT_ARR, cut_on, H, IND, counter,  lmn, sym_dat1, sym_dat2)

    #    sym_dat1[1] = 1.0
    #    sym_dat1[2] = 0.0
    #    sym_dat1[3] = 0.0

    #    sym_dat2[1] = 1.0
    #    sym_dat2[2] = 0.0
    #    sym_dat2[3] = 0.0
    
    
    @inbounds @simd    for o2x = 1:norb[a1]
        o2 = orbs_arr[a1,o2x,1]
        s2 = orbs_arr[a1,o2x,2]
        sum2 = orbs_arr[a1,o2x,3]
        for o1x = 1:norb[a1]

            #            sym_dat1[1] = 1.0
            #            sym_dat2[1] = 1.0

            o1 = orbs_arr[a1,o1x,1]
            s1 = orbs_arr[a1,o1x,2]
            sum1 = orbs_arr[a1,o1x,3]

            #            println("t1 $t2 t2 $t2 sum1 $sum1 sum2 $sum2 s1 $s1 s2 $s2 o1 $o1 o2 $o2 o1x $o1x o2x $o2x")

            
            nind = DAT_IND_ARR_O[t1,t2,sum1,sum2,1]
            #DAT_IND_ARR[t1,t2,1,orbs_arr[a1,o1x,2],orbs_arr[a2,o2x,2],1]
            

            temp1 = 0.0
            if o1x == o2x
                for n = 1:5 #DAT_IND_ARR[t1,t2,1,orbs_arr[a1,o1x,2],orbs_arr[a2,o2x,2],1]
                    temp1 +=  lag_arr[n]*DAT_ARR[t1,t2,1,DAT_IND_ARR_O[t1,t2,sum1,sum2,n+1]  ]
                end
            end

            temp2 = 0.0
            if sum1 != sum2
                for n = 1:5 #DAT_IND_ARR[t1,t2,1,orbs_arr[a1,o1x,2],orbs_arr[a2,o2x,2],1]
                    temp2 +=  lag_arr[n]*DAT_ARR[t1,t2,1,DAT_IND_ARR_O[t1,t2,sum1,sum2,n+1]  ]
                end
                temp2= temp2* symmetry_factor_int(s1, 1, lmn, one)*symmetry_factor_int(s2, 1, lmn, one)
            end

            temp3 = 0.0
            if (sum1 == 2 && sum2 == 2) || (sum1 == 3 && sum2 == 3)
                for n = 1:5 #DAT_IND_ARR[t1,t2,1,orbs_arr[a1,o1x,2],orbs_arr[a2,o2x,2],1]
                    temp3 +=  lag_arr[n]*DAT_ARR[t1,t2,1,DAT_IND_ARR_O[t1,t2,sum1,sum2,n+1+5]  ]
                end
                temp3 = temp3 * symmetry_factor_int(s1, 1, lmn, one)*symmetry_factor_int(s2, 1, lmn, one)

            end

            #            H[c_zero][ o1, o2] += (temp1 + temp2 + temp3)  * cut_on

            H[c_zero][counter[c_zero]] = (temp1 + temp2 + temp3)  * cut_on
            
            IND[c_zero][counter[c_zero],1] = o1
            IND[c_zero][counter[c_zero],2] = o2
            counter[c_zero] += 1
            
            #if abs((temp1 + temp2)) > 1e-5
            #    println("$c_zero, $a1, $a2, $t1, $t2, $(temp1 + temp2)   $cut_on")
            #end
            
        end
    end
end

"""
    function calc_tb_LV_sparse(crys::crystal)

Main sparse matrix tight binding Hamiltonian calculator.
"""
function calc_tb_LV_sparse(crys::crystal, database=missing; reference_tbc=missing, verbose=false, var_type=missing, use_threebody=true, use_threebody_onsite=true, gamma=missing,background_charge_correction=0.0,  screening=1.0, set_maxmin=false, check_frontier=true, check_only=false, repel = true, DIST=missing, tot_charge=0.0, retmat=false, atom = -1)


    #        verbose=true
    #    println("test")
    #    @time twobody(10)

    
    At = crys.A'
    
    if verbose
        println()
        println("-----")
        println("Construct tight-binding model from crystal structure")
        println()
    end

    if ismissing(var_type)
        var_type=Float64
    end

    if verbose  println("LV $var_type")  end
    
    if ismissing(database)
        println("missing database, creating empty tbc")
        repel = false
        #    else
        #        println(keys(database))
    end
    
    if ismissing(reference_tbc)
        prepare_for_fitting = false
    else
        prepare_for_fitting = true
    end
    

    
    ind2orb, orb2ind, etotal, nval = orbital_index(crys)
    if verbose println("distances") end
    begin
        
        use_dist_arr = true
        if !ismissing(DIST)
            use_dist_arr = false
            R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3 = DIST

        else
            if (use_threebody || use_threebody_onsite ) && !ismissing(database)
                R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel_LV(crys,cutoff2X, cutoff3bX,var_type=var_type, return_floats=false)
                DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3

            else
                R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel_LV(crys,cutoff2X, 0.0,var_type=var_type, return_floats=false)
                DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3             
            end
        end

        
        
        within_fit = true
        
        if !ismissing(database)
            for key in keys(dmin_types)
                for key2 in keys(database)
                    if key == Set(key2)
                        if dmin_types[key] < database[key2].min_dist*1.0199 && length(key2) == 2 && var_type == Float64
                            println("WARNING : structure has 2body distances less than or close to min fitting data distances, may result in errors")
                            println(key," " ,key2, " : ", dmin_types[key], " <~ ", database[key2].min_dist, "   key $key key2 $key2")
                            #println(key," " ,key2, " : ", dmin_types[key], " <~ ", database[key2].min_dist)
                            within_fit = false
                        end
                    end
                end
            end
            
            c_zero_ref=1
            if !(ismissing(reference_tbc))
                if size(reference_tbc.tb.ind_arr)[1] > 1
                    c_zero_ref = reference_tbc.tb.r_dict[[0,0,0]]
                end
            end
        end
    end

    
    if verbose println("check_frontier  $check_frontier") end
    if !ismissing(database) && check_frontier
        #    if false
        #        println("check ---------------------------------")
        violation_list, vio_bool, repel_vals = calc_frontier(crys, database, test_frontier=true, diststuff=DIST, verbose=verbose, var_type=var_type, use_threebody=use_threebody)
        #        println("repel ", repel_vals)
        if vio_bool == false 
            within_fit = false
        end
    else
        repel_vals = zeros(var_type, crys.nat)
    end

#    println("stuff")
    begin
    
        
        #    println("repel_vals ", repel_vals)
        
        if check_only==true
            #        println("repel_vals ", repel_vals)
            return within_fit , sum(abs.(repel_vals)) < 1e-12
        end
        
        nwan = length(keys(ind2orb))

        nkeep=size(R_keep)[1]
        #    nkeep2=size(R_keep2)[1]    
        #    println("nkeep, $nkeep, nkeep2, $nkeep2")
        

        if verbose; println("memory"); end
        begin 
            H = SparseMatrixCSC{var_type, Int64}[]
            S = SparseMatrixCSC{var_type, Int64}[]
            
            for k in 1:nkeep
                push!(H, spzeros(nwan, nwan))
                push!(S, spzeros(nwan, nwan))
            end
            Hzero=H[c_zero]
            Szero=S[c_zero]
        end
        
        ind_arr = zeros(Int64, nkeep, 3)
        ind_arr[:,:] = R_keep[:,2:4]

        #    println("ind_arr")
        #    for c in 1:size(ind_arr)[1]
        #        println("$c ", ind_arr[c,:])
        #    end
        
        r_dict = make_rdict(ind_arr)
        
        norb = zeros(UInt16, crys.nat)
        orbs = zeros(UInt16, crys.nat, 1+3+5+7)
        sorbs = zeros(UInt16, crys.nat, 1+3+5+7)
        sumorbs = zeros(UInt16, crys.nat, 1+3+5+7)
        for a = 1:crys.nat
            ox = orb2ind[a]
            norb[a] = length(ox)            
            for (co,o) = enumerate(orb2ind[a])
                a1,t,s = ind2orb[o]
                sumO = summarize_orb_num(s)
                orbs[a,co] = o
                
                sorbs[a,co] = orb_num(s)
                sumorbs[a,co] = sumO + 1
                #            println("first ", [a1,t,s,sumO,o])
            end
        end

        
        warned = false
        warned_onsite = false

        nkeep_ab = size(R_keep_ab)[1]






#        println("not ismissing")
        if !ismissing(database)

            
            begin
                
                types_dict = Dict()
                types_dict_reverse = Dict()
                types_counter = 0
                types_arr = zeros(Int64, crys.nat)
                for a = 1:crys.nat
                    if !(crys.stypes[a]  in keys(types_dict))
                        types_counter += 1
                        types_dict[crys.stypes[a]] = types_counter
                        types_dict_reverse[types_counter] = crys.stypes[a]
                    end
                    types_arr[a] = types_dict[crys.stypes[a]]
                end

                begin 
                    orbs_arr = zeros(Int64, crys.nat, 9, 3)
                    DAT_IND_ARR = zeros(Int64, types_counter, types_counter, 2,4,4, 33 )
                    DAT_ARR = zeros(Float64, types_counter, types_counter, 2, 164)

                    DAT_IND_ARR_O = zeros(Int64, types_counter, types_counter, 4,4, 33 )


                    cutoff_arr = zeros(crys.nat, crys.nat, 2)
#                    cutoff_arr3 = zeros(crys.nat, crys.nat, crys.nat)
                    get_cutoff_pre = Dict()       
                    s = Set(crys.stypes)
                end
                for s1 in s
                    for s2 in s
                        get_cutoff_pre[(s1,s2)] = get_cutoff(s1,s2)[:] 
                        for s3 in s 
                            get_cutoff_pre[(s1,s2,s3)] = get_cutoff(s1,s2,s3)[1] 
                        end
                    end
                end
                
                for a = 1:crys.nat
                    ta = crys.stypes[a]
                    for b = 1:crys.nat
                        tb = crys.stypes[b]            
                        cutoff_arr[a,b,:] = get_cutoff_pre[(ta,tb)][:]
#                        for c = 1:crys.nat
#                            tc = crys.stypes[c]            
#                            cutoff_arr3[a,b,c] =  get_cutoff_pre[(ta,tb,tc)]
#                        end
                    end
                end

                badlist = Set()
                
                for c1 = 1:types_counter
                    for c2 = 1:types_counter
                        t1 = types_dict_reverse[c1]
                        t2 = types_dict_reverse[c2]
                        if (t1,t2) in keys(database)
                            coef = database[(t1,t2)]
                            indH, indS, inH, inS = coef.inds_int[[t1,t2]]
                            DAT_IND_ARR[c1,c2,1,:,:,2:33] = indH
                            DAT_IND_ARR[c1,c2,2,:,:,2:33] = indS
                            DAT_IND_ARR[c1,c2,1,:,:,1] = inH
                            DAT_IND_ARR[c1,c2,2,:,:,1] = inS
                            DAT_ARR[c1,c2,1,1:coef.sizeH] = coef.datH
                            DAT_ARR[c1,c2,2,1:coef.sizeS] = coef.datS
                            indO, inO = coef.inds_int[[t1]]
                            DAT_IND_ARR_O[c1,c2,:,:,1] = inO
                            DAT_IND_ARR_O[c1,c2,:,:,2:33] = indO
                        else
                            println("WARNING, ",(t1,t2), " database not found")
                            within_fit = false
                            push!(badlist, (t1,t2))
                            println("badlist ", badlist)
                        end
                    end
                end

                if use_threebody || use_threebody_onsite
                    
                    begin 
                        DAT_IND_ARR_3 = zeros(Int64, types_counter, types_counter, types_counter, 4,4, 33 )
                        DAT_ARR_3 = zeros(Float64, types_counter, types_counter, types_counter, 224)
                    end
                    
                    for c1 = 1:types_counter
                        t1 = types_dict_reverse[c1]
                        for c2 = 1:types_counter
                            t2 = types_dict_reverse[c2]
                            for c3 = 1:types_counter
                                t3 = types_dict_reverse[c3]
                                if (t1,t2,t3) in keys(database)
                                    cdat = database[(t1,t2,t3)]
                                    (cindX, nindX) = cdat.inds_int[[t1,t2,t3]]
                                    DAT_IND_ARR_3[c1,c2,c3,:,:,1] = nindX
                                    DAT_IND_ARR_3[c1,c2,c3,:,:,2:33] = cindX
                                    DAT_ARR_3[c1,c2,c3, 1:length(cdat.datH)] = cdat.datH
                                else
                                    println("WARNING, ",(t1,t2,t3), " database not found")
                                    within_fit = false
                                    push!(badlist, (t1,t2,t3))
                                    println("badlist ", badlist)
                                    
                                end
                            end
                        end
                    end
                end
                if use_threebody_onsite
                    DAT_IND_ARR_onsite_3 = zeros(Int64, types_counter, types_counter, types_counter, 4, 33 )
                    #            DAT_ARR_onsite_3 = zeros(Float64, types_counter, types_counter, types_counter, 32)
                    for c1 = 1:types_counter
                        t1 = types_dict_reverse[c1]
                        for c2 = 1:types_counter
                            t2 = types_dict_reverse[c2]
                            for c3 = 1:types_counter
                                t3 = types_dict_reverse[c3]
                                if !((t1,t2,t3) in badlist)
                                    cdat = database[(t1,t2,t3)]
                                    (cindX, nindX) = cdat.inds_int[[t1,t2,t3,:O]]
                                    DAT_IND_ARR_onsite_3[c1,c2,c3,:,1] = nindX
                                    DAT_IND_ARR_onsite_3[c1,c2,c3,:,2:33] = cindX
                                    #                        DAT_ARR_onsite_3[c1,c2,c3, 1:length(cdat.datO)] = cdat.datH
                                end
                            end
                        end
                    end
                end

                
                for a = 1:crys.nat
                    t1 = crys.stypes[a]
                    for o1x = 1:norb[a]
                        o1 = orbs[a,o1x]
                        s1 = sorbs[a,o1x]
                        sum1 = sumorbs[a,o1x]
                        orbs_arr[a,o1x,1] = o1
                        orbs_arr[a,o1x,2] = s1
                        orbs_arr[a,o1x,3] = sum1
                        #                println("$t1 norb[a] ", norb[a], ", orbs arr $a $o1x | $o1 $s1 $sum1")
                    end
                end

                begin 
                    lag_arr_TH = zeros(var_type, 6, nthreads())
                    lmn_arr_TH = zeros(var_type, 3, nthreads())
                    sym_arr_TH = zeros(var_type, 3, nthreads())
                    sym_arrS_TH = zeros(var_type, 3, nthreads())
                end
                
            end
#            println("end")
            
        end
    end    
    if verbose println("2body LV") end


    twobody_LV = begin

#        println("nonsense")
        begin
            #                println("nkeep $nkeep")
            counter_arr = zeros(Int64, nkeep)
        #    println("loop")
            for c = 1:nkeep_ab
                A1 = R_keep_ab[c,2]
                A2 = R_keep_ab[c,3]
                CHAM = R_keep_ab[c,7]
                counter_arr[CHAM] += norb[A1]*norb[A2]
            end
            counter_arr[c_zero] +=  nkeep_ab * 81  + nwan*2
            
            
            #                for k = 1:nkeep
            #                   println("nz $(counter_arr[k]) ")
            #                end
            Harr = Array{var_type}[]
            Sarr = Array{var_type}[]
            INDarr = Array{Int64,2}[]
#            println("push")
            for k in 1:nkeep
                push!(Harr, zeros(var_type, counter_arr[k]))
                push!(Sarr, zeros(var_type, counter_arr[k]))
                push!(INDarr, zeros(Int64, counter_arr[k],2))
            end
            
        end            
        
        
        #println("nkeep_ab $nkeep_ab")
 #       println("main twobody")
        begin
            #@time twobody(nkeep_ab)
            counter_arr .= 1
            for c = 1:nkeep_ab
                begin
                    id = threadid()
                    lmn_arr = lmn_arr_TH[:,id]
                    sym_arr = sym_arr_TH[:,id]
                    sym_arrS = sym_arrS_TH[:,id]
                    lag_arr = lag_arr_TH[:,id]
                    
                    cind = R_keep_ab[c,1]
                    cham = R_keep_ab[c,7]
                    #                        Hc = H[cham]
                    #                        Sc = S[cham]
                    a1 = R_keep_ab[c,2]
                    a2 = R_keep_ab[c,3]

                    if atom > 0 && !(a1 == atom || a2 == atom)
                        continue
                    end
                    #                        println("atom $atom a1 $a1 a2 $a2")
                    rind1 = ind_arr[cham,1:3]

                    t1 = types_arr[a1]
                    t2 = types_arr[a2]

                    if use_dist_arr
                        dist_a = dist_arr[c,1]
                        lmn_arr[1] = dist_arr[c,2]
                        lmn_arr[2] = dist_arr[c,3]
                        lmn_arr[3] = dist_arr[c,4]
                        cut_a = dist_arr[c,5]
                        cut_on = dist_arr[c,6]
                    else
                        dist_a, lmn_arr[:] = get_dist(a1,a2, R_keep_ab[c,4:6], crys, At)
                        cutoff2, cutoff2on = cutoff_arr[a1,a2,:]
                        cut_a = cutoff_fn_fast(dist_a, cutoff2 - cutoff_length, cutoff2)
                        cut_on = cutoff_fn_fast(dist_a, cutoff2on - cutoff_length, cutoff2on)                            
                    end
                    
                    
                    if dist_a <  1e-5    # true onsite

                        for o1x = 1:norb[a1]
                            o1 = orbs_arr[a1,o1x,1]
                            sum1 = orbs_arr[a1,o1x,3]
                            t1symbol =   crys.stypes[a1]
                            a1a,t1a,s1 = ind2orb[o1x]
                            
                            (h,s) = calc_onsite(t1symbol,sum1,sum1, database)
                            #                                S[o1, o1, c_zero] += s 
                            #                                H[o1, o1, c_zero] += h

                            #                                Szero[o1,o1] += s
                            #                                Hzero[o1,o1] += h


                            Harr[cham][counter_arr[cham]] = h

                            Sarr[cham][counter_arr[cham]] = s
                            INDarr[cham][counter_arr[cham],1] = o1
                            INDarr[cham][counter_arr[cham],2] = o1
                            counter_arr[cham] += 1
                            if repel
#                                println("repel")
                                #H[o1, o1, c_zero] += repel_vals[a1] * 0.1
                                #H[c_zero][o1, o1] += repel_vals[a1] * 0.1
                                Harr[cham][counter_arr[cham]] += repel_vals[a1] * 0.1
                                #println("add repel $a1 $o1 ", repel_vals[a1] * 0.1)
                            end
                            
                        end
                        
                        
                    else #normal case
                        

                        laguerre_fast!(dist_a, lag_arr)
                        core_sparse!(cham, a1, a2, t1, t2, norb, orbs_arr, DAT_IND_ARR, lag_arr, DAT_ARR, cut_a, Harr, Sarr, INDarr, counter_arr, lmn_arr, sym_arr, sym_arrS)
                        core_onsite_sparse!(c_zero, a1, a2, t1, t2, norb, orbs_arr, DAT_IND_ARR_O, lag_arr, DAT_ARR, cut_on, Harr, INDarr, counter_arr, lmn_arr, sym_arr, sym_arrS)

                        
                    end
                    
                end
                
            end
        end
    end

    

    #println("sparsify")
    for k in  1:nkeep
        #            println(size(INDarr[k][:,1]), " " , size(INDarr[k][:,2]), size(Harr[k]))
        #            println("k $k   max ", maximum(INDarr[k]),  " min ", minimum(INDarr[k]), " vs nwan ", nwan)
        if retmat
            H[k] += sparse(INDarr[k][1:counter_arr[k]-1,1], INDarr[k][1:counter_arr[k]-1,2], Harr[k][1:counter_arr[k]-1] .+ 1.234, nwan, nwan ) #this adding 1.234 is a trick to prevent the sparse array from deciding that the dual number is zero accidently if the real part is zero
            S[k] += sparse(INDarr[k][1:counter_arr[k]-1,1], INDarr[k][1:counter_arr[k]-1,2], Sarr[k][1:counter_arr[k]-1] .+ 1.234, nwan , nwan)
        else
            H[k] += sparse(INDarr[k][1:counter_arr[k]-1,1], INDarr[k][1:counter_arr[k]-1,2], Harr[k][1:counter_arr[k]-1], nwan, nwan )
            S[k] += sparse(INDarr[k][1:counter_arr[k]-1,1], INDarr[k][1:counter_arr[k]-1,2], Sarr[k][1:counter_arr[k]-1], nwan , nwan)
        end
    end


    
    
        #                dist_a, lmn = get_dist(a1,a2, R_keep_ab[c,4:6], crys, At)

        #            dist_a_old = dist_arr[a1,a2,cind,1]
        
        
        #                if (dist_a > cutoff2Xt || dist_a < 1e-5)
        #                    continue
        #                end
        
        
        
        


    threebdy_LV = begin

        #println("threebody prepare")
        if use_threebody || use_threebody_onsite

            #            println("nonsense3")
            counter_arr .= 0

            begin
#                println("loop")
                if use_threebody
                    
                    for  counter = 1: (size(array_ind3)[1] ) #add threads back
                        a1 = array_ind3[counter,1]
                        a2 = array_ind3[counter,2]
                        a3 = array_ind3[counter,3]
                        cind1 = array_ind3[counter,4]
                        #                println("add $cind1  $a1 $a2 $a3")
                        counter_arr[cind1] += norb[a1]*norb[a2]
                    end
                end
                if use_threebody_onsite
                    counter_arr[c_zero] += 9 * size(array_ind3)[1] + nwan*2
                end
                
                #            for i = 1:nkeep
                #                println("ca $i ", counter_arr[i])
                #            end
                
                Harr3 = Array{var_type}[]
                INDarr3 = Array{Int64,2}[]
#                println("push")
                for k in 1:nkeep
                    push!(Harr3, zeros(var_type, counter_arr[k]))
                    push!(INDarr3, zeros(Int64, counter_arr[k],2))
                end
                counter_arr .= 1
            end
        end
        
        begin
            #                @time H_thread = zeros(var_type,   nwan, nwan,  nkeep,  nthreads() )
            #                @time H_thread3 = zeros(var_type,   nwan, nwan,  nthreads() )
            #                H_thread2 = zeros(var_type,   9, 9,  nthreads() )

            hh = zeros(var_type, nthreads(), 3,3)

            skip = 0
            keep = 0

            lmn_arr_TH_1 = zeros(var_type, 3, nthreads())
            lmn_arr_TH_2 = zeros(var_type, 3, nthreads())
            lmn_arr_TH_3 = zeros(var_type, 3, nthreads())

            #                lag_arr_TH = zeros(var_type, 2, nthreads())
            #                lag_arr_TH = zeros(var_type, 2, nthreads())
            #                lag_arr_TH = zeros(var_type, 2, nthreads())

            memory_TH = zeros(var_type, 8, nthreads())
        end
        
        if use_threebody || use_threebody_onsite
            if verbose println("3body LV") end
            meta_count = []
            old = 1
            for  counter = 1: (size(array_ind3)[1]-1 ) #add threads back
                if array_ind3[counter,1] != array_ind3[counter+1,1]
                    push!(meta_count, old:counter)
                    #                        println("mc ", old:counter)
                    old = counter+1
                end
            end
            push!(meta_count, old:size(array_ind3)[1])
            

            
            #                for mc in meta_count #@threads 
            #                    for counter in mc
            #println("sorttime")
            sortind = sortperm(array_ind3[:,4])
            
            for  counterT = 1: (size(array_ind3)[1] ) #add threads back
                counter = sortind[counterT]
                id = threadid()


                
                a1 = array_ind3[counter,1]
                a2 = array_ind3[counter,2]
                a3 = array_ind3[counter,3]

                if atom > 0 && !(a1 == atom || a2 == atom || a3 == atom)
                    continue
                end

                
                t1s = crys.stypes[a1]
                t2s = crys.stypes[a2]
                t3s = crys.stypes[a3]

                if !haskey(database, (t1s, t2s, t3s))
                    continue
                end

                t1 = types_arr[a1]
                t2 = types_arr[a2]
                t3 = types_arr[a3]
                
                if (t1,t2,t3) in badlist
                    continue
                end
                
                lmn12 = lmn_arr_TH_1[:,id]
                lmn13 = lmn_arr_TH_2[:,id]
                lmn23 = lmn_arr_TH_3[:,id]

                sym_arr1 = sym_arr_TH[:,id]
                sym_arr2 = sym_arrS_TH[:,id]
                

                cind1 = array_ind3[counter,4]
                cind2 = array_ind3[counter,5]

                rind1 = ind_arr[cind1,1:3]
                rind2 = ind_arr[cind2,1:3]

                
                if use_dist_arr
                    dist12 = array_floats3[counter, 1]
                    dist13 = array_floats3[counter, 2]
                    dist23 = array_floats3[counter, 3]                    

                    for i = 1:3
                        lmn12[i] = array_floats3[counter, 3+i]
                        lmn13[i] = array_floats3[counter, 6+i]
                        lmn23[i] = array_floats3[counter, 9+i]
                    end
                    cut_h = array_floats3[counter,13]
                    cut_o = array_floats3[counter,14]
                else
                    dist12, lmn12[:] = get_dist(a1,a2, rind1, crys, At)
                    dist13, lmn13[:] = get_dist(a1,a3, rind2, crys, At)
                    dist23, lmn23[:] = get_dist(a2,a3, -rind1+rind2, crys, At)
                    
                    cutoffZZ = cutoff_arr[a1,a2,1]
                    #                    cutoff3 = cutoff_arr3[a1,a2,a3]
                    cutoff3 = get_cutoff_pre[(crys.stypes[a1],crys.stypes[a2],crys.stypes[a3])]

                    cut_ab = cutoff_fn_fast(dist12, cutoffZZ - cutoff_length, cutoffZZ)
                    cut_ab2 = cutoff_fn_fast(dist12, cutoff3 - cutoff_length, cutoff3)
                    cut_ac = cutoff_fn_fast(dist13, cutoff3 - cutoff_length, cutoff3)
                    cut_bc = cutoff_fn_fast(dist23, cutoff3 - cutoff_length, cutoff3)
                    
                    cut_h = cut_ab*cut_ac*cut_bc
                    cut_o = cut_ab2*cut_ac*cut_bc   
                    
                end

                #                    println("counter $counter $dist12 $dist13 $dist23 $cut_h $cut_o")

                memory = memory_TH[:,id]

                
                if use_threebody
                    laguerre_fast_threebdy!(dist12,dist13,dist23, t1==t2, t1 !=t2 && t1 != t3 && t2 != t3, memory)
                    #$core3!(cind1,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_3, memory, DAT_ARR_3, cut_h, H_thread, id,sym_arr1, sym_arr2, lmn13, lmn23)
                    #                            println("cuth $cut_h ",  [dist12,dist13,dist23], " ", exp(-1.0*dist13)*exp(-1.0*dist23)*1000 )
                    core3b_sparse!(cind1,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_3, memory, DAT_ARR_3, cut_h, Harr3[cind1], INDarr3[cind1] ,counter_arr, sym_arr1, sym_arr2, lmn13, lmn23)

                    
                    #                        core3a!(cind1,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_3, memory, DAT_ARR_3, cut_h, H_thread3,id, sym_arr1, sym_arr2, lmn13, lmn23)
                    #                        println("size ", size(H), size(H_thread2))
                    #                        println("x $a1 $a2 ", orbs_arr[a1,1,1]:orbs_arr[a1,norb[a1],1], " " , orbs_arr[a2,1,1]:orbs_arr[a2,norb[a2],1])

                    #    H[orbs_arr[a1,1,1]:orbs_arr[a1,norb[a1],1],orbs_arr[a2,1,1]:orbs_arr[a2,norb[a2],1], cind1] += @view H_thread2[1:norb[a1],1:norb[a2],1]

                end
                

                if use_threebody_onsite
                    laguerre_fast_threebdy_onsite!(dist12,dist13,dist23, t1==t2 && t1 == t3, memory)
                        core_onsite3b_sparse!(c_zero,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_onsite_3, memory, DAT_ARR_3, cut_o, Harr3[c_zero], INDarr3[c_zero], counter_arr )

#                    core_onsite3_sparse!(c_zero,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_onsite_3, memory, DAT_ARR_3, cut_o, Harr3, INDarr3, counter_arr, id)

                    #println("cuto $cut_o ",  [dist12,dist13,dist23], " ", exp(-1.0*dist13)*exp(-1.0*dist23)*1000, " ", exp(-1.0*dist13)*exp(-1.0*dist23)*1000*exp(-1.0*dist12))

                    #                    core_onsite3b_sparse!(c_zero,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_onsite_3, memory, DAT_ARR_3, cut_o, H)
                end
            end
        end
        #H[:,:,:] .+= sum(H_thread, dims=4)[ :,:,:]

        #           for id = 1:nthreads()
        #                println("id $id ", sum(H_thread[:,:,:,id]))
        #                println(H_thread[:,:,1,id])
        #                println()
        #            end
    end


#    println(INDarr3)
    

    
    if use_threebody || use_threebody_onsite
        #println("sparsify3")
        for k in  1:nkeep
            if retmat
                H[k] += sparse(INDarr3[k][1:counter_arr[k]-1,1], INDarr3[k][1:counter_arr[k]-1,2], Harr3[k][1:counter_arr[k]-1] .+ 0.1234, nwan, nwan )
            else
                H[k] += sparse(INDarr3[k][1:counter_arr[k]-1,1], INDarr3[k][1:counter_arr[k]-1,2], Harr3[k][1:counter_arr[k]-1], nwan, nwan )
            end
        end
    end    


    if retmat
        #return Harr, Sarr #, INDarr, Harr3, INDarr3
        #return

        return H, S
    end
    

    if verbose println("make") end
    if true
        #        println("typeof H ", typeof(H), " " , size(H), " S ", typeof(S), " " , size(S))
        #println("maketb")
        tb = make_tb_sparse( H , ind_arr, r_dict,  S)
        if !ismissing(database) && (haskey(database, "scf") || haskey(database, "SCF"))
            scf = database["scf"]
        else
            scf = false
        end
        #println("make")
        tbc = make_tb_crys_sparse(tb, crys, nval, 0.0, scf=scf, gamma=gamma, background_charge_correction=background_charge_correction, within_fit=within_fit, screening=screening)
        tbc.tot_charge = tot_charge
        tbc.nelec = tbc.nelec - tot_charge
    end
    if verbose 
        println("-----")
        println()
    end


    return tbc

end

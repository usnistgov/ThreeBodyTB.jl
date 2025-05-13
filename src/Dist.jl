
using LoopVectorization
using Base.Threads
using ..AtomicMod:atom
using ..Utility:cutoff_fn
using ..Utility:cutoff_fn_fast
using ..Atomdata:atoms
using ..Atomdata:get_cutoff


#const cutoff2X = 20.51 
#const cutoff3bX = 20.01 
#const cutoff_onX = 21.01

const cutoff2X = 18.51
const cutoff3bX = 13.01
const cutoff_onX = 18.01

const cutoff_length = 1.5


const aX2 = 1.0

const cutoff4X = 9.0

function  get_dist(a1,a2,Rvec, crys, At)
    
    v = At*( (@view crys.coords[a1,:])  - (@view crys.coords[a2,:]) + Rvec)
    dist = norm(v)
    
    return dist, v./(dist + 1e-30)
    
end

function  get_dist_only(a1,a2,Rvec, crys, At)
    
#    v = At*( (@view crys.coords[a1,:])  - (@view crys.coords[a2,:]) + Rvec)
#    dist = norm(v)
#    
#    return dist, v./(dist + 1e-30)
    return norm(At*( (@view crys.coords[a1,:])  - (@view crys.coords[a2,:]) + Rvec))
    
end


"""
    function distances_etc(crys, cutoff, cutoff2=missing)

Old distance function, no threebody only twobody.
Code primarly uses `distances_etc_3bdy_parallel`
"""
function distances_etc(crys, cutoff, cutoff2=missing)
    if ismissing(cutoff2)
        cutoff2=cutoff
    end
    
    R = get_grid(crys, 35.0)

    nr = (R[1]*2+1)*(R[2]*2+1)*(R[3]*2+1)
    
    dist_arr = zeros(crys.nat,crys.nat,nr,4)
                     
    R_f = zeros(3)
    lmn = zeros(3)
    dist = 0.0
    
    c=0

    c_zero = 0
    
    R_keep = []
    R_keep2 = []    
    for r1 = -R[1]:R[1]
        for r2 = -R[2]:R[2]
            for r3 = -R[3]:R[3]
                c+=1

                
                R_f[:] = Float64.([r1,r2,r3])
                found = false
                found2 = false                
                for a = 1:crys.nat
                    ta = crys.stypes[a]
                    for b = 1:crys.nat
                        tb = crys.stypes[b]            
                        dR = ( -crys.coords[a,:] .+ crys.coords[b,:] .+ R_f[:])'*crys.A
                        dist = sum(dR.^2)^0.5
#                        if dist < 7.0 
#                            println("calc $a $b [$r1 $r2 $r3] $dist $dR ")
#                        end
                        if dist > 1e-7
                            lmn[:] = dR/(dist )
                        else
                            lmn[:] .= 0.0
                        end
                        dist_arr[a,b,c,1] = dist
                        dist_arr[a,b,c,2:4] = lmn[:]

                        if dist < cutoff
                            found = true
                        end
                        if dist < cutoff2
                            found2 = true
                        end
                    end
                end
                if found
                    push!(R_keep, [c, r1,r2,r3])
                    if r1 == 0 && r2 == 0 && r3 == 0
                        c_zero = size(R_keep)[1]
                    end
                    
                end
                if found2
                    push!(R_keep2, [c, r1,r2,r3])
                end
            end
        end
    end

    return R_keep, R_keep2, dist_arr, c_zero
    
end
########################################################################################################################################################

"""
    function distances_etc_3bdy_parallel(crys, cutoff=missing, cutoff2=missing; var_type=Float64)

Finds atoms-atom and atom-atom-atom distances that are less than the cutoffs.

`    return R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3`

Where

- `R_keep` has the supercells that are necessary
- `R_keep_ab` has the 2body indexes that are less than the `cutoff`
- `array_ind3` has the 3body indexes
- `array_floats3` has the 3body actual distances and lmn.
- `dist_arr` has the actual 2body distances
- `c_zero` the index of cell `[0,0,0]`, the onsite cell
- `dmin_types` minimum distance between types of atoms.
- `dmin_types3` minimum 3body distances between types of atoms.

"""
function distances_etc_3bdy_parallel_old(crys, cutoff=missing, cutoff2=missing; var_type=Float64)

#    println("cutoff $cutoff $cutoff2")
    
    if ismissing(cutoff)
        cutoff = cutoff2X
    end
    
    if ismissing(cutoff2)
        cutoff2=cutoff3bX
    end
    threebody=true
    if cutoff2 < 1e-5
        threebody= false
    end

    dmin_types = Dict()
    dmin_types3 = Dict()
    for t1 = crys.stypes
        for t2 in crys.stypes
            dmin_types[Set((t1,t2))] = get_cutoff(t1,t2)[1]
            for t3 in crys.stypes
                dmin_types3[Set((t1,t2,t3))] = get_cutoff(t1,t2,t3)
            end
        end
    end

    

    R = get_grid(crys, 35.0)

    nr = (R[1]*2+1)*(R[2]*2+1)*(R[3]*2+1)
    
    dist_arr = zeros(var_type, crys.nat,crys.nat,nr,4)
                     
    R_f = zeros(var_type, 3, nthreads())
    
    lmn = zeros(var_type, 3)

    coords_ab = zeros(var_type, 3, crys.nat, crys.nat)
    for a = 1:crys.nat
        for b = 1:crys.nat
            coords_ab[:,a,b] = (crys.coords[a,:] .- crys.coords[b,:])' * crys.A 
        end
    end
    
#    c=0
    At = crys.A'

    S3 = (2*R[3]+1)
    S23 = (2*R[3]+1)*(2*R[2]+1)

    c_zero = 0

    Rind = zeros(Int64, nr, 3)

    found_arr = zeros(Bool, nr)
    found_arr[:] .= false
    
    @threads for r1 = -R[1]:R[1]
    #for r1 = -R[1]:R[1]    
        c1 = r1 + R[1] + 1
        id = threadid()
        R_f[1, id] = r1
        for (c2,r2) = enumerate(-R[2]:R[2])
            R_f[2, id] = r2
            for (c3,r3) = enumerate(-R[3]:R[3])
                c = c3 + (c2-1) * S3 + (c1-1)*S23

                Rind[c, :] .= [r1,r2,r3]

                
                R_f[3, id] = r3
                R_f2 = At*R_f

                found = false
                for a = 1:crys.nat
                    ta = crys.stypes[a]
                    for b = 1:crys.nat
                        tb = crys.stypes[b]            
#                        println(R_f2)
#                        println(coords_ab[:,a,b])
                        dR = coords_ab[:,a,b] + R_f2[:,id]
                        dist = (dR'*dR)^0.5

#                        println(typeof(dR))
#                        println(typeof(dist))
#                        println(size(dR))
#                        println(size(dist))
                        dist_arr[a,b,c,1] = dist
                        
                        if dist > 1e-7
                            dist_arr[a,b,c,2:4].= dR/(dist )
                        end
                        cutoff = get_cutoff(ta,tb)[1]
                        if dist < cutoff
                            found = true
                        end
                        
                    end
                end
                if r1 == 0 && r2 == 0 && r3 == 0
                    c_zero = c
                end
                if found
                    found_arr[c] = true
                end
            end
        end
    end

    R_reverse = Dict()
    for key in 1:size(Rind)[1]
#        println(key , " ", Rind[key,:])
        R_reverse[Rind[key,:]] = key
    end
        
    

    
    R_keep = zeros(Int64, 0, 4)
    R_dict = Dict()
    fcount =0
    for i in 1:nr
        if found_arr[i]
            c = Rind[i,:]
            R_keep = [R_keep; [0 c']]
            fcount += 1
            R_dict[c] = fcount
            if i == c_zero
                c_zero = fcount
            end
        end
    end

    R_keep_ab = zeros(Int64, crys.nat*crys.nat*nr, 7)
    
    ind_cutoff = Dict()
    ind_cutoff3bX = Dict()

    keep_counter = 0
    
    for a = 1:crys.nat
        ta = crys.stypes[a]
        for b = 1:crys.nat
            tb = crys.stypes[b]            
            ind = dist_arr[a,b,:,1] .> 1e-7
            
            dmin = minimum( dist_arr[a,b,ind,1])
            if dmin < dmin_types[Set((ta,tb))]
                dmin_types[Set((ta,tb))] = dmin
            end

            cutoff = get_cutoff(ta,tb)[1]
            ind2 = findall(dist_arr[a,b,:,1] .< cutoff)
            ind_cutoff[(a,b)] = deepcopy(ind2)
            for i in ind2
                keep_counter += 1
                r = Rind[i,:]
                ikeep = R_dict[r]
                
                R_keep_ab[keep_counter,:] = [i, a, b, r[1], r[2], r[3], ikeep ]
            end
            
            ind3 = findall(dist_arr[a,b,:,1] .< cutoff3bX)
            ind_cutoff3bX[(a,b)] = deepcopy(ind3)

        end
    end

    
    R_keep_ab = R_keep_ab[1:keep_counter,:]
    
    ############

    MEMCHUNK = min(nr*nr * crys.nat^3, 10000)
    TOTMEM = MEMCHUNK * ones(Int64, nthreads())

    AI3 = []
    AF3 = []
    COUNTER = zeros(Int64, nthreads())
    for i = 1:nthreads()
        array_ind3 = zeros(Int64, MEMCHUNK, 5)
        array_floats3 = zeros(var_type, MEMCHUNK , 14)
        push!(AI3, array_ind3)
        push!(AF3, array_floats3)
    end


        

    dmat = dist_arr[:,:,:,2:4] .* dist_arr[:,:,:, 1]
    if threebody
#        println("THREE")
        #        @threads for a = 1:crys.nat
        for a = 1:crys.nat        
        ta = crys.stypes[a]
            #        id = threadid()
            id = 1
        array_ind3X = AI3[id]
        array_floats3X = AF3[id]

        for b = 1:crys.nat
            tb = crys.stypes[b]
            cutoff = get_cutoff(ta,tb)[1]
            for c1 = ind_cutoff[(a,b)]

                dist_ab = dist_arr[a,b, c1 , 1]
                lmn_ab = dist_arr[a,b, c1  , 2:4]

                d_ab = dmat[a,b, c1,:]
                
                cut_ab = cutoff_fn(dist_ab, cutoff - cutoff_length, cutoff)



                for c = 1:crys.nat
                    
                    tc = crys.stypes[c]
                    cutoff3 = get_cutoff(ta,tb,tc)

                    cut_ab2 = cutoff_fn(dist_ab, cutoff3 - cutoff_length, cutoff3)
                    for c2 in ind_cutoff[(a,c)]

                        dist_ac = dist_arr[a,c, c2, 1]
                        if dist_ac > cutoff3
                            continue
                        end

#                        r1 = Rind[c1, :]
#                        r2 = Rind[c2, :]
#                        r3 = r1 - r2
#                        c12 = R_reverse[r3]

                        
                        Rdiff = @view(Rind[c1, :]) .- @view (Rind[c2, :])
                        if Rdiff in keys(R_reverse)
                            c12 = R_reverse[ Rdiff ]
                        else
                            continue
                        end
                            
#                        temp = ( d_ab .- dmat[a,c,c2, :])
                        
#                        dist_bc = (temp'*temp)
#                        if dist_bc > cutoff3_2
#                            continue
#                        end
#                        dist_bc = dist_bc^0.5

                        dist_bc = dist_arr[c, b, c12, 1]
                        if dist_bc > cutoff3
                            continue
                        end
                        lmn_bc = - (@view dmat[c,b,c12, :]) ./ dist_bc
                        
                        
#                        println(dist_bc_alt, " " , dist_bc)
                        
                        if dist_bc < cutoff3 && dist_ab > 1e-5 && dist_bc > 1e-5 && dist_ac > 1e-5
                        
                            lmn_ac = dist_arr[a,c, c2, 2:4]
#                            lmn_bc = - temp ./ dist_bc

                            

                            d0=dmin_types3[Set((ta,tb,tc))]
                            if (dist_bc + dist_ab + dist_ac)/3.0 < d0
                                dmin_types3[Set((ta,tb,tc))]  = (dist_bc + dist_ab + dist_ac)/3.0
                            end

                            cut_ac = cutoff_fn(dist_ac, cutoff3 - cutoff_length, cutoff3)
                            cut_bc = cutoff_fn(dist_bc, cutoff3 - cutoff_length, cutoff3)

                            COUNTER[id] += 1
                            if COUNTER[id] > TOTMEM[id] #need more memory
                                TOTMEM[id] += MEMCHUNK
                                array_ind3X = [array_ind3X;zeros(Int64, MEMCHUNK, 5)]
                                AI3[id] = array_ind3X
                                array_floats3X = [array_floats3X; zeros(var_type, MEMCHUNK , 14)]
                                AF3[id] = array_floats3X
                                
                            end

                            r = Rind[c1,:]
                            ikeep = R_dict[r]
                            
                            array_ind3X[COUNTER[id],:] = [a,b,c,ikeep,c2]
                            array_floats3X[COUNTER[id], :] = [dist_ab, dist_ac, dist_bc, lmn_ab[1],lmn_ab[2], lmn_ab[3], lmn_ac[1],lmn_ac[2], lmn_ac[3], lmn_bc[1],lmn_bc[2], lmn_bc[3], cut_ab*cut_bc*cut_ac, cut_ab2*cut_bc*cut_ac]

                                if COUNTER[id] == 761  ||  COUNTER[id] == 762 || COUNTER[id] == 763 || COUNTER[id] == 764

                                println(array_ind3X[COUNTER[id],:])
                            end
                            
#                            catch
#                                println(id, " " , COUNTER[id], " ", TOTMEM[id])
#                            end
                        end
                    end
                end
            end
        end
    end
    end
                
    array_ind3 = zeros(Int64, 0, 5)
    array_floats3 = zeros(var_type, 0 , 14)
    for (counter, i,f) in zip(COUNTER, AI3, AF3)
        array_ind3 = [array_ind3; i[1:counter, :]]
        array_floats3 = [array_floats3; f[1:counter,:]]
    end
    
    return R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3
#    return R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3

end





function distances_etc_3bdy_parallel(crys, cutoff=missing, cutoff2=missing; var_type=Float64, return_floats=true, shrink = 1.0)
    #    println("cutoff $cutoff $cutoff2")

    if ismissing(cutoff)
        cutoff = cutoff2X
    end
    
    if ismissing(cutoff2)
        cutoff2=cutoff3bX
    end
    threebody=true
    if cutoff2 < 1e-5
        threebody= false
    end

    dmin_types = Dict()
    dmin_types3 = Dict()
    for t1 = crys.stypes
        for t2 in crys.stypes
            dmin_types[Set((t1,t2))] = get_cutoff(t1,t2)[1] * shrink
            for t3 in crys.stypes
                dmin_types3[Set((t1,t2,t3))] = get_cutoff(t1,t2,t3) * shrink
            end
        end
    end

    

    R = get_grid(crys, 35.0)
    nr = (R[1]*2+1)*(R[2]*2+1)*(R[3]*2+1)

    #    println("R grid ", R, " " , nr)
    
    dist_arr = zeros(var_type, crys.nat,crys.nat,nr,4)
    
    R_f = zeros(var_type, 3, nthreads())
    
    lmn = zeros(var_type, 3)

    coords_ab = zeros(var_type, 3, crys.nat, crys.nat)
    for a = 1:crys.nat
        for b = 1:crys.nat
            coords_ab[:,a,b] = (crys.coords[a,:] .- crys.coords[b,:])' * crys.A 
        end
    end
    
    #    c=0
    At = crys.A'

    S3 = (2*R[3]+1)
    S23 = (2*R[3]+1)*(2*R[2]+1)

    c_zero = 0

    Rind = zeros(Int64, nr, 3)

    found_arr = zeros(Bool, nr)
    found_arr[:] .= false

    #threads
    for c = 1: (R[1]*2+1) * (R[2]*2+1) * (R[3]*2+1)
        
        r3 = mod(c-1 , R[3]*2+1 ) - R[3]
        r2 = mod((c-1) ÷ (R[3]*2+1), (R[2]*2+1)) - R[2]
        r1 = (c-1) ÷ ( (R[2]*2+1)*(R[3]*2+1)) - R[1]

        rf = [r1,r2,r3]
        
        Rind[c,:] .= rf

        RF = At*rf

        found = false
        for a = 1:crys.nat
            ta = crys.stypes[a]
            for b = 1:crys.nat
                tb = crys.stypes[b]            
                #                        println(R_f2)
                #                        println(coords_ab[:,a,b])
                dR = (@view coords_ab[:,a,b]) + RF
                dist = (dR'*dR)^0.5

                dist_arr[a,b,c,1] = dist
                
                if dist > 1e-7
                    dist_arr[a,b,c,2:4].= dR/(dist )
                end
                cutoffXX = get_cutoff(ta,tb)[1] * shrink
                if dist < cutoffXX
                    found = true
                end
                
                if r1 == 0 && r2 == 0 && r3 == 0
                    c_zero = c
                end
                if found
                    found_arr[c] = true
                end
            end
        end
    end

#    println("sum found_arr ", sum(found_arr))
    
    #    R_reverse = Dict()
    #    for key in 1:size(Rind)[1]
    ##        println(key , " ", Rind[key,:])
    #        R_reverse[Rind[key,:]] = key
    #    end

    begin
        Rdiff = zeros(UInt16, size(Rind)[1], size(Rind)[1])
        #Rdiff = zeros(UInt8, size(Rind)[1], size(Rind)[1])
        #Rdiff = spzeros(Int32, size(Rind)[1], size(Rind)[1])

        NR1 = (R[1]*2+1)
        NR2 = (R[2]*2+1)
        NR3 = (R[3]*2+1)

        rr1 = zeros(size(Rind)[1])
        rr2 = zeros(size(Rind)[1])
        rr3 = zeros(size(Rind)[1])
    end
    
    #        test = zeros(Bool, size(Rind[1]))

    for c1 = 1:size(Rind)[1]
        if found_arr[c1]

            rr1[:] .=  -(@view Rind[:, 1]) .+ Rind[c1, 1]
            rr2[:] .=  -(@view Rind[:, 2]) .+ Rind[c1, 2]
            rr3[:] .=  -(@view Rind[:, 3]) .+ Rind[c1, 3]
            
            
            #                test = (rr1 .<= R[1]) .& (rr1 .>= -R[1]) .& (rr2 .<= R[2]) .& (rr2 .>= -R[2]) .& (rr3 .<= R[3]) .& (rr3 .>= -R[3])
            #                Rdiff[c1,:] = test .* ( (rr1 .+ R[1])*NR3*NR2    .+  (rr2.+R[2])*NR3  .+    rr3 .+ 1 .+ R[3])

            #                println(size(test))

            for c2 = 1:size(Rind)[1]
                if found_arr[c2]
                    if rr1[c2] <= R[1] && rr1[c2] >= -R[1] && rr2[c2] <= R[2] && rr2[c2] >= -R[2] && rr3[c2] <= R[3] && rr3[c2] >= -R[3]
                        
                        Rdiff[c1,c2] =  (rr1[c2]+R[1])*NR3*NR2    +  (rr2[c2]+R[2])*NR3  +    rr3[c2] + 1 + R[3]
                        
                    end
                    
                    
                end
            end
        end
    end
    
    

    
    
    #            continue
    #            for c2 = 1:size(Rind)[1]
    #                if Rdiff in keys(R_reverse)
    #                    c12 = R_reverse[ [r1[c2],r2[c2],r3[c2]] ]
    #                    Rdiff[c1,c2] = c12
    #                end
    #            end
    #        end
    
    #    end

    begin
        
        R_keep = zeros(Int64, 0, 4)
        R_dict = Dict()
        fcount =0
        for i in 1:nr
            if found_arr[i]
                c = Rind[i,:]
                R_keep = [R_keep; [0 c']]
                fcount += 1
                R_dict[c] = fcount
                if i == c_zero
                    c_zero = fcount
                end
            end
        end

        R_keep_ab = zeros(Int64, crys.nat*crys.nat*nr, 7)
        
        ind_cutoff = Dict()
        ind_cutoff3bX = Dict()

        keep_counter = 0

    end

    for a = 1:crys.nat
        ta = crys.stypes[a]
        for b = 1:crys.nat
            tb = crys.stypes[b]            
            ind = dist_arr[a,b,:,1] .> 1e-7
            
            dmin = minimum( dist_arr[a,b,ind,1])
            if dmin < dmin_types[Set((ta,tb))]
                dmin_types[Set((ta,tb))] = dmin
            end

            cutoffYY = get_cutoff(ta,tb)[1] * shrink
            ind2 = findall(dist_arr[a,b,:,1] .< cutoffYY)
            ind_cutoff[(a,b)] = deepcopy(ind2)
            for i in ind2
                keep_counter += 1
                r = Rind[i,:]
                ikeep = R_dict[r]
                
                R_keep_ab[keep_counter,:] = [i, a, b, r[1], r[2], r[3], ikeep ]
            end
            
            ind3 = findall(dist_arr[a,b,:,1] .< cutoff3bX)
            ind_cutoff3bX[(a,b)] = deepcopy(ind3)

        end
    end
    
    begin
        
        R_keep_ab = R_keep_ab[1:keep_counter,:]
        
        ############

        MEMCHUNK = min(nr*nr * crys.nat^3, 2000)
        if !threebody || return_floats
            MEMCHUNK = 10
        end
        TOTMEM = MEMCHUNK * ones(Int32, nthreads())

        AI3 = []
        AF3 = []
        COUNTER = zeros(Int32, nthreads())
        for i = 1:nthreads()
            array_ind3 = zeros(Int32, MEMCHUNK, 5)
            #        array_floats3 = zeros(var_type, MEMCHUNK , 14)
            array_floats3 = zeros(var_type, MEMCHUNK , 11)
            
            push!(AI3, array_ind3)
            push!(AF3, array_floats3)
        end

        
        

        dmat = dist_arr[:,:,:,2:4] .* dist_arr[:,:,:, 1]

    end
    

    

    if threebody
        #        for a = 1:crys.nat        
        #                       ta = crys.stypes[a]
        #id = 1
        
        #             for b = 1:crys.nat

        #threads
        @threads for ab = 1:crys.nat^2
            b = mod(ab-1, crys.nat)+1
            a = (ab-1) ÷ crys.nat   +1
            
            ta = crys.stypes[a]

            tb = crys.stypes[b]
            cutoffZZ = get_cutoff(ta,tb)[1] * shrink

            id = threadid()
            #id = 1
            array_ind3X = AI3[id]
            array_floats3X = AF3[id]

            
            for c1 = ind_cutoff[(a,b)]

                dist_ab = dist_arr[a,b, c1 , 1]
                lmn_ab = dist_arr[a,b, c1  , 2:4]
                
                d_ab = dmat[a,b, c1,:]
                
                cut_ab = cutoff_fn(dist_ab, cutoffZZ - cutoff_length, cutoffZZ)
                

                
                for c = 1:crys.nat
                    
                    tc = crys.stypes[c]
                    cutoff3 = get_cutoff(ta,tb,tc) * shrink
                    
                    cut_ab2 = cutoff_fn(dist_ab, cutoff3 - cutoff_length, cutoff3)
                    
                    for c2 in ind_cutoff[(a,c)]
                        
                        dist_ac = dist_arr[a,c, c2, 1]
                        
                        if dist_ac > cutoff3
                            continue
                        end
                        
                        c12 = Rdiff[c1,c2]
                        
                        
                        if c12 == 0
                            continue
                        end
                        

                        
                        #                            
                        #                            
                        #                            if rdiff in keys(R_reverse)
                        #                                c12 = R_reverse[ rdiff ]
                        #                            else
                        #                                continue
                        #                            end
                        
                        dist_bc = dist_arr[c, b, c12, 1]
                        
                        if ( exp(-aX2*dist_ac)*exp(-aX2*dist_bc)*1000 < 0.2e-6 && exp(-aX2*dist_ac)*exp(-aX2*dist_ab)*1000 < 0.2e-6)
                            continue
                        end

                        
                        if dist_bc < cutoff3 && dist_ab > 1e-5 && dist_bc > 1e-5 && dist_ac > 1e-5
                            lmn_bc = - (@view  dist_arr[c,b,c12, 2:4])
                            
                            lmn_ac = (@view dist_arr[a,c, c2, 2:4])
                            #                            lmn_bc = - temp ./ dist_bc
                            
                            
                            
                            d0=dmin_types3[Set((ta,tb,tc))]
                            if (dist_bc + dist_ab + dist_ac)/3.0 < d0
                                dmin_types3[Set((ta,tb,tc))]  = (dist_bc + dist_ab + dist_ac)/3.0
                            end
                            
                            cut_ac = cutoff_fn(dist_ac, cutoff3 - cutoff_length, cutoff3)
                            cut_bc = cutoff_fn(dist_bc, cutoff3 - cutoff_length, cutoff3)

                            
                                COUNTER[id] += 1
                                if COUNTER[id] > TOTMEM[id] #need more memory
                                    TOTMEM[id] += MEMCHUNK
                                    array_ind3X = [array_ind3X;zeros(Int32, MEMCHUNK, 5)]
                                    AI3[id] = array_ind3X
                                    if return_floats
                                        array_floats3X = [array_floats3X; zeros(var_type, MEMCHUNK , 11)]
                                        AF3[id] = array_floats3X
                                    end
                                end

                            r = Rind[c1,:]
                            ikeep = R_dict[r]
#                            r = Rind[c1,:]
#                            ikeep = R_dict[r]
                            r2 = Rind[c2,:]
                            ikeep2 = R_dict[r2]                            
                            
                            
                            array_ind3X[COUNTER[id],:] .= [a,b,c,ikeep,ikeep2]
                            #                            array_floats3X[COUNTER[id], :] .= [dist_ab, dist_ac, dist_bc, lmn_ab[1],lmn_ab[2], lmn_ab[3], lmn_ac[1],lmn_ac[2], lmn_ac[3], lmn_bc[1],lmn_bc[2], lmn_bc[3], cut_ab*cut_bc*cut_ac, cut_ab2*cut_bc*cut_ac]
                            if return_floats
                                array_floats3X[COUNTER[id], :] .= [dist_ab, dist_ac, dist_bc,  lmn_ac[1],lmn_ac[2], lmn_ac[3], lmn_bc[1],lmn_bc[2], lmn_bc[3], cut_ab*cut_bc*cut_ac, cut_ab2*cut_bc*cut_ac]                            
                            end
                            
                            
                            #                            if COUNTER[id] == 761  ||  COUNTER[id] == 762 || COUNTER[id] == 763 || COUNTER[id] == 764
                            #                                println(array_ind3X[COUNTER[id],:])
                            #                            end

                            
                            #                            catch
                            #                                println(id, " " , COUNTER[id], " ", TOTMEM[id])
                            #                            end
                        end
                    end
                end
            end
        end
    end
    #    end



    
    begin
#        array_ind3 = zeros(Int64, 0, 5)

        #        array_ind3 = zeros(Int64, sum(COUNTER), 5)
        array_ind3 = zeros(UInt16, sum(COUNTER), 5)

        array_floats3 = zeros(var_type, 0 , 11)

        cx = 1
        for (counter, i,f) in zip(COUNTER, AI3, AF3)
#            array_ind3 = [array_ind3; i[1:counter, :]]
            array_ind3[cx:cx+counter-1,:] = i[1:counter, :]
            cx += counter

            if return_floats
                array_floats3 = [array_floats3; f[1:counter,:]]
            end
        end
    end
    
    return R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3, Rind
    #    return R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3

end


function distances_etc_3bdy_parallel2(crys, cutoff=missing, cutoff2=missing; var_type=Float64)

#    println("cutoff $cutoff $cutoff2")
    
    if ismissing(cutoff)
        cutoff = cutoff2X
    end
    
    if ismissing(cutoff2)
        cutoff2=cutoff3bX
    end
    threebody=true
    if cutoff2 < 1e-5
        threebody= false
    end

    dmin_types = Dict()
    dmin_types3 = Dict()
    for t1 = crys.stypes
        for t2 in crys.stypes
            dmin_types[Set((t1,t2))] = get_cutoff(t1,t2)[1]
            for t3 in crys.stypes
                dmin_types3[Set((t1,t2,t3))] = get_cutoff(t1,t2,t3)
            end
        end
    end

    

    R = get_grid(crys, 35.0)
    nr = (R[1]*2+1)*(R[2]*2+1)*(R[3]*2+1)

#    println("R grid ", R, " " , nr)
    
    dist_arr = zeros(var_type, crys.nat,crys.nat,nr,4)
                     
    R_f = zeros(var_type, 3, nthreads())
    
    lmn = zeros(var_type, 3)

    coords_ab = zeros(var_type, 3, crys.nat, crys.nat)
    for a = 1:crys.nat
        for b = 1:crys.nat
            coords_ab[:,a,b] = (crys.coords[a,:] .- crys.coords[b,:])' * crys.A 
        end
    end
    
#    c=0
    At = crys.A'

    S3 = (2*R[3]+1)
    S23 = (2*R[3]+1)*(2*R[2]+1)

    c_zero = 0

    Rind = zeros(Int64, nr, 3)

    found_arr = zeros(Bool, nr)
    found_arr[:] .= false


    @threads for c = 1: (R[1]*2+1) * (R[2]*2+1) * (R[3]*2+1)
        
        r3 = mod(c-1 , R[3]*2+1 ) - R[3]
        r2 = mod((c-1) ÷ (R[3]*2+1), (R[2]*2+1)) - R[2]
        r1 = (c-1) ÷ ( (R[2]*2+1)*(R[3]*2+1)) - R[1]

        rf = [r1,r2,r3]
        
        Rind[c,:] .= rf

        RF = At*rf

        found = false
        for a = 1:crys.nat
            ta = crys.stypes[a]
            for b = 1:crys.nat
                tb = crys.stypes[b]            
                #                        println(R_f2)
                #                        println(coords_ab[:,a,b])
                dR = (@view coords_ab[:,a,b]) + RF
                dist = (dR'*dR)^0.5

                dist_arr[a,b,c,1] = dist
                        
                if dist > 1e-7
                    dist_arr[a,b,c,2:4].= dR/(dist )
                end
                cutoffXX = get_cutoff(ta,tb)[1]
                if dist < cutoffXX
                    found = true
                end
                        
                if r1 == 0 && r2 == 0 && r3 == 0
                    c_zero = c
                end
                if found
                    found_arr[c] = true
                end
            end
        end
    end

    
#    R_reverse = Dict()
#    for key in 1:size(Rind)[1]
##        println(key , " ", Rind[key,:])
#        R_reverse[Rind[key,:]] = key
#    end

    begin
            Rdiff = zeros(UInt16, size(Rind)[1], size(Rind)[1])
            #Rdiff = spzeros(Int32, size(Rind)[1], size(Rind)[1])

        NR1 = (R[1]*2+1)
        NR2 = (R[2]*2+1)
        NR3 = (R[3]*2+1)

        rr1 = zeros(size(Rind)[1])
        rr2 = zeros(size(Rind)[1])
        rr3 = zeros(size(Rind)[1])
    end
        
#        test = zeros(Bool, size(Rind[1]))

    for c1 = 1:size(Rind)[1]
        if found_arr[c1]

            rr1[:] .=  -(@view Rind[:, 1]) .+ Rind[c1, 1]
            rr2[:] .=  -(@view Rind[:, 2]) .+ Rind[c1, 2]
            rr3[:] .=  -(@view Rind[:, 3]) .+ Rind[c1, 3]
                    
                    
                    #                test = (rr1 .<= R[1]) .& (rr1 .>= -R[1]) .& (rr2 .<= R[2]) .& (rr2 .>= -R[2]) .& (rr3 .<= R[3]) .& (rr3 .>= -R[3])
                    #                Rdiff[c1,:] = test .* ( (rr1 .+ R[1])*NR3*NR2    .+  (rr2.+R[2])*NR3  .+    rr3 .+ 1 .+ R[3])

                    #                println(size(test))

            for c2 = 1:size(Rind)[1]
                if found_arr[c2]
                    if rr1[c2] <= R[1] && rr1[c2] >= -R[1] && rr2[c2] <= R[2] && rr2[c2] >= -R[2] && rr3[c2] <= R[3] && rr3[c2] >= -R[3]
                        
                        Rdiff[c1,c2] =  (rr1[c2]+R[1])*NR3*NR2    +  (rr2[c2]+R[2])*NR3  +    rr3[c2] + 1 + R[3]
                        
                    end
                    
                    
                end
            end
        end
    end
    
    

                        
                
                #            continue
#            for c2 = 1:size(Rind)[1]
#                if Rdiff in keys(R_reverse)
#                    c12 = R_reverse[ [r1[c2],r2[c2],r3[c2]] ]
#                    Rdiff[c1,c2] = c12
#                end
#            end
#        end
        
#    end


    
    R_keep = zeros(Int64, 0, 4)
    R_dict = Dict()
    fcount =0
    for i in 1:nr
        if found_arr[i]
            c = Rind[i,:]
            R_keep = [R_keep; [0 c']]
            fcount += 1
            R_dict[c] = fcount
            if i == c_zero
                c_zero = fcount
            end
        end
    end

    R_keep_ab = zeros(Int64, crys.nat*crys.nat*nr, 7)
    
    ind_cutoff = Dict()
    ind_cutoff3bX = Dict()

    keep_counter = 0
    
    for a = 1:crys.nat
        ta = crys.stypes[a]
        for b = 1:crys.nat
            tb = crys.stypes[b]            
            ind = dist_arr[a,b,:,1] .> 1e-7
            
            dmin = minimum( dist_arr[a,b,ind,1])
            if dmin < dmin_types[Set((ta,tb))]
                dmin_types[Set((ta,tb))] = dmin
            end

            cutoffYY = get_cutoff(ta,tb)[1]
            ind2 = findall(dist_arr[a,b,:,1] .< cutoffYY)
            ind_cutoff[(a,b)] = deepcopy(ind2)
            for i in ind2
                keep_counter += 1
                r = Rind[i,:]
                ikeep = R_dict[r]
                
                R_keep_ab[keep_counter,:] = [i, a, b, r[1], r[2], r[3], ikeep ]
            end
            
            ind3 = findall(dist_arr[a,b,:,1] .< cutoff3bX)
            ind_cutoff3bX[(a,b)] = deepcopy(ind3)

        end
    end

    
    R_keep_ab = R_keep_ab[1:keep_counter,:]
    
    ############

    MEMCHUNK = min(nr*nr * crys.nat^3, 20000)
    TOTMEM = MEMCHUNK * ones(Int64, nthreads())

    AI3 = []
    AF3 = []
    COUNTER = zeros(Int64, nthreads())
    for i = 1:nthreads()
        array_ind3 = zeros(Int64, MEMCHUNK, 5)
        array_floats3 = zeros(var_type, MEMCHUNK , 14)
        push!(AI3, array_ind3)
        push!(AF3, array_floats3)
    end

    
        

    dmat = dist_arr[:,:,:,2:4] .* dist_arr[:,:,:, 1]


    

    
#    println("threebody")
    if threebody
        #        for a = 1:crys.nat        
        #                       ta = crys.stypes[a]
        #id = 1
        
        #             for b = 1:crys.nat

        @threads for ab = 1:crys.nat^2
            b = mod(ab-1, crys.nat)+1
            a = (ab-1) ÷ crys.nat   +1
            
            ta = crys.stypes[a]

            tb = crys.stypes[b]
            cutoffZZ = get_cutoff(ta,tb)[1]

            id = threadid()
            #id = 1
            array_ind3X = AI3[id]
            array_floats3X = AF3[id]

            
            for c1 = ind_cutoff[(a,b)]

                dist_ab = dist_arr[a,b, c1 , 1]
                lmn_ab = dist_arr[a,b, c1  , 2:4]
                
                d_ab = dmat[a,b, c1,:]
                
                cut_ab = cutoff_fn(dist_ab, cutoffZZ - cutoff_length, cutoffZZ)
                

                
                for c = 1:crys.nat
                    
                    tc = crys.stypes[c]
                    cutoff3 = get_cutoff(ta,tb,tc)
                    
                    cut_ab2 = cutoff_fn(dist_ab, cutoff3 - cutoff_length, cutoff3)
                    
                    for c2 in ind_cutoff[(a,c)]
                        
                        dist_ac = dist_arr[a,c, c2, 1]
                        
                        if dist_ac > cutoff3
                            continue
                        end
                        
                        c12 = Rdiff[c1,c2]
                        
                        
                        if c12 == 0
                            continue
                        end
                        

                        
                        #                            
                        #                            
                        #                            if rdiff in keys(R_reverse)
                        #                                c12 = R_reverse[ rdiff ]
                        #                            else
                        #                                continue
                        #                            end
                        
                        dist_bc = dist_arr[c, b, c12, 1]
                        

                        
                        if dist_bc < cutoff3 && dist_ab > 1e-5 && dist_bc > 1e-5 && dist_ac > 1e-5
                            lmn_bc = - (@view  dist_arr[c,b,c12, 2:4])
                            
                            lmn_ac = (@view dist_arr[a,c, c2, 2:4])
                            #                            lmn_bc = - temp ./ dist_bc
                            
                            
                            
                            d0=dmin_types3[Set((ta,tb,tc))]
                            if (dist_bc + dist_ab + dist_ac)/3.0 < d0
                                dmin_types3[Set((ta,tb,tc))]  = (dist_bc + dist_ab + dist_ac)/3.0
                            end
                            
                            cut_ac = cutoff_fn(dist_ac, cutoff3 - cutoff_length, cutoff3)
                            cut_bc = cutoff_fn(dist_bc, cutoff3 - cutoff_length, cutoff3)
                            
                            COUNTER[id] += 1
                            if COUNTER[id] > TOTMEM[id] #need more memory
                                TOTMEM[id] += MEMCHUNK
                                array_ind3X = [array_ind3X;zeros(Int64, MEMCHUNK, 5)]
                                AI3[id] = array_ind3X
                                array_floats3X = [array_floats3X; zeros(var_type, MEMCHUNK , 14)]
                                AF3[id] = array_floats3X
                                
                            end

                            r = Rind[c1,:]
                            ikeep = R_dict[r]
                            
                            
                            array_ind3X[COUNTER[id],:] .= [a,b,c,ikeep,c2]
                            array_floats3X[COUNTER[id], :] .= [dist_ab, dist_ac, dist_bc, lmn_ab[1],lmn_ab[2], lmn_ab[3], lmn_ac[1],lmn_ac[2], lmn_ac[3], lmn_bc[1],lmn_bc[2], lmn_bc[3], cut_ab*cut_bc*cut_ac, cut_ab2*cut_bc*cut_ac]

                            
#                            if COUNTER[id] == 761  ||  COUNTER[id] == 762 || COUNTER[id] == 763 || COUNTER[id] == 764
#                                println(array_ind3X[COUNTER[id],:])
#                            end

                            
                            #                            catch
                            #                                println(id, " " , COUNTER[id], " ", TOTMEM[id])
                            #                            end
                        end
                    end
                end
            end
        end
    end
#    end



                            
                
    array_ind3 = zeros(Int64, 0, 5)
    array_floats3 = zeros(var_type, 0 , 14)
    for (counter, i,f) in zip(COUNTER, AI3, AF3)
        array_ind3 = [array_ind3; i[1:counter, :]]
        array_floats3 = [array_floats3; f[1:counter,:]]
    end
    
    return R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3, Rind
#    return R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3

end

#function get_subgrid(f, subgrid)
#    for i = 1:3
#        for g = 1:subgrid[i]
#            if crys.coords[a,i] <= 1/g
#                subgrid[a,i] = g
#                break
#            end
#        end
#    end

#end

function distances_subgrid(crys; cutoff=missing, threebody=true, keep_extra= false)

    
    if ismissing(cutoff)
        cutoff = cutoff2X
    end

    cutoff_more = cutoff * 1.000001

    subgrid = [1,1,1]
    #adjust subgrid

    #adjust subgrid so that there are at least several atoms in the subgrids, if possible.
    for factor = 1.0:0.2:20.0
    
        subgrid1 = Int64(floor(1/factor * sqrt(sum(crys.A[1,:].^2)) / cutoff_more) + 1)
        subgrid2 = Int64(floor(1/factor * sqrt(sum(crys.A[2,:].^2)) / cutoff_more) + 1)
        subgrid3 = Int64(floor(1/factor * sqrt(sum(crys.A[3,:].^2)) / cutoff_more) + 1)
        
        subgrid = [subgrid1, subgrid2, subgrid3]
        
        #println("gridlength $gridlength")
        #look for fake dimensions (all atoms in same subgrid)
        for i = 1:3
            if sum(abs.(crys.coords[:,i] .- crys.coords[1,i])) < 1e-7
                subgrid[i] = min(2, subgrid[i])
            end
        end

        if crys.nat / prod(subgrid) > 10
            break
        end
#        break
        
    end


    
#    println("subgrid $subgrid")
    gridlength = [sqrt(sum(crys.A[1,:].^2 )), sqrt(sum(crys.A[2,:].^2 )),sqrt(sum(crys.A[3,:].^2 ))] ./ subgrid
    gridlength = maximum(gridlength)

#    return subgrid
    
    subgrid_dict = Dict()
    for g1 = 1:subgrid[1]
        for g2 = 1:subgrid[2]
            for g3 = 1:subgrid[3]
                subgrid_dict[g1,g2,g3] = Int64[]
            end
        end
    end
    
    crys.coords = mod.(crys.coords, 1.0)
    subgrid_atoms = zeros(Int64, crys.nat, 3)
    for a = 1:crys.nat
#        println("a $a")
        for i = 1:3
            for g = 1:subgrid[i]
#                println("$i $a compare $(crys.coords[a,i] + 1e-10)  $(g / subgrid[i])")
                if crys.coords[a,i] + 1e-10 <= g / subgrid[i]
                    subgrid_atoms[a,i] = g
                    break
                end
            end
        end
        
        push!(subgrid_dict[subgrid_atoms[a,1], subgrid_atoms[a,2],subgrid_atoms[a,3]], a)
    end
    maxgrid = 0
    for k in keys(subgrid_dict)
        maxgrid = max(maxgrid, length(subgrid_dict[k]))
    end

    
    
    R = get_grid(crys, 35.0)

    
    nr = (R[1]*2+1)*(R[2]*2+1)*(R[3]*2+1)
    grid_coords = zeros(prod(subgrid), 3)
    c = 0
    subgrid_num = zeros(Int64,prod(subgrid), 3)
    for g1 = 1:subgrid[1]
        for g2 = 1:subgrid[2]
            for g3 = 1:subgrid[3]
                c += 1
                grid_coords[c,:] = [ (g1-1) / subgrid[1] , (g2-1) / subgrid[2], (g3-1) / subgrid[3]]  + 0.5 * 1 ./ subgrid
               
                subgrid_num[c,:] = [g1,g2,g3]
            end
        end
    end
    #println("grid_coords")
    #println(grid_coords)
#    println("subgrid_num")
#    println(subgrid_num)
    
    fake_crys = makecrys(crys.A, grid_coords, repeat([:H], prod(subgrid)))
    #println("fake_crys")
    #println(fake_crys)
    #println("fake")
    nz_inds, R_keep_ab_TT, nz_ind3, dist3_nonzero, dist_arr_TT, c_zero_tt, dmin_types_TT, dmin_types3_TT = distances_etc_3bdy_parallel_LV_old(fake_crys, max(gridlength*2.001, 2.001*cutoff ), -0.1; var_type=Float64, return_floats=false, keep_extra=true)

#    println("R_keep_ab_TT")
#    println(R_keep_ab_TT)
    
#    R_keep_ab_TT_big = zeros(Int64, 0,7)
    counter = 0
    At = collect((crys.A)')
    #println("final")
    R_keep_ab_TT_big = zeros(Int64, 100000, 7)
    dist2_list = zeros(Float64, 100000)
    dist = 0.0

    crys_cart =zeros(  crys.nat,3)
    crys_cart[:,:] =    collect(crys.coords * crys.A)

    nz_inds = Set()
    dmin_types_TT = Dict()
    dmin_types_TT = Dict()
    s = Set(crys.stypes)
    for s1 in s
        for s2 in s
            ss = Set([s1,s2])
            dmin_types_TT[ss] = 1000000000000.0
        end
    end

    Ar  = zeros(3)
    Ar_t  = zeros(3)
    for c_ab = 1:size(R_keep_ab_TT)[1]
        (n,sg1,sg2, r1,r2,r3,q) = R_keep_ab_TT[c_ab,:]
#        println("sg1 $sg1 sg2 $sg2 [$r1 $r2 $r3]")

        Ar[:] = At * [r1,r2,r3]
        atoms1 = subgrid_dict[subgrid_num[sg1,1], subgrid_num[sg1,2], subgrid_num[sg1,3]]
        atoms2 = subgrid_dict[subgrid_num[sg2,1], subgrid_num[sg2,2], subgrid_num[sg2,3]]
#        println("atoms1 $atoms1 atoms2 $atoms2")
        @fastmath @inbounds @simd for a1 in atoms1
            ta = crys.stypes[a1]
            Ar_t .= Ar
            for i = 1:3
                Ar_t[i] += crys_cart[a1,i]
            end
            for a2 in atoms2
                #sum((At * (crys.coords[a1,:] - crys.coords[a2,:] + [r1,r2,r3] ) ).^2)
                dist = 0.0
                for i = 1:3
                    dist += (Ar_t[i] - crys_cart[a2,i] )^2
                end
                tb = crys.stypes[a2]
               
                

                if (keep_extra && dist < cutoff^2) || (!keep_extra && dist <  get_cutoff(ta,tb)[1]^2)

                    counter += 1
                
#                    println("x $counter   a1 $a1 a2 $a2 dist $(sqrt(dist)) r1 r2 r3 $r1 $r2 $r3")

                    #                    println("add $([counter a1 a2 r1 r2 r3 q]) !!!!!!!!!!!!!!!!!!!!!!!!!")
                    #R_keep_ab_TT_big = vcat(R_keep_ab_TT_big,[counter a1 a2 r1 r2 r3 q])
                    R_keep_ab_TT_big[counter,:] = [counter, a1 ,a2 ,r1, r2, r3, q]
                    dist2_list[counter] = dist
                    push!(nz_inds, [r1,r2,r3])
                    if mod(counter, 100000) == 0
                        R_keep_ab_TT_big = [R_keep_ab_TT_big;zeros(Int64, 100000, 7)]
                        dist2_list = [dist2_list ;zeros(Float64, 100000)]
                    end
                    if dist > 1e-5
                        dmin_types_TT[Set([crys.stypes[a1], crys.stypes[a2]])] = min(dmin_types_TT[Set([crys.stypes[a1], crys.stypes[a2]])], sqrt(dist))
                    end
                        
                end

            end
        end

    end
    R_keep_ab_TT_big = R_keep_ab_TT_big[1:counter,:]
    nz_inds = sort(collect(nz_inds))
    nz_inds2 = zeros(Int64, length(nz_inds), 4)
    c_zero = -1
    R_dict = Dict()
    for (c, nz) in enumerate(nz_inds)
#        println("c $c nz $nz")
        nz_inds2[c,2:4] = nz
        R_dict[nz[1], nz[2], nz[3]]= c
        
        if nz[1] ==0 && nz[2] ==0 && nz[3] ==0
#            println("c_zero $c")
            c_zero = c
        end
    end



    

    s = Set(crys.stypes)
    dmin_types_TT = Dict()
    for s1 in s
        for s2 in s
            ss = Set([s1,s2])
            dmin_types_TT[ss] = 1000000000000.0
            for c in 1:size(dist_arr_TT,1)
                a = R_keep_ab_TT[c,2]
                b = R_keep_ab_TT[c,3]
                if crys.stypes[a] == s1 && crys.stypes[b] == s2 && dist_arr_TT[c,1] > 1e-5
                    dmin_types_TT[ss] = min(dmin_types_TT[ss], dist_arr_TT[c,1])
                end
            end
        end
    end

    for nz =  1:size(R_keep_ab_TT_big)[1]
        R_keep_ab_TT_big[nz,7] = R_dict[R_keep_ab_TT_big[nz,4], R_keep_ab_TT_big[nz,5],R_keep_ab_TT_big[nz,6]]
    end
    
    if threebody
        #println("threebody")
        nz_inds3 = threebody_dist_subgrid(crys, R_keep_ab_TT_big, R_dict, dist2_list, c_zero)
    else
        nz_inds3 = missing
    end
    return nz_inds2, R_keep_ab_TT_big, nz_inds3, missing, sqrt.(dist2_list), c_zero, dmin_types_TT, missing, missing
    
end

function threebody_dist_subgrid(crys, R_keep_ab, R_dict, dist2_list, c_zero)
#    println("inside 3")
#    println("c_zero $c_zero")
    nab = size(R_keep_ab)[1]

    inds_by_atom= Dict()
    for a = 1:crys.nat
        inds_by_atom[a] = Int64[]
    end
    
    for ind1 = 1:nab
#        a1 = R_keep_ab[ind1,2]
        #        a2 = R_keep_ab[ind1,3]
        c_temp, a1 ,a2 ,r1, r2, r3, q = R_keep_ab[ind1,:]            
        if a1 == a2 && R_dict[r1,r2,r3] == c_zero
            continue
        end
        
        push!(inds_by_atom[a1], ind1)
    end

    twobody_dict = Dict{Tuple{Int64,Int64, Int64, Int64, Int64}, Float64}()
    for ind1 = 1:nab
        counter, a1 ,a2 ,r1, r2, r3, q = R_keep_ab[ind1,:]            
        twobody_dict[a1,a2,r1,r2,r3] = dist2_list[ind1]
    end


    
    nz_ind3 = zeros(Int64, 100000, 5)
#    println("main")
    counter = 0
    for a1 = 1:crys.nat
        ta1 = crys.stypes[a1]
        for ind12 = inds_by_atom[a1]
            c_temp, a1a ,a2 ,r1, r2, r3, q = R_keep_ab[ind12,:]            

#            if twobody_dict[a1,a2,r1,r2,r3] > get_cutoff(ta1, crys.stypes[a2] )[1]^2
#                continue
#            end
                
#            if a1 == a2 && R_dict[r1,r2,r3] == c_zero
#                println("continue")
#                continue
#            end

            for ind13 = inds_by_atom[a1]
                if ind12 == ind13
                    continue
                end
                c_temp_3, a1b ,a3 ,r1_3, r2_3, r3_3, q_3 = R_keep_ab[ind13,:]
#                if a1 == a3 && R_dict[r1_3,r2_3,r3_3] == c_zero
#                    continue
#                end

                #                if (a2,a3,r1_3 - r1,r2_3 - r2,r3_3 - r3) in keys(twobody_dict)
                if haskey(twobody_dict, (a2,a3,r1_3 - r1,r2_3 - r2,r3_3 - r3))
#                    x = true
#                end
                #if false

#                    println(" twobody_dict[a1,a3,r1_3,r2_3,r3_3] ", twobody_dict[a1,a3,r1_3,r2_3,r3_3])
#                    println(" get_cutoff[ta1, crys.stypes[a3]][2]^2 ",  get_cutoff(ta1, crys.stypes[a3])[2]^2 )
                    d13 = twobody_dict[a1,a3,r1_3,r2_3,r3_3]
                    if d13 > get_cutoff(ta1, crys.stypes[a2], crys.stypes[a3])[1]^2
                        continue
                    end
                    d23 = twobody_dict[a2,a3,r1_3 - r1,r2_3 - r2,r3_3-r3]
                    if d23 > get_cutoff(ta1, crys.stypes[a2], crys.stypes[a3])[1]^2
                        continue
                    end
                    
                    if  exp(-aX2*sqrt(d13))*exp(-aX2*sqrt(d23))*1000 < 0.2e-6
                        continue
                    end
                    
                    counter += 1
#                    println("add $counter ; $([a1, a2, a3, R_dict[r1,r2,r3], R_dict[r1_3,r2_3,r3_3]])")
                    nz_ind3[counter,1] = a1
                    nz_ind3[counter,2] = a2
                    nz_ind3[counter,3] = a3
                    nz_ind3[counter,4] = R_dict[r1,r2,r3]
                    nz_ind3[counter,5] = R_dict[r1_3,r2_3,r3_3]

                    if mod(counter, 100000) == 0
                        nz_ind3 = [nz_ind3; zeros(Int64, 100000, 5)]
                    end
                end
#                end
            end
        end
    end
                
        #counter, a1 ,a2 ,r1, r2, r3, q = R_keep_ab[ind1,:]
#    println("done")
    return nz_ind3[1:counter,:]

end


function distances_etc_3bdy_parallel_LV_new(crys, cutoff=missing, cutoff2=missing; var_type=Float64, return_floats=true, shrink = 1.0, R=missing, cutoff4 = -1.0, keep_extra = false)

    if crys.nat < 10 || cutoff4 > 1e-5
        return distances_etc_3bdy_parallel_LV_old(crys, cutoff, cutoff2; var_type=var_type, return_floats=return_floats, shrink = shrink, R=R, cutoff4 = cutoff4, keep_extra = keep_extra)
    else
        if ismissing(cutoff2)
            cutoff3 = true
        elseif cutoff2 > 1e-5
            cutoff3 = true
        else
            cutoff3 = false
        end
        println("distances_subgrid ", [cutoff, cutoff3, keep_extra])
        return distances_subgrid(crys; cutoff=cutoff, threebody=cutoff3, keep_extra= keep_extra)
    end


end

function distances_etc_3bdy_parallel_LV(crys, cutoff=missing, cutoff2=missing; var_type=Float64, return_floats=true, shrink = 1.0, R=missing, cutoff4 = -1.0, keep_extra = false)
    #        println("cutoff $cutoff $cutoff2")

    begin
        
        if ismissing(cutoff)
            cutoff = cutoff2X
        end
        
        if ismissing(cutoff2)
            cutoff2=cutoff3bX
        end
        threebody=true
        if cutoff2 < 1e-5
            threebody= false
        end

        fourbody=true
        if cutoff4 < 1e-5
            fourbody=false
        end
        
        if ismissing(R)
            R = get_grid(crys, 35.0)
        end
        
        nr = (R[1]*2+1)*(R[2]*2+1)*(R[3]*2+1)

        dist_arr_TT = zeros(var_type, crys.nat,crys.nat,nr,4)
        
        coords_ab_TT = zeros(var_type, 3 , crys.nat, crys.nat)
        for a = 1:crys.nat
            for b = 1:crys.nat
                coords_ab_TT[:,a,b] = (crys.coords[a,:] .- crys.coords[b,:])' * crys.A 
            end
        end
        
        At = deepcopy(crys.A')

        #Rind_TT = zeros(Int64, nr, 3)

        cutoff_arr = zeros(crys.nat, crys.nat, 2)


        #        get_cutoff_pre = Dict()
        get_cutoff_pre = Dict()       
        s = Set(crys.stypes)
        for s1 in s
            for s2 in s
                get_cutoff_pre[(s1,s2)] = get_cutoff(s1,s2)[:] * shrink
                for s3 in s 
                    get_cutoff_pre[(s1,s2,s3)] = get_cutoff(s1,s2,s3)[1] * shrink
                end
            end
        end

#        println("get_cutoff_pre")
#        println(get_cutoff_pre)
        
        for a = 1:crys.nat
            ta = crys.stypes[a]
            for b = 1:crys.nat
                tb = crys.stypes[b]            
                cutoff_arr[a,b,:] = get_cutoff_pre[(ta,tb)][:]
#                for c = 1:crys.nat
#                    tc = crys.stypes[c]            
#                    cutoff_arr3[a,b,c] =  get_cutoff_pre[(ta,tb,tc)]
#                end
            end
        end

            
        
        
#        for c = 1:nr
#            
#            r3 = mod(c-1 , R[3]*2+1 ) - R[3]
#            r2 = mod((c-1) ÷ (R[3]*2+1), (R[2]*2+1)) - R[2]
#            r1 = (c-1) ÷ ( (R[2]*2+1)*(R[3]*2+1)) - R[1]
#
#            Rind_TT[c,1] = r1
#            Rind_TT[c,2] = r2
#            Rind_TT[c,3] = r3
#        end

        nat = Int64(crys.nat)

       rf1 = Float64.(-R[1]:R[1])
        rf2 = Float64.(-R[2]:R[2])
        rf3 = Float64.(-R[3]:R[3])
        nr1 = (2*R[1]+1)
        nr2 = (2*R[2]+1)
        nr3 = (2*R[3]+1)

        A = crys.A

#        println(" R ", R)
#        println("size coords ", size(coords_ab_TT))
        
        dist_TT = zeros(var_type, nr1,nr2,nr3,nat, nat,4)
        @turbo for r1 = eachindex(rf1) 
            for r2 = eachindex(rf2)
                for r3 = eachindex(rf3)
                    for a = 1:nat
                        for b = 1:nat
                            for i = 1:3
                                dist_TT[r1,r2,r3,a,b,1] +=  (coords_ab_TT[i,a,b] +At[i,1]*rf1[r1] + At[i,2]*rf2[r2]  + At[i,3]*rf3[r3])^2
                                dist_TT[r1,r2,r3,a,b,i+1] +=  (coords_ab_TT[i,a,b] +At[i,1]*rf1[r1] + At[i,2]*rf2[r2]  + At[i,3]*rf3[r3])
                            end
                        end
                    end
                end
            end
        end

        #println("slow")
        @fastmath @inbounds @simd   for r1 = eachindex(rf1) #1
            for r2 = eachindex(rf2)
                for r3 = eachindex(rf3)
                    for a = 1:nat
                        for b = 1:nat
#                            dtemp = dist_TT[r1,r2,r3,a,b,1]^0.5
                            dist_TT[r1,r2,r3,a,b,1] = sqrt(dist_TT[r1,r2,r3,a,b,1])
                        end
                    end
                end
            end
        end
        @turbo  for r1 = eachindex(rf1) #1
            for r2 = eachindex(rf2)
                for r3 = eachindex(rf3)
                    for a = 1:nat
                        for b = 1:nat
                            dist_TT[r1,r2,r3,a,b,2] = dist_TT[r1,r2,r3,a,b,2]/(dist_TT[r1,r2,r3,a,b,1] + 1e-20)
                            dist_TT[r1,r2,r3,a,b,3] = dist_TT[r1,r2,r3,a,b,3]/(dist_TT[r1,r2,r3,a,b,1] + 1e-20)
                            dist_TT[r1,r2,r3,a,b,4] = dist_TT[r1,r2,r3,a,b,4]/(dist_TT[r1,r2,r3,a,b,1] + 1e-20)
                        end
                    end
                end
            end
        end

        found_arr_TT = zeros(Bool, nr1, nr2, nr3)
        #println("mem")
        found_arr_TT_ab = zeros(Bool, nr1, nr2, nr3,nat,nat)
        #println("decide")
        if keep_extra
            for a = 1:crys.nat
                for b = 1:crys.nat
                    found_arr_TT[:,:,:] = found_arr_TT[:,:,:]  .|     (dist_TT[:,:,:,a,b,1] .< max(cutoff_arr[a,b,1], cutoff) )
                    for r1 = eachindex(rf1)
                        for r2 = eachindex(rf2)
                            for r3 = eachindex(rf3)
                                found_arr_TT_ab[r1,r2,r3,a,b] = dist_TT[r1,r2,r3,a,b,1] < max(cutoff_arr[a,b,1], cutoff)
                            end
                        end
                    end
                end
            end

        else
            
            for a = 1:crys.nat
                for b = 1:crys.nat
                    found_arr_TT[:,:,:] = found_arr_TT[:,:,:]  .|     (dist_TT[:,:,:,a,b,1] .< cutoff_arr[a,b,1])
                    for r1 = eachindex(rf1)
                        for r2 = eachindex(rf2)
                            for r3 = eachindex(rf3)
                                found_arr_TT_ab[r1,r2,r3,a,b] = dist_TT[r1,r2,r3,a,b,1] < cutoff_arr[a,b,1]
                            end
                        end
                    end
                end
            end
        end
#        println("sum found_arr_TT ", sum(found_arr_TT))
        
        nz_ind = zeros(Bool, nr)
        nz_ints = zeros(Int64, nr,3)
        nz_inds = zeros(Int64, nr,4)
        R_dict_tt = Dict()

        c=0
        cfound = 0
        for r1 = eachindex(rf1)
            for r2 = eachindex(rf2)
                for r3 = eachindex(rf3)
                    c += 1

                    nz_ints[c,1] = r1
                    nz_ints[c,2] = r2
                    nz_ints[c,3] = r3
                    nz_inds[c,2] = r1 - R[1] - 1
                    nz_inds[c,3] = r2 - R[2] - 1
                    nz_inds[c,4] = r3 - R[3] - 1
                    nz_ind[c] = found_arr_TT[r1,r2,r3]
                    if found_arr_TT[r1,r2,r3]
                        cfound += 1
                        R_dict_tt[[r1 - R[1] - 1, r2 - R[2] - 1, r3 - R[3] - 1 ]] = cfound
                    end
                end
            end
        end

#        println("R_dict [-1, 2, 0]")
#        println(R_dict_tt[[-1,2,0]])
#        println()
        
        nz_ints = nz_ints[nz_ind,:]
        nz_inds = nz_inds[nz_ind,:]

        #println("nz_inds ")
        #for c = 1:size(nz_inds)[1]
        #    println("$c nz ", nz_inds[c,:])
        #end
        
        c_zero_tt = argmin(sum(abs.(nz_inds), dims=2))[1]


        ctt = sum(found_arr_TT_ab)
        R_keep_ab_TT = zeros(Int64, ctt, 7)
        dist_arr_TT = zeros(var_type, ctt,6)    
        c = 0

        for a = 1:crys.nat
            for b = 1:crys.nat
                for r1 = eachindex(rf1)
                    for r2 = eachindex(rf2)
                        for r3 = eachindex(rf3)
                            if found_arr_TT_ab[r1,r2,r3,a,b]
                                c += 1
                                R_keep_ab_TT[c,1] = c
                                R_keep_ab_TT[c,2] = a
                                R_keep_ab_TT[c,3] = b
                                R_keep_ab_TT[c,4] = r1 - R[1] - 1
                                R_keep_ab_TT[c,5] = r2 - R[2] - 1
                                R_keep_ab_TT[c,6] = r3 - R[3] - 1
                                R_keep_ab_TT[c,7] = R_dict_tt[[r1 - R[1] - 1,r2 - R[2] - 1,r3 - R[3] - 1]]

                                dist_arr_TT[c,1] = dist_TT[r1,r2,r3,a,b,1]
                                dist_arr_TT[c,2] = dist_TT[r1,r2,r3,a,b,2]
                                dist_arr_TT[c,3] = dist_TT[r1,r2,r3,a,b,3]
                                dist_arr_TT[c,4] = dist_TT[r1,r2,r3,a,b,4]
                                dist_arr_TT[c,5] = cutoff_fn_fast(dist_TT[r1,r2,r3,a,b,1], cutoff_arr[a, b, 1] - cutoff_length, cutoff_arr[a, b, 1])
                                dist_arr_TT[c,6] = cutoff_fn_fast(dist_TT[r1,r2,r3,a,b,1], cutoff_arr[a, b, 2] - cutoff_length, cutoff_arr[a, b, 2])
                                
                            end
                        end
                    end
                end
            end
        end


        s = Set(crys.stypes)
        dmin_types_TT = Dict()
        for s1 in s
            for s2 in s
                ss = Set([s1,s2])
                dmin_types_TT[ss] = 1000000000000.0
                for c in 1:size(dist_arr_TT,1)
                    a = R_keep_ab_TT[c,2]
                    b = R_keep_ab_TT[c,3]
                    if crys.stypes[a] == s1 && crys.stypes[b] == s2 && dist_arr_TT[c,1] > 1e-5
                        dmin_types_TT[ss] = min(dmin_types_TT[ss], dist_arr_TT[c,1])
                    end
                end
            end
        end
        
        dmin_types3_TT = Dict()
    end
    #    println("threebody dist")
    
    
    if threebody
        
        nz_ind3, dist3_nonzero, dmin_types3_TT, nz_ind4 = dist_threebody(crys, nz_ints, nat, cutoff_arr, coords_ab_TT, At, dist_TT, rf1, rf2, rf3, get_cutoff_pre,R_dict_tt,fourbody ; cutoff2=cutoff2, var_type=var_type,return_floats=return_floats, shrink = 1.0, R=R, cutoff4 = cutoff4, keep_extra = keep_extra)
        
        return nz_inds, R_keep_ab_TT, nz_ind3, dist3_nonzero, dist_arr_TT, c_zero_tt, dmin_types_TT, dmin_types3_TT, nz_ind4
    else
        return nz_inds, R_keep_ab_TT, missing, missing, dist_arr_TT, c_zero_tt, dmin_types_TT, missing, missing
    end        
    
end


function dist_threebody(crys, nz_ints, nat, cutoff_arr, coords_ab_TT, At, dist_TT, rf1, rf2, rf3, get_cutoff_pre,R_dict_tt,fourbody ; cutoff2=-0.1, var_type=Float64,return_floats=true, shrink = 1.0, R=missing, cutoff4 = -1.0, keep_extra = false)
    
    nz_ab = ones(Int64, nat, size(nz_ints,1)*nat,4)
    max_ind_a = zeros(UInt16, nat)

    cutoff_arr3 = zeros(crys.nat, crys.nat, crys.nat)

    for a = 1:crys.nat
        ta = crys.stypes[a]
        for b = 1:crys.nat
            tb = crys.stypes[b]            
            cutoff_arr[a,b,:] = get_cutoff_pre[(ta,tb)][:]
            for c = 1:crys.nat
                tc = crys.stypes[c]            
                cutoff_arr3[a,b,c] =  get_cutoff_pre[(ta,tb,tc)]
            end
        end
    end
    
    for a = 1:nat
        c_ab = 0
        for i = 1:size(nz_ints,1)
            for b = 1:nat
                if dist_TT[nz_ints[i,1],nz_ints[i,2],nz_ints[i,3],a,b,1] < cutoff_arr[a,b,1]
                    c_ab += 1
                    nz_ab[a, c_ab,1] = nz_ints[i,1]
                    nz_ab[a, c_ab,2] = nz_ints[i,2]
                    nz_ab[a, c_ab,3] = nz_ints[i,3]
                    nz_ab[a, c_ab,4] = b
                    max_ind_a[a] = c_ab
                end
            end
        end
    end             

    max_a = maximum(max_ind_a)


    dist3a = zeros(var_type, nat,max_a, max_a)
    @turbo for a = 1:nat
        for i_b = 1:max_a
            for i_c = 1:max_a
                for i = 1:3
                    temp = -coords_ab_TT[i,nz_ab[a,i_b,4], nz_ab[a,i_c,4]] +At[i,1]*(rf1[nz_ab[a,i_b,1]]-rf1[nz_ab[a,i_c,1]]) + At[i,2]*(rf2[nz_ab[a,i_b,2]]-rf2[nz_ab[a,i_c,2]] )  + At[i,3]*(rf3[nz_ab[a,i_b,3]]-rf3[nz_ab[a,i_c,3]])
                    dist3a[a,i_b,i_c] += temp^2
                end
                dist3a[a,i_b,i_c] = dist3a[a,i_b,i_c]^0.5
                
            end
        end
    end

    counter = 0
    NZ = zeros(UInt16, nat*max_a*max_a, 3)
    @inbounds @fastmath @simd for a = 1:nat
        for i_b = 1:max_a
            d12 = dist_TT[nz_ab[a,i_b,1],nz_ab[a,i_b,2],nz_ab[a,i_b,3],a,nz_ab[a,i_b,4],1]
            if i_b <= max_ind_a[a] && d12 > 1e-3 && d12 < cutoff_arr[a, nz_ab[a,i_b,4],1]
                for i_c = 1:max_a 
                    
                    cut3 = cutoff_arr3[a, nz_ab[a,i_b,4],nz_ab[a,i_c,4]]
                    d13 = dist_TT[nz_ab[a,i_c,1],nz_ab[a,i_c,2],nz_ab[a,i_c,3],a,nz_ab[a,i_c,4],1]
                    d23 = dist3a[a,i_b,i_c]
                    
                    #                        if i_c <= max_ind_a[a] && d12 > 1e-3 && d13 > 1e-3 && d23  > 1e-3 && d12 < cutoff_arr[a, nz_ab[a,i_b,4],1] && d13 < cut3 && d23 < cut3
                    #if i_c <= max_ind_a[a] && d13 > 1e-3 && d23  > 1e-3  && d13 < cut3 && d23 < cut3 && ( exp(-aX2*d13)*exp(-aX2*d23)*1000 > 0.2e-6 || exp(-aX2*d13)*exp(-aX2*d12)*1000 > 0.2e-6)   #fixme
                    if i_c <= max_ind_a[a] && d13 > 1e-3 && d23  > 1e-3  && d13 < cut3 && d23 < cut3 &&  exp(-aX2*d13)*exp(-aX2*d23)*1000 > 0.2e-6 
                        counter += 1
                        NZ[counter,1] = a
                        NZ[counter,2] = i_b
                        NZ[counter,3] = i_c
                    end
                end
            end
        end
    end
    #        println("D3 NZ $counter xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    begin
        dist3_nonzero = zeros(var_type, counter,3+3*3+2)
        nz_ind3 = zeros(UInt16 ,counter,5)
    end
    cx = [0,0,0]
    for nz = 1:counter
        a = NZ[nz,1]
        i_b = NZ[nz,2]
        i_c = NZ[nz,3]
        
        nz_ind3[nz, 1] = a
        nz_ind3[nz, 2] = nz_ab[a,i_b,4]
        nz_ind3[nz, 3] = nz_ab[a,i_c,4]
        cx[1] = nz_ab[a,i_b,1]-R[1]-1
        cx[2] = nz_ab[a,i_b,2]-R[2]-1
        cx[3] = nz_ab[a,i_b,3]-R[3]-1
        nz_ind3[nz, 4] = R_dict_tt[cx]
        cx[1] = nz_ab[a,i_c,1]-R[1]-1
        cx[2] = nz_ab[a,i_c,2]-R[2]-1
        cx[3] = nz_ab[a,i_c,3]-R[3]-1
        nz_ind3[nz, 5] = R_dict_tt[cx]
        
    end

    #println("asdf")
    @turbo for nz = 1:counter
        a = NZ[nz,1]
        i_b = NZ[nz,2]
        i_c = NZ[nz,3]

        dist3_nonzero[nz,1] = dist_TT[nz_ab[a,i_b,1],nz_ab[a,i_b,2],nz_ab[a,i_b,3],a,nz_ab[a,i_b,4],1]
        dist3_nonzero[nz,2] = dist_TT[nz_ab[a,i_c,1],nz_ab[a,i_c,2],nz_ab[a,i_c,3],a,nz_ab[a,i_c,4],1]
        
        dist3_nonzero[nz,4] = dist_TT[nz_ab[a,i_b,1],nz_ab[a,i_b,2],nz_ab[a,i_b,3],a,nz_ab[a,i_b,4],2]
        dist3_nonzero[nz,5] = dist_TT[nz_ab[a,i_b,1],nz_ab[a,i_b,2],nz_ab[a,i_b,3],a,nz_ab[a,i_b,4],3]
        dist3_nonzero[nz,6] = dist_TT[nz_ab[a,i_b,1],nz_ab[a,i_b,2],nz_ab[a,i_b,3],a,nz_ab[a,i_b,4],4]
        
        dist3_nonzero[nz,7] = dist_TT[nz_ab[a,i_c,1],nz_ab[a,i_c,2],nz_ab[a,i_c,3],a,nz_ab[a,i_c,4],2]
        dist3_nonzero[nz,8] = dist_TT[nz_ab[a,i_c,1],nz_ab[a,i_c,2],nz_ab[a,i_c,3],a,nz_ab[a,i_c,4],3]
        dist3_nonzero[nz,9] = dist_TT[nz_ab[a,i_c,1],nz_ab[a,i_c,2],nz_ab[a,i_c,3],a,nz_ab[a,i_c,4],4]
        
        for i = 1:3
            temp = -coords_ab_TT[i,nz_ab[a,i_b,4], nz_ab[a,i_c,4]] +At[i,1]*(rf1[nz_ab[a,i_b,1]]-rf1[nz_ab[a,i_c,1]]) + At[i,2]*(rf2[nz_ab[a,i_b,2]]-rf2[nz_ab[a,i_c,2]] )  + At[i,3]*(rf3[nz_ab[a,i_b,3]]-rf3[nz_ab[a,i_c,3]])
            dist3_nonzero[nz,3] += temp^2
            dist3_nonzero[nz,9+i] += temp
            #                        println("temp $temp")
            
        end
        dist3_nonzero[nz,3] = dist3_nonzero[nz,3]^0.5
        
    end
    
    @turbo for nz = 1:counter
        a = NZ[nz,1]
        i_b = NZ[nz,2]
        i_c = NZ[nz,3]

        for i = 1:3
            dist3_nonzero[nz,9+i] = -dist3_nonzero[nz,9+i] / (dist3_nonzero[nz,3] + 1e-20)
        end
        
        cut3 = cutoff_arr3[a, nz_ab[a,i_b,4],nz_ab[a,i_c,4]]
        
        cut_ab = cutoff_fn_fast(dist3_nonzero[nz,1], cutoff_arr[a, nz_ab[a,i_b,4], 1] - cutoff_length, cutoff_arr[a, nz_ab[a,i_b,4], 1])
        cut_ac = cutoff_fn_fast(dist3_nonzero[nz,2], cut3 - cutoff_length, cut3)
        cut_bc = cutoff_fn_fast(dist3_nonzero[nz,3], cut3 - cutoff_length, cut3)
        
        
        cut_ab2 = cutoff_fn_fast(dist3_nonzero[nz,1], cut3 - cutoff_length, cut3)
        
        dist3_nonzero[nz,13] = cut_ab*cut_bc*cut_ac
        dist3_nonzero[nz,14] = cut_ab2*cut_bc*cut_ac
        
    end            
    
    s = Set(crys.stypes)
    dmin_types3_TT = Dict()
    for s1 in s
        for s2 in s
            for s3 in s
                ss = Set([s1,s2,s3])
                dmin_types3_TT[ss] = 1000000000000.0
                for count in 1:size(dist3_nonzero,1)
                    a = nz_ind3[count,1]
                    b = nz_ind3[count,2]
                    c = nz_ind3[count,3]
                    if crys.stypes[a] == s1 && crys.stypes[b] == s2 && crys.stypes[c] == s3 && (dist3_nonzero[count,1]+dist3_nonzero[count,2]+dist3_nonzero[count,3])/3.0 > 1e-3
                        dmin_types3_TT[ss] = min(dmin_types3_TT[ss], (dist3_nonzero[count,1]+dist3_nonzero[count,2]+dist3_nonzero[count,3])/3.0)
                    end
                end
            end
        end
    end

    nz_ind4 = zeros(UInt16 ,counter,7)

    if fourbody

        nz_ab = ones(Int64, nat, size(nz_ints,1)*nat,4)
        max_ind_a = zeros(UInt16, nat)
        
        for a = 1:nat
            c_ab = 0
            for i = 1:size(nz_ints,1)
                for b = 1:nat
                    if dist_TT[nz_ints[i,1],nz_ints[i,2],nz_ints[i,3],a,b,1] < cutoff4
                        c_ab += 1
                        nz_ab[a, c_ab,1] = nz_ints[i,1]
                        nz_ab[a, c_ab,2] = nz_ints[i,2]
                        nz_ab[a, c_ab,3] = nz_ints[i,3]
                        nz_ab[a, c_ab,4] = b
                        max_ind_a[a] = c_ab
                    end
                end
            end
        end             

        max_a = maximum(max_ind_a)

        #recalculate
        dist3a = zeros(var_type, nat,max_a, max_a)
        @turbo for a = 1:nat
            for i_b = 1:max_a
                for i_c = 1:max_a
                    for i = 1:3
                        temp = -coords_ab_TT[i,nz_ab[a,i_b,4], nz_ab[a,i_c,4]] +At[i,1]*(rf1[nz_ab[a,i_b,1]]-rf1[nz_ab[a,i_c,1]]) + At[i,2]*(rf2[nz_ab[a,i_b,2]]-rf2[nz_ab[a,i_c,2]] )  + At[i,3]*(rf3[nz_ab[a,i_b,3]]-rf3[nz_ab[a,i_c,3]])
                        dist3a[a,i_b,i_c] += temp^2
                    end
                    dist3a[a,i_b,i_c] = dist3a[a,i_b,i_c]^0.5
                    
                end
            end
        end

        
        counter = 0
        NZ4 = zeros(UInt16, nat*max_a*max_a*max_a, 4)
        @inbounds @fastmath @simd for a = 1:nat
            for i_b = 1:max_a
                d12 = dist_TT[nz_ab[a,i_b,1],nz_ab[a,i_b,2],nz_ab[a,i_b,3],a,nz_ab[a,i_b,4],1]
                if i_b <= max_ind_a[a] && d12 > 1e-3 && d12 < cutoff4
                    for i_c = 1:max_a 
                        
                        d13 = dist_TT[nz_ab[a,i_c,1],nz_ab[a,i_c,2],nz_ab[a,i_c,3],a,nz_ab[a,i_c,4],1]
                        d23 = dist3a[a,i_b,i_c]
                        
                        #                        if i_c <= max_ind_a[a] && d12 > 1e-3 && d13 > 1e-3 && d23  > 1e-3 && d12 < cutoff_arr[a, nz_ab[a,i_b,4],1] && d13 < cut3 && d23 < cut3
                        if i_c <= max_ind_a[a] && d13 > 1e-3 && d23  > 1e-3  && d13 < cutoff4 && d23 < cutoff4
                            # && ( exp(-aX2*d13)*exp(-aX2*d23)*1000 > 0.2e-6 || exp(-aX2*d13)*exp(-aX2*d12)*1000 > 0.2e-6)
                            for i_d = 1:max_a 
                                d14 = dist_TT[nz_ab[a,i_d,1],nz_ab[a,i_d,2],nz_ab[a,i_d,3],a,nz_ab[a,i_d,4],1]
                                d24 = dist3a[a,i_b,i_d]
                                d34 = dist3a[a,i_c,i_d]
                                if i_d <= max_ind_a[a] && d14 > 1e-3 && d24  > 1e-3  && d34  > 1e-3  &&  d14  < cutoff4  && d24 < cutoff4 && d34 < cutoff4
                                    counter += 1
                                    NZ4[counter,1] = a
                                    NZ4[counter,2] = i_b
                                    NZ4[counter,3] = i_c
                                    NZ4[counter,4] = i_d
                                end
                            end
                        end
                    end
                end
            end
        end

        begin
            dist4_nonzero = zeros(var_type, counter,3+3*4+2)
            nz_ind4 = zeros(UInt16 ,counter,7)
        end
        cx = [0,0,0]
        for nz = 1:counter
            a = NZ4[nz,1]
            i_b = NZ4[nz,2]
            i_c = NZ4[nz,3]
            i_d = NZ4[nz,4]
            
            nz_ind4[nz, 1] = a
            nz_ind4[nz, 2] = nz_ab[a,i_b,4]
            nz_ind4[nz, 3] = nz_ab[a,i_c,4]
            nz_ind4[nz, 4] = nz_ab[a,i_d,4]
            cx[1] = nz_ab[a,i_b,1]-R[1]-1
            cx[2] = nz_ab[a,i_b,2]-R[2]-1
            cx[3] = nz_ab[a,i_b,3]-R[3]-1
            nz_ind4[nz, 5] = R_dict_tt[cx]
            cx[1] = nz_ab[a,i_c,1]-R[1]-1
            cx[2] = nz_ab[a,i_c,2]-R[2]-1
            cx[3] = nz_ab[a,i_c,3]-R[3]-1
            nz_ind4[nz, 6] = R_dict_tt[cx]

            cx[1] = nz_ab[a,i_d,1]-R[1]-1
            cx[2] = nz_ab[a,i_d,2]-R[2]-1
            cx[3] = nz_ab[a,i_d,3]-R[3]-1
            nz_ind4[nz, 7] = R_dict_tt[cx]
            
        end
        

        #       return nz_inds, R_keep_ab_TT, nz_ind3, dist3_nonzero, dist_arr_TT, c_zero_tt, dmin_types_TT, dmin_types3_TT, nz_ind4

        
    end #end fourbody

    return nz_ind3, dist3_nonzero, dmin_types3_TT, nz_ind4 

end
#return nz_inds, R_keep_ab_TT, nz_ind3, dist3_nonzero, dist_arr_TT, c_zero_tt, dmin_types_TT, dmin_types3_TT


    #    return R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3


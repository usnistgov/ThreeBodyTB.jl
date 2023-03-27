###include("Crystal.jl")
###include("DFToutMod.jl")
#using XMLDict



"""
    module CalcTB

Create TB matrix sets from coefficients, or prepare to fit the coefficients.
"""
module Classical


using LinearAlgebra
using SpecialFunctions
using GZip
using EzXML
using XMLDict
using Printf
using SparseArrays
using LoopVectorization
using ..CrystalMod:distances_etc_3bdy_parallel_LV
using ..CalcTB:laguerre_fast!

#using Statistics

using Base.Threads
#using Base.lock

using ..AtomicMod:atom
using ..CrystalMod:crystal
#using ..TB:orbital_index
using ..CrystalMod:orbital_index

using ..CrystalMod:get_dist

using ..Utility:cutoff_fn
using ..Utility:cutoff_fn_fast
using ForwardDiff
using ..CrystalMod:makecrys

using Random


const n_2body_cl = 6
const n_3body_cl = 11
const n_3body_cl_same = 7



using ..CrystalMod:cutoff2X
using ..CrystalMod:cutoff3bX
using ..CrystalMod:cutoff_onX
using ..CrystalMod:cutoff_length

#const cutoff2X = 18.51 
#const cutoff3bX = 13.01 
#const cutoff_onX = 18.01 

#const cutoff_length = 1.0


"""
    struct coefs

Hold the TB coefficients for 2body or 3body interactions

- `dim::Int64` - dimension (2 or 3)
- `datH::Array{Float64,1}` Hamiltonian parameters
- `datS::Array{Float64,1}`  Overlap parameters, if any
- `sizeH::Int64` 
- `sizeS::Int64`
- `inds::Dict{Array{Symbol}, Array{Int64,1}}` Holds inds that tell which coeffiecients correspond to which atoms/orbitals
- `inds_int::Dict{Array{Symbol}, Any}` Holds inds that tell which coeffiecients correspond to which atoms/orbitals, but integer array
- `ninds_int::Dict{Array{Symbol}, Array{UInt16,2}}` Holds number of coefs
- `names::Set` Names of atoms
- `orbs::Array{Any,1}` Orbitals 
- `cutoff::Float64` cutoff distance
- `min_dist::Float64` minimum atom-atom distance in fitting data
- `maxmin_val_train::Dict` max and min value of matrix elements within fitting data
- `dist_frontier::Dict` dictionary of pareto frontier of shortest fitting distances
- `version::Int64` version number
"""
struct coefs_cl

    dim::Int64
    datH::Array{Float64,1}
    sizeH::Int64
    names::Set
    min_dist::Float64
    dist_frontier::Dict
    version::Int64
    inds::Dict{Array{Symbol}, Array{Int64,1}}
end


"""
    function write_coefs(filename, co::coefs; compress=true)

Write `coefs` to a file. Compress uses gzip. See `read_coefs`
"""
function write_coefs_cl(filename, co::coefs_cl; compress=true)
    """
    write xml coefs object
    """

    doc = XMLDocument()
    root = ElementNode("root")
    setroot!(doc, root)

    c = ElementNode("coefs")
    link!(root, c)

    addelement!(c, "version", string(co.version))

    addelement!(c, "dim", string(co.dim))
    addelement!(c, "datH", arr2str(co.datH))
    addelement!(c, "sizeH", string(co.sizeH))

    addelement!(c, "names", str_w_spaces(co.names))
#    addelement!(c, "orbs", str_w_spaces(co.orbs))
    addelement!(c, "cutoff", string(co.cutoff))
    addelement!(c, "min_dist", string(co.min_dist))
    addelement!(c, "dist_frontier", dict2str(co.dist_frontier))
    

    if compress
        io=gzopen(filename*".gz", "w")
    else
        io=open(filename, "w")
    end
    prettyprint(io, doc);
    close(io)

    return Nothing
    
end


"""
    function read_coefs(filename, directory = missing)

Read `coefs` from filename. Can read gzip directly.
"""
function read_coefs(filename, directory = missing)
    """
    read xml coefs object
    """
    if !ismissing(directory)
        filename=directory*"/"*filename
    end
    if !isfile(filename)
        if isfile(filename*".xml")
            filename=filename*".xml"
        elseif isfile(filename*".gz")
            filename=filename*".gz"
        elseif isfile(filename*".xml.gz")
            filename=filename*".xml.gz"
        end
    end
    if !isfile(filename)
        println("warning error read_coefs $filename $directory not found")
    else
        println("found $filename")
    end
    f = gzopen(filename, "r")
    
    fs = read(f, String)
    close(f)

    d = xml_dict(fs)["root"]
    
    dim = parse(Int64, (d["coefs"]["dim"]))
    sizeH = parse(Int64, d["coefs"]["sizeH"])

    if sizeH > 0
        datH = parse_str_ARR_float(d["coefs"]["datH"])
    else
        datH = Float64[]
    end
    
    names = Set(String.(split(d["coefs"]["names"])))
    min_dist = parse(Float64, d["coefs"]["min_dist"])
    dist_frontier = str2tuplesdict(eval(d["coefs"]["dist_frontier"]))

    println("version $version")
    
    co = make_coefs_cl(names,dim, datH=datH, datS=datS, min_dist=min_dist, dist_frontier = dist_frontier, version=version)

    return co
    
end


function test(at_list, dim)
    println( at_list)
    println(dim)
end

"""
    function make_coefs(at_list, dim; datH=missing, datS=missing, cutoff=18.01, min_dist = 3.0, fillzeros=false, maxmin_val_train=missing, dist_frontier=missing)

Constructor for `coefs`. Can create coefs filled with ones for testing purposes.

See `coefs` to understand arguments.
"""
function make_coefs_cl(at_list, dim; datH=missing, min_dist = 3.0, fillzeros=false, dist_frontier=missing, version=3)

    println("make coefs")
    if dim == 2
        totH = 6
    else
        println("asdf")
    end
    
    if ismissing(datH)
        if fillzeros
            datH = zeros(totH)
        else
            datH = ones(totH) 
        end
    end

    
    dist_frontier2 = Dict()
    if !ismissing(dist_frontier)
        
        for key in keys(dist_frontier)
#            println("key ", key, " at_list ", at_list)
            if dim == length(key)
                dist_frontier2[key] = dist_frontier[key]
                dist_frontier2[Symbol.(key)] = dist_frontier[key]
            end
        end
  end


    return coefs_cl(dim, datH, totH, Set(at_list), min_dist, dist_frontier2, version)
    
    
end
    

function get_cl_info(at_list::Set, dim)

    at_list = Symbol.([i for i in at_set])
    sort!(at_list)
    inds_dict = Dict{Any, Array{Int64,1}}()
    if dim == 2
        if length(at_list) == 1
            totH = n_2body_cl
            inds_dict[(at_list[1], at_list[2])] = collect(1:totH)
        elseif length(at_list) == 2
            totH = n_2body_cl
            inds_dict[(at_list[1], at_list[2])] = collect(1:totH)
            inds_dict[(at_list[2], at_list[1])] = collect(1:totH)
        else
            println("$dim error length(at_list) == ", length(at_list))
        end
    elseif dim == 3
        if length(at_list) == 1
            totH = n_3body_cl_same
            inds_dict[(at_list[1], at_list[1], at_list[1])] = collect(1:totH)
        elseif length(at_list) == 2
            permutations = [[at_list[1], at_list[2], at_list[2]],
                            [at_list[2], at_list[1], at_list[2]],
                            [at_list[2], at_list[2], at_list[1]],
                            [at_list[1], at_list[1], at_list[2]],
                            [at_list[1], at_list[2], at_list[1]],
                            [at_list[2], at_list[1], at_list[1]]]
            totH = 0
#            for p in permutations
#                inds_dict[

            #inds_dict[(at_list[1], at_list[1], at_list[2])] = totH = n_3body_cl_same

        elseif length(at_list) == 3
            totH = n_3body_cl
        else
            println("$dim error length(at_list) == ", length(at_list))
        end
        
    else
        println("error get_cl_info dim $dim")
    end

    return totH, inds_dict
    

end


Base.show(io::IO, d::coefs_cl) = begin

    println(io, "coeffs classical ", d.names)
    println(io, "min dist ", d.min_dist)
    println(io, "dist frontier ")
    for t in keys(d.dist_frontier)
        println(io, t, "    ", d.dist_frontier[t])
    end

    println(io)
end


function get_fit_cl(CRYS; use_threebody=true)

    types_dict = Dict()
    types_dict_reverse = Dict()
    types_counter = 0


    vars = Set()
    at_types_set = Set()
    for crys in CRYS
        st = Set(crys.stypes)
        for s1 in st
            push!(at_types_set, s1)
            for s2 in st
                push!(vars,[Set([s1,s2]),2])
                if use_threebody
                    for s3 in st
                        push!(vars,[Set([s1,s2]), 3])
                    end
                end
            end
        end
    end
    vars_list = collect(vars)
    at_types = collect(at_types_set)
    ntot = 0
    at_ind = []
    for v in vars_list
        if v[2] == 2
            ntot += 6
        end
    end

#    for (ind, crys) in enumerate(CRYS)
#        en = calc_energy_cl(crys, dat_vars=ones(ntot), at_types=at_types, vars_list= vars_list, use_threebody=use_threebody)
#        println("ind $ind en $en")
#    end
    
    V = zeros(length(CRYS), ntot)
    for (ind, crys) in enumerate(CRYS)
        V[ind, :] = ForwardDiff.gradient(x->calc_energy_cl(crys, dat_vars=x, at_types=at_types, vars_list= vars_list, use_threebody=use_threebody), ones(ntot) )
    end

    return V
    
end



function calc_energy_cl(crys::crystal;  database=missing, dat_vars=missing, at_types = missing, vars_list = missing,  DIST=missing, verbose=true, use_threebody=false)


    At = crys.A'
    
    if verbose
        println()
        println("-----")
        println("Get classical energy from crystal structure")
        println()
    end

    if ismissing(dat_vars)
        var_type=Float64
    else
        var_type = eltype(dat_vars)

    end

    if verbose  println("LV $var_type")  end
    
    
    

    
    ind2orb, orb2ind, etotal, nval = orbital_index(crys)
    if verbose println("distances") end
    begin
        
        use_dist_arr = true
        if !ismissing(DIST)
            use_dist_arr = false
            R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3 = DIST

        else
            if (use_threebody ) && !ismissing(database)
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
                            println(key," " ,key2, " : ", dmin_types[key], " <~ ", database[key2].min_dist)
                            within_fit = false
                        end
                    end
                end
            end
            
#            c_zero_ref=1
#            if !(ismissing(reference_tbc))
#                if size(reference_tbc.tb.ind_arr)[1] > 1
#                    c_zero_ref = reference_tbc.tb.r_dict[[0,0,0]]
#                end
#            end
        end
    end

    if !ismissing(database)
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
    else
        types_dict = Dict()
        types_dict_reverse = Dict()
        types_counter = 0
        types_arr = zeros(Int64, crys.nat)
        for a = at_types
            if !(Symbol(a) in keys(types_dict))
                types_counter += 1
                types_dict[Symbol(a)] = types_counter
                types_dict_reverse[types_counter] = Symbol(a)
            end
        end

        println("types_counter ", types_counter)
        
        for (ind, a) in enumerate(crys.stypes)
            types_arr[ind] = types_dict[a]
        end
        println("types_arr, ", types_arr)
        
    end
    
    
    DAT_IND_ARR = zeros(var_type, types_counter, types_counter,1:n_2body_cl )

    badlist = Set()
    
    if !ismissing(database)
        
        for c1 = 1:types_counter
            for c2 = 1:types_counter
                t1 = types_dict_reverse[c1]
                t2 = types_dict_reverse[c2]
                if (t1,t2) in keys(database)
                    coef = database[(t1,t2)]
                    DAT_IND_ARR[c1,c2,1:n_2body_cl] = coef.datH[:]
                else
                    println("WARNING, ",(t1,t2), " database not found")
                    within_fit = false
                    push!(badlist, (t1,t2))
                    println("badlist ", badlist)
                end
            end
        end
    else
        counter = 0
        for v in vars_list
            println("v $v")
            for c1 = 1:types_counter
                for c2 = 1:types_counter
                    t1 = types_dict_reverse[c1]
                    t2 = types_dict_reverse[c2]
                    println("t1 $t1 t2 $t2 ", typeof(t1))
                    if t1 in v[1] && t2 in v[1] && length(Set([t1,t2])) == length(v[1])
                        DAT_IND_ARR[c1,c2,1:n_2body_cl] = dat_vars[counter.+(1:n_2body_cl)]
                    end
                end
            end
            counter += n_2body_cl
        end
    end
        
    println("DAT_IND_ARR")
    println(DAT_IND_ARR)

    if verbose println("2body CL") end

    nkeep_ab = size(R_keep_ab)[1]
    
    lag_arr_TH = zeros(var_type, 6, nthreads())
    energy_2bdy_TH = zeros(var_type,nthreads())
    twobody_Cl = begin
        for c = 1:nkeep_ab
            id = threadid()
            lag_arr = lag_arr_TH[:,id]

            cind = R_keep_ab[c,1]
            cham = R_keep_ab[c,7]
            a1 = R_keep_ab[c,2]
            a2 = R_keep_ab[c,3]

            t1 = types_arr[a1]
            t2 = types_arr[a2]
            
            if use_dist_arr
                dist_a = dist_arr[c,1]
                cut_a = dist_arr[c,5]
            else
                dist_a, _ = get_dist(a1,a2, R_keep_ab[c,4:6], crys, At)
                cutoff2, cutoff2on = cutoff_arr[a1,a2,:]
                cut_a = cutoff_fn_fast(dist_a, cutoff2 - cutoff_length, cutoff2)
            end
            
            if dist_a <  1e-5    # true onsite
                continue
            end
            
            laguerre_fast!(dist_a, lag_arr)
            energy_2bdy_TH[id] += core_cl(t1, t2, lag_arr, DAT_IND_ARR, var_type) * cut_a
        end
        
    end
        
    threebdy_LV = begin
        
        if use_threebody 

            meta_count = []
            old = 1
            for  counter = 1: (size(array_ind3)[1]-1 ) 
                if array_ind3[counter,1] != array_ind3[counter+1,1]
                    push!(meta_count, old:counter)
                    old = counter+1
                end
            end
            push!(meta_count, old:size(array_ind3)[1])
            
            if verbose println("3body LV") end
            
            for mc in meta_count #@threads 
                for counter in mc
                    
                    id = threadid()

                    a1 = array_ind3[counter,1]
                    a2 = array_ind3[counter,2]
                    a3 = array_ind3[counter,3]

                    t1s = crys.stypes[a1]
                    t2s = crys.stypes[a2]
                    t3s = crys.stypes[a3]

                    t1 = types_arr[a1]
                    t2 = types_arr[a2]
                    t3 = types_arr[a3]
                    
                    if (t1,t2,t3) in badlist
                        continue
                    end
                    
                    cind1 = array_ind3[counter,4]
                    cind2 = array_ind3[counter,5]

                    rind1 = ind_arr[cind1,1:3]
                    rind2 = ind_arr[cind2,1:3]

                    
                    if use_dist_arr
                        dist12 = array_floats3[counter, 1]
                        dist13 = array_floats3[counter, 2]
                        dist23 = array_floats3[counter, 3]                    

                        cut_h = array_floats3[counter,13]
                    else
                        dist12, _ = get_dist(a1,a2, rind1, crys, At)
                        dist13, _ = get_dist(a1,a3, rind2, crys, At)
                        dist23, _= get_dist(a2,a3, -rind1+rind2, crys, At)
                        
                        cutoff3 = cutoff_arr3[a1,a2,a3]

                        cut_ab = cutoff_fn_fast(dist12, cutoffZZ - cutoff_length, cutoff3)
                        cut_ac = cutoff_fn_fast(dist13, cutoff3 - cutoff_length, cutoff3)
                        cut_bc = cutoff_fn_fast(dist23, cutoff3 - cutoff_length, cutoff3)
                        
                        cut_h = cut_ab*cut_ac*cut_bc
                        
                    end

                    

                    memory = memory_TH[:,id]

                    
                    if use_threebody
                        laguerre_fast_threebdy_classical!(dist12,dist13,dist23, t1==t2, t1 !=t2 && t1 != t3 && t2 != t3, memory)
                        core3b!(cind1,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_3, memory, DAT_ARR_3, cut_h, H,sym_arr1, sym_arr2, lmn13, lmn23)
                        
                        
                        #                        core3a!(cind1,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_3, memory, DAT_ARR_3, cut_h, H_thread3,id, sym_arr1, sym_arr2, lmn13, lmn23)
                        #                        println("size ", size(H), size(H_thread2))
                        #                        println("x $a1 $a2 ", orbs_arr[a1,1,1]:orbs_arr[a1,norb[a1],1], " " , orbs_arr[a2,1,1]:orbs_arr[a2,norb[a2],1])

                        #    H[orbs_arr[a1,1,1]:orbs_arr[a1,norb[a1],1],orbs_arr[a2,1,1]:orbs_arr[a2,norb[a2],1], cind1] += @view H_thread2[1:norb[a1],1:norb[a2],1]

                    end
                    

                    if use_threebody_onsite
                        laguerre_fast_threebdy_onsite!(dist12,dist13,dist23, t1==t2 && t1 == t3, memory)
                        #core_onsite3!(c_zero,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_onsite_3, memory, DAT_ARR_3, cut_o, H_thread, id)

                        #println("cuto $cut_o ",  [dist12,dist13,dist23], " ", exp(-1.0*dist13)*exp(-1.0*dist23)*1000, " ", exp(-1.0*dist13)*exp(-1.0*dist23)*1000*exp(-1.0*dist12))

                        core_onsite3b!(c_zero,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_onsite_3, memory, DAT_ARR_3, cut_o, H)
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
    end

    if retmat
        #return H, S
        return 
    end
    if verbose println("make") end
    if true
        #        println("typeof H ", typeof(H), " " , size(H), " S ", typeof(S), " " , size(S))
        #println("maketb")
        tb = make_tb( reshape(H, 1,size(H)[1], size(H)[2], size(H)[3])  , ind_arr, S)
        if !ismissing(database) && (haskey(database, "scf") || haskey(database, "SCF"))
            scf = database["scf"]
        else
            scf = false
        end
        #println("make")
        tbc = make_tb_crys(tb, crys, nval-tot_charge, 0.0, scf=scf, gamma=gamma, background_charge_correction=background_charge_correction, within_fit=within_fit, screening=screening)
    end
    if verbose 
        println("-----")
        println()
    end


    return tbc

end

function laguerre_fast_cl!(dist, memory)

    a=2.2
    ad = a*dist
    expa=exp.(-0.5*ad)
    memory[1] = 1.0 * expa
    memory[2] = (1.0 .- ad) .* expa
    memory[3]= 0.5*(ad.^2 .- 4.0*ad .+ 2) .* expa
    memory[4] = 1.0/6.0*(-ad.^3 .+ 9.0*ad.^2 .- 18.0*ad .+ 6.0) .* expa
    memory[5] = 1.0/24.0*(ad.^4 .- 16.0 * ad.^3 .+ 72.0*ad.^2 .- 96.0*ad .+ 24.0) .* expa
    memory[6] = 1.0/120*(-ad.^5 .+ 25*ad.^4 .- 200 * ad.^3 .+ 600.0*ad.^2 .- 600.0*ad .+ 120.0) .* expa

end


function core_cl(t1, t2, lag_arr, DAT_ARR, var_type )
    energy = zero(var_type)
    for n = 1:n_2body_cl
        energy +=  lag_arr[n]*DAT_ARR[t1,t2,n]
    end
    return energy
end

function core_onsite!(c_zero, a1, a2, t1, t2, norb, orbs_arr, DAT_IND_ARR_O, lag_arr, DAT_ARR, cut_on, H,  lmn, sym_dat1, sym_dat2)

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

            if (sum1 == 2 && sum2 == 2) || (sum1 == 3 && sum2 == 3)
                for n = 1:5 #DAT_IND_ARR[t1,t2,1,orbs_arr[a1,o1x,2],orbs_arr[a2,o2x,2],1]
                    temp2 +=  lag_arr[n]*DAT_ARR[t1,t2,1,DAT_IND_ARR_O[t1,t2,sum1,sum2,n+1+5]  ]
                end
                temp2 = temp2* symmetry_factor_int(s1, 1, lmn, one)*symmetry_factor_int(s2, 1, lmn, one)

            end

            H[ o1, o2, c_zero] += (temp1 + temp2)  * cut_on

            
        end
    end
end

function core3!(cind1,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_3, memory, DAT_ARR_3, cut_h, H_thread, id, sym_dat1, sym_dat2, lmn31, lmn32)

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
            H_thread[o1,o2,cind1,id] += temp * sym31*sym32*10^3 * cut_h
            #temp = temp * sym31*sym32*10^3 * cut_h
            
        end
    end
            
end

function core3b!(cind1,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_3, memory, DAT_ARR_3, cut_h, H_thread, sym_dat1, sym_dat2, lmn31, lmn32)

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
                H_thread[o1,o2,cind1] += temp * sym31*sym32*10^3 * cut_h
           # end
            #println("add $o1 $o2 $cind1   ", temp * sym31*sym32*10^3 * cut_h)
            
            #temp = temp * sym31*sym32*10^3 * cut_h
            
        end
    end
            
end

function core3a!(cind1,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_3, memory, DAT_ARR_3, cut_h, H_thread, id,sym_dat1, sym_dat2, lmn31, lmn32)
    H_thread .= 0.0
    for o2x = 1:norb[a2]
        id = threadid()
        o2 = orbs_arr[a2,o2x,1]
        s2 = orbs_arr[a2,o2x,2]
        sum2 = orbs_arr[a2,o2x,3]

        sym32 = symmetry_factor_int(s2,1,lmn32, one ) 

        @simd for o1x = 1:norb[a1]

            o1 = orbs_arr[a1,o1x,1]
            s1 = orbs_arr[a1,o1x,2]
            sum1 = orbs_arr[a1,o1x,3]

            
            sym31 = symmetry_factor_int(s1,1,lmn31, one )    

                      
            temp = 0.0
            for i = 1:DAT_IND_ARR_3[t1,t2,t3, sum1, sum2,1]
                temp += memory[i] * DAT_ARR_3[t1,t2,t3, DAT_IND_ARR_3[t1,t2,t3, sum1, sum2,1+i]]
            end
            H_thread[o1x,o2x,1] += temp * sym31*sym32*10^3 * cut_h
            #temp = temp * sym31*sym32*10^3 * cut_h
            
        end
    end
    
end

function core_onsite3!(c_zero,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_onsite_3, memory, DAT_ARR_3, cut_o, H_thread, id)

    #@inbounds @fastmath @simd     
    @inbounds @simd for o1x = 1:norb[a1]

        o1 = orbs_arr[a1,o1x,1]
        s1 = orbs_arr[a1,o1x,2]
        sum1 = orbs_arr[a1,o1x,3]
        
        
        temp = 0.0
        for i = 1:DAT_IND_ARR_onsite_3[t1,t2,t3, sum1,1 ]
            temp += memory[i] * DAT_ARR_3[t1,t2,t3, DAT_IND_ARR_onsite_3[t1,t2,t3, sum1, 1+i]]
        end
#            println(DAT_IND_ARR_3[t1,t2,t3, sum1, sum2,1], ", $a1 $a2 $a3 | $o1x $o2x | temp $temp | sum(mem) ", sum(memory), " sum DAT_ARR " , sum(DAT_ARR_3[t1,t2,t3,sum1, sum2, :]))
#        println("$a1 $o1 o  $(temp *10^3 * cut_o)   $cut_o")
        H_thread[o1,o1,c_zero,id] += temp *10^3 * cut_o
    end
    
end

function core_onsite3b!(c_zero,  a1, a2, a3, t1, t2, t3, norb, orbs_arr, DAT_IND_ARR_onsite_3, memory, DAT_ARR_3, cut_o, H_thread)

    #@inbounds @fastmath @simd     
    @inbounds @simd for o1x = 1:norb[a1]

        o1 = orbs_arr[a1,o1x,1]
        s1 = orbs_arr[a1,o1x,2]
        sum1 = orbs_arr[a1,o1x,3]
        
        
        temp = 0.0
        for i = 1:DAT_IND_ARR_onsite_3[t1,t2,t3, sum1,1 ]
            temp += memory[i] * DAT_ARR_3[t1,t2,t3, DAT_IND_ARR_onsite_3[t1,t2,t3, sum1, 1+i]]
        end
#            println(DAT_IND_ARR_3[t1,t2,t3, sum1, sum2,1], ", $a1 $a2 $a3 | $o1x $o2x | temp $temp | sum(mem) ", sum(memory), " sum DAT_ARR " , sum(DAT_ARR_3[t1,t2,t3,sum1, sum2, :]))
#        if abs(temp *10^3 * cut_o) > 1e-7
#            println("$a1 $o1 o  $(temp *10^3 * cut_o)   $cut_o")
            H_thread[o1,o1,c_zero] += temp *10^3 * cut_o
        #        end
    end
    
end


#function twobody(nkeep_ab::Int64)
#    @threads for c = 1:nkeep_ab
#        begin
#        end
#    end#
#
#    
#end

#include("CalcTB_laguerre_deriv.jl")



end #end module




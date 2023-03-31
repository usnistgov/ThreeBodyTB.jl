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
using ..CalcTB:calc_frontier_list
using Suppressor

#using Statistics

using Base.Threads
#using Base.lock

using ..AtomicMod:atom
using ..CrystalMod:crystal
#using ..TB:orbital_index
using ..CrystalMod:orbital_index

using ..CrystalMod:get_dist
using ..CrystalMod:get_dist_only

using ..Utility:cutoff_fn
using ..Utility:cutoff_fn_fast
using ForwardDiff
using ..CrystalMod:makecrys

using Random


const n_2body_cl = 6
#const n_3body_cl_same = 7
#const n_3body_cl_pair = 13
#const n_3body_cl_diff = 20

const n_3body_cl_same = 9
const n_3body_cl_pair = 13
const n_3body_cl_diff = 29

#const n_3body_cl_same = 4
#const n_3body_cl_pair = 13
#const n_3body_cl_diff = 8



using ..CrystalMod:cutoff2X
using ..CrystalMod:cutoff3bX
using ..CrystalMod:cutoff_onX
using ..CrystalMod:cutoff_length

using ..Force_Stress:reshape_vec
using ..Atomdata:get_cutoff

using ..DFToutMod:dftout

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



"""
    function make_coefs(at_list, dim; datH=missing, datS=missing, cutoff=18.01, min_dist = 3.0, fillzeros=false, maxmin_val_train=missing, dist_frontier=missing)

Constructor for `coefs`. Can create coefs filled with ones for testing purposes.

See `coefs` to understand arguments.
"""
function make_coefs_cl(at_list, dim; datH=missing, min_dist = 3.0, fillzeros=false, dist_frontier=missing, version=3)

    println("make coefs")
    if dim == 2
        totH = n_2body_cl
    elseif dim == 3

        totH = n_3body_cl_diff
        #        t1,t2,t3 = at_list[:]
#        if t1 == t2 && t2 == t3
#            totH = n_3body_cl_same
#        elseif t1 == t2 && t1 != t3 || t1 == t3 && t1 != t2 || t2 == t3 && t2 != t1
#            totH = n_3body_cl_pair
#        elseif t1 != t2 && t2 != t3 && t1 != t3
#            totH = n_3body_cl_diff
#        else
#            println("something wrong make coefs cl $t1 t2 $t3 ", [ t1 == t2 && t2 == t3,t1 == t2 && t1 != t3 || t1 == t3 && t1 !##= t2 || t2 == t3 && t2 != t1,t1 != t2 && t2 != t3 && t1 != t3])
#        end
        
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
    

Base.show(io::IO, d::coefs_cl) = begin

    println(io, "coeffs classical ", d.names)
    if d.dim == 2
        println(io, "min dist ", round.(d.min_dist, digits=4))
    end
    println(io, "dist frontier ")
    for t in keys(d.dist_frontier)
        println(io, t)
        for temp in d.dist_frontier[t]
            println(io,  round.(temp, digits=2))
        end
    end
    println(io)
    println(io, "datH ", round.(d.datH[:], digits=3))
    println(io)
end

function atom_perm_to_dist_perm(perm)

    d12 = sort([perm[1], perm[2]])
    d13 = sort([perm[1], perm[3]])
    d23 = sort([perm[2], perm[3]])

    normal_order = [[1,2], [1,3], [2,3]]

    new_order = [findfirst(x->x == normal_order[1], [d12,d13,d23]),
                 findfirst(x->x == normal_order[2], [d12,d13,d23]),
                 findfirst(x->x == normal_order[3], [d12,d13,d23])]

    
#    new_ind = Int64[]
#    for n in new_order
#        push!(findfirst(n, normal_order), new_ind)
#    end

    return new_order
    
    
end


function do_fit_cl(DFT::Array{Any,1}; use_threebody=true, energy_weight=1.0, use_energy=true, use_force=true, use_stress = true , database=missing, lambda = -1.0)

    if !use_energy && !use_force && !use_stress
        println("everything is false")
        return Dict()
    end
    
    CRYS = crystal[]
    if use_energy
        ENERGIES = Float64[]
    else
        ENERGIES=missing
    end
    if use_force
        FORCES = Float64[]
    else
        FORCES = missing
    end
    if use_stress
        STRESSES = Float64[]
    else
        STRESSES = missing
    end
    
    for dft in DFT
        push!(CRYS, dft.crys)

        #subtract old forces
        if !ismissing(database)
            energyT, forceT, stressT = energy_force_stress(dft.crys, database=database, use_threebody=use_threebody, verbose=false)
        else
            energyT = 0.0
            forceT = zeros(dft.crys.nat, 3)
            stressT= zeros(3,3)
        end
        
        if use_energy
            ENERGIES = vcat(ENERGIES, [dft.atomize_energy/dft.crys.nat - energyT/dft.crys.nat])
        end
        if use_force
            FORCES = vcat(FORCES, dft.forces[:] - forceT[:])
        end
        if use_stress
            stress = dft.stress - stressT
            STRESSES = vcat(STRESSES, [stress[1,1], stress[1,2],stress[1,3],stress[2,2],stress[2,3],stress[3,3]])
        end
    end
    ENERGIES = ENERGIES * energy_weight

    return do_fit_cl(CRYS, use_threebody=use_threebody, energy_weight = energy_weight, ENERGIES=ENERGIES, FORCES=FORCES, STRESSES=STRESSES, database=database, lambda=lambda)
    
end

function do_fit_cl(CRYS::Array{crystal,1}; use_threebody=true,ENERGIES=missing, FORCES=missing, STRESSES=missing, energy_weight = 1.0, database=missing, lambda = -1.0)

    if !ismissing(FORCES) || !ismissing(STRESSES)
        get_force=true
    else
        get_force = false
    end

    Ve,Vf,Vs, vars, at_types, ind_set = prepare_fit_cl(CRYS; use_threebody=use_threebody, get_force=get_force, database=database)

    ntot = max(size(Ve)[2], size(Vf)[2], size(Vs)[2])
    
    if !ismissing(ENERGIES) && !ismissing(FORCES) && !ismissing(STRESSES)
        println("Use energy, force and stress; energy weight $energy_weight")
        Vtot = vcat(energy_weight*Ve,Vf,Vs)
        Rtot = vcat(ENERGIES, FORCES, STRESSES)
    elseif !ismissing(ENERGIES) && !ismissing(FORCES)
        println("Use energy and force; energy weight $energy_weight")
        Vtot = vcat(energy_weight*Ve,Vf)
        Rtot = vcat(ENERGIES, FORCES)
    elseif !ismissing(FORCES)
        println("Use forces only")
        Vtot = Vf
        Rtot = FORCES
    elseif !ismissing(ENERGIES)
        println("Use energies only")
        Vtot = Ve
        Rtot = ENERGIES
    else
        println("something wrong do_fit_cl ", [ismissing(ENERGIES), ismissing(FORCES), ismissing(STRESSES)])
        return Dict()
    end

    if lambda > 0
        LAM = I(ntot) * lambda
        Vtot = vcat(Vtot, LAM)
        Rtot = vcat(Rtot, zeros(ntot))
    end
    
    x = Vtot \ Rtot
    
    frontier = calc_frontier_list(CRYS)
    if ismissing(database)
        database = Dict()
    end
    println("put stuff in database")
    for ind in keys(ind_set)
        println("ind $ind")
        if length(ind) == 2
            t1,t2 = ind
#            println("t1 t2 $t1 $t2 ind ", ind_set[ind])
#            println("t1 t2 $t1 $t2 x   ", x[ind_set[ind]])
            if (t1,t2) in keys(frontier)
                c = make_coefs_cl([t1, t2], 2, datH = x[ind_set[ind]], dist_frontier = frontier, version = 1, min_dist = frontier[(t1,t2)])
            else
                c = make_coefs_cl([t1, t2], 2, datH = x[ind_set[ind]], version = 1)
            end                
            database[ind] = c
        elseif length(ind) == 3
            t1,t2,t3 = ind
            if (t1,t2,t3) in keys(frontier)
                c = make_coefs_cl([t1, t2, t3], 3, datH = x[ind_set[ind]], dist_frontier = frontier, version = 1)
            else
                c = make_coefs_cl([t1, t2, t3], 3, datH = x[ind_set[ind]], version = 1)
            end
            database[ind] = c
            
        else
            println("something wrong do_fit_cl $ind")            
        end
            
    end
    return database
end    

function prepare_fit_cl(CRYS; use_threebody=true, get_force=true, database=missing)

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
                if !ismissing(database)
                    if ! ( (s1,s2) in keys(database))
                        push!(vars,[Set([s1,s2]),2])
                    end 
                else
                    push!(vars,[Set([s1,s2]),2])
                end
                if use_threebody
                    for s3 in st
                        if !ismissing(database)
                            if ! ( (s1,s2,s3) in keys(database))
                                 push!(vars,[Set([s1,s2,s3]), 3])
                            end
                        else
                            push!(vars,[Set([s1,s2,s3]), 3])
                        end
                    end
                end
            end
        end
    end
    println("vars ", vars)
    
    vars_list = collect(vars)
    at_types = collect(at_types_set)
    ntot = 0
    at_ind = []

    ind_set = Dict()

    pattern = [ [1,1,1], #these are distances 12, 13, 23

                [2,1,1],
                [1,2,1],
                [1,1,2],

                [3,1,1],
                [1,3,1],
                [1,1,3],

                [2,2,2],

                [1,2,2],
                [2,1,2],
                [2,2,1], 

                [0,1,1], 
                [1,0,1], 
                [1,1,0], 

                [0,1,2], 
                [0,2,1], 
                [2,0,1], 
                [1,0,2], 
                [1,2,0], 
                [2,1,0]]


    
    for v in vars_list
        println("v $v")
        println()
        if v[2] == 2
            if length(v[1]) == 1
                t1 = collect(v[1])[1]
                t2 = t1
            elseif length(v[1]) == 2
                t1,t2 = collect(v[1])[1:2]
            else
                println("something wrong v length(v[1]) ", length(v[1]) )
            end
            ind_set[(t1,t2)] = collect(1:n_2body_cl) .+ ntot
            ind_set[(t2,t1)] = collect(1:n_2body_cl) .+ ntot
            ntot += n_2body_cl
        elseif v[2] == 3
            println("three ", v[1])            
            if length(v[1]) == 1
                println("elemental")
                t1 = collect(v[1])[1]
                ind_set[(t1,t1,t1)] = ntot .+ [1,2,2,2,3,3,3,4,5,5,5,6,6,6,7,7,7,7,7,7, 8, 8, 8, 9,9,9,9,9,9]
                #ind_set[(t1,t1,t1)] = ntot .+ [1,2,2,2,3,4,4,4]
                ntot += n_3body_cl_same
                
            elseif length(v[1]) == 2
                t1,t2 = collect(v[1])[1:2]

                #1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20
                ind_set[(t1,t2,t2)  ] = ntot .+ [1, 2, 2, 3, 4, 4, 5, 6, 7, 7, 8, 9, 9,10,13,12,12,13,11,11] #3 special
                ind_set[(t2,t1,t2)  ] = ntot .+ [1, 2, 3, 2, 4, 5, 4, 6, 7, 8, 7, 9,10, 9,12,13,11,11,13,12] #2 special
                ind_set[(t2,t2,t1)  ] = ntot .+ [1, 3, 2, 2, 5, 4, 4, 6, 8, 7, 7,10, 9, 9,11,11,13,12,12,13] #1 special

                ind_set[(t2,t1,t1)  ] = ntot .+ n_3body_cl_pair .+ [1, 2, 2, 3, 4, 4, 5, 6, 7, 7, 8, 9, 9,10,11,12,12,11,13,13] #3 special
                ind_set[(t1,t2,t1)  ] = ntot .+ n_3body_cl_pair .+ [1, 2, 3, 2, 4, 5, 4, 6, 7, 8, 7, 9,10, 9,12,13,11,11,13,12] #2 special
                ind_set[(t1,t1,t2)  ] = ntot .+ n_3body_cl_pair .+ [1, 3, 2, 2, 5, 4, 4, 6, 8, 7, 7,10, 9, 9,11,11,13,12,12,13] #1 special
                
                            
                ntot += n_3body_cl_pair * 2 #two possible pairs
                            
                
                
            elseif length(v[1]) == 3
                ntot += n_3body_cl_diff
                
                t1,t2,t3 = collect(v[1])[1:3]
                t123 = [t1,t2,t3]
                ind_set[(t1,t2,t3)] = ntot .+ collect(1:20)
                
                PERM = [ [1,2,3],[2,1,3],[3,1,2],[3,2,1],[1,3,2],[2,3,1]] #atom permutations
                for perm in PERM
                    #
                    dist_perm = atom_perm_to_dist_perm(perm)
                    
                    ind_set[ (t123[perm[1]], t123[perm[2]],t123[perm[3]])  ] = Int64[]
                    for p in pattern
                        pnew = p[dist_perm]
                        ind = findfirst(x->x==pnew,pattern)
                        push!(ind_set[ (t123[perm[1]], t123[perm[2]],t123[perm[3]])  ], ind)
                    end
                end
                    
                ntot += n_3body_cl_diff
                    
                    #                for p in pattern
                    #                    
                    #                ind_set[(t2,t1,t3)] = 
                    
            else
                println("something wrong get_fit_cl v[2] == 3  $v")
            end                
        else
                println("something wrong get_fit_cl $v")
        end
        println("ntot $ntot ")
        
    end

    println("at_types $at_types")
    println("ind_set")
    println(ind_set)
    println()
    #return ntot, ind_set
    
        #    for (ind, crys) in enumerate(CRYS)
        #        en = calc_energy_cl(crys, dat_vars=ones(ntot), at_types=at_types, vars_list= vars_list, use_threebody=use_threebody)
        #        println("ind $ind en $en")
        #    end
        
    V = zeros(length(CRYS), ntot)
    for (ind, crys) in enumerate(CRYS)
#        println("ind crys")
#        println(crys)
        begin
            V[ind, :] = ForwardDiff.gradient(x->calc_energy_cl(crys, dat_vars=x, at_types=at_types, ind_set=ind_set,vars_list= vars_list, use_threebody=use_threebody), ones(ntot) ) / crys.nat
            #V[ind, :] = ForwardDiff.gradient(x->efs(crys, x, at_types, ind_set,vars_list, use_threebody), ones(ntot))
        end
    end


    if get_force 
        NAT = 0
        for c in CRYS
            NAT += c.nat
        end
        Vf = zeros(NAT*3, ntot)
        Vs = zeros(length(CRYS) * 6, ntot)
        
        counter_f = 0
        counter_s = 0
        for (ind, crys) in enumerate(CRYS)
#            println("ind crys")
#            println(crys)
            begin
                ret = ForwardDiff.jacobian(x->efs(crys, x, at_types, ind_set,vars_list, use_threebody), ones(ntot) )
                Vf[counter_f .+ (1:crys.nat*3),:] = ret[1:crys.nat*3,:]
                Vs[counter_s .+ (1:6),:] = ret[(crys.nat*3+1):(crys.nat*3+6),:]
                counter_f += 3*crys.nat
                counter_s += 6
            end
        end
    else
        Vf = zeros(0, ntot)
        Vs = zeros(0, ntot)
    end        
        
    
    return V, Vf, Vs, vars, at_types, ind_set
    
end

function efs(crys, dat_vars, at_types, ind_set,vars_list, use_threebody)
    var_type = eltype(dat_vars)
#    println("EFS ----------------------------------------------------------------------------------------")
    #energy = calc_energy_cl(crys, dat_vars=dat_vars, at_types=at_types, ind_set=ind_set,vars_list= vars_list, use_threebody=use_threebody)
    energy, force, stress = energy_force_stress(crys, dat_vars=dat_vars, at_types=at_types, ind_set=ind_set,vars_list= vars_list, use_threebody=use_threebody, var_type=var_type, verbose=false)
    #return energy 
    v = vcat(force[:], [stress[1,1], stress[1,2],stress[1,3],stress[2,2],stress[2,3],stress[3,3]])
    #    return v
    #return energy 
end

function calc_energy_cl(crys::crystal;  database=missing, dat_vars=missing, at_types = missing, vars_list = missing,  DIST=missing, verbose=false, use_threebody=true, ind_set=missing, var_type=missing, turn_off_warn = false)


    At = crys.A'
    
    if verbose
        println()
        println("-----")
        println("Get classical energy from crystal structure")
        println()
    end

    if ismissing(dat_vars) && ismissing(var_type)
        var_type=Float64
    elseif ismissing(var_type)
        var_type = eltype(dat_vars)
    end

#    if verbose  println("LV $var_type")  end
    
    
    

    
    ind2orb, orb2ind, etotal, nval = orbital_index(crys)
    if verbose println("distances") end
     begin
        
        use_dist_arr = true
        if !ismissing(DIST)
            use_dist_arr = false
            R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3 = DIST

        else
            if (use_threebody ) 
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
                            if !turn_off_warn; println("WARNING : structure has 2body distances less than or close to min fitting data distances, may result in errors"); end
                            if !turn_off_warn; println(key," " ,key2, " : ", dmin_types[key], " <~ ", database[key2].min_dist); end
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

    nkeep=size(R_keep)[1]
    ind_arr = zeros(Int64, nkeep, 3)
    ind_arr[:,:] = R_keep[:,2:4]

    if verbose; println("types stuff "); end
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

#        println("types_counter ", types_counter)
        
        for (ind, a) in enumerate(crys.stypes)
            types_arr[ind] = types_dict[a]
        end
#        println("types_arr, ", types_arr)
        
    end

    if verbose; println("cutoff stuff"); end

    
    get_cutoff_pre = Dict()       
    s = Set(crys.stypes)
    for s1 in s
        for s2 in s
            ind1 = types_dict[s1]
            ind2 = types_dict[s2]
            get_cutoff_pre[(ind1, ind2)] = get_cutoff(s1,s2)[:] 
            if use_threebody
                for s3 in s 
                    ind3 = types_dict[s3] 
                    get_cutoff_pre[(ind1, ind2, ind3)] = get_cutoff(s1,s2,s3)[1] 
                end
            end
        end
    end
        
    if false
        cutoff_arr = zeros(crys.nat, crys.nat, 2)
        cutoff_arr3 = zeros(crys.nat, crys.nat, crys.nat)
        for a = 1:crys.nat
            ta = crys.stypes[a]
            for b = 1:crys.nat
                tb = crys.stypes[b]            
                cutoff_arr[a,b,:] = get_cutoff_pre[(ta,tb)][:]
                if use_threebody
                    for c = 1:crys.nat
                        tc = crys.stypes[c]            
                        cutoff_arr3[a,b,c] =  get_cutoff_pre[(ta,tb,tc)]
                    end
                end
            end
        end
    end
    

    badlist = Set()

    if verbose; println("setup database stuff"); end
    #println("ismissing database ", ismissing(database))
     begin
        DAT_IND_ARR = zeros(var_type, types_counter, types_counter,1:n_2body_cl )
        #        DAT_IND_ARR3 = zeros(var_type, types_counter, types_counter, types_counter, 1: max(n_3body_cl_same,n_3body_cl_pair,n_3body_cl_diff) )
                DAT_IND_ARR3 = zeros(var_type, max(n_3body_cl_same,n_3body_cl_pair,n_3body_cl_diff), types_counter, types_counter, types_counter )

        
        if !ismissing(database)
#            println("inside not ismissing database")
            for c1 = 1:types_counter
                for c2 = 1:types_counter
                    t1 = types_dict_reverse[c1]
                    t2 = types_dict_reverse[c2]
                    if (t1,t2) in keys(database)
                        coef = database[(t1,t2)]
                        DAT_IND_ARR[c1,c2,1:n_2body_cl] = coef.datH[:]
                    else
                        if !turn_off_warn; println("WARNING, ",(t1,t2), " database not found"); end
                        within_fit = false
                        push!(badlist, (t1,t2))
                        if !turn_off_warn; println("badlist ", badlist); end
                    end
                    if use_threebody
                        for c3 = 1:types_counter
                            t3 = types_dict_reverse[c3]
                            if (t1,t2,t3) in keys(database)
                                coef = database[(t1,t2,t3)]
                                #                                DAT_IND_ARR3[c1,c2,c3,1:coef.sizeH] = coef.datH[1:coef.sizeH]
                                DAT_IND_ARR3[1:coef.sizeH, c1,c2,c3] = coef.datH[1:coef.sizeH]
                            else
                                if !turn_off_warn;                             println("WARNING, ",(t1,t2,t3), " database not found"); end
                                within_fit = false
                                push!(badlist, (t1,t2,t3))
                                if !turn_off_warn; println("badlist ", badlist); end
                            end
                        end
                    end
                    
                end
            end
        else
            counter = 0
            for v in vars_list
                if v[2] == 2
                    for c1 = 1:types_counter
                        for c2 = 1:types_counter
                            t1 = types_dict_reverse[c1]
                            t2 = types_dict_reverse[c2]
                            if (t1,t2) in keys(ind_set)
                                ind = ind_set[(t1,t2)] 
                                if t1 in v[1] && t2 in v[1] && length(Set([t1,t2])) == length(v[1])
                                #                            println("v $v")
                                #                            println("DAT_IND_ARR $c1 $c2 $n_2body_cl ", dat_vars[ind])
                                #                            println("ind $ind ", typeof(ind))
                                #                            println("dat_vars  $ind ", typeof(dat_vars[ind]))
                                    DAT_IND_ARR[c1,c2,1:n_2body_cl] = dat_vars[ind]
                                #                            println("after")
                                end
                            end
                        end
                    end
                    counter += n_2body_cl
                elseif v[2] == 3
                    for c1 = 1:types_counter
                        for c2 = 1:types_counter
                            for c3 = 1:types_counter
                                t1 = types_dict_reverse[c1]
                                t2 = types_dict_reverse[c2]
                                t3 = types_dict_reverse[c3]
                                if t1 in v[1] && t2 in v[1] && t3 in v[1] && length(Set([t1,t2,t3])) == length(v[1])
                                    if (t1,t2,t3) in keys(ind_set)
                                        ind = ind_set[(t1,t2,t3)] 
                                        #DAT_IND_ARR3[c1,c2,c3,1:length(ind)] = dat_vars[ind]
                                        DAT_IND_ARR3[1:length(ind), c1,c2,c3] = dat_vars[ind]
                                    end
                                end
                            end
                        end
                    end
                    
                else
                    println("something wrong energy_cl v $v ")
                end
                
            end
        end
    end
        #    println("DAT_IND_ARR")
#    println(DAT_IND_ARR)

    warned = false

    nkeep_ab = size(R_keep_ab)[1]
    
    lag_arr_TH = zeros(var_type, n_2body_cl, nthreads())
    energy_2bdy_TH = zeros(var_type,nthreads())

    if verbose println("2body CL") end
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

            if !turn_off_warn && (t1,t2) in badlist
                if warned == false
                    println("WARNING missing $( (t1,t2) ) ")
                    warned = true
                end
                continue
            end
            
            if use_dist_arr
                dist_a = dist_arr[c,1]
                cut_a = dist_arr[c,5]
            else
                dist_a = get_dist_only(a1,a2, R_keep_ab[c,4:6], crys, At)
                #cutoff2, cutoff2on = cutoff_arr[a1,a2,:]
                cutoff2, cutoff2on = get_cutoff_pre[(t1,t2)][:]
                cut_a = cutoff_fn_fast(dist_a, cutoff2 - cutoff_length, cutoff2)
            end
            
            if dist_a <  1e-5    # true onsite
                continue
            end
            
            laguerre_fast!(dist_a, lag_arr)
            energy_2bdy_TH[id] += core_cl(t1, t2, lag_arr, DAT_IND_ARR, var_type) * cut_a 
        end
        
    end
    
    lag_arr3_TH = zeros(var_type, max(n_3body_cl_same, n_3body_cl_pair, n_3body_cl_diff), nthreads())
    energy_3bdy_TH = zeros(var_type,nthreads())

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
            #            println("asdf")
            
            #println("meta_count $meta_count")
            
             @threads for mc in meta_count 
                for counter in mc
                    
                    id = threadid()

                    a1 = array_ind3[counter,1]
                    a2 = array_ind3[counter,2]
                    a3 = array_ind3[counter,3]

            #        println("a1 a2 a3 $a1 $a2 $a3 ", crys.stypes)

                    t1s = crys.stypes[a1]
                    t2s = crys.stypes[a2]
                    t3s = crys.stypes[a3]

                    t1 = types_arr[a1]
                    t2 = types_arr[a2]
                    t3 = types_arr[a3]
                    
                    if !turn_off_warn && (t1,t2,t3) in badlist
                        if warned == false
                            println("WARNING missing $( (t1,t2,t3) ) ")
                            warned = true
                        end
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
                        
                        cut_h = array_floats3[counter,14]
                    else
                        dist12 = get_dist_only(a1,a2, rind1, crys, At)
                        dist13 = get_dist_only(a1,a3, rind2, crys, At)
                        dist23 = get_dist_only(a2,a3, -rind1+rind2, crys, At)

                        cutoff3 = get_cutoff_pre[(t1,t2,t3)]
                            
                        cut_ab = cutoff_fn_fast(dist12, cutoff3 - cutoff_length, cutoff3)
                        cut_ac = cutoff_fn_fast(dist13, cutoff3 - cutoff_length, cutoff3)
                        cut_bc = cutoff_fn_fast(dist23, cutoff3 - cutoff_length, cutoff3)
                        
                        cut_h = cut_ab*cut_ac*cut_bc
                    end
                    
                    
                    #lag_arr3 = lag_arr3_TH[:,id]
                    #ntot = 20
#                    ntot = laguerre_fast_threebdy_cl!(dist12,dist13,dist23, t1,t2,t3, lag_arr3)
                    #                    energy_3bdy_TH[id] += core_3_cl(t1,t2,t3,ntot,lag_arr3,DAT_IND_ARR3, var_type) * cut_h * 100.0

                     energy_3bdy_TH[id] += core_3_cl_laguerre_fast_threebdy_cl!(dist12, dist13, dist23, t1,t2,t3, DAT_IND_ARR3, var_type)  * cut_h * 100.0
                    
                end
            end
        end
    end

    return sum(energy_2bdy_TH) + sum(energy_3bdy_TH)
    
    
end

function ham(x :: Vector, ct, database, donecheck, DIST, FloatX, use_threebody, dat_vars, at_types, vars_list, ind_set, turn_off_warn, verbose)
    T=eltype(x)

    x_r, x_r_strain = reshape_vec(x, ct.nat, strain_mode=true)
    A = FloatX.(ct.A) * (I(3) + x_r_strain)
    crys_dual = makecrys( A , ct.coords + x_r, ct.types, units="Bohr")

    energy = calc_energy_cl(crys_dual;  database=database, DIST=DIST, verbose=verbose, use_threebody=use_threebody, var_type=T, dat_vars=dat_vars, at_types=at_types, vars_list=vars_list, ind_set=ind_set, turn_off_warn=turn_off_warn )
    

end

function energy_force_stress(crys::crystal;  database=missing, verbose=false, use_threebody=true, dat_vars=missing, at_types = missing, vars_list = missing,  DIST=missing, ind_set=missing, var_type=Float64, turn_off_warn = false)

    if verbose; println("dist energy_force_stress");end
     if ismissing(DIST)
        if use_threebody
            R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel_LV(crys,cutoff2X,cutoff3bX,var_type=var_type, return_floats=false)
            DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3
        else
            R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel_LV(crys,cutoff2X,0.0,var_type=var_type, return_floats=false)
            DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3
        end
    else
        R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3 = DIST
    end
    
#    energy = calc_energy_cl(crys, dat_vars=dat_vars, at_types=at_types, ind_set=ind_set,vars_list= vars_list, use_threebody=use_threebody)
#    return energy


#    energy = calc_energy_cl(crys, dat_vars=dat_vars, at_types=at_types, ind_set=ind_set,vars_list= vars_list, use_threebody=use_threebody, DIST=missing)
#    return energy
    #    energy = calc_energy_cl(crys,  database=missing,  DIST=missing, verbose=false,  use_threebody=use_threebody,  dat_vars=dat_vars, at_types=at_types)

    FN_ham = x->ham(x,crys,database, true, DIST, var_type, use_threebody, dat_vars, at_types, vars_list, ind_set, turn_off_warn, verbose)

    if verbose; println("energy test "); end
     energy_test = FN_ham(zeros(var_type, 3*crys.nat + 6))

    #    return energy_test
#    println("energy $energy energy_test $energy_test")
                         
    chunksize=min(15, 3*crys.nat + 6)

    if verbose; println("config");end
     cfg = ForwardDiff.GradientConfig(FN_ham, zeros(var_type, 3*crys.nat + 6), ForwardDiff.Chunk{chunksize}())
#    println("jacham")
#    println("typeof ", typeof(zeros(Float64, 3*crys.nat + 6)))

    if verbose; println("garr"); end
     garr = ForwardDiff.gradient(FN_ham, zeros(var_type, 3*crys.nat + 6) , cfg ) 
    if verbose; println("end garr"); end

    x, stress = reshape_vec(garr, crys.nat)
    f_cart = -1.0 * x
    f_cart = f_cart * inv(crys.A)' 
    stress = -stress / abs(det(crys.A)) 

    #neaten                                                                                                                                                                                                              
    for i = 1:3
	for j = 1:3
            if abs(stress[i,j]) < 1e-9
                stress[i,j] = 0.0
            end
        end
    end
    for i = 1:crys.nat
        for j = 1:3
            if abs(f_cart[i,j]) < 1e-7
                f_cart[i,j] = 0.0
            end
        end
    end

    
    return energy_test, f_cart, stress


end

#function laguerre_fast_cl!(dist, memory)
#
#    a=2.2
#    ad = a*dist
#    expa=exp.(-0.5*ad)
#    memory[1] = 1.0 * expa
#    memory[2] = (1.0 .- ad) .* expa
#    memory[3] = 0.5*(ad.^2 .- 4.0*ad .+ 2) .* expa
#    memory[4] = 1.0/6.0*(-ad.^3 .+ 9.0*ad.^2 .- 18.0*ad .+ 6.0) .* expa
#    memory[5] = 1.0/24.0*(ad.^4 .- 16.0 * ad.^3 .+ 72.0*ad.^2 .- 96.0*ad .+ 24.0) .* expa
#    memory[6] = 1.0/120*(-ad.^5 .+ 25*ad.^4 .- 200 * ad.^3 .+ 600.0*ad.^2 .- 600.0*ad .+ 120.0) .* expa
#
#end


function core_cl(t1, t2, lag_arr, DAT_ARR, var_type )
    energy = zero(var_type)
    for n = 1:n_2body_cl
        energy +=  lag_arr[n]*DAT_ARR[t1,t2,n]
    end
    return energy
end

function core_3_cl(t1, t2, t3,ntot,lag_arr, DAT_ARR3, var_type )
    energy = zero(var_type)
    for n = 1:ntot
        energy +=  lag_arr[n]*DAT_ARR3[t1,t2,t3,n]
    end
    return energy
end

function core_3_cl_laguerre_fast_threebdy_cl!(dist_1, dist_2, dist_3, t1,t2,t3, DAT_ARR3, var_type)

#    println("a")
    begin
        a=2.2

        ad_1 = a*dist_1
        expa_1 =exp.(-0.5*ad_1) 

        L_A_0 = 1.0
        L_A_1 = expa_1
        L_A_2 = (1.0 .- ad_1)*expa_1
        L_A_3 = 0.5*(ad_1.^2 .- 4.0*ad_1 .+ 2.0) * expa_1
        
        ad_2 = a*dist_2
        expa_2 =exp.(-0.5*ad_2) 

        L_B_0 = 1.0
        L_B_1 = expa_2
        L_B_2 = (1.0 .- ad_2)*expa_2
        L_B_3 = 0.5*(ad_2.^2 .- 4.0*ad_2 .+ 2.0) * expa_2
        
        ad_3 = a*dist_3
        expa_3 =exp.(-0.5*ad_3) 

        L_C_0 = 1.0
        L_C_1 = expa_3
        L_C_2 = (1.0 .- ad_3)*expa_3
        L_C_3 = 0.5*(ad_3.^2 .- 4.0*ad_3 .+ 2.0) * expa_3
    end    

#    e3 = expa_1 * expa_2 * expa_3
    
#    if t1==t2 && t1 == t3 #all same
#        memory[1] = L_A_1 * L_B_1 * L_C_1
#        memory[2] = (L_A_2 * L_B_1 * L_C_1 + L_A_1 * L_B_2 * L_C_1 + L_A_1 * L_B_1 * L_C_2)
#        memory[3] = (L_A_3 * L_B_1 * L_C_1 + L_A_1 * L_B_3 * L_C_1 + L_A_1 * L_B_1 * L_C_3)
#        memory[4] = (L_A_2 * L_B_2 * L_C_2)
#        memory[5] = (L_A_1 * L_B_2 * L_C_2 + L_A_2 * L_B_1 * L_C_2 + L_A_2 * L_B_2 * L_C_1)
#        memory[6] = (L_A_0 * L_B_1 * L_C_1 + L_A_1 * L_B_0 * L_C_1 + L_A_1 * L_B_1 * L_C_0)
#        memory[7] = (L_A_0 * L_B_1 * L_C_2 +
#                     L_A_0 * L_B_2 * L_C_1 +
#                     L_A_2 * L_B_0 * L_C_1 +
#                     L_A_1 * L_B_0 * L_C_2 +
#                     L_A_1 * L_B_2 * L_C_0 +
#                     L_A_2 * L_B_1 * L_C_0)
#        return 7
#
#        elseif t1 == t2
#        
#    elseif t1 != t2 && t1 != t3 && t1 != t2

    #println("b")
    @inbounds @fastmath begin
        energy = zero(var_type)

#=        energy += DAT_ARR3[1,t1,t2,t3] * L_A_1 * L_B_1 * L_C_1 #1       #a

        energy += DAT_ARR3[2,t1,t2,t3] * L_A_2 * L_B_1 * L_C_1 #2-4  #b
        energy += DAT_ARR3[3,t1,t2,t3] * L_A_1 * L_B_2 * L_C_1       #b
        energy += DAT_ARR3[4,t1,t2,t3] * L_A_1 * L_B_1 * L_C_2       #c

        energy += DAT_ARR3[5,t1,t2,t3] * L_A_2 * L_B_2 * L_C_2 #8    #f

        energy += DAT_ARR3[6,t1,t2,t3] * L_A_0 * L_B_1 * L_C_1 #12-14  #j
        energy += DAT_ARR3[7,t1,t2,t3] * L_A_1 * L_B_0 * L_C_1         #j
        energy += DAT_ARR3[8,t1,t2,t3] * L_A_1 * L_B_1 * L_C_0         #i
  =#      
        
        energy += DAT_ARR3[1,t1,t2,t3] * L_A_1 * L_B_1 * L_C_1 #1       #a

        energy += DAT_ARR3[2,t1,t2,t3] * L_A_2 * L_B_1 * L_C_1 #2-4  #b
        energy += DAT_ARR3[3,t1,t2,t3] * L_A_1 * L_B_2 * L_C_1       #b
        energy += DAT_ARR3[4,t1,t2,t3] * L_A_1 * L_B_1 * L_C_2       #c
        
        energy += DAT_ARR3[5,t1,t2,t3] * L_A_3 * L_B_1 * L_C_1 #5-7  #d
        energy += DAT_ARR3[6,t1,t2,t3] * L_A_1 * L_B_3 * L_C_1       #d
        energy += DAT_ARR3[7,t1,t2,t3] * L_A_1 * L_B_1 * L_C_3       #e
        
        energy += DAT_ARR3[8,t1,t2,t3] * L_A_2 * L_B_2 * L_C_2 #8    #f

        energy += DAT_ARR3[9,t1,t2,t3] * L_A_1 * L_B_2 * L_C_2 #9-11  #g
        energy += DAT_ARR3[10,t1,t2,t3] * L_A_2 * L_B_1 * L_C_2      #g
        energy += DAT_ARR3[11,t1,t2,t3] * L_A_2 * L_B_2 * L_C_1      #h

        energy += DAT_ARR3[12,t1,t2,t3] * L_A_0 * L_B_1 * L_C_1 #12-14  #j
        energy += DAT_ARR3[13,t1,t2,t3] * L_A_1 * L_B_0 * L_C_1         #j
        energy += DAT_ARR3[14,t1,t2,t3] * L_A_1 * L_B_1 * L_C_0         #i

        energy += DAT_ARR3[15,t1,t2,t3] * L_A_0 * L_B_1 * L_C_2 #15-20  #k 
        energy += DAT_ARR3[16,t1,t2,t3] * L_A_0 * L_B_2 * L_C_1         #l
        energy += DAT_ARR3[17,t1,t2,t3] * L_A_2 * L_B_0 * L_C_1         #l
        energy += DAT_ARR3[18,t1,t2,t3] * L_A_1 * L_B_0 * L_C_2         #k
        energy += DAT_ARR3[19,t1,t2,t3] * L_A_1 * L_B_2 * L_C_0         #m
        energy += DAT_ARR3[20,t1,t2,t3] * L_A_2 * L_B_1 * L_C_0         #m

        ###
        energy += DAT_ARR3[21,t1,t2,t3] * L_A_0 * L_B_2 * L_C_2         #m
        energy += DAT_ARR3[22,t1,t2,t3] * L_A_2 * L_B_0 * L_C_2         #m
        energy += DAT_ARR3[23,t1,t2,t3] * L_A_2 * L_B_2 * L_C_0         #m

        energy += DAT_ARR3[24,t1,t2,t3] * L_A_0 * L_B_1 * L_C_3 #15-20  #k 
        energy += DAT_ARR3[25,t1,t2,t3] * L_A_0 * L_B_3 * L_C_1         #l
        energy += DAT_ARR3[26,t1,t2,t3] * L_A_3 * L_B_0 * L_C_1         #l
        energy += DAT_ARR3[27,t1,t2,t3] * L_A_1 * L_B_0 * L_C_3         #k
        energy += DAT_ARR3[28,t1,t2,t3] * L_A_1 * L_B_3 * L_C_0         #m
        energy += DAT_ARR3[29,t1,t2,t3] * L_A_3 * L_B_1 * L_C_0         #m
        
    end
    
    return energy
    
#=    memory[1] = L_A_1 * L_B_1 * L_C_1 #1       #a

    memory[2] = L_A_2 * L_B_1 * L_C_1 #2-4  #b
    memory[3] = L_A_1 * L_B_2 * L_C_1       #b
    memory[4] = L_A_1 * L_B_1 * L_C_2       #c
    
    memory[5] = L_A_3 * L_B_1 * L_C_1 #5-7  #d
    memory[6] = L_A_1 * L_B_3 * L_C_1       #d
    memory[7] = L_A_1 * L_B_1 * L_C_3       #e
 
    memory[8] = L_A_2 * L_B_2 * L_C_2 #8    #f

    memory[9] = L_A_1 * L_B_2 * L_C_2 #9-11  #g
    memory[10] = L_A_2 * L_B_1 * L_C_2      #g
    memory[11] = L_A_2 * L_B_2 * L_C_1      #h

    memory[12] = L_A_0 * L_B_1 * L_C_1 #12-14  #j
    memory[13] = L_A_1 * L_B_0 * L_C_1         #j
    memory[14] = L_A_1 * L_B_1 * L_C_0         #i

    memory[15] = L_A_0 * L_B_1 * L_C_2 #15-20  #k 
    memory[16] = L_A_0 * L_B_2 * L_C_1         #l
    memory[17] = L_A_2 * L_B_0 * L_C_1         #l
    memory[18] = L_A_1 * L_B_0 * L_C_2         #k
    memory[19] = L_A_1 * L_B_2 * L_C_0         #m
    memory[20] = L_A_2 * L_B_1 * L_C_0         #m

#    memory[2:20] .= 0.0
    
    return 20
#    end
=#

end


function laguerre_fast_threebdy_cl!(dist_1, dist_2, dist_3, t1,t2,t3, memory)

    a=2.2

    ad_1 = a*dist_1
    expa_1 =exp.(-0.5*ad_1) 

    L_A_0 = 1.0
    L_A_1 = expa_1
    L_A_2 = (1.0 .- ad_1)*expa_1
    L_A_3 = 0.5*(ad_1.^2 .- 4.0*ad_1 .+ 2.0) * expa_1
    
    ad_2 = a*dist_2
    expa_2 =exp.(-0.5*ad_2) 

    L_B_0 = 1.0
    L_B_1 = expa_2
    L_B_2 = (1.0 .- ad_2)*expa_2
    L_B_3 = 0.5*(ad_2.^2 .- 4.0*ad_2 .+ 2.0) * expa_2
    
    ad_3 = a*dist_3
    expa_3 =exp.(-0.5*ad_3) 

    L_C_0 = 1.0
    L_C_1 = expa_3
    L_C_2 = (1.0 .- ad_3)*expa_3
    L_C_3 = 0.5*(ad_3.^2 .- 4.0*ad_3 .+ 2.0) * expa_3
    

#    e3 = expa_1 * expa_2 * expa_3
    
#    if t1==t2 && t1 == t3 #all same
#        memory[1] = L_A_1 * L_B_1 * L_C_1
#        memory[2] = (L_A_2 * L_B_1 * L_C_1 + L_A_1 * L_B_2 * L_C_1 + L_A_1 * L_B_1 * L_C_2)
#        memory[3] = (L_A_3 * L_B_1 * L_C_1 + L_A_1 * L_B_3 * L_C_1 + L_A_1 * L_B_1 * L_C_3)
#        memory[4] = (L_A_2 * L_B_2 * L_C_2)
#        memory[5] = (L_A_1 * L_B_2 * L_C_2 + L_A_2 * L_B_1 * L_C_2 + L_A_2 * L_B_2 * L_C_1)
#        memory[6] = (L_A_0 * L_B_1 * L_C_1 + L_A_1 * L_B_0 * L_C_1 + L_A_1 * L_B_1 * L_C_0)
#        memory[7] = (L_A_0 * L_B_1 * L_C_2 +
#                     L_A_0 * L_B_2 * L_C_1 +
#                     L_A_2 * L_B_0 * L_C_1 +
#                     L_A_1 * L_B_0 * L_C_2 +
#                     L_A_1 * L_B_2 * L_C_0 +
#                     L_A_2 * L_B_1 * L_C_0)
#        return 7
#
#        elseif t1 == t2
#        
#    elseif t1 != t2 && t1 != t3 && t1 != t2

    memory[1] = L_A_1 * L_B_1 * L_C_1 #1       #a

    memory[2] = L_A_2 * L_B_1 * L_C_1 #2-4  #b
    memory[3] = L_A_1 * L_B_2 * L_C_1       #b
    memory[4] = L_A_1 * L_B_1 * L_C_2       #c
    
    memory[5] = L_A_3 * L_B_1 * L_C_1 #5-7  #d
    memory[6] = L_A_1 * L_B_3 * L_C_1       #d
    memory[7] = L_A_1 * L_B_1 * L_C_3       #e
 
    memory[8] = L_A_2 * L_B_2 * L_C_2 #8    #f

    memory[9] = L_A_1 * L_B_2 * L_C_2 #9-11  #g
    memory[10] = L_A_2 * L_B_1 * L_C_2      #g
    memory[11] = L_A_2 * L_B_2 * L_C_1      #h

    memory[12] = L_A_0 * L_B_1 * L_C_1 #12-14  #j
    memory[13] = L_A_1 * L_B_0 * L_C_1         #j
    memory[14] = L_A_1 * L_B_1 * L_C_0         #i

    memory[15] = L_A_0 * L_B_1 * L_C_2 #15-20  #k 
    memory[16] = L_A_0 * L_B_2 * L_C_1         #l
    memory[17] = L_A_2 * L_B_0 * L_C_1         #l
    memory[18] = L_A_1 * L_B_0 * L_C_2         #k
    memory[19] = L_A_1 * L_B_2 * L_C_0         #m
    memory[20] = L_A_2 * L_B_1 * L_C_0         #m

#    memory[2:20] .= 0.0
    
    return 20
#    end


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




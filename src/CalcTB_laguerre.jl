###include("Crystal.jl")
###include("DFToutMod.jl")
#using XMLDict



"""
    module CalcTB

Create TB matrix sets from coefficients, or prepare to fit the coefficients.
"""
module CalcTB


using LinearAlgebra
using SpecialFunctions
using GZip
using EzXML
using XMLDict
using Printf
using SparseArrays
using LoopVectorization



#using Statistics

using Base.Threads
#using Base.lock

using ..AtomicMod:atom
using ..CrystalMod:crystal
#using ..TB:orbital_index
using ..CrystalMod:orbital_index

using ..CrystalMod:distances_etc
using ..CrystalMod:distances_etc_3bdy_parallel_old
using ..CrystalMod:distances_etc_3bdy_parallel
using ..CrystalMod:distances_etc_3bdy_parallel2
using ..CrystalMod:distances_etc_3bdy_parallel_LV
using ..CrystalMod:get_dist

using ..TB:make_tb_crys
using ..TB:make_tb
using ..TB:make_kgrid
using ..TB:orb_num
using ..TB:Hk
using ..TB:summarize_orb
using ..TB:summarize_orb_num
using ..TB:calc_energy
using ..CrystalMod:get_grid
using ..TB:myfft
using ..TB:trim
using ..TB:tb_crys
using ..Utility:cutoff_fn
using ..Utility:cutoff_fn_fast
using ..Utility:arr2str
using ..Utility:str_w_spaces
using ..Utility:parse_str_ARR_float
using ..TB:make_rdict

using ..Utility:dict2str
using ..Utility:str2tuplesdict
using ForwardDiff
using ..CrystalMod:makecrys

using Random

#include("Coef_format_convert.jl")

const one = [1.0]

#using PyPlot
using Plots


#include("Atomdata.jl")
using ..Atomdata:atoms
using ..Atomdata:get_cutoff


#constants

const sqrt3 = sqrt(3)
const sqrt3d2 = sqrt(3)/2.0

#this set up the number of terms in parts of the model

#const n_2body = 7
#const n_2body_onsite = 6
#const n_2body_S = 7

const n_2body = 6
const n_2body_onsite = 5
const n_2body_S = 6

const n_3body = 8
const n_3body_same = 6
const n_3body_triple = 4

const n_3body_onsite = 2
#const n_3body_onsite_same = 4
const n_3body_onsite_same = 5

EXP_a = [2.5]

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
mutable struct coefs

    dim::Int64
    datH::Array{Float64,1}
    datS::Array{Float64,1}
    sizeH::Int64
    sizeS::Int64
    #    inds::Dict{Tuple, Array{Int64,1}}
    inds::Dict{Array{Symbol}, Array{Int64,1}}
    inds_int::Dict{Array{Symbol}, Any}

#    names::Array{String,1}
    names::Set
    orbs::Array{Any,1}
    cutoff::Float64
    min_dist::Float64

#    maxmin_val_train::Dict
    dist_frontier::Dict
    version::Int64
    lim::Dict
    repval::Dict
    
end

function construct_coef_string(co)
    inds = co.inds
    st = String[]
    for sym = [:H, :O, :S]
        if sym == :H
            push!(st, "# Hopping Hamiltonian index\n")
        elseif sym == :O
            push!(st, "# Onsite Hamiltonian index\n")
        elseif sym == :S
            push!(st, "# Overlap (S) index\n")
        end
        for k in keys(inds)
            if k[end] == sym
                push!(st, "$k  ;  "*arr2str(inds[k])*" ; " *arr2str(co.datH[inds[k]])*"\n")
            end
        end
    end
    str = ""
    for line in st
        str *= line
    end

    return str

end


"""
    function write_coefs(filename, co::coefs; compress=true)

Write `coefs` to a file. Compress uses gzip. See `read_coefs`
"""
function write_coefs(filename, co::coefs; compress=true)
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
    addelement!(c, "datS", arr2str(co.datS))
    addelement!(c, "sizeH", string(co.sizeH))
    addelement!(c, "sizeS", string(co.sizeS))

    indsstr = construct_coef_string(co)

    addelement!(c, "inds", indsstr)        #string(co.inds))
    addelement!(c, "names", str_w_spaces(co.names))
#    addelement!(c, "orbs", str_w_spaces(co.orbs))
    addelement!(c, "cutoff", string(co.cutoff))
    addelement!(c, "min_dist", string(co.min_dist))
#    addelement!(c, "maxmin_val_train", string(co.maxmin_val_train))
#    addelement!(c, "dist_frontier", string(co.dist_frontier))

#    println(dict2str(co.maxmin_val_train))
#    println("type ", typeof(co.maxmin_val_train))
    
#    addelement!(c, "maxmin_val_train", dict2str(co.maxmin_val_train))
    if !isempty(co.dist_frontier)
        addelement!(c, "dist_frontier", dict2str(co.dist_frontier))
    end
    if !isempty(co.lim)
        addelement!(c, "lim", dict2str(co.lim))
    end
    if !isempty(co.repval)
        addelement!(c, "repval", dict2str(co.repval))
    end
    

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
#    try
#        println("asdfasdfasdf $filename")
        f = gzopen(filename, "r")
#    catch
#        println("error opening $filename")
#    end

    fs = read(f, String)
    close(f)

    d = xml_dict(fs)["root"]
    
    dim = parse(Int64, (d["coefs"]["dim"]))
    sizeH = parse(Int64, d["coefs"]["sizeH"])
    sizeS = parse(Int64, d["coefs"]["sizeS"])

    if sizeH > 0
        datH = parse_str_ARR_float(d["coefs"]["datH"])
    else
        datH = Float64[]
    end
    if sizeS > 0
        datS = parse_str_ARR_float(d["coefs"]["datS"])
    else
        datS = Float64[]
    end
    
#    addelement!(c, "inds", string(co.inds))
    
    names = Set(String.(split(d["coefs"]["names"])))
#    orbs = Symbol.(split(d["coefs"]["orbs"]))
    cutoff = parse(Float64, d["coefs"]["cutoff"])
    min_dist = parse(Float64, d["coefs"]["min_dist"])

#    println(d["coefs"]["maxmin_val_train"])
    
#    maxmin_val_train = str2tuplesdict(d["coefs"]["maxmin_val_train"])

#    println("dist")
#    println(d["coefs"]["dist_frontier"])
#    println(eval(d["coefs"]["dist_frontier"]))
#    println()
    if haskey(d["coefs"], "dist_frontier")
        dist_frontier = str2tuplesdict(eval(d["coefs"]["dist_frontier"]))
    else
        dist_frontier = missing
    end
    
    version = 1
    if "version" in keys(d["coefs"])
        version = parse(Int64, d["coefs"]["version"])
    end
#    println("version $version")

    if haskey(d["coefs"], "lim")
        lim = str2tuplesdict(eval(d["coefs"]["lim"]))
    else
        lim = missing
    end
    if haskey(d["coefs"], "repval")
        repval = str2tuplesdict(eval(d["coefs"]["repval"]))
    else
        repval = missing
    end
        
    co = make_coefs(names,dim, datH=datH, datS=datS, min_dist=min_dist, dist_frontier = dist_frontier, version=version, lim=lim, repval=repval)

    return co
    
end



"""
    function make_coefs(at_list, dim; datH=missing, datS=missing, cutoff=18.01, min_dist = 3.0, fillzeros=false, maxmin_val_train=missing, dist_frontier=missing)

Constructor for `coefs`. Can create coefs filled with ones for testing purposes.

See `coefs` to understand arguments.
"""
function make_coefs(at_list, dim; datH=missing, datS=missing, cutoff=18.01, min_dist = 3.0, fillzeros=false, dist_frontier=missing, version=3, lim=missing, repval=missing)

#    println("make coefs")
#    sort!(at_list)

    if version == 2 || version == 3
        totH,totS, data_info, orbs = get_data_info_v2(at_list, dim)
    elseif version == 1
        totH,totS, data_info, orbs = get_data_info_v1(at_list, dim)
    else
        println("warning, bad version $version")
        version = 2
        totH,totS, data_info, orbs = get_data_info_v2(at_list, dim)
    end

#    println("make_coefs ", at_list)
#    for key in keys(data_info)
#        println(key, " data_info ", minimum(data_info[key]), " " ,maximum(data_info[key]))
#    end

#    println("make coeffs $totH $totS")

#    if ismissing(maxmin_val_train)
#        maxmin_val_train = Dict()
#        for d in keys(data_info)
#            maxmin_val_train[d] = (1e7, -1e7)
#        end
#    end

    if ismissing(datH)
        if fillzeros
            datH = zeros(totH)
        else
            datH = ones(totH)
        end
    end
    if ismissing(datS)
        if fillzeros
            datS = zeros(totS)
        else
            datS = ones(totS) * 0.1
        end
    end

    #    if length(datH) < totH
    #        println("dat convert ", at_list, dim)
    #        datH = fix_format_change(datH, totH, dim, at_list, data_info)
    #    end
    
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


    inds_int = Dict()
    ninds_int = Dict()

    
    
    if dim == 2
        ninds_int = zeros(UInt16, 4,4)

        for k in keys(data_info)

#            println("k $k")

            if length(k) == 5
                t1,s1,t2,s2,HH = k
                if !haskey(inds_int, [t1,t2])
                    inds_int[[t1,t2]] = (zeros(UInt16, 4,4,32),zeros(UInt16, 4,4,32), zeros(UInt16, 4,4), zeros(UInt16, 4,4))
                end

                if HH == :H
                    i1 = summarize_orb_num(s1) + 1
                    i2 = summarize_orb_num(s2) + 1
                    n = length(data_info[[t1,s1,t2,s2,HH]])
                    inds_int[[t1,t2]][3][i1,i2] = n
                    inds_int[[t1,t2]][1][i1,i2,1:n] = data_info[[t1,s1,t2,s2,HH]]
                end
                if HH == :S
                    i1 = summarize_orb_num(s1) + 1
                    i2 = summarize_orb_num(s2) + 1
                    n = length(data_info[[t1,s1,t2,s2,HH]])
                    inds_int[[t1,t2]][4][i1,i2] = n
                    inds_int[[t1,t2]][2][i1,i2,1:n] = data_info[[t1,s1,t2,s2,HH]]
                end
            elseif length(k) == 4
                t1,s1,s2,HH = k
                if !haskey(inds_int, [t1])
                    inds_int[[t1]] = (zeros(UInt16, 4,4,32),zeros(UInt16, 4,4))
                end
                if HH == :O
                    i1 = summarize_orb_num(s1) + 1
                    i2 = summarize_orb_num(s2) + 1
                    n = length(data_info[[t1,s1,s2,HH]])
                    inds_int[[t1]][2][i1,i2] = n
                    inds_int[[t1]][1][i1,i2,1:n] = data_info[[t1,s1,s2,HH]]
                end
            end
                
        end
    end


    if dim == 3
        ninds_int = zeros(UInt16, 4,4)

        for k in keys(data_info)

            if length(k) == 6

                t1,s1,t2,s2,t3,HH = k

                if !haskey(inds_int, [t1,t2,t3])
                    inds_int[[t1,t2,t3]] = (zeros(UInt16, 4,4,32), zeros(UInt16, 4,4))
                end

                if HH == :H
                    i1 = summarize_orb_num(s1) + 1
                    i2 = summarize_orb_num(s2) + 1
                    n = length(data_info[[t1,s1,t2,s2,t3,HH]])
                    inds_int[[t1,t2,t3]][2][i1,i2] = n
                    inds_int[[t1,t2,t3]][1][i1,i2,1:n] = data_info[[t1,s1,t2,s2,t3,HH]]
                end
            elseif length(k) == 5
                t1,t2,t3,s1,HH = k

                if !haskey(inds_int, [t1,t2,t3,:O])
                    inds_int[[t1,t2,t3,:O]] = (zeros(UInt16, 4,32), zeros(UInt16, 4))
                end

                if HH == :O
                    i1 = summarize_orb_num(s1) + 1
                    n = length(data_info[[t1,t2,t3,s1,HH]])

                    inds_int[[t1,t2,t3,:O]][2][i1] = n
                    inds_int[[t1,t2,t3,:O]][1][i1,1:n] = data_info[[t1,t2,t3,s1,HH]]
                    
                end
            end                
                

        end
    end

    if ismissing(lim)
        lim = Dict()
    end
    if ismissing(repval)
        repval = Dict()
    end
    
    if dim == 2
        at_arr = Symbol.([i for i in at_list])
        if length(at_list) == 1
            at_arr = [(at_arr[1], at_arr[1])]
        else
            at_arr = [(at_arr[1], at_arr[2]), (at_arr[2], at_arr[1])]
        end

        for a in at_arr
            #println("a $a")
#            println(typeof(a))
            if !haskey(lim, a)
                lim[a] = 0.02
            end
            if !haskey(repval, a)
                repval[a] = 0.01
            end
        end
    elseif dim == 3
        at_arr = Symbol.([i for i in at_list])
        if length(at_list) == 1
            at_arr = [(at_arr[1], at_arr[1], at_arr[1])]
        elseif length(at_list) == 2
            at_arr = [(at_arr[1], at_arr[1], at_arr[2]),
                      (at_arr[1], at_arr[2], at_arr[1]),
                      (at_arr[2], at_arr[1], at_arr[1]),
                      (at_arr[1], at_arr[2], at_arr[2]),
                      (at_arr[2], at_arr[1], at_arr[2]),
                      (at_arr[2], at_arr[2], at_arr[1])]
        elseif length(at_list) == 3
            at_arr = [(at_arr[1], at_arr[2], at_arr[3]), 
                      (at_arr[1], at_arr[3], at_arr[2]),
                      (at_arr[2], at_arr[1], at_arr[3]),
                      (at_arr[2], at_arr[3], at_arr[1]),
                      (at_arr[3], at_arr[1], at_arr[2]),
                      (at_arr[3], at_arr[2], at_arr[1]),
                      ]
        end

        for a in at_arr
            if !haskey(lim,a)
                #println("add $a")
                lim[a] = 0.04
            end
            if !haskey(repval,a)
                repval[a] = 0.001
            end
        end
    end

#    for x in [dim, datH, datS, totH, totS, data_info, inds_int, at_list, orbs, cutoff, min_dist, dist_frontier2, version, lim, repval]
#        println(typeof(x))
#        end
    
    return coefs(dim, datH, datS, totH, totS, data_info, inds_int, at_list, orbs, cutoff, min_dist, dist_frontier2, version, lim, repval)

end
    

"""
    function plot_database(database, entry, t=missing)

Plot some data from `coefs`. Needs to be updated to work with Plots, probably doesn't work now.
"""
function plot_database(database, entry, t=missing)#::coefs)
#    title("dat.names")

    dat = database[entry]
    function sym_tup_2_str(tup)
        st = ""
        c=0
        for s in tup
            if c >= 1
                st = st *"_"*String(s)
            else
                st = String(s)
            end
            c+=1
        end
        return st
    end

    dist = collect(4.0:0.1:21)

    cut2 = cutoff_fn.(dist, cutoff2X-cutoff_length, cutoff2X)
    cut3 = cutoff_fn.(dist, cutoff3bX-cutoff_length, cutoff3bX)
    cut_on = cutoff_fn.(dist, cutoff_onX-cutoff_length, cutoff_onX)

#    plot([dist[1]-5, dist[end]+5], zeros(2,1), "--", color="0.5", LineWidth=0.5)
    xlims!(dist[1]-0.5, dist[end]+0.5)
#    styles = ["-b", "-g", "-r", "-c", "-m", "-k", "-y","--b", "--g", "--r", "--c", "--m", "--k", "--y", ":b", ":g", ":r", ":c", ":m", ":k", ":y"]
#    ns = length(styles)
    linew = [2.5, 2.5,  2.5,  2.5, 2.0, 2.0,	2.0,  2.0, 1.5, 1.5,	1.5,  1.5, 1.0, 1.0,	1.0,  1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    
    cs=1
    handels = String[]
    plts = []

    println(dat.names)

    plot()

    if length(entry) == 2 #twobody
        for key in keys(dat.inds)
            if !ismissing(t)
                if key[end] != t
                    continue
                end
            end
            println("key ", key)

            ind = dat.inds[key]
            legendS = sym_tup_2_str(key)
            if key[end] == :H 

                d = dat.datH
                a = two_body_H(dist, d[ind[1:n_2body]])

 #               s = styles[cs%ns + 1]
                cs += 1

                plot!(dist, (a.*cut2),  label=legendS, LineWidth=linew[cs])


                if key[2] == :p && key[4] == :p
                    b = two_body_H(dist, d[ind[1+n_2body:n_2body*2]])
#                    s = styles[cs%ns + 1]
                    cs += 1
                    plot!(dist, b.*cut2,  label=legendS*"_π", LineWidth=linew[cs])

                end
            end

            if key[end] == :S 

                d = dat.datS
                a = two_body_H(dist, d[ind[1:n_2body_S]])

#                s = styles[cs%ns + 1]
                cs += 1
                plot!(dist, a.*cut2,  label=legendS, LineWidth=linew[cs])

                if key[2] == :p && key[4] == :p
                    b = two_body_H(dist, d[ind[1+n_2body_S:n_2body_S*2]])
#                    s = styles[cs%ns + 1]
                    cs += 1
                    plot!(dist, b.*cut2,  label=legendS*"_π", LineWidth=linew[cs])

                end
            end

            if key[3] == :O 
#            if false

                d = dat.datH
                if key[1] == key[2]

                    a = two_body_O(dist, d[ind[1:n_2body_onsite]])
#                    s = styles[cs%ns + 1]
                    cs += 1
                    plot!(dist, a.*cut_on,  label=legendS, LineWidth=linew[cs])
                end
                if (key[1] == :s && key[2] == :p) ||    (key[2] == :s && key[1] == :p)
                    b = two_body_O(dist, d[1:n_2body_onsite])
#                    s = styles[cs%ns + 1]
                    cs += 1
                    plot!(dist, b.*cut_on,  label=legendS*"_π", LineWidth=linew[cs])
                    
                elseif (key[1] == :p && key[2] == :p) 
                    b = two_body_O(dist, d[ind[1+n_2body_onsite:n_2body_onsite*2]])
#                    s = styles[cs%ns + 1]
                    cs += 1
                    plot!(dist, b.*cut_on,  label=legendS*"_π", LineWidth=linew[cs])
                end

            end

        end
    else
        for key in keys(dat.inds)
            if !ismissing(t)
                if key[end] != t
                    continue
                end
            end
            println("key ", key)
            ind = dat.inds[key]
            legendS = sym_tup_2_str(key)
            if key[end] == :H 

                d = dat.datH
                o = ones(size(dist))

                sameat = (key[1] == key[3])

                a = three_body_H(6.0*o, dist, dist,sameat,  d[ind]) .*cut3
                b = three_body_H(6.0*o, dist, 6.0*o, sameat, d[ind]) .*cut3
                c = three_body_H( dist, dist, dist, sameat,d[ind]) .*cut3 #.*cut3 .*cut3

#                print([size(a), size(b), size(c)])

#                s = styles[cs%ns + 1]
#                cs += 1
                plot!(dist, a,  label=legendS*"_1")
#                s = styles[cs%ns + 1]
#                cs += 1
#                plot!(dist, b,  label=legendS*"_2")
#                s = styles[cs%ns + 1]
#                cs += 1
#                plot!(dist, c,  label=legendS*"_3")

            end
            if key[end] == :O

                d = dat.datH
                o = ones(size(dist))
                
#                sameat = (dat.names[2] == dat.names[3])
#                sameat = (entry[2] == entry[3])
                sameat = (key[3] == key[4])


#                a = three_body_O(dist, 5.0*o,5.0*o,sameat, d[ind]) .*cut3
#                b = three_body_O(5.0*o, dist,5.0*o,sameat, d[ind]) .*cut3
#                c = three_body_O(dist, dist, dist,sameat, d[ind]) .*cut3 .*cut3 .*cut3

#                print([size(a), size(b), size(c)])

#                s = styles[cs%ns + 1]
#                cs += 1
#                plot(dist, a, s, label=legendS*"_1")
#                s = styles[cs%ns + 1]
#                cs += 1
#                plot(dist, b, s, label=legendS*"_2")
#                s = styles[cs%ns + 1]
#                cs += 1
#                plot!(dist, c, label=legendS*"_3")

            end

        end
    end
#    println(handels)
#    legend(plts, handels)

    plot!(legend=true, box=true)
end

"""
    function get_data_info(at_set, dim)

Figure out the arrangement of data in a `coefs` file.

Loops over various combinations of orbitals and atoms and assigns them places in `datH` and `datS`, depending on the terms included in the model and the dimensionaly.
"""
function get_data_info_v2(at_set, dim)
    
    
    data_info = Dict{Any, Array{Int64,1}}()
    orbs = []
    if dim == 0
        at_list = Symbol.(collect(at_set))
        data_info[[:eam, at_list[1]]] = [1,2,3]
        return 3 ,0, data_info, orbs

        #return totHO ,totS, data_info, orbs
        
    elseif dim == 2 #2body
        
        at_list = Symbol.([i for i in at_set])
        #        println(at_list)
        if length(at_list) == 1
            at_list = [at_list[1], at_list[1]]
        end
        sort!(at_list)
        #        println(at_list)
        
        orbs1 = atoms[at_list[1]].orbitals
        orbs2 = atoms[at_list[2]].orbitals

        at1 = at_list[1]
        at2 = at_list[2]

        if at1 == at2
            same_at = true
        else
            same_at = false
        end

        #        orbs = []

        #2body part
        function get2bdy(n, symb)
            tot=0
            for o1 in orbs1
                for o2 in orbs2
                    if same_at && ((o2 == :s && o1 == :p) || (o2 == :s && o1 == :d) || (o2 == :p && o1 == :d))
                        continue
                    end
                    #                    push!(orbs, (o1, o2, symb))

                    if [at1, o1, at2, o2, symb] in keys(data_info)
                        continue
                    end

                    if o1 == :s && o2 == :s
                        data_info[[at1, o1, at2, o2, symb]] = tot+1:tot+n
                        data_info[[at2, o2, at1, o1, symb]] = tot+1:tot+n
                        tot += n
                        
                    elseif (o1 == :s && o2 == :p ) || (o1 == :p && o2 == :s )
                        data_info[[at1, o1, at2, o2, symb]] = tot+1:tot+n
                        data_info[[at2, o2, at1, o1, symb]] = tot+1:tot+n
                        tot += n
                        #                        if same_at
                        #                            data_info[[o2, o1, symb]] = data_info[(o1, o2, symb)]
                        #                        end
                        
                    elseif (o1 == :p && o2 == :p )
                        data_info[[at1, o1, at2, o2, symb]] = tot+1:tot+n*2
                        data_info[[at2, o2, at1, o1, symb]] = tot+1:tot+n*2
                        tot += n*2

                        #                    elseif (o1 == :p && o2 == :p )
                        #                        data_info[[at1, o1, at2, o2, symb]] = tot+1:tot+n*2
                        #                        data_info[[at2, o2, at1, o1, symb]] = tot+1:tot+n*2
                        #                        tot += n*2

                    elseif (o1 == :s && o2 == :d ) || (o1 == :d && o2 == :s )
                        data_info[[at1, o1, at2, o2, symb]] = tot+1:tot+n
                        data_info[[at2, o2, at1, o1, symb]] = tot+1:tot+n
                        tot += n

                    elseif (o1 == :p && o2 == :d ) || (o1 == :d && o2 == :p )
                        data_info[[at1, o1, at2, o2, symb]] = tot+1:tot+n*2
                        data_info[[at2, o2, at1, o1, symb]] = tot+1:tot+n*2
                        tot += n*2

                    elseif (o1 == :d && o2 == :d ) 
                        data_info[[at1, o1, at2, o2, symb]] = tot+1:tot+n*3
                        data_info[[at2, o2, at1, o1, symb]] = tot+1:tot+n*3
                        tot += n*3


                        
                    end
                end
            end
            return tot
        end

        totH = get2bdy(n_2body, :H)
        totS = get2bdy(n_2body_S, :S)
        #        println("totH $totH totS $totS")
        #onsite part
        function getonsite(atX,orbsX, tot, n)
            for o1 in orbsX
                for o2 in orbsX
                    if (o2 == :s && o1 == :p) || (o2 == :s && o1 == :d) || (o2 == :p && o1 == :d)
                        continue
                    end

                    if [atX, o1, o2, :O] in keys(data_info)
                        continue
                    end


                    #                    push!(orbs, (o1, o2, :O))
                    if o1 == :s && o2 == :s
                        data_info[[atX, o1, o2, :O]] = tot+1:tot+n
                        #                        println("data_info"[, (atX, o1, o2, :O], tot+1:tot+n)
                        tot += n
                    elseif (o1 == :s && o2 == :p )
                        data_info[[atX, o1, o2, :O]] = tot+1:tot+n
                        data_info[[atX, o2, o1, :O]] = data_info[[atX, o1, o2, :O]]
                        tot += n
                    elseif (o1 == :p && o2 == :p )
                        data_info[[atX, o1, o2, :O]] = tot+1:tot+n*2
                        tot += n*2

                    elseif o1 == :s && o2 == :d
                        data_info[[atX, o1, o2, :O]] = tot+1:tot+n
                        data_info[[atX, o2, o1, :O]] = data_info[[atX, o1, o2, :O]]
                        tot += n

                    elseif o1 == :p && o2 == :d
                        data_info[[atX, o1, o2, :O]] = tot+1:tot+n
                        data_info[[atX, o2, o1, :O]] = data_info[[atX, o1, o2, :O]]
                        tot += n

                    elseif o1 == :d && o2 == :d
                        data_info[[atX, o1, o2, :O]] = tot+1:tot+n*2
                        data_info[[atX, o2, o1, :O]] = data_info[[atX, o1, o2, :O]]
                        tot += n*2
                        
                    end
                end
            end
            return tot
        end


        #PUREONSITE
        #        if same_at #true onsite terms
        #           for o in orbs1
        ##                println("true onsite ", o)
        #                data_info[[at1, o, :A]] = [totH+1]
        #                totH += 1
        #            end
        #        end


        totHO = getonsite(at1, orbs1, totH, n_2body_onsite)

        if !(same_at) #need reverse if not same atom
            totHO = getonsite(at2, orbs2, totHO, n_2body_onsite)
        end

    elseif dim == 3 #3body
        
        totS = 0 #no 3body overlap terms

        at_list = Symbol.([i for i in at_set])
        sort!(at_list)
        
        if length(at_list) == 1


            #permutations are trivial
            perm_ij = [[at_list[1], at_list[1], at_list[1]]]
            perm_on = [[at_list[1], at_list[1], at_list[1]]]
        elseif length(at_list) == 2


            #unique permutations
            perm_ij = [[at_list[1], at_list[1], at_list[2]] ,
                       [at_list[2], at_list[2], at_list[1]] ,
                       [at_list[1], at_list[2], at_list[1]] ,
                       [at_list[1], at_list[2], at_list[2]] ]

            perm_on = [[at_list[1], at_list[2], at_list[2]] ,
                       [at_list[1], at_list[1], at_list[2]] ,
                       [at_list[2], at_list[1], at_list[1]] ,
                       [at_list[2], at_list[1], at_list[2]] ]
            
        elseif length(at_list) == 3


            #all permutations exist hij
            perm_ij = [[at_list[1], at_list[2], at_list[3]] ,
                       [at_list[1], at_list[3], at_list[2]] ,
                       [at_list[2], at_list[1], at_list[3]] ,
                       [at_list[2], at_list[3], at_list[1]] ,
                       [at_list[3], at_list[1], at_list[2]] ,
                       [at_list[3], at_list[2], at_list[1]] ]

            #onsite can flip last 2 atoms
            perm_on = [[at_list[1], at_list[2], at_list[3]] ,
                       [at_list[2], at_list[1], at_list[3]] ,
                       [at_list[3], at_list[1], at_list[2]]]

            
        else
            println("ERROR  get_data_info dim $at_set $at_list")
        end
        



        function get3bdy(n, symb, start, at1, at2, at3)
            tot=start

            orbs1 = atoms[at1].orbitals
            orbs2 = atoms[at2].orbitals

            if at1 == at2
                same_at = true
            else
                same_at = false
            end
            if at1 != at2 && at1 != at3 && at2 != at3
                triple = true
            else
                triple = false
            end
            
            for o1 in orbs1
                for o2 in orbs2


                    
                    if same_at && ((o2 == :s && o1 == :p) || (o2 == :s && o1 == :d) || (o2 == :p && o1 == :d))
                        continue
                    end
                    
                    #                    push!(orbs, (o1, o2, symb))

                    if [at1, o1, at2, o2, at3,  symb] in keys(data_info)
                        continue
                    end

                    if triple
                        data_info[[at1, o1, at2, o2, at3,  symb]] = collect(tot+1:tot+n_3body_triple)
                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1, 3, 2, 4]
                        tot += n_3body_triple
                        continue
                    end

                    if same_at
                        data_info[[at1, o1, at2, o2, at3,  symb]] = collect(tot+1:tot+n)
                        data_info[[at2, o2, at1, o1, at3,  symb]] = collect(tot+1:tot+n)
                    else                                                  #[1 2 3 4 5 6 7 8  9  10 11 12 13 14 15 16 17 18]
                        data_info[[at1, o1, at2, o2, at3,  symb]] = collect(tot+1:tot+n)
                        #                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1 4 6 2 5 3 7 10 12 8  11 9  13 16 18 14 17 15]' #switch 2 4 and 3 6
                        #                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1 4 6 2 5 3 7 10 12 8  11 9  ]' #switch 2 4 and 3 6
                        #                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1 3 2 4 5 7 6 8 9]' #switch 2 4 and 3 6
                        #                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1 3 2 4 5 7 6]' #switch 2 4 and 3 6
                        #                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1 3 2  4 6 5]' #switch 2 4 and 3 6
                        #                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1 3 2 4 6 5  7 9 8 ]' #switch 2 4 and 3 6
                        #                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1, 3, 2, 4, 5, 7, 6  ] #switch 2 4 and 3 6
#                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1, 3, 2, 4 ] #switch 2 4 and 3 6
                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1, 3, 2, 4,5,6, 8, 7 ] #switch 2 4 and 3 6
                    end
                    
                    #                    println([at1, o1, at2, o2, at3,  symb], tot, " ",  n, " " , data_info[[at1, o1, at2, o2, at3,  symb]] )
                    tot += n

                    #                    if same_at
                    #                        data_info[[o2, o1, symb]] = data_info[[o1, o2, symb]]
                    #                    end
                    
                    
                end
            end
            return tot
        end

        #        if at_list[2] == at_list[3]
        #            same_at_on = true
        #        else
        #            same_at_on = false
        #       end
        

        function get3bdy_onsite(n, same_at,symb, start, at1, at2, at3)
            #            if at2 == at3  #|| at1 == at2 || at1 == at3
            #                same_at = true
            #            else
            #                same_at = false
            #            end

            #            orbs1 = atoms[at1].orbitals

            #            println("get3bdy_onsite $at1 $at2 $at3 $n")

            tot=start
            #            for o1 in orbs1
            #                data_info[[at1, o1,at2, at3,  symb]] = collect(tot+1:tot+n]
            #                data_info[[at1, o1,at3, at2,  symb]] = collect(tot+1:tot+n)

            #                push!(orbs, (at1, o1,at2, at3,  symb))
            #                push!(orbs, (at1, o1,at3, at2,  symb))

            orbs1 = atoms[at1].orbitals
            if same_at
                for o1 in orbs1
                    
                    if [at1, at2, at3, o1,  symb] in keys(data_info)
                        continue
                    end

                    data_info[[at1, at2, at3,o1,  symb]] = collect(tot+1:tot+n)
                    data_info[[at1, at3, at2,o1,  symb]] = collect(tot+1:tot+n)
                    tot += n                               #       1 2 3 4 5 6 7 8
                end
            else
                for o1 in orbs1
                    if [at1, at2, at3, o1,  symb] in keys(data_info)
                        continue
                    end
                    data_info[[at1, at2, at3,o1,  symb]] = collect(tot+1:tot+n)
                    data_info[[at1, at3, at2,o1,  symb]] = collect(tot+1:tot+n)
                end
                tot += n
            end
#                    data_info[[at1, at3, at2,o1,  symb]] = tot .+ [1, 3, 2, 4]
                    
                    #                    data_info[[at1, o1,at2, at3,  symb]] = collect(tot+1:tot+n)
                    #                    data_info[[at1, o1,at3, at2,  symb]] = tot .+ [1 3 2 4]'
#                end
#                tot += n                               #       1 2 3 4 5 6 7 8
#            end                


            return tot
        end
        
        #this is not efficient storage, we are reassigning permutations multiple times.
        tot_size = 0
        for p in perm_ij 
            if p[1] == p[2]
                tot_size = get3bdy(n_3body_same, :H, tot_size, p[1], p[2], p[3])
            else
                tot_size = get3bdy(n_3body, :H, tot_size, p[1], p[2], p[3])
            end
        end

        #        println("data_info between")
        #        println(data_info)

        for p in perm_on
            #            if  (p[1] == p[2] ||  p[2] == p[3] || p[1] == p[3])
            if p[2] == p[3] && p[2] == p[1]
                tot_size = get3bdy_onsite(n_3body_onsite_same,true, :O, tot_size, p[1], p[2], p[3]) #all diff
            else
                tot_size = get3bdy_onsite(n_3body_onsite,false, :O, tot_size, p[1], p[2], p[3]) #
            end
        end
        
        totHO = tot_size

    else
        println("error, only 2 or 3 body terms, you gave me : ", at_list)
    end
    return totHO ,totS, data_info, orbs

end            


"""
    function get_data_info(at_set, dim)

Figure out the arrangement of data in a `coefs` file.

Loops over various combinations of orbitals and atoms and assigns them places in `datH` and `datS`, depending on the terms included in the model and the dimensionaly.
"""
function get_data_info_v1(at_set, dim)

    
    data_info = Dict{Any, Array{Int64,1}}()
    orbs = []
    if dim == 2 #2body

        at_list = Symbol.([i for i in at_set])
#        println(at_list)
        if length(at_list) == 1
            at_list = [at_list[1], at_list[1]]
        end
        sort!(at_list)
#        println(at_list)

        orbs1 = atoms[at_list[1]].orbitals
        orbs2 = atoms[at_list[2]].orbitals

        at1 = at_list[1]
        at2 = at_list[2]

        if at1 == at2
            same_at = true
        else
            same_at = false
        end

#        orbs = []

        #2body part
        function get2bdy(n, symb)
            tot=0
            for o1 in orbs1
                for o2 in orbs2
                    if same_at && ((o2 == :s && o1 == :p) || (o2 == :s && o1 == :d) || (o2 == :p && o1 == :d))
                        continue
                    end
#                    push!(orbs, (o1, o2, symb))
                    if o1 == :s && o2 == :s
                        data_info[[at1, o1, at2, o2, symb]] = tot+1:tot+n
                        data_info[[at2, o2, at1, o1, symb]] = tot+1:tot+n
                        tot += n
                        
                    elseif (o1 == :s && o2 == :p ) || (o1 == :p && o2 == :s )
                        data_info[[at1, o1, at2, o2, symb]] = tot+1:tot+n
                        data_info[[at2, o2, at1, o1, symb]] = tot+1:tot+n
                        tot += n
#                        if same_at
#                            data_info[[o2, o1, symb]] = data_info[(o1, o2, symb)]
#                        end
                        
                    elseif (o1 == :p && o2 == :p )
                        data_info[[at1, o1, at2, o2, symb]] = tot+1:tot+n*2
                        data_info[[at2, o2, at1, o1, symb]] = tot+1:tot+n*2
                        tot += n*2

#                    elseif (o1 == :p && o2 == :p )
#                        data_info[[at1, o1, at2, o2, symb]] = tot+1:tot+n*2
#                        data_info[[at2, o2, at1, o1, symb]] = tot+1:tot+n*2
#                        tot += n*2

                    elseif (o1 == :s && o2 == :d ) || (o1 == :d && o2 == :s )
                        data_info[[at1, o1, at2, o2, symb]] = tot+1:tot+n
                        data_info[[at2, o2, at1, o1, symb]] = tot+1:tot+n
                        tot += n

                    elseif (o1 == :p && o2 == :d ) || (o1 == :d && o2 == :p )
                        data_info[[at1, o1, at2, o2, symb]] = tot+1:tot+n*2
                        data_info[[at2, o2, at1, o1, symb]] = tot+1:tot+n*2
                        tot += n*2

                    elseif (o1 == :d && o2 == :d ) 
                        data_info[[at1, o1, at2, o2, symb]] = tot+1:tot+n*3
                        data_info[[at2, o2, at1, o1, symb]] = tot+1:tot+n*3
                        tot += n*3


                        
                    end
                end
            end
            return tot
        end

        totH = get2bdy(n_2body, :H)
        totS = get2bdy(n_2body_S, :S)
#        println("totH $totH totS $totS")
        #onsite part
        function getonsite(atX,orbsX, tot, n)
            for o1 in orbsX
                for o2 in orbsX
                    if (o2 == :s && o1 == :p) || (o2 == :s && o1 == :d) || (o2 == :p && o1 == :d)
                        continue
                    end

#                    push!(orbs, (o1, o2, :O))
                    if o1 == :s && o2 == :s
                        data_info[[atX, o1, o2, :O]] = tot+1:tot+n
#                        println("data_info"[, (atX, o1, o2, :O], tot+1:tot+n)
                        tot += n
                    elseif (o1 == :s && o2 == :p )
                        data_info[[atX, o1, o2, :O]] = tot+1:tot+n
                        data_info[[atX, o2, o1, :O]] = data_info[[atX, o1, o2, :O]]
                        tot += n
                    elseif (o1 == :p && o2 == :p )
                        data_info[[atX, o1, o2, :O]] = tot+1:tot+n*2
                        tot += n*2

                    elseif o1 == :s && o2 == :d
                        data_info[[atX, o1, o2, :O]] = tot+1:tot+n
                        data_info[[atX, o2, o1, :O]] = data_info[[atX, o1, o2, :O]]
                        tot += n

                    elseif o1 == :p && o2 == :d
                        data_info[[atX, o1, o2, :O]] = tot+1:tot+n
                        data_info[[atX, o2, o1, :O]] = data_info[[atX, o1, o2, :O]]
                        tot += n

                    elseif o1 == :d && o2 == :d
                        data_info[[atX, o1, o2, :O]] = tot+1:tot+n*2
                        data_info[[atX, o2, o1, :O]] = data_info[[atX, o1, o2, :O]]
                        tot += n*2

                    end
                end
            end
            return tot
        end


#PUREONSITE
#        if same_at #true onsite terms
#           for o in orbs1
##                println("true onsite ", o)
#                data_info[[at1, o, :A]] = [totH+1]
#                totH += 1
#            end
#        end


        totHO = getonsite(at1, orbs1, totH, n_2body_onsite)

        if !(same_at) #need reverse if not same atom
            totHO = getonsite(at2, orbs2, totHO, n_2body_onsite)
        end

    elseif dim == 3 #3body
        
        totS = 0 #no 3body overlap terms

        at_list = Symbol.([i for i in at_set])
        sort!(at_list)
        
        if length(at_list) == 1


            #permutations are trivial
            perm_ij = [[at_list[1], at_list[1], at_list[1]]]
            perm_on = [[at_list[1], at_list[1], at_list[1]]]
        elseif length(at_list) == 2


            #unique permutations
            perm_ij = [[at_list[1], at_list[1], at_list[2]] ,
                       [at_list[2], at_list[2], at_list[1]] ,
                       [at_list[1], at_list[2], at_list[1]] ,
                       [at_list[1], at_list[2], at_list[2]] ]

            perm_on = [[at_list[1], at_list[2], at_list[2]] ,
                        [at_list[1], at_list[1], at_list[2]] ,
                        [at_list[2], at_list[1], at_list[1]] ,
                        [at_list[2], at_list[1], at_list[2]] ]
            
        elseif length(at_list) == 3


            #all permutations exist hij
            perm_ij = [[at_list[1], at_list[2], at_list[3]] ,
                       [at_list[1], at_list[3], at_list[2]] ,
                       [at_list[2], at_list[1], at_list[3]] ,
                       [at_list[2], at_list[3], at_list[1]] ,
                       [at_list[3], at_list[1], at_list[2]] ,
                       [at_list[3], at_list[2], at_list[1]] ]

            #onsite can flip last 2 atoms
            perm_on = [[at_list[1], at_list[2], at_list[3]] ,
                       [at_list[2], at_list[1], at_list[3]] ,
                       [at_list[3], at_list[1], at_list[2]]]

            
        else
            println("ERROR  get_data_info dim $at_set $at_list")
        end
        



        function get3bdy(n, symb, start, at1, at2, at3)
            tot=start

            orbs1 = atoms[at1].orbitals
            orbs2 = atoms[at2].orbitals

            if at1 == at2
                same_at = true
            else
                same_at = false
            end

            for o1 in orbs1
                for o2 in orbs2
                    if same_at && ((o2 == :s && o1 == :p) || (o2 == :s && o1 == :d) || (o2 == :p && o1 == :d))
                        continue
                    end
                    
#                    push!(orbs, (o1, o2, symb))

                    if same_at
                        data_info[[at1, o1, at2, o2, at3,  symb]] = collect(tot+1:tot+n)
                        data_info[[at2, o2, at1, o1, at3,  symb]] = collect(tot+1:tot+n)
                    else                                                  #[1 2 3 4 5 6 7 8  9  10 11 12 13 14 15 16 17 18]
                        data_info[[at1, o1, at2, o2, at3,  symb]] = collect(tot+1:tot+n)
#                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1 4 6 2 5 3 7 10 12 8  11 9  13 16 18 14 17 15]' #switch 2 4 and 3 6
#                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1 4 6 2 5 3 7 10 12 8  11 9  ]' #switch 2 4 and 3 6
#                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1 3 2 4 5 7 6 8 9]' #switch 2 4 and 3 6
#                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1 3 2 4 5 7 6]' #switch 2 4 and 3 6
#                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1 3 2  4 6 5]' #switch 2 4 and 3 6
#                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1 3 2 4 6 5  7 9 8 ]' #switch 2 4 and 3 6
                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1, 3, 2, 4, 5, 6, 8, 7  ] #switch 2 4 and 3 6
                    end
                    
#                    println([at1, o1, at2, o2, at3,  symb], tot, " ",  n, " " , data_info[[at1, o1, at2, o2, at3,  symb]] )
                    tot += n

#                    if same_at
#                        data_info[[o2, o1, symb]] = data_info[[o1, o2, symb]]
#                    end
                        
                    
                end
            end
            return tot
        end

#        if at_list[2] == at_list[3]
#            same_at_on = true
#        else
#            same_at_on = false
 #       end
        

        function get3bdy_onsite(n, same_at,symb, start, at1, at2, at3)
#            if at2 == at3  #|| at1 == at2 || at1 == at3
#                same_at = true
#            else
#                same_at = false
#            end

            orbs1 = atoms[at1].orbitals

#            println("get3bdy_onsite $at1 $at2 $at3 $n")

            tot=start
            for o1 in orbs1
#                data_info[[at1, o1,at2, at3,  symb]] = collect(tot+1:tot+n]
#                data_info[[at1, o1,at3, at2,  symb]] = collect(tot+1:tot+n)

#                push!(orbs, (at1, o1,at2, at3,  symb))
#                push!(orbs, (at1, o1,at3, at2,  symb))

                if same_at
                    data_info[[at1, o1,at2, at3,  symb]] = collect(tot+1:tot+n)
                    data_info[[at1, o1,at3, at2,  symb]] = collect(tot+1:tot+n)
                else
                    data_info[[at1, o1,at2, at3,  symb]] = collect(tot+1:tot+n)
#                    data_info[[at1, o1,at3, at2,  symb]] = collect(tot+1:tot+n)
                    data_info[[at1, o1,at3, at2,  symb]] = tot .+ [1, 3, 2, 4]

#                    data_info[[at1, o1,at2, at3,  symb]] = collect(tot+1:tot+n)
#                    data_info[[at1, o1,at3, at2,  symb]] = tot .+ [1 3 2 4]'
                end
                tot += n                               #       1 2 3 4 5 6 7 8

            end
            return tot
        end
        
        #this is not efficient storage, we are reassigning permutations multiple times.
        tot_size = 0
        for p in perm_ij 
            if p[1] == p[2]
                tot_size = get3bdy(n_3body_same, :H, tot_size, p[1], p[2], p[3])
            else
                tot_size = get3bdy(n_3body, :H, tot_size, p[1], p[2], p[3])
            end
        end

#        println("data_info between")
#        println(data_info)

        for p in perm_on
#            if  (p[1] == p[2] ||  p[2] == p[3] || p[1] == p[3])
            if p[2] == p[3]
                tot_size = get3bdy_onsite(n_3body_onsite_same,true, :O, tot_size, p[1], p[2], p[3]) #all diff
            else
                tot_size = get3bdy_onsite(n_3body_onsite,false, :O, tot_size, p[1], p[2], p[3]) #
            end
        end
    
        totHO = tot_size

    else
        println("error, only 2 or 3 body terms, you gave me : ", at_list)
    end
    return totHO ,totS, data_info, orbs

end            
                




Base.show(io::IO, d::coefs) = begin

    println(io, "coeffs ", d.names)
    for key in keys(d.inds)
        i = d.inds[key]
        if key[end] == :S
            println(io, key, ": " , d.datS[i])
        else
            println(io, key, ": " , d.datH[i])
        end
    end
    println(io, "min dist ", round(d.min_dist, digits=4))
    println(io, "dist frontier ")
    for t in keys(d.dist_frontier)
        println(io, t, "    ", d.dist_frontier[t])
    end

    println(io)
end




#=
#old version, slower.
function distances_etc_3bdy(crys, cutoff=missing, cutoff2=missing; var_type=Float64)

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
    
    R = get_grid(crys, 35.0)

    nr = (R[1]*2+1)*(R[2]*2+1)*(R[3]*2+1)
    
    dist_arr = zeros(var_type, crys.nat,crys.nat,nr,4)
                     
    R_f = zeros(var_type, 3)
    lmn = zeros(var_type, 3)
    dist = 0.0
    
    c=0

    c_zero = 0
    
    R_keep = zeros(Int64, 0,4)
    R_keep_ab = zeros(Int64, 0,7)

    R_keep2 = Dict()
    for a = 1:crys.nat
        for b = 1:crys.nat
            R_keep2[a,b] = [zeros(Int64,0, 7),zeros(Int64,0, 7)]
            #R_keep2[a,b] = [zeros(Int64,crys.nat^2*nr, 7),zeros(Int64,crys.nat^2*nr, 7)]
        end
    end

    dmin_types = Dict()
    dmin_types3 = Dict()
    for t1 = crys.stypes
        for t2 in crys.stypes
#            dmin_types[Set((t1,t2))] = cutoff
            dmin_types[Set((t1,t2))] = get_cutoff(t1,t2)[1]
            for t3 in crys.stypes
#                dmin_types3[Set((t1,t2,t3))] = cutoff2
                dmin_types3[Set((t1,t2,t3))] = get_cutoff(t1,t2,t3)
            end
        end
    end

    keep_counter = 1
    dR = zeros(var_type, 3)
    dist = 0.0

    coords_ab = zeros(var_type, 3, crys.nat, crys.nat)
    for a = 1:crys.nat
        for b = 1:crys.nat
            coords_ab[:,a,b] = crys.coords[a,:] .- crys.coords[b,:]
        end
    end
#    println("two body")
    for r1 = -R[1]:R[1]
        R_f[1] = r1
        for r2 = -R[2]:R[2]
            R_f[2] = r2
            for r3 = -R[3]:R[3]
                c+=1
                R_f[3] = r3
#                R_f[:] = [r1,r2,r3]
                found = false
                found2 = false                
                for a = 1:crys.nat
                    ta = crys.stypes[a]
                    for b = 1:crys.nat
                        tb = crys.stypes[b]            

                        dR[:] = ( coords_ab[:,a,b]  .+ R_f[:])'*crys.A

                       
                

                        dist = (dR'*dR)^0.5
#                        dist = sum(dR.^2)^0.5

                        if dist > 1e-7
                            lmn[:] .= dR/(dist )
                        else
                            lmn[:] .= 0.0
                        end

                        dist_arr[a,b,c,1] = dist
                        dist_arr[a,b,c,2:4] .= lmn[:]

#                        if dist < cutoff

                        cutoff = get_cutoff(ta,tb)[1]
                        
                        if dist < cutoff
                            found = true
                            R_keep_ab = [R_keep_ab; [c, a, b, r1,r2,r3, keep_counter]']
                            #R_keep_ab[keep_counter,:] .= [c, a, b, r1,r2,r3, keep_counter]
                            
                            if dist < dmin_types[Set((ta,tb))] && dist > 1e-7
                                dmin_types[Set((ta,tb))] = dist
                            end
                            
                        end
                        if dist < cutoff
                            R_keep2[a,b][1] = [R_keep2[a,b][1]; [keep_counter, a, b, r1,r2,r3, c]']
                        end
                        if dist < cutoff3bX
                            R_keep2[a,b][2] = [R_keep2[a,b][2]; [keep_counter, a, b, r1,r2,r3, c]']
                        end
                        
                    end
                end
                if found
#                    push!(R_keep, [c, r1,r2,r3])
                    R_keep = [R_keep; [c, r1,r2,r3]']
                    if r1 == 0 && r2 == 0 && r3 == 0
#                        c_zero = size(R_keep)[1]
                        c_zero = keep_counter
                    end

                    keep_counter += 1
                    
                end

            end
        end
    end

#    R_keep_ab = R_keep_ab[1:keep_counter, :]
    
    #    println("total c $c")
#    println("keep_counter $keep_counter")
    
#    array_ind3 = zeros(Int64, 0, 5)
#    array_floats3 = zeros(Float64, 0, 14)
    

    MEMCHUNK = min(keep_counter^2 * crys.nat^3, 5000)

    array_ind3 = zeros(Int64, MEMCHUNK, 5)
    array_floats3 = zeros(var_type, MEMCHUNK , 14)

    TOTMEM = MEMCHUNK


#    println("nkeep1 = ", size(R_keep2[1,1][1])[1])
#    println("nkeep2 = ", size(R_keep2[1,1][2])[1])

    counter = 0
#    println("three body")

    if threebody
    for a = 1:crys.nat
        ta = crys.stypes[a]
        for b = 1:crys.nat
            tb = crys.stypes[b]
            cutoff = get_cutoff(ta,tb)[1]


            R2 = R_keep2[a,b][1]
            nkeep2 = size(R2)[1]
            for c1 = 1:nkeep2
                cind1 = R2[c1][1]
                
                c_dist_ind1 = R2[c1,7]
                dist_ab = dist_arr[a,b,c_dist_ind1,1]
                if dist_ab < 1e-5
                    continue
                end
                
                lmn_ab = dist_arr[a,b,c_dist_ind1,2:4]
                

                cut_ab = cutoff_fn(dist_ab, cutoff - cutoff_length, cutoff)
        
                for c = 1:crys.nat
                    tc = crys.stypes[c]
                    cutoff3 = get_cutoff(ta,tb,tc)
                    cut_ab2 = cutoff_fn(dist_ab, cutoff3 - cutoff_length, cutoff3)

                    R22 = R_keep2[a,c][2]
                    nkeep22 = size(R22)[1]
                    for c2 = 1:nkeep22
                        cind2 = R22[c2][1]
                        c_dist_ind2 = R22[c2,7]
                        dist_ac = dist_arr[a,c,c_dist_ind2, 1]
                        
                        if dist_ac > cutoff3
                            continue
                        end
                        
                        lmn_ac = dist_arr[a,c,c_dist_ind2,2:4]

                        dist_bc = sum( (dist_ab*lmn_ab[:] - dist_ac*lmn_ac).^2)^0.5

                        if dist_bc > cutoff3
                            continue
                        end


                        if dist_bc > 1e-5
                            lmn_bc = -(dist_ab*lmn_ab[:] .- dist_ac*lmn_ac) ./ dist_bc
                        end
                        if dist_bc < cutoff3 && dist_ab > 1e-5 && dist_bc > 1e-5 && dist_ac > 1e-5

                            d0=dmin_types3[Set((ta,tb,tc))]
                            if (dist_bc + dist_ab + dist_ac)/3.0 < d0
                                dmin_types3[Set((ta,tb,tc))]  = (dist_bc + dist_ab + dist_ac)/3.0
                            end



                            cut_ac = cutoff_fn(dist_ac, cutoff3 - cutoff_length, cutoff3)
                            cut_bc = cutoff_fn(dist_bc, cutoff3 - cutoff_length, cutoff3)
                            
                            counter += 1
                            if counter > TOTMEM #need more memory
                                TOTMEM += MEMCHUNK
                                array_ind3 = [array_ind3;zeros(Int64, MEMCHUNK, 5)]
                                array_floats3 = [array_floats3; zeros(var_type, MEMCHUNK , 14)]
                            end
                                
                            array_ind3[counter,:] = [a,b,c,cind1,cind2]
                            array_floats3[counter, :] = [dist_ab, dist_ac, dist_bc, lmn_ab[1],lmn_ab[2], lmn_ab[3], lmn_ac[1],lmn_ac[2], lmn_ac[3], lmn_bc[1],lmn_bc[2], lmn_bc[3], cut_ab*cut_bc*cut_ac, cut_ab2*cut_bc*cut_ac]
#                            array_ind3=[array_ind3; [a,b,c,cind1, cind2]' ]
#                            array_floats3 = [array_floats3; [dist_ab, dist_ac, dist_bc, lmn_ab[1],lmn_ab[2], lmn_ab[3], lmn_ac[1],lmn_ac[2], lmn_ac[3], lmn_bc[1],lmn_bc[2], lmn_bc[3], cut_ab*cut_bc*cut_ac, cut_ab2*cut_bc*cut_ac]']
                        end
                    end
                end
            end
        end
    end
    end

    array_ind3 = array_ind3[1:counter,:]
    array_floats3 = array_floats3[1:counter,:]

    return R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3
    
end
=#

########################################################################################################################################################

"""
    function trim_dist(tbc, cutoff=18.0001)

Reduce the atom-atom hamiltonian terms longer than cutoff.

Not used in typical code, but can make `tbc` run faster / reduce memory.
"""
function trim_dist(tbc, cutoff=18.0001)

    R_keep,R_keep2, dist_arr, c_zero = distances_etc(tbc.crys, cutoff)
    nkeep=size(R_keep)[1]

    ind_arr = zeros(Int64, nkeep, 3)
    R_keep_dict = Dict()
    for c= 1:nkeep
        ind_arr[c,:] = R_keep[c][2:4]
        cind = R_keep[c][1]

        R_keep_dict[R_keep[c][2:4]] = (c, cind)
    end


    ind2orb, orb2ind, etotal, nval = orbital_index(tbc.crys)
    
    nwan = tbc.tb.nwan


    var_type = typeof(real(tbc.tb.H[1,1,1]))
    
    H = zeros(var_type, nwan, nwan, nkeep)
    S = zeros(var_type, nwan, nwan, nkeep)    

    i = 0
    for c = 1:tbc.tb.nr
        ind_old = tbc.tb.ind_arr[c,:]

        if haskey(R_keep_dict, ind_old)
            found = true
            (cnew, cind) = R_keep_dict[ind_old]
        else
            found = false
        end
#        ind_arr[c,:] = R_keep[c][2:4]
#        cind = R_keep[c][1]
        for o1 = 1:nwan
            a1,t1,s1 = ind2orb[o1]
            for o2 = 1:nwan
                a2, t2,s2 = ind2orb[o2]
                if found 
                    dist = dist_arr[a1,a2,cind,1]

                    if (dist <= cutoff )
                        H[o1,o2,cnew] += real(tbc.tb.H[o1,o2,c])
                        S[o1,o2,cnew] += real(tbc.tb.S[o1,o2,c])
                        i += 1
                    else
                        found = false
                    end
                end                        
                if found == false
                    H[o1,o2,c_zero] += real(tbc.tb.H[o1,o2,c])
                    S[o1,o2,c_zero] += real(tbc.tb.S[o1,o2,c])
                end
            end
        end
    end

    println("dist_trim keep $i out of ", nwan*nwan*tbc.tb.nr)
    
    nelec = tbc.nelec
    dftenergy = tbc.dftenergy
    tb = make_tb(H, ind_arr, S)
    tbc_new = make_tb_crys(tb, tbc.crys, nval, dftenergy)

    return tbc_new

end    

########################################################################################################################################################
"""
    calc_tb_fast(crys::crystal; database=missing, use_threebody=true, use_threebody_onsite=true)

Construct `tb_crys` from crystal stucture, but does not solve. This is
usually called internally by functions like `scf_energy`, but you can
use it directly if you want. Until you do a SCF energy calculation, the
electron density and Fermi level will be wrong.

# Arguments
- `crys::crystal` - Required crystal structure
- `database=missing` - Source of coefficients. Will load from default source if not specified.
- `use_threebody=true` - Use three-body off-site interactions. Only turn off for testing purposes.
- `use_threebody_onsite=true` - Use three-body on-site interactions. Only turn off for testing purposes.
- `verbose=true` - set to false for less output.
- `var_type=missing` - variable type of `tb_crys`. Default is `Float64`.
"""
function calc_tb_fast(crys::crystal, database=missing; reference_tbc=missing, verbose=true, var_type=missing, use_threebody=true, use_threebody_onsite=true, gamma=missing, u3=missing,background_charge_correction=0.0, screening=1.0, set_maxmin=false, check_frontier=true, check_only=false, repel = true)

    #    use_threebody= false
    #    use_threebody_onsite=false
    
####    verbose = true

#    println("repel $repel -------------------------------")
    
    if verbose
        println()
        println("-----")
        println("Construct tight-binding model from crystal structure")
        println()
    end

    if ismissing(var_type)
        var_type=Float64
    end

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

    if (use_threebody || use_threebody_onsite ) && !ismissing(database)
        #        parallel =true
        #        if parallel

        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3, Rind = distances_etc_3bdy_parallel(crys,cutoff2X, cutoff3bX,var_type=var_type)
#        println("done dist")
        #        else
        #            R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy(crys,cutoff2X, cutoff3bX,var_type=var_type)
        #        end        
    else
        #        parallel = true
        #        if parallel

        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3, Rind = distances_etc_3bdy_parallel(crys,cutoff2X, 0.0,var_type=var_type)

        #else
        #        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy(crys,cutoff2X, 0.0,var_type=var_type)
        #    end
        
    end
#    println("end d")
    #    R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero = distances_etc_3bdy(crys,cutoff2X, 7.0)

    #    println("size R_keep_ab ", size(R_keep_ab))
    #    println("size array_ind3 ", size(array_ind3))

    #    if ismissing(dist_etc)
    #    else
    #        (R_keep, R_keep2, dist_arr, c_zero) = dist_etc
    #    end
    
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
        #        for key in keys(dmin_types3)
        #            for key2 in keys(database)
        #                if key == Set(key2)
        #                    if dmin_types3[key] < database[key2].min_dist*1.04 && length(key2) == 3
        #                        println("WARNING : structure has 3body distances less than or close to min fitting data distances, may result in errors")
        #                        println(key, "  ",key2, " : ", dmin_types3[key], " <~ ", database[key2].min_dist)
        #                        within_fit = false
        #                    end
        #                end
        #            end
        #        end
        
        c_zero_ref=1
        if !(ismissing(reference_tbc))
            if size(reference_tbc.tb.ind_arr)[1] > 1
                c_zero_ref = reference_tbc.tb.r_dict[[0,0,0]]
            end
        end
    end

    if verbose println("check_frontier") end
    if !ismissing(database) && check_frontier
        #    if false
        #        diststuff = (R_keep, R_keep_ab, dist_arr, c_zero, dmin_types, dmin_types3)
        diststuff = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3, Rind
        violation_list, vio_bool, repel_vals = calc_frontier(crys, database, test_frontier=true, diststuff=diststuff, verbose=verbose, var_type=var_type)
        if vio_bool == false
            within_fit = false
        end
    else
        repel_vals = zeros(var_type, crys.nat)
    end
    
    if check_only==true
        return within_fit, sum(abs.(repel_vals)) < 1e-12
    end
    
    nwan = length(keys(ind2orb))

    nkeep=size(R_keep)[1]
    #    nkeep2=size(R_keep2)[1]    
    #    println("nkeep, $nkeep, nkeep2, $nkeep2")
    

    H = zeros(var_type, 1, nwan, nwan, nkeep)
    S = zeros(var_type, nwan, nwan, nkeep)    


    ind_arr = zeros(Int64, nkeep, 3)
    ind_arr[:,:] = R_keep[:,2:4]

    #    lmn = zeros(var_type, 3)
    #    dist = 0.0
    #    lmn31 = zeros(var_type, 3)
    #    dist31 = 0.0

    #    lmn32 = zeros(var_type, 3)
    #    dist32 = 0.0

    #    lmn41 = zeros(var_type, 3)
    #    dist14 = 0.0
    #    dist43 = 0.0


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
        end
    end

    
    warned = false
    warned_onsite = false

    nkeep_ab = size(R_keep_ab)[1]

    if !ismissing(database)

        if verbose println("2body") end
        LMN = zeros(var_type, 3, nthreads())

        @threads for c = 1:nkeep_ab
            id = threadid()

            #        ind_arr[c,:] = R_keep_ab[c][4:6]
            cind = R_keep_ab[c,1]
            cham = R_keep_ab[c,7]
            a1 = R_keep_ab[c,2]
            a2 = R_keep_ab[c,3]

            t1 = crys.stypes[a1]
            t2 = crys.stypes[a2]
            
            coef = database[(t1,t2)]

            indH, indS, inH, inS = coef.inds_int[[t1,t2]]

            cutoff2Xt = get_cutoff(t1,t2)[1]
            dist_a = dist_arr[a1,a2,cind,1]
            if (dist_a > cutoff2Xt || dist_a < 1e-5)
                continue
            end
            LMN[:,id] .= (dist_arr[a1,a2,cind,2:4])

            lag = two_body_S(dist_a)

            if dist_a < cutoff2Xt - cutoff_length
                cut = 1.0
            else
                cut = cutoff_fn(dist_a, cutoff2Xt - cutoff_length, cutoff2Xt)
            end
            cutoff2Xa = get_cutoff(t1,t2)[1]

            ####            for o1 = orb2ind[a1]
            ####                a1a,t1,s1 = ind2orb[o1]
            ####                sum1 = summarize_orb(s1)
            ####                for o2 = orb2ind[a2]
            ####                    a2a,t2,s2 = ind2orb[o2]
            ####                    at_set = Set((t1,t2))
            ####
            ####                    cutoff2Xa = get_cutoff(t1,t2)[1]
            ####
            ####                    #                println("asdf $c $cind $a1 $a2 $o1 $o2 $s1 $s2")
            ####                    #        for o1 = 1:nwan
            ####                    #            a1,t1,s1 = ind2orb[o1]
            ####                    #            for o2 = 1:nwan
            ####                    #                a2, t2,s2 = ind2orb[o2]
            ####                    
            ####                    dist = dist_arr[a1,a2,cind,1]
            ####
            ####
            ####                    if (dist > cutoff2Xa || dist < 1e-5)
            ####                        continue
            ####                    end
            ####                    
            ####                    lmn = dist_arr[a1,a2,cind,2:4]
            ####                    #               println("jkl $t1 $t2 $s1 $s2 $dist $lmn")
            ####
            ####
            ####                    (h,s) = calc_twobody(t1,t2,s1,s2,dist,lmn, database)
            ####
            ####                    if dist < cutoff2Xa - cutoff_length
            ####                        cut = 1.0
            ####                    else
            ####                        cut = cutoff_fn(dist, cutoff2Xa - cutoff_length, cutoff2Xa)
            ####                    end
            ####
            ####                    H[o1, o2, cham] += h  *cut
            ####                    S[o1, o2, cham] += s  *cut

            
            #            for o1 = orb2ind[a1]
            #                a1a,t1,s1 = ind2orb[o1]
            #                sum1 = summarize_orb(s1)
            #                for o2 = orb2ind[a2]
            #                    a2a,t2,s2 = ind2orb[o2]
            #                    at_set = Set((t1,t2))

            
            for o1x = 1:norb[a1]
                o1 = orbs[a1,o1x]
                s1 = sorbs[a1,o1x]
                sum1 = sumorbs[a1,o1x]
                
                for o2x = 1:norb[a2]
                    o2 = orbs[a2,o2x]
                    s2 = sorbs[a2,o2x]
                    sum2 = sumorbs[a2,o2x] 

                    (hw,sw) = calc_twobody_faster(t1,t2,s1,s2,sum1, sum2, dist_a,LMN[:,id], coef, indH[sum1,sum2,:], indS[sum1,sum2,:],lag)

                    H[1, o1, o2, cham] += hw  *cut
                    S[o1, o2, cham] += sw  *cut


                end
            end
        end

        #############
        #threebody

        #        memory0=zeros(var_type, 3)
        #        memory1=zeros(var_type, 3)
        #        memory2=zeros(var_type, 3)
        #        memoryV=zeros(var_type, n_3body)

        #        lmn = zeros(var_type, 3)

        H_thread = zeros(var_type,  nwan, nwan,  nkeep,  nthreads() )


        memory0_th=zeros(var_type, maximum([n_3body, n_3body_onsite, n_3body_onsite_same, n_3body_same]) , nthreads())
        memory1_th=zeros(var_type, maximum([n_3body, n_3body_onsite, n_3body_onsite_same, n_3body_same]) , nthreads())
        memory2_th=zeros(var_type, maximum([n_3body, n_3body_onsite, n_3body_onsite_same, n_3body_same]) , nthreads())
        memoryV_th=zeros(var_type, maximum([n_3body, n_3body_onsite, n_3body_onsite_same, n_3body_same]) , nthreads())

        hh = zeros(var_type, 3,3, nthreads())
        Htemp = zeros(var_type, 16,16, nthreads())
        
        #        v =Array{Symbol}(undef, 6, nthreads())
        

        #        sdict = Dict()
        #        sumdict = Dict()
        #        for a1 = 1:crys.nat
        #            sdict[a1] = Symbol[]
        #            sumdict[a1] = Symbol[]
        #            
        #            for o1 = orb2ind[a1]
        #                a1a,t1,s1 = ind2orb[o1]
        #                sum1 = summarize_orb(s1)
        #                push!(sdict[a1], s1)
        #                push!(sumdict[a1], sum1)
        #            end
        #            
        #        end
        




        if verbose println("3body") end
        if use_threebody || use_threebody_onsite
            #        if false
            @time @threads for counter = 1:size(array_ind3)[1]
                #            for counter = 1:size(array_ind3)[1]
                id = threadid()
                #id = 1
                a1 = array_ind3[counter,1]
                a2 = array_ind3[counter,2]
                a3 = array_ind3[counter,3]

                cind1 = array_ind3[counter,4]
                #                cind2 = array_ind3[counter,5]

                dist = array_floats3[counter, 1]
                dist31 = array_floats3[counter, 2]
                dist32 = array_floats3[counter, 3]
                
                #                lmn[:] = array_floats3[counter, 4:6]
                #                lmn31[:] = array_floats3[counter, 7:9]
                #                lmn32[:] = array_floats3[counter, 10:12]

#                lmn = @view array_floats3[counter, 4:6]
                lmn31 = @view array_floats3[counter, 4:6]
                lmn32 = @view array_floats3[counter, 7:9]

                memory0= @view memory0_th[:,id]
                memory1= @view memory1_th[:,id]
                memory2= @view memory2_th[:,id]
                memoryV= @view memoryV_th[:,id]


                cut = array_floats3[counter, 10]

                t1 = crys.stypes[a1]
                t2 = crys.stypes[a2]
                t3 = crys.stypes[a3]

                #                v = [t1,:s,t2,:s,t3,:H]


                
                if haskey(database, (t1, t2, t3))
                    cdat = database[(t1,t2,t3)]
                    (cindX, nindX) = cdat.inds_int[[t1,t2,t3]]
                    if use_threebody

                        
                        
                        three_body_H(dist, dist31, dist32,t1==t2, t1 !=t2 && t1 != t3 && t2 != t3, memory0=memory0, memory1=memory1, memory2=memory2, memoryV=memoryV)

                        #puts what we need in memoryV



                        for sum1 = 1:3
                            for sum2 = 1:3
                                @inbounds hh[sum1,sum2,id] = ( (@view memoryV[1:nindX[sum1, sum2]])'* (@view cdat.datH[ (@view cindX[sum1, sum2, 1:nindX[sum1, sum2]])   ]))[1]
                            end
                        end

                        sym31 = 1.0
                        sym32 = 1.0                        
                        
                        #                        Htemp = zeros(norb[a1], norb[a2])
                        Htemp[:,:, id] .= 0.0


                        for o1x = 1:norb[a1]
                            o1 = orbs[a1,o1x]
                            s1 = sorbs[a1,o1x]
                            sum1 = sumorbs[a1,o1x]
                            
                            sym31 = symmetry_factor_int(s1,1,lmn31,one ) 

                            for o2x = 1:norb[a2]
                                o2 = orbs[a2,o2x]
                                s2 = sorbs[a2,o2x]
                                sum2 = sumorbs[a2,o2x] 

                                sym32 = symmetry_factor_int(s2,1,lmn32, one)    

                                #                                @inbounds Htemp[o1x,o2x,id] += hh[sum1,sum2, id]  * sym31 * sym32
                                @inbounds Htemp[o1x,o2x,id] += hh[sum1,sum2, id] * sym31 * sym32 
                                
                                #                                @inbounds H_thread[o1, o2, cind1, id ] += hh[sum1,sum2, id]  * sym31 * sym32


                                
                                
                            end
                        end

                        @inbounds Htemp[1:norb[a1],1:norb[a2], id] .*= cut * 10^3
                        
                        @inbounds H_thread[orbs[a1,1:norb[a1]] , orbs[a2,1:norb[a2]]  , cind1, id ] .+= (@view Htemp[1:norb[a1],1:norb[a2], id])
                        
                    end
                    ############################################
                    if use_threebody_onsite
                        #        if false
                        cut2 = array_floats3[counter, 11]
                        for o1 = orb2ind[a1]
                            a1a,t1,s1 = ind2orb[o1]
                            o = calc_threebody_onsite(t1,t2,t3,s1,dist,dist31,dist32, cdat, set_maxmin=set_maxmin, memory=memoryV)

                            
                            @inbounds H_thread[ o1, o1,c_zero, id] += o  * cut2
                        end
                    elseif !warned_onsite && use_threebody_onsite
                        println("WARNING, missing 3bdy onsite ", (t1, t2, t3))
                        warned_onsite = true
                        within_fit = false
                    end
                    ###########################################

                elseif !warned
                    println("WARNING, missing 3bdy ", (t1, t2, t3))
                    within_fit = false

                    warned = true
                end
            end
        end

        #        println("thread sum")
        #        println(size(H))
        #        println(size(H_thread))

        
        H[1,:,:,:] .+= sum(H_thread, dims=4)[:,:,:]
        #        H += sum(H_thread, dims=4)[:, :,:]



        lmn = zeros(var_type, 3)

        ############ONSITE

        Hon = zeros(var_type, nwan,nwan, nthreads())
        Son = zeros(var_type, nwan,nwan, nthreads())
        
        if verbose println("onsite") end
        @threads for c = 1:nkeep_ab
            id = threadid()

            #        ind_arr[c,:] = R_keep_ab[c][4:6]
            cind = R_keep_ab[c,1]
            a1 = R_keep_ab[c,2]
            a2 = R_keep_ab[c,3]
            t1 = crys.stypes[a1]
            t2 = crys.stypes[a2]

            dist = dist_arr[a1,a2,cind,1]
            LMN[:, id] = dist_arr[a1,a2,cind,2:4]

            cutoff_onXa = get_cutoff(t1,t2)[2]

            if (dist > cutoff_onXa)
                continue
            end

            if dist < cutoff_onXa - cutoff_length
                cut = 1.0
            else
                cut = cutoff_fn(dist, cutoff_onXa - cutoff_length, cutoff_onXa)
            end

            
            for o1 = orb2ind[a1]
                a1a,t1a,s1 = ind2orb[o1]
                sum1 = summarize_orb(s1)
                for o2 = orb2ind[a1]
                    a2a,t2a,s2 = ind2orb[o2]
                    
                    if dist < 1e-5 #true onsite
                        (h,s) = calc_onsite(t1,s1,s2, database)
                        #                        S[o1, o2, c_zero] += s 
                        #                        H[o1, o2, c_zero] += h
                        Son[o1, o2, id] += s 
                        Hon[o1, o2, id] += h
                        if repel
                            if o1 == o2
                                Hon[o1, o1, id] += repel_vals[a1a] * 0.1
                            end
                        end
                        
                    else
                        o = calc_twobody_onsite(t1,t2, s1,s2,dist,LMN[:,id], database)
                        Hon[o1, o2, id] += o * cut
                        #                        H[o1, o2, c_zero] += o  * cut
                    end
                    



                end

            end
        end
        H[1, :,:,c_zero] += sum(Hon, dims=3)[:,:]
        S[:,:,c_zero] += sum(Son, dims=3)[:,:]

    end

    

    if verbose println("make") end
    if true
#        println("typeof H ", typeof(H), " " , size(H), " S ", typeof(S), " " , size(S))
        tb = make_tb(H, ind_arr, S)
        if !ismissing(database) && (haskey(database, "scf") || haskey(database, "SCF"))
            scf = database["scf"]
        else
            scf = false
        end
        tbc = make_tb_crys(tb, crys, nval, 0.0, scf=scf, gamma=gamma, u3=u3,background_charge_correction=background_charge_correction, within_fit=within_fit, screening=screening)
    end
    if verbose 
        println("-----")
        println()
    end

    return tbc

end


function calc_tb_fast_old(crys::crystal, database=missing; reference_tbc=missing, verbose=true, var_type=missing, use_threebody=true, use_threebody_onsite=true, gamma=missing, u3=missing, background_charge_correction=0.0, screening=1.0, set_maxmin=false, check_frontier=true, check_only=false)

#    use_threebody= false
#    use_threebody_onsite=false


    
    if verbose
        println()
        println("-----")
        println("Construct tight-binding model from crystal structure")
        println()
    end

    if ismissing(var_type)
        var_type=Float64
    end

    if ismissing(database)
        println("missing database, creating empty tbc")
    end

    
    if ismissing(reference_tbc)
        prepare_for_fitting = false
    else
        prepare_for_fitting = true
    end
    

    
    ind2orb, orb2ind, etotal, nval = orbital_index(crys)

    if verbose println("distances") end

    if (use_threebody || use_threebody_onsite ) && !ismissing(database)
#        parallel =true
#        if parallel

        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3,Rind = distances_etc_3bdy_parallel(crys,cutoff2X, cutoff3bX,var_type=var_type)

        #        else
#            R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy(crys,cutoff2X, cutoff3bX,var_type=var_type)
#        end        
    else
#        parallel = true
#        if parallel

        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3,Rind = distances_etc_3bdy_parallel(crys,cutoff2X, 0.0,var_type=var_type)

    #else
    #        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy(crys,cutoff2X, 0.0,var_type=var_type)
        #    end
        
    end

    #    R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero = distances_etc_3bdy(crys,cutoff2X, 7.0)

#    println("size R_keep_ab ", size(R_keep_ab))
#    println("size array_ind3 ", size(array_ind3))

#    if ismissing(dist_etc)
#    else
#        (R_keep, R_keep2, dist_arr, c_zero) = dist_etc
#    end
    
    within_fit = true
    
    if !ismissing(database)
        for key in keys(dmin_types)
            for key2 in keys(database)
                if key == Set(key2)
                    if dmin_types[key] < database[key2].min_dist*1.0199 && length(key2) == 2
                        println("WARNING : structure has 2body distances less than or close to min fitting data distances, may result in errors")
                        println(key," " ,key2, " : ", dmin_types[key], " <~ ", database[key2].min_dist)
                        within_fit = false
                    end
                end
            end
        end
#        for key in keys(dmin_types3)
#            for key2 in keys(database)
#                if key == Set(key2)
#                    if dmin_types3[key] < database[key2].min_dist*1.04 && length(key2) == 3
#                        println("WARNING : structure has 3body distances less than or close to min fitting data distances, may result in errors")
#                        println(key, "  ",key2, " : ", dmin_types3[key], " <~ ", database[key2].min_dist)
#                        within_fit = false
#                    end
#                end
#            end
#        end
        
        c_zero_ref=1
        if !(ismissing(reference_tbc))
            if size(reference_tbc.tb.ind_arr)[1] > 1
                c_zero_ref = reference_tbc.tb.r_dict[[0,0,0]]
            end
        end
    end

    if !ismissing(database) && check_frontier
#    if false
        diststuff = (R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3, Rind)
        violation_list, vio_bool = calc_frontier(crys, database, test_frontier=true, diststuff=diststuff, verbose=verbose)
        if vio_bool == false
            within_fit = false
        end
    end
    if check_only==true
        return within_fit, sum(abs.(repel_vals)) < 1e-12
    end
    
    nwan = length(keys(ind2orb))

    nkeep=size(R_keep)[1]
#    nkeep2=size(R_keep2)[1]    
#    println("nkeep, $nkeep, nkeep2, $nkeep2")
    

    H = zeros(var_type, nwan, nwan, nkeep)
    S = zeros(var_type, nwan, nwan, nkeep)    


    ind_arr = zeros(Int64, nkeep, 3)
    ind_arr[:,:] = R_keep[:,2:4]

#    lmn = zeros(var_type, 3)
#    dist = 0.0
#    lmn31 = zeros(var_type, 3)
#    dist31 = 0.0

#    lmn32 = zeros(var_type, 3)
#    dist32 = 0.0

#    lmn41 = zeros(var_type, 3)
#    dist14 = 0.0
#    dist43 = 0.0

    
    warned = false
    warned_onsite = false

    nkeep_ab = size(R_keep_ab)[1]

    if !ismissing(database)

        if verbose println("2body") end
        @time @threads for c = 1:nkeep_ab
            #        ind_arr[c,:] = R_keep_ab[c][4:6]
            cind = R_keep_ab[c,1]
            cham = R_keep_ab[c,7]
            a1 = R_keep_ab[c,2]
            a2 = R_keep_ab[c,3]
            for o1 = orb2ind[a1]
                a1a,t1,s1 = ind2orb[o1]
                sum1 = summarize_orb(s1)
                for o2 = orb2ind[a2]
                    a2a,t2,s2 = ind2orb[o2]
                    at_set = Set((t1,t2))

                    cutoff2Xa = get_cutoff(t1,t2)[1]

                    #                println("asdf $c $cind $a1 $a2 $o1 $o2 $s1 $s2")
                    #        for o1 = 1:nwan
                    #            a1,t1,s1 = ind2orb[o1]
                    #            for o2 = 1:nwan
                    #                a2, t2,s2 = ind2orb[o2]
                    
                    dist = dist_arr[a1,a2,cind,1]


                    if (dist > cutoff2Xa || dist < 1e-5)
                        continue
                    end
                    
                    lmn = dist_arr[a1,a2,cind,2:4]
                    #               println("jkl $t1 $t2 $s1 $s2 $dist $lmn")


                    (h,s) = calc_twobody(t1,t2,s1,s2,dist,lmn, database)

                    if dist < cutoff2Xa - cutoff_length
                        cut = 1.0
                    else
                        cut = cutoff_fn(dist, cutoff2Xa - cutoff_length, cutoff2Xa)
                    end

                    H[1, o1, o2, cham] += h  *cut
                    S[o1, o2, cham] += s  *cut


                end
            end
        end

        #############
        #threebody

#        memory0=zeros(var_type, 3)
#        memory1=zeros(var_type, 3)
#        memory2=zeros(var_type, 3)
#        memoryV=zeros(var_type, n_3body)

#        lmn = zeros(var_type, 3)

        H_thread = zeros(var_type, nwan, nwan, nkeep, nthreads() )


        memory0_th=zeros(var_type, 3, nthreads())
        memory1_th=zeros(var_type, 3, nthreads())
        memory2_th=zeros(var_type, 3, nthreads())
        memoryV_th=zeros(var_type, n_3body, nthreads())

#        v =Array{Symbol}(undef, 6, nthreads())
        

#        sdict = Dict()
#        sumdict = Dict()
#        for a1 = 1:crys.nat
#            sdict[a1] = Symbol[]
#            sumdict[a1] = Symbol[]
#            
#            for o1 = orb2ind[a1]
#                a1a,t1,s1 = ind2orb[o1]
#                sum1 = summarize_orb(s1)
#                push!(sdict[a1], s1)
#                push!(sumdict[a1], sum1)
#            end
#            
#        end
                
        

        if verbose println("3body") end
        if use_threebody || use_threebody_onsite
            #        if false
            for counter = 1:size(array_ind3)[1]
                #            for counter = 1:size(array_ind3)[1]
                #id = threadid()
                id = 1
                a1 = array_ind3[counter,1]
                a2 = array_ind3[counter,2]
                a3 = array_ind3[counter,3]

                cind1 = array_ind3[counter,4]
#                cind2 = array_ind3[counter,5]

                dist = array_floats3[counter, 1]
                dist31 = array_floats3[counter, 2]
                dist32 = array_floats3[counter, 3]
                
#                lmn[:] = array_floats3[counter, 4:6]
#                lmn31[:] = array_floats3[counter, 7:9]
#                lmn32[:] = array_floats3[counter, 10:12]

                lmn = array_floats3[counter, 4:6]
                lmn31 = array_floats3[counter, 7:9]
                lmn32 = array_floats3[counter, 10:12]

                memory0=memory0_th[:,id]
                memory1=memory1_th[:,id]
                memory2=memory2_th[:,id]
                memoryV=memoryV_th[:,id]


                cut = array_floats3[counter, 13]

                t1 = crys.stypes[a1]
                t2 = crys.stypes[a2]
                t3 = crys.stypes[a3]

                v = [t1,:s,t2,:s,t3,:H]

                
                if haskey(database, (t1, t2, t3))
                    cdat = database[(t1,t2,t3)]
                    if use_threebody

                        
                        three_body_H(dist, dist31, dist32,t1==t2, t1 !=t2 && t1 != t3 && t2 != t3, memory0=memory0, memory1=memory1, memory2=memory2, memoryV=memoryV) #puts what we need in memoryV
                        
                        for o1 = orb2ind[a1]
                        #    id = threadid()

                            a1a,t1,s1 = ind2orb[o1]
                            sum1 = summarize_orb(s1)

                            sym31 = symmetry_factor(s1,:s,lmn31, [1.0])

                            
                            v[2] = sum1
                            for o2 = orb2ind[a2]

                                a2a,t2,s2 = ind2orb[o2]
                                
                                sum2 = summarize_orb(s2)
                                v[4] = sum2

                                sym32 = symmetry_factor(s2,:s,lmn32, [1.0])    
                                ind = cdat.inds[v]
                                s = length(ind)
                                h = ( (@view memoryV[1:s])'* (@view cdat.datH[ind]))[1] * 10^3        

                                @inbounds H_thread[o1, o2, cind1, id] += h  * cut * sym31 * sym32
                                
#                                @inbounds H_thread[o1, o2, cind1, id] +=  cut * sym31 * sym32 
                                

                                #                                
                                #o1 = summarize_orb(orb1)
                                #o2 = summarize_orb(orb2)    

#                                v = (t1,sum1,t2,sum2,t3,:H)
                                #                                ind = cdat.inds[[t1,sum1,t2,sum2,t3,:H]]
                                #                                ind = cdat.inds[v]

                                #                                ind = cdat.inds[v]
#                                ind = cdat.inds[[t1,sum1,t2,sum2,t3,:H]]
#                                ind = cdat.inds[v]
#                        if true


                                #                            h=0.0

#                                h = calc_threebody(t1,t2,t3,s1,s2,dist,dist31,dist32,lmn, lmn31,lmn32, database, memory0, memory1, memory2, memoryV, set_maxmin=set_maxmin)

#                                h = calc_threebody(cdat, t1,t2,t3,s1,s2,dist,dist31,dist32,lmn, lmn31,lmn32, database, memory0, memory1, memory2, memoryV, precalc=true, set_maxmin=set_maxmin)
 #                               h = calc_threebody(cdat,ind  , t1,t2,t3,s1,s2,dist,dist31,dist32,lmn, lmn31,lmn32, memoryV=memoryV, precalc=true, set_maxmin=set_maxmin)
#                                h = 0.0
#                              @inbounds H_thread[o1, o2, cind1, id] += h  * cut

  

                        end
                        end
                        end
                    ############################################
                    if use_threebody_onsite
                    #if false
                        cut2 = array_floats3[counter, 14]
                        o = calc_threebody_onsite(t1,t2,t3,dist,dist31,dist32, database, set_maxmin=set_maxmin, memory=memoryV)
                        for o1 = orb2ind[a1]
#                            a1a,t1,s1 = ind2orb[o1]
                            @inbounds H_thread[o1, o1, c_zero, id] += o  * cut2
                        end
                    elseif !warned_onsite && use_threebody_onsite
                        println("WARNING, missing 3bdy onsite ", (t1, t2, t3))
                        warned_onsite = true
                        within_fit = false
                    end
                    ###########################################

                elseif !warned
                    println("WARNING, missing 3bdy ", (t1, t2, t3))
                    within_fit = false

                    warned = true
                end
            end
        end

#        println("thread sum")
        H += sum(H_thread, dims=4)[:,:,:]



        lmn = zeros(var_type, 3)

        ############ONSITE
        if verbose println("onsite") end
        for c = 1:nkeep_ab
            #        ind_arr[c,:] = R_keep_ab[c][4:6]
            cind = R_keep_ab[c,1]
            a1 = R_keep_ab[c,2]
            a2 = R_keep_ab[c,3]
            t1 = crys.stypes[a1]
            t2 = crys.stypes[a2]

            dist = dist_arr[a1,a2,cind,1]
            lmn[:] = dist_arr[a1,a2,cind,2:4]

            cutoff_onX = get_cutoff(t1,t2)[2]

            if (dist > cutoff_onX)
                continue
            end

            for o1 = orb2ind[a1]
                a1a,t1a,s1 = ind2orb[o1]
                sum1 = summarize_orb(s1)
                for o2 = orb2ind[a1]
                    a2a,t2a,s2 = ind2orb[o2]
                    
                    if dist < 1e-5 #true onsite
                        (h,s) = calc_onsite(t1,s1,s2, database)
                        S[o1, o2, c_zero] += s 
                        H[1, o1, o2, c_zero] += h 
                    else

                        if dist < cutoff_onX - cutoff_length
                            cut = 1.0
                        else
                            cut = cutoff_fn(dist, cutoff_onX - cutoff_length, cutoff_onX)
                        end
                        o = calc_twobody_onsite(t1,t2, s1,s2,dist,lmn, database)
                        H[1, o1, o2, c_zero] += o  * cut
                    end
                    



                end

            end
        end

    end

    #println("make")
    if true
        tb = make_tb(H, ind_arr, S)
        if !ismissing(database) && (haskey(database, "scf") || haskey(database, "SCF"))
            scf = database["scf"]
        else
            scf = false
        end
        tbc = make_tb_crys(tb, crys, nval, 0.0, scf=scf, gamma=gamma, u3=u3, background_charge_correction=background_charge_correction, within_fit=within_fit, screening=screening)
    end
    if verbose 
        println("-----")
        println()
    end

    return tbc

end

function calc_frontier_list(crys_list, frontier=missing)

    if ismissing(frontier)
        frontier = Dict()
    end
    for c in crys_list

        if typeof(c) == crystal{Float64}
            frontier = calc_frontier(c, frontier)
        else
            frontier = calc_frontier(c.crys, frontier)
        end

    end
    return frontier
end


function repel_short_dist_fn(dist, dref, lim)

    if dist >= dref * (1.0 + lim); return 0.0 * dist; end
    
    x = (dref*(1.0 + lim) - dist) / (dref * lim)
    #println("repel_short_dist_fn $dref $lim $dist $x")
    return x^3

end


"""
    function calc_frontier(crys::crystal, frontier; var_type=Float64, test_frontier=false, diststuff=missing, verbose=true)

Calculate a pareto frontier of short distances.
For 2body interactions, this is just a single distance, which is easy.
For 3body interactions, there are 3 distances, so multiple points are on the frontier.
Useful for deciding if old fitting data applies to a new structure with atoms that are close together.

`test_frontier=true` is used to check if new structure is in the frontier.
othersise, see if new structure changes frontier.
"""
function calc_frontier(crys::crystal, frontier; var_type=Float64, test_frontier=false, diststuff=missing, verbose=true, use_threebody=true)

    #println("calc_frontier use_threebody $use_threebody")
    
    if ismissing(var_type)
        var_type=Float64
    end

    lim = 0.04
    repval = 0.001
        
    ind2orb, orb2ind, etotal, nval = orbital_index(crys)
    #use_threebody=true

    if ismissing(diststuff)
#        println("distances")
        #        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy(crys,cutoff2X, cutoff3bX,var_type=var_type)
        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel(crys,cutoff2X, cutoff3bX,var_type=var_type, return_floats=false)
    else
        #        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = diststuff
        R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3 = diststuff
        
    end


#    println("size R_keep_ab ", size(R_keep_ab))
#    println("size array_ind3 ", size(array_ind3))

    #twobody is simple
    nkeep=size(R_keep)[1]
    nkeep_ab = size(R_keep_ab)[1]

    twobody_test = true
    threebody_test = true
    violation_list = []

    nowarn = true

    repel_vals = zeros(var_type, crys.nat)
    At = (crys.A)'
    
    
    #println("2bd")
    for c = 1:nkeep_ab

        cind = R_keep_ab[c,1]
        cham = R_keep_ab[c,7]
        a1 = R_keep_ab[c,2]
        a2 = R_keep_ab[c,3]
        t1= crys.stypes[a1]
        t2= crys.stypes[a2]
        #dist = dist_arr[a1,a2,cind,1]


        dist, lmn = get_dist(a1,a2, R_keep_ab[c,4:6], crys, At)

#        println("dist $dist")
        
        if dist > 6.0
            continue
        end
        
        if test_frontier
            if haskey(frontier, (t1,t2))
                if !isa(frontier[(t1,t2)], Number) && !ismissing(frontier[(t1,t2)]) 
#                if typeof(frontier[(t1,t2)]) == coefs
                    
                    if haskey(frontier[(t1,t2)].dist_frontier, (t1,t2))
                        d = frontier[(t1,t2)].dist_frontier[(t1,t2)]
                    else
                        if nowarn
                            println("warning missing frontier")
                            nowarn = false
                        end
                        d = 0.0
                    end
                else
                    d = frontier[(t1,t2)]
                end
                if !ismissing(d) && dist < d && dist > 0.1
                    twobody_test = false

                    if !((t1,t2,dist) in violation_list)
                        push!(violation_list, (t1,t2,dist))
                    end
                end

                if frontier[(t1,t2)].lim[(t1,t2)] > -1e-12
                    lim = frontier[(t1,t2)].lim[(t1,t2)]
                else
                    lim = 0.04
                end

                if frontier[(t1,t2)].repval[(t1,t2)] > -1e-12
                    repval = frontier[(t1,t2)].repval[(t1,t2)]
                else
                    repval = 0.001
                end
                
                if !ismissing(d) && dist <= d*1.005  && dist > 0.1
#                    println("type of repel_vals[a1] ", typeof(repel_vals[a1]))
#                    println("type fo repel ", typeof(repel_short_dist_fn(dist, d, lim) * 0.1))

                    repel_vals[a1] += repel_short_dist_fn(dist, d, lim) * repval
                    repel_vals[a2] += repel_short_dist_fn(dist, d, lim) * repval

#                    println("repel_2 ", repel_vals)
                end
                
            else
                if !( (t1,t2,0) in violation_list)
                    push!(violation_list, (t1,t2,0))
                end
            end
        else
                    
            add = false
            if haskey(frontier, (t1,t2))
                if frontier[(t1,t2)] > dist && dist > 0.1
                    frontier[(t1,t2)] = dist
                end
            elseif dist > 0.1
                frontier[(t1,t2)] = dist
            end
        end

    end
    #threebody

    if use_threebody
        #println("3bd $use_threebody")
        for counter = 1:size(array_ind3)[1]
            a1 = array_ind3[counter,1]
            a2 = array_ind3[counter,2]
            a3 = array_ind3[counter,3]

            #        cind1 = array_ind3[counter,4]
            #        cind2 = array_ind3[counter,5]

            #dist = array_floats3[counter, 1]
            #dist31 = array_floats3[counter, 2]
            #dist32 = array_floats3[counter, 3]
            
            cind1 = array_ind3[counter,4]
            cind2 = array_ind3[counter,5]

            rind1 = R_keep[cind1,2:4]
            rind2 = R_keep[cind2,2:4]
            
            #        rind1 = R_keep_ab[cind1,4:6]
            #        rind2 = R_keep_ab[cind2,4:6]
            #        rind2 = Rind[cind2,1:3]
            
            dist, lmn = get_dist(a1,a2, rind1, crys, At)
            dist31, lmn31 = get_dist(a1,a3, rind2, crys, At)
            dist32, lmn32 = get_dist(a2,a3, -rind1+rind2, crys, At)

            t1 = crys.stypes[a1]
            t2 = crys.stypes[a2]
            t3 = crys.stypes[a3]
            
            if t1 == t2
                tmp1 = min(dist31, dist32)
                tmp2 = max(dist31, dist32)
                
                dist31=tmp1
                dist32=tmp2
            end


            #        if dist < 6.0 && dist31 < 6.0 && dist32 < 6.0
            #            println("dist $t1 $t2 $t3 ", [dist,dist31,dist32])
            #        end

            if test_frontier                 ##############

                if dist > 6.25 || dist31 > 6.25 || dist32 > 6.25
                    continue
                end
#                println("check $t1 $t2 $t3 $dist $dist31 $dist32")
                
                if haskey(frontier, (t1,t2,t3))

                   
                    if !isa(frontier[(t1,t2,t3)], Array) ####&& !ismissing(frontier[(t1,t2)])
                        #                if typeof(frontier[(t1,t2,t3)]) == coefs
                        if ismissing(frontier[(t1,t2,t3)].dist_frontier) || !( (t1,t2,t3) in keys(frontier[(t1,t2,t3)].dist_frontier))
                            continue
                        end
                        vals = frontier[(t1,t2,t3)].dist_frontier[(t1,t2,t3)]
                        #                    vals = frontier[(t1,t2,t3)].dist_frontier
                    else
                        vals = frontier[(t1,t2,t3)]
                    end

                    if frontier[(t1,t2,t3)].lim[(t1,t2,t3)] > -1e-12
                        lim = frontier[(t1,t2,t3)].lim[(t1,t2,t3)]
                    else
                        lim = 0.04
                    end

                    if frontier[(t1,t2,t3)].repval[(t1,t2,t3)] > -1e-12
                        repval = frontier[(t1,t2,t3)].repval[(t1,t2,t3)]
                    else
                        repval = 0.001
                    end

                    
                    vio = true
                    vio_lim = true
                    for f in vals

                        if dist >= f[1]-1e-5 && dist31 >= f[2]-1e-5 && dist32 >= f[3]-1e-5
                            vio = false
                        end
                        if sum(abs.([dist, dist31, dist32] - f)) < 1e-5
                            vio = false
                        end
                        if dist >= f[1] && dist31 >= f[2] && dist32 >= f[3]
                            vio_lim = false
                        end
                        
                    end

#=                    vio = false
                    vio_lim = false
                    for f in vals

                        if dist <= f[1]-1e-5 && dist31 <= f[2]-1e-5 && dist32 <= f[3]-1e-5
                            vio = true
                        end
                        if dist <= f[1]*(1+lim) && dist31 <= f[2]*(1+lim) && dist32 <= f[3]*(1+lim)
                            vio_lim = true
                        end
                    end
=#


                    
                    if vio_lim == true
                        
                        rsum = 10000000.0
                        rvals = zeros(var_type, 3)
                        for f in vals
                            if dist <= f[1] && dist31 <= f[2] && dist32 <= f[3] 
#                                println("$lim dist $dist $dist31 $dist32 " , f, " ", dist <= f[1]*(1+lim) && dist31 <= f[2]*(1+lim) && dist32 <= f[3]*(1+lim), " ", [dist <= f[1]*(1+lim) , dist31 <= f[2]*(1+lim) , dist32 <= f[3]*(1+lim)])


                                rvals_t = [repel_short_dist_fn(dist, f[1], lim),repel_short_dist_fn(dist31, f[2], lim),repel_short_dist_fn(dist32, f[3], lim)]
                                if sum(rvals_t) < rsum
                                    rsum = sum(rvals_t)
                                    rvals[:] = rvals_t[:]
                                end
                            end
                        end
                        repel_vals[a1] += rvals[1] * repval
                        repel_vals[a2] += rvals[1] * repval
                        
                        repel_vals[a1] += rvals[2] * repval
                        repel_vals[a3] += rvals[2] * repval
                        
                        repel_vals[a2] += rvals[3] * repval
                        repel_vals[a3] += rvals[3] * repval

#                        println("repel_3 ", repel_vals)
                        
                    end
                    

                    
                    
                    if vio
                        threebody_test = false
                        
                        need = true
                        for v in violation_list
                            if t1 == v[1] && t2 == v[2] && t3 == v[3] && abs(dist - v[4]) < 1e-5 && abs(dist31 - v[5]) < 1e-5 && abs(dist32 -v[6]) < 1e-5
                                need = false
                                break
                            end
                        end
                        if need
                            push!(violation_list, (t1,t2,t3,dist, dist31, dist32))
                        end
                        

                        
                    end
                else
                    if !( (t1,t2,t3,0,0,0) in violation_list)
                        push!(violation_list, (t1,t2,t3,0,0,0))
                    end
                end

            else #add to frontier if necessary #############

                if haskey(frontier, (t1,t2,t3))
                    need = true
                    
                    for f in frontier[(t1,t2,t3)]
                        if dist >= f[1]-1e-5 && dist31 >= f[2]-1e-5 && dist32 >= f[3]-1e-5
                            need = false
                        end
                        if sum(abs.([dist, dist31, dist32] - f)) < 1e-5
                            need = false
                        end
                    end
                    
                    if need
                        push!(frontier[(t1,t2,t3)], [dist, dist31,dist32])
                        for iter = 1:5
                            del = -1
                            for (c1,f1) in enumerate(frontier[(t1,t2,t3)])
                                need2 = true
                                for (c2,f2) in enumerate(frontier[(t1,t2,t3)])
                                    if c1 != c2
                                        
                                        if f2[1]-1e-5 < f1[1] && f2[2]-1e-5 < f1[2] && f2[3]-1e-5 < f1[3] 
                                            need2 = false
                                        end
                                        if sum(abs.(f1-f2)) < 1e-5
                                            need2 = false
                                        end
                                    end
                                end
                                if need2 == false
                                    del = c1
                                    break
                                end
                            end
                            if del > 0
                                deleteat!(frontier[(t1,t2,t3)], del)
                            end
                        end
                    end
                else
                    frontier[(t1,t2,t3)] = [[dist, dist31,dist32]]
                end
            end
        end
    end
    
#    println("repel_vals , ", repel_vals)
    
    if test_frontier   
        if twobody_test && threebody_test
            if verbose println("CHECK FRONTIER - everything fine") end
        else
            if var_type == Float64
                println("CHECK FRONTIER WARNING- twobody $twobody_test threebody $threebody_test")
                for v in violation_list
                    println(v)
                end
            end

        end
        return violation_list, twobody_test && threebody_test, repel_vals
    else
        return frontier
    end

end


########################################################################################################################################################
#old slow version
#=
function calc_tb(crys::crystal, database; reference_tbc=missing, verbose=false, var_type=missing, use_threebody=true, use_threebody_onsite=true)

    if ismissing(var_type)
        var_type=Float64
    end

    
    if ismissing(reference_tbc)
        prepare_for_fitting = false
    else
        prepare_for_fitting = true
    end
    

    
    ind2orb, orb2ind, etotal, nval = orbital_index(crys)

    R_keep, R_keep2, dist_arr, c_zero = distances_etc(crys,cutoff2X, cutoff3bX)

#    if ismissing(dist_etc)
#    else
#        (R_keep, R_keep2, dist_arr, c_zero) = dist_etc
#    end
    
        
    c_zero_ref=1
    if !(ismissing(reference_tbc))
        if size(reference_tbc.tb.ind_arr)[1] > 1
            c_zero_ref = reference_tbc.tb.r_dict[[0,0,0]]
        end
    end
    
    nwan = length(keys(ind2orb))

    nkeep=size(R_keep)[1]
    nkeep2=size(R_keep2)[1]    
    println("nkeep, $nkeep, nkeep2, $nkeep2")
    

    H = zeros(var_type, nwan, nwan, nkeep)
    S = zeros(var_type, nwan, nwan, nkeep)    


    ind_arr = zeros(Int64, nkeep, 3)

    lmn = zeros(3)
    dist = 0.0
    lmn31 = zeros(3)
    dist31 = 0.0

    lmn32 = zeros(3)
    dist32 = 0.0

    lmn41 = zeros(3)
    dist14 = 0.0
    dist43 = 0.0

    
    warned = false
    warned_onsite = false

    for c = 1:nkeep
        ind_arr[c,:] = R_keep[c][2:4]
        cind = R_keep[c][1]
        for o1 = 1:nwan
            a1,t1,s1 = ind2orb[o1]
            for o2 = 1:nwan
                a2, t2,s2 = ind2orb[o2]
                
                dist = dist_arr[a1,a2,cind,1]

                if (dist > cutoff2X || dist < 1e-5)
                    continue
                end
                
                lmn[:] = dist_arr[a1,a2,cind,2:4]

                

                (h,s) = calc_twobody(t1,t2,s1,s2,dist,lmn, database)

                if dist < cutoff2X - cutoff_length
                    cut = 1.0
                else
                    cut = cutoff_fn(dist, cutoff2X - cutoff_length, cutoff2X)
                end

                H[o1, o2, c] += h  *cut
                S[o1, o2, c] += s  *cut

                if use_threebody && dist < cutoff2X  && dist > 1e-4
                    for a3 = 1:crys.nat
                        t3 = crys.types[a3]
                        if haskey(database, (t1, t2, t3))
                            for c3 = 1:nkeep2


                                cind3 = R_keep2[c3][1]

                                dist31 = dist_arr[a1,a3,cind3,1]
                                lmn31[:] = dist_arr[a1,a3,cind3,2:4]
                                
                                dist32 = sum( (dist*lmn[:] - dist31*lmn31).^2)^0.5
                                
                                      
                                #no onsite terms or beyond cutoff terms
                                if dist31 > cutoff3bX || dist32 > cutoff3bX || dist > cutoff2X || dist < 1e-4 || dist31 < 1e-4 || dist32 < 1e-4
                                    continue
                                end
                                
                                #cutoffs
                                if dist31 < cutoff3bX - cutoff_length
                                    cut31 = 1.0
                                else
                                    cut31 = cutoff_fn(dist31, cutoff3bX - cutoff_length, cutoff3bX)
                                end

                                if dist32 < cutoff3bX - cutoff_length
                                    cut32 = 1.0
                                else
                                    cut32 = cutoff_fn(dist32, cutoff3bX - cutoff_length, cutoff3bX)
                                end
                                #
                                      
                                lmn32[:] = -(dist*lmn[:] .- dist31*lmn31) ./ dist32
                                
                                h = calc_threebody(t1,t2,t3,s1,s2,dist,dist31,dist32,lmn, lmn31,lmn32, database)


                                H[o1, o2, c] += h  *cut31 * cut32 * cut

                            end
                        elseif !warned
                            println("WARNING, missing 3bdy ", (t1, t2, t3))
                            warned = true
                        end
                    end
                end #end 3bdy

            end
        end
    end

############ONSITE

    for a1 = 1:crys.nat
        orbind = orb2ind[a1]
#        println("orbind ", orbind, " a1 ", a1)

        for o1 = orbind
            a1a,t1,s1 = ind2orb[o1]
            for o2 = orbind
                a2a, t2,s2 = ind2orb[o2]

                if a1a != a2a || a1a != a1 # some checks
                    println("ERROR orb2ind ", [a1a, a2a, a1])
                    continue
                end

                for a3 = 1:crys.nat
                    t3 = crys.types[a3]
                
                    for c = 1:nkeep

                        cind = R_keep[c][1]
                        dist = dist_arr[a1,a3,cind,1]

                        if (dist > cutoff_onX)
                            continue
                        end

                                #


                    
                        lmn[:] = dist_arr[a1,a3,cind,2:4]
                        
                        if dist < 1e-5 #true onsite
                            (h,s) = calc_onsite(t1,s1,s2)
                            S[o1, o2, c_zero] += s 
                            H[o1, o2, c_zero] += h 
                        else

                            if dist < cutoff_onX - cutoff_length
                                cut = 1.0
                            else
                                cut = cutoff_fn(dist, cutoff_onX - cutoff_length, cutoff_onX)
                            end
                            o = calc_twobody_onsite(t1,t3, s1,s2,dist,lmn, database)
                            H[o1, o2, c_zero] += o  * cut
                        end
                        
                        if dist > 1e-5 && dist < cutoff3bX && o1 == o2 && use_threebody_onsite
                            
                            for a4 = 1:crys.nat
                                t4 = crys.types[a4]
                                
                                if haskey(database, (t1, t3, t4))

                                    for c4 = 1:nkeep2
                                        cind4 = R_keep2[c4][1]
                                        dist14 = dist_arr[a1,a4,cind4,1]
                                        
                                        lmn41[:] = dist_arr[a1,a4,cind4,2:4]
                                        
                                        dist43 = sum( (dist*lmn[:] - dist14*lmn41).^2)^0.5
                                        
                                        if (dist14 < 1e-5 || dist43 < 1e-5 || dist14 > cutoff3bX || dist43 > cutoff3bX || dist > cutoff3bX)
                                            continue
                                        end

                                        cut13 = cutoff_fn(dist, cutoff3bX - cutoff_length, cutoff3bX)
                                        cut14 = cutoff_fn(dist14, cutoff3bX - cutoff_length, cutoff3bX)
                                        cut43 = cutoff_fn(dist43, cutoff3bX - cutoff_length, cutoff3bX)
                                        
                                        o = calc_threebody_onsite(t1,t3,t4,s1,dist,dist14,dist43, database)
                                        
                                        H[o1, o1, c_zero] += o  * (cut13 * cut14 * cut43)
                                    end
                                elseif !warned_onsite
                                    println("WARNING, missing 3bdy onsite ", (t1, t3, t4))
                                    warned_onsite = true
                                end
                            end
                        end



                    end

                end
            end
        end
    end

     
    tb = make_tb(H, ind_arr, S)
    tbc = make_tb_crys(tb, crys, nval, 0.0)
    return tbc

end
=#

################################################################################################################################ppp
"""
    function calc_tb_prepare_fast(reference_tbc::tb_crys; use_threebody=false, use_threebody_onsite=false)

Take a `tbc` from DFT and rearrange it for use in fitting code. Basically set it up
so that it is ready for a least squares linear fit of coefficients.

`    return twobody_arrays, threebody_arrays, hvec, svec, Rvec, INDvec, h_onsite, ind_conversion, dmin_types, dmin_types3`

Where
- `twobody_arrays` - info for fitting twobody coefficients
- `threebody_array` - info for fitting twobody coefficients
- `hvec` - vector of reference TB matrix els, arranged for fitting
- `svec` - vector of reference overlap matrix els, arranged for fitting
- `Rvec` - displaments
- `INDvec` info on which matrix el goes with which row
- `h_onsite` info for subtracting atomic terms, I think
- `ind_conversion` - dict to convert between place in hamiltonian and overall counter, which removes duplicates.
- `dmin_types` - shortest 2body distances
- `dmin_types` - shortest 3body distances
"""

function calc_tb_prepare_fast(reference_tbc::tb_crys; use_threebody=false, use_threebody_onsite=false, spin=1, factor_dict = missing, use_eam=false)

    #    println("calc_tb_prepare_fast 3bdy $use_threebody    3bdy-onsite $use_threebody_onsite")
    #    println(reference_tbc.crys)
    #    println()
    
    crys = reference_tbc.crys
    
    var_type=typeof(crys.coords[1,1])


    if ismissing(factor_dict)
        use_factor_dict = false
    else
        use_factor_dict = true
    end
    #    if ismissing(var_type)
#        var_type=Float64
#    end

    #    if ismissing(var_type)
    #        var_type=Float64
    #    end

    
    ind2orb, orb2ind, etotal, nval = orbital_index(crys)

    if use_threebody || use_threebody_onsite
        #println("distances")
        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3,Rind = distances_etc_3bdy_parallel(crys,cutoff2X, cutoff3bX)
        #R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel_LV(crys,cutoff2X, cutoff3bX,var_type=var_type, return_floats=false)
        
    else
        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3,Rind = distances_etc_3bdy_parallel(crys,cutoff2X, 0.0)
        #R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel_LV(crys,cutoff2X, cutoff3bX,var_type=var_type, return_floats=false)        
    end

    #    R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel(crys,cutoff2X, 0.0)

    


    
    #    println("size(array_ind3) ", size(array_ind3))

    #    R_keepP, R_keep_abP, array_ind3P, array_floats3P, dist_arrP, c_zeroP, dmin_typesP, dmin_types3P = distances_etc_3bdy_parallel(crys,cutoff2X, cutoff3bX)

    #    R_keepP, R_keep_abP, array_ind3P, array_floats3P, dist_arrP, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel(crys,cutoff2X, 0.0)

    
    
    
    
    #    println("c_zero $c_zero $c_zeroP")
    #    println(sum(R_keep, dims=1), " R_keep " , sum(R_keepP, dims=1))
    #    println(sum(R_keep_ab, dims=1), " R_keep_ab " , sum(R_keep_abP, dims=1))
    #    println(sum(array_floats3P, dims=1), " array_floats3P ", sum(array_floats3, dims=1))
    #    println(sum(dist_arr), " dist_arr " , sum(dist_arrP))

    #    for key in keys(dmin_types)
    #        println("dmin ", dmin_types[key], " ", dmin_typesP[key])
    #    end

    #    for key in keys(dmin_types3)
    #        println("dmin ", dmin_types3[key], " ", dmin_types3P[key])
    #    end

    #    println("xxxxxxxxxxx")
    
    
    if size(reference_tbc.tb.ind_arr)[1] > 1
        c_zero_ref = reference_tbc.tb.r_dict[[0,0,0]]
    end
    
    nwan = length(keys(ind2orb))

    nkeep=size(R_keep)[1]
    #    nkeep2=size(R_keep2)[1]    
    #    println("nkeep, $nkeep, nkeep2, $nkeep2")
    
    #nkeep_nwan_nwan = size(R_keep_ab)[1]

    twobody_arrays = Dict()

    for c in crys.stypes
        for c2 in crys.stypes
            at_set = Set((c, c2))

            if !haskey(twobody_arrays, at_set)

                coef = make_coefs(at_set, 2)

                #                println("2bdy $at_set")                
                #                println("atset ", at_set)
                #                println("coef")
                #                println(coef)
                #                println()

                hmat = zeros(var_type, nkeep*nwan*nwan, coef.sizeH)
                smat = zeros(var_type, nkeep*nwan*nwan, coef.sizeS)

                #                hmat = spzeros(var_type, nkeep*nwan*nwan, coef.sizeH)
                #                smat = spzeros(var_type, nkeep*nwan*nwan, coef.sizeS)


                twobody_arrays[at_set] = [hmat, smat, coef]

            end
        end
    end

    threebody_arrays = Dict()

    if use_threebody || use_threebody_onsite

        for c in crys.stypes
            for c2 in crys.stypes
                for c3 in crys.stypes
                    at_set = Set((c, c2, c3))
                    if !haskey(threebody_arrays, at_set)
                        #                        println("3bdy $at_set")
                        coef = make_coefs(at_set, 3)
                        hmat = zeros(var_type, nkeep*nwan*nwan, coef.sizeH)
                        #                        hmat = spzeros(var_type, nkeep*nwan*nwan, coef.sizeH)

                        threebody_arrays[at_set] = [hmat, coef]

                    end

                end
            end
        end
    end


    eam_arrays = Dict()

    if use_eam

        for c in crys.stypes
            at_set = Set((c, ))
            if !haskey(eam_arrays, at_set)
                coef = make_coefs(at_set, 0)
                hmat = zeros(var_type, nkeep*nwan*nwan, coef.sizeH)
                eam_arrays[at_set] = [hmat, coef]
            end
            
        end
    end
    


    
    hvec = zeros(var_type, nkeep*nwan*nwan)
    svec = zeros(var_type, nkeep*nwan*nwan)
    dist2 = zeros(var_type,  nkeep*nwan*nwan)


    h_onsite = zeros(var_type, nwan, nwan)

    Rvec = zeros(var_type, nkeep*nwan*nwan,3)
    INDvec = zeros(Int64, nkeep*nwan*nwan,2)


    lmn = zeros(3)
    dist = 0.0
    lmn31 = zeros(3)
    dist31 = 0.0

    lmn32 = zeros(3)
    dist32 = 0.0

    lmn41 = zeros(3)
    dist14 = 0.0

    lmn43 = zeros(3)
    dist43 = 0.0

    
    ind_conversion = Dict()

    counter = 0

    #    println("c_zero $c_zero")

    nkeep_ab = size(R_keep_ab)[1]

    
    rho = zeros(var_type, crys.nat, 3)
    
    println("assign twobody")
    @time for c = 1:nkeep_ab
        #        ind_arr[c,:] = R_keep_ab[c][4:6]
        cind = R_keep_ab[c,1]
        cham = R_keep_ab[c,7]
        a1 = R_keep_ab[c,2]
        a2 = R_keep_ab[c,3]

        if !(R_keep_ab[c,4:6] in keys(reference_tbc.tb.r_dict))
            continue
        end

        c_ref = reference_tbc.tb.r_dict[R_keep_ab[c,4:6]]

        dist = dist_arr[a1,a2,cind,1]
        lmn[:] = dist_arr[a1,a2,cind,2:4]


        for o1 = orb2ind[a1]
            a1a,t1,s1 = ind2orb[o1]
            sum1 = summarize_orb(s1)
            for o2 = orb2ind[a2]
                a2a,t2,s2 = ind2orb[o2]
                cutoff2X = get_cutoff(t1,t2)[1]
                
                sum2 = summarize_orb(s2)
                at_set = Set((t1, t2))
                #                println("$c $o1 $o2 $t1 $t2 $s1 $s2 $sum1 $sum2 $at_set")


                #                if s1 != :s || s2 != :s || dist > 5.0
                #                    continue
                #               end

                if !haskey(ind_conversion, (o1, o2, cham))
                    counter += 1
                    ind_conversion[(o1,o2,cham)] = counter
                    Rvec[counter,:] = R_keep_ab[c,4:6]
                    INDvec[counter,1] = o1
                    INDvec[counter,2] = o2

                end
                
                ind = ind_conversion[(o1,o2,cham)]

                hvec[ind] = real(reference_tbc.tb.H[spin, o1,o2,c_ref])
                if use_factor_dict
                    svec[ind] = real(reference_tbc.tb.S[o1,o2,c_ref]) * (factor_dict[t1]*factor_dict[t2])^0.5
#                    println("apply factor_dict $t1 $t2 ", (factor_dict[t1]*factor_dict[t2])^0.5)
                else
                    svec[ind] = real(reference_tbc.tb.S[o1,o2,c_ref]) 
                end
                dist2[ind] = dist
                #                if abs(svec[ind]) > 0.001
                #                    println(" svec $ind $o1 $o2 $c_ref ",   svec[ind] )
                #                end

                if (dist > cutoff2X || dist < 1e-5)
                    continue

                end

                if dist < cutoff2X - cutoff_length
                    cut = 1.0
                else
                    cut = cutoff_fn(dist, cutoff2X - cutoff_length, cutoff2X)
                end


                (h,s) = fit_twobody(s1,s2,dist,lmn)

                coef = twobody_arrays[at_set][3]

                ih = coef.inds[[t1,sum1,t2,sum2,:H]]
                is = coef.inds[[t1,sum1,t2,sum2,:S]]

                twobody_arrays[at_set][1][ind,ih] += h[:] * cut
                twobody_arrays[at_set][2][ind,is] += s[:] * cut


                
                

                #                if (s1 == :px) && (s2 == :dxz) && dist < 7.0
                #                    println(c," q ", ind,  " " , (t1,s1,t2,s2,:S), " " , is, " " , s[:] * cut)
                #                end

            end
        end
    end

    #    println("counter: ", counter)

    #resize arrays
    Rvec = Rvec[1:counter,:] 
    INDvec = INDvec[1:counter,:]
    hvec = hvec[1:counter]
    svec = svec[1:counter]

    #    for i = 1:counter
    #        if abs(svec[i]) > 0.001
    #            println("SVEC ", Rvec[i,:], " " ,   INDvec[i,:], " ", svec[i])
    #        end
    #    end

    for key in keys(twobody_arrays)
        hmat = twobody_arrays[key][1]
        smat = twobody_arrays[key][2]
        coef = twobody_arrays[key][3]
        hmat = hmat[1:counter,:]
        smat = smat[1:counter,:]

        twobody_arrays[key] = [hmat, smat, coef]
    end

    if use_threebody || use_threebody_onsite
        for key in keys(threebody_arrays)
            hmat = threebody_arrays[key][1] 
            hmat = hmat[1:counter,:]
            coef = threebody_arrays[key][2] 
            threebody_arrays[key] = [hmat, coef]
        end
    end

    ############3body

    println("assign three body")
    @time if use_threebody  || use_threebody_onsite
        for counter = 1:size(array_ind3)[1]
            a1 = array_ind3[counter,1]
            a2 = array_ind3[counter,2]
            a3 = array_ind3[counter,3]

            cind1 = array_ind3[counter,4]
            cind2 = array_ind3[counter,5]

            dist = array_floats3[counter, 1]
            dist31 = array_floats3[counter, 2]
            dist32 = array_floats3[counter, 3]
            
            #            lmn[:] = array_floats3[counter, 4:6]
            lmn31[:] = array_floats3[counter, 4:6]
            lmn32[:] = array_floats3[counter, 7:9]

            cut = array_floats3[counter, 10]
            cut2 = array_floats3[counter, 11]

            t1 = crys.stypes[a1]
            t2 = crys.stypes[a2]
            t3 = crys.stypes[a3]

            at_set3 = Set((t1, t2, t3))

            if use_threebody

                for o1 = orb2ind[a1]
                    a1a,t1a,s1 = ind2orb[o1]
                    #                if t1a != t1
                    #                    println("error t1 t1a $t1 $t1a")
                    #                end
                    sum1 = summarize_orb(s1)

                    for o2 = orb2ind[a2]
                        a2a,t2a,s2 = ind2orb[o2]

                        #                    if t2a != t2
                        #                        println("error t2 t2a $t2 $t2a")
                        #                    end
                        

                        sum2 = summarize_orb(s2)

                        #                    if s1 != :s || s2 != :s
                        #                        continue
                        #                    end


                        if haskey(ind_conversion, (o1, o2, cind1))
                            
                            ind = ind_conversion[(o1,o2,cind1)]
                            
                            h = fit_threebody(t1,t2,t3,s1,s2,dist,dist31,dist32,lmn, lmn31,lmn32)
                            
                            
                            coef = threebody_arrays[at_set3][2]
                            if !([t1, sum1,t2, sum2,t3, :H] in keys(coef.inds))
                                println("err")
                                println((t1, sum1,t2, sum2,t3, :H))
                                println(at_set3)
                            end
                            
                            ih = coef.inds[[t1, sum1,t2, sum2,t3, :H]]



                            #println(at_set3, " zzzzzzzzzz ",(t1, sum1,t2, sum2,t3, :H), " " ,  size(h[:]) , " " , size(threebody_arrays[at_set3][1][ind,ih]))
                            #                            println(ih)
                            #                            println(ind)
                            #                            println(h)
                            #                            println("asdf", at_set3,[t1,t2,t3], size(h), size(ih))
                            threebody_arrays[at_set3][1][ind,ih] += h[1:size(ih)[1]] * cut * 1000
                        end
                        
                    end
                end
            end
            if use_threebody_onsite  #3bdy onsite!!!!!!!!!!!
                


                for o1 = orb2ind[a1]
                    a1a,t1,s1 = ind2orb[o1]
                    sum1 = summarize_orb(s1)
                    h = fit_threebody_onsite(t1,t2,t3,s1,dist,dist31,dist32)

                    #                    if s1 != :s 
                    #                        continue
                    #                    end

                    ind = ind_conversion[(o1,o1,c_zero)]

                    
                    coef = threebody_arrays[at_set3][2]
                    ih = coef.inds[[t1,t2,t3,sum1,:O]]
                    threebody_arrays[at_set3][1][ind,ih] += h[:] * cut2
                    #threebody_arrays[at_set3][1][ind,ih] += [h] * cut2

                end
            end
        end #end 3bdy
    end



    ############
    ############ONSITE

    if true
        for c = 1:nkeep_ab
            #        ind_arr[c,:] = R_keep_ab[c][4:6]
            cind = R_keep_ab[c,1]
            cham = R_keep_ab[c,7]
            a1 = R_keep_ab[c,2]
            a2 = R_keep_ab[c,3]

            t1 = crys.stypes[a1]
            t2 = crys.stypes[a2]
            cutoff_onX = get_cutoff(t1,t2)[2]

            at_set = Set((t1,t2))
            

            dist = dist_arr[a1,a2,cind,1]
            lmn[:] = dist_arr[a1,a2,cind,2:4]
            if (dist > cutoff_onX)
                continue
            end

            
            if dist < cutoff_onX - cutoff_length
                cut = 1.0
            else
                cut = cutoff_fn(dist, cutoff_onX - cutoff_length, cutoff_onX)
            end

            if dist > 1e-5
                ad = EXP_a[1]*dist
                expa=exp.(-0.5*ad)
            
                rho[a1, 1] += (1.0 * expa) * cut
                rho[a1, 2] += (1.0 .- ad) * expa * cut
            #    println("ADD RHO FIT $a1 dist $dist 1 $((1.0 * expa) * cut)   2  $( (1.0 .- ad) * expa * cut)     cut $cut")
            end
            
            for o1 = orb2ind[a1]
                a1a,t1a,s1 = ind2orb[o1]
                sum1 = summarize_orb(s1)
                for o2 = orb2ind[a1]
                    a2a,t2a,s2 = ind2orb[o2]
                    sum2 = summarize_orb(s2)

                    ind = ind_conversion[(o1,o2,c_zero)]


                    if dist < 1e-5 # subtract true onsite from variables to fit

                        #PUREONSITE
                        #                    if s1 == s2
                        #                        coef = twobody_arrays[at_set][3]
                        #                        io = coef.inds[(t1, sum1,:A)][1]
                        #                        twobody_arrays[at_set][1][ind,io] += 1.0
                        #                    end



                        (h,s) = calc_onsite(t1,s1,s2)

                        h_onsite[o1,o2] = h

                        hvec[ind] = hvec[ind] - h
                        svec[ind] = svec[ind] - s
                        #                    if o1 == 1 && o2 == 1
                        #                        println(" onsite $o1 $o2 $a1a $a2a $ind $s1 $s2 $s   ", svec[ind])
                        #                    end
                    else
                        
                        #                    if dist < cutoff_onX - cutoff_length
                        #                        cut = 1.0
                        #                    else
                        #                        cut = cutoff_fn(dist, cutoff_onX - cutoff_length, cutoff_onX)
                        #                    end
                        
                        o = fit_twobody_onsite(t1,t2, s1,s2,dist,lmn)
                        
                        coef = twobody_arrays[at_set][3]

                        io = coef.inds[[t1, sum1,sum2,:O]]
                        twobody_arrays[at_set][1][ind,io] += o[:] * cut

                        #                    println("cath $t1 $t2 $sum1 $sum2 $at_set $t1a $t2a $a1 $a2 $o1 $o2")

                        #                    twobody_arrays[at_set][1][ind,io] += o[:] * cut
                        
                    end

                end
            end
        end
    end

    if use_eam
        for a in 1:crys.nat
            for o = orb2ind[a]
                aa,t,s = ind2orb[o]
                ind = ind_conversion[(o,o,c_zero)]
#                println("t $t ", typeof(t))
                at_set = Set((t,))
#                println(at_set)
#                println("at_set ", at_set)
#                println(keys(eam_arrays))
#                println("RHO $(rho[a,:])    vals $([rho[a,1]^2, rho[a,2]^2, rho[a,1]*rho[a,2]])")
                eam_arrays[at_set][1][ind,:] += [rho[a,1]^2, rho[a,2]^2, rho[a,1]*rho[a,2]]
            end
        end
    end    
    
    ########Check for duplicates. For fitting purposes, we don't need symmetrically equivalent entries. Memory saver.
    println("time duplicates")
    @time if true
        already_found = Dict()
        keep = Int64[]
        nkeep = 0
        ind_conversion = zeros(Int64, length(hvec))
        for i  = 1:length(hvec)

            hval = hvec[i]
            sval = svec[i]
            dval = dist2[i]

            hint = Int64(round(hval*10000000))
            sint = Int64(round(sval*10000000))
            dint = Int64(round(dval*10000000))
            
            sumval = 0.0
            for key in keys(twobody_arrays)
                hmat_new = twobody_arrays[key][1][i,:]
                sumval += sum(hmat_new)
            end
            sumint = Int64(round(sumval*10000000))
            


            #            if i < 100
            #                println("hs $hint $sint $dint")
            #            end

            need = true
            
            if (hint, sint, dint) in  keys(already_found)
                
                counters_old = already_found[(hint, sint, dint, sumint)]
                
                haves = Bool[]
                found_num = -1
                for c in counters_old
                    have = true
                    for key in keys(twobody_arrays)
                        hmat_old = twobody_arrays[key][1][c,:]
                        hmat_new = twobody_arrays[key][1][i,:]
                        
                        if sum(abs.(hmat_old - hmat_new)) > 1e-5
                            have = false
                            break
                        end

                    end
                    
                    if have && (use_threebody || use_threebody_onsite)
                        for key in keys(threebody_arrays)
                            hmat_old = threebody_arrays[key][1][c,:]
                            hmat_new = threebody_arrays[key][1][i,:]
                            if sum(abs.(hmat_old - hmat_new)) > 1e-5
                                #                            println("key", key)#
                                #
                                #                           println(Rvec[c,:], INDvec[c,1], INDvec[c,2])
                                #                            println(Rvec[i,:], INDvec[i,1], INDvec[i,2])

                                #                            println("hmat_old")
                                #                            println(hmat_old)
                                #                            println("hmat_new")
                                #                            println(hmat_new)
                                #                            println()
                                have = false
                                break
                            end
                        end
                    end
                    push!(haves, have)
                    if have == true
                        break
                    end

                end
                #            println("haves $haves")
                if any(haves .== true)
                    need = false
                    for (c,h) in zip(counters_old, haves)
                        if h == true
                            ind_conversion[i] = ind_conversion[c]
                        end
                    end
                else
                    need = true
                    #                println("counters_old $counters_old i $i")
                    #                for c in counters_old
                    #                    println("old ", Rvec[c,:]," ",  INDvec[c,1]," ", INDvec[c,2])
                    #                end
                    #                println("new ", Rvec[i,:]," ",  INDvec[i,1]," ", INDvec[i,2])



                    push!( already_found[(hint, sint, dint, sumint)], i)
                end
                
            else
                need = true
                already_found[(hint, sint, dint, sumint)] = [i]
            end

            if need
                push!(keep, i)
                nkeep += 1 
                ind_conversion[i] = nkeep
            end
            
        end

        #        println("nkeep: ", nkeep)
        
        for key in keys(twobody_arrays)

            th =  deepcopy(twobody_arrays[key][1][keep,:])
            ts =  deepcopy(twobody_arrays[key][2][keep,:])
            coef = deepcopy(twobody_arrays[key][3])
            twobody_arrays[key] = [th, ts, coef]

        end

        for key in keys(threebody_arrays)

            th =  deepcopy(threebody_arrays[key][1][keep,:])
            coef = deepcopy(threebody_arrays[key][2])

            threebody_arrays[key] = [th, coef]

        end

        for key in keys(eam_arrays)

            th =  deepcopy(eam_arrays[key][1][keep,:])
            coef = deepcopy(eam_arrays[key][2])

            eam_arrays[key] = [th, coef]

        end
        
        hvec = hvec[keep]
        svec = svec[keep]
    else
        ind_conversion = collect(1:length(hvec))
    end

    #    for i = 1:counter
    #        if abs(svec[i]) > 0.001
    #            println("SVEC2  ", Rvec[i,:], " " ,   INDvec[i,:], " ", svec[i])
    #        end
    #    end

    #this saves memory during fitting, as these arrays are rather sparse under normal circumstances.
    println("make sparse")
    for k in keys(twobody_arrays)
        twobody_arrays[k][1] = sparse(twobody_arrays[k][1])
        twobody_arrays[k][2] = sparse(twobody_arrays[k][2])
    end
    for k in keys(threebody_arrays)
        threebody_arrays[k][1] = sparse(threebody_arrays[k][1])
    end
    for k in keys(eam_arrays)
        eam_arrays[k][1] = sparse(eam_arrays[k][1])
    end


    #    println("test correlation")
    #    for key in keys(threebody_arrays)
    #        println("key ", key)
    #        X = threebody_arrays[key][1]
    #        cols = size(X)[2]
    #        cmat = cor(X,dims=2)
    #        for c in 1:cols
    #            for c2 in (c+1):cols
    #                if cmat[c,c2] > 0.9999 || cmat[c,c2] < -0.9999
    #                    println("$c $c2 ", cmat[c,c2])
    #                end
    #            end
    #        end
    #    end

    return twobody_arrays, threebody_arrays, hvec, svec, Rvec, INDvec, h_onsite, ind_conversion, dmin_types, dmin_types3, eam_arrays

end
################################################################################################################################ppp
#old slow version
#=
function calc_tb_prepare(reference_tbc::tb_crys; var_type=missing, use_threebody=false, use_threebody_onsite=false)
    
    crys = reference_tbc.crys

    if ismissing(var_type)
        var_type=Float64
    end
    
    ind2orb, orb2ind, etotal, nval = orbital_index(crys)

    R_keep, R_keep2, dist_arr, c_zero = distances_etc(crys,cutoff2X, cutoff3bX)

    
    if size(reference_tbc.tb.ind_arr)[1] > 1
        c_zero_ref = reference_tbc.tb.r_dict[[0,0,0]]
    end
    
    nwan = length(keys(ind2orb))

    nkeep=size(R_keep)[1]
    nkeep2=size(R_keep2)[1]    
#    println("nkeep, $nkeep, nkeep2, $nkeep2")
    

    twobody_arrays = Dict()

    for c in crys.types
        for c2 in crys.types
            if !haskey(twobody_arrays, (c,c2))

#                at_set = Set((c, c2))

                coef = make_coefs(at_set, 2)

                hmat = zeros(var_type, nwan*nwan*nkeep, coef.sizeH)
                smat = zeros(var_type, nwan*nwan*nkeep, coef.sizeS)

                twobody_arrays[(c, c2)] = (hmat, smat, coef)

            end
        end
    end

    threebody_arrays = Dict()

    if use_threebody || use_threebody_onsite

        for c in crys.types
            for c2 in crys.types
                for c3 in crys.types
                    if !haskey(threebody_arrays, (c,c2,c3))
                        
                        coef = make_coefs([c, c2, c3], 3)
                        hmat = zeros(var_type, nwan*nwan*nkeep, coef.sizeH)
                        threebody_arrays[(c, c2, c3)] = (hmat, coef)

                    end
                end
            end
        end
    end

    hvec = zeros(var_type, nwan*nwan*nkeep)
    svec = zeros(var_type, nwan*nwan*nkeep)

    h_onsite = zeros(var_type, nwan, nwan)

    Rvec = zeros(var_type, nwan*nwan*nkeep,3)
    INDvec = zeros(Int64, nwan*nwan*nkeep,2)


    lmn = zeros(3)
    dist = 0.0
    lmn31 = zeros(3)
    dist31 = 0.0

    lmn32 = zeros(3)
    dist32 = 0.0

    lmn41 = zeros(3)
    dist14 = 0.0

    lmn43 = zeros(3)
    dist43 = 0.0

    
    ind_conversion = Dict()

    counter = 0

#    println("c_zero $c_zero")

    for c = 1:nkeep

        cind = R_keep[c][1]
        for o1 = 1:nwan
            a1,t1,s1 = ind2orb[o1]
            sum1 = summarize_orb(s1)
            for o2 = 1:nwan
                a2, t2,s2 = ind2orb[o2]
                sum2 = summarize_orb(s2)
                if !haskey(ind_conversion, (o1, o2, c))
                    counter += 1
                    ind_conversion[(o1,o2,c)] = counter
                end
                
                ind = ind_conversion[(o1,o2,c)]

                lmn[:] = dist_arr[a1,a2,cind,2:4]

                
                c_ref = reference_tbc.tb.r_dict[R_keep[c][2:4]]
                
                hvec[ind] = real(reference_tbc.tb.H[o1,o2,c_ref])
                svec[ind] = real(reference_tbc.tb.S[o1,o2,c_ref])
                Rvec[ind,:] = R_keep[c][2:4]
                INDvec[ind,1] = o1
                INDvec[ind,2] = o2
               
                dist = dist_arr[a1,a2,cind,1]


                if (dist > cutoff2X || dist < 1e-5)
                    continue

                end

                if dist < cutoff2X - cutoff_length
                    cut = 1.0
                else
                    cut = cutoff_fn(dist, cutoff2X - cutoff_length, cutoff2X)
                end


                (h,s) = fit_twobody(t1,t2,s1,s2,dist,lmn)

                coef = twobody_arrays[(t1, t2)][3]

                ih = coef.inds[(sum1,sum2,:H)]
                is = coef.inds[(sum1,sum2,:S)]

                twobody_arrays[(t1, t2)][1][ind,ih] += h[:] * cut
                twobody_arrays[(t1, t2)][2][ind,is] += s[:] * cut

                if use_threebody && dist < cutoff2X  && dist > 1e-4
                    for c3 = 1:nkeep2
                        cind3 = R_keep2[c3][1]
                        for a3 = 1:crys.nat

                            t3 = crys.types[a3]
                            
                            dist31 = dist_arr[a1,a3,cind3,1]
                            lmn31[:] = dist_arr[a1,a3,cind3,2:4]
                                
                            dist32 = sum( (dist*lmn[:] - dist31*lmn31).^2)^0.5

                                
                            #no onsite terms or beyond cutoff terms
                            if dist31 > cutoff3bX || dist32 > cutoff3bX || dist > cutoff2X || dist < 1e-4 || dist31 < 1e-4 || dist32 < 1e-4
                                continue
                            end
                            #cutoffs
                            if dist31 < cutoff3bX - cutoff_length
                                cut31 = 1.0
                            else
                                cut31 = cutoff_fn(dist31, cutoff3bX - cutoff_length, cutoff3bX)
                            end
                            
                            if dist32 < cutoff3bX - cutoff_length
                                cut32 = 1.0
                            else
                                cut32 = cutoff_fn(dist32, cutoff3bX - cutoff_length, cutoff3bX)
                            end
                            #

#                            if dist < cutoff3bX - cutoff_length
#                                cut = 1.0
#                            else
#                                cut = cutoff_fn(dist, cutoff3bX - cutoff_length, cutoff3bX)
#                            end
                            
                            
                            lmn32[:] = -(dist*lmn[:] .- dist31*lmn31) ./ dist32
                            
                            h = fit_threebody(t1,t2,t3,s1,s2,dist,dist31,dist32,lmn, lmn31,lmn32)

                            coef = threebody_arrays[(t1, t2, t3)][2]
                            ih = coef.inds[(sum1,sum2,:H)]


                            threebody_arrays[(t1, t2, t3)][1][ind,ih] += h[:] * cut31 * cut32 * cut
                            
                        end
                    end
                end #end 3bdy


            end
        end
    end

############ONSITE

    for a1 = 1:crys.nat
        orbind = orb2ind[a1]
        for o1 = orbind
            a1a,t1,s1 = ind2orb[o1]
            sum1 = summarize_orb(s1)
            for o2 = orbind
                a2a, t2,s2 = ind2orb[o2]
                sum2 = summarize_orb(s2)

                if a1a != a2a || a1a != a1 # some checks
                    println("ERROR orb2ind ", [a1a, a2a, a1])
                    continue
                end

                for a3 = 1:crys.nat
                    t3 = crys.types[a3]
                
                    for c = 1:nkeep

                        cind = R_keep[c][1]
                        dist = dist_arr[a1,a3,cind,1]

                        ind = ind_conversion[(o1,o2,c_zero)]

                        if (dist > cutoff_onX)
                            continue
                        end
                  
                        lmn[:] = dist_arr[a1,a3,cind,2:4]


                        if dist < 1e-5 # subtract true onsite from variables to fit
                            (h,s) = calc_onsite(t1,s1,s2)

                            h_onsite[o1,o2] = h

                            hvec[ind] = hvec[ind] - h
                            svec[ind] = svec[ind] - s
                            if o1 == 1 && o2 == 1
                                println(" onsite $o1 $o2 $a1a $a2a $ind $s1 $s2 $s   ", svec[ind])
                            end
                        else

                            if dist < cutoff_onX - cutoff_length
                                cut = 1.0
                            else
                                cut = cutoff_fn(dist, cutoff_onX - cutoff_length, cutoff_onX)
                            end

                            o = fit_twobody_onsite(t1,t3, s1,s2,dist,lmn)

                            coef = twobody_arrays[(t1, t3)][3]
                            io = coef.inds[(sum1,sum2,:O)]
                            twobody_arrays[(t1, t3)][1][ind,io] += o[:] * cut

                        end

                        #threebody part
                        if dist > 1e-5 && dist < cutoff3bX && o1 == o2 && use_threebody_onsite

                            for a4 = 1:crys.nat
                                t4 = crys.types[a4]
                                for c4 = 1:nkeep2
                                    cind4 = R_keep2[c4][1]
                                    dist14 = dist_arr[a1,a4,cind4,1]
                                    
                                    ind = ind_conversion[(o1,o1,c_zero)]
                                    lmn41[:] = dist_arr[a1,a4,cind4,2:4]
                                    
                                    dist43 = sum( (dist*lmn[:] - dist14*lmn41).^2)^0.5
                                    
                                    if (dist14 < 1e-5 || dist43 < 1e-5 || dist14 > cutoff3bX || dist43 > cutoff3bX || dist > cutoff3bX)
                                        continue
                                    end

                                    cut13 = cutoff_fn(dist, cutoff3bX - cutoff_length, cutoff3bX)
                                    cut14 = cutoff_fn(dist14, cutoff3bX - cutoff_length, cutoff3bX)
                                    cut43 = cutoff_fn(dist43, cutoff3bX - cutoff_length, cutoff3bX)
                            
                                    h = fit_threebody_onsite(t1,t3,t4,s1,dist,dist14,dist43)
                                    
                                    coef = threebody_arrays[(t1, t3, t4)][2]
                                    ih = coef.inds[(sum1,:O)]
                                    threebody_arrays[(t1, t3, t4)][1][ind,ih] += h[:] * (cut14 * cut43 * cut13)
                                    

                                end
                            end
                        end
                    end
                end
            end
        end
    end


#    println("test correlation")
#    for key in keys(threebody_arrays)
#        println("key ", key)
#        X = threebody_arrays[key][1]
#        cols = size(X)[2]
#        cmat = cor(X,dims=2)
#        for c in 1:cols
#            for c2 in (c+1):cols
#                if cmat[c,c2] > 0.9999 || cmat[c,c2] < -0.9999
#                    println("$c $c2 ", cmat[c,c2])
#                end
#            end
#        end
#    end

    return twobody_arrays, threebody_arrays, hvec, svec, Rvec, INDvec, h_onsite

end
=#



"""
    function calc_onsite(t1,s1,s2, database=missing)

Handles atomic matrix els. We do not currently fit onsite matrix els
"""
function calc_onsite(t1,s1,s2, database=missing)

    #S
    if s1 == s2
        S = 1.0
    else
        S = 0.0
    end

    #H
#    println("t1 $t1 s1 $s1 ")
#    println(atoms[t1].eigs)
    if s1 == s2
        if typeof(s1) <: Int
            H = atoms[t1].eigs[s1]
        else
            H = atoms[t1].eigs[summarize_orb(s1)]
        end

#        if !ismissing(database)
#            c=database[(t1,t1)]
#PUREONSITE
#            if (t1, summarize_orb(s1), :A) in keys(c.inds)
#                ind = c.inds[(t1, summarize_orb(s1), :A)] 
#                H += c.datH[ind[1]]
#            end
#        end


    else
        H = 0.0
    end

    return H,S
    
end

function laguerre_fast!(dist, memory; a = EXP_a[1])

#    a=2.0
    ad = a*dist
    expa=exp.(-0.5*ad)
    memory[1] = 1.0 * expa
    memory[2] = (1.0 .- ad) .* expa
    memory[3]= 0.5*(ad.^2 .- 4.0*ad .+ 2) .* expa
    memory[4] = 1.0/6.0*(-ad.^3 .+ 9.0*ad.^2 .- 18.0*ad .+ 6.0) .* expa
    memory[5] = 1.0/24.0*(ad.^4 .- 16.0 * ad.^3 .+ 72.0*ad.^2 .- 96.0*ad .+ 24.0) .* expa
    memory[6] = 1.0/120*(-ad.^5 .+ 25*ad.^4 .- 200 * ad.^3 .+ 600.0*ad.^2 .- 600.0*ad .+ 120.0) .* expa
#    memory[7] = 1.0/720*(ad.^6  .- 36*ad.^5 .+ 450*ad.^4 .- 2400 * ad.^3 .+ 5400.0*ad.^2 .- 4320*ad .+ 720.0) .* expa

    
end

function laguerre_fast_threebdy!(dist_0, dist_a, dist_b, same_atom, triple, memory)

    #a=2.0
    a=EXP_a[1]

    ad_0 = a*dist_0
    expa_0 =exp.(-0.5*ad_0) #* 10.0

    ad_a = a*dist_a
    expa_a =exp.(-0.5*ad_a)

    ad_b = a*dist_b
    expa_b =exp.(-0.5*ad_b)

    exp_ab = expa_a * expa_b
    
    if triple
        memory[1] = exp_ab
        memory[2] = exp_ab * (1 - ad_b)
        memory[3] = (1 - ad_a) * exp_ab
        memory[4] = expa_0 *exp_ab
    elseif same_atom
        memory[1] = exp_ab
        memory[2] = exp_ab* (  (1 - ad_b) + (1 - ad_a))
        memory[3] = expa_0 * exp_ab
        memory[4] = (1 - ad_0)*expa_0 * exp_ab
        memory[5] = (1 - ad_a)*(1 - ad_b)*exp_ab
        memory[6] = expa_0*exp_ab * ( (1 - ad_b) + (1 - ad_a))
    else
        memory[1] = exp_ab
        memory[2] = exp_ab *(1 - ad_b)
        memory[3] = (1 - ad_a)*exp_ab
        memory[4] = expa_0 * exp_ab
        memory[5] = (1 - ad_0)*expa_0 * exp_ab
        memory[6] = (1 - ad_a)*(1 - ad_b)*exp_ab
        memory[7] = expa_0 * (1 - ad_b) * exp_ab
        memory[8] = expa_0 * (1 - ad_a) * exp_ab
    end


end


function laguerre_fast_threebdy_onsite!(dist_0, dist_a, dist_b, same_atom, memory)

    #a=2.0
    a=EXP_a[1]

    ad_0 = a*dist_b
    expa_0 =exp.(-0.5*ad_0) #* 10.0

    ad_a = a*dist_0
    expa_a =exp.(-0.5*ad_a)

    ad_b = a*dist_a
    expa_b =exp.(-0.5*ad_b)

    expa_ab = expa_a * expa_b
    expa_ab0 = expa_ab * expa_0

    if same_atom
        memory[1] = expa_ab0
        memory[2] = expa_ab0*( (1 - ad_a) + (1 - ad_b) )
        memory[3] = expa_ab0 * (1 - ad_0)
        memory[4] = expa_ab 
        memory[5] = expa_ab * ( (1 - ad_a) + (1 - ad_b) )

    else
        memory[1] = expa_ab
        memory[2] = expa_ab0
        #println("memory ", memory[1:2])
    end

end

"""
    function laguerre(dist, ind=missing; nmax=6, memory=missing)

Calculate laguerre polynomials up to order `nmax`
"""
function laguerre(dist, ind=missing; nmax=6, memory=missing)

    #    a=2.0
    a=EXP_a[1]


#    a=3.0
#    a=1.0


    ad = a*dist
    expa=exp.(-0.5*ad)

    if length(dist) > 1
        l_0 = 1.0 * expa
        l_1 = (1.0 .- ad) .* expa
        l_2 = 0.5*(ad.^2 .- 4.0*ad .+ 2) .* expa
        l_3 = 1.0/6.0*(-ad.^3 .+ 9.0*ad.^2 .- 18.0*ad .+ 6.0) .* expa
        l_4 = 1.0/24.0*(ad.^4 .- 16.0 * ad.^3 .+ 72.0*ad.^2 .- 96.0*ad .+ 24.0) .* expa
        l_5 = 1.0/120*(-ad.^5 .+ 25*ad.^4 .- 200 * ad.^3 .+ 600.0*ad.^2 .- 600.0*ad .+ 120.0) .* expa
        l_6 = 1.0/720*(ad.^6  .- 36*ad.^5 .+ 450*ad.^4 .- 2400 * ad.^3 .+ 5400.0*ad.^2 .- 4320*ad .+ 720.0) .* expa
        L = [l_0 l_1 l_2 l_3 l_4 l_5 l_6]

        if !ismissing(ind)
            s = length(ind)
            if isa(dist, Array)
                return (L[:,1:s]*ind)
            else
                return (L[:,1:s]*ind)[1]
            end

        else
            return L
        end
        return L

    end

    if ismissing(memory)
        memory=zeros(typeof(dist),nmax+1)
    end

    memory[1] = 1.0 * expa
    memory[2] = (1.0 .- ad) .* expa  
    if nmax >= 2
        memory[3] = 0.5*(ad.^2 .- 4.0*ad .+ 2) .* expa 
        if nmax >= 3
            memory[4] = 1.0/6.0*(-ad.^3 .+ 9.0*ad.^2 .- 18.0*ad .+ 6.0) .* expa
            if nmax >= 4
                memory[5] = 1.0/24.0*(ad.^4 .- 16.0 * ad.^3 .+ 72.0*ad.^2 .- 96.0*ad .+ 24.0) .* expa  
                if nmax >= 5
                    memory[6] = 1.0/120*(-ad.^5 .+ 25*ad.^4 .- 200 * ad.^3 .+ 600.0*ad.^2 .- 600.0*ad .+ 120.0) .* expa 
                    if nmax >= 6
                        memory[7] = 1.0/720*(ad.^6  .- 36*ad.^5 .+ 450*ad.^4 .- 2400 * ad.^3 .+ 5400.0*ad.^2 .- 4320*ad .+ 720.0) .* expa
                    end
                end
            end
        end
    end
    

#    L = [l_0 l_1 l_2 l_3 l_4 l_5 l_6]
#    L = [l_0 l_1 l_2 l_3 ]

    if !ismissing(ind)
        s = length(ind)
        if isa(dist, Array)
            return (memory[1:s]'*ind)
        else
            return (memory[1:s]'*ind)[1]
        end

    else
        return memory
    end

#    if !ismissing(ind)
#        s = length(ind)
#        if isa(dist, Array)
#            return (L[:,1:s]*ind)
#        else
#            return (L[:,1:s]*ind)[1]
#        end

#    else
#        return L
#    end
end



"""
    function two_body_O(dist, ind=missing)

Two body onsite
"""
function two_body_O(dist, ind=missing)
    if ismissing(ind)
        return laguerre(dist, ind)[1:n_2body_onsite]
    else
        return  laguerre(dist, ind)
    end
end

"""
    function two_body_O(dist, ind=missing)

three body onsite.
"""
function three_body_O(dist1, dist2, dist3, same_atom, ind=missing; memoryV = missing)

    d1 = laguerre(dist1, missing, nmax=1)
    d2 = laguerre(dist2, missing, nmax=1)
    d3 = laguerre(dist3, missing, nmax=1)
        
    

    if same_atom
    
        if  isa(dist1, Array)
            
#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2] d1[:,2].*d2[:,2].*d3[:,2] ]

#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2] d1[:,2].*d2[:,2].*d3[:,2] (d1[:,3].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,3].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,3]  ]
#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2] d1[:,2].*d2[:,2].*d3[:,2] (d1[:,3].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,3].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,3] (d1[:,2].*d2[:,1].*d3[:,2] + d1[:,1].*d2[:,2].*d3[:,2]) d1[:,2].*d2[:,2].*d3[:,1]  ]
#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2] d1[:,2].*d2[:,2].*d3[:,2] (d1[:,3].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,3].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,3] (d1[:,2].*d2[:,1].*d3[:,2] + d1[:,1].*d2[:,2].*d3[:,2]) d1[:,2].*d2[:,2].*d3[:,1]  ]

            #            V = [10.0*d1[:,1].*d2[:,1].*d3[:,1] 10.0*(d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) 10.0*d1[:,1].*d2[:,1].*d3[:,2] d1[:,1].*d2[:,1] (d1[:,1].*d2[:,2] + d1[:,2].*d2[:,1]) ]
            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2] d1[:,1].*d2[:,1] (d1[:,1].*d2[:,2] + d1[:,2].*d2[:,1]) ]


#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2]  ]
#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1] + d1[:,1].*d2[:,1].*d3[:,2]) ]


        else

            if ismissing(memoryV)
                V = zeros(typeof(d1[1]), n_3body_onsite_same) 
            else
                V = memoryV
            end
            
            #V[1:4] .= [d1[1].*d2[1].*d3[1], (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]), d1[1].*d2[1].*d3[2], d1[2].*d2[2].*d3[2]]
#            V[1:8] .= [d1[1].*d2[1].*d3[1], (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]), d1[1].*d2[1].*d3[2], d1[2].*d2[2].*d3[2], (d1[3].*d2[1].*d3[1] + d1[1].*d2[3].*d3[1]), d1[1].*d2[1].*d3[3], (d1[2].*d2[1].*d3[2] + d1[1].*d2[2].*d3[2]), d1[2].*d2[2].*d3[1]]
#            V[1:8] .= [d1[1].*d2[1].*d3[1], (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]), d1[1].*d2[1].*d3[2], d1[2].*d2[2].*d3[2], (d1[3].*d2[1].*d3[1] + d1[1].*d2[3].*d3[1]), d1[1].*d2[1].*d3[3], (d1[2].*d2[1].*d3[2] + d1[1].*d2[2].*d3[2]), d1[2].*d2[2].*d3[1]]

            #V[1:n_3body_onsite_same] = [d1[1].*d2[1].*d3[1], (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]), d1[1].*d2[1].*d3[2], d1[1].*d2[1],  (d1[1].*d2[2] + d1[2].*d2[1]) ] 

#            V[1] = d1[1].*d2[1].*d3[1] *10.0
#            V[2] = (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1])*10.0
#            V[3] = d1[1].*d2[1].*d3[2]*10.0
#            V[4] = d1[1].*d2[1]
#            V[5] = (d1[1].*d2[2] + d1[2].*d2[1])

            V[1] = d1[1].*d2[1].*d3[1] 
            V[2] = (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1])
            V[3] = d1[1].*d2[1].*d3[2]
            V[4] = d1[1].*d2[1]
            V[5] = (d1[1].*d2[2] + d1[2].*d2[1])
            

            
#            V = [d1[1].*d2[1].*d3[1] (d1[1].*d2[1].*d3[2]+d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1])]


        end
    else
        if  isa(dist1, Array)
#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2]  ]

#            V = [d1[:,1].*d2[:,1].*d3[:,1]  (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]+ d1[:,1].*d2[:,1].*d3[:,2]) ]

            #            V = [d1[:,1].*d2[:,1] d1[:,1].*d2[:,1].*d3[:,1]*10.0]
                        V = [d1[:,1].*d2[:,1] d1[:,1].*d2[:,1].*d3[:,1]]
#            V = [d1[:,1].*d2[:,1].*d3[:,1] d1[:,2].*d2[:,1].*d3[:,1] d1[:,1].*d2[:,2].*d3[:,1] d1[:,1].*d2[:,1].*d3[:,2]]

        else

#            V = [d1[1].*d2[1].*d3[1] (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]) d1[1].*d2[1].*d3[2]]
#            V = d1[1].*d2[1].*d3[1] #(d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]) d1[1].*d2[1].*d3[2]]

#            V = [d1[1].*d2[1].*d3[1] (d1[1].*d2[1].*d3[2]+d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1])]

            if ismissing(memoryV)
                V = zeros(typeof(d1[1]), n_3body_onsite) 
            else
                V = memoryV
            end

#            V[1:4] .= [d1[1].*d2[1].*d3[1]
            
#            V[1:n_3body_onsite] .= [d1[1].*d2[1].*d3[1],       d1[2].*d2[1].*d3[1] ,      d1[1].*d2[2].*d3[1] ,       d1[1].*d2[1].*d3[2]]
#            V = [d1[1].*d2[1] d1[1].*d2[1].*d3[1]]

            V[1] = d1[1].*d2[1]
            #            V[2] = d1[1].*d2[1].*d3[1]*10.0
            V[2] = d1[1].*d2[1].*d3[1]
            
#            V = [d1[1].*d2[1].*d3[1]       d1[2].*d2[1].*d3[1]       d1[1].*d2[2].*d3[1]        d1[1].*d2[1].*d3[2]]
            
        end
    end

    if !ismissing(ind)
        if  isa(dist1, Array)

            return (V*ind) *10^3
        else
            s=size(ind)[1]
            return ((@view V[1:s])'*ind)[1] * 10^3

        end
    else
        return V * 10^3
    end
end    



function three_body_O_lag(d1,d2,d3, same_atom, ind=missing; memoryV = missing)

#    d1 = laguerre(dist1, missing, nmax=1)
#    d2 = laguerre(dist2, missing, nmax=1)
#    d3 = laguerre(dist3, missing, nmax=1)
        
    if same_atom
    
        if  false
            
#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2] d1[:,2].*d2[:,2].*d3[:,2] ]

#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2] d1[:,2].*d2[:,2].*d3[:,2] (d1[:,3].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,3].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,3]  ]
#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2] d1[:,2].*d2[:,2].*d3[:,2] (d1[:,3].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,3].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,3] (d1[:,2].*d2[:,1].*d3[:,2] + d1[:,1].*d2[:,2].*d3[:,2]) d1[:,2].*d2[:,2].*d3[:,1]  ]
#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2] d1[:,2].*d2[:,2].*d3[:,2] (d1[:,3].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,3].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,3] (d1[:,2].*d2[:,1].*d3[:,2] + d1[:,1].*d2[:,2].*d3[:,2]) d1[:,2].*d2[:,2].*d3[:,1]  ]

            #            V = [d1[:,1].*d2[:,1].*d3[:,1]*10.0 (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1])*10.0 d1[:,1].*d2[:,1].*d3[:,2]*10.0 d1[:,1].*d2[:,1] (d1[:,1].*d2[:,2] + d1[:,2].*d2[:,1]) ]
            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2] d1[:,1].*d2[:,1] (d1[:,1].*d2[:,2] + d1[:,2].*d2[:,1]) ]


#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2]  ]
#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1] + d1[:,1].*d2[:,1].*d3[:,2]) ]


        else

            if ismissing(memoryV)
                V = zeros(typeof(d1[1]), n_3body_onsite_same) 
            else
                V = memoryV
            end
            
            #V[1:4] .= [d1[1].*d2[1].*d3[1], (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]), d1[1].*d2[1].*d3[2], d1[2].*d2[2].*d3[2]]
#            V[1:8] .= [d1[1].*d2[1].*d3[1], (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]), d1[1].*d2[1].*d3[2], d1[2].*d2[2].*d3[2], (d1[3].*d2[1].*d3[1] + d1[1].*d2[3].*d3[1]), d1[1].*d2[1].*d3[3], (d1[2].*d2[1].*d3[2] + d1[1].*d2[2].*d3[2]), d1[2].*d2[2].*d3[1]]
#            V[1:8] .= [d1[1].*d2[1].*d3[1], (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]), d1[1].*d2[1].*d3[2], d1[2].*d2[2].*d3[2], (d1[3].*d2[1].*d3[1] + d1[1].*d2[3].*d3[1]), d1[1].*d2[1].*d3[3], (d1[2].*d2[1].*d3[2] + d1[1].*d2[2].*d3[2]), d1[2].*d2[2].*d3[1]]

            #V[1:n_3body_onsite_same] = [d1[1].*d2[1].*d3[1], (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]), d1[1].*d2[1].*d3[2], d1[1].*d2[1],  (d1[1].*d2[2] + d1[2].*d2[1]) ] 

#            V[1] = d1[1].*d2[1].*d3[1] * 10.0
#            V[2] = (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]) * 10.0
#            V[3] = d1[1].*d2[1].*d3[2] * 10.0
#            V[4] = d1[1].*d2[1]
#            V[5] = (d1[1].*d2[2] + d1[2].*d2[1])

            V[1] = d1[1].*d2[1].*d3[1] 
            V[2] = (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]) 
            V[3] = d1[1].*d2[1].*d3[2] 
            V[4] = d1[1].*d2[1]
            V[5] = (d1[1].*d2[2] + d1[2].*d2[1])
            

            
#            V = [d1[1].*d2[1].*d3[1] (d1[1].*d2[1].*d3[2]+d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1])]


        end
    else
        if  false
#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2]  ]

#            V = [d1[:,1].*d2[:,1].*d3[:,1]  (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]+ d1[:,1].*d2[:,1].*d3[:,2]) ]

            #            V = [d1[:,1].*d2[:,1] d1[:,1].*d2[:,1].*d3[:,1] * 10.0]
            V = [d1[:,1].*d2[:,1] d1[:,1].*d2[:,1].*d3[:,1] ]
#            V = [d1[:,1].*d2[:,1].*d3[:,1] d1[:,2].*d2[:,1].*d3[:,1] d1[:,1].*d2[:,2].*d3[:,1] d1[:,1].*d2[:,1].*d3[:,2]]

        else

#            V = [d1[1].*d2[1].*d3[1] (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]) d1[1].*d2[1].*d3[2]]
#            V = d1[1].*d2[1].*d3[1] #(d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]) d1[1].*d2[1].*d3[2]]

#            V = [d1[1].*d2[1].*d3[1] (d1[1].*d2[1].*d3[2]+d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1])]

            if ismissing(memoryV)
                V = zeros(typeof(d1[1]), n_3body_onsite) 
            else
                V = memoryV
            end

#            V[1:4] .= [d1[1].*d2[1].*d3[1]
            
#            V[1:n_3body_onsite] .= [d1[1].*d2[1].*d3[1],       d1[2].*d2[1].*d3[1] ,      d1[1].*d2[2].*d3[1] ,       d1[1].*d2[1].*d3[2]]
#            V = [d1[1].*d2[1] d1[1].*d2[1].*d3[1]]

            V[1] = d1[1].*d2[1]
            #            V[2] = d1[1].*d2[1].*d3[1] * 10.0
            V[2] = d1[1].*d2[1].*d3[1] 
            
#            V = [d1[1].*d2[1].*d3[1]       d1[2].*d2[1].*d3[1]       d1[1].*d2[2].*d3[1]        d1[1].*d2[1].*d3[2]]
            
        end
    end

    if !ismissing(ind)
        if  false

            return (V*ind) *10^3
        else

#            println("V ", V)
            s=size(ind)[1]
            return ((@view V[1:s])'*ind)[1] * 10^3

        end
    else
        return V * 10^3
    end
end    

#=
    d = (dist1+dist2+dist3)/3.0
#    d = sqrt(dist1^2 + dist2^2 + dist3^2)/3.0
    if ismissing(ind)
        return laguerre(d, missing)[:,1:n_3body_onsite]
    else
        return  laguerre(d, ind)
    end
=#


#=

    d1 = laguerre(dist1, missing)
    d2 = laguerre(dist2, missing)
    d3 = laguerre(dist3, missing)
    
    if same_atom
    
        if  isa(dist1, Array)
            
#            Vt = [d2[:,1] .* d3[:,1:2] d2[:,2] .* d3[:,1:2]]
            Vt = [d1[:,1].*d2[:,1] (d1[:,1].*d2[:,2] + d1[:,2].*d2[:,1]) d1[:,2].*d2[:,2]]
            V = [d3[:,1] .* Vt d3[:,2] .* Vt]
        
#return laguerre(dist, ind)[:,1:n_2body_onsite]
        else
#            Vt = [d2[1] .* d3[1:2]' d2[2] .* d3[1:2]']
#            V = [ (d1[1] .* Vt) (d1[2] .* Vt) ]

            Vt = [d1[1].*d2[1] (d1[1].*d2[2] + d1[2].*d2[1]) d1[2].*d2[2]]
            V = [d3[1] .* Vt d3[2] .* Vt]

        end
    else
        if  isa(dist1, Array)

            Vt = [d1[:,1] .* d2[:,1] d1[:,1] .* d2[:,2] d1[:,2] .* d2[:,1] d1[:,2] .* d2[:,2] ]
            V = [d3[:,1] .* Vt d3[:,2] .* Vt]
        
        
        else
            Vt = [d1[1] .* d2[1] d1[1] .* d2[2] d1[2] .* d2[1] d1[2] .* d2[2] ]
            V = [ (d3[1] .* Vt) (d3[2] .* Vt) ]

        end
    end

    if !ismissing(ind)
        if  isa(dist1, Array)

            return (V*ind)
        else
            return (V*ind)[1]
        end
    else
        return V
    end
=#
#end

"""
    function two_body_H(dist, ind=missing)

Two body intersite Hamiltonian.
"""
function two_body_H(dist, ind=missing)
    if ismissing(ind)
        return laguerre(dist, ind, nmax=n_2body-1)[1:n_2body]
    else
        return  laguerre(dist, ind, nmax=n_2body-1)
    end
end

"""
    function two_body_H(dist, ind=missing)

Two body intersite overlap.
"""
function two_body_S(dist, ind=missing)
    if ismissing(ind)
        return laguerre(dist, ind,nmax=n_2body_S-1)[1:n_2body_S]
    else
        return laguerre(dist, ind,nmax=n_2body_S-1)
    end
end

"""
    function two_body_H(dist, ind=missing)

get 3body hamiltonian terms together.
"""
function three_body_H(dist0, dist1, dist2, same_atom, triple, ind=missing; memory0=missing, memory1=missing, memory2=missing, memoryV=missing)

#    return 0.0

    #    zero =  laguerre(dist0,missing, nmax=1, memory=memory0) * 10.0
    zero =  laguerre(dist0,missing, nmax=1, memory=memory0) 
    a = laguerre(dist1,missing, nmax=1, memory=memory1)
    b = laguerre(dist2,missing, nmax=1, memory=memory2)

#    println(same_atom, "three_body_H ", zero, a,b)

#    zero = memory0
#    a = memory1
#    b = memory2
    
    if triple
        if  isa(dist1, Array)
            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]  a[:,2].*b[:,1] zero[:,1].*a[:,1].*b[:,1]]
        else
            if ismissing(memoryV)
                memoryV=zeros(typeof(dist0), max(n_3body, n_3body_same))
            end
            memoryV[1] =  a[1].*b[1]
            memoryV[2] =  a[1].*b[2]
            memoryV[3] =  a[2].*b[1]
            memoryV[4] =  zero[1].*a[1].*b[1]
            
        end

    elseif same_atom
        if  isa(dist1, Array)
#            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  a[:,2].*b[:,2] (a[:,1].*b[:,3]+ a[:,3].*b[:,1])  zero[:,1].*a[:,1].*b[:,1]  zero[:,1].*(a[:,1].*b[:,2]+a[:,2].*b[:,1]) ]
#            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  a[:,2].*b[:,2] (a[:,1].*b[:,3]+ a[:,3].*b[:,1]) zero[:,1].*a[:,1].*b[:,1]   ]
#            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  a[:,2].*b[:,2] (a[:,1].*b[:,3]+ a[:,3].*b[:,1]) zero[:,1].*a[:,1].*b[:,1]      ]
#            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  a[:,2].*b[:,2] (a[:,1].*b[:,3]+ a[:,3].*b[:,1]) zero[:,1].*a[:,1].*b[:,1]   zero[:,1].*(a[:,1].*b[:,2]+a[:,2].*b[:,1]) zero[:,2].*(a[:,1].*b[:,1])   ]

#            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  zero[:,1].*a[:,1].*b[:,1]  ]
            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  zero[:,1].*a[:,1].*b[:,1]  zero[:,2].*a[:,1].*b[:,1]   a[:,2].*b[:,2]  zero[:,1].*(a[:,1].*b[:,2]+ a[:,2].*b[:,1])  ]
            #Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  zero[:,1].*a[:,1].*b[:,1]  zero[:,2].*a[:,1].*b[:,1]    zero[:,1].*(a[:,1].*b[:,2]+ a[:,2].*b[:,1]) ]
#            Vt = [a[:,1].*b[:,1]   zero[:,1].*a[:,1].*b[:,1]  zero[:,2].*a[:,1].*b[:,1]    zero[:,1].*(a[:,1].*b[:,2]+ a[:,2].*b[:,1]) ]
            V = Vt
#            println("V ", V)
        else 
#           try
#                Vt = [a[1].*b[1]  (a[1].*b[2]+a[2].*b[1])  a[2].*b[2] (a[1].*b[3]+ a[3].*b[1])  zero[1].*a[1].*b[1] zero[1].*(a[1].*b[2]+a[2].*b[1]) ]
#                V = Vt
#               println("case ")
                if ismissing(memoryV)
                    memoryV=zeros(typeof(dist0), max(n_3body, n_3body_same))
                end
                memoryV[1] = a[1].*b[1]
                memoryV[2] =  (a[1].*b[2]+a[2].*b[1])
                memoryV[3] =  zero[1].*a[1].*b[1]


               memoryV[4] = zero[2].*a[1].*b[1]   
               memoryV[5] = a[2].*b[2]  
               memoryV[6] = zero[1].*(a[1].*b[2]+ a[2].*b[1]) 
#               memoryV[7] = zero[1].*a[2].*b[2]
#               println("mem ", memoryV)
#                memoryV[3] =   a[2].*b[2] 
#                memoryV[4] =  (a[1].*b[3]+ a[3].*b[1])
#                memoryV[6] = zero[1].*(a[1].*b[2]+a[2].*b[1])
#                memoryV[7] = zero[2].*(a[1].*b[1])


#            catch
#                println("asdf ",size(a), " " , size(b))
#            end

        end
    else
        if  isa(dist1, Array)
#            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]    a[:,2].*b[:,1] a[:,2].*b[:,2]  a[:,1].*b[:,3]    a[:,3].*b[:,1]   zero[:,1].*a[:,1].*b[:,1]   zero[:,1].*a[:,1].*b[:,2]  zero[:,1].*a[:,2].*b[:,1]]

#            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]    a[:,2].*b[:,1] a[:,2].*b[:,2]   zero[:,1].*a[:,1].*b[:,1] a[:,3].*b[:,1] a[:,1].*b[:,3]]
            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]    a[:,2].*b[:,1]    zero[:,1].*a[:,1].*b[:,1]  zero[:,2].*a[:,1].*b[:,1]  a[:,2].*b[:,2]   zero[:,1].*a[:,1].*b[:,2] zero[:,1].*a[:,2].*b[:,1]  ]

#            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]    a[:,2].*b[:,1] a[:,2].*b[:,2]   zero[:,1].*a[:,1].*b[:,1] ]
            V = Vt
        else
#            try
#                Vt = [a[1].*b[1]  a[1].*b[2]   a[2].*b[1]   a[2].*b[2]  a[1].*b[3]    a[3].*b[1]   zero[1].*a[1].*b[1]  zero[1].*a[1].*b[2]  zero[1].*a[2].*b[1] ]
#                V = Vt
                if ismissing(memoryV)
                    memoryV=zeros(typeof(dist0), n_3body)
                end
                memoryV[1] =  a[1].*b[1]
                memoryV[2] =  a[1].*b[2]
                memoryV[3] =  a[2].*b[1]
                memoryV[4] =  zero[1].*a[1].*b[1]

                memoryV[5] =  zero[2].*a[1].*b[1]
                memoryV[6] =  a[2].*b[2]

                memoryV[7] =  zero[1].*a[1].*b[2]
                memoryV[8] =  zero[1].*a[2].*b[1]
                

#                memoryV[6] =  

#                memoryV[6] =  a[3].*b[1]
#                memoryV[7] =  a[1].*b[3]

#                memoryV[7] =  zero[1].*a[1].*b[1]
#                memoryV[8] =  zero[1].*a[1].*b[2]
#                memoryV[9] =  zero[1].*a[2].*b[1]

#            catch
#                println("asdf ",size(a), " " , size(b))
#            end
        end
    end

    if !ismissing(ind)
        if  isa(dist1, Array)

            return (V* (ind * 10^3) )
        else
#            println("three_body_H ", same_atom, " " , size(V), " ", size(ind))
            s=size(ind)[1]
            return (memoryV[1:s]'* (ind*10^3))[1]
        end
    else
        return memoryV #* 10^3
    end

end



function three_body_H_lag(zero,a,b, same_atom, triple, ind=missing; memory0=missing, memory1=missing, memory2=missing, memoryV=missing)

#    return 0.0

#    zero =  laguerre(dist0,missing, nmax=1, memory=memory0)
#    a = laguerre(dist1,missing, nmax=1, memory=memory1)
#    b = laguerre(dist2,missing, nmax=1, memory=memory2)

#    println(same_atom, "three_body_H ", zero, a,b)

#    zero = memory0
#    a = memory1
#    b = memory2

    if triple
        if  false
            #            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]  a[:,2].*b[:,1] zero[:,1].*a[:,1].*b[:,1] * 10.0]
            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]  a[:,2].*b[:,1] zero[:,1].*a[:,1].*b[:,1] ]
        else
            if ismissing(memoryV)
                memoryV=zeros(typeof(zero[1]), n_3body_triple)
            end
            memoryV[1] =  a[1].*b[1]
            memoryV[2] =  a[1].*b[2]
            memoryV[3] =  a[2].*b[1]
            #            memoryV[4] =  zero[1].*a[1].*b[1] * 10.0
            memoryV[4] =  zero[1].*a[1].*b[1] 
        end

        
    elseif same_atom
        if  false
#            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  a[:,2].*b[:,2] (a[:,1].*b[:,3]+ a[:,3].*b[:,1])  zero[:,1].*a[:,1].*b[:,1]  zero[:,1].*(a[:,1].*b[:,2]+a[:,2].*b[:,1]) ]
#            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  a[:,2].*b[:,2] (a[:,1].*b[:,3]+ a[:,3].*b[:,1]) zero[:,1].*a[:,1].*b[:,1]   ]
#            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  a[:,2].*b[:,2] (a[:,1].*b[:,3]+ a[:,3].*b[:,1]) zero[:,1].*a[:,1].*b[:,1]      ]
#            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  a[:,2].*b[:,2] (a[:,1].*b[:,3]+ a[:,3].*b[:,1]) zero[:,1].*a[:,1].*b[:,1]   zero[:,1].*(a[:,1].*b[:,2]+a[:,2].*b[:,1]) zero[:,2].*(a[:,1].*b[:,1])   ]

#            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  zero[:,1].*a[:,1].*b[:,1]  ]
#            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  10.0*zero[:,1].*a[:,1].*b[:,1]  10.0*zero[:,2].*a[:,1].*b[:,1]   a[:,2].*b[:,2]  10.0*zero[:,1].*(a[:,1].*b[:,2]+ a[:,2].*b[:,1])  ]
            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  zero[:,1].*a[:,1].*b[:,1]  zero[:,2].*a[:,1].*b[:,1]   a[:,2].*b[:,2]  zero[:,1].*(a[:,1].*b[:,2]+ a[:,2].*b[:,1])  ]
            #Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  zero[:,1].*a[:,1].*b[:,1]  zero[:,2].*a[:,1].*b[:,1]    zero[:,1].*(a[:,1].*b[:,2]+ a[:,2].*b[:,1]) ]
#            Vt = [a[:,1].*b[:,1]   zero[:,1].*a[:,1].*b[:,1]  zero[:,2].*a[:,1].*b[:,1]    zero[:,1].*(a[:,1].*b[:,2]+ a[:,2].*b[:,1]) ]
            V = Vt
#            println("V ", V)
        else 
#           try
#                Vt = [a[1].*b[1]  (a[1].*b[2]+a[2].*b[1])  a[2].*b[2] (a[1].*b[3]+ a[3].*b[1])  zero[1].*a[1].*b[1] zero[1].*(a[1].*b[2]+a[2].*b[1]) ]
#                V = Vt
#               println("case ")
                if ismissing(memoryV)
                    memoryV=zeros(typeof(zero[1]), max(n_3body, n_3body_same))
                end
                memoryV[1] = a[1].*b[1]
                memoryV[2] =  (a[1].*b[2]+a[2].*b[1])
            #                memoryV[3] =  10.0*zero[1].*a[1].*b[1]
            memoryV[3] =  zero[1].*a[1].*b[1]


            #               memoryV[4] = 10.0*zero[2].*a[1].*b[1]
            memoryV[4] = zero[2].*a[1].*b[1]   
               memoryV[5] = a[2].*b[2]  
            #               memoryV[6] = 10.0*zero[1].*(a[1].*b[2]+ a[2].*b[1])
            memoryV[6] = zero[1].*(a[1].*b[2]+ a[2].*b[1]) 
#               memoryV[7] = zero[1].*a[2].*b[2]
#               println("mem ", memoryV)
#                memoryV[3] =   a[2].*b[2] 
#                memoryV[4] =  (a[1].*b[3]+ a[3].*b[1])
#                memoryV[6] = zero[1].*(a[1].*b[2]+a[2].*b[1])
#                memoryV[7] = zero[2].*(a[1].*b[1])


#            catch
#                println("asdf ",size(a), " " , size(b))
#            end

        end
    else
        if false
#            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]    a[:,2].*b[:,1] a[:,2].*b[:,2]  a[:,1].*b[:,3]    a[:,3].*b[:,1]   zero[:,1].*a[:,1].*b[:,1]   zero[:,1].*a[:,1].*b[:,2]  zero[:,1].*a[:,2].*b[:,1]]

#            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]    a[:,2].*b[:,1] a[:,2].*b[:,2]   zero[:,1].*a[:,1].*b[:,1] a[:,3].*b[:,1] a[:,1].*b[:,3]]
            #            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]    a[:,2].*b[:,1]    10.0*zero[:,1].*a[:,1].*b[:,1]  10.0*zero[:,2].*a[:,1].*b[:,1]  a[:,2].*b[:,2]   10.0*zero[:,1].*a[:,1].*b[:,2] 10.0*zero[:,1].*a[:,2].*b[:,1]  ]
            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]    a[:,2].*b[:,1]    zero[:,1].*a[:,1].*b[:,1]  zero[:,2].*a[:,1].*b[:,1]  a[:,2].*b[:,2]   zero[:,1].*a[:,1].*b[:,2] zero[:,1].*a[:,2].*b[:,1]  ]

#            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]    a[:,2].*b[:,1] a[:,2].*b[:,2]   zero[:,1].*a[:,1].*b[:,1] ]
            V = Vt
        else
#            try
#                Vt = [a[1].*b[1]  a[1].*b[2]   a[2].*b[1]   a[2].*b[2]  a[1].*b[3]    a[3].*b[1]   zero[1].*a[1].*b[1]  zero[1].*a[1].*b[2]  zero[1].*a[2].*b[1] ]
#                V = Vt
                if ismissing(memoryV)
                    memoryV=zeros(typeof(a[1]), n_3body)
                end
                memoryV[1] =  a[1].*b[1]
                memoryV[2] =  a[1].*b[2]
                memoryV[3] =  a[2].*b[1]
            #memoryV[4] =  10.0*zero[1].*a[1].*b[1]
            memoryV[4] =  zero[1].*a[1].*b[1]

            #                memoryV[5] =  10.0*zero[2].*a[1].*b[1]
                            memoryV[5] =  zero[2].*a[1].*b[1]
                memoryV[6] =  a[2].*b[2]

#                memoryV[7] =  10.0*zero[1].*a[1].*b[2]
#                memoryV[8] =  10.0*zero[1].*a[2].*b[1]
                memoryV[7] =  zero[1].*a[1].*b[2]
                memoryV[8] =  zero[1].*a[2].*b[1]
                

#                memoryV[6] =  

#                memoryV[6] =  a[3].*b[1]
#                memoryV[7] =  a[1].*b[3]

#                memoryV[7] =  zero[1].*a[1].*b[1]
#                memoryV[8] =  zero[1].*a[1].*b[2]
#                memoryV[9] =  zero[1].*a[2].*b[1]

#            catch
#                println("asdf ",size(a), " " , size(b))
#            end
        end
    end

    if !ismissing(ind)
        if  false

            return (V* (ind * 10^3) )
        else
#            println("three_body_H ", same_atom, " " , size(V), " ", size(ind))
            s=size(ind)[1]
            return (memoryV[1:s]'* (ind*10^3))[1]
        end
    else
        return memoryV #* 10^3
    end

end


#function two_body(dist, d1,d2,d3,d4)
    #    return    exp(-d1 * dist) * (d2  + dist * d3 + dist^2 * d4)
#    return    exp(-d1 * dist) * (d2  + dist * d3 )    
#end


function symmetry_factor(s1,s2,lmn, dat)
    i1 = orb_num(s1)
    i2 = orb_num(s2)
    return symmetry_factor_int(i1,i2,lmn,dat)
end
    

"""
    function symmetry_factor(s1,s2,lmn, dat)

All of the spd Slater-Koster matrix elements. `dat` is preallocated memory.
"""
function symmetry_factor_int(s1,s2,lmn, dat)
"""
Slater-Koster factors
"""
    
    if s1 == 1 && s2 == 1    #    if s1 == :s && s2 == :s
        return dat[1]

    elseif s1 == 1 && s2 == 3    #    elseif s1 == :s && s2 == :px
        return lmn[1]*dat[1]
    elseif s1 == 1 && s2 == 4    #    elseif s1 == :s && s2 == :py
        return lmn[2]*dat[1]
    elseif s1 == 1 && s2 == 2    #    elseif s1 == :s && s2 == :pz
        return lmn[3]*dat[1]
    elseif s1 == 3 && s2 == 1    #    elseif s1 == :px && s2 == :s
        return -lmn[1]*dat[1]
    elseif s1 == 4 && s2 == 1    #    elseif s1 == :py && s2 == :s
        return -lmn[2]*dat[1]
    elseif s1 == 2 && s2 == 1    #    elseif s1 == :pz && s2 == :s
        return -lmn[3]*dat[1]

    elseif s1 == 3 && s2 == 3    #    elseif s1 == :px && s2 == :px
        return lmn[1]^2*dat[1] + (1.0-lmn[1]^2)*dat[2]
    elseif s1 == 4 && s2 == 4    #    elseif s1 == :py && s2 == :py
        return lmn[2]^2*dat[1] + (1.0-lmn[2]^2)*dat[2]
    elseif s1 == 2 && s2 == 2    #    elseif s1 == :pz && s2 == :pz
        return lmn[3]^2*dat[1] + (1.0-lmn[3]^2)*dat[2]

    elseif s1 == 3 && s2 == 4    #    elseif s1 == :px && s2 == :py
        return lmn[1]*lmn[2]*dat[1] - lmn[1]*lmn[2]*dat[2]
    elseif s1 == 3 && s2 == 2    #    elseif s1 == :px && s2 == :pz
        return lmn[1]*lmn[3]*dat[1] - lmn[1]*lmn[3]*dat[2]
    elseif s1 == 4 && s2 == 2    #    elseif s1 == :py && s2 == :pz
        return lmn[2]*lmn[3]*dat[1] - lmn[2]*lmn[3]*dat[2]
    elseif s1 == 4 && s2 == 3    #    elseif s1 == :py && s2 == :px
        return lmn[2]*lmn[1]*dat[1] - lmn[2]*lmn[1]*dat[2]
    elseif s1 == 2 && s2 == 3    #    elseif s1 == :pz && s2 == :px
        return lmn[3]*lmn[1]*dat[1] - lmn[3]*lmn[1]*dat[2]
    elseif s1 == 2 && s2 == 4    #    elseif s1 == :pz && s2 == :py
        return lmn[3]*lmn[2]*dat[1] - lmn[3]*lmn[2]*dat[2]


#s d
    elseif s1 == 1 && s2 == 5    #    elseif s1 == :s && s2 == :dz2
        return (lmn[3]^2 - (lmn[1]^2 + lmn[2]^2)/2.0) * dat[1]

    elseif s1 == 1 && s2 == 6    #    elseif s1 == :s && s2 == :dxz
        return sqrt3 * lmn[1] * lmn[3] * dat[1]

    elseif s1 == 1 && s2 == 7    #    elseif s1 == :s && s2 == :dyz
        return sqrt3 * lmn[2] * lmn[3] * dat[1]

    elseif s1 == 1 && s2 == 9    #    elseif s1 == :s && s2 == :dxy
        return sqrt3 * lmn[1] * lmn[2] * dat[1]

    elseif s1 == 1 && s2 == 8    #    elseif s1 == :s && s2 == :dx2_y2
        return sqrt3d2 * (lmn[1]^2 - lmn[2]^2) * dat[1]
#d s are same
    elseif s1 == 5 && s2 == 1    #    elseif s1 == :dz2 && s2 == :s
        return (lmn[3]^2 - (lmn[1]^2 + lmn[2]^2)/2.0) * dat[1]

    elseif s1 == 6 && s2 == 1    #    elseif s1 == :dxz && s2 == :s
        return sqrt3 * lmn[1] * lmn[3] * dat[1]

    elseif s1 == 7 && s2 == 1    #    elseif s1 == :dyz && s2 == :s
        return sqrt3 * lmn[2] * lmn[3] * dat[1]

    elseif s1 == 9 && s2 == 1    #    elseif s1 == :dxy && s2 == :s
        return sqrt3 * lmn[1] * lmn[2] * dat[1]

    elseif s1 == 8 && s2 == 1    #    elseif s1 == :dx2_y2 && s2 == :s
        return sqrt3d2 * (lmn[1]^2 - lmn[2]^2) * dat[1]
# p d
    elseif s1 == 3 && s2 == 9    #    elseif s1 == :px && s2 == :dxy
        return sqrt3 * lmn[1]^2 * lmn[2] * dat[1] + lmn[2] * ( 1.0 - 2.0 * lmn[1]^2) * dat[2]

    elseif s1 == 3 && s2 == 7    #    elseif s1 == :px && s2 == :dyz
        return sqrt3 * lmn[1] * lmn[2] * lmn[3] * dat[1] - 2.0 * lmn[1] * lmn[2] * lmn[3] * dat[2]

    elseif s1 == 3 && s2 == 6    #    elseif s1 == :px && s2 == :dxz
        return sqrt3 * lmn[1]^2 * lmn[3] * dat[1] + lmn[3] * (1.0 - 2.0 * lmn[1]^2) * dat[2]

    elseif s1 == 3 && s2 == 8    #    elseif s1 == :px && s2 == :dx2_y2
        return sqrt3d2 * lmn[1] * (lmn[1]^2 - lmn[2]^2) * dat[1] + lmn[1] * (1.0 - lmn[1]^2 + lmn[2]^2) * dat[2]

    elseif s1 == 3 && s2 == 5    #    elseif s1 == :px && s2 == :dz2
        return lmn[1] * (lmn[3]^2 - (lmn[1]^2 + lmn[2]^2)/2.0) * dat[1] -  sqrt3 * lmn[1] * lmn[3]^2 * dat[2]
##
    elseif s1 == 4 && s2 == 9    #    elseif s1 == :py && s2 == :dxy
        return sqrt3 * lmn[2]^2 * lmn[1] * dat[1] + lmn[1] * ( 1.0 - 2.0 * lmn[2]^2) * dat[2]

    elseif s1 == 4 && s2 == 7    #    elseif s1 == :py && s2 == :dyz
        return sqrt3 * lmn[2]^2 * lmn[3] * dat[1] + lmn[3] * ( 1.0 - 2.0 * lmn[2]^2) * dat[2]

    elseif s1 == 4 && s2 == 6    #    elseif s1 == :py && s2 == :dxz
        return sqrt3 * lmn[1] * lmn[2] * lmn[3] * dat[1] - 2.0 * lmn[1] * lmn[2] * lmn[3] * dat[2]


    elseif s1 == 4 && s2 == 8    #    elseif s1 == :py && s2 == :dx2_y2
        return sqrt3d2 * lmn[2] * (lmn[1]^2 - lmn[2]^2) * dat[1] - lmn[2] * (1.0 - lmn[2]^2 + lmn[1]^2) * dat[2]
    elseif s1 == 4 && s2 == 5    #    elseif s1 == :py && s2 == :dz2
        return lmn[2] * (lmn[3]^2 - (lmn[1]^2 + lmn[2]^2)/2.0) * dat[1] -  sqrt3 * lmn[2] * lmn[3]^2 * dat[2]

##
    elseif s1 == 2 && s2 == 9    #    elseif s1 == :pz && s2 == :dxy
        return sqrt3 * lmn[1] * lmn[2] * lmn[3] * dat[1] - 2.0 * lmn[1] * lmn[2] * lmn[3] * dat[2]


    elseif s1 == 2 && s2 == 7    #    elseif s1 == :pz && s2 == :dyz
        return sqrt3 * lmn[3]^2 * lmn[2] * dat[1] + lmn[2] * ( 1.0 - 2.0 * lmn[3]^2) * dat[2]

    elseif s1 == 2 && s2 == 6    #    elseif s1 == :pz && s2 == :dxz
        return sqrt3 * lmn[3]^2 * lmn[1] * dat[1] + lmn[1] * ( 1.0 - 2.0 * lmn[3]^2) * dat[2]

    elseif s1 == 2 && s2 == 8    #    elseif s1 == :pz && s2 == :dx2_y2
        return sqrt3d2 * lmn[3] * (lmn[1]^2 - lmn[2]^2) * dat[1] - lmn[3] * (lmn[1]^2 - lmn[2]^2) * dat[2]
    elseif s1 == 2 && s2 == 5    #    elseif s1 == :pz && s2 == :dz2
        return lmn[3] * (lmn[3]^2 - 0.5*(lmn[1]^2 + lmn[2]^2)) * dat[1] +  sqrt3 * lmn[3] *( lmn[1]^2 +  lmn[2]^2) * dat[2]
## d  p picks up negative signs
##
    elseif  s1 == 9 && s2 == 3     #    elseif  s1 == :dxy && s2 == :px 
        return -sqrt3 * lmn[1]^2 * lmn[2] * dat[1] - lmn[2] * ( 1.0 - 2.0 * lmn[1]^2) * dat[2]

    elseif  s1 == 7 && s2 == 3    #    elseif  s1 == :dyz && s2 == :px
        return -sqrt3 * lmn[1] * lmn[2] * lmn[3] * dat[1] + 2.0 * lmn[1] * lmn[2] * lmn[3] * dat[2]

    elseif  s1 == 6 && s2 == 3    #    elseif  s1 == :dxz && s2 == :px
        return -sqrt3 * lmn[1]^2 * lmn[3] * dat[1] - lmn[3] * (1.0 - 2.0 * lmn[1]^2) * dat[2]

    elseif  s1 == 8 && s2 == 3    #    elseif  s1 == :dx2_y2 && s2 == :px
        return -sqrt3d2 * lmn[1] * (lmn[1]^2 - lmn[2]^2) * dat[1] - lmn[1] * (1.0 - lmn[1]^2 + lmn[2]^2) * dat[2]

    elseif  s1 == 5 && s2 == 3    #    elseif  s1 == :dz2 && s2 == :px
        return -lmn[1] * (lmn[3]^2 - (lmn[1]^2 + lmn[2]^2)/2.0) * dat[1] +  sqrt3 * lmn[1] * lmn[3]^2 * dat[2]
##
    elseif  s1 == 9 && s2 == 4    #    elseif  s1 == :dxy && s2 == :py
        return -sqrt3 * lmn[2]^2 * lmn[1] * dat[1] - lmn[1] * ( 1.0 - 2.0 * lmn[2]^2) * dat[2]

    elseif  s1 == 7 && s2 == 4    #    elseif  s1 == :dyz && s2 == :py
        return -sqrt3 * lmn[2]^2 * lmn[3] * dat[1] - lmn[3] * ( 1.0 - 2.0 * lmn[2]^2) * dat[2]

    elseif  s1 == 6 && s2 == 4    #    elseif  s1 == :dxz && s2 == :py
        return -sqrt3 * lmn[1] * lmn[2] * lmn[3] * dat[1] + 2.0 * lmn[1] * lmn[2] * lmn[3] * dat[2]


    elseif  s1 == 8 && s2 == 4    #    elseif  s1 == :dx2_y2 && s2 == :py
        return -sqrt3d2 * lmn[2] * (lmn[1]^2 - lmn[2]^2) * dat[1] + lmn[2] * (1.0 - lmn[2]^2 + lmn[1]^2) * dat[2]
    elseif  s1 == 5 && s2 == 4    #    elseif  s1 == :dz2 && s2 == :py
        return -lmn[2] * (lmn[3]^2 - (lmn[1]^2 + lmn[2]^2)/2.0) * dat[1] +  sqrt3 * lmn[2] * lmn[3]^2 * dat[2]

##
    elseif  s1 == 9 && s2 == 2    #    elseif  s1 == :dxy && s2 == :pz
        return -sqrt3 * lmn[1] * lmn[2] * lmn[3] * dat[1] + 2.0 * lmn[1] * lmn[2] * lmn[3] * dat[2]


    elseif  s1 == 7 && s2 == 2    #    elseif  s1 == :dyz && s2 == :pz
        return -sqrt3 * lmn[3]^2 * lmn[2] * dat[1] - lmn[2] * ( 1.0 - 2.0 * lmn[3]^2) * dat[2]

    elseif  s1 == 6 && s2 == 2    #    elseif  s1 == :dxz && s2 == :pz
        return -sqrt3 * lmn[3]^2 * lmn[1] * dat[1] - lmn[1] * ( 1.0 - 2.0 * lmn[3]^2) * dat[2]

    elseif  s1 == 8 && s2 == 2    #    elseif  s1 == :dx2_y2 && s2 == :pz
        return -sqrt3d2 * lmn[3] * (lmn[1]^2 - lmn[2]^2) * dat[1] + lmn[3] * (lmn[1]^2 - lmn[2]^2) * dat[2]
    elseif  s1 == 5 && s2 == 2    #    elseif  s1 == :dz2 && s2 == :pz
        return -lmn[3] * (lmn[3]^2 - (lmn[1]^2 + lmn[2]^2)/2.0) * dat[1] -  sqrt3 * lmn[3] *( lmn[1]^2 +  lmn[2]^2) * dat[2]


########## d d

    elseif  (s1 == 9 && s2 == 9)     #    elseif  (s1 == :dxy && s2 == :dxy) 
        return 3.0 * lmn[1]^2 * lmn[2]^2 * dat[1] + (lmn[1]^2 + lmn[2]^2 - 4.0 * lmn[1]^2 * lmn[2]^2) * dat[2] + (lmn[3]^2 + lmn[1]^2 * lmn[2]^2) * dat[3]
    elseif  (s1 == 9 && s2 == 7) || (s2 == 9 && s1 == 7)    #    elseif  (s1 == :dxy && s2 == :dyz) || (s2 == :dxy && s1 == :dyz)
        return 3.0 * lmn[1] * lmn[2]^2 * lmn[3] * dat[1] + lmn[1] * lmn[3] * (1.0 - 4.0 * lmn[2]^2) * dat[2] + lmn[1] * lmn[3] *(lmn[2]^2 - 1.0) * dat[3]
    elseif  (s1 == 9 && s2 == 6) || (s2 == 9 && s1 == 6)    #    elseif  (s1 == :dxy && s2 == :dxz) || (s2 == :dxy && s1 == :dxz)
        return 3.0 * lmn[1]^2 * lmn[2] * lmn[3] * dat[1] + lmn[2] * lmn[3] * (1.0 - 4.0 * lmn[1]^2) * dat[2] + lmn[2] * lmn[3] *(lmn[1]^2 - 1.0) * dat[3]
    elseif  (s1 == 9 && s2 == 8) || (s2 == 9 && s1 == 8)    #    elseif  (s1 == :dxy && s2 == :dx2_y2) || (s2 == :dxy && s1 == :dx2_y2)
        return 3.0/2.0 * lmn[1] * lmn[2] * (lmn[1]^2 - lmn[2]^2)* dat[1] + 2.0 * lmn[1] * lmn[2] *(lmn[2]^2 - lmn[1]^2) * dat[2] + 0.5* lmn[1] * lmn[2] * (lmn[1]^2 - lmn[2]^2) * dat[3]
    elseif  (s1 == 9 && s2 == 5) || (s2 == 9 && s1 == 5)    #    elseif  (s1 == :dxy && s2 == :dz2) || (s2 == :dxy && s1 == :dz2)
        return sqrt3 * lmn[1] * lmn[2] * (lmn[3]^2 - 0.5*(lmn[1]^2 + lmn[2]^2)) * dat[1] - 2.0*sqrt3 * lmn[1] * lmn[2] * lmn[3]^2 * dat[2] + sqrt3d2* lmn[1] * lmn[2] * (1.0 + lmn[3]^2) * dat[3]

    elseif  (s1 == 7 && s2 == 7)     #    elseif  (s1 == :dyz && s2 == :dyz) 
        return  3.0 * lmn[2]^2 * lmn[3]^2 * dat[1] + (lmn[2]^2 + lmn[3]^2 - 4.0 * lmn[2]^2 * lmn[3]^2) * dat[2] + (lmn[1]^2 + lmn[2]^2 * lmn[3]^2) * dat[3]
    elseif  (s1 == 6 && s2 == 6)     #    elseif  (s1 == :dxz && s2 == :dxz) 
        return  3.0 * lmn[1]^2 * lmn[3]^2 * dat[1] + (lmn[1]^2 + lmn[3]^2 - 4.0 * lmn[1]^2 * lmn[3]^2) * dat[2] + (lmn[2]^2 + lmn[1]^2 * lmn[3]^2) * dat[3]
    elseif  (s1 == 6 && s2 == 7) || (s2 == 6 && s1 == 7)    #    elseif  (s1 == :dxz && s2 == :dyz) || (s2 == :dxz && s1 == :dyz)
        return 3.0 * lmn[1] * lmn[3]^2 * lmn[2] * dat[1] + lmn[1] * lmn[2] * (1.0 - 4.0 * lmn[3]^2) * dat[2] + lmn[1] * lmn[2] *(lmn[3]^2 - 1.0) * dat[3]


    elseif  (s1 == 7 && s2 == 8) || (s2 == 7 && s1 == 8)    #    elseif  (s1 == :dyz && s2 == :dx2_y2) || (s2 == :dyz && s1 == :dx2_y2)
        return 3.0/2.0 * lmn[2] * lmn[3] * (lmn[1]^2 - lmn[2]^2)* dat[1] -  lmn[2] * lmn[3] *(1.0 + 2.0* (lmn[1]^2 - lmn[2]^2)) * dat[2] + lmn[2] * lmn[3] * (1.0 + 0.5*(lmn[1]^2 - lmn[2]^2)) * dat[3]
    elseif  (s1 == 7 && s2 == 5) || (s2 == 7 && s1 == 5)    #    elseif  (s1 == :dyz && s2 == :dz2) || (s2 == :dyz && s1 == :dz2)
        return sqrt3 * lmn[2] * lmn[3] * (lmn[3]^2 - 0.5*(lmn[1]^2 + lmn[2]^2)) * dat[1] + sqrt3 * lmn[2] * lmn[3] * (lmn[1]^2 + lmn[2]^2 - lmn[3]^2) * dat[2] - sqrt3d2* lmn[2] * lmn[3] * (lmn[1]^2 + lmn[2]^2) * dat[3]

    elseif  (s1 == 6 && s2 == 8) || (s2 == 6 && s1 == 8)    #    elseif  (s1 == :dxz && s2 == :dx2_y2) || (s2 == :dxz && s1 == :dx2_y2)
        return 3.0/2.0 * lmn[3] * lmn[1] * (lmn[1]^2 - lmn[2]^2)* dat[1] +  lmn[3] * lmn[1] *(1.0 - 2.0* (lmn[1]^2 - lmn[2]^2)) * dat[2] - lmn[3] * lmn[1] * (1.0 - 0.5*(lmn[1]^2 - lmn[2]^2)) * dat[3]


    elseif  (s1 == 6 && s2 == 5) || (s2 == 6 && s1 == 5)    #    elseif  (s1 == :dxz && s2 == :dz2) || (s2 == :dxz && s1 == :dz2)
        return sqrt3 * lmn[1] * lmn[3] * (lmn[3]^2 - 0.5*(lmn[1]^2 + lmn[2]^2)) * dat[1] + sqrt3 * lmn[1] * lmn[3] * (lmn[1]^2 + lmn[2]^2 - lmn[3]^2) * dat[2] - sqrt3d2* lmn[1] * lmn[3] * (lmn[1]^2 + lmn[2]^2) * dat[3]

    elseif  (s1 == 8 && s2 == 8)    #    elseif  (s1 == :dx2_y2 && s2 == :dx2_y2)
        return 0.75*(lmn[1]^2 - lmn[2]^2)^2 * dat[1] + (lmn[1]^2 + lmn[2]^2 - (lmn[1]^2 - lmn[2]^2)^2) * dat[2] + (lmn[3]^2 + 0.25*(lmn[1]^2 - lmn[2]^2)^2) * dat[3]
    elseif  (s1 == 5 && s2 == 5)    #    elseif  (s1 == :dz2 && s2 == :dz2)
        return (lmn[3]^2 - 0.5*(lmn[1]^2 + lmn[2]^2))^2 * dat[1] + 3.0*lmn[3]^2*(lmn[1]^2 + lmn[2]^2) * dat[2] + 0.75 * (lmn[1]^2 + lmn[2]^2)^2 * dat[3]
    elseif  (s1 == 8 && s2 == 5) || (s2 == 8 && s1 == 5)    #    elseif  (s1 == :dx2_y2 && s2 == :dz2) || (s2 == :dx2_y2 && s1 == :dz2)
        return sqrt3d2 * (lmn[1]^2 - lmn[2]^2)*(lmn[3]^2 - 0.5* (lmn[1]^2 + lmn[2]^2)) * dat[1] + sqrt3 * lmn[3]^2 * (lmn[2]^2 - lmn[1]^2) * dat[2] + 0.25*sqrt3*(1 + lmn[3]^2)*(lmn[1]^2 - lmn[2]^2) * dat[3]


    else
        println("symmetry_factor not coded yet: $s1, $s2  !!!!!!!!!!!")
    end
end

"""
    function symmetry_factor_fit(s1,s2,lmn)

All the Slater-Koster factors, for fitting.
"""
function symmetry_factor_fit(s1,s2,lmn)
"""
Slater-Koster factors
"""
    if s1 == :s && s2 == :s
        return [1.0, 0.0, 0.0]

    elseif s1 == :s && s2 == :px
        return [lmn[1], 0.0, 0.0]
    elseif s1 == :s && s2 == :py
        return [lmn[2], 0.0, 0.0]
    elseif s1 == :s && s2 == :pz
        return [lmn[3], 0.0, 0.0]
    elseif s1 == :px && s2 == :s
        return [-lmn[1], 0.0, 0.0]
    elseif s1 == :py && s2 == :s
        return [-lmn[2], 0.0, 0.0]
    elseif s1 == :pz && s2 == :s
        return [-lmn[3], 0.0, 0.0]
    elseif s1 == :px && s2 == :px
        return [lmn[1]^2, (1.0-lmn[1]^2), 0.0]
    elseif s1 == :py && s2 == :py
        return [lmn[2]^2,  (1.0-lmn[2]^2), 0.0]
    elseif s1 == :pz && s2 == :pz
        return [lmn[3]^2,  (1.0-lmn[3]^2), 0.0]
    elseif s1 == :px && s2 == :py
        return [lmn[1]*lmn[2], -lmn[1]*lmn[2], 0.0]
    elseif s1 == :px && s2 == :pz
        return [lmn[1]*lmn[3],  -lmn[1]*lmn[3], 0.0]
    elseif s1 == :py && s2 == :pz
        return [lmn[2]*lmn[3],  -lmn[2]*lmn[3], 0.0]
    elseif s1 == :py && s2 == :px
        return [lmn[2]*lmn[1],  -lmn[2]*lmn[1], 0.0]
    elseif s1 == :pz && s2 == :px
        return [lmn[3]*lmn[1],  -lmn[3]*lmn[1], 0.0]
    elseif s1 == :pz && s2 == :py
        return [lmn[3]*lmn[2],  -lmn[3]*lmn[2], 0.0]
#new
    elseif s1 == :s && s2 == :dz2
        return [(lmn[3]^2 - (lmn[1]^2 + lmn[2]^2)/2.0) , 0.0,0.0]

    elseif s1 == :s && s2 == :dxz
        return [sqrt3 * lmn[1] * lmn[3]  , 0.0, 0.0]

    elseif s1 == :s && s2 == :dyz
        return [sqrt3 * lmn[2] * lmn[3] , 0.0,0.0]

    elseif s1 == :s && s2 == :dxy
        return [sqrt3 * lmn[1] * lmn[2] , 0.0,0.0]

    elseif s1 == :s && s2 == :dx2_y2
        return [sqrt3d2 * (lmn[1]^2 - lmn[2]^2) , 0.0, 0.0]
#d s are same
    elseif s1 == :dz2 && s2 == :s
        return [(lmn[3]^2 - (lmn[1]^2 + lmn[2]^2)/2.0), 0.0,0.0]

    elseif s1 == :dxz && s2 == :s
        return [sqrt3 * lmn[1] * lmn[3] , 0.0, 0.0]

    elseif s1 == :dyz && s2 == :s
        return [sqrt3 * lmn[2] * lmn[3]  , 0.0, 0.0]

    elseif s1 == :dxy && s2 == :s
        return [sqrt3 * lmn[1] * lmn[2]  , 0.0, 0.0]

    elseif s1 == :dx2_y2 && s2 == :s
        return [sqrt3d2 * (lmn[1]^2 - lmn[2]^2)  , 0.0, 0.0]
# p d
    elseif s1 == :px && s2 == :dxy
        return [ sqrt3 * lmn[1]^2 * lmn[2]  ,  lmn[2] * ( 1.0 - 2.0 * lmn[1]^2) , 0.0]

    elseif s1 == :px && s2 == :dyz
        return [ sqrt3 * lmn[1] * lmn[2] * lmn[3] , -2.0 * lmn[1] * lmn[2] * lmn[3] , 0.0]

    elseif s1 == :px && s2 == :dxz
        return [ sqrt3 * lmn[1]^2 * lmn[3] ,  lmn[3] * (1.0 - 2.0 * lmn[1]^2) , 0.0]

    elseif s1 == :px && s2 == :dx2_y2
        return [ sqrt3d2 * lmn[1] * (lmn[1]^2 - lmn[2]^2) ,  lmn[1] * (1.0 - lmn[1]^2 + lmn[2]^2) , 0.0]

    elseif s1 == :px && s2 == :dz2
        return [ lmn[1] * (lmn[3]^2 - (lmn[1]^2 + lmn[2]^2)/2.0) ,  -sqrt3 * lmn[1] * lmn[3]^2 , 0.0]
##
    elseif s1 == :py && s2 == :dxy
        return [ sqrt3 * lmn[2]^2 * lmn[1] ,  lmn[1] * ( 1.0 - 2.0 * lmn[2]^2), 0.0]

    elseif s1 == :py && s2 == :dyz
        return [ sqrt3 * lmn[2]^2 * lmn[3] , lmn[3] * ( 1.0 - 2.0 * lmn[2]^2) , 0.0]

    elseif s1 == :py && s2 == :dxz
        return [ sqrt3 * lmn[1] * lmn[2] * lmn[3] ,  -2.0 * lmn[1] * lmn[2] * lmn[3] , 0.0]


    elseif s1 == :py && s2 == :dx2_y2
        return [ sqrt3d2 * lmn[2] * (lmn[1]^2 - lmn[2]^2) , -lmn[2] * (1.0 - lmn[2]^2 + lmn[1]^2) , 0.0]
    elseif s1 == :py && s2 == :dz2
        return [ lmn[2] * (lmn[3]^2 - (lmn[1]^2 + lmn[2]^2)/2.0) , -sqrt3 * lmn[2] * lmn[3]^2 , 0.0]

##
    elseif s1 == :pz && s2 == :dxy
        return [ sqrt3 * lmn[1] * lmn[2] * lmn[3] , -2.0 * lmn[1] * lmn[2] * lmn[3] , 0.0]


    elseif s1 == :pz && s2 == :dyz
        return [ sqrt3 * lmn[3]^2 * lmn[2] , lmn[2] * ( 1.0 - 2.0 * lmn[3]^2) , 0.0]

    elseif s1 == :pz && s2 == :dxz
        return [ sqrt3 * lmn[3]^2 * lmn[1] , lmn[1] * ( 1.0 - 2.0 * lmn[3]^2) , 0.0]

    elseif s1 == :pz && s2 == :dx2_y2
        return [ sqrt3d2 * lmn[3] * (lmn[1]^2 - lmn[2]^2) , -lmn[3] * (lmn[1]^2 - lmn[2]^2) , 0.0]
    elseif s1 == :pz && s2 == :dz2
        return [ lmn[3] * (lmn[3]^2 - 0.5*(lmn[1]^2 + lmn[2]^2)) ,   sqrt3 * lmn[3] *( lmn[1]^2 +  lmn[2]^2) , 0.0]
## d  p picks up negative signs
##
    elseif  s1 == :dxy && s2 == :px 
        return [ -sqrt3 * lmn[1]^2 * lmn[2] ,  -lmn[2] * ( 1.0 - 2.0 * lmn[1]^2) , 0.0]

    elseif  s1 == :dyz && s2 == :px
        return [ -sqrt3 * lmn[1] * lmn[2] * lmn[3] ,  2.0 * lmn[1] * lmn[2] * lmn[3] , 0.0]

    elseif  s1 == :dxz && s2 == :px
        return [ -sqrt3 * lmn[1]^2 * lmn[3] ,  -lmn[3] * (1.0 - 2.0 * lmn[1]^2) , 0.0]

    elseif  s1 == :dx2_y2 && s2 == :px
        return [ -sqrt3d2 * lmn[1] * (lmn[1]^2 - lmn[2]^2) ,  -lmn[1] * (1.0 - lmn[1]^2 + lmn[2]^2) , 0.0]

    elseif  s1 == :dz2 && s2 == :px
        return [ -lmn[1] * (lmn[3]^2 - (lmn[1]^2 + lmn[2]^2)/2.0) , sqrt3 * lmn[1] * lmn[3]^2 , 0.0]
##
    elseif  s1 == :dxy && s2 == :py
        return [ -sqrt3 * lmn[2]^2 * lmn[1] ,  -lmn[1] * ( 1.0 - 2.0 * lmn[2]^2) , 0.0]

    elseif  s1 == :dyz && s2 == :py
        return [ -sqrt3 * lmn[2]^2 * lmn[3] ,  -lmn[3] * ( 1.0 - 2.0 * lmn[2]^2) , 0.0]

    elseif  s1 == :dxz && s2 == :py
        return [ -sqrt3 * lmn[1] * lmn[2] * lmn[3] ,  2.0 * lmn[1] * lmn[2] * lmn[3] , 0.0]


    elseif  s1 == :dx2_y2 && s2 == :py
        return [ -sqrt3d2 * lmn[2] * (lmn[1]^2 - lmn[2]^2) ,  lmn[2] * (1.0 - lmn[2]^2 + lmn[1]^2) , 0.0]
    elseif  s1 == :dz2 && s2 == :py
        return [ -lmn[2] * (lmn[3]^2 - (lmn[1]^2 + lmn[2]^2)/2.0) ,   sqrt3 * lmn[2] * lmn[3]^2 , 0.0]

##
    elseif  s1 == :dxy && s2 == :pz
        return [ -sqrt3 * lmn[1] * lmn[2] * lmn[3] , 2.0 * lmn[1] * lmn[2] * lmn[3] , 0.0]


    elseif  s1 == :dyz && s2 == :pz
        return [ -sqrt3 * lmn[3]^2 * lmn[2] ,  -lmn[2] * ( 1.0 - 2.0 * lmn[3]^2) , 0.0]

    elseif  s1 == :dxz && s2 == :pz
        return [ -sqrt3 * lmn[3]^2 * lmn[1],  -lmn[1] * ( 1.0 - 2.0 * lmn[3]^2) , 0.0]

    elseif  s1 == :dx2_y2 && s2 == :pz
        return [ -sqrt3d2 * lmn[3] * (lmn[1]^2 - lmn[2]^2) ,  lmn[3] * (lmn[1]^2 - lmn[2]^2) , 0.0]
    elseif  s1 == :dz2 && s2 == :pz
        return [ -lmn[3] * (lmn[3]^2 - (lmn[1]^2 + lmn[2]^2)/2.0) ,  -sqrt3 * lmn[3] *( lmn[1]^2 +  lmn[2]^2) , 0.0]


########## d d

    elseif  (s1 == :dxy && s2 == :dxy) 
        return [ 3.0 * lmn[1]^2 * lmn[2]^2 ,  (lmn[1]^2 + lmn[2]^2 - 4.0 * lmn[1]^2 * lmn[2]^2) ,  (lmn[3]^2 + lmn[1]^2 * lmn[2]^2) ]
    elseif  (s1 == :dxy && s2 == :dyz) || (s2 == :dxy && s1 == :dyz)
        return [ 3.0 * lmn[1] * lmn[2]^2 * lmn[3] , lmn[1] * lmn[3] * (1.0 - 4.0 * lmn[2]^2) ,lmn[1] * lmn[3] *(lmn[2]^2 - 1.0) ]
    elseif  (s1 == :dxy && s2 == :dxz) || (s2 == :dxy && s1 == :dxz)
        return [ 3.0 * lmn[1]^2 * lmn[2] * lmn[3] ,  lmn[2] * lmn[3] * (1.0 - 4.0 * lmn[1]^2) ,  lmn[2] * lmn[3] *(lmn[1]^2 - 1.0) ]
    elseif  (s1 == :dxy && s2 == :dx2_y2) || (s2 == :dxy && s1 == :dx2_y2)
        return [ 3.0/2.0 * lmn[1] * lmn[2] * (lmn[1]^2 - lmn[2]^2) ,  2.0 * lmn[1] * lmn[2] *(lmn[2]^2 - lmn[1]^2) , 0.5*lmn[1] * lmn[2] * (lmn[1]^2 - lmn[2]^2) ]
    elseif  (s1 == :dxy && s2 == :dz2) || (s2 == :dxy && s1 == :dz2)
        return [ sqrt3 * lmn[1] * lmn[2] * (lmn[3]^2 - 0.5*(lmn[1]^2 + lmn[2]^2)) , -2.0*sqrt3 * lmn[1] * lmn[2] * lmn[3]^2 , sqrt3d2*lmn[1] * lmn[2] * (1.0 + lmn[3]^2) ]

    elseif  (s1 == :dyz && s2 == :dyz) 
        return [  3.0 * lmn[2]^2 * lmn[3]^2 ,  (lmn[2]^2 + lmn[3]^2 - 4.0 * lmn[2]^2 * lmn[3]^2) ,  (lmn[1]^2 + lmn[2]^2 * lmn[3]^2) ]
    elseif  (s1 == :dxz && s2 == :dxz) 
        return [  3.0 * lmn[1]^2 * lmn[3]^2 , (lmn[1]^2 + lmn[3]^2 - 4.0 * lmn[1]^2 * lmn[3]^2) , (lmn[2]^2 + lmn[1]^2 * lmn[3]^2) ]
    elseif  (s1 == :dxz && s2 == :dyz) || (s2 == :dxz && s1 == :dyz)
        return [ 3.0 * lmn[1] * lmn[3]^2 * lmn[2] , lmn[1] * lmn[2] * (1.0 - 4.0 * lmn[3]^2) , lmn[1] * lmn[2] *(lmn[3]^2 - 1.0) ]


    elseif  (s1 == :dyz && s2 == :dx2_y2) || (s2 == :dyz && s1 == :dx2_y2)
        return [ 3.0/2.0 * lmn[2] * lmn[3] * (lmn[1]^2 - lmn[2]^2) , -lmn[2] * lmn[3] *(1.0 + 2.0* (lmn[1]^2 - lmn[2]^2)) , lmn[2] * lmn[3] * (1.0 + 0.5*(lmn[1]^2 - lmn[2]^2)) ]
    elseif  (s1 == :dyz && s2 == :dz2) || (s2 == :dyz && s1 == :dz2)
        return [ sqrt3 * lmn[2] * lmn[3] * (lmn[3]^2 - 0.5*(lmn[1]^2 + lmn[2]^2)) , sqrt3 * lmn[2] * lmn[3] * (lmn[1]^2 + lmn[2]^2 - lmn[3]^2) , -sqrt3d2* lmn[2] * lmn[3] * (lmn[1]^2 + lmn[2]^2) ]

    elseif  (s1 == :dxz && s2 == :dx2_y2) || (s2 == :dxz && s1 == :dx2_y2)
        return [ 3.0/2.0 * lmn[3] * lmn[1] * (lmn[1]^2 - lmn[2]^2) ,  lmn[3] * lmn[1] *(1.0 - 2.0* (lmn[1]^2 - lmn[2]^2)) , -lmn[3] * lmn[1] * (1.0 - 0.5*(lmn[1]^2 - lmn[2]^2)) ]


    elseif  (s1 == :dxz && s2 == :dz2) || (s2 == :dxz && s1 == :dz2)
        return [ sqrt3 * lmn[1] * lmn[3] * (lmn[3]^2 - 0.5*(lmn[1]^2 + lmn[2]^2)) , sqrt3 * lmn[1] * lmn[3] * (lmn[1]^2 + lmn[2]^2 - lmn[3]^2) , -sqrt3d2*lmn[1] * lmn[3] * (lmn[1]^2 + lmn[2]^2) ]

    elseif  (s1 == :dx2_y2 && s2 == :dx2_y2)
        return [ 0.75*(lmn[1]^2 - lmn[2]^2)^2 ,  (lmn[1]^2 + lmn[2]^2 - (lmn[1]^2 - lmn[2]^2)^2) , (lmn[3]^2 + 0.25*(lmn[1]^2 - lmn[2]^2)^2) ]
    elseif  (s1 == :dz2 && s2 == :dz2)
        return [ (lmn[3]^2 - 0.5*(lmn[1]^2 + lmn[2]^2))^2 , 3.0*lmn[3]^2*(lmn[1]^2 + lmn[2]^2) , 0.75 * (lmn[1]^2 + lmn[2]^2)^2 ]
    elseif  (s1 == :dx2_y2 && s2 == :dz2) || (s2 == :dx2_y2 && s1 == :dz2)
        return [ sqrt3d2 * (lmn[1]^2 - lmn[2]^2)*(lmn[3]^2 - 0.5* (lmn[1]^2 + lmn[2]^2)) , sqrt3 * lmn[3]^2 * (lmn[2]^2 - lmn[1]^2) , 0.25*sqrt3*(1 + lmn[3]^2)*(lmn[1]^2 - lmn[2]^2) ]



    else
        println("symmetry_factor_fit not coded yet: $s1, $s2  !!!!!!!!!!!")
    end
end



"""
    function calc_twobody(t1,t2,orb1,orb2,dist,lmn, database) 

Two body intersite Hamiltonian and overlap matrix els.
"""
function calc_twobody(t1,t2,orb1,orb2,dist,lmn, database)

    
    c = database[(t1,t2)]
    
    if dist < c.min_dist * 0.95
        dist = c.min_dist * 0.95
    end

    s1=summarize_orb(orb1)
    s2=summarize_orb(orb2)

    indH = c.inds[[t1,s1,t2,s2,:H]]
    indS = c.inds[[t1,s1,t2,s2,:S]]
    
    H1 =  two_body_H(dist, c.datH[indH[1:n_2body]])
    S1 =  two_body_S(dist, c.datS[indS[1:n_2body_S]])

    if (s1 == :p && s2 == :p) 

        H2 =  two_body_H(dist, c.datH[indH[1+n_2body:n_2body*2]])
        S2 =  two_body_S(dist, c.datS[indS[1+n_2body_S:n_2body_S*2]])

        Htot = symmetry_factor(orb1,orb2,lmn, [H1 H2 ])
        Stot = symmetry_factor(orb1,orb2,lmn, [S1 S2 ])

    elseif (s1 == :p && s2 == :d) || ( s2 == :p && s1 == :d)

        H2 =  two_body_H(dist, c.datH[indH[1+n_2body:n_2body*2]])
        S2 =  two_body_S(dist, c.datS[indS[1+n_2body_S:n_2body_S*2]])

        Htot = symmetry_factor(orb1,orb2,lmn, [H1 H2 ])
        Stot = symmetry_factor(orb1,orb2,lmn, [S1 S2 ])

    elseif (s1 == :d && s2 == :d)

        H2 =  two_body_H(dist, c.datH[indH[1+n_2body:n_2body*2]])
        S2 =  two_body_S(dist, c.datS[indS[1+n_2body_S:n_2body_S*2]])

        H3 =  two_body_H(dist, c.datH[indH[1+2*n_2body:n_2body*3]])
        S3 =  two_body_S(dist, c.datS[indS[1+2*n_2body_S:n_2body_S*3]])

        Htot = symmetry_factor(orb1,orb2,lmn, [H1 H2 H3 ])
        Stot = symmetry_factor(orb1,orb2,lmn, [S1 S2 S3 ])
        
    else
        
        Htot = symmetry_factor(orb1,orb2,lmn, [H1 ])
        Stot = symmetry_factor(orb1,orb2,lmn, [S1 ])

    end
    
    return Htot, Stot
    
end


"""
    function calc_twobody(t1,t2,orb1,orb2,dist,lmn, database) 

Two body intersite Hamiltonian and overlap matrix els.
"""
function calc_twobody_faster(t1,t2,orb1,orb2,s1, s2, dist,lmn, c, indH, indS, lag)

#    c = database[(t1,t2)]
    
    if dist < c.min_dist * 0.95
        dist = c.min_dist * 0.95
    end

#    println("faster", [t1,t2,orb1,orb2,s1, s2])
    
#    s1=summarize_orb(orb1)
#    s2=summarize_orb(orb2)

#    indH = c.inds[[t1,s1,t2,s2,:H]]
#    indS = c.inds[[t1,s1,t2,s2,:S]]
    
#    H1 =  two_body_H(dist, c.datH[indH[1:n_2body]])
#    S1 =  two_body_S(dist, c.datS[indS[1:n_2body_S]])

    H1 = sum((@view lag[1:n_2body]) .* c.datH[@view indH[1:n_2body]])
    S1 = sum((@view lag[1:n_2body_S]) .* c.datS[@view indS[1:n_2body_S]])

    
    if (s1 == 2 && s2 == 2)  #:p :p

#        H2 =  two_body_H(dist, c.datH[indH[1+n_2body:n_2body*2]])
#        S2 =  two_body_S(dist, c.datS[indS[1+n_2body_S:n_2body_S*2]])

        H2 = sum((@view lag[1:n_2body]) .* c.datH[@view indH[1+n_2body:n_2body*2]])
        S2 = sum((@view lag[1:n_2body_S]) .* c.datS[@view indS[1+n_2body_S:n_2body_S*2]])

        
        Htot = symmetry_factor_int(orb1,orb2,lmn, [H1 H2 ])
        Stot = symmetry_factor_int(orb1,orb2,lmn, [S1 S2 ])

    elseif (s1 == 2 && s2 == 3) || ( s2 == 2 && s1 == 3) #p d

#        H2 =  two_body_H(dist, c.datH[indH[1+n_2body:n_2body*2]])
#        S2 =  two_body_S(dist, c.datS[indS[1+n_2body_S:n_2body_S*2]])

        H2 = sum((@view lag[1:n_2body]) .* c.datH[@view indH[1+n_2body:n_2body*2]])
        S2 = sum((@view lag[1:n_2body_S]) .* c.datS[@view indS[1+n_2body_S:n_2body_S*2]])

        
        
        Htot = symmetry_factor_int(orb1,orb2,lmn, [H1 H2 ])
        Stot = symmetry_factor_int(orb1,orb2,lmn, [S1 S2 ])

    elseif (s1 == 3 && s2 == 3) #dd

#        H2 =  two_body_H(dist, c.datH[indH[1+n_2body:n_2body*2]])
#        S2 =  two_body_S(dist, c.datS[indS[1+n_2body_S:n_2body_S*2]])

#        H3 =  two_body_H(dist, c.datH[indH[1+2*n_2body:n_2body*3]])
#        S3 =  two_body_S(dist, c.datS[indS[1+2*n_2body_S:n_2body_S*3]])

        H2 = sum((@view lag[1:n_2body]) .* c.datH[@view indH[1+n_2body:n_2body*2]])
        S2 = sum((@view lag[1:n_2body_S]) .* c.datS[@view indS[1+n_2body_S:n_2body_S*2]])

        H3 = sum((@view lag[1:n_2body]) .* c.datH[@view indH[1+n_2body*2:n_2body*3]])
        S3 = sum((@view lag[1:n_2body_S]) .* c.datS[@view indS[1+2*n_2body_S:n_2body_S*3]])
        
        
        Htot = symmetry_factor_int(orb1,orb2,lmn, [H1 H2 H3 ])
        Stot = symmetry_factor_int(orb1,orb2,lmn, [S1 S2 S3 ])
        
    else
        
        Htot = symmetry_factor_int(orb1,orb2,lmn, [H1 ])
        Stot = symmetry_factor_int(orb1,orb2,lmn, [S1 ])

    end
    
    return Htot, Stot
    
end


"""
    function fit_twobody(orb1,orb2,dist,lmn)

Fit Two body intersite Hamiltonian and overlap matrix els.
"""
function fit_twobody(orb1,orb2,dist,lmn)

    
    s1=summarize_orb(orb1)
    s2=summarize_orb(orb2)

    
    H1 =  two_body_H(dist, missing)
    S1 =  two_body_S(dist, missing)

    if (s1 == :p && s2 == :p) 

        H2 =  two_body_H(dist, missing)
        S2 =  two_body_S(dist, missing)

        sym = symmetry_factor_fit(orb1,orb2,lmn)
        
        H = hcat(H1*sym[1], H2*sym[2])
        S = hcat(S1*sym[1], S2*sym[2])

    elseif (s1 == :p && s2 == :d)  || (s2 == :p && s1 == :d)


        H2 =  two_body_H(dist, missing)
        S2 =  two_body_S(dist, missing)

        sym = symmetry_factor_fit(orb1,orb2,lmn)
        
        H = hcat(H1*sym[1], H2*sym[2])
        S = hcat(S1*sym[1], S2*sym[2])

    elseif (s1 == :d && s2 == :d)  

        H2 =  two_body_H(dist, missing)
        S2 =  two_body_S(dist, missing)

        H3 =  two_body_H(dist, missing)
        S3 =  two_body_S(dist, missing)

        sym = symmetry_factor_fit(orb1,orb2,lmn)
        
        H = hcat(H1*sym[1], H2*sym[2], H3*sym[3])
        S = hcat(S1*sym[1], S2*sym[2], S3*sym[3])

        
    else

        sym = symmetry_factor_fit(orb1,orb2,lmn)
        
        H = H1 * sym[1]
        S = S1 * sym[1]


    end
    
    return H, S
    
end


"""
    function calc_twobody_onsite(t1,t2,orb1,orb2, dist,lmn, database) 

Calculate 2body onsite interactions.
"""
function calc_twobody_onsite(t1,t2,orb1,orb2, dist,lmn, database)

    c = database[(t1,t2)]

#    if dist < c.min_dist * 0.95
#        dist = c.min_dist * 0.95
#    end


    o1 = summarize_orb(orb1)
    o2 = summarize_orb(orb2)    

    indO = c.inds[[t1,o1,o2,:O]]


    if orb1 == orb2
        Otot = two_body_O(dist, c.datH[indO[1:n_2body_onsite]])
    else
        Otot = 0.0
    end


    if (o1 == :s && o2 == :p) || (o1 == :p && o2 == :s) || (o1 == :s && o2 == :d) || (o1 == :d && o2 == :s) || (o1 == :p && o2 == :d) || (o1 == :d && o2 == :p)
        O_split = two_body_O(dist, c.datH[indO[1:n_2body_onsite]])
        Otot += O_split * symmetry_factor(orb1,:s,lmn, [1.0 ]) * symmetry_factor(orb2,:s,lmn, [1.0 ])
    elseif (o1 == :p && o2 == :p) || (o1 == :d && o2 == :d)
        O_split = two_body_O(dist, c.datH[indO[1+n_2body_onsite:n_2body_onsite*2]])
        Otot += O_split * symmetry_factor(orb1,:s,lmn, [1.0 ]) * symmetry_factor(orb2,:s,lmn, [1.0 ])

    end

#    if o1 == :s && o2 == :s
#        println("calc_twobody_onsite ",(t1,t2), (t1,o1,o2,:O), " ", orb1, " " , orb2, " " , dist , " lmn ", lmn, " Otot " , Otot, " " , indO, " ", c.datH[indO[1:n_2body_onsite]] )
#    end
    return Otot
        
end

"""
    function fit_twobody_onsite(t1,t2,orb1,orb2, dist,lmn)

Fit 2body onsite interactions.
"""
function fit_twobody_onsite(t1,t2,orb1,orb2, dist,lmn)

    o1 = summarize_orb(orb1)
    o2 = summarize_orb(orb2)    

    if orb1 == orb2
        O = two_body_O(dist)
    else
        O = zeros(n_2body_onsite)
    end


    if (o1 == :s && o2 == :p) || (o1 == :p && o2 == :s) || (o1 == :s && o2 == :d) || (o1 == :d && o2 == :s) || (o1 == :p && o2 == :d) || (o1 == :d && o2 == :p)

        O_split = two_body_O(dist)
        sym = symmetry_factor(orb1,:s,lmn, [1.0 ]) * symmetry_factor(orb2,:s,lmn, [1.0 ])
        O = O_split*sym #only split
#    elseif o1 == :p && o2 == :p
    elseif (o1 == :p && o2 == :p) || (o1 == :d && o2 == :d)
        
        O_split = two_body_O(dist)
        sym = symmetry_factor(orb1,:s,lmn, [1.0 ]) * symmetry_factor(orb2,:s,lmn, [1.0 ])
#        println("so ", size(O), size(O_split))
        O = hcat(O', O_split'*sym) #two terms
    end

    return O
        
end

########
"""
    function calc_threebody_onsite(t1,t2,t3,orb1,dist12,dist13,dist23, database; set_maxmin=false, memory=missing)

calculate threebody onsite interactions
"""
function calc_threebody_onsite(t1,t2,t3,orb1,dist12,dist13,dist23, cdat; set_maxmin=false, memory=missing)



    o1 = summarize_orb(orb1)

#    o2 = summarize_orb(orb2)    

    sameat = (t1 == t2 && t1 == t3 )


    #    sameat = (t2 == t3 || t1 == t2 || t1 == t3)
    #sameat = (t2 == t3)
    #sameat = (t1 == t2 && t1 == t3 )

    indO = cdat.inds[[t1,t2,t3,o1,:O]]
    Otot = three_body_O(dist12, dist13, dist23, sameat, cdat.datH[indO], memoryV=memory)
#    Otot = 0.0
    

    return Otot
        
end


function calc_threebody_onsite_lag(t1,t2,t3,orb1,d1,d2,d3, cdat; set_maxmin=false, memory=missing)



    o1 = summarize_orb(orb1)

#    o2 = summarize_orb(orb2)    

    sameat = (t1 == t2 && t1 == t3 )


    #    sameat = (t2 == t3 || t1 == t2 || t1 == t3)
    #sameat = (t2 == t3)
    #sameat = (t1 == t2 && t1 == t3 )

    indO = cdat.inds[[t1,t2,t3,o1,:O]]
    Otot = three_body_O_lag( d1,d2,d3, sameat, cdat.datH[indO], memoryV=memory)
#    Otot = 0.0
    

    return Otot
        
end

"""
    function fit_threebody_onsite(t1,t2,t3,orb1,dist12,dist13,dist23) 

Fit three body onsite interactions.
"""
function fit_threebody_onsite(t1,t2,t3,orb1,dist12,dist13,dist23)

    o1 = summarize_orb(orb1)
#    o2 = summarize_orb(orb2)    

#    println("f3o $t1 $t2 $t3 $orb1  ", sum(dist12-dist13), " ", sum(dist12-dist23), " ",sum(dist13-dist23))

#    sameat = (t2 == t3 || t1 == t2 || t1 == t3)
#    sameat = (t2 == t3)
    
    sameat = (t1 == t2 && t1 == t3)

    Otot = three_body_O(dist12, dist13, dist23, sameat)

    return Otot
        
end


########
"""
    function calc_threebody(c,ind, t1,t2,t3,orb1,orb2,dist,dist31,dist32,lmn12, lmn31,lmn32; database=missing, memory0=missing, memory1=missing, memory2=missing, memoryV=missing, precalc=false, set_maxmin=false)

Calculate 3body intersite hamiltonian interactions.

`memory1`, etc have preallocated memory. This function is important for performance.
"""
function calc_threebody(c,ind, t1,t2,t3,orb1,orb2,dist,dist31,dist32,lmn12, lmn31,lmn32; database=missing, memory0=missing, memory1=missing, memory2=missing, memoryV=missing, precalc=false, set_maxmin=false)
#function calc_threebody(t1,t2,t3,orb1,orb2,dist,dist31,dist32,lmn12, lmn31,lmn32, database; set_maxmin=false)



    if precalc
        s = length(ind)
#        H = sum(memoryV[1:s] .* c.datH[ind])*10^3
#        s=size(c.datH[ind])[1]

        #H = (memoryV[1:s]'* (c.datH[ind]))[1] * 10^3
        H = ( (@view memoryV[1:s])'* (@view c.datH[ind]))[1] * 10^3        

    else
        H =  three_body_H(dist, dist31, dist32,t1==t2, t1 !=t2 && t1 != t3 && t2 != t3,  c.datH[ind], memory0=memory0, memory1=memory1, memory2=memory2, memoryV=memoryV)
    end

    sym31 = symmetry_factor(orb1,:s,lmn31, [1.0])
    sym32 = symmetry_factor(orb2,:s,lmn32, [1.0])    

    H1 = H * (sym31 * sym32   )


    return  H1


end

"""
    function fit_threebody(t1,t2,t3,orb1,orb2,dist,dist31,dist32,lmn12, lmn31,lmn32)

Fit threebody intersite interactions
"""
function fit_threebody(t1,t2,t3,orb1,orb2,dist,dist31,dist32,lmn12, lmn31,lmn32)

    #currently only for offsite Hamiltonian 3body terms, treat third atom with :s symmetry
    
    o1 = summarize_orb(orb1)
    o2 = summarize_orb(orb2)    
    
    H =  three_body_H(dist, dist31, dist32, t1==t2, t1 !=t2 && t1 != t3 && t2 != t3 )
#    println("H ", H)

    sym31 = symmetry_factor(orb1,:s,lmn31, [1.0])
    sym32 = symmetry_factor(orb2,:s,lmn32, [1.0])    

    H1 = H * (sym31 * sym32  )

    return H1


end



"""
    function renormalize_S(tbc, database, cutoff=17.99)

Function for changing S but keeping the band structure fixed. Not currently used.
"""
function renormalize_S(tbc, database, cutoff=17.99)
    #steps: 1) calculate new S_k
    #       2) change tbc to use new S_k
    #       3) trim new tbc in real space
    #       4) readjust onsite elements to match new energy


    c = tbc.crys
    tbc_new = calc_tb(c, database)

    grid= get_grid(c)

    kpts, kweights = make_kgrid(grid)

    nk = size(kpts)[1]
    nw = tbc.tb.nwan

    thetype=typeof(real(tbc.tb.H[1,1,1]))
    
    Hk_arr = zeros(Complex{thetype}, nw,nw,nk)
    Sk_arr = zeros(Complex{thetype}, nw,nw,nk)
    
    for i = 1:nk
        
        vect_old, vals_old, hk_old, sk_old = Hk(tbc.tb,     kpts[i,:])
        vect_new, vals_new, hk_new, sk_new = Hk(tbc_new.tb, kpts[i,:])

        Sk_arr[:,:,i] = (sk_new + sk_new')/2.0
        
        sk_old_invsqrt = sqrt(sk_old)^-1
        hk_orth = sk_old_invsqrt * hk_old * sk_old_invsqrt

        sk_new_sqrt = sqrt(sk_new)
        hk_new = sk_new_sqrt * hk_orth * sk_new_sqrt

        Hk_arr[:,:,i] = (hk_new + hk_new')/2.0

        if i == 1
            F_old = eigen(Hermitian(hk_old),Hermitian( sk_old))
            F_orth = eigen(Hermitian(hk_orth))
            F_new = eigen(Hermitian(hk_new), Hermitian(sk_new))
#            println("vals old ", F_old.values)
#            println("vals orth ", real(F_orth.values))
#            println("vals new ", real(F_new.values)                        )
#            println("sk_old ", diag(sk_old))
#            println("sk_new ", diag(sk_new))            
#            println("hk_new ", diag(hk_new))            
        end
    end
        
    ham_r, S_r, r_dict, ind_arr = myfft(c, true, grid, kpts,Hk_arr, Sk_arr)

    tb = make_tb(ham_r, ind_arr, r_dict, S_r )        
    nelec = tbc.nelec
    dftenergy = tbc.dftenergy

    tbc_renorm = make_tb_crys(tb, c, nelec, dftenergy)

    vect_old, vals_old, hk_old, sk_old = Hk(tbc_renorm.tb,  [0,0,0])
    F = eigen(Hermitian(hk_old), Hermitian(sk_old))
#    println("vals renorm ", real(F.values)                        )
#    println("sk_renorm ", diag(sk_old))
#    println("hk_renorm ", diag(hk_old))

    
#    trim_dist(tbc_renorm.tb, cutoff)

    tbc_renorm = trim_dist(tbc_renorm)

    vect_old, vals_old, hk_old, sk_old = Hk(tbc_renorm.tb,  [0,0,0])
    F = eigen(Hermitian(hk_old), Hermitian(sk_old))
#    println("vals trim ", real(F.values) )
#    println("sk_trim ", diag(sk_old))
#    println("hk_trim ", diag(hk_old))

    
    energy_orig = tbc.dftenergy
    energy_new = calc_energy(tbc_renorm)
    
    shift = (energy_orig - energy_new)/tbc_renorm.nelec
    c_zero = tbc_renorm.tb.r_dict[[0,0,0]]

    vects, vals, hk, sk = Hk(tbc_renorm.tb, [0.0,0.0,0.0])

    sk2 = sqrt(sk)

    toadd = sk2 * (I(nw) * shift) * sk2

    tbc_renorm.tb.H[:,:,c_zero] += toadd

    energy_after_shift = calc_energy(tbc_renorm)

    println("ENERGY SHIFTING DFT: $dftenergy , new $energy_new shift $shift final $energy_after_shift")

    vect_old, vals_old, hk_old, sk_old = Hk(tbc_renorm.tb,  [0,0,0])
    F = eigen(Hermitian(hk_old), Hermitian(sk_old))
#    println("vals shift ", real(F.values) )

    
    return tbc_renorm
    
end




function calc_tb_lowmem(crys::crystal, database=missing; reference_tbc=missing, verbose=true, var_type=missing, use_threebody=true, use_threebody_onsite=true, gamma=missing, u3=missing, background_charge_correction=0.0, screening=1.0, set_maxmin=false, check_frontier=true, check_only=false, repel = true, DIST=missing)

    #    use_threebody= false
    #    use_threebody_onsite=false
    
####    verbose = true

#    println("repel $repel -------------------------------")

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

#R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3, Rind = DIST    

#    R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3, Rind = distances_etc_3bdy_parallel(crys,cutoff2X, cutoff3bX,var_type=var_type)

    
    @time if !ismissing(DIST)
    
        R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3, Rind = DIST

    else
        if (use_threebody || use_threebody_onsite ) && !ismissing(database)
            R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3, Rind = distances_etc_3bdy_parallel(crys,cutoff2X, cutoff3bX,var_type=var_type, return_floats=false)
            DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3, Rind 

        else
            R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3, Rind = distances_etc_3bdy_parallel(crys,cutoff2X, 0.0,var_type=var_type, return_floats=false)
            DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3, Rind             
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


    
    if verbose println("check_frontier") end
    @time if !ismissing(database) && check_frontier
        #    if false

        violation_list, vio_bool, repel_vals = calc_frontier(crys, database, test_frontier=true, diststuff=DIST, verbose=verbose, var_type=var_type)
        if vio_bool == false 
            within_fit = false
        end
    else
        repel_vals = zeros(var_type, crys.nat)
    end
    
    
    
    if check_only==true
        return within_fit, sum(abs.(repel_vals)) < 1e-12
    end
    
    nwan = length(keys(ind2orb))

    nkeep=size(R_keep)[1]
    #    nkeep2=size(R_keep2)[1]    
    #    println("nkeep, $nkeep, nkeep2, $nkeep2")
    

    H = zeros(var_type, 1, nwan, nwan, nkeep)
    S = zeros(var_type, nwan, nwan, nkeep)    


    ind_arr = zeros(Int64, nkeep, 3)
    ind_arr[:,:] = R_keep[:,2:4]

    #    lmn = zeros(var_type, 3)
    #    dist = 0.0
    #    lmn31 = zeros(var_type, 3)
    #    dist31 = 0.0

    #    lmn32 = zeros(var_type, 3)
    #    dist32 = 0.0

    #    lmn41 = zeros(var_type, 3)
    #    dist14 = 0.0
    #    dist43 = 0.0


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
        end
    end

    
    warned = false
    warned_onsite = false

    nkeep_ab = size(R_keep_ab)[1]

    if !ismissing(database)


        if verbose println("2body") end
        LMN = zeros(var_type, 3, nthreads())

        @time @threads for c = 1:nkeep_ab #@threads
            id = threadid()

            #        ind_arr[c,:] = R_keep_ab[c][4:6]
            cind = R_keep_ab[c,1]
            cham = R_keep_ab[c,7]
            a1 = R_keep_ab[c,2]
            a2 = R_keep_ab[c,3]

            t1 = crys.stypes[a1]
            t2 = crys.stypes[a2]
            
            coef = database[(t1,t2)]

            indH, indS, inH, inS = coef.inds_int[[t1,t2]]

            cutoff2Xt = get_cutoff(t1,t2)[1]

            dist_a, lmn = get_dist(a1,a2, R_keep_ab[c,4:6], crys, At)

#            dist_a_old = dist_arr[a1,a2,cind,1]
            
            
            if (dist_a > cutoff2Xt || dist_a < 1e-5)
                continue
            end
            #            LMN[:,id] .= (dist_arr[a1,a2,cind,2:4])

#            LMN[:,id] .= (dist_arr[a1,a2,cind,2:4])

#            println("dist_a $dist_a $dist_a_old")

            
            lag = two_body_S(dist_a)

            if dist_a < cutoff2Xt - cutoff_length
                cut = 1.0
            else
                cut = cutoff_fn(dist_a, cutoff2Xt - cutoff_length, cutoff2Xt)
            end
            cutoff2Xa = get_cutoff(t1,t2)[1]

            ####            for o1 = orb2ind[a1]
            ####                a1a,t1,s1 = ind2orb[o1]
            ####                sum1 = summarize_orb(s1)
            ####                for o2 = orb2ind[a2]
            ####                    a2a,t2,s2 = ind2orb[o2]
            ####                    at_set = Set((t1,t2))
            ####
            ####                    cutoff2Xa = get_cutoff(t1,t2)[1]
            ####
            ####                    #                println("asdf $c $cind $a1 $a2 $o1 $o2 $s1 $s2")
            ####                    #        for o1 = 1:nwan
            ####                    #            a1,t1,s1 = ind2orb[o1]
            ####                    #            for o2 = 1:nwan
            ####                    #                a2, t2,s2 = ind2orb[o2]
            ####                    
            ####                    dist = dist_arr[a1,a2,cind,1]
            ####
            ####
            ####                    if (dist > cutoff2Xa || dist < 1e-5)
            ####                        continue
            ####                    end
            ####                    
            ####                    lmn = dist_arr[a1,a2,cind,2:4]
            ####                    #               println("jkl $t1 $t2 $s1 $s2 $dist $lmn")
            ####
            ####
            ####                    (h,s) = calc_twobody(t1,t2,s1,s2,dist,lmn, database)
            ####
            ####                    if dist < cutoff2Xa - cutoff_length
            ####                        cut = 1.0
            ####                    else
            ####                        cut = cutoff_fn(dist, cutoff2Xa - cutoff_length, cutoff2Xa)
            ####                    end
            ####
            ####                    H[o1, o2, cham] += h  *cut
            ####                    S[o1, o2, cham] += s  *cut

            
            #            for o1 = orb2ind[a1]
            #                a1a,t1,s1 = ind2orb[o1]
            #                sum1 = summarize_orb(s1)
            #                for o2 = orb2ind[a2]
            #                    a2a,t2,s2 = ind2orb[o2]
            #                    at_set = Set((t1,t2))

            
            for o1x = 1:norb[a1]
                o1 = orbs[a1,o1x]
                s1 = sorbs[a1,o1x]
                sum1 = sumorbs[a1,o1x]
                
                for o2x = 1:norb[a2]
                    o2 = orbs[a2,o2x]
                    s2 = sorbs[a2,o2x]
                    sum2 = sumorbs[a2,o2x] 

                    (hw,sw) = calc_twobody_faster(t1,t2,s1,s2,sum1, sum2, dist_a,lmn, coef, indH[sum1,sum2,:], indS[sum1,sum2,:],lag)

                    H[1, o1, o2, cham] += hw  *cut
                    S[o1, o2, cham] += sw  *cut


                end
            end
        end

            #############
        #threebody

        #        memory0=zeros(var_type, 3)
        #        memory1=zeros(var_type, 3)
        #        memory2=zeros(var_type, 3)
        #        memoryV=zeros(var_type, n_3body)

        #        lmn = zeros(var_type, 3)

        H_thread = zeros(var_type,  nwan, nwan,  nkeep,  nthreads() )


#        memory0_th=zeros(var_type, maximum([n_3body, n_3body_onsite, n_3body_onsite_same, n_3body_same]) , nthreads())
#        memory1_th=zeros(var_type, maximum([n_3body, n_3body_onsite, n_3body_onsite_same, n_3body_same]) , nthreads())
#        memory2_th=zeros(var_type, maximum([n_3body, n_3body_onsite, n_3body_onsite_same, n_3body_same]) , nthreads())
        memoryV_th=zeros(var_type, maximum([n_3body, n_3body_onsite, n_3body_onsite_same, n_3body_same]) , nthreads())

        hh = zeros(var_type, 3,3, nthreads())
#        Htemp = zeros(var_type, 16,16, nthreads())
        
        #        v =Array{Symbol}(undef, 6, nthreads())
        

        #        sdict = Dict()
        #        sumdict = Dict()
        #        for a1 = 1:crys.nat
        #            sdict[a1] = Symbol[]
        #            sumdict[a1] = Symbol[]
        #            
        #            for o1 = orb2ind[a1]
        #                a1a,t1,s1 = ind2orb[o1]
        #                sum1 = summarize_orb(s1)
        #                push!(sdict[a1], s1)
        #                push!(sumdict[a1], sum1)
        #            end
        #            
        #        end
        




        if verbose println("3body") end
#        println("3bdy")
        @time if use_threebody || use_threebody_onsite
            #        if false
            @threads for counter = 1:size(array_ind3)[1] #@threads
                #            for counter = 1:size(array_ind3)[1]
                id = threadid()
                #id = 1
                a1 = array_ind3[counter,1]
                a2 = array_ind3[counter,2]
                a3 = array_ind3[counter,3]

#                println("3bdy $a1 $a2 $a3")
                
                cind1 = array_ind3[counter,4]
                cind2 = array_ind3[counter,5]

#                println("$a1 $a2 $a3 $cind1 $cind2")
#                r_dict
                rind1 = ind_arr[cind1,1:3]
                rind2 = ind_arr[cind2,1:3]
                #                rind2 = Rind[cind2,1:3]

#                println("rind1 $rind1  rind2 $rind2")
                

                dist, lmn = get_dist(a1,a2, rind1, crys, At)
                dist31, lmn31 = get_dist(a1,a3, rind2, crys, At)
                dist32, lmn32 = get_dist(a2,a3, -rind1+rind2, crys, At)

                
#                distO = array_floats3[counter, 1]
#                dist31O = array_floats3[counter, 2]
#                dist32O = array_floats3[counter, 3]
                
                #                lmn[:] = array_floats3[counter, 4:6]
                #                lmn31[:] = array_floats3[counter, 7:9]
                #                lmn32[:] = array_floats3[counter, 10:12]

#                lmn = @view array_floats3[counter, 4:6]
#                lmn31 = @view array_floats3[counter, 4:6]
#                lmn32 = @view array_floats3[counter, 7:9]

#                println("dist3 $dist $distO  $dist31 $dist31O   $dist32 $dist32O")
#                println("lmn $lmn31_a_new $lmn31     $lmn32_a_new  $lmn32")

                
#                memory0= @view memory0_th[:,id]
#                memory1= @view memory1_th[:,id]
#                memory2= @view memory2_th[:,id]
#                memoryV= @view memoryV_th[:,id]


                t1 = crys.stypes[a1]
                t2 = crys.stypes[a2]
                t3 = crys.stypes[a3]

                cutoff3 = get_cutoff(t1,t2,t3)
                cutoffZZ = get_cutoff(t1,t2)[1]

                cut_ab = cutoff_fn(dist, cutoffZZ - cutoff_length, cutoffZZ)
                cut_ab2 = cutoff_fn(dist, cutoff3 - cutoff_length, cutoff3)
                cut_ac = cutoff_fn(dist31, cutoff3 - cutoff_length, cutoff3)
                cut_bc = cutoff_fn(dist32, cutoff3 - cutoff_length, cutoff3)

                cut = cut_ab*cut_ac*cut_bc
                cut2 = cut_ab2*cut_ac*cut_bc   
                #cut = array_floats3[counter, 10]

                #println("cut $cut $(cut_ab*cut_ac*cut_bc)")
                #                v = [t1,:s,t2,:s,t3,:H]

                
                d1 = laguerre(dist, missing, nmax=1)
                d2 = laguerre(dist31, missing, nmax=1)
                d3 = laguerre(dist32, missing, nmax=1)
                
                
                if haskey(database, (t1, t2, t3))
                    cdat = database[(t1,t2,t3)]
                    (cindX, nindX) = cdat.inds_int[[t1,t2,t3]]
                    if use_threebody  

                        
                        
                        #                        three_body_H(dist, dist31, dist32,t1==t2, memory0=memory0, memory1=memory1, memory2=memory2, memoryV=memoryV)
                        #        memoryV = three_body_H(dist, dist31, dist32,t1==t2)
                        memoryV = three_body_H_lag(d1,d2,d3,t1==t2, t1 !=t2 && t1 != t3 && t2 != t3 )

                        #puts what we need in memoryV

                        
                        for sum1 = 1:maximum(sumorbs[a1,:])
                            for sum2 = 1:maximum(sumorbs[a2,:])
                                @inbounds hh[sum1,sum2,id] = ( (@view memoryV[1:nindX[sum1, sum2]])'* (@view cdat.datH[ (@view cindX[sum1, sum2, 1:nindX[sum1, sum2]])   ]))[1]
                            end
                        end

                            
                        sym31 = 1.0
                        sym32 = 1.0                        
                        
#                        Htemp = zeros(norb[a1], norb[a2])
#                        Htemp[:,:, id] .= 0.0
                        

                        for o1x = 1:norb[a1]
                            o1 = orbs[a1,o1x]
                            s1 = sorbs[a1,o1x]
                            sum1 = sumorbs[a1,o1x]
                            
                            sym31 = symmetry_factor_int(s1,1,lmn31,one ) 

                            for o2x = 1:norb[a2]
                                o2 = orbs[a2,o2x]
                                s2 = sorbs[a2,o2x]
                                sum2 = sumorbs[a2,o2x] 

#                                hht = ( (@view memoryV[1:nindX[sum1, sum2]])'* (@view cdat.datH[ (@view cindX[sum1, sum2, 1:nindX[sum1, sum2]])   ]))[1]
                                
                                sym32 = symmetry_factor_int(s2,1,lmn32, one)    

#                                @inbounds Htemp[o1x,o2x,id] += hh[sum1,sum2, id] * sym31 * sym32

                                @inbounds H_thread[orbs[a1,o1x] , orbs[a2,o2x]  , cind1, id ] += (hh[sum1,sum2, id]) * (sym31 * sym32  * cut * 10^3)
                                
                                
                                
                            end
                        end

#                        @inbounds Htemp[1:norb[a1],1:norb[a2], id] .*= cut * 10^3
                        
#                                @inbounds H_thread[orbs[a1,1:norb[a1]] , orbs[a2,1:norb[a2]]  , cind1, id ] .+= (@view Htemp[1:norb[a1],1:norb[a2], id]) * cut * 10^3


                        
                    end
                    ############################################
                    if use_threebody_onsite 

#                        cut2 = array_floats3[counter, 11]

                        
                        for o1 = orb2ind[a1]
                            a1a,t1,s1 = ind2orb[o1]
                            o = calc_threebody_onsite_lag(t1,t2,t3,s1,d1,d2,d3, cdat, set_maxmin=set_maxmin, memory=memoryV) * cut2

                            
                            @inbounds H_thread[ o1, o1,c_zero, id] += o  ####* cut2
                        end
                    elseif !warned_onsite && use_threebody_onsite
                        println("WARNING, missing 3bdy onsite ", (t1, t2, t3))
                        warned_onsite = true
                        within_fit = false
                    end
                    ###########################################

                elseif !warned
                    println("WARNING, missing 3bdy ", (t1, t2, t3))
                    within_fit = false

                    warned = true
                end
            end
        end
        #        println("thread sum")
        #        println(size(H))
        #        println(size(H_thread))

        
        H[1,:,:,:] .+= sum(H_thread, dims=4)[:,:,:]
        #        H += sum(H_thread, dims=4)[:, :,:]



#        lmn = zeros(var_type, 3)

        ############ONSITE

        Hon = zeros(var_type, nwan,nwan, nthreads())
        Son = zeros(var_type, nwan,nwan, nthreads())
        
        if verbose println("onsite") end

        @time @threads for c = 1:nkeep_ab #@threads 
            id = threadid()

            #        ind_arr[c,:] = R_keep_ab[c][4:6]
            cind = R_keep_ab[c,1]
            a1 = R_keep_ab[c,2]
            a2 = R_keep_ab[c,3]
            t1 = crys.stypes[a1]
            t2 = crys.stypes[a2]

            dist, lmn = get_dist(a1,a2, R_keep_ab[c,4:6], crys, At)

#            dist = dist_arr[a1,a2,cind,1]
#            LMN[:, id] = dist_arr[a1,a2,cind,2:4]

            cutoff_onXa = get_cutoff(t1,t2)[2]

            if (dist > cutoff_onXa)
                continue
            end

            if dist < cutoff_onXa - cutoff_length
                cut = 1.0
            else
                cut = cutoff_fn(dist, cutoff_onXa - cutoff_length, cutoff_onXa)
            end

            
            for o1 = orb2ind[a1]
                a1a,t1a,s1 = ind2orb[o1]
                sum1 = summarize_orb(s1)
                for o2 = orb2ind[a1]
                    a2a,t2a,s2 = ind2orb[o2]
                    
                    if dist < 1e-5 #true onsite
                        (h,s) = calc_onsite(t1,s1,s2, database)
                        #                        S[o1, o2, c_zero] += s 
                        #                        H[o1, o2, c_zero] += h
                        Son[o1, o2, id] += s 
                        Hon[o1, o2, id] += h
                        if repel
                            if o1 == o2
                                Hon[o1, o1, id] += repel_vals[a1a] * 0.1
                            end
                        end
                        
                    else
                        o = calc_twobody_onsite(t1,t2, s1,s2,dist,lmn, database)
                        Hon[o1, o2, id] += o * cut
                        #                        H[o1, o2, c_zero] += o  * cut
                    end
                    



                end

            end
        end
        end
        H[1, :,:,c_zero] += sum(Hon, dims=3)[:,:]
        S[:,:,c_zero] += sum(Son, dims=3)[:,:]


    

    if verbose println("make") end
    if true
#        println("typeof H ", typeof(H), " " , size(H), " S ", typeof(S), " " , size(S))
        tb = make_tb(H, ind_arr, S)
        if !ismissing(database) && (haskey(database, "scf") || haskey(database, "SCF"))
            scf = database["scf"]
        else
            scf = false
        end
        tbc = make_tb_crys(tb, crys, nval, 0.0, scf=scf, gamma=gamma, u3=u3, background_charge_correction=background_charge_correction, within_fit=within_fit, screening=screening)
    end
    if verbose 
        println("-----")
        println()
    end

    return tbc

end















function calc_tb_lowmem2(crys::crystal, database=missing; reference_tbc=missing, verbose=true, var_type=missing, use_threebody=true, use_threebody_onsite=true, gamma=missing, u3=missing, background_charge_correction=0.0,  screening=1.0, set_maxmin=false, check_frontier=true, check_only=false, repel = true, DIST=missing)

    println("low mem 2")
    
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

#R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3, Rind = DIST    

#    R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3, Rind = distances_etc_3bdy_parallel(crys,cutoff2X, cutoff3bX,var_type=var_type)


    @time if !ismissing(DIST)
        
        R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3 = DIST

    else
        if (use_threebody || use_threebody_onsite ) && !ismissing(database)
            R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel(crys,cutoff2X, cutoff3bX,var_type=var_type, return_floats=false)
            DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3 

        else
            R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel(crys,cutoff2X, 0.0,var_type=var_type, return_floats=false)
            DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3             
        end
    end

    #println("stuff")
    begin
        
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
            
            c_zero_ref=1
            if !(ismissing(reference_tbc))
                if size(reference_tbc.tb.ind_arr)[1] > 1
                    c_zero_ref = reference_tbc.tb.r_dict[[0,0,0]]
                end
            end
        end


        
        if verbose println("check_frontier") end
        if !ismissing(database) && check_frontier
            #    if false
            #println("check ---------------------------------")
            violation_list, vio_bool, repel_vals = calc_frontier(crys, database, test_frontier=true, diststuff=DIST, verbose=verbose, var_type=var_type)

            if vio_bool == false 
                within_fit = false
            end
        else
            repel_vals = zeros(var_type, crys.nat)
        end
        
        #    println("repel_vals ", repel_vals)
        
        if check_only==true
            return within_fit, sum(abs.(repel_vals)) < 1e-12
        end
        
        nwan = length(keys(ind2orb))

        nkeep=size(R_keep)[1]
        #    nkeep2=size(R_keep2)[1]    
        #    println("nkeep, $nkeep, nkeep2, $nkeep2")
        

        H = zeros(var_type, 1, nwan, nwan, nkeep)
        S = zeros(var_type, nwan, nwan, nkeep)    


        ind_arr = zeros(Int64, nkeep, 3)
        ind_arr[:,:] = R_keep[:,2:4]

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
            end
        end

        
        warned = false
        warned_onsite = false

        nkeep_ab = size(R_keep_ab)[1]
    end
        
    if !ismissing(database)

        twobdy =  begin
            if verbose println("2body") end
            LMN = zeros(var_type, 3, nthreads())

            println("nkeep_ab $nkeep_ab")
            
            @time @threads for c = 1:nkeep_ab #add back threads
                id = threadid()

                #        ind_arr[c,:] = R_keep_ab[c][4:6]
                cind = R_keep_ab[c,1]
                cham = R_keep_ab[c,7]
                a1 = R_keep_ab[c,2]
                a2 = R_keep_ab[c,3]


                if a2 > a1
                    continue
                end

                rind1 = ind_arr[cham,1:3]
                cham_reverse = r_dict[-rind1]

                
                t1 = crys.stypes[a1]
                t2 = crys.stypes[a2]
                
                coef = database[(t1,t2)]

                indH, indS, inH, inS = coef.inds_int[[t1,t2]]

                cutoff2Xt = get_cutoff(t1,t2)[1]

                dist_a, lmn = get_dist(a1,a2, R_keep_ab[c,4:6], crys, At)

                #            dist_a_old = dist_arr[a1,a2,cind,1]
                
                
                if (dist_a > cutoff2Xt || dist_a < 1e-5)
                    continue
                end
                #            LMN[:,id] .= (dist_arr[a1,a2,cind,2:4])

                #            LMN[:,id] .= (dist_arr[a1,a2,cind,2:4])

                #            println("dist_a $dist_a $dist_a_old")

                
                lag = two_body_S(dist_a)

                if dist_a < cutoff2Xt - cutoff_length
                    cut = 1.0
                else
                    cut = cutoff_fn(dist_a, cutoff2Xt - cutoff_length, cutoff2Xt)
                end
                cutoff2Xa = get_cutoff(t1,t2)[1]


                
                for o1x = 1:norb[a1]
                    o1 = orbs[a1,o1x]
                    s1 = sorbs[a1,o1x]
                    sum1 = sumorbs[a1,o1x]
                    
                    for o2x = 1:norb[a2]
                        o2 = orbs[a2,o2x]
                        s2 = sorbs[a2,o2x]
                        sum2 = sumorbs[a2,o2x] 

                        (hw,sw) = calc_twobody_faster(t1,t2,s1,s2,sum1, sum2, dist_a,lmn, coef, indH[sum1,sum2,:], indS[sum1,sum2,:],lag)

                        
                        H[1, o1, o2, cham] += hw  *cut
                        S[o1, o2, cham] += sw  *cut

                        if (a1 != a2) 
                            
                            H[1, o2, o1, cham_reverse] += hw  *cut
                            S[o2, o1, cham_reverse] += sw  *cut
                        end

                    end
                end
            end
        end

        #############
        #threebody

        #H_thread = zeros(var_type,  nwan, nwan,  nkeep )
        H_thread = zeros(var_type,   nwan, nwan,  nkeep,  nthreads() )
        
        threebdy = begin
        


            memoryV_th=zeros(var_type, maximum([n_3body, n_3body_onsite, n_3body_onsite_same, n_3body_same]) , nthreads())

            hh = zeros(var_type, nthreads(), 3,3)

            skip = 0
            keep = 0

#            d1 = zeros(var_type, 2)
#            d2 = zeros(var_type, 2)
#            d3 = zeros(var_type, 2)
            
            
            Htemp = zeros(var_type,  9,9,nthreads())
            
            if verbose println("3body") end
            #        println("3bdy")
            @time if use_threebody || use_threebody_onsite
                #        if false
                
                @threads for  counter = 1: (size(array_ind3)[1] ) 
                    #            for counter = 1:size(array_ind3)[1]
                    a1 = array_ind3[counter,1]
                    a2 = array_ind3[counter,2]

#                    if a1 > a2
#                        continue
#                    end
                    
                    id = threadid()
                    #id = 1
                    
                    a3 = array_ind3[counter,3]

                    cind1 = array_ind3[counter,4]
                    cind2 = array_ind3[counter,5]

                    
                    rind1 = ind_arr[cind1,1:3]
                    cind1_reverse = r_dict[-rind1]

                    
                    
                    #                println("$a1 $a2 $a3 $cind1 $cind2")
                    #                r_dict
                    rind1 = ind_arr[cind1,1:3]
                    rind2 = ind_arr[cind2,1:3]
                    #                rind2 = Rind[cind2,1:3]

                    #                println("rind1 $rind1  rind2 $rind2")
                    

                    dist, lmn = get_dist(a1,a2, rind1, crys, At)
                    dist31, lmn31 = get_dist(a1,a3, rind2, crys, At)
                    dist32, lmn32 = get_dist(a2,a3, -rind1+rind2, crys, At)

                    


                    t1 = crys.stypes[a1]
                    t2 = crys.stypes[a2]
                    t3 = crys.stypes[a3]

                    cutoff3 = get_cutoff(t1,t2,t3)
                    cutoffZZ = get_cutoff(t1,t2)[1]

                    cut_ab = cutoff_fn(dist, cutoffZZ - cutoff_length, cutoffZZ)
                    cut_ab2 = cutoff_fn(dist, cutoff3 - cutoff_length, cutoff3)
                    cut_ac = cutoff_fn(dist31, cutoff3 - cutoff_length, cutoff3)
                    cut_bc = cutoff_fn(dist32, cutoff3 - cutoff_length, cutoff3)

                    cut = cut_ab*cut_ac*cut_bc
                    cut2 = cut_ab2*cut_ac*cut_bc   
                    #cut = array_floats3[counter, 10]

                    #println("cut $cut $(cut_ab*cut_ac*cut_bc)")
                    #                v = [t1,:s,t2,:s,t3,:H]

#                    laguerre(dist, missing, nmax=1, memory=d1)
#                    laguerre(dist31, missing, nmax=1, memory = d2)
#                    laguerre(dist32, missing, nmax=1, memory = d3)

                    d1 = laguerre(dist, missing, nmax=1)
                    d2 = laguerre(dist31, missing, nmax=1)
                    d3 = laguerre(dist32, missing, nmax=1)
                    
                    
                    if haskey(database, (t1, t2, t3))
                        cdat = database[(t1,t2,t3)]
                        (cindX, nindX) = cdat.inds_int[[t1,t2,t3]]
                        memoryV = three_body_H_lag(d1,d2,d3,t1==t2, t1 !=t2 && t1 != t3 && t2 != t3 )

                        if use_threebody && a1 <=  a2

                            #memoryV = three_body_H_lag(d1,d2,d3,t1==t2)
                            

                            
                            #                        three_body_H(dist, dist31, dist32,t1==t2, memory0=memory0, memory1=memory1, memory2=memory2, memoryV=memoryV)
                            #        memoryV = three_body_H(dist, dist31, dist32,t1==t2)

                            #puts what we need in memoryV

                            
                            for sum1 = 1:maximum(sumorbs[a1,:])
                                for sum2 = 1:maximum(sumorbs[a2,:])
                                    @inbounds hh[id, sum1,sum2] = ( (@view memoryV[1:nindX[sum1, sum2]])'* (@view cdat.datH[ (@view cindX[sum1, sum2, 1:nindX[sum1, sum2]])   ]))[1]
                                end
                            end

                            sym31 = 1.0
                            sym32 = 1.0                        
                            
                            
                            
                            #                        Htemp = zeros(norb[a1], norb[a2])
                            #                        Htemp[:,:, id] .= 0.0
                            
                            
                            for o1x = 1:norb[a1]
#                                o1 = orbs[a1,o1x]
                                s1 = sorbs[a1,o1x]
                                sum1 = sumorbs[a1,o1x]
                                
                                sym31 = symmetry_factor_int(s1,1,lmn31,one ) 
                                
                                for o2x = 1:norb[a2]
#                                    o2 = orbs[a2,o2x]
                                    s2 = sorbs[a2,o2x]
                                    sum2 = sumorbs[a2,o2x] 
                                    
                                    #                                hht = ( (@view memoryV[1:nindX[sum1, sum2]])'* (@view cdat.datH[ (@view cindX[sum1, sum2, 1:nindX[sum1, sum2]])   ]))[1]
                                    
                                    sym32 = symmetry_factor_int(s2,1,lmn32, one)    
                                    
                                    #@inbounds Htemp[o1x,o2x,id] += hh[sum1,sum2, id] * sym31 * sym32
                                    
                                    #vtemp = (hh[id, sum1,sum2]) * (sym31 * sym32  * cut * 10^3)

                                    Htemp[ o1x, o2x, id] = (hh[id, sum1,sum2]) * (sym31 * sym32  * cut * 10^3)

                                    #                                    @inbounds H_thread[id, orbs[a1,o1x] , orbs[a2,o2x]  , cind1 ] += vtemp
#                                    if (a1 != a2) 
#                                        @inbounds H_thread[id, orbs[a2,o2x], orbs[a1,o1x] , cind1_reverse] += vtemp##
#                                    end
                                    
                                    
                                end
                            end

                            #@inbounds Htemp[1:norb[a1],1:norb[a2], id] .*= cut * 10^3
                            
                            #                                @inbounds H_thread[orbs[a1,1:norb[a1]] , orbs[a2,1:norb[a2]]  , cind1, id ] .+= (@view Htemp[1:norb[a1],1:norb[a2], id]) * cut * 10^3




                            @inbounds     H_thread[ orbs[a1,1]:orbs[a1,norb[a1]], orbs[a2,1]:orbs[a2,norb[a2]], cind1, id]              .+= (@view Htemp[1:norb[a1],1:norb[a2], id])
                            if (a1 != a2) 
                             @inbounds H_thread[ orbs[a2,1]:orbs[a2,norb[a2]], orbs[a1,1]:orbs[a1,norb[a1]] , cind1_reverse, id] .+= (@view Htemp[1:norb[a1],1:norb[a2], id])'
                            end
                            
                        end
                        
                        ############################################
                        if use_threebody_onsite 

                            #                        cut2 = array_floats3[counter, 11]


                            for o1 = orb2ind[a1]
                                a1a,t1,s1 = ind2orb[o1]
                                o = calc_threebody_onsite_lag(t1,t2,t3,s1,d1,d2,d3, cdat, set_maxmin=set_maxmin) * cut2 
#                                println("$a1 $o1 o $(o*cut2) $cut2")
                                
                                @inbounds H_thread[  o1, o1,c_zero, id] += o  ####* cut2
                            end
                        elseif !warned_onsite && use_threebody_onsite
                            println("WARNING, missing 3bdy onsite ", (t1, t2, t3))
                            warned_onsite = true
                            within_fit = false
                        end
                        ###########################################

                    elseif !warned
                        println("WARNING, missing 3bdy ", (t1, t2, t3))
                        within_fit = false

                        warned = true
                    end
                end
            end


        end
            
        #        println("thread sum")
        #        println(size(H))
        #        println(size(H_thread))

        
#        H[1,:,:,:] .+= sum(H_thread, dims=4)[:,:,:]


        
        #        H += sum(H_thread, dims=4)[:, :,:]



#        lmn = zeros(var_type, 3)

        ############ONSITE

        
        Hon = zeros(var_type, nwan,nwan, nthreads())
        Son = zeros(var_type, nwan,nwan, nthreads())
        
        if verbose println("onsite") end

        @threads for c = 1:nkeep_ab #@threads 
            id = threadid()

            #        ind_arr[c,:] = R_keep_ab[c][4:6]
            cind = R_keep_ab[c,1]
            a1 = R_keep_ab[c,2]
            a2 = R_keep_ab[c,3]
            t1 = crys.stypes[a1]
            t2 = crys.stypes[a2]

            dist, lmn = get_dist(a1,a2, R_keep_ab[c,4:6], crys, At)

#            dist = dist_arr[a1,a2,cind,1]
#            LMN[:, id] = dist_arr[a1,a2,cind,2:4]

            cutoff_onXa = get_cutoff(t1,t2)[2]

            if (dist > cutoff_onXa)
                continue
            end

            if dist < cutoff_onXa - cutoff_length
                cut = 1.0
            else
                cut = cutoff_fn(dist, cutoff_onXa - cutoff_length, cutoff_onXa)
            end

            
            for o1 = orb2ind[a1]
                a1a,t1a,s1 = ind2orb[o1]
                sum1 = summarize_orb(s1)
                for o2 = orb2ind[a1]
                    a2a,t2a,s2 = ind2orb[o2]
                    
                    if dist < 1e-5 #true onsite
                        (h,s) = calc_onsite(t1,s1,s2, database)
                        #                        S[o1, o2, c_zero] += s 
                        #                        H[o1, o2, c_zero] += h
                        Son[o1, o2, id] += s 
                        Hon[o1, o2, id] += h
                        if repel
                            if o1 == o2
                                Hon[o1, o1, id] += repel_vals[a1a] * 0.1
                            end
                        end
                        
                    else
                        o = calc_twobody_onsite(t1,t2, s1,s2,dist,lmn, database)
                        Hon[o1, o2, id] += o * cut
                        #                        H[o1, o2, c_zero] += o  * cut
                    end
                    



                end

            end
        end
        end
#        H[1, :,:,c_zero] += sum(Hon, dims=3)[:,:]
#        S[:,:,c_zero] += sum(Son, dims=3)[:,:]

        
#    wait(twobdy)
#    wait(twobdy_on)
#    wait(threebdy)

    H[1,:,:,:] .+= sum(H_thread, dims=4)[ :,:,:]
    H[1, :,:,c_zero] += sum(Hon, dims=3)[:,:]
    S[:,:,c_zero] += sum(Son, dims=3)[:,:]
        


    if verbose println("make") end
    if true
#        println("typeof H ", typeof(H), " " , size(H), " S ", typeof(S), " " , size(S))
        tb = make_tb(H, ind_arr, S)
        if !ismissing(database) && (haskey(database, "scf") || haskey(database, "SCF"))
            scf = database["scf"]
        else
            scf = false
        end
        tbc = make_tb_crys(tb, crys, nval, 0.0, scf=scf, gamma=gamma, u3=u3,background_charge_correction=background_charge_correction, within_fit=within_fit, screening=screening)
    end
    if verbose 
        println("-----")
        println()
    end

    return tbc

end





#-----------------------------------------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

function calc_tb_LV(crys::crystal, database=missing; reference_tbc=missing, verbose=true, var_type=missing, use_threebody=true, use_threebody_onsite=true,use_eam=false, gamma=missing,u3=missing, background_charge_correction=0.0,  screening=1.0, set_maxmin=false, check_frontier=true, check_only=false, repel = true, DIST=missing, tot_charge=0.0, retmat=false, Hin=missing, Sin=missing, atom = -1)


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

    if verbose println("check_frontier") end
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

    #println("stuff")
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
            if ismissing(Hin)
                H = zeros(var_type, nwan, nwan, nkeep)
                S = zeros(var_type, nwan, nwan, nkeep)
            else
                H = Hin
                S = Sin
                H .= 0.0
                S .= 0.0
            end
        end
        #println("done memory")
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
                
                orbs_arr = zeros(Int64, crys.nat, 9, 3)
                DAT_IND_ARR = zeros(Int64, types_counter, types_counter, 2,4,4, 33 )
                DAT_ARR = zeros(Float64, types_counter, types_counter, 2, 164)

                DAT_IND_ARR_O = zeros(Int64, types_counter, types_counter, 4,4, 33 )


                cutoff_arr = zeros(crys.nat, crys.nat, 2)
                cutoff_arr3 = zeros(crys.nat, crys.nat, crys.nat)
                get_cutoff_pre = Dict()       
                s = Set(crys.stypes)
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
                        for c = 1:crys.nat
                            tc = crys.stypes[c]            
                            cutoff_arr3[a,b,c] =  get_cutoff_pre[(ta,tb,tc)]
                        end
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
                    
                    DAT_IND_ARR_3 = zeros(Int64, types_counter, types_counter, types_counter, 4,4, 33 )
                    DAT_ARR_3 = zeros(Float64, types_counter, types_counter, types_counter, 224)

                    
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
                
                lag_arr_TH = zeros(var_type, maximum([n_2body, n_2body_onsite, n_2body_S]), nthreads())
                lmn_arr_TH = zeros(var_type, 3, nthreads())
                sym_arr_TH = zeros(var_type, 3, nthreads())
                sym_arrS_TH = zeros(var_type, 3, nthreads())

            end
            #; println("end")
            

        end
    end
    
    if verbose println("2body LV") end


    twobody_LV = begin
        
        rho_th = zeros(var_type, crys.nat, 3, nthreads())
        #println("nkeep_ab $nkeep_ab")
        begin
            #@time twobody(nkeep_ab)
            for c = 1:nkeep_ab
                begin
                    id = threadid()
                    lmn_arr = lmn_arr_TH[:,id]
                    sym_arr = sym_arr_TH[:,id]
                    sym_arrS = sym_arrS_TH[:,id]
                    lag_arr = lag_arr_TH[:,id]
                    
                    cind = R_keep_ab[c,1]
                    cham = R_keep_ab[c,7]
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
                            S[o1, o1, c_zero] += s 
                            H[o1, o1, c_zero] += h

                            
                            if repel
                                H[o1, o1, c_zero] += repel_vals[a1] * 0.1
                                #println("add repel $a1 $o1 ", repel_vals[a1] * 0.1)
                            end
                        end
                        
                        
                    else #normal case
                        

                        laguerre_fast!(dist_a, lag_arr)

                        rho_th[a1, 1, id] += lag_arr[1] * cut_on
                        rho_th[a1, 2, id] += lag_arr[2] * cut_on

#                        println("add rho $a1 dist $dist_a 1 $(lag_arr[1] * cut_on) 2 $(lag_arr[2] * cut_on)")
                        
                        core!(cham, a1, a2, t1, t2, norb, orbs_arr, DAT_IND_ARR, lag_arr, DAT_ARR, cut_a, H, S, lmn_arr, sym_arr, sym_arrS)
                        core_onsite!(c_zero, a1, a2, t1, t2, norb, orbs_arr, DAT_IND_ARR_O, lag_arr, DAT_ARR, cut_on, H, lmn_arr, sym_arr, sym_arrS)
                        
                    end
                    
                end
                
            end
        end
    end
    
    if use_eam
        rho = sum(rho_th, dims=3)
#        println("rho $rho")
        for a = 1:crys.nat
            t = crys.stypes[a]
            d = database[(:eam, t)].datH
            temp = d[1]*rho[a,1]^2 +  d[2]*rho[a,2]^2 + d[3]*rho[a,1]*rho[a,2]
#            println("a $a d $d rho[a,:] $(rho[a,:])  temp $temp     vals $([rho[a,1]^2, rho[a,2]^2, rho[a,1]*rho[a,2]])  ")
            for ox = 1:norb[a]
                oxx = orbs_arr[a,ox,1]
                H[ oxx, oxx, c_zero] += temp
            end
        end
    end
        #o = orbs_arr[a,ox,1]
            #s = orbs_arr[a,ox,2]
            #sum1 = orbs_arr[a,ox,3]
            
        
        #                dist_a, lmn = get_dist(a1,a2, R_keep_ab[c,4:6], crys, At)

        #            dist_a_old = dist_arr[a1,a2,cind,1]
        
        
        #                if (dist_a > cutoff2Xt || dist_a < 1e-5)
        #                    continue
        #                end
        
        
        
        


    threebdy_LV = begin
        
        
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
            
            if verbose println("3body LV") end
            
            #                for mc in meta_count #@threads 
            #                    for counter in mc
            for  counter = 1: (size(array_ind3)[1] ) #add threads back
                
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

                #                    println("$counter ref $a1 $a2 $a3 , $cind1 $cind2")
                
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
                    cutoff3 = cutoff_arr3[a1,a2,a3]

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



    if retmat
        #return H, S
        return 
    end
    if verbose println("make") end
    begin
        #        println("typeof H ", typeof(H), " " , size(H), " S ", typeof(S), " " , size(S))
        #println("maketb")
        tb = make_tb( reshape(H, 1,size(H)[1], size(H)[2], size(H)[3])  , ind_arr, S)
        if !ismissing(database) && (haskey(database, "scf") || haskey(database, "SCF"))
            scf = database["scf"]
        else
            scf = false
        end
        #println("make")
        tbc = make_tb_crys(tb, crys, nval, 0.0, scf=scf, gamma=gamma, u3=u3, background_charge_correction=background_charge_correction, within_fit=within_fit, screening=screening)
        tbc.tot_charge = tot_charge
        tbc.nelec = tbc.nelec - tot_charge
    end
    if verbose 
        println("-----")
        println()
    end

    
    return tbc

end


function core!(cham, a1, a2, t1, t2, norb, orbs_arr, DAT_IND_ARR, lag_arr, DAT_ARR, cut_a, H, S, lmn, sym_dat, sym_datS)
    
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

            
            H[ o1, o2, cham] = symmetry_factor_int(s1, s2, lmn, sym_dat)  * cut_a
            S[ o1, o2, cham] = symmetry_factor_int(s1, s2, lmn, sym_datS)  * cut_a
            
        end
    end
end

function core_onsite!(c_zero, a1, a2, t1, t2, norb, orbs_arr, DAT_IND_ARR_O, lag_arr, DAT_ARR, cut_on, H,  lmn, sym_dat1, sym_dat2)

    
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

            H[ o1, o2, c_zero] += (temp1 + temp2 + temp3)  * cut_on

            #if abs((temp1 + temp2)) > 1e-5
            #    println("$c_zero, $a1, $a2, $t1, $t2, $(temp1 + temp2)   $cut_on")
            #end
            
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
#            println("add $o1 $o2 $cind1")
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


include("CalcTB_laguerre_sparse.jl")

end #end module




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
#using Statistics

using Base.Threads
#using Base.lock

using ..AtomicMod:atom
using ..CrystalMod:crystal
#using ..TB:orbital_index
using ..CrystalMod:orbital_index

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
using ..Utility:arr2str
using ..Utility:str_w_spaces
using ..Utility:parse_str_ARR_float

using ..Utility:dict2str
using ..Utility:str2tuplesdict

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



#old

const n_2body = 5
const n_2body_onsite = 4

const n_2body_S = 6

const n_3body = 4
const n_3body_same = 3

const n_3body_onsite = 4
const n_3body_onsite_same = 4



#=
const n_2body = 6

const n_2body_onsite = 5

const n_2body_S = 6

#const n_3body = 8
const n_3body = 4
const n_3body_same = 6

#const n_3body_onsite = 2
const n_3body_onsite = 4

#const n_3body_onsite_same = 4
const n_3body_onsite_same = 5
=#


const cutoff2X = 18.51 
const cutoff3bX = 13.01 
const cutoff_onX = 18.01 

const cutoff_length = 0.5


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
struct coefs

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
    dist_frontier = str2tuplesdict(eval(d["coefs"]["dist_frontier"]))

    version = 1
    if "version" in keys(d["coefs"])
        version = parse(Int64, d["coefs"]["version"])
    end
    println("version $version")
    
    co = make_coefs(names,dim, datH=datH, datS=datS, min_dist=min_dist, dist_frontier = dist_frontier, version=version)

    return co
    
end



"""
    function make_coefs(at_list, dim; datH=missing, datS=missing, cutoff=18.01, min_dist = 3.0, fillzeros=false, maxmin_val_train=missing, dist_frontier=missing)

Constructor for `coefs`. Can create coefs filled with ones for testing purposes.

See `coefs` to understand arguments.
"""
function make_coefs(at_list, dim; datH=missing, datS=missing, cutoff=18.01, min_dist = 3.0, fillzeros=false, dist_frontier=missing, version=3)

    println("make coefs")
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
            if length(k) == 4
                continue
            end

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
        end
    end


    if dim == 3
        ninds_int = zeros(UInt16, 4,4)

        for k in keys(data_info)
            if length(k) == 4
                continue
            end

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
        end
    end

    return coefs(dim, datH, datS, totH, totS, data_info, inds_int, at_list, orbs, cutoff, min_dist, dist_frontier2, version)

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
            
            for o1 in orbs1
                for o2 in orbs2
                    if same_at && ((o2 == :s && o1 == :p) || (o2 == :s && o1 == :d) || (o2 == :p && o1 == :d))
                        continue
                    end
                    
                    #                    push!(orbs, (o1, o2, symb))

                    if [at1, o1, at2, o2, at3,  symb] in keys(data_info)
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
                        data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1, 3, 2, 4 ] #switch 2 4 and 3 6

                        #data_info[[at2, o2, at1, o1, at3,  symb]] = tot .+ [1, 3, 2, 4,5,6, 8, 7 ] #switch 2 4 and 3 6

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
#                for o1 in orbs1
                    
                    if [at1, at2, at3,  symb] in keys(data_info)
                        return tot
                    end

                    data_info[[at1, at2, at3,  symb]] = collect(tot+1:tot+n)
                    data_info[[at1, at3, at2,  symb]] = collect(tot+1:tot+n)
                    tot += n                               #       1 2 3 4 5 6 7 8
#                end
            else
                data_info[[at1, at2, at3,  symb]] = collect(tot+1:tot+n)
                data_info[[at1, at3, at2,  symb]] = tot .+ [1, 3, 2, 4]

                #               for o1 in orbs1
 #                   if [at1, at2, at3, o1,  symb] in keys(data_info)
 #                       continue
 #                   end
 #                   data_info[[at1, at2, at3,o1,  symb]] = collect(tot+1:tot+n)
 #                   data_info[[at1, at3, at2,o1,  symb]] = collect(tot+1:tot+n)
     #           end
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
            if p[2] == p[3]
#            if p[2] == p[3] && p[2] == p[1]
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
    println(io, "min dist ", d.min_dist)
    println(io, "dist frontier ")
    for t in keys(d.dist_frontier)
        println(io, t, "    ", d.dist_frontier[t])
    end

    println(io)
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





function distances_etc_3bdy_parallel(crys, cutoff=missing, cutoff2=missing; var_type=Float64)
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
        
        begin
    
    R_keep_ab = R_keep_ab[1:keep_counter,:]
    
    ############

    MEMCHUNK = min(nr*nr * crys.nat^3, 2000)
            if !threebody
MEMCHUNK = 10
            end
    TOTMEM = MEMCHUNK * ones(Int32, nthreads())

    AI3 = []
    AF3 = []
    COUNTER = zeros(Int32, nthreads())
    for i = 1:nthreads()
        array_ind3 = zeros(Int32, MEMCHUNK, 5)
        array_floats3 = zeros(var_type, MEMCHUNK , 14)
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
                                array_ind3X = [array_ind3X;zeros(Int32, MEMCHUNK, 5)]
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

    return R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3
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
    
    return R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3
#    return R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3

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
function calc_tb_fast(crys::crystal, database=missing; reference_tbc=missing, verbose=true, var_type=missing, use_threebody=true, use_threebody_onsite=true, gamma=missing, screening=1.0, set_maxmin=false, check_frontier=true, check_only=false, repel = true)

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

        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel(crys,cutoff2X, cutoff3bX,var_type=var_type)
#        println("done dist")
        #        else
        #            R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy(crys,cutoff2X, cutoff3bX,var_type=var_type)
        #        end        
    else
        #        parallel = true
        #        if parallel

        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel(crys,cutoff2X, 0.0,var_type=var_type)

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

    if verbose println("check_frountier") end
    if !ismissing(database) && check_frontier
        #    if false
        diststuff = (R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3)
        violation_list, vio_bool, repel_vals = calc_frontier(crys, database, test_frontier=true, diststuff=diststuff, verbose=verbose, var_type=var_type)
        if vio_bool == false
            within_fit = false
        end
    else
        repel_vals = zeros(var_type, crys.nat)
    end
    
    if check_only==true
        return within_fit
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
            @threads for counter = 1:size(array_ind3)[1]
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

                lmn = @view array_floats3[counter, 4:6]
                lmn31 = @view array_floats3[counter, 7:9]
                lmn32 = @view array_floats3[counter, 10:12]

                memory0= @view memory0_th[:,id]
                memory1= @view memory1_th[:,id]
                memory2= @view memory2_th[:,id]
                memoryV= @view memoryV_th[:,id]


                cut = array_floats3[counter, 13]

                t1 = crys.stypes[a1]
                t2 = crys.stypes[a2]
                t3 = crys.stypes[a3]

                #                v = [t1,:s,t2,:s,t3,:H]


                
                if haskey(database, (t1, t2, t3))
                    cdat = database[(t1,t2,t3)]
                    (cindX, nindX) = cdat.inds_int[[t1,t2,t3]]
                    if use_threebody

                        
                        
                        three_body_H(dist, dist31, dist32,t1==t2, memory0=memory0, memory1=memory1, memory2=memory2, memoryV=memoryV)

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

                        cut2 = array_floats3[counter, 14]
                        o = calc_threebody_onsite(t1,t2,t3,dist,dist31,dist32, database, set_maxmin=set_maxmin, memory=memoryV)
                        for o1 = orb2ind[a1]
                            #                            a1a,t1,s1 = ind2orb[o1]
                            @inbounds H_thread[ o1, o1,c_zero, id] += o  * cut2
                        end


#                        cut2 = array_floats3[counter, 14]
#                        for o1 = orb2ind[a1]
#                            a1a,t1,s1 = ind2orb[o1]
#                            o = calc_threebody_onsite(t1,t2,t3,s1,dist,dist31,dist32, database, set_maxmin=set_maxmin, memory=memoryV)

 #                           
 #                           @inbounds H_thread[ o1, o1,c_zero, id] += o  * cut2
 #                       end

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
                            Hon[o1, o2, id] += repel_vals[a1a]
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
        tbc = make_tb_crys(tb, crys, nval, 0.0, scf=scf, gamma=gamma, within_fit=within_fit, screening=screening)
    end
    if verbose 
        println("-----")
        println()
    end

    return tbc

end


function calc_tb_fast_old(crys::crystal, database=missing; reference_tbc=missing, verbose=true, var_type=missing, use_threebody=true, use_threebody_onsite=true, gamma=missing, screening=1.0, set_maxmin=false, check_frontier=true, check_only=false)

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

        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel(crys,cutoff2X, cutoff3bX,var_type=var_type)

        #        else
#            R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy(crys,cutoff2X, cutoff3bX,var_type=var_type)
#        end        
    else
#        parallel = true
#        if parallel

        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel(crys,cutoff2X, 0.0,var_type=var_type)

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
        diststuff = (R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3)
        violation_list, vio_bool = calc_frontier(crys, database, test_frontier=true, diststuff=diststuff, verbose=verbose)
        if vio_bool == false
            within_fit = false
        end
    end
    if check_only==true
        return within_fit
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

                        
                        three_body_H(dist, dist31, dist32,t1==t2, memory0=memory0, memory1=memory1, memory2=memory2, memoryV=memoryV) #puts what we need in memoryV
                        
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
        tbc = make_tb_crys(tb, crys, nval, 0.0, scf=scf, gamma=gamma, within_fit=within_fit, screening=screening)
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
function calc_frontier(crys::crystal, frontier; var_type=Float64, test_frontier=false, diststuff=missing, verbose=true)

    if ismissing(var_type)
        var_type=Float64
    end

    lim = 0.05
    
    ind2orb, orb2ind, etotal, nval = orbital_index(crys)
    use_threebody=true

    if ismissing(diststuff)
#        println("distances")
        #        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy(crys,cutoff2X, cutoff3bX,var_type=var_type)
        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel(crys,cutoff2X, cutoff3bX,var_type=var_type)
    else
        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = diststuff
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

    
#    println("a")
    for c = 1:nkeep_ab

        cind = R_keep_ab[c,1]
        cham = R_keep_ab[c,7]
        a1 = R_keep_ab[c,2]
        a2 = R_keep_ab[c,3]
        t1= crys.stypes[a1]
        t2= crys.stypes[a2]
        dist = dist_arr[a1,a2,cind,1]

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

                if !ismissing(d) && dist < d*(1+lim) && dist > 0.1
#                    println("type of repel_vals[a1] ", typeof(repel_vals[a1]))
#                    println("type fo repel ", typeof(repel_short_dist_fn(dist, d, lim) * 0.1))

                    repel_vals[a1] += repel_short_dist_fn(dist, d, lim) * 0.1
                    repel_vals[a2] += repel_short_dist_fn(dist, d, lim) * 0.1
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
            else
                frontier[(t1,t2)] = dist
            end
        end

    end

#    println("repel 2 ", repel_vals)
    #threebody


    for counter = 1:size(array_ind3)[1]
        a1 = array_ind3[counter,1]
        a2 = array_ind3[counter,2]
        a3 = array_ind3[counter,3]

#        cind1 = array_ind3[counter,4]
#        cind2 = array_ind3[counter,5]

        dist = array_floats3[counter, 1]
        dist31 = array_floats3[counter, 2]
        dist32 = array_floats3[counter, 3]
                
        t1 = crys.stypes[a1]
        t2 = crys.stypes[a2]
        t3 = crys.stypes[a3]
        
        if t1 == t2
            tmp1 = min(dist31, dist32)
            tmp2 = max(dist31, dist32)
            
            dist31=tmp1
            dist32=tmp2
        end



#        println("dist $t1 $t2 $t3 ", [dist,dist31,dist32])


        if test_frontier                 ##############
            if dist > 9.5 || dist31 > 9.5 || dist32 > 9.5
                continue
            end

            if haskey(frontier, (t1,t2,t3))

                if !isa(frontier[(t1,t2,t3)], Array) ####&& !ismissing(frontier[(t1,t2)])
#                if typeof(frontier[(t1,t2,t3)]) == coefs
                    if ismissing(frontier[(t1,t2,t3)].dist_frontier)
                        continue
                    end
                    vals = frontier[(t1,t2,t3)].dist_frontier[(t1,t2,t3)]
#                    vals = frontier[(t1,t2,t3)].dist_frontier
                else
                    vals = frontier[(t1,t2,t3)]
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
                    if dist >= f[1]*(1+lim) && dist31 >= f[2]*(1+lim) && dist32 >= f[3]*(1+lim)
                        vio_lim = false
                    end
                end

                if vio_lim == true
#                    println("vio_lim true")
                    rsum = 10000000.0
                    rvals = zeros(var_type, 3)
                    for f in vals
                        if dist <= f[1]*(1+lim) || dist31 <= f[2]*(1+lim) || dist32 <= f[3]*(1+lim) 
#                            println("dist $dist $dist31 $dist32 " , f)

                            rvals_t = [repel_short_dist_fn(dist, f[1], lim),repel_short_dist_fn(dist31, f[2], lim),repel_short_dist_fn(dist32, f[3], lim)]
                            if sum(rvals_t) < rsum
                                rsum = sum(rvals_t)
                                rvals[:] = rvals_t[:]
                            end
                        end
                    end
                    repel_vals[a1] += rvals[1] * 0.1
                    repel_vals[a2] += rvals[1] * 0.1
                    
                    repel_vals[a1] += rvals[2] * 0.1
                    repel_vals[a3] += rvals[2] * 0.1
                    
                    repel_vals[a2] += rvals[3] * 0.1
                    repel_vals[a3] += rvals[3] * 0.1
                    
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
function calc_tb_prepare_fast(reference_tbc::tb_crys; use_threebody=false, use_threebody_onsite=false, spin=1)

#    println("calc_tb_prepare_fast 3bdy $use_threebody    3bdy-onsite $use_threebody_onsite")
#    println(reference_tbc.crys)
#    println()
    
    crys = reference_tbc.crys
    
    var_type=typeof(crys.coords[1,1])

#    if ismissing(var_type)
#        var_type=Float64
#    end
    
    ind2orb, orb2ind, etotal, nval = orbital_index(crys)

    if use_threebody || use_threebody_onsite
        #println("distances")
        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel(crys,cutoff2X, cutoff3bX)
    else
        R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel(crys,cutoff2X, 0.0)
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
                svec[ind] = real(reference_tbc.tb.S[o1,o2,c_ref])
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
            
            lmn[:] = array_floats3[counter, 4:6]
            lmn31[:] = array_floats3[counter, 7:9]
            lmn32[:] = array_floats3[counter, 10:12]

            cut = array_floats3[counter, 13]
            cut2 = array_floats3[counter, 14]

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
                            threebody_arrays[at_set3][1][ind,ih] += h[1:size(ih)[1]] * cut
                        end
                        
                    end
                end
            end
            if use_threebody_onsite  #3bdy onsite!!!!!!!!!!!
                

                h = fit_threebody_onsite(t1,t2,t3,dist,dist31,dist32)

                for o1 = orb2ind[a1]
                    a1a,t1,s1 = ind2orb[o1]
                    sum1 = summarize_orb(s1)

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
                    
                    if dist < cutoff_onX - cutoff_length
                        cut = 1.0
                    else
                        cut = cutoff_fn(dist, cutoff_onX - cutoff_length, cutoff_onX)
                    end
                    
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

    return twobody_arrays, threebody_arrays, hvec, svec, Rvec, INDvec, h_onsite, ind_conversion, dmin_types, dmin_types3

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
    if s1 == s2
        H = atoms[t1].eigs[summarize_orb(s1)]

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

"""
    function laguerre(dist, ind=missing; nmax=6, memory=missing)

Calculate laguerre polynomials up to order `nmax`
"""
function laguerre(dist, ind=missing; nmax=6, memory=missing)

    a=2.0


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
    if nmax >= 1
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
    d1 = laguerre(dist1, missing, nmax=2 )
    d2 = laguerre(dist2, missing, nmax=2)
    d3 = laguerre(dist3, missing, nmax=2)

        
    if same_atom
    
        if  isa(dist1, Array)
            
            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2] d1[:,2].*d2[:,2].*d3[:,2] ]


#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2]  ]
#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1] + d1[:,1].*d2[:,1].*d3[:,2]) ]


        else

            if ismissing(memoryV)
                V = zeros(typeof(d1[1]), 4) 
            else
                V = memoryV
            end
            
            V[1:4] .= [d1[1].*d2[1].*d3[1], (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]), d1[1].*d2[1].*d3[2], d1[2].*d2[2].*d3[2]]

#            V = [d1[1].*d2[1].*d3[1] (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]) d1[1].*d2[1].*d3[2]]

#            V = [d1[1].*d2[1].*d3[1] (d1[1].*d2[1].*d3[2]+d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1])]


        end
    else
        if  isa(dist1, Array)
#            V = [d1[:,1].*d2[:,1].*d3[:,1] (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]) d1[:,1].*d2[:,1].*d3[:,2]  ]

#            V = [d1[:,1].*d2[:,1].*d3[:,1]  (d1[:,2].*d2[:,1].*d3[:,1] + d1[:,1].*d2[:,2].*d3[:,1]+ d1[:,1].*d2[:,1].*d3[:,2]) ]
            V = [d1[:,1].*d2[:,1].*d3[:,1] d1[:,2].*d2[:,1].*d3[:,1] d1[:,1].*d2[:,2].*d3[:,1] d1[:,1].*d2[:,1].*d3[:,2]]

        else

#            V = [d1[1].*d2[1].*d3[1] (d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]) d1[1].*d2[1].*d3[2]]
#            V = d1[1].*d2[1].*d3[1] #(d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1]) d1[1].*d2[1].*d3[2]]

#            V = [d1[1].*d2[1].*d3[1] (d1[1].*d2[1].*d3[2]+d1[2].*d2[1].*d3[1] + d1[1].*d2[2].*d3[1])]

            if ismissing(memoryV)
                V = zeros(typeof(d1[1]), 4) 
            else
                V = memoryV
            end

            
            V[1:4] .= [d1[1].*d2[1].*d3[1],       d1[2].*d2[1].*d3[1] ,      d1[1].*d2[2].*d3[1] ,       d1[1].*d2[1].*d3[2]]

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
function three_body_H(dist0, dist1, dist2, same_atom, ind=missing; memory0=missing, memory1=missing, memory2=missing, memoryV=missing)

#    return 0.0

    zero =  laguerre(dist0,missing, nmax=2, memory=memory0)
    a = laguerre(dist1,missing, nmax=2, memory=memory1)
    b = laguerre(dist2,missing, nmax=2, memory=memory2)



#    zero = memory0
#    a = memory1
#    b = memory2
    
    if same_atom
        if  isa(dist1, Array)
#            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  a[:,2].*b[:,2] (a[:,1].*b[:,3]+ a[:,3].*b[:,1])  zero[:,1].*a[:,1].*b[:,1]  zero[:,1].*(a[:,1].*b[:,2]+a[:,2].*b[:,1]) ]
#            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  a[:,2].*b[:,2] (a[:,1].*b[:,3]+ a[:,3].*b[:,1]) zero[:,1].*a[:,1].*b[:,1]   ]
#            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  a[:,2].*b[:,2] (a[:,1].*b[:,3]+ a[:,3].*b[:,1]) zero[:,1].*a[:,1].*b[:,1]      ]
#            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  a[:,2].*b[:,2] (a[:,1].*b[:,3]+ a[:,3].*b[:,1]) zero[:,1].*a[:,1].*b[:,1]   zero[:,1].*(a[:,1].*b[:,2]+a[:,2].*b[:,1]) zero[:,2].*(a[:,1].*b[:,1])   ]

            Vt = [a[:,1].*b[:,1]  (a[:,1].*b[:,2]+ a[:,2].*b[:,1])  zero[:,1].*a[:,1].*b[:,1]      ]
            V = Vt
        else 
           try
#                Vt = [a[1].*b[1]  (a[1].*b[2]+a[2].*b[1])  a[2].*b[2] (a[1].*b[3]+ a[3].*b[1])  zero[1].*a[1].*b[1] zero[1].*(a[1].*b[2]+a[2].*b[1]) ]
#                V = Vt
                if ismissing(memoryV)
                    memoryV=zeros(typeof(dist0), n_3body)
                end
                memoryV[1] = a[1].*b[1]
                memoryV[2] =  (a[1].*b[2]+a[2].*b[1])
                memoryV[3] =  zero[1].*a[1].*b[1]

#                memoryV[3] =   a[2].*b[2] 
#                memoryV[4] =  (a[1].*b[3]+ a[3].*b[1])
#                memoryV[6] = zero[1].*(a[1].*b[2]+a[2].*b[1])
#                memoryV[7] = zero[2].*(a[1].*b[1])


            catch
                println("asdf ",size(a), " " , size(b))
            end
        end
    else
        if  isa(dist1, Array)
#            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]    a[:,2].*b[:,1] a[:,2].*b[:,2]  a[:,1].*b[:,3]    a[:,3].*b[:,1]   zero[:,1].*a[:,1].*b[:,1]   zero[:,1].*a[:,1].*b[:,2]  zero[:,1].*a[:,2].*b[:,1]]

#            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]    a[:,2].*b[:,1] a[:,2].*b[:,2]   zero[:,1].*a[:,1].*b[:,1] a[:,3].*b[:,1] a[:,1].*b[:,3]]
            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]    a[:,2].*b[:,1]    zero[:,1].*a[:,1].*b[:,1] ]

#            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]    a[:,2].*b[:,1] a[:,2].*b[:,2]   zero[:,1].*a[:,1].*b[:,1] ]
            V = Vt
        else
            try
#                Vt = [a[1].*b[1]  a[1].*b[2]   a[2].*b[1]   a[2].*b[2]  a[1].*b[3]    a[3].*b[1]   zero[1].*a[1].*b[1]  zero[1].*a[1].*b[2]  zero[1].*a[2].*b[1] ]
#                V = Vt
                if ismissing(memoryV)
                    memoryV=zeros(typeof(dist0), n_3body)
                end
                memoryV[1] =  a[1].*b[1]
                memoryV[2] =  a[1].*b[2]
                memoryV[3] =  a[2].*b[1]
                memoryV[4] =  zero[1].*a[1].*b[1]

#                memoryV[4] =  a[2].*b[2]
                
#                memoryV[6] =  

#                memoryV[6] =  a[3].*b[1]
#                memoryV[7] =  a[1].*b[3]

#                memoryV[7] =  zero[1].*a[1].*b[1]
#                memoryV[8] =  zero[1].*a[1].*b[2]
#                memoryV[9] =  zero[1].*a[2].*b[1]
            catch
                println("asdf ",size(a), " " , size(b))
            end
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
        return memoryV * 10^3
    end

#=
#    return 0.0

    zero =  laguerre(dist0,missing, nmax=2, memory=memory0)
    a = laguerre(dist1,missing, nmax=2, memory=memory1)
    b = laguerre(dist2,missing, nmax=2, memory=memory2)

#    println(same_atom, "three_body_H ", zero, a,b)

#    zero = memory0
#    a = memory1
#    b = memory2
    
    if same_atom
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
#            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]    a[:,2].*b[:,1]    zero[:,1].*a[:,1].*b[:,1]  zero[:,2].*a[:,1].*b[:,1]  a[:,2].*b[:,2]   zero[:,1].*a[:,1].*b[:,2] zero[:,1].*a[:,2].*b[:,1]  ]

#            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]    a[:,2].*b[:,1] a[:,2].*b[:,2]   zero[:,1].*a[:,1].*b[:,1] ]

            Vt = [a[:,1].*b[:,1]  a[:,1].*b[:,2]    a[:,2].*b[:,1]    zero[:,1].*a[:,1].*b[:,1] ]

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

#                memoryV[5] =  zero[2].*a[1].*b[1]
#                memoryV[6] =  a[2].*b[2]

#                memoryV[7] =  zero[1].*a[1].*b[2]
#                memoryV[8] =  zero[1].*a[2].*b[1]
                

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
        return memoryV * 10^3
    end
=#
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
function calc_threebody_onsite(t1,t2,t3,dist12,dist13,dist23, database; set_maxmin=false, memory=missing)

    c = database[(t1,t2,t3)]


    #o1 = summarize_orb(orb1)

#    o2 = summarize_orb(orb2)    

#    sameat = (t1 == t2 && t1 == t3 )
    sameat = (t2 == t3)

#    indO = c.inds[[t1,t2,t3,o1,:O]]
    indO = c.inds[[t1,t2,t3,:O]]

    #    sameat = (t2 == t3 || t1 == t2 || t1 == t3)
    #sameat = (t2 == t3)
    #sameat = (t1 == t2 && t1 == t3 )

    Otot = three_body_O(dist12, dist13, dist23, sameat, c.datH[indO], memoryV=memory)


    return Otot
        
end

"""
    function fit_threebody_onsite(t1,t2,t3,orb1,dist12,dist13,dist23) 

Fit three body onsite interactions.
"""
function fit_threebody_onsite(t1,t2,t3,dist12,dist13,dist23)

#    o1 = summarize_orb(orb1)
#    o2 = summarize_orb(orb2)    

#    println("f3o $t1 $t2 $t3 $orb1  ", sum(dist12-dist13), " ", sum(dist12-dist23), " ",sum(dist13-dist23))

#    sameat = (t2 == t3 || t1 == t2 || t1 == t3)
#    sameat = (t2 == t3)
    
#    sameat = (t1 == t2 && t1 == t3)
    sameat = (t2 == t3)

    Otot = three_body_O(dist12, dist13, dist23, sameat)

#    if !sameat
#        Otot[1] = 0.0
#    end
    
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
        H =  three_body_H(dist, dist31, dist32,t1==t2, c.datH[ind], memory0=memory0, memory1=memory1, memory2=memory2, memoryV=memoryV)
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
    
    H =  three_body_H(dist, dist31, dist32, t1==t2)
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
            F_old = eigen(hk_old, sk_old)
            F_orth = eigen(hk_orth)
            F_new = eigen(hk_new, sk_new)
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
    F = eigen(hk_old, sk_old)
#    println("vals renorm ", real(F.values)                        )
#    println("sk_renorm ", diag(sk_old))
#    println("hk_renorm ", diag(hk_old))

    
#    trim_dist(tbc_renorm.tb, cutoff)

    tbc_renorm = trim_dist(tbc_renorm)

    vect_old, vals_old, hk_old, sk_old = Hk(tbc_renorm.tb,  [0,0,0])
    F = eigen(hk_old, sk_old)
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
    F = eigen(hk_old, sk_old)
#    println("vals shift ", real(F.values) )

    
    return tbc_renorm
    
end





end #end module


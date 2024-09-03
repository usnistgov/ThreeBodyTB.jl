###include("Crystal.jl")
###include("DFToutMod.jl")
#using XMLDict



"""
    module Classical

Run a CLASSICAL model based on similar computational machinery to the TB3 model.
"""
module Classical

using ..Utility:arr2str
using ..Utility:str_w_spaces
using ..Utility:parse_str_ARR_float
using ..Utility:dict2str
using ..Utility:str2tuplesdict

using Random
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
using ..CalcTB:calc_frontier
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

using ..TB:ewald_guess
using ..Ewald:estimate_best_kappa
#using ..ThreeBodyTB:scf_energy_force_stress
#using ..ThreeBodyTB:scf_energy

using Random


const n_2body_cl = 6

const n_2body_em = 6

const n_2body_charges = 7

#const n_2body_em = 12

#const n_3body_cl_same = 7
#const n_3body_cl_pair = 13
#const n_3body_cl_diff = 20

const n_3body_cl_same = 8
const n_3body_cl_pair = 15
const n_3body_cl_diff = 23


const a_var = 2.0
const L_A_0 = 1.0
const L_B_0 = 1.0
const L_C_0 = 1.0


#const n_3body_cl_same = 4
#const n_3body_cl_pair = 13
#const n_3body_cl_diff = 8



using ..CrystalMod:cutoff2X
using ..CrystalMod:cutoff3bX
using ..CrystalMod:cutoff_onX
using ..CrystalMod:cutoff_length
using ..CrystalMod:cutoff4X

using ..Utility:reshape_vec
using ..Atomdata:get_cutoff

using ..DFToutMod:dftout

#const cutoff2X = 18.51 
#const cutoff3bX = 13.01 
#const cutoff_onX = 18.01 

#const cutoff_length = 1.0


"""
    struct coefs_cl

Hold the classical coefficients for 2body or 3body interactions and EAM interactions.

"""
struct coefs_cl

    dim::Int64
    datH::Array{Float64,1}
    sizeH::Int64
    names::Set
    min_dist::Float64
    dist_frontier::Dict
    version::Int64
    em::Bool
    charges::Bool
    lim::Dict
    repval::Dict
end



"""
    function write_coefs_cl(filename, co::coefs; compress=true)

Write `coefs_cl` to an file xml. Compress uses gzip. See `read_coefs_cl`
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
#    addelement!(c, "cutoff", string(co.cutoff))
    addelement!(c, "min_dist", string(co.min_dist))
    addelement!(c, "dist_frontier", dict2str(co.dist_frontier))
    addelement!(c, "em", string(co.em))
    addelement!(c, "charges", string(co.charges))

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
    function read_coefs_cl(filename, directory = missing)

Read Classical `coefs_cl` from filename. Can read gzipped files directly.
"""
function read_coefs_cl(filename, directory = missing)
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
    em = parse(Bool, d["coefs"]["em"])
    charges = parse(Bool, d["coefs"]["charges"])

    version = parse(Int64, d["coefs"]["version"])
    
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


#    println("read ----------------")
#    println("names $names")
#    println("dim $dim")
#    println("datH $datH")
#    println("at_list, $at_list)")
#    println("min_dist $min_dist")
#    println("version $version")
#    println("em $em")
#    println("dist_frontier ", dist_frontier)
#    println("-")
    
    co = make_coefs_cl(names,dim, datH=datH, min_dist=min_dist, dist_frontier = dist_frontier, version=version, em=em,charges=charges, lim=lim, repval=repval)

    return co
    
end



"""
    function make_coefs_cl(at_list, dim)

Constructor for `coefs`. Can create coefs filled with ones for testing purposes.

See `coefs_cl` to understand arguments.
"""
function make_coefs_cl(at_list, dim; datH=missing, min_dist = 3.0, fillzeros=false, dist_frontier=missing, version=3, em=false, charges=false, lim=missing, repval=missing)

    if ismissing(lim)
        lim = Dict()
    end
    if ismissing(repval)
        repval = Dict()
    end


    println("make coefs ", at_list, " ", typeof(at_list))
    if em == true
        totH = n_2body_em
    elseif charges == true
        totH = n_2body_charges
    else
        if dim == 2
            totH = n_2body_cl

            
        elseif dim == 3
            totH = n_3body_cl_diff

            
        elseif dim == 4
            totH = 2
        else
            println("WARNING, make_coefs_cl only implemented up to 4")
            toth = 0
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
    end
    
    if ismissing(datH)
        if fillzeros
            datH = zeros(totH)
        else
            datH = ones(totH) 
        end
    end

    at_list = Symbol.(at_list)
    dist_frontier2 = Dict()
    if !ismissing(dist_frontier)
        println("not missing ")
        for key in keys(dist_frontier)
#            println("key ", key, " at_list ", at_list, " ", dim == length(key), " " , Set(at_list) == Set(key) )
            #println([dim, length(key), Set(Symbol.(at_list)),Set(key)])
            if dim == length(key) && Set(at_list) == Set(key)
                dist_frontier2[String.(key)] = dist_frontier[key]
                dist_frontier2[Symbol.(key)] = dist_frontier[key]
            end
        end
    end
    println("dim $dim")
    
    if dim == 2
        at_arr = Symbol.([i for i in at_list])
        if length(at_list) == 1
            at_arr = [(at_arr[1], at_arr[1])]
        else
            at_arr = [(at_arr[1], at_arr[2]), (at_arr[2], at_arr[1])]
        end


        for a in at_arr
            #println("a $a")
            println(typeof(a))
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
    
    println("dim $dim")
    println("datH $datH")
    println("totH $totH")
    println("at_list, $(Set(at_list))")
    println("min_dist $min_dist")
    println("version $version")
    println("em $em")
    println(typeof.([dim, datH, totH, Set(at_list), min_dist, dist_frontier2, version, em]))
    println("lim ", lim)
    println("repval ", repval)
    return coefs_cl(dim, datH, totH, Set(at_list), min_dist, dist_frontier2, version, em, charges, lim, repval)
    
end
    

Base.show(io::IO, d::coefs_cl) = begin

    println(io, "coeffs classical ", d.names, " embedding: ", d.em, " charges: ", d.charges)
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




function ham(x :: Vector, ct, database, donecheck, DIST, FloatX, use_threebody, dat_vars, at_types, vars_list, ind_set, turn_off_warn, verbose, use_fourbody, use_em, use_charges, fixed_charges, kappa, factor)
    T=eltype(x)

    x_r, x_r_strain = reshape_vec(x, ct.nat, strain_mode=true)
    A = FloatX.(ct.A) * (I(3) + x_r_strain)
    crys_dual = makecrys( A , ct.coords + x_r, ct.types, units="Bohr")

    energy, _ = calc_energy_cl(crys_dual;  database=database, DIST=DIST, verbose=verbose, use_threebody=use_threebody, var_type=T, dat_vars=dat_vars, at_types=at_types, vars_list=vars_list, ind_set=ind_set, turn_off_warn=turn_off_warn, use_fourbody=use_fourbody , use_em=use_em, use_charges=use_charges, check_frontier=false, fixed_charges=fixed_charges, kappa=kappa, factor=factor)

    return energy
    
end

"""
    function energy_force_stress_cl(crys::crystal)

Main function for running classical model. Returns energy/force/stress from crystal structure.
"""
function energy_force_stress_cl(crys::crystal;  database=missing, verbose=false, use_threebody=true, dat_vars=missing, at_types = missing, vars_list = missing,  DIST=missing, ind_set=missing, var_type=Float64, turn_off_warn = false, use_fourbody=false, use_em = true, use_charges=true, fixed_charges = missing, factor=1.0)

    if verbose; println("dist energy_force_stress");end
    if ismissing(DIST)
        if use_threebody
            if use_fourbody
                R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3, array_ind4 = distances_etc_3bdy_parallel_LV(crys,cutoff2X,cutoff3bX,var_type=var_type, return_floats=false, cutoff4 = cutoff4X)
                DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3, array_ind4
            else
                R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel_LV(crys,cutoff2X,cutoff3bX,var_type=var_type, return_floats=false)
                DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3
            end
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

    kappa = estimate_best_kappa(crys.A)

    FN_ham = x->ham(x,crys,database, true, DIST, var_type, use_threebody, dat_vars, at_types, vars_list, ind_set, turn_off_warn, verbose, use_fourbody, use_em, use_charges, fixed_charges, kappa, factor)

    if verbose; println("energy test "); end
    energy_test = FN_ham(zeros(var_type, 3*crys.nat + 6))

    #    return energy_test
#    println("energy $energy energy_test $energy_test")
                         
    chunksize=min(100, 3*crys.nat + 6)
    #println("chunksize $chunksize")
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

    #neaten (this breaks the fitting in very rare cases, do not use)
    #=            
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
    =#
    
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

function core_cl_charges(t1, t2, lag_arr, DAT_ARR_CHARGES, var_type )
    energy = zero(var_type)
    for n = 1:(n_2body_charges-1)
        energy +=  lag_arr[n]*DAT_ARR_CHARGES[t1,t2,n+1] #1 reserved for onsite
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
    @inbounds @fastmath begin
#        a_var=2.0

        ad_1 = a_var*dist_1
        expa_1 =exp.(-0.5*ad_1) 

#        L_A_0 = 1.0
        L_A_1 = expa_1
        L_A_2 = (1.0 .- ad_1).*expa_1
        L_A_3 = 0.5*(ad_1.^2 .- 4.0*ad_1 .+ 2.0) .* expa_1
        
        ad_2 = a_var*dist_2
        expa_2 =exp.(-0.5*ad_2) 

#        L_B_0 = 1.0
        L_B_1 = expa_2
        L_B_2 = (1.0 .- ad_2).*expa_2
        L_B_3 = 0.5*(ad_2.^2 .- 4.0*ad_2 .+ 2.0) .* expa_2
        
        ad_3 = a_var*dist_3
        expa_3 =exp.(-0.5*ad_3) 

#        L_C_0 = 1.0
        L_C_1 = expa_3
        L_C_2 = (1.0 .- ad_3).*expa_3
        L_C_3 = 0.5*(ad_3.^2 .- 4.0*ad_3 .+ 2.0) .* expa_3
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

#        energy += DAT_ARR3[24,t1,t2,t3] * L_A_0 * L_B_1 * L_C_3 #15-20  #k 
#        energy += DAT_ARR3[25,t1,t2,t3] * L_A_0 * L_B_3 * L_C_1         #l
#        energy += DAT_ARR3[26,t1,t2,t3] * L_A_3 * L_B_0 * L_C_1         #l
#        energy += DAT_ARR3[27,t1,t2,t3] * L_A_1 * L_B_0 * L_C_3         #k
#        energy += DAT_ARR3[28,t1,t2,t3] * L_A_1 * L_B_3 * L_C_0         #m
#        energy += DAT_ARR3[29,t1,t2,t3] * L_A_3 * L_B_1 * L_C_0         #m
        
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

    a=2.0

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

"""
Four body terms are experimental, may not be supported
"""
function core_4_cl_laguerre_fast_fourbdy_cl!(dist_1, dist_2, dist_3, dist_4, dist_5, dist_6, t1,t2,t3,t4, DAT_ARR4, var_type)


    begin
        a=2.0

        ad_1 = a*dist_1
        ad_2 = a*dist_2
        ad_3 = a*dist_3
        ad_4 = a*dist_4
        ad_5 = a*dist_5
        ad_6 = a*dist_6
        
        expa_1 =exp.(-0.5*ad_1)
        expa_2 =exp.(-0.5*ad_2)
        expa_3 =exp.(-0.5*ad_3)
        expa_4 =exp.(-0.5*ad_4)
        expa_5 =exp.(-0.5*ad_5)
        expa_6 =exp.(-0.5*ad_6)


#        L_A_0 = 1.0
#        L_B_0 = 1.0
#        L_C_0 = 1.0
#        L_D_0 = 1.0
#        L_E_0 = 1.0
#        L_F_0 = 1.0

        L0 = ones(6)
        
#        L_A_1 = expa_1
#        L_B_1 = expa_2
#        L_C_1 = expa_3
#        L_D_1 = expa_4
#        L_E_1 = expa_5
#        L_F_1 = expa_6

        L1 = [expa_1, expa_2, expa_3, expa_4, expa_5, expa_6]
        
#        L_A_2 = (1.0 .- ad_1)*expa_1
#        L_B_2 = (1.0 .- ad_2)*expa_2
#        L_C_2 = (1.0 .- ad_3)*expa_3
#        L_D_2 = (1.0 .- ad_4)*expa_4
#        L_E_2 = (1.0 .- ad_5)*expa_5
        #        L_F_2 = (1.0 .- ad_6)*expa_6

#        L2 = [(1.0 .- ad_1)*expa_1,(1.0 .- ad_2)*expa_2,(1.0 .- ad_3)*expa_3,(1.0 .- ad_4)*expa_4,(1.0 .- ad_5)*expa_5,(1.0 .- ad_6)*expa_6]
        
    end    

    @inbounds @fastmath begin
        energy = zero(var_type)
        
        #        for c in collect(combinations(1:6,4))
        for c in [ [1, 2, 3, 4], [1, 2, 3, 5], [1, 2, 3, 6], [1, 2, 4, 5], [1, 2, 4, 6], [1, 2, 5, 6], [1, 3, 4, 5], [1, 3, 4, 6], [1, 3, 5, 6], [1, 4, 5, 6], [2, 3, 4, 5], [2, 3, 4, 6], [2, 3, 5, 6], [2, 4, 5, 6], [3, 4, 5, 6]]
            energy += DAT_ARR4[1,t1,t2,t3,t4] * prod(L1[c])
        end

#        for c2 in collect(combinations(1:6,5))
        #            for c1 in collect(combinations(1:6,1))
        for c2 in [ [1, 2, 3, 4, 5], [1, 2, 3, 4, 6], [1, 2, 3, 5, 6], [1, 2, 4, 5, 6], [1, 3, 4, 5, 6], [2, 3, 4, 5, 6]]
            for c1 in 1:6
                energy += DAT_ARR4[2,t1,t2,t3,t4] * L0[c1] * prod(L1[c2])
            end
        end
        
    #    @inbounds @fastmath begin
#        energy = zero(var_type)
#        energy += DAT_ARR4[1,t1,t2,t3,t4] * L_A_1 * L_B_1 * L_C_1 * L_D_1
#        energy += DAT_ARR4[1,t1,t2,t3,t4] * L_A_1 * L_B_1 * L_C_1 * L_E_1
 #       energy += DAT_ARR4[1,t1,t2,t3,t4] * L_A_1 * L_B_1 * L_C_1 * L_F_1
 ##       energy += DAT_ARR4[1,t1,t2,t3,t4] * L_A_1 * L_B_1 * L_C_1 * L_D_1
 #       energy += DAT_ARR4[1,t1,t2,t3,t4] * L_A_1 * L_B_1 * L_C_1 * L_D_1
#
#        energy += DAT_ARR4[1,t1,t2,t3,t4] * L_A_2 * L_B_1 * L_C_1 * L_D_1
#        energy += DAT_ARR4[1,t1,t2,t3,t4] * L_A_1 * L_B_2 * L_C_1 * L_D_1
#        energy += DAT_ARR4[1,t1,t2,t3,t4] * L_A_1 * L_B_1 * L_C_2 * L_D_1
#        energy += DAT_ARR4[1,t1,t2,t3,t4] * L_A_1 * L_B_1 * L_C_1 * L_D_2
#
#        energy += DAT_ARR4[1,t1,t2,t3,t4] * L_A_0 * L_B_1 * L_C_1 * L_D_1
#        energy += DAT_ARR4[1,t1,t2,t3,t4] * L_A_1 * L_B_0 * L_C_1 * L_D_1
#        energy += DAT_ARR4[1,t1,t2,t3,t4] * L_A_1 * L_B_1 * L_C_0 * L_D_1
#        energy += DAT_ARR4[1,t1,t2,t3,t4] * L_A_1 * L_B_1 * L_C_1 * L_D_0
        
        
    end
    
    return energy
    

end


#------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
    function calc_energy_cl(crys::crystal

Main function for classical energy calculation from crystal structure.
"""
function calc_energy_cl(crys::crystal;  database=missing, dat_vars=missing, at_types = missing, vars_list = missing,  DIST=missing, verbose=false, use_threebody=true, ind_set=missing, var_type=missing, turn_off_warn = false, use_fourbody = false, use_em = true, use_charges= true, check_only=false, check_frontier=true, fixed_charges=missing, kappa=missing, factor=1.0)

    no_errors = true
    #fixed_charges .= [100.0, 100.0]
#    println("calc_energy_cl use_em $use_em use_threebody $use_threebody use_charges $use_charges fixed_charges $fixed_charges")
#    verbose=true


    #    println("CRYS.nat ", crys.nat)
    
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
    CHARGES = zeros(var_type, crys.nat)

    if use_charges
        if ismissing(fixed_charges) 
            CHARGES .= ewald_guess(crys, kappa=kappa, factor=factor)
        else
            CHARGES .= fixed_charges
        end
    end    
    #    if verbose  println("LV $var_type")  end
    
    
    

    
    ind2orb, orb2ind, etotal, nval = orbital_index(crys)
    if verbose println("distances") end
     begin
        
        use_dist_arr = true
        if !ismissing(DIST)
            use_dist_arr = false
            if use_fourbody
                R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3,array_ind4 = DIST
            else
                R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3 = DIST
            end                

        else
            if (use_threebody ) 
                if (use_fourbody)
                    R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3,array_ind4 = distances_etc_3bdy_parallel_LV(crys,cutoff2X, cutoff3bX,var_type=var_type, return_floats=false, cutoff4 = cutoff4X)
                    DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3, array_ind4
                else
                
                    R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel_LV(crys,cutoff2X, cutoff3bX,var_type=var_type, return_floats=false)
                    DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3
                end

            else
                R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel_LV(crys,cutoff2X, 0.0,var_type=var_type, return_floats=false)
                DIST = R_keep, R_keep_ab, array_ind3, c_zero, dmin_types, dmin_types3             
            end
        end

        
        
        within_fit = true
        
        if !ismissing(database)
            for key in keys(dmin_types)
                for key2 in keys(database)
                    if !(typeof(key) <: Set) || !(typeof(key2) <: Set)
                        continue
                    end
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

    if check_frontier
#        println("keys ", keys(database))
        violation_list, vio_bool, repel_vals = calc_frontier(crys, database, test_frontier=true, diststuff=DIST, verbose=verbose, var_type=var_type, use_threebody=use_threebody)
        if vio_bool
            no_errors = false
        end
        if check_only
            return vio_bool
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
         if use_em
            DAT_EM_ARR = zeros(var_type, types_counter, types_counter,1:n_2body_em )
         end
         if use_charges
             DAT_CHARGES_ARR = zeros(var_type, types_counter, types_counter,1:n_2body_charges )
         end
         
        #        DAT_IND_ARR3 = zeros(var_type, types_counter, types_counter, types_counter, 1: max(n_3body_cl_same,n_3body_cl_pair,n_3body_cl_diff) )
         if use_threebody
             DAT_IND_ARR3 = zeros(var_type, max(n_3body_cl_same,n_3body_cl_pair,n_3body_cl_diff), types_counter, types_counter, types_counter )
             if use_fourbody
                 DAT_IND_ARR4 = zeros(var_type, 2, types_counter, types_counter, types_counter, types_counter )
             end
         end
        
        if !ismissing(database)
#            println("not is missing ")
            #            if use_em
#                if :em in keys(database)
#                    em = database[:em]                    
#                else
#                    if use_em
#                        if !turn_off_warn; println("WARNING, no embed key :em "); end
#                        use_em = false
#                    end
#                end
#            end

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

                    if use_em
                        found_em=false
                        if (t1,t2,:em) in keys(database)
                            found_em = true
                            DAT_EM_ARR[c1,c2,1:n_2body_em] = database[(t1,t2, :em)].datH
                        end
                        if !found_em
                            use_em=false
                        end
                    end

                    if use_charges
#                        println("use ", keys(database))
                        found_charge = false
                        if (t1,t2,:charges) in keys(database)
#                            println("add ")
                            found_charge=true
                            DAT_CHARGES_ARR[c1,c2,1:n_2body_charges] = database[(t1,t2, :charges)].datH
                        end
                        if !found_charge
                            use_charges = false
                        end

                        
#                        println("DAT_CHARGES_ARR ", DAT_CHARGES_ARR)
#                        println("CHARGES $CHARGES")
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
                            if use_fourbody
                                for c4 = 1:types_counter
                                    t4 = types_dict_reverse[c4]
                                    if (t1,t2,t3,t4) in keys(database)
                                        coef = database[(t1,t2,t3,t4)]
                                        #                                DAT_IND_ARR3[c1,c2,c3,1:coef.sizeH] = coef.datH[1:coef.sizeH]
                                        DAT_IND_ARR4[1:coef.sizeH, c1,c2,c3,c4] = coef.datH[1:coef.sizeH]
                                    else
                                        if !turn_off_warn;                             println("WARNING, ",(t1,t2,t3,t4), " database not found"); end
                                        within_fit = false
                                        push!(badlist, (t1,t2,t3,t4))
                                        if !turn_off_warn; println("badlist ", badlist); end
                                    end
                                end
                            end
                        end
                    end
                    
                end
            end
        else

            for v in vars_list
#                println("v $v")
                if v[2] == 2
                    for c1 = 1:types_counter
                        for c2 = 1:types_counter
                            t1 = types_dict_reverse[c1]
                            t2 = types_dict_reverse[c2]
                            if (t1,t2) in keys(ind_set)
                                ind = ind_set[(t1,t2)] 
#                                println("ind $ind")
                                if t1 in v[1] && t2 in v[1] && length(Set([t1,t2])) == length(v[1])
                                #                            println("v $v")
                                #                            println("DAT_IND_ARR $c1 $c2 $n_2body_cl ", dat_vars[ind])
                                #                            println("ind $ind ", typeof(ind))
                                #                            println("dat_vars  $ind ", typeof(dat_vars[ind]))
                                    DAT_IND_ARR[c1,c2,1:n_2body_cl] = dat_vars[ind]
                                #                            println("after")
                                end
                            end

                            if use_em
#                                println("use ------------------------------------------ $((t1,t2,:em)) ")
#                                println(keys(ind_set), " " , (t1,t2,:em) in keys(ind_set))
                                if (t1,t2,:em) in keys(ind_set)
#                                    println("key")
                                    ind = ind_set[(t1,t2,:em)] 

                                    #                                    println(v , " " , v[1])
                                    if t1 in v[1] && t2 in v[1]
                                        DAT_EM_ARR[c1,c2,1:n_2body_em] = dat_vars[ind]
#                                        println("add $t1 $t2 xxxxxxxxxxxxxxxxxxxxxxxxxxx")
                                    end
                                end
                            end

                            if use_charges
#                                println("use charges")
                                if (t1,t2,:charges) in keys(ind_set)
#                                    println("key")
                                    ind = ind_set[(t1,t2,:charges)] 
#                                    println("ind charges $ind")
#                                    println(v , " " , v[1])
                                    if t1 in v[1] && t2 in v[1]
                                        DAT_CHARGES_ARR[c1,c2,1:length(dat_vars[ind])] = dat_vars[ind]
#                                        println("blah $((t1,t2,:charges)) ", dat_vars[ind])
                                    end
                                end
                            end

                            
                        end
                    end
#                    counter += n_2body_cl
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
#                                        println("dat_vars[ind] ", dat_vars[ind], " size ", size(DAT_IND_ARR3))
                                        DAT_IND_ARR3[1:length(ind), c1,c2,c3] = dat_vars[ind]
                                    end
                                end
                            end
                        end
                    end

                elseif v[2] == 4
                    for c1 = 1:types_counter
                        for c2 = 1:types_counter
                            for c3 = 1:types_counter
                                for c4 = 1:types_counter
                                    t1 = types_dict_reverse[c1]
                                    t2 = types_dict_reverse[c2]
                                    t3 = types_dict_reverse[c3]
                                    t4 = types_dict_reverse[c4]
                                    if t1 in v[1] && t2 in v[1] && t3 in v[1] && t4 in v[1] && length(Set([t1,t2,t3,t4])) == length(v[1])
                                        if (t1,t2,t3,t4) in keys(ind_set)
                                            ind = ind_set[(t1,t2,t3,t4)] 
                                            DAT_IND_ARR4[1:length(ind), c1,c2,c3,c4] = dat_vars[ind]
                                        end
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

    rho = zeros(var_type, crys.nat, types_counter)
    rho2 = zeros(var_type, crys.nat, types_counter)
#    rho3 = zeros(var_type, crys.nat, types_counter)
    if verbose println("2body CL") end
     twobody_Cl = begin
        for c = 1:nkeep_ab #rho not thread safe
            #id = threadid()
            id = 1
            lag_arr = lag_arr_TH[:,id]

            cind = R_keep_ab[c,1]
            cham = R_keep_ab[c,7]
            a1 = R_keep_ab[c,2]
            a2 = R_keep_ab[c,3]

#            println("c $c a1 $a1 a2 $a2 types_arr $types_arr")
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
                if use_charges
                    energy_2bdy_TH[id] += 0.5 * CHARGES[a1]^2 * DAT_CHARGES_ARR[t1,t2,1] #onsite charges
                end
                continue
            end
            
            laguerre_fast!(dist_a, lag_arr, a=2.0)
            energy_2bdy_TH[id] += core_cl(t1, t2, lag_arr, DAT_IND_ARR, var_type) * cut_a
            #            println("ADD ", core_cl(t1, t2, lag_arr, DAT_IND_ARR, var_type) * cut_a)
            if use_charges
                energy_2bdy_TH[id] += core_cl_charges(t1, t2, lag_arr, DAT_CHARGES_ARR, var_type) * cut_a *  CHARGES[a1] * CHARGES[a2] 
            end
            rho[a1,t2] += exp(-1.0*dist_a)
            rho2[a1,t2] += (1.0 .- dist_a * 2.0)* exp(-1.0*dist_a)
#            rho3[a1, t2] +=0.5*((dist_a * 2.0).^2 .- 4.0*dist_a * 2.0 .+ 2)*exp(-1.0*dist_a)
        end
        
    end
    
    energy_embed = zero(var_type)
    if verbose; println("do em"); end
    if use_em
        for a1 = 1:crys.nat
#            if var_type == Float64
#                println("rho $a1 $(rho[a1])  $( em[1]*rho[a1] + em[2]*rho[a1]^2 + em[3]*rho[a1]^3)")
#            end
            t1 = types_arr[a1]
            for t=1:types_counter
                em = DAT_EM_ARR[t1,t,1:n_2body_em]
                energy_embed +=  em[1]*rho[a1,t]^2 + em[2]*rho[a1,t]^3  #+ em[6] * rho[a1,t]^1
                energy_embed +=  em[3]*rho2[a1,t]^2 + em[4]*rho2[a1,t]^3  #+ em[7] * rho2[a1,t]^1
                energy_embed +=  em[5]*rho[a1,t]*rho2[a1,t] #+ em[6]*(rho[a1,t]*rho2[a1,t])^2
                if t1 != t
                    energy_embed +=  em[6]*rho[a1,t]*rho[a1,t1] 
                end

#                energy_embed +=  em[6]*rho3[a1,t]^2 + em[7]*rho3[a1,t]^3  + em[12]* rho3[a1,t]^4
#                energy_embed +=  em[8]*rho[a1,t]*rho3[a1,t] + em[9]*rho2[a1,t]*rho3[a1,t]
                
            end
            
            #            energy_embed +=  em[7]*rho[a1]*rho2[a1] + em[8]*rho2[a1]*rho3[a1]

            
        end
    end
        

    if verbose println("3body LV") end
     begin
        lag_arr3_TH = zeros(var_type, max(n_3body_cl_same, n_3body_cl_pair, n_3body_cl_diff), nthreads())
        energy_3bdy_TH = zeros(var_type,nthreads())
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
            
            #            println("asdf")
            
            #println("meta_count $meta_count")
            
            @threads for mc in meta_count 
                for counter in mc
                    

                    a1 = array_ind3[counter,1]
                    a2 = array_ind3[counter,2]
                    a3 = array_ind3[counter,3]


                    factor = 1.0
                    if a1 != a2 && a1 != a3 && a2 != a3
                        if a1 < a2 && a1 < a3 && a2 < a3
                            factor = 6.0
                        else
                            continue
                        end
                    end

                    id = threadid()
                    
#                    elseif a1 < a2 && a2 == a3
#                        factor = 2.0
#                    elseif a1 == a2 && a2 < a3
#                        factor = 2.0
#                    elseif a1 == a2 && a1 == a3
#                        factor = 1.0
#                    else
#                        continue
#                    end
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

                    energy_3bdy_TH[id] += core_3_cl_laguerre_fast_threebdy_cl!(dist12, dist13, dist23, t1,t2,t3, DAT_IND_ARR3, var_type)  * cut_h * 100.0 * factor
                    
                end
            end
        end
    end

    energy_4bdy_TH = zeros(var_type,nthreads())
    if verbose println("4body LV") end

     fourbody_LV = begin
        
        if use_fourbody

            meta_count = []
            old = 1
            for  counter = 1: (size(array_ind4)[1]-1 ) 
                if array_ind4[counter,1] != array_ind4[counter+1,1]
                    push!(meta_count, old:counter)
                    old = counter+1
                end
            end
            push!(meta_count, old:size(array_ind4)[1])
            
            #            println("asdf")
            
            #println("meta_count $meta_count")
            
            @threads for mc in meta_count 
                for counter in mc
                    
                    id = threadid()

                    a1 = array_ind4[counter,1]
                    a2 = array_ind4[counter,2]
                    a3 = array_ind4[counter,3]
                    a4 = array_ind4[counter,4]

            #        println("a1 a2 a3 $a1 $a2 $a3 ", crys.stypes)

                    t1s = crys.stypes[a1]
                    t2s = crys.stypes[a2]
                    t3s = crys.stypes[a3]
                    t4s = crys.stypes[a4]

                    t1 = types_arr[a1]
                    t2 = types_arr[a2]
                    t3 = types_arr[a3]
                    t4 = types_arr[a4]
                    
                    if !turn_off_warn && (t1,t2,t3,t4) in badlist
                        if warned == false
                            println("WARNING missing $( (t1,t2,t3,t4) ) ")
                            warned = true
                        end
                        continue
                    end
                    
                    cind1 = array_ind4[counter,5]
                    cind2 = array_ind4[counter,6]
                    cind3 = array_ind4[counter,7]

                    rind1 = ind_arr[cind1,1:3]
                    rind2 = ind_arr[cind2,1:3]
                    rind3 = ind_arr[cind3,1:3]

                    
                    
                    if false
                        dist12 = array_floats3[counter, 1]
                        dist13 = array_floats3[counter, 2]
                        dist23 = array_floats3[counter, 3]                    
                        
                        cut_h = array_floats3[counter,14]
                    else
                        dist12 = get_dist_only(a1,a2, rind1, crys, At)
                        dist13 = get_dist_only(a1,a3, rind2, crys, At)
                        dist14 = get_dist_only(a1,a4, rind3, crys, At)
                        
                        dist23 = get_dist_only(a2,a3, -rind1+rind2, crys, At)
                        dist24 = get_dist_only(a2,a4, -rind1+rind3, crys, At)

                        dist34 = get_dist_only(a3,a4, -rind2+rind3, crys, At)
                        
                        cutoff4 = cutoff4X
                            
#                        cut_ab = cutoff_fn_fast(dist12, cutoff3 - cutoff_length, cutoff3)
#                        cut_ac = cutoff_fn_fast(dist13, cutoff3 - cutoff_length, cutoff3)
#                        cut_bc = cutoff_fn_fast(dist23, cutoff3 - cutoff_length, cutoff3)
                        #                        cut_h = cut_ab*cut_ac*cut_bc
                        
                    end
                    
                    
                    #lag_arr3 = lag_arr3_TH[:,id]
                    #ntot = 20
#                    ntot = laguerre_fast_threebdy_cl!(dist12,dist13,dist23, t1,t2,t3, lag_arr3)
                    #                    energy_3bdy_TH[id] += core_3_cl(t1,t2,t3,ntot,lag_arr3,DAT_IND_ARR3, var_type) * cut_h * 100.0

                    energy_4bdy_TH[id] += core_4_cl_laguerre_fast_fourbdy_cl!(dist12, dist13, dist14, dist23,dist24, dist34, t1,t2,t3,t4, DAT_IND_ARR4, var_type) * 1000.0
                    
                end
            end
        end
    end

#    if var_type == Float64
#        println(" energy comps  $(sum(energy_2bdy_TH))    $(sum(energy_3bdy_TH))    $(sum(energy_4bdy_TH))    $energy_embed")
#    end
    if verbose; println("done e cl");end
    return sum(energy_2bdy_TH) + sum(energy_3bdy_TH) + sum(energy_4bdy_TH) + energy_embed, no_errors
    
    
end


end #end module




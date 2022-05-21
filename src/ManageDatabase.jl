"""
    module ManageDatabase

Module for reading `coefs` from files and making database as needed for calculations
"""
module ManageDatabase

#using FileIO
#using JLD2

using ..CrystalMod:crystal
using ..CalcTB:coefs

using ..ThreeBodyTB:DATSDIR1
using ..ThreeBodyTB:DATSDIR2
using ..CalcTB:read_coefs

datdir1 = DATSDIR1 #primary directory
datdir2 = DATSDIR2 #backup directory

database_list = Set()

database_cached = Dict()
database_cached["SCF"] = true
database_cached["scf"] = true

"""
    function prepare_database(c::crystal)

Get ready database of precalculated `coefs` for `crystal`
"""
function prepare_database(c::crystal; directory = missing)
#    println("prepare c ", c.types)
    prepare_database(c.types, directory = directory)
end

"""
    function clear_database()

Remove loaded databasea
"""
function clear_database()

    if length(database_list) > 0
        for a in database_list
            delete!(database_list, a)
        end
    end

    for a in keys(database_cached)
        delete!(database_cached, a)
    end
    database_cached["SCF"] = true
    database_cached["scf"] = true


end


"""
    function prepare_database(at_list)
"""
function prepare_database(at_list; directory = missing, verbose=true)
    
    if verbose; println("prepare atoms ", at_list); end
    at_list = Symbol.(at_list)
    #    println("database_list ", database_list)
    s = Set(at_list)
    for s1 in s
        for s2 in s
            for s3 in s
                if ! ( Set([s1,s2,s3]) in database_list)
                    add_to_database(Set([s1,s2,s3]), directory = directory, verbose=verbose)
                end
            end
        end
    end
end

"""
    function add_to_database(s::Set)

Load elements or twobody terms from precalcuated `coefs` from files.
"""
function add_to_database(s::Set; directory = missing, verbose=true)#

    if verbose; println("add_to_database $s"); end

    dirlist = []
    if !ismissing(directory)
        push!(dirlist, directory)
    end
    push!(dirlist, datdir1)
    push!(dirlist, datdir2)
    

    at_arr = collect(s)
    if length(s) == 1
        a1 = Symbol(at_arr[1])
        if !haskey(database_cached , (a1, a1))
            loaded = false
            for d in dirlist
                f =  "$d/els/coef.el.2bdy.$a1.xml.gz"

                if isfile(f) || isfile(f*".gz")
                    try
                        dat = read_coefs(f)
                        database_cached[(a1, a1)] = dat
                        println("added to cache ", (a1, a1))
                        loaded=true
                    catch
                        println("WARNING - error loading $f")
                    end
                end
                if loaded == true
                    break
                end

            end
            if loaded==false
                println("WARNING, FAILED TO LOAD 2bdy ", a1)
            end

        end

        if !haskey(database_cached , (a1, a1, a1))
#            f = "$defaultdatdir/v0.1_dat_3body_scf_pbesol_el.$a1.xml"
#            f =  "$datdir1/coef.el.3bdy.$a1.xml.gz"
#            f2 = "$datdir2/coef.el.3bdy.$a1.xml.gz"

            loaded = false
            for d in dirlist
                f = "$d/els/coef.el.3bdy.$a1.xml.gz"

                if isfile(f) || isfile(f*".gz")
                    try
                        #                    jldopen(f)                    
                        dat = read_coefs(f)
                        database_cached[(a1, a1, a1)] = dat
                        println("added to cache ", (a1, a1, a1))
                        loaded == true
                    catch
                        println("WARNING - error loading $f")
                    end
                    loaded=true
                end
                if loaded == true
                    break
                end

            end
            if loaded == false
                println("WARNING, FAILED LOADING 3bdy ", a1)
            end

#                elseif isfile(f2) || isfile(f2*".gz")
#                try
#                    jldopen(f)                    
#                    dat = read_coefs(f2)
#                    database_cached[(a1, a1, a1)] = dat
#                    println("added to cache ", (a1, a1, a1))
#                catch
#                    println("WARNING - error loading $f2")
#                end
#            else
#                println("WARNING, no file for 3bdy database $a1")
#            end
        end
           
 
    elseif length(s) == 2

        a1 = Symbol(at_arr[1])
        a2 = Symbol(at_arr[2])
        if !haskey(database_cached , (a1, a2))
            
#            fab = "$defaultdatdir/v0.1_dat_2body_scf_pbesol_binary.$a1.$a2.xml"
#            fba = "$defaultdatdir/v0.1_dat_2body_scf_pbesol_binary.$a2.$a1.xml"

#            fab =  "$datdir1/coef.el.2bdy.$a1.$a2.xml.gz"
#            fba =  "$datdir1/coef.el.2bdy.$a2.$a1.xml.gz"

#            fab2 = "$datdir2/coef.el.2bdy.$a1.$a2.xml.gz"
#            fba2 = "$datdir2/coef.el.2bdy.$a2.$a1.xml.gz"

            loaded = false
            for d in dirlist

                fab =  "$d/binary/coef.el.2bdy.$a1.$a2.xml.gz"
                fba =  "$d/binary/coef.el.2bdy.$a2.$a1.xml.gz"

                f=missing
                if isfile(fab) || isfile(fab*".gz")
                    f = fab
                elseif isfile(fba) || isfile(fba*".gz")
                    f = fba
                end
#            elseif isfile(fab2) || isfile(fab2*".gz")
#                f = fab2
#            elseif isfile(fba2) || isfile(fba2*".gz")
#                f = fba2
#            else
#                f = missing
#                println("WARNING - binary file missing ")
#            end

                if !ismissing(f)
                    try
                        dat = read_coefs(f)
                        
                        
                        database_cached[(a1, a2)] = dat
                        database_cached[(a2, a1)] = dat
                        println("added to cache ", (a1, a2), " twobody ")
                        loaded = true
                    catch
                        println("WARNING - error loading binary file $f")
                    end
                end
                if loaded == true
                    break
                end
            end
            if loaded == false
                println("WARNING, FAILED LOAD 2BDY ", a1, " ", a2)
            end

#################

#            fab = "$defaultdatdir/v0.1_dat_3body_scf_pbesol_binary.$a1.$a2.xml"
#            fba = "$defaultdatdir/v0.1_dat_3body_scf_pbesol_binary.$a2.$a1.xml"

#            fab =  "$datdir1/coef.el.3bdy.$a1.$a2.xml.gz"
#            fba =  "$datdir1/coef.el.3bdy.$a2.$a1.xml.gz"

#            fab2 = "$datdir2/coef.el.3bdy.$a1.$a2.xml.gz"
#            fba2 = "$datdir2/coef.el.3bdy.$a2.$a1.xml.gz"

            loaded = false
            for d in dirlist
                fab =  "$d/binary/coef.el.3bdy.$a1.$a2.xml.gz"
                fba =  "$d/binary/coef.el.3bdy.$a2.$a1.xml.gz"

                f=missing
                if isfile(fab) || isfile(fab*".gz")
                    f = fab
                elseif isfile(fba) || isfile(fba*".gz")
                    f = fba
                end
#                elseif isfile(fab2) || isfile(fab2*".gz")
#                    f = fab2
#                elseif isfile(fba2) || isfile(fba2*".gz")
#                f = fba2

#            else
#                f = missing
 #               println("WARNING - binary file missing ")
 #           end
                
                if !ismissing(f)
                    try
#                    jldopen(f)
                        dat = read_coefs(f)
                        
                        database_cached[(a1,a1,a2)] = dat
                        database_cached[(a1,a2,a1)] = dat
                        database_cached[(a2,a1,a1)] = dat
                        
                        database_cached[(a1,a2,a2)] = dat
                        database_cached[(a2,a1,a2)] = dat
                        database_cached[(a2,a2,a1)] = dat
                        println("added to cache ", (a1, a2), " threebody ")
                        loaded=true
                    catch
                        println(" MISSING FILE - WARNING loading binary 3body $f ")
                    end
                end
                if loaded == true
                    break
                end

            end
            if loaded==false
                println("WARNING, FAILED TO LOAD 3bdy ", a1, " " , a2)
            end
            
        end
    elseif length(s) == 3

        println("warning, ternary not currently supported")
        return

        a1 = Symbol(at_arr[1])
        a2 = Symbol(at_arr[2])
        a3 = Symbol(at_arr[3])

        ats  = sort([a1,a2,a3])

        a1 = ats[1]
        a2 = ats[2]
        a3 = ats[3]

        if !haskey(database_cached , Tuple(ats))

            loaded = false
            for d in dirlist
            
                f =  "$d/tern_$a1/coef.el.3bdy.$a1.$a2.$a3.xml.gz"
                if  isfile(f) || isfile(f*".gz")
                    try
                        dat = read_coefs(f)
                        database_cached[(a1,a2,a3)] = dat
                        database_cached[(a1,a3,a2)] = dat
                        database_cached[(a2,a1,a3)] = dat
                        database_cached[(a2,a3,a1)] = dat
                        database_cached[(a3,a1,a2)] = dat
                        database_cached[(a3,a2,a1)] = dat
                        println("added to cache ", (a1, a2, a3), " threebody ")
                        loaded=true
                    catch
                        println("failed to load $f")
                    end
                end
                if loaded == true
                    break
                end
            end
            if loaded == false
                println("WARNING, FAILED TO LOAD ", [a1, a2, a3])
            end
        end

    else
        println("WARNING - not setup for quaternary+ ", s)
    end

    
    push!(database_list, s)

    
end





end #end module

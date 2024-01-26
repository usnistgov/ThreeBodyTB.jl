
module CrystalMod

#using AtomsBase
#using Unitful
#using UnitfulAtomic

using LinearAlgebra
using Printf
using ..Atomdata:atoms
using ..Atomdata:atom_radius
using ..Atomdata:cutoff_dist
using ..Atomdata:charge_cell_missing_atom_correction
using GZip

using ..ThreeBodyTB:global_length_units
using ..ThreeBodyTB:convert_length
using ..ThreeBodyTB:convert_force
using ..ThreeBodyTB:convert_stress
using ..ThreeBodyTB:Ang

import Plots
import Plots:plot

export crystal
export makecrys
export generate_supercell
export generate_random_distortion
export write_poscar
export write_efs
export plot

#holds data
"""
    mutable struct crystal{T}

Holds basic crystal structure information, type T. Use `makecrys` to easily construct.

Note: you can create supercells like

```julia-repl
julia> c = makecrys([5.0 0 0; 0 5.0 0; 0 0 5.0], [0.0 0.0 0.0], ["H"])
A1=     5.00000  0.00000  0.00000
A2=     0.00000  5.00000  0.00000
A3=     0.00000  0.00000  5.00000

H    0.00000  0.00000  0.00000


julia> c*[2,2,2]
A1=     10.00000  0.00000  0.00000
A2=     0.00000  10.00000  0.00000
A3=     0.00000  0.00000  10.00000

H    0.00000  0.00000  0.00000
H    0.00000  0.00000  0.50000
H    0.00000  0.50000  0.00000
H    0.00000  0.50000  0.50000
H    0.50000  0.00000  0.00000
H    0.50000  0.00000  0.50000
H    0.50000  0.50000  0.00000
H    0.50000  0.50000  0.50000
```

# Holds
- `A::Array{T,2}` 3 × 3 lattice vectors, Bohr (atomic units) internally.
- `coords::Array{T,2}` num_atoms × 3  atomic positions, fractional units.
- `types::Array{String,1}` atomic names, like `"H"` or `"Zn"`.
- `types::Array{Symbol,1}` atomic names, but julia Symbols like `:H` or `:Zn`, for nominally faster internal evaluation.
- `nat::Int64` number of atoms.
"""
mutable struct crystal{T}

    A::Array{T,2}
    coords::Array{T,2}
    types::Array{String,1}
    stypes::Array{Symbol,1}
    nat::Int64

end

include("Dist.jl")


#printing
Base.show(io::IO, c::crystal) = begin

    println(io, "Units: "* global_length_units)
    println(io)
    if global_length_units == "Å"
        At = c.A * Ang
    else
        At = c.A
    end
        
    for i in 1:3
            @printf(io, "A%.1i=     %.5f  %.5f  %.5f\n", i, At[i,1], At[i,2], At[i,3])
    end
    @printf(io, "\n")
    for i in 1:c.nat
        @printf(io, "%-3s  %.5f  %.5f  %.5f\n", c.types[i], c.coords[i,1], c.coords[i,2], c.coords[i,3])        
    end
    
end   

"""
    function print_with_force_stress(c::crystal, force, stress)

pretty io function
"""
function print_with_force_stress(c::crystal, force, stress)

#    println(io, "Units: "* global_length_units)
    #    println(io)

    force = convert_force(force)
    stress = convert_stress(stress)
    
    if global_length_units == "Å"
        At = c.A * Ang
    else
        At = c.A
    end
        
    println("Lattice Vectors                       | Stress")
    println("-----------------------------------------------------------------------")
    for i in 1:3
            @printf( "A%.1i=     %+.5f  %+.5f  %+.5f  |  %+.5f  %+.5f  %+.5f\n", i, At[i,1], At[i,2], At[i,3], stress[i,1],stress[i,2], stress[i,3] )
    end
    @printf( "\n")
    println("Crystal coords                        | Force (Cartesian)")
    println("-----------------------------------------------------------------------")
    for i in 1:c.nat
        @printf( "%-3s     %+.5f  %+.5f  %+.5f  |  %+.5f  %+.5f  %+.5f\n", c.types[i], c.coords[i,1], c.coords[i,2], c.coords[i,3], force[i,1],force[i,2],force[i,3])        
    end
    println("-----------------------------------------------------------------------")
    println()
end   


Base.:*(c::crystal,strain::String) = begin

    println("string version")
    
end

Base.:*(c::crystal,strain::Array{<:Real,1}) = begin
    if size(strain) == (3,)
        A_new = c.A*(I(3) + diagm(strain))
        return makecrys(A_new, c.coords, c.types, units=:Bohr)
    else
        println("we cannot interpret as strain")
        println(strain)
        return c
    end
end

Base.:*(strain::Array{<:Real,1}, c::crystal) = begin
    return c*strain 
end

Base.:*(c::crystal,strain::Array{<:Real,2}) = begin
    if size(strain) == (3,3)
        A_new = c.A*(I(3) + strain)
        return makecrys(A_new, c.coords, c.types, units=:Bohr)
    elseif size(strain) == (1,3)
        A_new = c.A*(I(3) + diagm(strain[:]))
        return makecrys(A_new, c.coords, c.types, units=:Bohr)
    else
        println("we cannot interpret as strain")
        println(strain)
        return c
    end
end

Base.:*(strain::Array{<:Real, 2}, c::crystal) = begin
    return c*strain 
end


Base.:+(c1::crystal, c2::crystal) = begin

    if c1.nat != c2.nat
        println("nat not same ", c1.nat, c2.nat)
        return c1
    end
    A = c1.A .+ c2.A
    coords = c1.coords .+ c2.coords

    types = deepcopy(c1.types)
    
    return makecrys(A, coords, types, units="Bohr")
end   

Base.:(==)(c1::crystal, c2::crystal) = begin

    if c1.nat != c2.nat
        return false
    end
    
    if c1.types != c2.types
        return false
    end
    
    if sum(abs.(c1.A - c2.A)) > 1e-11
        return false
    end
    
    if sum(abs.(c1.coords - c2.coords)) > 1e-11
        return false
    end
    
    return true

end   


Base.:*(a::Real, c::crystal) = begin

    A = a * c.A 
    coords = deepcopy( c.coords)

    types = deepcopy(c.types)
    
    return makecrys(A, coords, types, units="Bohr")
end   

Base.:*(c::crystal,a::Real) = begin
    return a*c
end

Base.:*(a::Vector{<:Int}, c::crystal) = begin

    return generate_supercell(c, a)
    
#=    A_new = zeros(3,3)
    A_new[1,:] = c.A[1,:] *a[1]
    A_new[2,:] = c.A[2,:] *a[2]
    A_new[3,:] = c.A[3,:] *a[3]

    n = prod(a)
    coords_new = zeros(n*c.nat, 3)
    types_new = String[]

    counter = 0
    for x = 1:a[1]
        for y = 1:a[2]
            for z = 1:a[3]
                for i = 1:c.nat
                    counter +=1
                    coords_new[counter,:] = (c.coords[i,:] + [x-1, y-1, z-1]) ./ [a[1], a[2], a[3]]
                    push!(types_new, c.types[i])
                end
            end
        end
    end

    return makecrys(A_new, coords_new, types_new, units="Bohr")
    =#
    
end   

Base.:*(c::crystal,  a::Vector{<:Int} ) = begin
    return a*c
end


#setup crystal
"""
    makecrys(A,coords,types; units=missing)

Return a crystal object from 3×3 lattice, nat × 3 coords in crystal units, and nat element strings.

Will use units set by `set_units`, with default to Ang. Can override with `units`.

Note: also export-ed directly from ThreeBodyTB for convenience


```julia-repl
julia> makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0.0 0.0 0.0], ["H"])
A1=     10.00000  0.00000  0.00000
A2=     0.00000  10.00000  0.00000
A3=     0.00000  0.00000  10.00000

H    0.00000  0.00000  0.00000
```
"""
function makecrys(A,coords,types; units=missing, type=missing)
    if ismissing(type)
        T = typeof(A[1,1])
    else
        T = type
    end

    
    if !ismissing(units)
        if units == "A" || units == "Ang" || units == "Ang." || units == "Angstrom" || units == "a" || units == "ang" || units == "ang." || units == "angstrom" || units == "Å"
            factor = 1.0 / Ang
        else
            factor = 1.0
        end
    else
        if global_length_units == "Å"
            factor = 1.0 / Ang
        else
            factor = 1.0
        end
    end

    if isa(types, Tuple)
        types = collect(types)
    end
    
    #entering integer coords or A messes everything up later
    if T == Int64
        coords = Float64.(coords)
        T = typeof(coords[1,1])
    end
    if typeof(A[1,1]) == Int64
        A = Float64.(A) 
    end

    A = A * factor

    if T == Float64 || T == Float32
        if rank(A) != 3
            println("WARNING, possible linearly degenerate lattice vectors ....")
        end
    end
    
    if isa(types,Array) && length(size(types)) == 2  # if are accidently given a 1 x nat types array 
        types = types[:]
    end

    if typeof(types[1]) == Symbol
        types = String.(types)
    end

    for t in types
        if !(t in keys(atoms))
            println()
            println("WARNING, atom $t not currently supported --------------")
            periodictable()
            break
        end
    end
    
    nat = length(types)

    if nat != size(coords,1) #try to fix
        coords = coords'
    end

    if nat != size(coords,1)

        println(types)
        println(coords)
        error("Error initilzing crystal, nat doesn't match, nat from types: ", nat," size(coords,1): ", size(coords,1)  )
        
        return 
    end

    if (3,3) != size(A)
        println(A)
        error("Error initilzing crystal, A wrong ", size(A))
        
        return 
    end

    if 3 != size(coords,2)
        println(coords)
        error("Error initilzing crystal, coords wrong ", size(coords))
        
        return 
    end
    
    A = convert(Array{T,2}, A)
    coords = convert(Array{T,2}, coords)


#    coords = mod.(coords, 1.0)
    
    if size(coords)[1] != length(types)
        error("Error making crys, types and coords sizes don't match")
        return 
    end

    nat = size(coords)[1]
    stypes = Symbol.(types)
    return crystal{T}(A, coords, types, stypes, nat)
end

"""
   function set_types(c::crystal, types)

Helper function to keep crystal structure but change atoms on sites.
cd Types is a new set of atom types like `["Si", "Ge"]` or `[:Si, :Ge]`
Must match number of atoms in `c`
"""
function set_types(c::crystal, types)
    makecrys(c.A, c.coords, types)
end


#=
"""
    function makecrys(atomsb::FlexibleSystem)

Attempt to make a crystal from an AtomsBase object. This may not always work. Must be periodic
"""
function makecrys(atomsb::FlexibleSystem)


    thetype = typeof(atomsb.box[1][1].val)
    
    A = zeros(thetype, 3,3)
    for i = 1:3
        for j = 1:3
            x = atomsb.box[i][j]
            try
                A[i,j] = UnitfulAtomic.auconvert(u"Å",x )
            catch
                A[i,j] = x.val * 0.529177
            end
        end
    end

    stypes = []
    pos = zeros(Float64, 0,3)
    parts = atomsb.particles
    for p in parts
        push!(stypes , p.atomic_symbol)
        pos = vcat(pos, [0.0  0.0  0.0])
        try
            
            pos[end,1] = ustrip(Unitful.uconvert(u"Å",p.position[1]))
            pos[end,2] = ustrip(Unitful.uconvert(u"Å",p.position[2]))
            pos[end,3] = ustrip(Unitful.uconvert(u"Å",p.position[3]))
        catch
            pos[end,1] = ustrip(p.position[1])*0.529177
            pos[end,2] = ustrip(p.position[2])*0.529177
            pos[end,3] = ustrip(p.position[3])*0.529177
        end
            
    end

    pos = pos*inv(A)
    
    return makecrys(A, pos, stypes, units="Å")

end
=#            

                  



"""
    makecrys(filename::String)

Read filename, return crystal object.
File can be POSCAR, or simple quantum espresso inputfile

The entire file can be in the string instead of the filename. The code decides which is which by looking for newlines"
"""
function makecrys(filename::String)

    if !occursin("\n", filename)
        if !isfile(filename)
            println("tried to open $filename to load crystal, file does not exist")
            throw(error("noFile"))
        end

        f = gzopen(filename, "r")
        lines = readlines(f)
        close(f)
    
        c = makecrys(lines)
        return c
    else
        println("newlines found, interpret as string containing file instead of filename")
        lines = String.(split(filename, "\n"))

        c = makecrys(lines)
    end
    
end


"""
    makecrys(lines::Array{String,1})

Read string array, return crystal object.
string array can be POSCAR, or simple quantum espresso inputfile
"""
function makecrys(lines::Array{String,1})

    intype="POSCAR"
    for line in lines
        sp = split(line)
        if size(sp)[1] >= 1 && sp[1] == "&control"
            intype="QE"
        end
        if  size(sp)[1] >= 1 && sp[1] == "A1="
            intype="MYFORMAT"
        end
    end

    if intype == "POSCAR"

#        println("parse POSCAR")
        A, coords, types = parsePOSCAR(lines)

    elseif intype == "QE"

        A, coords, types = parseQEinput(lines)

    elseif intype == "MYFORMAT"

        A, coords, types = parseMYFORMATinput(lines)
        if global_length_units=="Å"
            A = A / Ang
        end

    end



    coords = neaten_coords(coords)
    coords = neaten_coords(coords)
    A = neaten_lattice(A)
    A = neaten_lattice(A)
    A = neaten_lattice(A)

    
    c = makecrys(A, coords,types, units="Bohr")
    return c
    
    end

function parseARRfloat(sp)

    return map(x->parse(Float64,x),sp)

end

function parseARRint(sp)

    return map(x->parse(Int64,x),sp)

end

function parseMYFORMATinput(lines)

    A = zeros(3,3)

    coords = zeros(1000, 3)
    types = []
    nat = 0 
    for line in lines
        sp = split(line)
        if length(sp) >= 1
            if sp[1] == "A1="
                A[1,:] = parseARRfloat(sp[2:4])
            elseif sp[1] == "A2="
                A[2,:] = parseARRfloat(sp[2:4])
            elseif sp[1] == "A3="
                A[3,:] = parseARRfloat(sp[2:4])
            elseif sp[1] == "Crystal" || sp[1] == "Lattice" 
                continue
            elseif length(sp) >= 4
                push!(types, String(sp[1]))
                nat += 1
                coords[nat,:] = parseARRfloat(sp[2:4])
            end
        end
    end

    coords = coords[1:nat, :]
    
    return A, coords, types

end

"""
    function parsePOSCAR(lines)

Parse a POSCAR from VASP

Called by makecrys, doesn't need to be called directly.
"""
function parsePOSCAR(lines)

    
    title = lines[1]

    
    a = parse(Float64, split(lines[2])[1])

    A = zeros(3,3)
    A[1,:] = parseARRfloat(split(lines[3]))
    A[2,:] = parseARRfloat(split(lines[4]))
    A[3,:] = parseARRfloat(split(lines[5]))

    A = A * a / 0.529177 

    thetypes = split(lines[6])

    

    thenumbers = parseARRint(split(lines[7]))


    
    
    if size(thetypes)[1] != size(thenumbers)[1]

        return -1,-1,-1
    end
    
    nat = sum(thenumbers)


    dc = split(lines[8])[1][1]
    cart = false
    if dc == 'c' || dc == 'C' 
        cart = true
    end

    
    types = String[]
    for (t,n) in zip(thetypes,thenumbers)

        for i = 1:n
            push!(types, t)
        end
    end
    
    coords = zeros(nat,3)
    for i = 1:nat

        coords[i,:] = parseARRfloat(split(lines[8+i])[1:3])
    end

    if cart
        coords = (coords / 0.529177 )* inv(A)
    end
    coords = neaten_coords(coords)
    coords = neaten_coords(coords)
    A = neaten_lattice(A)
    A = neaten_lattice(A)
    A = neaten_lattice(A)
    
    return A, coords, types
end

"""
    function parseQEinput(lines)

Parse a quantum espresso inputfile. Can only handle simple cases with explict CELL_PARAMETERS
Cannot handle nonzero ibrav. or celldm

Called by makecrys, doesn't need to be called directly.
"""
function parseQEinput(lines)

#    println("parseQE")

    posind = -1
    Aind = -1
    types = String[]
    A=zeros(3,3)
    coords = Float64[]
    nat = -1
    units = 1.0
    for line in lines
        liner = replace(replace(line, "," => " "), "=" => " ")

        sp = split(liner)
        
        if size(sp)[1] >= 1

            if posind >= 0
                posind += 1
                coords[posind,:] = parseARRfloat(sp[2:4])
                push!(types, sp[1])
                if posind >= nat
                    posind = -1
                end
            end

            if Aind >= 0
                Aind += 1
                A[Aind,:] = parseARRfloat(sp[1:3])
                if Aind >= 3
                    Aind = -1
                end
            end
            
            if sp[1] == "nat"
                nat = parse(Int64, sp[2])
                coords = zeros(nat,3)
                
            elseif  sp[1] == "ATOMIC_POSITIONS"
                posind = 0
                if length(sp) > 1
                    if sp[2] != "crystal" && sp[2] != "(crystal)"
                        println("warning, only crystal coords supported !!!!!!!!!!!!!!!!!!")
                    end
                end
            elseif sp[1] == "CELL_PARAMETERS"
                Aind = 0
                if length(sp) > 1
                    if sp[2] == "angstrom" || sp[2] == "Angstrom" || sp[2] == "(angstrom)"
                        units = Ang
                    else 
                        println("warning alat or other CELL_PARAMETERS not supported !!!!!!!!!!!!!!!!!!!!!!!")
                    end
                end
            end
                
                   
        end 
    end

#    println("A")
#    println(A)
    
    A = A / units

    coords = neaten_coords(coords)
    coords = neaten_coords(coords)
    A = neaten_lattice(A)
    A = neaten_lattice(A)
    A = neaten_lattice(A)
    

    return A, coords, types

    
end

"""
    function generate_supercell(crys, cell)

Generate supercell. cell is `[1,1,2]`, etc. 

Note, perfered notation is to use syntax:
`c * [1,1,2]` , where `c` is a crystal, thus 
using overloading of the `*` operator, 
rather than calling directly.
"""
function generate_supercell(c, cell)

    if typeof(cell) == Int64
        cell = [cell,cell,cell]
    end

    A_new = zeros(3,3)
    A_new[1,:] = c.A[1,:] *cell[1]
    A_new[2,:] = c.A[2,:] *cell[2]
    A_new[3,:] = c.A[3,:] *cell[3]

    n = prod(cell)
    coords_new = zeros(n*c.nat, 3)
    types_new = String[]

    counter = 0
    for x = 1:cell[1]
        for y = 1:cell[2]
            for z = 1:cell[3]
                for i = 1:c.nat
                    counter +=1
                    coords_new[counter,:] = (c.coords[i,:] + [x-1, y-1, z-1]) ./ [cell[1], cell[2], cell[3]]
                    push!(types_new, c.types[i])
                end
            end
        end
    end

    return makecrys(A_new, coords_new, types_new, units="Bohr")

#=    
    A = copy(crys.A)
    A[1,:] *= cell[1]
    A[2,:] *= cell[2]
    A[3,:] *= cell[3]

    cells = prod(cell)

#    println("cells ", cell, " " , cells)
    
    coords = zeros(crys.nat*cells, 3)

    types=String[]
    
    c = 0
    for i in 0.0:cell[1]-1
        for j in 0.0:cell[2]-1
            for k in 0.0:cell[3]-1
                for at in 1:crys.nat
                    c+= 1

                    coords_new[counter,:] = (c.coords[i,:] + [x-1, y-1, z-1]) ./ [a[1], a[2], a[3]]
                    coords[c,:] = (crys.coords[at,:] + [i j k]) ./ cell
                    push!(types, crys.types[at])
                end
            end
        end
    end

    csuper = makecrys(A, coords, types, units="Bohr")
    return csuper

    =#
    
end


"""
    function generate_optimum_supercell(c::crystal, dist)

Generate a supercell of crystal `c` where the periodic copies of all atoms
are at least `dist` apart. Will consider (some) linear combinations of the
initial lattice vectors and look for the "best" cell.  
"""
function generate_optimum_supercell(c::crystal, dist)

    #dist = convert_length(dist)
    if global_length_units == "Å"
        dist = dist / 0.529177210903
    end
    

    check_list1 = [1,-1,2,-2]
    check_list2 = [0,1,-1,2,-2, 3,-3, 4, -4]
    check_list3 = [0,1,-1,2,-2, 3,-3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8]
    check_list_small = [0,1,-1, 2, -2]

    cutoff_dist[(:Hx, :Hx)] = [dist, dist]

    #this function does the search. We first do a quicker search, then more thorough if the first one fails.
    good_list, A_new = gen(c, dist, check_list1, check_list_small)
    if length(good_list) == 0
        good_list, A_new = gen(c, dist, check_list2, check_list_small)
        if length(good_list) == 0
            good_list, A_new = gen(c, dist, check_list3, check_list_small)
            if length(good_list) == 0
                println("Sorry we did not find an appropriate cell.")
                return false
            end
        end
    end
    vol_factor = Int64(round(abs(det(A_new))/abs(det(c.A))))
    println("vol_factor ", vol_factor)

    coords_new = zeros(c.nat * vol_factor, 3)

    #generate new crystal
    found = 0
    types_new = []
    @time for x = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8]
        for y = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8]
            for z = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8]
                for n = 1:c.nat
                    coords = ((c.coords[n,:]+[x, y, z])' * c.A)*inv(A_new)
                    if coords[1] > -1e-10 && coords[1] < (1.0-1e-10) && coords[2] > -1e-10 && coords[2] < (1.0-1e-10) && coords[3] > -1e-10 && coords[3] < (1.0-1e-10)
                        found += 1
                        coords_new[found,:] = coords[:]
                        push!(types_new, c.stypes[n])
                        if found == c.nat * vol_factor
                            break
                        end
                    end
                end
            end
        end
    end

    
    
    cnew = makecrys(A_new, coords_new, types_new, units=:Bohr)

    
     
    return cnew
    
end

"""
    function gen(c, dist, check_list, check_list_small)

Helper function to search for optimum supercells 
"""
function gen(c, dist, check_list, check_list_small)
    
    println("dist $dist")
    good_list = []
    A_new = zeros(3,3)

    vol = zeros(length(check_list)^3 * length(check_list_small)^6)
    alength =  zeros(length(check_list)^3 * length(check_list_small)^6)
    counter = 0
    inds = zeros(Int64, length(check_list)^3 * length(check_list_small)^6, 9)
    dist_sq = dist^2
    
    @time for x1 in check_list
        for y1 in check_list_small
            for z1 in check_list_small
                A_new[1,:] = c.A[1,:] * x1 + c.A[2,:] * y1 + c.A[3,:] * z1
                if sum(A_new[1,:].^2) < dist_sq
                    continue
                end
                for x2 in check_list_small
                    for y2 in check_list
                        for z2 in check_list_small
                            A_new[2,:] = c.A[1,:] * x2 + c.A[2,:] * y2 + c.A[3,:] * z2
                            if sum(A_new[2,:].^2) < dist_sq
                                continue
                            end
                            for x3 in check_list_small
                                for y3 in check_list_small
                                    for z3 in check_list
#                                        println([x1 y1 z1 x2 y2 z2 x3 y3 z3])
                                        A_new[3,:] = c.A[1,:] * x3 + c.A[2,:] * y3 + c.A[3,:] * z3

                                        v = det(A_new)
                                        if v < 1e-3
                                            continue
                                        end
                                        if sum(A_new[3,:].^2) < dist_sq
                                            continue
                                        end
                                        
                                        counter +=1
                                        vol[counter] = (det(A_new))
                                        alength[counter] = sqrt(sum(A_new.^2))
                                        inds[counter, :] = [x1,x2,x3,y1,y2,y3,z1,z2,z3]
                                        
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end    

    
    bv = vol .> 1e-3 #we only consider right handed cells, rest are redundant
    vol = vol[bv]
    inds = inds[bv,:]
    alength = alength[bv]
    
    #find the smallest vol cells the meet requirements
    perm = sortperm(vol + 0.0000001 * alength) #we use heuristic that best candidates will have short total lattice vector lengths in addition to small volumes
    good_vec = []
    good_vol = 1000000000.0
    ctemp = makecrys([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], [0.0 0.0 0.0], [:Hx])
    println("perm")
    counter = 0

    tried = Dict()
    
    @time for p in perm

        
        key = Int64(round(vol[p]*10000000))
        if ! (key in keys(tried))
            tried[key] = 0
        end
#        println("key $key tried $(tried[key])")
        if tried[key] > 301
            continue
        end
            
        if vol[p] > good_vol + 1e-3
            break
        end
        x1,x2,x3,y1,y2,y3,z1,z2,z3 = inds[p,:]


        A_new[1,:] = c.A[1,:] * x1 + c.A[2,:] * y1 + c.A[3,:] * z1
        A_new[2,:] = c.A[1,:] * x2 + c.A[2,:] * y2 + c.A[3,:] * z2
        A_new[3,:] = c.A[1,:] * x3 + c.A[2,:] * y3 + c.A[3,:] * z3
        ctemp.A[:,:] = A_new[:,:]
        
        nz_inds, R_keep_ab_TT, nz_ind3, dist3_nonzero, dist_arr_TT, c_zero_tt, dmin_types_TT, dmin_types3_TT = distances_etc_3bdy_parallel_LV(ctemp, missing, 0.0, R=[2,2,2])

        if length(nz_inds) > 1
            tried[key] += 1
        end
        
        if size(nz_inds)[1] == 1
            nz_inds, R_keep_ab_TT, nz_ind3, dist3_nonzero, dist_arr_TT, c_zero_tt, dmin_types_TT, dmin_types3_TT = distances_etc_3bdy_parallel_LV(ctemp, missing, 0.0, R=[3,3,3])
            if size(nz_inds)[1] == 1
                nz_inds, R_keep_ab_TT, nz_ind3, dist3_nonzero, dist_arr_TT, c_zero_tt, dmin_types_TT, dmin_types3_TT = distances_etc_3bdy_parallel_LV(ctemp, missing, 0.0, R=[6,6,6])
                if size(nz_inds)[1] == 1
                    push!(good_list, p)
                    good_vol = vol[p]
                end
            end
        end
        if length(good_list) > 300 #enough candidates
            break
        end
    end



    if length(good_list) == 0
        return good_list, A_new
    end

    #now pick "best" lattice vectors from the small vol set.
    best_vec = 10000000000.0
    current_best = good_list[1]
    @time for p in good_list

#        vec_new = sum(inds[p,:].^2)

        x1,x2,x3,y1,y2,y3,z1,z2,z3 = inds[p,:]
        
        A_new[1,:] = c.A[1,:] * x1 + c.A[2,:] * y1 + c.A[3,:] * z1
        A_new[2,:] = c.A[1,:] * x2 + c.A[2,:] * y2 + c.A[3,:] * z2
        A_new[3,:] = c.A[1,:] * x3 + c.A[2,:] * y3 + c.A[3,:] * z3

        vec_new = sum(A_new[1,:].^2) + sum(A_new[2,:].^2) + sum(A_new[3,:].^2)
        
        if vec_new < best_vec
            current_best = p
            best_vec = vec_new
            
        end
    end

    x1,x2,x3,y1,y2,y3,z1,z2,z3 = inds[current_best,:]
    
    A_new[1,:] = c.A[1,:] * x1 + c.A[2,:] * y1 + c.A[3,:] * z1
    A_new[2,:] = c.A[1,:] * x2 + c.A[2,:] * y2 + c.A[3,:] * z2
    A_new[3,:] = c.A[1,:] * x3 + c.A[2,:] * y3 + c.A[3,:] * z3
    
    return good_list, A_new
        
end




"""
    function generate_random_distortion(crys, amag, strain_mag)

Randomly distort a crystal. `amag` is the atom distance, `strain_mag` is the strain magnitude
"""
function generate_random_distortion(crys, amag, strain_mag)

    amag = convert_length(amag)
    
    st = (rand(3,3) .- 0.5) * strain_mag
    st = (st + st')/2.0
    
    A = crys.A * (I + st)
    coords_real = (crys.coords * A) + (rand(crys.nat, 3) .- 0.5)*amag
    coords = coords_real * inv(A)

    return makecrys(A, coords, copy(crys.types), units="Bohr")
    
end

"""
    function write_poscar(crys, filename)

Write a `crystal` to a POSCAR.
"""
function write_poscar(crys, filename)

    fil = open(filename, "w")
    write(fil, "title kfg\n")
    write(fil, "1.0000000\n")

    Aang = crys.A * Ang
    
    write(fil, "$(Aang[1,1])  $(Aang[1,2])  $(Aang[1,3]) \n")
    write(fil, "$(Aang[2,1])  $(Aang[2,2])  $(Aang[2,3]) \n")
    write(fil, "$(Aang[3,1])  $(Aang[3,2])  $(Aang[3,3]) \n")

    tstr = ""
    numstr = ""

    ut = Set(crys.types)
    for t in ut
        tstr = tstr * t * " " 
        count = sum( t .== crys.types)
        numstr = numstr * "$count "
    end

#    tstr = ""
#    numstr = ""
 #   for t in crys.types
 #       tstr = tstr * t * " " 
 #       numstr = numstr * "1 "
 #   end


    write(fil,tstr*'\n')
    write(fil,numstr*'\n')
 
    write(fil,"Direct\n")
    for t in ut
        for at in 1:crys.nat
            if crys.types[at] == t
                write(fil, "$(crys.coords[at,1])  $(crys.coords[at,2]) $(crys.coords[at,3])\n")
            end
        end
    end
    close(fil)
    
    return 0
    
end

"""
    function write_xsf(c::crystal, filename="t.xsf", force=missing)

Write file for xcrysden. Forces optional
"""
function write_xsf(c::crystal; filename="t.xsf", force=missing)

    outstr = " INFO
nunit      1    1    1
unit   cell
celltype   primcell
shape   parapipedal
 END_INFO
 DIM-GROUP
           3           1
 PRIMVEC
"
    As = string.(c.A * Ang)
    for i = 1:3
        outstr *= As[i,1] * "   " *  As[i,2] * "   " * As[i,3]* "\n"
    end

    outstr *= " CONVVEC
"
    for i = 1:3
        outstr *= As[i,1] * "   " *  As[i,2] * "   " * As[i,3]* "\n"
    end

    outstr *= " PRIMCOORD
"                       
    nat = c.nat
    outstr *= "          $nat       1\n"

    pos = string.(c.coords * c.A * Ang)
    if !ismissing(force)
        fstring = string.(forces)
    else
        fstring = string.(zeros(size(c.coords)))
    end
    
    for i = 1:nat
#        println(atoms[c.types[i]].Z)
        outstr *= string(atoms[c.types[i]].Z)*"   "* pos[i,1] * "   " *  pos[i,2] * "   " * pos[i,3]*"   " * fstring[i,1]*"   " * fstring[i,2]*"   " * fstring[i,3]*"\n"
    end

    write(filename, outstr)
    
end

"""
    function write_axsf(CRYSTAL; filename="t.axsf", FORCES=missing)

write animated axsf file for xcrysden
Takes in array of crystal and optionally forces
"""
function write_axsf(CRYSTAL; filename="t.axsf", FORCES=missing)

    if typeof(CRYSTAL) == crystal
        write_xsf(CRYSTAL, filename=filename, force=FORCES)
    end
        
    nstep = length(CRYSTAL)
    outstr= "ANIMSTEPS  $nstep
CRYSTAL
"
    for i in 1:nstep
        c = CRYSTAL[i]
        nat = c.nat
        if !ismissing(FORCES)
            f = FORCES[i]
        else
            f = zeros(size(c.coords))
        end
        outstr *= "PRIMVEC $i\n"

        As = string.(c.A * Ang)
        for i = 1:3
            outstr *= As[i,1] * "   " *  As[i,2] * "   " * As[i,3]* "\n"
        end
        outstr *= "PRIMCOORD $i\n"
        outstr *= "$nat        1\n"
        pos = string.(c.coords * c.A * Ang)
        fstring = string.(f)
        
        for i = 1:nat
            outstr *= string(atoms[c.types[i]].Z)*"   "* pos[i,1] * "   " *  pos[i,2] * "   " * pos[i,3]*"   " * fstring[i,1]*"   " * fstring[i,2]*"   " * fstring[i,3]*"\n"
        end
    end
     
    write(filename, outstr)

end
    
"""
    function function write_efs(crys, energy, forces, stress, filename)

Write `crystal`, energy, force, stress to fake quantum espresso output file.
This is for testing purposes only.
"""
function write_efs(crys, energy, forces, stress, filename)

    fil = open(filename, "w")

    write(fil,"Program PWSCF fake\n")
    write(fil,"number of atoms/cell = $(crys.nat) \n")
#    ntypes = pos.shape[0]
    ntypes = 1
    write(fil,"number of types = 1\n")
    write(fil,"celldm(1)= 1.00\n")

    write(fil,"a(1) = ( $(crys.A[1,1]) $(crys.A[1,2]) $(crys.A[1,3]) ) \n")
    write(fil,"a(2) = ( $(crys.A[2,1]) $(crys.A[2,2]) $(crys.A[2,3]) ) \n")
    write(fil,"a(3) = ( $(crys.A[3,1]) $(crys.A[3,2]) $(crys.A[3,3]) ) \n")

    write(fil,"     site n.     atom                  positions (cryst. coord.)\n")
    for na in 1:crys.nat
        write(fil,"         1       $(crys.types[na])   tau(   $na) = (  $(crys.coords[na,1]) $(crys.coords[na,2]) $(crys.coords[na,3])  )\n")
    end
    write(fil,"\n")

    write(fil,"!    total energy              =     $energy Ry\n")

                   

    write(fil,"     Forces acting on atoms (Ry/au):\n")
    for na in 1:crys.nat
        write(fil,"     atom    $na type  1   force =   $(forces[na,1]) $(forces[na,2]) $(forces[na,3]) \n")
    end
        
    write(fil,"The non-local\n")

    write(fil,"          total   stress  (Ry/bohr**3)                   (kbar)     P=  ???\n")
    write(fil,"$(stress[1,1]) \t $(stress[1,2]) \t  $(stress[1,3])   0 0 0 \n")
    write(fil,"$(stress[2,1]) \t $(stress[2,2]) \t  $(stress[2,3])   0 0 0 \n")
    write(fil,"$(stress[3,1]) \t $(stress[3,2]) \t  $(stress[3,3])   0 0 0 \n")


    write(fil,"JOB DONE.\n")

    
    
    close(fil)

end

"""
    function get_grid(c, kden=55.0)

Get a default k-point grid size with `kden` density.
"""
function get_grid(c, kden=55.0)

    B = transpose(inv(c.A))
#    kden = 55.0 

#    b1 = 1.0/norm(c.A[1,:])
#    b2 = 1.0/norm(c.A[2,:])
#    b3 = 1.0/norm(c.A[3,:])

    b1 = norm(B[1,:])
    b2 = norm(B[2,:])
    b3 = norm(B[3,:])


    k1 = convert(Int, round(kden * b1))
    k2 = convert(Int, round(kden * b2))
    k3 = convert(Int, round(kden * b3))

    if k1%2 == 1
        k1 = k1 + 1
    end
    if k2%2 == 1
        k2 = k2 + 1
    end
    if k3%2 == 1
        k3 = k3 + 1
    end
    k1 = max(k1, 2)
    k2 = max(k2, 2)
    k3 = max(k3, 2)

    k1 = min(k1, 14)
    k2 = min(k2, 14)
    k3 = min(k3, 14)

#    k1 = min(k1, 30)
#    k2 = min(k2, 30)
#    k3 = min(k3, 30)

    
    kpoints = [k1, k2, k3]
#    println("get grid ", kpoints)
    return kpoints
    
#    R = [0,0,0]
#    for i = 1:3
#        R[i] = Int64(round(20.0/sum(c.A[1,:].^2)^0.5))
#        R[i] = max(R[i], 2)
#    end
#    println("get grid ", R)
#    return R

end    

"""
    function orbital_index(c::crystal)

Get correspondence between `crystal` and the TB orbital numbers.

`return ind2orb, orb2ind, etotal, nval`

- `ind2orb` dictionary which gives `[atom_number,atom_type, :orbital_symbol]` from TB index integer.
- `orb2ind` dictionary which gives TB index integer from `[atom_number,atom_type, :orbital_symbol]`
- `etotal` total DFT energy of atoms, for calculation atomization energy.
- `nval` number of valence orbitals.
"""
function orbital_index(c::crystal)

    ind2orb = Dict()
    orb2ind = Dict()

    atomtypes = []
    ntypes=0
    for t in c.stypes
        if !(t in atomtypes)
            ntypes += 1
            push!(atomtypes, t)
        end
    end

    wan_counter = 1

    nsemi = 0
    nval = 0
    nwan = 0
    ind = 0

    etotal = 0.0
    
#    for at in atomtypes
#        for (i, t) in enumerate(c.types)
#            if t == at
    for (i,t) in enumerate(c.stypes)
        atom = atoms[t]
        nsemi += atom.nsemicore
        nval += atom.nval
        nwan += atom.nwan

        etotal += atom.total_energy
        
        for o in atom.orbitals
            if o == :s
                ind += 1
                ind2orb[ind] = [i,t, :s]
            elseif o == :p
                ind += 1
                ind2orb[ind] = [i,t, :pz]
                ind += 1
                ind2orb[ind] = [i,t, :px]
                ind += 1
                ind2orb[ind] = [i,t, :py]
            elseif o == :d
                ind += 1
                ind2orb[ind] = [i,t, :dz2]
                ind += 1
                ind2orb[ind] = [i,t, :dxz]
                ind += 1
                ind2orb[ind] = [i,t, :dyz]
                ind += 1
                ind2orb[ind] = [i,t, :dx2_y2]
                ind += 1
                ind2orb[ind] = [i,t, :dxy]
            elseif o == :f
                ind += 1
                ind2orb[ind] = [i,t, :fz3]
                ind += 1
                ind2orb[ind] = [i,t, :fxz2]
                ind += 1
                ind2orb[ind] = [i,t, :fyz2]
                ind += 1
                ind2orb[ind] = [i,t, :fz_x2y2]
                ind += 1
                ind2orb[ind] = [i,t, :fxyz]
                ind += 1
                ind2orb[ind] = [i,t, :fx_x2_3y2]
                ind += 1
                ind2orb[ind] = [i,t, :fy_3x2_y2]
            else
                error("orbitals ????")
            end

        end
        if ind != round(nwan/2)
            error("counting error", ind, " " , nwan)
        end
        
    end
            
    nwan = 0
    for (i,t) in enumerate(c.types)
        atom = atoms[t]
        orb2ind[i] = 1+nwan : nwan + Int(round(atom.nwan/2))
        nwan += Int(round(atom.nwan/2))
    end

    #reverse dictionary
#    for ind = 1:convert(Int, round(nwan/2))
#        orb2ind[ind2orb[ind]] = ind
#        println(ind, " ", ind2orb[ind][1], " ", ind2orb[ind][2])
#        
#    end

    return ind2orb, orb2ind, etotal, nval
    
end


"""
    function interpolate(c1::crystal, c2::crystal; n = 5)

Linearly interpolate between two crystal structures. 

`return C`, a crystal array

- `c1, c2` start and ending crystal structures
- `n=5` number of interpolating structures, counting endpoints
"""
function interpolate(c1::crystal, c2::crystal; n = 5)

    if c1.nat != c2.nat
        println("error interplate nat $(c1.nat) != $(c2.nat)")
        return
    end
    if n < 3
        println("error n must be >= 3")
        return
    end
    C = []
    for x = 0.0:(1.0 / (n-1) ):1.0
        c = deepcopy(c1)
        c.coords = c1.coords * (1-x) + c2.coords * x
        c.A = c1.A * (1-x) + c2.A * x
        push!(C, c)
    end
    return C
    
    
end

"""
    function periodictable( print_atoms = missing)

Silly function. Default no argument display periodic table and
currently allowed atoms for calculations.

With array of strings as argument: print the atoms in the list. You
can also set argument to `false` to just print a periodic table for fun!
"""
function periodictable( print_atoms = missing)

    st = "Input atoms: "
    
    if ismissing(print_atoms)
        print_atoms = keys(atoms)
        st = "Currently supported atoms: "
    end
    
    st1 = "H                                                  He"
    st2 = "Li Be                               B  C  N  O  F  Ne"
    st3 = "Na Mg                               Al Si P  S  Cl Ar"
    st4 = "K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr"
    st5 = "Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe"
    st6 = "Cs Ba La Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn"

    lan = "         Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
    println()
    println("Full table: ")
    println("...........")
    println(st1)
    println(st2)
    println(st3)
    println(st4)
    println(st5)
    println(st6)
    println()
    println(lan)
    println()
    if print_atoms != false
        println(st)
        for a in keys(atom_radius)
            if typeof(a) != String
                continue
            end
            if !(a in print_atoms)
                
                if length(a) == 1
                    ar = a*" "
                else
                    ar = a
                end
                
                st1 = replace(st1, ar => "--")
                st2 = replace(st2, ar => "--")
                st3 = replace(st3, ar => "--")
                st4 = replace(st4, ar => "--")
                st5 = replace(st5, ar => "--")
                st6 = replace(st6, ar => "--")
                lan = replace(lan, ar => "--")
            end
        end
        println("...........")
        println(st1)
        println(st2)
        println(st3)
        println(st4)
        println(st5)
        println(st6)
        println()
        println(lan)
        println()
    end    
    
    
end



"""
    function concatenate(c1::crystal, c2::crystal)

Average unit cells, put atoms from each cell into average of the two cells
"""
function concatenate(c1::crystal, c2::crystal)
    A = (c1.A + c2.A)/2.0

    tt = vcat(c1.types, c2.types)
    coords = vcat(c1.coords, c2.coords)

    return makecrys(A,coords, tt, units="Bohr")
end


function neaten_coords(mat, tol=1e-4)
#    println("neaten")
    mat = mod.(mat , 1.0)
    for i in 1:size(mat)[1]
        for j in 1:size(mat)[2]
            for t in [0.0, 1/4,1/2,3/4,1/8,3/8,5/8,7/8,1/3,2/3,1/6,5/6]
                if abs(mat[i,j] - t) < tol
#                    println("$i $j mat[i,j]  $t")
                    mat[i,j] = t
                end
            end
        end
    end
    for i in 1:size(mat)[1]
        for j in 1:size(mat)[2]
            for ii in 1:size(mat)[1]
                for jj in 1:size(mat)[2]
                    if abs(mat[ii,jj]) < tol
                        mat[ii,jj] = 0.0
                    elseif abs(mat[ii,jj] - mat[i,j]) < tol
                        mat[ii,jj] = mat[i,j]
                    elseif abs( mat[ii,jj] - (1-mat[i,j]) ) < tol
                        mat[ii,jj] = (1-mat[i,j])
                    elseif abs( mat[ii,jj] - (0.5-mat[i,j]) ) < tol
                        mat[ii,jj] = (0.5-mat[i,j])
                    elseif abs( mat[ii,jj] - (0.5+mat[i,j]) ) < tol
                        mat[ii,jj] = (0.5+mat[i,j])
                    elseif abs( mat[ii,jj] - (0.5*mat[i,j]) ) < tol
                        mat[ii,jj] = (0.5*mat[i,j])
                    end
                end
            end
        end
    end
    
    return mat
end
function neaten_lattice(mat, tol=1e-4)

    
    for i in 1:size(mat)[1]
        for j in 1:size(mat)[2]
            for ii in 1:size(mat)[1]
                for jj in 1:size(mat)[2]
                    if abs(mat[ii,jj]) < tol
                        mat[ii,jj] = 0.0
                    else
                        for t in [1.0, -1.0, 0.5,-0.5,0.25,-0.25,1/3,-1/3,2/3, -2/3, sqrt(3)/2,-sqrt(3)/2, 1/sqrt(3),-1/sqrt(3),1/2/sqrt(3),-1/2/sqrt(3),sqrt(8/9),sqrt(8/9)]
                            if abs(mat[ii,jj] - t*mat[i,j]) < tol
                                mat[ii,jj] = mat[i,j] * t
                            end
                        end
                    end

                #=elseif abs(mat[ii,jj] - mat[i,j]) < tol
                        mat[ii,jj] = mat[i,j]
                    elseif abs( mat[ii,jj] - (1-mat[i,j]) ) < tol
                        mat[ii,jj] = (1-mat[i,j])
                    elseif abs( mat[ii,jj] - (0.5-mat[i,j]) ) < tol
                        mat[ii,jj] = (0.5-mat[i,j])
                    elseif abs( mat[ii,jj] - (0.5+mat[i,j]) ) < tol
                        mat[ii,jj] = (0.5+mat[i,j])
                    elseif abs( mat[ii,jj] - (0.5*mat[i,j]) ) < tol
                        mat[ii,jj] = (0.5*mat[i,j])
                    elseif abs( mat[ii,jj] - (sqrt(3)/2)*mat[i,j]) < tol
                        mat[ii,jj] = (sqrt(3)/2)*mat[i,j]
                    elseif abs( mat[ii,jj] - (1/sqrt(3))*mat[i,j])  < tol
                        mat[ii,jj] = (1/sqrt(3))*mat[i,j]
                    elseif abs( mat[ii,jj] - (1/2/sqrt(3))*mat[i,j]) < tol
                        mat[ii,jj] = (1/2/sqrt(3))*mat[i,j]
                    elseif abs( mat[ii,jj] - (1/sqrt(2))*mat[i,j])  < tol
                        mat[ii,jj] = (1/sqrt(2))*mat[i,j]
                    elseif abs( mat[ii,jj] - (1/3)*mat[i,j])  < tol
                        mat[ii,jj] = (1/3)*mat[i,j]
                    elseif abs( mat[ii,jj] - sqrt(8/9)*mat[i,j])  < tol
                        mat[ii,jj] = sqrt(8/9)*mat[i,j]
                    end
=#
#                    println("$i $j $ii $jj ", abs( mat[ii,jj] - (1/3)*mat[i,j]) )
                    
                end
            end
        end
    end
    
    return mat
end

"""
    function delete_atom(c::crystal, n)

Delete atom `n` from crystal `c`. Returns new crystal with 1 fewer atom.
"""
function delete_atom(c::crystal, n)

    if c.nat == 1
        println("cannot delete atom from 1 atom cell. Consider making a supercell first")
        return c
    end
    if n > c.nat
        println("warning, only $(c.nat) atoms in this crystal, you wanted to delete num: $n.")
    end
    
    bv = (1:c.nat) .!= n

    coords_new = c.coords[bv, :]
    types_new = c.types[bv]

    return makecrys(c.A, coords_new, types_new, units=:Bohr)
    
end

"""
    function add_atom(c::crystal, coord, type)

Add atom `coord` and `type` to crystal `c`. Returns new crystal with 1 more atom.
"""
function add_atom(c::crystal, coord, type)

    coords_new = zeros(eltype(c.coords), c.nat+1, 3)
    coords_new[end,:] = coord
    
    types_new = vcat(c.types, String(type))

    return makecrys(c.A, coords_new, types_new, units=:Bohr)
    
end


function neaten_crys(c::crystal)

    coords = neaten_coords(c.coords)
    coords = neaten_coords(coords)
    A = neaten_lattice(c.A)
    A = neaten_lattice(A)
    A = neaten_lattice(A)

    return makecrys(A, coords, c.types, units=:Bohr)

end

#=
function sphere(r, C;  n=100)   # r: radius; C: center [cx,cy,cz]
          
           u = range(-π, π; length = n)
           v = range(0, π; length = n)
           x = C[1] .+ r*cos.(u) * sin.(v)'
           y = C[2] .+ r*sin.(u) * sin.(v)'
           z = C[3] .+ r*ones(n) * cos.(v)'
           return x, y, z
end

function cylinder(p0, p1; r=0.25,  n=100)   # r: radius; C: center [cx,cy,cz]

    p0n = p0 / sqrt(sum(p0.^2))
    p1n = p0 / sqrt(sum(p0.^2))

    c = p0n * p1n'

    vcross = [0 -p0n[3] p0n[2]; p0n[3] 0 -p0n[1]; -p0n[2] p0n[1] 0] * b'

    R = I(3) + [0 -vcross[3] vcross[2]; vcross[3] 0 -vcross[1]; -vcross[2] vcross[1] 0] + v*v * 1 / (1+c)
    

    h = 1.0
    len = sqrt((p1 - p0).^2)
    u = range(0, 2pi, length=n)
    
    v1 = range(p0[1], p1[1], length=n)
    v2 = range(p0[2], p1[2], length=n)
    v3 = range(p0[3], p1[3], length=n)
    
    
    
    v = range(0, len, length=n)

    
    us = ones(n)*u'
    vs = v*ones(n)'
    #Surface parameterization
    x = r*cos.(us)
    y = r*sin.(us)
    z = vs

    
    
    return x, y, z
end
=#

function Plots.plot(c::crystal; bondlength=5.0)

    p = Plots.plot3d(legend=false, aspect_ratio=:equal)

    z = [0.0,0.0,0.]
    A1 = c.A[1,:]
    A2 = c.A[2,:]
    A3 = c.A[3,:]
    
    ap = [[z,A1], [z,A2], [z,A3], [A1, A1+A2], [A1,A1+A3], [A2, A2+A3], [A1+A2, A1+A2+A3], [A1+A3, A1+A2+A3], [A2+A3, A1+A2+A3], [A2, A1+A2], [A3, A1+A3], [A3, A2+A3]]

    mmin = min(0, minimum(c.A))
    mmax = max(0, maximum(c.A))
    
    for pts in ap
        Plots.plot3d!( [pts[1][1],pts[2][1]], [pts[1][2],pts[2][2]],[pts[1][3],pts[2][3]], color=:grey)
        mmin = min(mmin, pts[1][1])
        mmin = min(mmin, pts[1][2])
        mmin = min(mmin, pts[1][3])
        mmin = min(mmin, pts[2][1])
        mmin = min(mmin, pts[2][1])
        mmin = min(mmin, pts[2][1])

        mmax = max(mmax, pts[1][1])
        mmax = max(mmax, pts[1][2])
        mmax = max(mmax, pts[1][3])
        mmax = max(mmax, pts[2][1])
        mmax = max(mmax, pts[2][1])
        mmax = max(mmax, pts[2][1])
        
    end
    
#    plot3d!( [0, c.A[2,1]], [0, c.A[2,2]], [0, c.A[2,3]], color="grey")
#    plot3d!( [0, c.A[3,1]], [0, c.A[3,2]], [0, c.A[3,3]], color="grey")

#    plot3d!( [c.A[1,1], c.A[1,2]], [0, c.A[1,2]], [0, c.A[1,3]], color="grey")
#    plot3d!( [0, c.A[2,1]], [0, c.A[2,2]], [0, c.A[2,3]], color="grey")
#    plot3d!( [0, c.A[3,1]], [0, c.A[3,2]], [0, c.A[3,3]], color="grey")
    


    Plots.xlims!(mmin - 2.5, mmax+2.5)
    Plots.ylims!(mmin - 2.5, mmax+2.5)
    Plots.zlims!(mmin - 2.5, mmax+2.5)
    ccart = c.coords * c.A

    for bl = [bondlength, 5,5.5,6]
        count = 0
        
        for at1 in 1:c.nat
            for at2 in at1+1:c.nat
                if sqrt(sum((ccart[at1,:] - ccart[at2,:]).^2)) < bl
                    Plots.plot3d!([ccart[at1,1], ccart[at2,1]], [ccart[at1,2], ccart[at2,2]],[ccart[at1,3], ccart[at2,3]], color=:darkviolet, linewidth=1, alpha = 0.8)
                    count+=1
                end
            end
        end
        if count > 0
            break
        end
    end
    
    COLORS = [:blue, :orange, :green, :red, :black, :yellow, :cyan, :white]

    st = Set(c.types)

    for (ind, s) in enumerate(st)
        color = COLORS[mod1(ind, 8)]
        println("s $s $color")
        bitv  = (s .== c.types)
        Plots.scatter3d!(ccart[bitv,1], ccart[bitv,2], ccart[bitv,3], markersize=4, color=color, showaxis=false, ticks=false, grid=false, alpha = 0.7)
    end
    #    for at in 1:c.nat
        
        #        x,y,z = sphere(1.0, ccart[at,:];  n=n)
        #        surface!(x,y,z, grid=false, showaxis=false, ticks=false)
        #scatter!([x], [y], [z], markersize=2, edgecolor=false)
#        count += 1
#    end

    
    
    
    display(p)

end

function correct_charged_cell_alignment(A, atom, charge)

    return -charge_cell_missing_atom_correction[atom] / abs(det(A)) * charge
    
end

function correct_charged_cell_alignment(crys::crystal, atom, charge)

    return correct_charged_cell_alignment(crys.A, atom, charge)
    
end

function correct_charged_cell_alignment(crys::crystal, atom::Array, charge)

    tot = 0.0
    for a in atom
        tot +=  correct_charged_cell_alignment(crys.A, a, charge)
    end

    return tot
    
end

function correct_charged_cell_alignment(A, atom::Array, charge)

    tot = 0.0
    for a in atom
        tot +=  correct_charged_cell_alignment(A, a, charge)
    end

    return tot
    
end


function generate_defect_structure(c_start ; defect_type=:vacancy, defect_atom = missing, sub_atom= missing, distortion_mag = 1.5, number_of_rand = 3, rand_dist = 5.1, center_defect=false)


    c_start = deepcopy(c_start)
    c_start.coords = mod.(c_start.coords , 1.0)

    defect_type = Symbol(defect_type)
    if !ismissing(defect_atom)
        defect_atom = Symbol(defect_atom)
    end
    if !ismissing(sub_atom)
        sub_atom = Symbol(sub_atom)
    end
    
    
    if defect_type == :vacancy || defect_type == :vac

        if ismissing(defect_atom )
            println("we use first atom as defect type, otherwise specify defect_atom : ", c_start.stypes)
            defect_atom = c_start.stypes[1]
        end

        
        stypes = deepcopy(c_start.stypes)
        #old version ind = findfirst(c_start.stypes .== Symbol(defect_atom))
        
        #search for defect location closest to 0.5 0.5 0.5 in fractional coords
        ind = 0
        best_dist = 10000000000000000.0
        for n = 1:c_start.nat
            if c_start.stypes[n] == Symbol(defect_atom)
                if center_defect == true
                    t = sum((c_start.coords[n,:] - [0.5,0.5,0.5]).^2)
                else
                    t = sum((c_start.coords[n,:] - [0.0,0.0,0.0]).^2)
                end                    
                if t < best_dist
                    best_dist = t
                    ind = n
                end
            end
        end

        if ind == 0
            println("warning, we could not find the atom you wanted to remove: $defect_atom")
            return
        end
        defect_location = c_start.coords[[ind],:]*c_start.A

        

        keep = (1:c_start.nat) .!= ind
        c_pristine = makecrys(c_start.A, c_start.coords[keep,:], c_start.stypes[keep], units=:Bohr)
        if distortion_mag < 1e-10
            C = [deepcopy(c_pristine)]
        else
            R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel_LV(c_start,6.0, 0.0,var_type=Float64, return_floats=true)
            nkeep_ab = size(R_keep_ab)[1]
            todistort = Set()
            for c = 1:nkeep_ab
                a1 = R_keep_ab[c,2]
                a2 = R_keep_ab[c,3]
                dist_a = dist_arr[c,1]
                if dist_a <= rand_dist
                    if a1 == ind 
                        todistort = push!(todistort, a2)
                    end
                    if a2 == ind
                        todistort = push!(todistort, a1)
                    end
                end
            end

            C = []
            for counter = 1:number_of_rand

                therand = rand(length(todistort)-1, 3)
                therand = therand - repeat(sum(therand, dims=1) / (length(todistort)-1), length(todistort)-1)
#                println("therand")
#                println(therand)

                
                coords = deepcopy(c_start.coords)
                x = 0
                for t in todistort
                    if keep[t] == true
                        x+=1
                        coords[t,:] = (coords[[t],:] * c_start.A + therand[[x],:]  * distortion_mag) * inv(c_start.A)
#                        stypes[t] = :Hg #for testing

                    end
                end

                
#                coords = deepcopy(c_start.coords)
#                for t in todistort
#                    coords[t,:] = (coords[[t],:] * c_start.A + (rand(1,3) .- 0.5) * distortion_mag) * inv(c_start.A)
#                end

 #               println("sum start ", sum(c_start.coords[keep,:], dims=1))
 #               coords = coords + repeat( (sum(c_start.coords[keep,:], dims=1)  - sum(coords[keep,:], dims=1)) / (c_start.nat-1), c_start.nat)
 #               println("sum coords ", sum(coords[keep,:], dims=1))
                
                c_temp = makecrys(c_start.A, coords[keep,:], stypes[keep],units=:Bohr)
                push!(C, c_temp)
                
            end
        end

    elseif defect_type == :substitution || defect_type == :sub

        if ismissing(defect_atom )
            println("we use first atom as defect type, otherwise specify defect_atom : ", c_start.stypes[1])
            defect_atom = c_start.stypes[1]
        end
        if ismissing(sub_atom )
            x = findfirst(c_start.stypes .!= defect_atom)
            if isnothing(x)
                println("this calculation makes no sense, you shouldn't substitue an atom with the same atom and expect anything different")
                x = 1
            end
            sub_atom = c_start.stypes[x]
            println("we use second atom as substituation type, otherwise specify sub_atom : ", sub_atom)
        end
        
        
        #ind = findfirst(c_start.stypes .== Symbol(defect_atom))
        ind = 0
        best_dist = 10000000000000000.0
        for n = 1:c_start.nat
            if c_start.stypes[n] == Symbol(defect_atom)
#                t = sum((c_start.coords[n,:] - [0.5,0.5,0.5]).^2)
                if center_defect == true
                    t = sum((c_start.coords[n,:] - [0.5,0.5,0.5]).^2)
                else
                    t = sum((c_start.coords[n,:] - [0.0,0.0,0.0]).^2)
                end                    

                if t < best_dist
                    best_dist = t
                    ind = n
                end
            end
        end
        if ind == 0
            println("warning, we could not find the atom you wanted to remove: $defect_atom")
            return
        end
        defect_location = c_start.coords[[ind],:]*c_start.A

        
        stypes = deepcopy(c_start.stypes)
        stypes[ind] = sub_atom

            
        c_pristine = makecrys(c_start.A, c_start.coords[:,:], stypes,units=:Bohr)
        if distortion_mag < 1e-10
            C = [deepcopy(c_pristine)]
        else
            R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel_LV(c_start,6.0, 0.0,var_type=Float64, return_floats=true)
            nkeep_ab = size(R_keep_ab)[1]
            todistort = Set()
            for c = 1:nkeep_ab
                a1 = R_keep_ab[c,2]
                a2 = R_keep_ab[c,3]
                dist_a = dist_arr[c,1]
                if dist_a <= rand_dist
                    if a1 == ind 
                        todistort = push!(todistort, a2)
                    end
                    if a2 == ind
                        todistort = push!(todistort, a1)
                    end
                end
            end

            C = []
            for counter = 1:number_of_rand

                therand = rand(length(todistort), 3)
                therand = therand - repeat(sum(therand, dims=1) / (length(todistort)), length(todistort))
#                println("sum therand , ", sum(therand, dims=1))

                
                coords = deepcopy(c_start.coords)
                for (i,t) in enumerate(todistort)
                    coords[t,:] = (coords[[t],:] * c_start.A + (therand[i,:]' ) * distortion_mag) * inv(c_start.A)

                    #                    stypes[t] = :Hg #for testing

                end
                stypes[ind] = sub_atom


                
                c_temp = makecrys(c_start.A, coords[:,:], stypes,units=:Bohr)
                push!(C, c_temp)
                
            end
        end

        
    else

        println("I didn't recognize  defect_type $defect_type, please try :vac or :sub")
        return c_start, [c_start], [0 0 0]
        
    end

    
    return c_pristine, C,         defect_location 

end

function get_dist(cart, c)
    dist = ones(c.nat) * 1000000000000.0
    c_cart = c.coords*c.A
    for at = 1:c.nat
        for x = -2:2
            for y = -2:2
                for z = -2:2
                    d = sum( ( c_cart[[at],:] + ([x y z] * c.A) - cart).^2)

                    if d < dist[at]
                        dist[at] = d
                    end
                end
            end
        end
    end
    return sqrt.(dist)
end

function align_crystal(c1, c2)

    if sum(abs.(c1.A - c2.A)) > 1e-6
        println("WARNING, trying align_crystal with different unit cells probably will not work")
    end
    
    if c1.nat > c2.nat
        ct = deepcopy(c1)
        c1 = deepcopy(c2)
        c2 = ct
    end

    c1cart = c1.coords * c1.A
    c2cart = c2.coords * c2.A
    
    best_dist = ones(c2.nat) * 1000000000.0
    best_match = zeros(Int64, c2.nat)
    for j in 1:c2.nat

        for at = 1:c1.nat
            for x = -2:2
                for y = -2:2
                    for z = -2:2
                        d = sum( (c2cart[j,:] + ([x y z] * c2.A)[:] - c1cart[at,:]).^2)
                        if d < best_dist[j]
                            best_match[j] = at
                            best_dist[j] = d
                        end
                    end
                end
            end
        end
    end    

    return sqrt.(best_dist), best_match

end



function merge_crystal(c_small, c_big, vacancy_location = missing)

    A = deepcopy(c_big.A)
    coords = deepcopy(c_big.coords)
    cart = coords  * A
    stypes = deepcopy(c_big.stypes)


    small_cart = c_small.coords * c_small.A
    
    if !ismissing(vacancy_location)

        todel = findmin(sum( (cart - repeat(vacancy_location, c_big.nat)).^2, dims=2))
        println("todel ", todel)
        coords = coords[ (1:c_big.nat) .!= todel, :]
        stypes = stypes[ (1:c_big.nat) .!= todel]

    end

    for j in 1:c_small.nat

        best_dist = 1000000000.0
        best_match = 0
        for at = 1:c_big.nat
            for x = -3:3
                for y = -3:3
                    for z = -3:3
                        d = sum( (cart[at,:] + ([x y z] * A)[:] - small_cart[j,:]).^2)
                        if d < best_dist
                            best_match = at
                            best_dist = d
                        end
                    end
                end
            end
        end
#        match = sortperm(sum( (cart - repeat(small_cart[[j],:], c_big.nat)).^2, dims=2)[:])[1]
#        println("match ", match  ,  " " , cart[[match], :], " ", small_cart[[j],:], " dist ", sum( (cart[[match], :] - small_cart[[j],:]).^2 )   )
        coords[best_match, :] = mod.(small_cart[[j],:] * inv(A), 1.0)
        stypes[best_match]  = c_small.stypes[j]
    end

    return makecrys(A, coords, stypes,units=:Bohr)

    
end


end #end module

#using .CrystalMod
#Base.:(==)(x::crystal, y::crystal) = (x.A == y.A && x.coords == y.coords && x.types == y.types && x.nat==y.nat)

"""
    module AtomicMod

Simple module for atom type.
"""
module AtomicMod

using ..Utility:str_w_spaces
using ..BandTools:band_energy
using ..BandTools:smearing_energy

"""
    struct atom

Hold basic atomic information

- `name::String` 
- `Z::Int64` Atomic number
- `row::Float64` Periodic table row
- `col::Float64` Periodic table col
- `mass::Float64` Mass in amu 
- `nval::Float64` Number of valence electrons in TB calculation.
- `nsemicore::Int64` Number of semicore electrons. Depends on pseudopotential choice
- `nwan::Int64` Number of TB orbitals
- `orbitals::Array{Symbol,1}` Names of the orbitals, like `[:s, :p]`
- `total_energy::Float64` DFT total energy, depends on pseudopotentials.
- `eigs::Dict` orbital eigenvalues of isolated non-spin-polarized nuetral atom.
- `energy_offset::Float64` energy to add to TB calculation to make isolated atoms have zero energy.
- `U::Float64` U value for Ewald correction

"""
struct atom
    name::String
    Z::Int64
    row::Float64
    col::Float64    
    mass::Float64
    nval::Float64
    nsemicore::Int64
    nwan::Int64
    orbitals::Array{Symbol,1}
    total_energy::Float64
    eigs::Dict
    energy_offset::Float64
    U::Float64
end




#printing
Base.show(io::IO, a::atom) = begin

    println(io,str_w_spaces([a.name, "Z=", a.Z, "Mass=",a.mass, "RC=",a.row, a.col]))
    println(io,str_w_spaces(["Nval=", a.nval, "Nsemicore=", a.nsemicore, "Orbitals=",a.orbitals]))
    println(io,"U = ",a.U)
    println(io,"total energy=", a.total_energy)
    for i in a.orbitals
        println(i, "  " , a.eigs[i])
    end
end




"""
    function makeatom(name, Z, row, col, mass, nval, nsemicore, orbitals, etot, eigs, vac_potential=0.0, U=0.0)

Constructor for atom.
"""
function makeatom(name, Z, row, col, mass, nval, nsemicore, orbitals, etot, eigs, vac_potential=0.0, U=0.0)
    #vac potential in ryd

    orb2 = map(x->convert(Symbol, x), orbitals)

    nwan = 0

    for o in orb2
        if o == :s
            nwan += 2
        elseif o == :p
            nwan += 6
        elseif o == :d
            nwan += 10
        elseif o == :f
            nwan += 14
        else
            exit("oribtal entry incorrect, good values are s p d f: ", orbitals, orb2)
        end
    end

    convert_ev_ryd = 1.0/13.605693122
    EIGS = zeros(Int64(nwan/2), 1)
    c=0
    for (i,o) in enumerate(orb2)
        if o == :s
            c+=1
            EIGS[c,1] = eigs[i]*convert_ev_ryd
        elseif o == :p
            c += 3
            EIGS[c-2,1] = eigs[i]*convert_ev_ryd
            EIGS[c-1,1] = eigs[i]*convert_ev_ryd
            EIGS[c  ,1] = eigs[i]*convert_ev_ryd            
        elseif o == :d
            c += 5
#            EIGS[c-5,1] = eigs[i]*convert_ev_ryd
            EIGS[c-4,1] = eigs[i]*convert_ev_ryd
            EIGS[c-3,1] = eigs[i]*convert_ev_ryd
            EIGS[c-2,1] = eigs[i]*convert_ev_ryd
            EIGS[c-1,1] = eigs[i]*convert_ev_ryd
            EIGS[c  ,1] = eigs[i]*convert_ev_ryd            
        elseif o == :f
            c += 7
#            EIGS[c-7,1] = eigs[i]*convert_ev_ryd
            EIGS[c-6,1] = eigs[i]*convert_ev_ryd
            EIGS[c-5,1] = eigs[i]*convert_ev_ryd
            EIGS[c-4,1] = eigs[i]*convert_ev_ryd
            EIGS[c-3,1] = eigs[i]*convert_ev_ryd
            EIGS[c-2,1] = eigs[i]*convert_ev_ryd
            EIGS[c-1,1] = eigs[i]*convert_ev_ryd
            EIGS[c  ,1] = eigs[i]*convert_ev_ryd            
        end
    end
    

    EIGS = EIGS .- vac_potential

#    println("EIGS after subtraction")
#    println(EIGS)

    band_energy0, efermi = band_energy(EIGS, [1.0], nval, returnef=true)
    smear_energy = smearing_energy(EIGS, [1.0], efermi)

#    if name=="Ta"
#        println("name $name")
#        println("EIGS before subtraction")
#        println(EIGS)
#        println("band_energy0 ", band_energy0, " ", efermi)
#    end




#    shift = (0.0 - band_energy0)/nval

    energy_offset = -(band_energy0+smear_energy)  #this is the atomic energy, no spin, after vac correction

    if name == "X" || name == "Xa"
        energy_offset = 0.0
    end

#    println("$name energy_offset $energy_offset    $band_energy0 $smear_energy  efermi $efermi")

    d = Dict()

    for (o, i) in zip(orb2, eigs)
#        d[o] = (i * convert_ev_ryd + shift)
        d[o] = (i * convert_ev_ryd - vac_potential)

#        println("d ", d[o], " ", o)
    end

#    band_energy_new = band_energy(EIGS .- band_energy0 / nval , [1.0], nval)#
#
#    println("Loading $name, new band energy $band_energy_new, ", d[:s])
    
    
    return atom(name, Z, row, col, mass, nval, nsemicore, nwan, orb2, etot, d, energy_offset, U)
    
end
    
end #end module




"""
    module AtomicMod

Simple module for atom type.
"""
module AtomicMod

using ..Utility:str_w_spaces
using ..BandTools:band_energy
using ..BandTools:smearing_energy
using ..Utility:parse_str_ARR_float

struct mag
    
    name::String
    Z::Int64
    Wss::Float64
    Wsp::Float64
    Wsd::Float64
    Wpp::Float64
    Wpd::Float64
    Wdd::Float64
    WpD::Float64
    WdD::Float64
    Xss::Float64
    Xsp::Float64
    Xsd::Float64
    Xpp::Float64
    Xpd::Float64
    Xdd::Float64
    XpD::Float64
    XdD::Float64
    W::Array{Float64,2}
    X::Array{Float64,2}
    
end

function make_mag(str)
    
    sp = split(str)
    num = parse(Int64, sp[1])
    name = sp[2]

    mag_arr = parse_str_ARR_float(sp[3:end])

    W = zeros(9,9)
    X = zeros(9,9)

    Wss = mag_arr[1]
    Wsp = mag_arr[2]
    Wsd = mag_arr[3]
    Wpp = mag_arr[4]
    Wpd = mag_arr[5]
    Wdd = mag_arr[6]
    WpD = mag_arr[7]
    WdD = mag_arr[8]
    Xss = mag_arr[9]
    Xsp = mag_arr[10]
    Xsd = mag_arr[11]    
    Xpp = mag_arr[12]    
    Xpd = mag_arr[13]    
    Xdd = mag_arr[14]    
    XpD = mag_arr[15]    
    XdD = mag_arr[16]    

    inds = [:s, :d, :d, :d, :d, :d, :p, :p, :p]
    for i = 1:9
        ind1 = inds[i]
        for j = 1:9
            ind2 = inds[j]            
            
            if ind1 == :s && ind2 == :s
                w = Wss
                x = Xss
            elseif ind1 == :s && ind2 == :p
                w = Wsp
                x = Xsp
            elseif ind1 == :s && ind2 == :d
                w = Wsd
                x = Xsd
            elseif ind1 == :p && ind2 == :s
                w = Wsp
                x = Xsp
            elseif ind1 == :d && ind2 == :s
                w = Wsd
                x = Xsd
            elseif ind1 == :p && ind2 == :d
                w = Wpd
                x = Xpd
            elseif ind1 == :d && ind2 == :p
                w = Wpd
                x = Xpd
            elseif ind1 == :d && ind2 == :d
                if ind1 == ind2
                    w = WdD
                    x = XdD
                else
                    w = Wdd
                    x = Xdd
                end
            elseif ind1 == :p && ind2 == :p
                if ind1 == ind2
                    w = WpD
                    x = XpD
                else
                    w = Wpp
                    x = Xpp
                end
            end
            W[i,j] = w
            X[i,j] = x

        end
    end


    return mag(name, num, mag_arr[1], mag_arr[2], mag_arr[3], mag_arr[4], mag_arr[5], mag_arr[6], mag_arr[7], mag_arr[8], mag_arr[9], mag_arr[10], mag_arr[11], mag_arr[12], mag_arr[13], mag_arr[14], mag_arr[15], mag_arr[16], W, X)

end


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
mutable struct atom
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
    U3::Float64
    efermi::Float64
    fullU::Bool
    Us::Float64
    Usp::Float64
    Up::Float64
    Upp::Float64
    Ud::Float64
    Usd::Float64
    Upd::Float64
    Udd::Float64
    Umat::Array{Float64,2}
end




#printing
Base.show(io::IO, a::atom) = begin

    println(io,str_w_spaces([a.name, "Z=", a.Z, "Mass=",a.mass, "RC=",a.row, a.col]))
    println(io,str_w_spaces(["Nval=", a.nval, "Nsemicore=", a.nsemicore, "Orbitals=",a.orbitals]))
    println(io,"U = ",a.U)
    println(io,"total energy=", a.total_energy)
    println(io,"efermi (isolated netural)=", a.efermi)
    for i in a.orbitals
        println(i, "  " , a.eigs[i])
    end
end




"""
    function makeatom(name, Z, row, col, mass, nval, nsemicore, orbitals, etot, eigs, vac_potential=0.0, U=0.0)

Constructor for atom.
"""
function makeatom(name, Z, row, col, mass, nval, nsemicore, orbitals, etot, eigs, vac_potential=0.0, U=0.0,U3=0.0, fullU=false, Uarr=zeros(8))
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
        if o == :s
            d[1] = (i * convert_ev_ryd - vac_potential)
        elseif o == :p
            d[2] = (i * convert_ev_ryd - vac_potential)
        elseif o == :d
            d[3] = (i * convert_ev_ryd - vac_potential)
        elseif o == :f
            d[4] = (i * convert_ev_ryd - vac_potential)
        end
        
#        println("d ", d[o], " ", o)
    end

#    band_energy_new = band_energy(EIGS .- band_energy0 / nval , [1.0], nval)#
#
#    println("Loading $name, new band energy $band_energy_new, ", d[:s])

    
    if fullU == false
        Umat = ones(nwan,nwan) * U
    else
        Umat = zeros(nwan,nwan)
        if orb2 == [:s]
            orblist = [:s]
        elseif orb2 == [:s, :p]
            orblist = [:s, :p, :p, :p]
        elseif orb2 == [:s, :d, :p]
            orblist = [:s, :d, :d, :d, :d, :d, :p, :p, :p]
        elseif orb2 == [:s, :p, :d]
            orblist = [:s, :p, :p, :p, :d, :d, :d, :d, :d]
        end
        for (c1,o1) in enumerate(orblist)
            for (c2,o2) in enumerate(orblist)
                if o1 == :s && o2 == :s
                    Umat[c1,c2] = Uarr[1]
                elseif (o1 == :s && o2 == :p) || (o1 == :p && o2 == :s)
                    Umat[c1,c2] = Uarr[2]
                elseif (o1 == :p && o2 == :p) && c1 == c2
                    Umat[c1,c2] = Uarr[3]
                elseif (o1 == :p && o2 == :p) && c1 != c2
                    Umat[c1,c2] = Uarr[4]
                elseif (o1 == :d && o2 == :d) && c1 == c2
                    Umat[c1,c2] = Uarr[5]
                elseif (o1 == :s && o2 == :d) || (o1 == :d && o2 == :s) 
                    Umat[c1,c2] = Uarr[6]
                elseif (o1 == :p && o2 == :d) || (o1 == :d && o2 == :p) 
                    Umat[c1,c2] = Uarr[7]
                elseif (o1 == :d && o2 == :d) && c1 != c2
                    Umat[c1,c2] = Uarr[8]
                end
            end
        end
    end

    return atom(name, Z, row, col, mass, nval, nsemicore, nwan, orb2, etot, d, energy_offset, U, U3,efermi, fullU, Uarr[1], Uarr[2], Uarr[3], Uarr[4], Uarr[5], Uarr[6], Uarr[7], Uarr[8], Umat)
    
end
    
end #end module




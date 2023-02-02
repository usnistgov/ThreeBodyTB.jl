
#no longer a stand alone module

#module Magnetic

#using ..TB:tb_crys
#using ..CrystalMod:crystal
#using ..CrystalMod:orbital_index
#using ..TB:summarize_orb
#using ..Atomdata:atoms
using ..AtomicMag:magnetic

function get_spin_h1(tbc::tb_crys)

     return get_spin_h1(tbc, tbc.eden)

end

function get_spin_h1(tbcK::tb_crys_kspace)

     return get_spin_h1(tbcK.crys, tbcK.eden)

end

function get_spin_h1(tbcK::tb_crys_kspace, eden)

     return get_spin_h1(tbcK.crys, eden)

end

function magnetic_energy(tbc::tb_crys, eden)
    return magnetic_energy(tbc.crys,eden)
end

function magnetic_energy(tbc::tb_crys)
    return magnetic_energy(tbc.crys,tbc.eden)
end

function magnetic_energy(tbcK::tb_crys_kspace)
    return magnetic_energy(tbcK.crys,tbcK.eden)
end

function magnetic_energy(tbcK::tb_crys_kspace, eden)
    return magnetic_energy(tbcK.crys, eden)
end


function magnetic_energy(crys::crystal, eden::Array{Float64,2})
    if size(eden)[1] == 1
        return 0.0
    else
        Wm, Xm, energy = WX(crys, eden[1,:] - eden[2,:])
#        println("magnetic energy $energy")
        return energy 
    end
end

function WX(tbc::tb_crys, spin)
    return WX(tbc.crys, spin)
end

function WX(tbc::tb_crys_kspace, spin)
    return WX(tbc.crys, spin)
end

function WX(crys::crystal, spin)

    c = 0
    spin = spin[:]

    nwan = size(spin)[1]

    Wm = zeros(nwan)
    Xm = zeros(nwan)
    energy = 0.0
    for (at,t) in enumerate(crys.types)
        atom = atoms[t]
        mag = magnetic[t]
        nwan = Int64(atom.nwan/2)
        spin_at = spin[c .+ (1:nwan)]
        Wm[c.+(1:nwan)] = mag.W * spin_at
        Xm[c.+(1:nwan)] = mag.X * abs.(spin_at)
        energy += 0.5*(sum(spin_at .* Wm[c.+(1:nwan)]) + sum(abs.(spin_at) .* Xm[c.+(1:nwan)]))
        c += nwan

#        println("W")
#        println(mag.W)
#        println("X")
#        println(mag.X)

    end

#    println(Wm)
#    println(Xm)
    
    return Wm, Xm, energy

end


function get_spin_h1(tbc::tb_crys, chargeden::Array{Float64,2})
    
    return get_spin_h1(tbc.crys, chargeden::Array{Float64,2})

end

function get_spin_h1(crys::crystal, chargeden::Array{Float64,2})
    

    spin = chargeden[1,:] - chargeden[2,:]

#    h1up = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)
#    h1dn = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)
    
    ind2orb, orb2ind, etotal, nval = orbital_index(crys)
    nwan = maximum(keys(ind2orb))

    h1spin = zeros(Complex{Float64}, 2, nwan, nwan)


    Wm, Xm, energy = WX(crys, spin)
    

    for i1 = 1:nwan
        a1,t1,orb1 = ind2orb[i1]
        sorb1 = summarize_orb(orb1)

        for i2 = 1:nwan

            a2,t2,orb2 = ind2orb[i2]
            sorb2 = summarize_orb(orb2)
            
#            h1up[i1,i2] = 0.5*Wm[i1] + 0.5*Xm[i1] + 0.5*Wm[i2] + 0.5*Xm[i2]
#            h1dn[i1,i2] = -0.5*Wm[i1] + 0.5*Xm[i1] + -0.5*Wm[i2] + 0.5*Xm[i2]

            h1spin[1,i1,i2] = 0.5*Wm[i1] + 0.5*Xm[i1] + 0.5*Wm[i2] + 0.5*Xm[i2]
            h1spin[2,i1,i2] = -0.5*Wm[i1] + 0.5*Xm[i1] + -0.5*Wm[i2] + 0.5*Xm[i2]

        end
    end

#    return h1up, h1dn
    return h1spin

end

function get_magmom(tbc::tb_crys)
    if size(tbc.eden)[1] == 1
        return zeros(tbc.crys.nat)
    else
        return get_magmom(tbc.crys, tbc.eden)
    end
end

function get_magmom(tbc::tb_crys, eden)
    return get_magmom(tbc.crys, eden)    
end

function get_magmom(c::crystal, eden)

    spin = eden[1,:] - eden[2,:]
    ind = 0
    mm = zeros(c.nat)
    for (counter, t) in enumerate(c.types)
        atom = atoms[t]
        mag = magnetic[t]
        nwan = Int64(atom.nwan/2)
        mm[counter] = sum(spin[ind .+ (1:nwan)])
        ind += nwan
    end
    return mm

end



#end #end module

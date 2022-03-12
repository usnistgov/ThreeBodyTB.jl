
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

function magnetic_energy(tbc::tb_crys)

    return magnetic_energy(tbc,tbc.eden)
end

function magnetic_energy(tbc::tb_crys, eden::Array{Float64,2})
    if size(eden)[1] == 1
        return 0.0
    else
        Wm, Xm, energy = WX(tbc, eden[1,:] - eden[2,:])
        return energy
    end
end

function WX(tbc::tb_crys, spin)

    c = 0
    spin = spin[:]

    Wm = zeros(tbc.tb.nwan)
    Xm = zeros(tbc.tb.nwan)
    energy = 0.0
    for (at,t) in enumerate(tbc.crys.types)
        atom = atoms[t]
        mag = magnetic[t]
        nwan = Int64(atom.nwan/2)
#        println("c $c at $at t $t nwan $nwan ", c .+ (1:nwan))

        spin_at = spin[c .+ (1:nwan)]

#        println("size W ", size(mag.W), " spin_at ", size(spin_at))
        Wm[c.+(1:nwan)] = mag.W * spin_at
        Xm[c.+(1:nwan)] = mag.X * abs.(spin_at)

#        println("Size Wm ", size(Wm), " size Wm2 ", size(Wm[c.+(1:nwan)]), " spin_at ", size(spin_at)

        energy += sum(spin_at .* Wm[c.+(1:nwan)]) + sum(abs.(spin_at) .* Xm[c.+(1:nwan)])
        c += nwan
    end
    
    return Wm, Xm, energy

end



function get_spin_h1(tbc::tb_crys, chargeden::Array{Float64,2})


    spin = chargeden[1,:] - chargeden[2,:]

    h1up = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)
    h1dn = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)

    Wm, Xm, energy = WX(tbc, spin)
    
    ind2orb, orb2ind, etotal, nval = orbital_index(tbc.crys)

    for i1 = 1:tbc.tb.nwan
        a1,t1,orb1 = ind2orb[i1]
        sorb1 = summarize_orb(orb1)

        for i2 = 1:tbc.tb.nwan

            a2,t2,orb2 = ind2orb[i2]
            sorb2 = summarize_orb(orb2)
            
            h1up[i1,i2] = 0.5*Wm[i1] + 0.5*Xm[i1] + 0.5*Wm[i2] + 0.5*Xm[i2]
            h1dn[i1,i2] = -0.5*Wm[i1] + 0.5*Xm[i1] + -0.5*Wm[i2] + 0.5*Xm[i2]

        end
    end

    return h1up, h1dn

end






#end #end module

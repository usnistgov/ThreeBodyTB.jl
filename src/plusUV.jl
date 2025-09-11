
module UV

using LinearAlgebra
using ..TB:summarize_orb

function addU(denmat, Umat, nwan, ind2orb, orb2ind, crys, nspin)

    potU = zeros(Float64, nwan, nwan,nspin)
    iden = I(nwan)
    energy = 0.0
    for spin = 1:nspin
        for a1 = 1:nwan
            at1,t1,o1 =  ind2orb[a1]
            s1 = summarize_orb(o1)
            for a2 = 1:nwan
                at2,t2,o2 =  ind2orb[a2]
                s2 = summarize_orb(o2)
                if at1 == at2 && s1 == s2
                    energy += 0.5 * Umat[a1] * (denmat[a1,a2, spin]*( iden[a2,a1] - denmat[a2,a1, spin]))
#                    println("add $a1 $a2 $([at1,t1,o1, at2,t2,o2]) $(Umat[a1] * (denmat[a1,a2, spin]*( iden[a2,a1] - denmat[a2,a1, spin])))  $([0.5 , Umat[a1] , denmat[a1,a2, spin], ( iden[a2,a1] - denmat[a2,a1, spin])])")
                end
            end
        end
    end
    if nspin == 1
        energy = energy * 2.0
    end
    for spin = 1:nspin
        for a1 = 1:nwan
            at1,t1,o1 =  ind2orb[a1]
            s1 = summarize_orb(o1)
            for a2 = 1:nwan
                at2,t2,o2 =  ind2orb[a2]
                s2 = summarize_orb(o2)
                if at1 == at2 && s1 == s2
                    potU[a1,a2,spin] = 0.5*Umat[a1] *( I[a1,a2] - 2.0 * denmat[a1,a2,spin])
                end
            end
        end
    end
    
    return energy, potU

end



end #end module

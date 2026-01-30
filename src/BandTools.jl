"""
    module BandTools

Utility functions for manipulating band structures. These need to be defined early in the code so other Modules have acccess to them.
"""
module BandTools

using LinearAlgebra
using SpecialFunctions

function gaussian(de, smearing=0.01)
    return 0.5 * (erfc.(de/smearing))
end    

function gaussian_derivative(de, smearing=0.01)
    return (-0.5 / sqrt(pi)) * exp.(-1.0*(de/smearing).^2)
end    


"""
    function calc_fermi(eigs, weights, nelec, smearing = 0.01)

calculate fermi energy using bisection
"""
function calc_fermi(eigs, weights, nelec, smearing = 0.01)


#    println("calc_fermi eigs ", size(eigs))

    efermi_max = maximum(eigs[:])+0.1
    efermi_min = minimum(eigs[:])-0.1    

    efermi = (efermi_max+efermi_min)/2.0

    norm = sum(weights)

    n2= Int64(round(nelec/2))

#    println("n2 $n2 size(eigs) $(size(eigs))")
    
   
    #handle case where all bands are full
    if length(size(eigs)) == 2 && size(eigs)[2] > 1
        if size(eigs)[2] <= (nelec/2.0 + 1e-10)
            efermi = efermi_max
            println("case 1")
            return maximum(eigs[:])
        end
    else
        println("lenght(eigs) $(length(eigs)) (nelec + 1e-10) $((nelec + 1e-10))")
        if length(eigs) <= (nelec/2.0 + 1e-10)
            efermi = efermi_max
            println("case 2")            
            return maximum(eigs[:])
        end
    end

    #handle case where all bands are empty
    if nelec <= 1e-10
        println("case 3")        
        return minimum(eigs[:])
    end
    
    
#    println(size(eigs), " n2 $n2")

    if  abs(nelec - n2*2) < 1e-7
        if length(size(eigs)) == 2 && size(eigs)[2] > 1
            gap  = minimum(eigs[:,n2+1]) - maximum(eigs[:,n2])   #large gaps become numerically unstable with small smearing
        else
            gap  = minimum(eigs[n2+1]) - maximum(eigs[n2])   #large gaps become numerically unstable with small smearing
        end
        
        if gap > 0.05  #if gap is positive/large and we have an integer number of electrons
            smearing = max(smearing, gap / 10.0)
        end
        #    println("gap $gap smearing $smearing n2 $n2")
        
    end

#    println("do iter ", smearing)
    
    for iter = 1:40
#        println("iter  $iter ")
        efermi = (efermi_max+efermi_min)/2.0

        occ = gaussian.(eigs.-efermi, smearing)

        n = sum(sum(occ, dims=2).*weights)/norm
#        println(iter, "  ", n,"  " ,  efermi)
        if n > (nelec/2.0)
            efermi_max = efermi
        else
            efermi_min = efermi
        end
    end

    return efermi
    
end



function calc_fermi_sp(eigs, weights, nelec, smearing = 0.01)
    nspin = 1
    if length(size(eigs)) == 3
        nspin = size(eigs)[3]
    end

#    println("FERMI nspin $nspin")
#    println("size eigs ", size(eigs))
    
    if nspin == 2

        nw = size(eigs)[2]
        nk = size(eigs)[1]
        eigs2 = zeros(nk, nw*nspin)
        for spin = 1:nspin
            eigs2[1:nk,1:nw] = eigs[:,:,1] 
            eigs2[1:nk,nw+1:2*nw] = eigs[:,:,2] 
        end
        eigs2 = sort(eigs2, dims=2)
#        println("eigs2 ")
#        println(eigs2)
        return calc_fermi(eigs2, weights, nelec*2, smearing )
    else
        return calc_fermi(eigs[:,:,1], weights, nelec, smearing )
    end        
end

"""
    function band_energy(eigs, weights, nelec, smearing = 0.01; returnk=false, returnocc=false, returnef=false, returnboth=false)

Calculate band energy. Has options for additional return variables. Calculates fermi energy internally.
"""
function band_energy(eigs, weights, nelec, smearing = 0.01; returnk=false, returnocc=false, returnef=false, returnboth=false)

    efermi = calc_fermi_sp(eigs, weights, nelec, smearing)
#    println("band_energy efermi $efermi")
    norm = sum(weights)
    
    occ = gaussian.(eigs.-efermi, smearing)

    nspin = 1
    if length(size(eigs)) == 3
        nspin = size(eigs)[3]
    end

    if nspin == 1
        ek = sum(eigs.*occ, dims=2).*weights * 2.0 / norm  
    elseif nspin == 2
        ek = sum(eigs.*occ, dims=[2,3]).*weights  / norm 
    else
        println("error nspin $nspin")
    end

    #    energy = sum(sum(eigs.*occ, dims=2).*weights) * 2.0 / norm
    energy=sum(ek)

    if returnk
        return energy,ek
    elseif returnocc
        return energy, occ
    elseif returnef
        return energy, efermi
    elseif returnboth
        return energy, efermi, occ
    else
        return energy
    end
    
end

"""
    function smearing_energy(eigs, weights, efermi, smearing = 0.01)

Smearing contribution to total energy. If you don't take this int account in metals or 
small gap semiconductors, your energy is not variational 
in the correct way, and things like forces and stresses become wrong.
"""
function smearing_energy(eigs, weights, efermi, smearing = 0.01)

#    efermi = calc_fermi(eigs, weights, nelec, smearing)
    norm = sum(weights)
    
    dgauss = gaussian_derivative(eigs.-efermi, smearing)

    nspin = 1
    if length(size(eigs)) == 3
        nspin = size(eigs)[3]
    end

    if nspin == 1
        ek = smearing * sum(dgauss, dims=2).*weights * 2.0 / norm  
    else
        ek = smearing * sum(dgauss, dims=[2,3]).*weights * 1.0 / norm  
    end
    
    energy=sum(ek)

    return energy

    
end


end 

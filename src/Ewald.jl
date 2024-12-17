###include("Crystal.jl")
###include("DFToutMod.jl")
#using XMLDict



####################### Wannier90 specific 
"""
    module Ewald

Module for electrostatic preperation.
"""
module Ewald
"""
Electrostatics
"""

#include("Atomdata.jl")
using ..Atomdata:atoms
using ..Atomdata:uniform_charge_interaction

using LinearAlgebra
using ..CrystalMod:crystal
###using ..CalcTB:distances_etc_3bdy
using SpecialFunctions
import Base.Threads.@spawn
using Base.Threads
using LoopVectorization


EWALD_PARAMS = [0.0,0.0,0.0]
EWALD_U = [0.0,0.0,0.0]
EWALD_FACTOR = [1.0]
"""
   function getU(types)

Get values of `U` from `Atomdata`.
"""
function getU(types, var_type=Float64)


    ntot = 0
    Usize = zeros(length(types))
    for (i,t) in enumerate(types)
        at1 = atoms[t ]
        nw1 = Int64(at1.nwan/2)
        ntot += nw1
        Usize[i] += at1.U
    end

    U = zeros(var_type, ntot, ntot)

    counter = 0
    for t1 in types
        at1 = atoms[t1 ]
        nw1 = Int64(at1.nwan/2)
        for i = 1:nw1
            for j = 1:nw1
                U[i+counter,j+counter] = at1.Umat[i,j]
            end
        end
        counter += nw1
    end
    

    return U, Usize

end

function getU3(types, var_type=Float64)

    U3 = zeros(var_type, length(types))
    for (c,t) in enumerate(types)
        U3[c] = atoms[t].U3
#        println("add U3 t ", U3[c], " " , atoms[t].U3 , " xxxxxxxxxxxxxxxxxxxxxxxxxxx")
    end

    return U3

end

"""
    function get_onsite(crys::crystal, U::Array{Float64,1})

return Onsite terms, Coulumb `U`.
"""
#function get_onsite(crys::crystal, U::Array{Float64,1})

#    return diagm(U)
#    gamma_onsite = zeros(Float64, crys.nat, crys.nat)
#
#    for i = 1:crys.nat
#        gamma_onsite[i,i] = U[i]
#    end##
#
#    return gamma_onsite

#end

"""
    function electrostatics_getgamma(crys::crystal;  kappa=missing, noU=false, onlyU=false, screening = 1.0)

Main function. Does Ewald calculation on `crys`. Returns `gamma`, which is used in scf calculation.
This is only run once for a given `tb_crys` object and stored.

- `kappa` is the splitting parameter between real/k-space in Ewald calculation. Will estimate best one if not provided.
- `noU=false` for testing only
- `onlyU=false` for testing only
- `screening=1.0` Not used. Purpose is to reduce U values for values < 1.
"""
function electrostatics_getgamma(crys::crystal;  kappa=missing, noU=false, onlyU=false, screening = 1.0, factor = 1.0)
#noU and onlyU are for testing purposes
    factor = EWALD_FACTOR[1]
    println("factor $factor")
#R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero = distances_etc_3bdy(crys,cutoff2X, 0.0)

#    println("EWALD")
#    println(crys)
#    println("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE")

#    noU = true
    
#    println("estimate")
    if ismissing(kappa)
#        kappa_default = 0.25
        kappa = estimate_best_kappa(crys.A)
#        println("kappa $kappa")
    end
    kappa = Float64(kappa)

#    println("kappa ", kappa)
    
    T = typeof(crys.coords[1,1])

    if noU
        println("noU - FOR TESTING")
        #U = zeros(Float64, crys.nat)
        U, Usize = getU(crys.types, T)
        U .= 0.0
    else
        U, Usize = getU(crys.types, T)
        U = U * screening
    end

    starting_size_rspace = 1
    starting_size_kspace = 1


    gamma_rs = zeros(T, crys.nat, crys.nat)
    gamma_U = zeros(T, crys.nat, crys.nat)
    gamma_k = zeros(T, crys.nat, crys.nat)
    gamma_self = zeros(T, crys.nat, crys.nat)


    #can do these at same time.
    # can run real space and k-space in parallel with asyncronous parallelization. Only a minor improvement.
#    println("rs")
    rs = begin 
        #gamma_rs, gamma_U, gamma_U_opt = real_space_LV(crys, kappa, [0.05, 0.05], starting_size_rspace, factor)
        gamma_rs, gamma_U, gamma_U_opt = real_space_LV(crys, kappa, max.(Usize, 0.2), starting_size_rspace, factor)
        
        #gamma_rs, gamma_U = real_space(crys, kappa, U, starting_size_rspace)
#        println("gamma_U_opt ", gamma_U_opt)
    end

    #gamma_onsiteU = get_onsite(crys, U)
    #for i = 1:crys.nat
    #    if crys.stypes[i] != :F
    #        println("$i old $(gamma_onsiteU[i]) add $(gamma_U_opt[i]) new $(gamma_onsiteU[i] + gamma_U_opt[i])")
    #        gamma_onsiteU[i] += gamma_U_opt[i]##
    #
    #    end
    #end
    
#    println("gamma_rs")
#    println(gamma_rs)
    #println("ks")
    ks = begin 
    #ks = begin
        gamma_k = k_space_LV(crys, kappa, starting_size_kspace)
        #gamma_k = k_space(crys, kappa, starting_size_kspace)
    end
#    println("gamma_k")
#    println(gamma_k)
    #println("self")
    self = begin 
        #self
        for i = 1:crys.nat
            gamma_self[i,i] -= kappa / sqrt(pi) * 2.0
        end
    end

    

    
#    wait(ks)
#    wait(self)
#    wait(rs)
    
    if false #for debugging
        println("gamma_rs")
        println(gamma_rs)
        println("gamma_k")
        println(gamma_k)
        println("gamma_self")
        println(gamma_self)
#        println("gamma_onsiteU")
#        println(gamma_onsiteU)
#        println("gamma_U")
#        println(gamma_U)
        println()
        println("only 1/r")
        println(gamma_rs + gamma_k + gamma_self)
    end
#    println("gamma_rs + gamma_k + gamma_self + gamma_U ", [gamma_rs ,  gamma_k ,  gamma_self ,  gamma_U])
    
    #rydberg units
    e2 = 2.0

    #    gamma_tot = e2*(gamma_rs + gamma_k + gamma_self + gamma_U) + gamma_onsiteU
    gamma_tot = e2*(gamma_rs + gamma_k + gamma_self + gamma_U) ####+ gamma_onsiteU      #xxyy

    println("gamma_rs" , gamma_rs)
    println("gamma_k", gamma_k)
    println("gamma_self", gamma_self)
    println("gamma_U", gamma_U)
   
    #interaction of net charge with uniform background
    background_charge_correction = -e2 * pi  / (2 * abs(det(crys.A)) * kappa^2) 

    gamma_bc = zeros(T, crys.nat, crys.nat)
    for i = 1:crys.nat
        for j = 1:crys.nat
            if i == j
                gamma_bc[i,j] += background_charge_correction
            else
                gamma_bc[i,j] += background_charge_correction
            end                
        end
    end

    gamma_tot += gamma_bc *2.0

    gamma_tot_expand = U
    println("gamma_tot_expand U ", gamma_tot_expand)
#    for a1 = 1:crys.nat
#        for a2 = 1:crys.nat
    if  !onlyU #for debugging

        o1 = 1
        for i = 1:crys.nat
            at1 = atoms[crys.types[i]  ]
            nw1 = Int64(at1.nwan/2)
            o2 = 1
            for j = 1:crys.nat
                at2 = atoms[crys.types[j]]
                nw2 = Int64(at2.nwan/2)
                for c1 = o1:o1+nw1-1
                    for c2 = o2:o2+nw2-1
                        gamma_tot_expand[c1,c2] += gamma_tot[i,j]
                    end
                end
                o2 += nw2
                
            end
            o1 += nw1
        end
    end

    
    background_charge_correction = 0.0
    for t = crys.stypes
        if t in keys(uniform_charge_interaction)
            background_charge_correction += uniform_charge_interaction[t]
        end
    end
    background_charge_correction = background_charge_correction / abs(det(crys.A)) / crys.nat
       

    


    #background_charge_correction = e2 * pi  / (2 * abs(det(crys.A)) * kappa^2)


    #onsite only U3
    U3 = zeros(T, crys.nat)
    U3 .= getU3(crys.stypes, T)
    
    #   return gamma_tot, background_charge_correction
    return gamma_tot_expand, background_charge_correction, U3

end

"""
    function real_space(crys::crystal, kappa::Float64, U::Array{Float64}, starting_size_rspace=2)

Real-space Ewald sum.
"""
function real_space(crys::crystal, kappa::Float64, U::Array{Float64}, starting_size_rspace=2)
    
    T = typeof(crys.coords[1,1])

    first_iter = true
    old_val = 0.0
    
    gamma_ij_tot = zeros(T, crys.nat, crys.nat)
    gamma_ij_new = zeros(T, crys.nat, crys.nat)

    gamma_U_tot = zeros(T, crys.nat, crys.nat)
    gamma_U_new = zeros(T, crys.nat, crys.nat)

    Uconst = zeros(T, crys.nat, crys.nat)

    if sum(abs.(U)) > 1e-5
        useU = true
        for i = 1:crys.nat
            Fi = sqrt( 8* log(2)/pi ) / (U[i]/2.0) #rydberg units, U in Ryd , we divide by 2, and multiply by e^2 layer. For Hartree formula, see koskinen comp mater sci 47 (2009) 237
            for j = 1:crys.nat
#                Uconst[i,j] = sqrt(pi/2 * (U[i]^2 * U[j]^2 / (U[i]^2 + U[j]^2)))

                Fj = sqrt( 8* log(2)/pi ) / (U[j]/2.0)
                Uconst[i,j] = sqrt(4 * log(2) / (Fi^2 + Fj^2))
            end
        end
    else
        useU = false
    end

    #    R = zeros(T,1, 3)
    R = zeros(T, 3)
    Ra = zeros(T, 3)
    Rb = zeros(T, 3)
    
    coords_cart = crys.coords * crys.A
    coords_cartT = coords_cart'

    coords_cartTij = zeros(T, 3, crys.nat, crys.nat)
    for i = 1:crys.nat
        for j = 1:crys.nat
            coords_cartTij[:,i,j] = (@view coords_cartT[:,i]) - (@view coords_cartT[:,j]) 
        end
    end
    
    
    converged = false
    converged_old = false
    R0 = false

    
    N_list = [[1,1,1]]
    a1 = sqrt(sum(crys.A[1,:].^2))
    a2 = sqrt(sum(crys.A[2,:].^2))
    a3 = sqrt(sum(crys.A[3,:].^2))
    a = (minimum([a1,a2,a3]) - 0.0010) * 0.95

#    println("a $a ", [a1, a2, a3])
    
    for N = 1:70

        n1 = max(Int64(ceil( a * N / a1)), 2)
        n2 = max(Int64(ceil( a * N / a2)), 2)
        n3 = max(Int64(ceil( a * N / a3)), 2)

        if [n1,n2,n3] != N_list[end]
            N_list = push!(N_list, [n1,n2,n3])
#            println("N_list ", (n1,n2,n3))
        end
    end

    
    newcontr = 0.0
    newcontrU = 0.0
    
#    for N = starting_size_rspace:20
#        Nx = N
#        Ny = N
#        Nz = N
    Nxold = 0
    Nyold = 0
    Nzold = 0

    At = crys.A'
    
    for (Nx, Ny, Nz) in N_list
#        println("N ", [Nx, Ny, Nz])
        gamma_ij_new[:,:] .= 0.0
        gamma_U_new[:,:] .= 0.0
        for x = -Nx:Nx
            Ra[1] = x
            for y = -Ny:Ny
                Ra[2] = y
                for z = -Nz:Nz
                    Ra[3] = z
                    if x != 0 || y != 0 || z != 0
                        R0 = false
                    else
                        R0 = true
                    end
                    R[:] .= At * Ra 
                    if first_iter == true || (abs(x) > Nxold || abs(y) > Nyold || abs(z) > Nzold )
                        for i = 1:crys.nat
                            for j = i:crys.nat
                                #                                Rb = (@view coords_cartT[:,i]) - (@view coords_cartT[:,j]) - (@view R[:])
                                Rb = (@view coords_cartTij[:,i, j])  - (@view R[:])
                                r = (Rb'*Rb).^0.5
#                                r = sum(( (@view coords_cartT[:,i]) - (@view coords_cartT[:,j]) - (@view R[:]) ).^2)^0.5
                                #r = sum(( (coords_cart[i,:]) - (coords_cart[j,:]) - ( R[1,:] ) ).^2)^0.5
                                if i != j ||  R0 == false
                                    @inbounds gamma_ij_new[i,j] += erfc( kappa * r) / r

                                    if useU
                                        @inbounds gamma_U_new[i,j] += -erfc( Uconst[i,j] * r) / r  #see eq 3 in prb 66 075212, or koskinen comp mater sci 47 (2009) 237

                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        gamma_U_new = gamma_U_new + gamma_U_new' - diagm(diag(gamma_U_new))
        gamma_ij_new = gamma_ij_new + gamma_ij_new' - diagm(diag(gamma_ij_new))
        
        newcontr = sum(abs.(gamma_ij_new))
        newcontrU = sum(abs.(gamma_U_new))

#        println("new $Nx $Ny $Nz real_space $newcontr $newcontrU")

        if newcontr < 1e-5 && newcontrU < 1e-5
            if converged_old == true || (newcontr < 1e-7 && newcontrU < 1e-7)
                converged = true
                break
            end
            converged_old = true
        end
        first_iter = false
        
        gamma_ij_tot += gamma_ij_new
        gamma_U_tot += gamma_U_new

        Nxold = Nx
        Nyold = Ny
        Nzold = Nz


    end

    if converged == false
        println("WARNING, real_space EWALD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        println("real_space NOT converged  : $newcontr    $newcontrU")
        println(crys)
        
    end

    return gamma_ij_tot, gamma_U_tot

end


function real_space_LV(crys::crystal, kappa::Float64, U::Array{Float64}, starting_size_rspace=2, factor=1.0)
    
    T = typeof(crys.coords[1,1])

    first_iter = true
    old_val = 0.0

    gamma_U_opt = zeros(T, crys.nat)
    
    gamma_ij_tot = zeros(T, crys.nat, crys.nat)
    gamma_ij_new = zeros(T, crys.nat, crys.nat)

    gamma_U_tot = zeros(T, crys.nat, crys.nat)
    gamma_U_new = zeros(T, crys.nat, crys.nat)
 
    gamma_ij_new_th = zeros(T, crys.nat, crys.nat, nthreads())
    gamma_U_new_th = zeros(T, crys.nat, crys.nat, nthreads())
   
    Uconst = zeros(T, crys.nat, crys.nat)

    if sum(abs.(U)) > 1e-5
        useU = true
        for i = 1:crys.nat
            Fi = sqrt( 8* log(2)/pi ) / (U[i]/2.0) #rydberg units, U in Ryd , we divide by 2, and multiply by e^2 layer. For Hartree formula, see koskinen comp mater sci 47 (2009) 237
            for j = 1:crys.nat
#                Uconst[i,j] = sqrt(pi/2 * (U[i]^2 * U[j]^2 / (U[i]^2 + U[j]^2)))

                Fj = sqrt( 8* log(2)/pi ) / (U[j]/2.0)
                Uconst[i,j] = sqrt(4 * log(2) / (Fi^2 + Fj^2)) * factor
            end
        end
    else
        useU = false
        Uconst .= 300.0 #for debug
    end

#    println("Uconst ")
#    println(Uconst)
    
    #    R = zeros(T,1, 3)
    R = zeros(T, 3)
    Ra = zeros(T, 3)
    Rb = zeros(T, 3)
    
    coords_cart = crys.coords * crys.A
    coords_cartT = coords_cart'

    coords_cartTij = zeros(T, 3, crys.nat, crys.nat)
    for i = 1:crys.nat
        for j = 1:crys.nat
            coords_cartTij[:,i,j] = (@view coords_cartT[:,i]) - (@view coords_cartT[:,j]) 
        end
    end
    
    
    converged = false
    converged_old = false
    R0 = false

    
    N_list = [[1,1,1]]
    a1 = sqrt(sum(crys.A[1,:].^2))
    a2 = sqrt(sum(crys.A[2,:].^2))
    a3 = sqrt(sum(crys.A[3,:].^2))
    a = (minimum([a1,a2,a3]) - 0.0010) * 0.95

#    println("a $a ", [a1, a2, a3])
    
    for N = 1:70

        n1 = max(Int64(ceil( a * N / a1)), 2)
        n2 = max(Int64(ceil( a * N / a2)), 2)
        n3 = max(Int64(ceil( a * N / a3)), 2)

        if [n1,n2,n3] != N_list[end]
            N_list = push!(N_list, [n1,n2,n3])
#            println("N_list ", (n1,n2,n3))
        end
    end

    
    newcontr = 0.0
    newcontrU = 0.0
    
#    for N = starting_size_rspace:20
#        Nx = N
#        Ny = N
#        Nz = N
    Nxold = 0
    Nyold = 0
    Nzold = 0

    At = crys.A'
    
    for (Nx, Ny, Nz) in N_list
#        println("N ", [Nx, Ny, Nz])
        gamma_ij_new[:,:] .= 0.0
        gamma_U_new[:,:] .= 0.0

        gamma_ij_new_th[:,:,:] .= 0.0
        gamma_U_new_th[:,:,:] .= 0.0

        for x = -Nx:Nx
            Ra[1] = x
            for y = -Ny:Ny
                Ra[2] = y
                for z = -Nz:Nz
                    Ra[3] = z
                    if x != 0 || y != 0 || z != 0
                        R0 = false
                    else
                        R0 = true
                    end
                    R[:] .= At * Ra 
                    if first_iter == true || (abs(x) > Nxold || abs(y) > Nyold || abs(z) > Nzold )
                        if R0
                            rs_core(crys, coords_cartTij, R, kappa, Uconst, gamma_ij_new, gamma_U_new, useU, R0, gamma_U_opt)
                        else 
                            rs_core_R0false_thread(crys.nat, coords_cartTij, R, kappa, Uconst, gamma_ij_new_th, gamma_U_new_th, useU, R0,gamma_U_opt)
                        end

                    end
                end
            end
        end
        gamma_U_new = gamma_U_new + sum(gamma_U_new_th, dims=3)[:,:]
        gamma_ij_new = gamma_ij_new + sum(gamma_ij_new_th, dims=3)[:,:]

        gamma_U_new = gamma_U_new + gamma_U_new' - diagm(diag(gamma_U_new))
        gamma_ij_new = gamma_ij_new + gamma_ij_new' - diagm(diag(gamma_ij_new))
        
        newcontr = sum(abs.(gamma_ij_new))
        newcontrU = sum(abs.(gamma_U_new))

#        println("new $Nx $Ny $Nz real_space $newcontr $newcontrU")

        if newcontr < 1e-5 && newcontrU < 1e-5
            if converged_old == true || (newcontr < 1e-7 && newcontrU < 1e-7)
                converged = true
                break
            end
            converged_old = true
        end
        first_iter = false
        
        gamma_ij_tot += gamma_ij_new
        gamma_U_tot += gamma_U_new

        Nxold = Nx
        Nyold = Ny
        Nzold = Nz


    end

    if !(useU)
        gamma_U_tot .= 0.0
    end
    
    if converged == false
        println("WARNING, real_space EWALD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        println("real_space NOT converged  : $newcontr    $newcontrU")
        println(crys)
        
    end

    return gamma_ij_tot, gamma_U_tot, gamma_U_opt 

end

function rs_core(crys, coords_cartTij, R, kappa, Uconst, gamma_ij_new, gamma_U_new, useU, R0, gamma_U_opt)

    @inbounds @simd for i = 1:crys.nat
        for j = i:crys.nat
            Rb = (@view coords_cartTij[:,i, j])  - (@view R[:])
            r = (Rb'*Rb).^0.5
            if i != j ||  R0 == false
                gamma_ij_new[i,j] += erfc( kappa * r) / r
#                gamma_ij_new[i,j] += 5*exp(-2.5*0.5 * r)* ( EWALD_PARAMS[1] + (1 - 2.5*r) * EWALD_PARAMS[2] + 0.5*((2.5*r).^2 .- 4.0*(2.5*r) .+ 2) *  EWALD_PARAMS[3])

                gamma_U_new[i,j] += -erfc( Uconst[i,j] * r) / r  #see eq 3 in prb 66 075212, or koskinen comp mater sci 47 (2009) 237

                #gamma_U_opt[i] += exp(-2.5*0.5 * r)* ( EWALD_U[1] + (1 - 2.5*r) * EWALD_U[2] + 0.5*((2.5*r).^2 .- 4.0*(2.5*r) .+ 2) *  EWALD_U[3])
                
            end
        end
    end

end

function rs_core_R0false(nat, coords_cartTij, R, kappa, Uconst, gamma_ij_new, gamma_U_new, useU, R0, gamma_U_opt)
    @inbounds @simd for i = 1:nat
        for j = i:nat
            r=0.0
            for a = 1:3
                r += (coords_cartTij[a,i, j]  - R[a])^2
            end
            r=r^0.5
            gamma_ij_new[i,j] += erfc( kappa * r) / r
#            gamma_ij_new[i,j] += exp(-2.5*0.5 * r)* ( EWALD_PARAMS[1] + (1 - 2.5*r) * EWALD_PARAMS[2] + 0.5*((2.5*r).^2 .- 4.0*(2.5*r) .+ 2) *  EWALD_PARAMS[3])

            gamma_U_new[i,j] += -erfc( Uconst[i,j] * r) / r  #see eq 3 in prb 66 075212, or koskinen comp mater sci 47 (2009) 237
#            gamma_U_opt[i] += exp(-2.5*0.5 * r)* ( EWALD_U[1] + (1 - 2.5*r) * EWALD_U[2] + 0.5*((2.5*r).^2 .- 4.0*(2.5*r) .+ 2) *  EWALD_U[3])

        end
    end
    
end

function rs_core_R0false_thread(nat, coords_cartTij, R, kappa, Uconst, gamma_ij_new, gamma_U_new, useU, R0, gamma_U_opt)
    @inbounds for i = 1:nat #@threads
        id = threadid()
        for j = i:nat
           r=0.0
           for a = 1:3
               r += (coords_cartTij[a,i, j]  - R[a])^2
           end
           r=r^0.5
           gamma_ij_new[i,j,id] += erfc( kappa * r) / r
            gamma_U_new[i,j,id] += -erfc( Uconst[i,j] * r) / r  #see eq 3 in prb 66 075212, or koskinen comp mater sci 47 (2009) 237

#            gamma_ij_new[i,j,id] += exp(-2.5*0.5 * r)* ( EWALD_PARAMS[1] + (1 - 2.5*r) * EWALD_PARAMS[2] + 0.5*((2.5*r).^2 .- 4.0*(2.5*r) .+ 2) *  EWALD_PARAMS[3])
#            gamma_U_opt[i] += exp(-2.5*0.5 * r)* ( EWALD_U[1] + (1 - 2.5*r) * EWALD_U[2] + 0.5*((2.5*r).^2 .- 4.0*(2.5*r) .+ 2) *  EWALD_U[3])
            
       end
    end
    
end


"""
    function estimate_best_kappa(A)

Estimate best value of `kappa` for Ewald sum.
Shouldn't effect final value, only calculation speed.
There is probably a better way to do this.
"""
function estimate_best_kappa(A)

    a1 = sqrt(sum(A[1,:].^2))
    a2 = sqrt(sum(A[2,:].^2))
    a3 = sqrt(sum(A[3,:].^2))

    #a = minimum([a1,a2,a3])
    a = (a1+a2+a3)/3

#    println("a ", a)
#    println(a1)
#    println(a2)
#    println(a3)

    B = inv(A)'

    b1 = sqrt(sum(B[1,:].^2))
    b2 = sqrt(sum(B[2,:].^2))
    b3 = sqrt(sum(B[3,:].^2))
    
#    b = minimum([b1,b2,b3])
    b = (b1+b2+b3)/3
    
#    println("b ", b)
#    println(b1)
#    println(b2)
#    println(b3)


    tot = Float64[]
    kappa_test = [0.00002 0.0001 0.0005 0.001 0.002 0.005 0.01 0.015 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.12 0.15 0.17 0.2 0.22 0.25 0.27  0.3 0.33 0.35 0.4 0.5 0.65 0.75 0.80 0.90 1.0 1.1 1.3 2.0 5.0 10.0 20.0 50.0 100.0]
    for kappa = kappa_test
        rs = erfc( kappa * a) / a * 200
        ks = exp(-b^2 / 4.0 / kappa^2  * (2 * pi)^2 ) / (2*pi*b)^2
        push!(tot, rs + ks)
#        println("$kappa rs $rs ks $ks  tot ",rs + ks)
    end
    i = argmin(tot)
    kappa = kappa_test[i]


    #println("bestkappa :", kappa)
    
    return kappa
end

"""
    function k_space(crys::crystal, kappa, starting_size_kspace=2)

K-space Ewald sum.
"""
function k_space(crys::crystal, kappa, starting_size_kspace=2)

    T = typeof(crys.coords[1,1])

    first_iter = true
    old_val = 0.0
    
    gamma_ij_tot = zeros(T, crys.nat, crys.nat)
    gamma_ij_new = zeros(T, crys.nat, crys.nat)
    K = zeros(T, 1,3)
    
    coords_cart = crys.coords * crys.A
    
    converged = false
    converged_old = false
    
    B = inv(crys.A)'
    vol = abs(det(crys.A))

    N_list = [[1,1,1]]
    b1 = sqrt(sum(B[1,:].^2))
    b2 = sqrt(sum(B[2,:].^2))
    b3 = sqrt(sum(B[3,:].^2))
    b = (minimum([b1,b2,b3]) - 0.0010) * 0.95

    for N = 1:30

        n1 = max(Int64(ceil( b * N / b1)), 2)
        n2 = max(Int64(ceil( b * N / b2)), 2)
        n3 = max(Int64(ceil( b * N / b3)), 2)

        if [n1,n2,n3] != N_list[end]
            N_list = push!(N_list, [n1,n2,n3])
#            println("k list ", [n1,n2,n3])
        end
    end


    newcontr = 0.0
    
    Nxold = 0
    Nyold = 0
    Nzold = 0

    for (Nx, Ny, Nz) in N_list
#        println("k ", [Nx, Ny, Nz])
#    for N = starting_size_kspace:25
        gamma_ij_new[:,:] .= 0.0
        for kx = -Nx:Nx
            for ky = -Ny:Ny
                for kz = -Nz:Nz
                    if kx == 0 && ky == 0 && kz == 0
                        continue
                    end
                    K[1,:] .= [kx,ky,kz] 
                    K = K * B
                    k2 = sum(K.^2) * (2*pi)^2
                    factor_k = (2*pi) * exp(-k2 / 4.0 / kappa^2  ) / k2

                    if first_iter == true || abs(kx) > Nxold || abs(ky) > Nyold || abs(kz) > Nzold
                        for i = 1:crys.nat
                            for j = i:crys.nat
                                kr = (K * ( (@view coords_cart[i,:]) - (@view coords_cart[j,:]) ) )[1]

                                exp_c = exp(2*pi*im*kr)
                                temp = real(exp_c )
                                gamma_ij_new[i,j] += factor_k * temp

                            end
                        end
                    end
                    
                end
            end
        end
        gamma_ij_new = gamma_ij_new + gamma_ij_new' - diagm(diag(gamma_ij_new))
        first_iter = false
        newcontr = sum(abs.(gamma_ij_new))

        if newcontr < 1e-5
            if converged_old == true || newcontr < 1e-7
                converged = true
                break
            end
            converged_old = true #need convergence twice

        end

        gamma_ij_tot += 2.0*gamma_ij_new /vol
        
        Nxold = Nx
        Nyold = Ny
        Nzold = Nz

    end

    if converged == false
        println("WARNING, k_space EWALD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        println("newcontr $newcontr")
        println(crys)
    end

    return gamma_ij_tot

end



"""
    function k_space(crys::crystal, kappa, starting_size_kspace=2)

K-space Ewald sum.
"""
function k_space_LV(crys::crystal, kappa, starting_size_kspace=2)

    T = typeof(crys.coords[1,1])

    first_iter = true
    old_val = 0.0
    
    gamma_ij_tot = zeros(T, crys.nat, crys.nat)
    gamma_ij_new = zeros(T, crys.nat, crys.nat)
    gamma_ij_new_th = zeros(T, crys.nat, crys.nat, nthreads())    
    K = zeros(T, 1,3)
    
    coords_cart = crys.coords * crys.A
    
    converged = false
    converged_old = false
    
    B = inv(crys.A)'
    vol = abs(det(crys.A))

    N_list = [[1,1,1]]
    b1 = sqrt(sum(B[1,:].^2))
    b2 = sqrt(sum(B[2,:].^2))
    b3 = sqrt(sum(B[3,:].^2))
    b = (minimum([b1,b2,b3]) - 0.0010) * 0.95

    for N = 1:30

        n1 = max(Int64(ceil( b * N / b1)), 2)
        n2 = max(Int64(ceil( b * N / b2)), 2)
        n3 = max(Int64(ceil( b * N / b3)), 2)

        if [n1,n2,n3] != N_list[end]
            N_list = push!(N_list, [n1,n2,n3])
#            println("k list ", [n1,n2,n3])
        end
    end


    newcontr = 0.0
    
    Nxold = 0
    Nyold = 0
    Nzold = 0

    for (Nx, Ny, Nz) in N_list
#        println("k ", [Nx, Ny, Nz])
#    for N = starting_size_kspace:25
        gamma_ij_new[:,:] .= 0.0
        gamma_ij_new_th[:,:,:] .= 0.0
        for kx = -Nx:Nx
            for ky = -Ny:Ny
                for kz = -Nz:Nz
                    if kx == 0 && ky == 0 && kz == 0
                        continue
                    end
                    K[1,:] .= [kx,ky,kz] 
                    K = K * B
                    k2 = sum(K.^2) * (2*pi)^2
                    factor_k = (2*pi) * exp(-k2 / 4.0 / kappa^2  ) / k2

                    if first_iter == true || abs(kx) > Nxold || abs(ky) > Nyold || abs(kz) > Nzold
                        ks_core(crys, coords_cart, K, factor_k, gamma_ij_new_th)
                    end
                    
                end
            end
        end
        gamma_ij_new = gamma_ij_new + sum(gamma_ij_new_th, dims=3)[:,:]
        
        gamma_ij_new = gamma_ij_new + gamma_ij_new' - diagm(diag(gamma_ij_new))
        first_iter = false
        newcontr = sum(abs.(gamma_ij_new))

        if newcontr < 1e-5
            if converged_old == true || newcontr < 1e-7
                converged = true
                break
            end
            converged_old = true #need convergence twice

        end

        gamma_ij_tot += 2.0*gamma_ij_new /vol
        
        Nxold = Nx
        Nyold = Ny
        Nzold = Nz

    end

    if converged == false
        println("WARNING, k_space EWALD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        println("newcontr $newcontr")
        println(crys)
    end

    return gamma_ij_tot

end


function ks_core(crys, coords_cart, K, factor_k, gamma_ij_new_th)

    @inbounds @threads for i = 1:crys.nat
        id = threadid()
        for j = i:crys.nat
            kr = 0.0
            for a = 1:3
                kr += K[a] * ( coords_cart[i,a] - coords_cart[j,a])
            end
                       
            #            exp_c = exp(2*pi*im*kr)
            #            temp = real(exp_c )
            #gamma_ij_new_th[i,j,id] += factor_k * real(  exp(2*pi*im*kr)  )
            gamma_ij_new_th[i,j,id] += factor_k * cos(2*pi*kr)
            
        end
    end
    
end

end #end module

"""
Module for classical monte carlo finite temperature sampling.
"""
module MonteCarlo

using ..CrystalMod:crystal
using ..ThreeBodyTB:scf_energy
using ..ThreeBodyTB:set_units
using Suppressor
using LinearAlgebra

function run_mc(c_start::crystal, tempK; step_size = 0.1, adjust_step = true, adjust_strain = true, nsteps = 1000, nsteps_thermal = 100, database = missing, smearing=0.01, grid = missing, conv_thr = 2e-5, iters = 100, mix = -1.0, mixing_mode=:simple, nspin=1, eden=missing, verbose=false, repel=true, tot_charge=0.0, use_sym=true, do_classical=true, do_tb=true, database_classical=missing, sparse=:auto)

    old_units = set_units()
    set_units(both="atomic")
    
    #temperature in K to atomic units energy
    temp = tempK * 8.617333262 * 10^-5 / 13.6057039763
    beta = 1/ temp
    println("temp (K) $tempK, temp Ryd $temp, beta $beta")

    println("thermalization")
    step_size_strain = 0.01
    
    energies_thermal, c_thermal, step_size_thermal,step_size_strain_thermal = mc_helper(c_start, beta, true, step_size, step_size_strain, adjust_strain = adjust_strain, nsteps = nsteps_thermal, database = database, smearing=smearing, grid = grid, conv_thr = conv_thr, iters = iters, mix = mix, mixing_mode=mixing_mode,  nspin=nspin, eden=eden, verbose=verbose, repel=repel, tot_charge=tot_charge, use_sym=use_sym, do_classical=do_classical, do_tb=do_tb, database_classical= database_classical, sparse=sparse)

    println("c_thermal")
    println(c_thermal)
    println("energy ", energies_thermal[end])
    
    println("final run")
    energies, c_final, step_size, step_size_strain = mc_helper(c_thermal, beta, false, step_size_thermal, step_size_strain_thermal, adjust_strain = adjust_strain, nsteps = nsteps, database = database, smearing=smearing, grid = grid, conv_thr = conv_thr, iters = iters, mix = mix, mixing_mode=mixing_mode,  nspin=nspin, eden=eden, verbose=verbose, repel=repel, tot_charge=tot_charge, use_sym=use_sym, do_classical=do_classical, do_tb=do_tb, database_classical= database_classical, sparse=sparse)

    set_units(energy=old_units[1], length=old_units[2])

    
    println("-----------------------")
    println("c_final ")
    println(c_final)
    println("energy ", energies[end])
    
    
    return energies, c_final
end


function mc_helper(c_start, beta, adjust_step, step_size, step_size_strain ; adjust_strain = false, nsteps = 100, database = missing, smearing=0.01, grid = missing, conv_thr = 2e-5, iters = 100, mix = -1.0, mixing_mode=:simple, nspin=1, eden=missing, verbose=false, repel=true, tot_charge=0.0, use_sym=true, do_classical=true, do_tb=true, database_classical=missing, sparse=:auto)

    c_current = deepcopy(c_start)
    c_work = deepcopy(c_start)

    tbc = missing
    flag =true

    en, tbc, flag = scf_energy(c_start, database = database, smearing=smearing, grid = grid, conv_thr = conv_thr, iters = iters, mix = mix, mixing_mode=mixing_mode,  nspin=nspin, eden=eden, verbose=verbose, repel=repel, tot_charge=tot_charge, use_sym=use_sym, do_classical=do_classical, do_tb=do_tb, database_classical= database_classical, sparse=sparse)
    println("c start ")
    println(c_start)
    println("en start $en")
    en_new = 0.0
    energies = zeros(nsteps)


    #0 means strain
    if adjust_strain
        #atoms = 0:c_current.nat
        atoms = [0,1]
    else
        #atoms= 1:c_current.nat
        atoms = [1]
    end
    
    for step = 1:nsteps

#        accept = 0
#        reject = 0

        accept = []
#        reject = []
        
        for atom = atoms
#            atom_step = ((rand(1,3) .- 0.5)*step_size ) * Ainv
#            c_work.coords[atom, :] += atom_step[:]
#            println("atom $atom --------------------------------------------------------------------- $en")
    #        println("c_curent before")
    #        println(c_current)
    #        println("c_work before")
     #       println(c_work)

            generate_guess(atom, c_work, step_size, step_size_strain)

      #      println("c_curent after")
      #      println(c_current)
      #      println("c_work after")
       #     println(c_work)
            
            @suppress begin
                if !ismissing(tbc)
                    eden = tbc.eden
                else
                    eden = missing
                end
                    
                en_new, tbc, flag = scf_energy(c_work, database = database, smearing=smearing, grid = grid, conv_thr = conv_thr, iters = iters, mix = mix, mixing_mode=mixing_mode,  nspin=nspin, verbose=verbose, repel=repel, tot_charge=tot_charge, use_sym=use_sym, do_classical=do_classical, do_tb=do_tb, database_classical= database_classical, sparse=sparse, eden = tbc.eden)
            end
#            println("en_new $en_new")
            
            rand_num = rand(1)[1]
            W = min(1.0, exp(-beta * (en_new - en)))

            if flag == false
                W = 0.0
            end
            
            if rand_num < W #accept

                en = en_new
                c_current.coords[:,:] = c_work.coords
                c_current.A[:,:]   = c_work.A

                #accept += 1
                push!(accept, true)
                
            else #reject

                c_work.coords[:,:] = mod.(c_current.coords, 1.0)
                c_work.A[:,:] = c_current.A

               
                #                reject += 1
                push!(accept, false)
                
            end

#            println("c_current end of loop")
#            println(c_current)
#            println("c_work end of loop")
#            println(c_work)
#            println("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        end
        println("step $step accept $accept en $en step_size $step_size")        
        if adjust_step
            if accept[2]
                step_size = step_size * 1.08
            else
                step_size = step_size * 0.90
            end
        end                
        if adjust_step && adjust_strain
            if accept[1]
                step_size_strain = step_size_strain * 1.08
            else
                step_size_strain = step_size_strain * 0.90
            end
            step_size_strain = min(step_size_strain, 0.015)
        end                

        energies[step] += en
            
    end
#    println("c_current end ")
#    println(c_current)
    println("end mc helper xyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy")
    println()
    return energies, c_current, step_size, step_size_strain

end

function generate_guess(atom, c_work, step_size, step_size_strain)

    if atom >= 1

        Ainv = inv(c_work.A)
        atom_step = ((rand(c_work.nat,3) .- 0.5)*step_size ) * Ainv
        c_work.coords[:, :] += atom_step[:,:]
        c_work.coords = mod.(c_work.coords, 1.0)

    else #strain case

        strain = (rand(3,3) .- 0.5)
        strain = 0.5*(strain + strain') * step_size_strain
        c_work.A[:,:] = c_work.A*(I(3) + strain)

    end
    
end

end  #end module


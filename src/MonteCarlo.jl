"""
Module for classical monte carlo finite temperature sampling.
"""
module MonteCarlo

using ..CrystalMod:crystal



function run_mc(c_start::crystal, tempK; step_size = 0.01, adjust_step = true, adjust_strain = true, nsteps = 100,  database = missing, smearing=0.01, grid = missing, conv_thr = 2e-5, iters = 100, mix = -1.0, mixing_mode=:simple, nspin=1, eden=missing, verbose=false, repel=true, tot_charge=0.0, use_sym=true, do_classical=true, do_tb=true, database_classical=missing, sparse=:auto)

    #temperature in K to atomic units energy
    temp = tempK * 8.617333262 * 10^-5 * 13.6057039763
    beta = 1/ temp
    
    c_current = deepcopy(c_start)
    c_work = deepcopy(c_start)

    en, tbc, flag = scf_energy(c_start, database = database, smearing=smearing, grid = grid, conv_thr = conv_thr, iters = iters, mix = mix, mixing_mode=mising_mode,  nspin=nspin, eden=eden, verbose=verbose, repel=repel, tot_charge=tot_charge, use_sym=use_sym, do_classical=do_classical, do_tb=do_tb, database_classical= database_classical, sparse=sparse)

    en_start = en
    
    accept = 0
    reject = 0
    for step = 1:nsteps
        Ainv = inv(c_current.A)
        for atom = 1:c_current.nat
            atom_step = ((rand(1,3) .- 0.5)*step_size ) * Ainv
            c_work.coords[atom, :] += atom_step

            en_new, tbc, flag = scf_energy(c_work, database = database, smearing=smearing, grid = grid, conv_thr = conv_thr, iters = iters, mix = mix, mixing_mode=mising_mode,  nspin=nspin, eden=eden, verbose=verbose, repel=repel, tot_charge=tot_charge, use_sym=use_sym, do_classical=do_classical, do_tb=do_tb, database_classical= database_classical, sparse=sparse)

            rand_num = rand(1)[1]
            W = min(1.0, exp(-beta * (en_new - en)))
                    
            if rand_num < W #accept

                en = en_new
                c_current.coords[atom,:] = c_work.coords[atom,:]
                accept += 1
                
            else #reject
                c_work.coords[atom,:] .= c_current.coords[atom,:]
                reject += 1
            end
            if adjust_step
                if accept > reject
                    step_size = step_size * 1.05
                elseif accept < reject
                    step_size = step_size * 0.95
                end
            end                
            println("step $step accept $accept reject $reject en $en")
        end

            
    end
    

end
    
end  #end module


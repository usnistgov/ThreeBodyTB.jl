    function construct_newXY(VECTS_FITTED::Dict{Int64, Array{Complex{Float64},4} }, OCCS_FITTED::Array{Float64,4}, ncalc::Int64, ncols::Int64, nlam::Int64, ERROR::Array{Int64,1}, EDEN_FITTED::Array{Float64,3}; leave_out=-1)


        #        nlam = 0

        #        NEWX = zeros(ncalc*nk_max*NWAN_MAX + ncalc + nlam, ncols)
        #        NEWY = zeros(ncalc*nk_max*NWAN_MAX + ncalc + nlam)

        counter = 0

        temp = 0.0+0.0*im
        
        nonzero_ind = Int64[]
        #count nonzero ahead of time
        for calc = 1:ncalc
            if calc == leave_out #skip this one
                continue
            end
            if ERROR[calc] == 1
                continue
            end
            row1, rowN, nw = ind_BIG[calc, 1:3]
            N = hermetian_indpt(nw)
            vals = ones(nw) * 100.0
            nk = size(KPOINTS[calc])[1]
            for spin = 1:SPIN[calc]
                for k = 1:nk
                    for i = 1:nw
                        counter += 1
                        push!(nonzero_ind, counter)
                    end
                end
            end
            counter += 1
            push!(nonzero_ind, counter)
        end

        NEWX = zeros(counter + nlam, ncols)
        NEWY = zeros(counter + nlam)
        
        counter = 0
        energy_counter = []
        for calc = 1:ncalc

            if calc == leave_out #skip this one
                continue
            end
            if ERROR[calc] == 1
                continue
            end

            row1, rowN, nw = ind_BIG[calc, 1:3]
            N = hermetian_indpt(nw)

            vals = ones(nw) * 100.0
            #            H = zeros(Complex{Float64}, nw, nw, ncols)
            #            H = zeros(Complex{Float64}, nw, nw)
            H_cols = zeros(Complex{Float64}, nw, nw, ncols)
            
            H_fixed = zeros(Complex{Float64}, nw, nw)
            
            
            #            H_cols = H_COLS[calc]
            
            VECTS = zeros(Complex{Float64}, nw, nw)
            #            S = zeros(Complex{Float64}, nw, nw)
            nk = size(KPOINTS[calc])[1]
            
            X_TOTEN = zeros(ncols)
            Y_TOTEN = ENERGIES[calc]
            
            if scf
                nat = list_of_tbcs[calc].crys.nat
                energy_charge, pot = ewald_energy(list_of_tbcs[calc], DQ[calc,1:nat])
            else
                energy_charge = 0.0
            end
            if list_of_tbcs[calc].tb.scfspin
                energy_magnetic = magnetic_energy(list_of_tbcs[calc], EDEN_FITTED[calc, :, 1:nw])
            else
                energy_magnetic = 0.0
            end
            
            etypes = types_energy(list_of_tbcs[calc].crys)
            
            energy_smear = ENERGY_SMEAR[calc]

            Y_TOTEN -= energy_charge + etypes + energy_smear + energy_magnetic
            
            VECTS = zeros(Complex{Float64}, nw, nw)
            VECTS_p = zeros(Complex{Float64}, nw, nw)

            vals_test_other = zeros(nw, ncols)
            vals_test_on = zeros(nw)

            for spin = 1:list_of_tbcs[calc].nspin
                for k = 1:nk
                    for i = 1:nw
                        for j = 1:nw
                            if i <= j
                                #                            H_cols[i,j,:] = X_Hnew_BIG[row1-1 + (k-1)*nw*nw + (i-1)*nw + j,:] + im*X_Hnew_BIG[row1-1 + (k-1)*nw*nw + (i-1)*nw + j + nk*nw*nw,:]
                                ind = hermetian_index(i,j, nw)
                                H_cols[i,j,:] = X_Hnew_BIG[row1-1 + (k-1)*N + ind,:] + im*X_Hnew_BIG[row1-1 + (k-1)*N + ind +  nk*N,:]
                                H_cols[j,i,:] = X_Hnew_BIG[row1-1 + (k-1)*N + ind,:] - im*X_Hnew_BIG[row1-1 + (k-1)*N + ind +  nk*N,:]                            
                            end
                            
                        end
                    end
                    for i = 1:nw
                        for j = 1:nw
                            if keep_bool
                                if i <= j
                                    ind = hermetian_index(i,j, nw)
                                    H_fixed[i,j] = Xc_Hnew_BIG[row1-1 + (k-1)*N + ind] + im*Xc_Hnew_BIG[row1-1 + (k-1)*N + ind +  nk*N]
                                    H_fixed[j,i] = Xc_Hnew_BIG[row1-1 + (k-1)*N + ind] - im*Xc_Hnew_BIG[row1-1 + (k-1)*N + ind +  nk*N]
                                    
                                end
                                #                            H_fixed[i,j] = Xc_keep[row1-1 + (k-1)*nw*nw + (i-1)*nw + j] + im*Xc_keep[row1-1 + (k-1)*nw*nw + (i-1)*nw + j + nk*nw*nw]
                            end
                        end
                    end


                    #                VECTS[:,:] = VECTS_FITTED[calc,k,1:nw,1:nw]
                    VECTS[:,:] = VECTS_FITTED[calc][spin,k,1:nw,1:nw]
                    VECTS_p[:,:] = VECTS'

                    if keep_bool
                        vals_test_on[:] = real(diag(VECTS_p * (h_on[calc] + H_fixed) * VECTS)) 
                    else
                        vals_test_on[:] = real(diag(VECTS_p * (h_on[calc] ) * VECTS)) 
                    end
                    
                    vals_test_other[:,:] .= 0.0

                    #                @time for i = 1:nw
                    #                    for j = 1:nw
                    #                        for kk = 1:nw
                    #                            for ii = 1:ncols
                    #                                vals_test_other[i,ii] += real(VECTS_p[i,j] .* H_cols[ j,kk,ii]  .* VECTS[kk,i])
                    #                            end
                    #                        end
                    #                    end
                    #                end

                    for ii = 1:ncols
                        for i = 1:nw
                            temp = 0.0+0.0im
                            for kk = 1:nw
                                for j = 1:nw
                                    temp += VECTS_p[i,j] * H_cols[ j,kk,ii]  * VECTS[kk,i]
                                end
                            end
                            vals_test_other[i,ii] += real(temp)
                        end
                    end
                    
                    for i = 1:nw
                        counter += 1


                        NEWX[counter, :] = vals_test_other[i,:] .* WEIGHTS[calc, k, i, spin]
                        X_TOTEN[:] +=   vals_test_other[i,:] .* (KWEIGHTS[calc][k] * OCCS_FITTED[calc,k,i, spin]) #* list_of_tbcs[calc].nspin

                        NEWY[counter] =  (VALS0[calc,k,i, spin] - vals_test_on[i]) .* WEIGHTS[calc, k, i, spin]
                        Y_TOTEN += -1.0 * vals_test_on[i] * (KWEIGHTS[calc][k] * OCCS_FITTED[calc,k,i, spin]) #*  list_of_tbcs[calc].nspin


                        push!(nonzero_ind, counter)

                        ###println("$calc $k $i : ",  vals_test[i], "\t" , VALS_FITTED[calc,k,i],"\t", vals_test_on[i] + vals_test_other[i,:]'*ch)
                    end
                    
                    
                end
            end
            counter += 1
            NEWX[counter, :] = X_TOTEN[:] * energy_weight * weights_list[calc]
            NEWY[counter] = Y_TOTEN * energy_weight * weights_list[calc] 
            push!(nonzero_ind, counter)
            push!(energy_counter, counter)
        end

        if lambda > 1e-10
            for ind3 = 1:nlam
                counter += 1
                NEWX[counter,ind3] = lambda
                push!(nonzero_ind, counter)

            end
        end

        #        println("len nonzero_ind ", length(nonzero_ind))
        #        println("size NEWX old ", size(NEWX))
        #        NEWX = NEWX[nonzero_ind,:]
        #        NEWY = NEWY[nonzero_ind]
        #        println("size NEWX new ", size(NEWX))
        #        println("size NEWY new ", size(NEWY))

        return NEWX, NEWY, energy_counter

    end

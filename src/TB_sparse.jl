
using SparseArrays


mutable struct tb_sparse{T}

    #    H::Array{Complex{Float64},3}
    H::Vector{SparseMatrixCSC{T, Int64}}
    ind_arr::Array{Int64,2}
    r_dict::Dict
    nwan::Int64
    nr::Int64
    nspin::Int64              
    nonorth::Bool
    #    S::Array{Complex{Float64},3}
    S::Vector{SparseMatrixCSC{Float64, Int64}}
    scf::Bool
    scfspin::Bool
    h1::Array{T,2} #scf term
    h1spin::Array{T,3} #scf term
end

Base.show(io::IO, h::tb_sparse) = begin
    println(io)
    nwan=h.nwan
    nr=h.nr
    nonorth = h.nonorth
    scf = h.scf
    scfspin = h.scfspin
    nspin=h.nspin
    println(io, "tight binding real space object (SPARSE); nwan = $nwan, nr = $nr, nonorth = $nonorth, scf = $scf, scfmagnetic = $scfspin, nspin = $nspin" )
    println(io)
    
end   

"""
        mutable struct tb_crys{T}

    Main tight-binding object, holds the tight-binding model `tb` and information about the `crystal`

    # Holds
    - `tb::tb` Has the key tb info (see above)
    - `crys::crystal` Has the crystal structure
    - `nelec::Float64` Number of electrons
    - `dftenergy::Float64` DFT energy for reference, only for fit to DFT cases.
    - `scf::Bool`  `true` if requires self-consistency.
    - `gamma::Array{T, 2}` has the Ewald calculation results, needed for self-consistency.
    - `eden::Array{Float64,2}` electron density, by orbital, if calculated by self-consistency.
    - `within_fit::Bool` is `true` if model is passes tests of being within the fitting parameter space, `false` for extrapolation
    - `energy::Float64` energy in Ryd, if calculated.
    - `efermi::Float64` Fermi energy in Ryd, if calculated.
    - `nspin::Int64` number of spins (2=magnetic)
    """
mutable struct tb_crys_sparse{T}

    tb::tb_sparse
    crys::crystal
    nelec::Float64
    dftenergy::Float64
    scf::Bool
    gamma::Array{T, 2}
    background_charge_correction::T
    eden::Array{Float64,2}
    within_fit::Bool
    energy::Float64
    efermi::Float64
    nspin::Int64
    tot_charge::Float64
    dq::Array{Float64,1}
end



Base.show(io::IO, x::tb_crys_sparse) = begin
    println(io)
    println(io, "tb_crys object (SPARSE)")
    println(io)    
    println(io, x.crys)
    println(io)
    println(io, "nelec: ", x.nelec, "; nspin (hoppings): ", x.nspin)
    ind2orb, orb2ind, etotal, nval = orbital_index(x.crys)
    if abs(nval - x.nelec) > 1e-10
        println(io, "tot_charge: $(nval-x.nelec)")
    end
    println(io, "within_fit: ", x.within_fit,"  ; scf: ", x.scf, "; scfspin: ", x.tb.scfspin)
    println(io, "calculated energy: ", round(convert_energy(x.energy)*1000)/1000, " $global_energy_units")
    println(io, "formation energy: ", round(convert_energy(get_formation_energy(x.energy, x.crys)), digits=3), " $global_energy_units")
    println(io, "efermi  : ", round(convert_energy(x.efermi)*1000)/1000, " $global_energy_units")
    dq = get_dq(x)
    println(io, "charges : ", round.(dq * 100)/100)
    if size(x.eden)[1] == 2
        mm = get_magmom(x)
        println(io, "mag mom.: ", round.(mm * 100)/100)
    end
    println(io)
    println(io, x.tb)    
    println(io)
    
end   

"""
        function make_tb(H, ind_arr, r_dict::Dict, S; h1=missing)

    Constructor function for `tb` with overlaps
    """
function make_tb_sparse(H, ind_arr, r_dict::Dict, S; h1=missing, h1spin = missing)
    nw=size(H[1],2)
    nr=length(H)
    nspin=1

    if size(ind_arr) != (nr,3) 
        error("make_tb ind_arr size ", size(ind_arr), size(H))
    end

    if length(S) != nr || size(S[1])[1] != nw
        error("make_tb size S doesn't match size H", size(S), size(H))
    end
    
    T=typeof(real(H[1][1,1]))

    if ismissing(h1) 
        h1 =zeros(T, nw, nw)
        scf = false
    else
        scf = true
    end

    if ismissing(h1spin)
        h1spin = zeros(T, 2,nw,nw)
        scfspin = false
    else
        scfspin = true
    end

    return tb_sparse{T}(H, ind_arr,  r_dict,nw, nr, nspin, true, S, scf,scfspin, h1, h1spin)
end

function make_tb_crys_sparse(ham::tb_sparse,crys::crystal, nelec::Float64, dftenergy::Float64; scf=false, eden = missing, gamma=missing, background_charge_correction=0.0, within_fit=true, screening=1.0, tb_energy=-999, fermi_energy=0.0 )

    T = typeof(crys.coords[1,1])
    nspin = ham.nspin
    if ismissing(eden)
#        if scf == false
#            eden = zeros(nspin,ham.nwan)
#        else
            eden = get_neutral_eden(crys, ham.nwan, nspin=nspin)
            bv = eden .> 1e-5
#            println("eden $eden sum $(sum(eden)) nelec $nelec")
            eden[bv] = eden[bv] .-  (sum(eden) - nelec / 2.0)/crys.nat
#            println("new ", eden)
 #       end
#        println("start eden ", eden)
    end

    dq = get_dq(crys, eden)
    tot_charge = -sum(dq)
    
    
    println("gamma")
    @time if ismissing(gamma) 
        #        println("ismissing gamma")
        gamma, background_charge_correction = electrostatics_getgamma(crys, screening=screening) #do this once and for all
    end

    nspin = ham.nspin
    return tb_crys_sparse{T}(ham,crys,nelec, dftenergy, scf, gamma, background_charge_correction, eden, within_fit, tb_energy, fermi_energy, nspin, tot_charge, dq)
end

function calc_energy_charge_fft_band2_sym_sparse(hk3, sk3, nelec; smearing=0.01, h1 = missing, h1spin=missing, VECTS=missing, DEN=missing, nk_red=nk_red, grid_ind=[1 1 1], kweights = [2.0] )

    #println("begin")
    begin
        thetype=typeof(real(sk3[1][1,1]))

        #     println("size ", size(hk3))
        nwan = size(sk3[1])[1]
        nspin = 1

        if ismissing(VECTS)
            VECTS = zeros(Complex{thetype}, nwan, nwan, nk_red, nspin)
        end

        if ismissing(DEN)
            DEN = zeros(Complex{thetype}, nwan, nwan, nk_red)
        end

        rDEN = zeros(thetype, nwan, nwan, nk_red)
        iDEN = zeros(thetype, nwan, nwan, nk_red)
        rv = zeros(thetype, nwan, nwan, nk_red)
        iv = zeros(thetype, nwan, nwan, nk_red)

#        println("h1 $h1")
        HK = zeros(Complex{thetype},nspin, size(h1)[1], size(h1)[1],nk_red)
        
        
        if true

            
            #grid = size(sk3)[3:5]
            #    print("calc_energy_charge_fft_band grid $grid")
            #nk = prod(grid)
            #nwan = size(hk3)[1]

            if !ismissing(h1spin)
                nspin = 2
            else
                nspin = 1
            end

            nspin_ham = 1


            VALS = zeros(Float64, nk_red, nwan, nspin)
            VALS0 = zeros(Float64, nk_red,nwan, nspin)
            #         c=0


            #         sk = zeros(Complex{thetype}, nwan, nwan)
            #         hk = zeros(Complex{thetype}, nwan, nwan)
            #         hk0 = zeros(Complex{thetype}, nwan, nwan)


            #         VECTS = zeros(Complex{thetype}, nk, nspin, nwan, nwan)
            #         SK = zeros(Complex{thetype}, nk, nwan, nwan)

            error_flag = false

            if ismissing(h1)
                h1 = zeros(nwan,nwan)
            else
                h1 = 0.5*(h1 + h1')
            end
            if ismissing(h1spin)
                h1spin = zeros(2,nwan,nwan)# , zeros(nwan,nwan)]
            else
                h1spin[1,:,:] .= 0.5*(h1spin[1,:,:] + h1spin[1,:,:]')
                h1spin[2,:,:] .= 0.5*(h1spin[2,:,:] + h1spin[2,:,:]')
            end
            
        end
    end


    println("go eig time")
    @time go_eig_sym_sparse(grid, nspin,nspin_ham, VALS, VALS0, VECTS, sk3, hk3, h1, h1spin, nk_red,grid_ind)

#    println("VALS ", VALS)
#    println("VALS0 ", VALS0)
    
    begin

        if nelec > 1e-10
#            println("VALS ", VALS[1,:,1])
            energy, efermi = band_energy(VALS, kweights, nelec, smearing, returnef=true)
            occ = gaussian.(VALS.-efermi, smearing)

#            println("nelec $nelec efermi $efermi sum(occ) $(sum(occ .* kweights)/2.0)   sym")
            
            max_occ = findlast(sum(occ, dims=[1,3]) .> 1e-8)[2]
            energy_smear = smearing_energy(VALS, kweights, efermi, smearing)

            energy0 = sum(occ .* VALS0 .* kweights) 

            energy0 += energy_smear * nspin#
        else

            energy_smear = 0.0
            energy0 = 0.0
            efermi = minimum(VALS)
            occ = zeros(size(VALS))
            
        end
        
            
        #         println("energy_smear , ", energy_smear * nspin, " energy0 ", sum(occ .* VALS0) / nk * 2.0)
        
    end

    #println("nelec $nelec")

    println("chargeden")
    @time if nelec > 1e-10
        chargeden = go_charge15_sym_sparse(VECTS, sk3, occ, nspin, max_occ, rDEN, iDEN, rv, iv, nk_red,grid_ind, kweights)
    else
        chargeden = zeros(nspin, nwan)
    end
    
    if nspin == 2
        energy0 = energy0 / 2.0
    end
    
    return energy0, efermi, chargeden, VECTS, VALS, error_flag #, denmat, HK


end 

function go_eig_sym_sparse(grid, nspin, nspin_ham, VALS, VALS0, VECTS, sk3, hk3, h1, h1spin, nk_red, grid_ind)

    #    max_num = nthreads()
    max_num = 1
    println("assemble memory")
    @time begin 
        hk = zeros(Complex{Float64}, size(h1)[1], size(h1)[1], max_num)
        sk = zeros(Complex{Float64}, size(h1)[1], size(h1)[1], max_num)
        #    hk0 = zeros(Complex{Float64}, size(h1)[1], size(h1)[1], max_num)
        vals = zeros(Complex{Float64}, size(h1)[1], max_num)
        vects = zeros(Complex{Float64}, size(h1)[1], size(h1)[1], max_num)
    end
    
#    hermH = Hermitian(zeros(Complex{Float64}, size(h1)[1], size(h1)[1]))
    #    hermS = Hermitian(zeros(Complex{Float64}, size(h1)[1], size(h1)[1]))

    
    
    @inbounds @fastmath for c = 1:nk_red
        id = threadid()
        
        sk[:,:,id] .= collect(sk3[c]) 
        for spin = 1:nspin
            spin_ind = min(spin, nspin_ham)
            
            #            hk0[:,:,id] .= ( (@view hk3[:,:,spin_ind, k1,k2,k3]) )
            #            hk0[:,:,id] = 0.5*(hk0[:,:,id]+hk0[:,:,id]')
            #            hk = hk0  .+ 0.5*sk .* (h1 + h1spin[spin,:,:] + h1' + h1spin[spin,:,:]')
            
            hk[:,:, id] .= hk3[c]  .+ sk[:,:,id] .* (h1 + (@view h1spin[spin,:,:] ))

#            HK[spin, :,:,c] = hk[:,:, id]
            #hk[:,:,id] .= 0.5*( (@view hk[:,:,id]) .+ (@view hk[:,:,id])')
            
            try
                #hermH[:,:] = (@view hk[:,:,id][:,:])
                #hermS[:,:] = (@view sk[:,:,id][:,:])
                @time vals[:,id], vects[:,:,id] = eigen( Hermitian(@view hk[:,:,id][:,:]), Hermitian(@view sk[:,:,id][:,:]))
#                if c == 1
#                    println()
#                    println("hk ")
#                    println(hk[:,:,id][:,:])
#                    println()
#                    println("vals ", vals[:,id])
#                    println()
#                end
                #vals[:,id], vects[:,:,id] = eigen( hermH, hermS)
            catch err
                typeof(err) == InterruptException && rethrow(err)
                vals[:,id], vects[:,:,id] = eigen( hk[:,:,id][:,:], sk[:,:,id][:,:])
            end
            
            if maximum(abs.(imag.(vals))) > 1e-10

                println("s ", eigvals(sk[:,:,id])[:,:])
                error_flag = true
            end
            
            VALS[c,:, spin] .= real.(vals[:,id])
            #            VALS0[c,:, spin] .= real.(diag(vects'*hk0[:,:,id]*vects))
            println("VALS0")
            @time VALS0[c,:, spin] .= real.(diag( (@view vects[:,:,id])'*(hk3[c])*(@view vects[:,:,id])))
            
#            if c == 1
#                println("vals0 ", VALS0[c,:, spin])
#                println()
#            end
            VECTS[:,:, c, spin] .= (@view vects[:,:,id])
            
        end
    end

end

function  myfft_R_to_K_sparse(tbc::tb_crys_sparse, nk_red, kpts)

    T = typeof(tbc.tb.H[1][1,1])
               
    Hk = SparseMatrixCSC{Complex{T}, Int64}[]
    Sk = SparseMatrixCSC{Complex{T}, Int64}[]

    twopi_i = 2*pi*(1.0 * im)
    
    for k in 1:nk_red

        push!(Hk, spzeros(Complex{T}, tbc.tb.nwan, tbc.tb.nwan))
        push!(Sk, spzeros(Complex{T}, tbc.tb.nwan, tbc.tb.nwan))
        
        kpoint = kpts[k,:]
        kmat = repeat(kpoint', tbc.tb.nr,1)
        exp_ikr = exp.(twopi_i * sum(kmat .* tbc.tb.ind_arr,dims=2))
        
        for nr = 1:tbc.tb.nr
            Hk[end] += tbc.tb.H[nr] * exp_ikr[nr]
            Sk[end] += tbc.tb.S[nr] * exp_ikr[nr]
        end
    end
    return  Hk, Sk


end

function get_dq(tbc::tb_crys_sparse)
    return get_dq(tbc.crys, tbc.eden)
end

function get_h1(tbc::tb_crys_sparse, chargeden::Array{Float64,2})
    dq = get_dq(tbc.crys, chargeden)
    h1 = get_h1_dq(tbc,dq)
    return h1, dq
end

function ewald_energy(tbc::tb_crys_sparse, delta_q=missing)

    background_charge_correction = tbc.background_charge_correction
    gamma = tbc.gamma 
    crys = tbc.crys

    if ismissing(delta_q)
        delta_q =  get_dq(crys , tbc.eden)
    end

    return ewald_energy(crys, gamma, background_charge_correction, delta_q)

end

function go_charge15_sym_sparse(VECTS, S, occ, nspin, max_occ, rDEN, iDEN, rv, iv,  nk_red, grid_ind,kweights )

    nw = size(S[1])[1]
    #    nk = size(S)[3]

    d = zeros(Complex{Float64}, nw,nw)
    charge = zeros(nspin, nw)

#    denmat = zeros(Complex{Float64}, nspin, nw, nw, nk_red)

#    println("size kw ", size(kweights))
#    println("size occ ", size(occ))
#    println("size rv ", size(rv))
#    println("nw $nw nk_red $nk_red max_occ $max_occ")
    for spin = 1:nspin
        rv[:,:,:] .= real.(VECTS[:,:,:,spin])
        iv[:,:,:] .= imag.(VECTS[:,:,:,spin])    

        rDEN .= 0.0
        iDEN .= 0.0
        @tturbo for n = 1:max_occ  #tturbo
            for k = 1:nk_red
                for b = 1:nw
                    for a = 1:nw
                        #DEN[a,b,k] += occ[k,n,spin].*conj(VECTS[a,n,k,spin]).*(VECTS[b,n,k,spin])

                        #                        DEN[a,b,k] += occ[k,n,spin].*  (real(VECTS[a,n,k,spin]) - im * imag(VECTS[a,n,k,spin])) .*(real(VECTS[b,n,k,spin]) + im * imag(VECTS[b,n,k,spin]))

#                        DEN[a,b,k] += occ[k,n,spin].*  ( real(VECTS[a,n,k,spin])*real(VECTS[b,n,k,spin]) + real(VECTS[a,n,k,spin])*im * imag(VECTS[b,n,k,spin]) + (-im)*imag(VECTS[a,n,k,spin])*real(VECTS[b,n,k,spin]) + (-im)*imag(VECTS[a,n,k,spin])*im*imag(VECTS[b,n,k,spin]))

                        rDEN[a,b,k] += kweights[k]*occ[k,n,spin]*(rv[a,n,k]*rv[b,n,k] + iv[a,n,k]*iv[b,n,k])
                        iDEN[a,b,k] += kweights[k]*occ[k,n,spin]*(rv[a,n,k]*iv[b,n,k] - iv[a,n,k]*rv[b,n,k])
                        
                        #DEN[a,b,k] += occ[k,n,spin]*(rv[a,n,k,spin]*rv[b,n,k,spin]-iv[b,n,k,spin]*iv[a,n,k,spin])

                    end
                end
            end
        end
#        for k = 1:nk_red
#            #            denmat[spin,:,:,k] = rDEN[:,:,k] + im*iDEN[:,:,k]
#            for n = 1:max_occ
#                denmat[spin,:,:,k] += occ[k,n,spin]*VECTS[:,n,k,spin] * VECTS[:,n,k,spin]'
#            end
#        end

        
#        println("size rDEN $(size(rDEN)) i $(size(iDEN)) s $(size(S)) d $(size(d))")
        d .= 0.0
        for k = 1:nk_red
            d += sum(((@view rDEN[:,:,k])+(@view iDEN[:,:,k])*im) .* S[k], dims=3)[:,:]
        end
        #        d .= sum(DEN .* S, dims=3)[:,:]
        charge[spin,:] = sum( real(0.5*(d + d')), dims=1)
    end
    return charge/2.0 #, denmat/2.0
end

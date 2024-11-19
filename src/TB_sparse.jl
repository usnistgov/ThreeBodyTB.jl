
using SparseArrays

#Contains sparse matrix implementations of severals structs and functions in parallel to the main TB.jl

"""
    mutable struct tb_sparse{T}

Holds the tb object, but with a sparse matrix implementation (see also `tb`). `H` and `S` are vectors that contain sparce matricies, as opposed to a higher-dimensional dense matrix that is implemented in the dense version.
"""
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
        mutable struct tb_crys_sparse{T}

    Main tight-binding object (SPARSE matrix version, see `tb_crys_dense`), holds the tight-binding model `tb` and information about the `crystal`

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
mutable struct tb_crys_sparse{T} <: tb_crys

    tb::tb_sparse
    crys::crystal
    nelec::Float64
    dftenergy::Float64
    scf::Bool
    gamma::Array{T, 2}
    u3::Array{T, 1}
    background_charge_correction::T
    eden::Array{Float64,2}
    within_fit::Bool
    energy::Float64
    efermi::Float64
    nspin::Int64
    tot_charge::Float64
    dq::Array{Float64,1}
    energy_band::Float64
    energy_smear::Float64
    energy_types::Float64
    energy_charge::Float64
    energy_mag::Float64
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
    dq, dq_eden = get_dq(x)
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
        function make_tb_sparse(H, ind_arr, r_dict::Dict, S; h1=missing)

Constructor function for `tb_sparse` with overlaps
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

"""
function make_tb_crys_sparse

Constructor function for `tb_crys_sparse`.
"""
function make_tb_crys_sparse(ham::tb_sparse,crys::crystal, nelec::Float64, dftenergy::Float64; scf=false, eden = missing, gamma=missing,u3=missing, background_charge_correction=0.0, within_fit=true, screening=1.0, tb_energy=-999, fermi_energy=0.0, energy_band= 0.0, energy_smear = 0.0, energy_types = 0.0, energy_charge = 0.0, energy_mag = 0.0 )

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

    dq, dq_eden = get_dq(crys, eden)
    tot_charge = -sum(dq)
    
    
    #println("gamma")
    if ismissing(gamma) 
        #        println("ismissing gamma")
        gamma, background_charge_correction,u3 = electrostatics_getgamma(crys, screening=screening) #do this once and for all
    end
    
    energy_types = types_energy(crys)
    nspin = ham.nspin
    println()
    println("start")
    println(ham)
    println(crys)
    println(nelec)
    println(dftenergy)
    println(scf)
    println(gamma)
    println(u3)
    println(background_charge_correction)
    println(eden)
    println(within_fit)
    println(tb_energy)
    println(fermi_energy)
    println(nspin)
    println(tot_charge)
    println(dq)
    println([energy_band, energy_smear , energy_types , energy_charge , energy_mag ])

    return tb_crys_sparse{T}(ham,crys,nelec, dftenergy, scf, gamma,u3, background_charge_correction, eden, within_fit, tb_energy, fermi_energy, nspin, tot_charge, dq, energy_band, energy_smear , energy_types , energy_charge , energy_mag )
end

function calc_energy_charge_fft_band2_sym_sparse(hk3, sk3, nelec; smearing=0.01, h1 = missing, h1spin=missing, VECTS=missing, DEN=missing, nk_red=nk_red, kweights = [2.0],SI=[], SJ=[], rSV=[], iSV=[], maxS=0 )

    if length(SI) == 0
        SI = []
        SJ = []
        rSV = []
        iSV = []
        maxS = 0
        for s in sk3
            I,J,V = findnz(s)
            push!(SI, I)
            push!(SJ, J)
            push!(rSV, real(V))
            push!(iSV, imag(V))
            maxS = max(maxS, length(V))
        end
    end
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


    #println("go eig time")
    go_eig_sym_sparse(grid, nspin,nspin_ham, VALS, VALS0, VECTS, sk3, hk3, h1, h1spin, nk_red)

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

#    println("chargeden")
    if nelec > 1e-10
#        @time chargeden = go_charge15_sym_sparse(VECTS, sk3, occ, nspin, max_occ, rDEN, iDEN, rv, iv, nk_red,grid_ind, kweights)
#        println("chargeden  ", sum(chargeden))
#        @time chargeden_opt1 = go_charge15_sym_sparse_opt(VECTS, sk3, occ, nspin, max_occ, rDEN, iDEN, rv, iv, nk_red,grid_ind, kweights)
#        println("chargedenO1 ", sum(chargeden_opt1))
        chargeden = go_charge15_sym_sparse_opt2(VECTS, sk3, occ, nspin, max_occ, rDEN, iDEN, rv, iv, nk_red,kweights, SI, SJ, rSV, iSV, maxS)
        
#        println("chargedenO2 ", sum(chargeden_opt2))
#        println("dchargeden ", sum(abs.(chargeden_opt1 - chargeden)))
#        println("dchargeden2 ", sum(abs.(chargeden_opt2 - chargeden)))
    else
        chargeden = zeros(nspin, nwan)
    end
    
    if nspin == 2
        energy0 = energy0 / 2.0
    end
    
    return energy0, efermi, chargeden, VECTS, VALS, error_flag #, denmat, HK


end 

"""
    function go_eig_sym_sparse

Helper function for calculating eigenvalues/vectors for sparse matricies. Uses dense linear algebra for the actual diagonalization.
Need to implement a linear-scaling version.
"""
function go_eig_sym_sparse(grid, nspin, nspin_ham, VALS, VALS0, VECTS, sk3, hk3, h1, h1spin, nk_red)

    #    max_num = nthreads()

    #println("assemble memory")
    begin 
        hk = zeros(Complex{Float64}, size(h1)[1], size(h1)[1])
        sk = zeros(Complex{Float64}, size(h1)[1], size(h1)[1])
        #    hk0 = zeros(Complex{Float64}, size(h1)[1], size(h1)[1], max_num)
        vals = zeros(Complex{Float64}, size(h1)[1])
        vects = zeros(Complex{Float64}, size(h1)[1], size(h1)[1])
    end
    
#    hermH = Hermitian(zeros(Complex{Float64}, size(h1)[1], size(h1)[1]))
    #    hermS = Hermitian(zeros(Complex{Float64}, size(h1)[1], size(h1)[1]))

    
    
    @inbounds @fastmath for c = 1:nk_red
        sk .= collect(sk3[c]) 
        for spin = 1:nspin
            spin_ind = min(spin, nspin_ham)
            
            #            hk0[:,:,id] .= ( (@view hk3[:,:,spin_ind, k1,k2,k3]) )
            #            hk0[:,:,id] = 0.5*(hk0[:,:,id]+hk0[:,:,id]')
            #            hk = hk0  .+ 0.5*sk .* (h1 + h1spin[spin,:,:] + h1' + h1spin[spin,:,:]')
            
            hk .= hk3[c]  .+ sk .* (h1 + (@view h1spin[spin,:,:] ))

#            HK[spin, :,:,c] = hk[:,:, id]
            #hk[:,:,id] .= 0.5*( (@view hk[:,:,id]) .+ (@view hk[:,:,id])')
            
            try
                #hermH[:,:] = (@view hk[:,:,id][:,:])
                #hermS[:,:] = (@view sk[:,:,id][:,:])
                #println("eig")
                vals, vects = eigen( Hermitian(hk), Hermitian(sk))
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
                println("eig failed")
                typeof(err) == InterruptException && rethrow(err)
                vals, vects = eigen(Hermitian( hk),Hermitian( sk))
            end
            
            if maximum(abs.(imag.(vals))) > 1e-10
                
                println("s ", eigvals(sk))
                error_flag = true
            end
            
            VALS[c,:, spin] .= real.(vals)
            #            VALS0[c,:, spin] .= real.(diag(vects'*hk0[:,:,id]*vects))
#            println("VALS0")
            VALS0[c,:, spin] .= real.(diag( ( vects)'*(hk3[c])*(vects)))
            
#            if c == 1
#                println("vals0 ", VALS0[c,:, spin])
#                println()
#            end
            VECTS[:,:, c, spin] .= vects
            
        end
    end

end

"""
    function  myfft_R_to_K_sparse(tbc::tb_crys_sparse, nk_red, kpts)

Sparse matrix fourier transform of Hamiltonian. Not actually an fft,
just a naive fourier transform. Becomes fast if the matricies are
actually sparse, so that matrix multiplication becomes linear in the
number of atoms instead of `nat^3`. In that case we only need a few
k-points anyways, so the fft is not needed.  
"""
function myfft_R_to_K_sparse(tbc::tb_crys_sparse, nk_red, kpts)

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



"""
    function go_charge15_sym_sparse_opt2

Helper function for calculating charge density using sparse matrix multiplcation. Becomes faster if matricies are actually sparse (large unit cells or low dimensions).
"""
function go_charge15_sym_sparse_opt2(VECTS, S, occ, nspin, max_occ, rDEN, iDEN, rp, ip,  nk_red, kweights, SI, SJ, rSV, iSV, maxS )

    nw = size(S[1])[1]
    #    nk = size(S)[3]
    #d = zeros(Complex{Float64}, nw,nw)
    rd = zeros(Float64, maxS)
    id = zeros(Float64, maxS)
    charge = zeros(nspin, nw)


    
    for spin = 1:nspin

        rp[:,:,:] .= real.(VECTS[:,:,:,spin])
        ip[:,:,:] .= imag.(VECTS[:,:,:,spin])    

        d = spzeros(Complex{Float64}, nw,nw)
        
        for k = 1:nk_red
            kw =kweights[k]
            #            I,J,V = findnz(S[k])
            I = SI[k]
            J = SJ[k]
            rV = rSV[k]
            iV = iSV[k]
            lv = length(rSV[k])
            #            rV = real(V)
            #            iV = imag(V)

            rd .= 0.0
            id .= 0.0

            @tturbo for c in 1:lv
                i = I[c]
                j = J[c]
                rv = rV[c]
                iv = iV[c]
                for a = 1:max_occ
                    #                    d[i,j] += conj(VECTS[i,a,k,spin])*VECTS[j,a,k,spin]*occ[k,a,spin] * kw * v
                    #                    rd[i,j] += (rv[i,a,k]*rv[j,a,k] + iv[i,a,k]*iv[j,a,k])*occ[k,a,spin] * kw * rV[c]
                    #                    id[i,j] += (rv[i,a,k]*iv[j,a,k] - iv[i,a,k]*rv[j,a,k])*occ[k,a,spin] * kw * rV[c]
                    #                    rd[i,j] += occ[k,a,spin] * kw*(-ip[j,a,k]*iv*rp[i,a,k] + ip[i,a,k]* ip[j,a,k]*rv + rp[i,a,k]*rp[j,a,k]*rv + ip[i,a,k]*ip[j,a,k]*iv)
                    #                    id[i,j] += occ[k,a,spin] * kw*(ip[i,a,k]* ip[j,a,k]*iv + iv*rp[i,a,k]*rp[j,a,k] + ip[j,a,k]*rp[i,a,k]*rv -  ip[i,a,k]*rp[j,a,k]*rv )


                    #                    rd[i,j] +=  occ[k,a,spin] * kw*( -ip[j,a,k]*iv*rp[i,a,k] + ip[i,a,k]*iv*rp[j,a,k] +  ip[i,a,k]*ip[j,a,k]*rv + rp[i,a,k]*rp[j,a,k]*rv)
                    #                    id[i,j] +=  occ[k,a,spin] * kw*(  ip[i,a,k]*ip[j,a,k]*iv + iv*rp[i,a,k]*rp[j,a,k] -  ip[i,a,k]*rp[j,a,k]*rv + ip[j,a,k]*rp[i,a,k]*rv)

#                    rd[i,j] +=  occ[k,a,spin] * kw*( -ip[j,a,k]*iv*rp[i,a,k] + ip[i,a,k]*iv*rp[j,a,k] +  ip[i,a,k]*ip[j,a,k]*rv + rp[i,a,k]*rp[j,a,k]*rv)
#                    id[i,j] +=  occ[k,a,spin] * kw*(  ip[i,a,k]*ip[j,a,k]*iv + iv*rp[i,a,k]*rp[j,a,k] -  ip[i,a,k]*rp[j,a,k]*rv + ip[j,a,k]*rp[i,a,k]*rv)

                    rd[c] +=  occ[k,a,spin] * kw*( -ip[j,a,k]*iv*rp[i,a,k] + ip[i,a,k]*iv*rp[j,a,k] +  ip[i,a,k]*ip[j,a,k]*rv + rp[i,a,k]*rp[j,a,k]*rv)
                    id[c] +=  occ[k,a,spin] * kw*(  ip[i,a,k]*ip[j,a,k]*iv + iv*rp[i,a,k]*rp[j,a,k] -  ip[i,a,k]*rp[j,a,k]*rv + ip[j,a,k]*rp[i,a,k]*rv)
                    
                    
                end
            end
            d += sparse(I, J,(@view  (rd+im*id)[1:lv]))
        end

        charge[spin,:] = 0.5*sum( real( d + d' ) , dims=1)
    end

    return charge/2.0  #, denmat/2.0
end

#sparse matrix version of arbitrary k-point diagonalization. Uses dense diagonalization
function Hk(h::tb_sparse, kpoint; spin=1)

#    println("Hk sparse")
    
    kpoint = vec(kpoint)
    hk = zeros(Complex{Float64}, h.nwan, h.nwan)
    sk = zeros(Complex{Float64}, h.nwan, h.nwan)
    hk0 = zeros(Complex{Float64}, size(hk))

    if h.nonorth  
        fill!(sk, zero(Float64))
    end

    twopi_i = -1.0im*2.0*pi

    #    println("repeat")
    kmat = repeat(kpoint', h.nr,1)

    #    println("exp")
    exp_ikr = exp.(twopi_i * sum(kmat .* h.ind_arr,dims=2))

    for m in 1:h.nwan
        for n in 1:h.nwan        
            for k = 1:h.nr
                if h.nspin == 2
                    hk0[m,n] += h.H[k][m,n]'*exp_ikr[k]
                else
                    hk0[m,n] += h.H[k][m,n]'*exp_ikr[k]
                end
                if h.nonorth
                    sk[m,n] += h.S[k][m,n]'*exp_ikr[k]
                end
            end
        end
    end
    hk0 = 0.5*(hk0 + hk0')

    hk .= hk0
    if h.scf
        hk .+= sk .* h.h1
    end
    if h.scfspin
        hk .+= sk .* h.h1spin[spin,:,:]
    end

    hk = 0.5*(hk[:,:] + hk[:,:]')

    
    nw=size(hk)[1]
    vects = zeros(nw,nw)
    vals = zeros(nw)
    vals0 = zeros(nw)

    try
        if h.nonorth
            sk = 0.5*(sk[:,:] + sk[:,:]')
            F=eigen(Hermitian(hk[:,:]),Hermitian( sk[:,:]))
        else
            #        println("orth")
            #        println(typeof(hk))
            hk = 0.5*(hk[:,:] + hk[:,:]')            
            F=eigen(Hermitian(hk[:,:])) #orthogonal
        end

        vects = F.vectors
        vals = real(F.values)
        vals0 = real.(diag(vects'*hk0*vects))

    catch
        println("warning eigen failed, ", kpoint)
        vects = collect(I(nw))
        vals = 1000.0 * ones(nw)
        vals0 = 1000.0 * ones(nw)

    end

    return vects, sort!(vals), hk, sk, vals0

end

#sparse version.
function calc_energy_fft(tbc::tb_crys_sparse; grid=missing, smearing=0.01, return_more_info=false, use_sym=false, scissors_shift = 0.0, scissors_shift_atoms = [])

    etypes = types_energy(tbc.crys)

    if tbc.nspin == 2 || tbc.tb.scfspin == true
        nspin = 2
    else
        nspin = 1
    end


    #     println("size(hk3) ", size(hk3))
    
    #     if tbc.tb.scfspin 
    #         h1 = tbc.tb.h1spin
    #     else
    #         h1 = missing
    #     end

    #     if tbc.nspin == 2 || tbc.tb.scfspin == true
    #         nspin = 2
    #     else
    #         nspin = 1
    #     end

    use_sym = true
    
    nk_red, grid_ind, kpts, kweights = get_kgrid_sym(tbc.crys, grid=grid)

    hk3, sk3 = myfft_R_to_K_sparse(tbc, nk_red, kpts)

    
    if abs(scissors_shift) > 1e-10
        return_more_info = true
    end
    
    ret =  calc_energy_fft_band_sparse(hk3, sk3, tbc.nelec, smearing=smearing, return_more_info=true, h1=tbc.tb.h1, h1spin = tbc.tb.h1spin , nspin=nspin, use_sym=use_sym, nk_red=nk_red, grid_ind = grid_ind, kweights=kweights)

#    println("ret ", ret)
    energy, efermi, vals, vects  = ret

    vals_old = deepcopy(vals)
    vects_old = deepcopy(vects)
    hk3 = deepcopy(hk3)
    sk3 = deepcopy(sk3)
    
    if abs(scissors_shift) > 1e-10
        println("apply scissors shift $scissors_shift (in ryd) to $scissors_shift_atoms")
        energy, efermi, vals, vects  = ret
        vals, vects = apply_scissors_shift(efermi, vals, vects, scissors_shift, scissors_shift_atoms, tbc.crys, hk3, sk3, tbc.nspin, tbc.tb.h1, tbc.tb.h1spin) 
        
    end
    
    
    if return_more_info

        #        println("PRE-precheck ", sum(vects[1,:,:]' * sk3[:,:,1,1,1] * vects[1,:,:]))

        etot = etypes + energy
        return etot, efermi, vals, vects, hk3, sk3

    end
    return energy + etypes

end

function calc_energy_fft_band_sparse(hk3, sk3, nelec; smearing=0.01, return_more_info=false, h1 = missing, h1spin = missing, nspin=1, use_sym=false, nk_red=missing, grid_ind = missing, kweights=missing)

    thetype=typeof(real(sk3[1][1]))
    nwan = size(hk3[1])[1]
    if use_sym
        println("use_sym")
        VALS = zeros(Float64, nk_red,nwan,nspin)
        VECTS = zeros(Complex{Float64}, nk_red,nspin, nwan, nwan)
        
        spin_size=1
        go_sym_sparse(grid, sk3, hk3, h1, h1spin, VALS, VECTS, nk_red, grid_ind, thetype, nwan, nspin, spin_size,return_more_info)
        
        band_en, efermi = band_energy(VALS, kweights, nelec, smearing, returnef=true)
        energy_smear = smearing_energy(VALS, kweights, efermi, smearing)
        
    end
    
    if return_more_info
        return  band_en + energy_smear, efermi, VALS, VECTS
    else
        return  band_en + energy_smear
    end
    
end 

function go_sym_sparse(grid, sk3, hk3, h1, h1spin, VALS, VECTS, nk_red, grid_ind, thetype, nwan, nspin, spin_size, return_more_info)
    c=0
    sk = zeros(Complex{thetype}, nwan, nwan)
    hk = zeros(Complex{thetype}, nwan, nwan)
    for c = 1:nk_red
        
        try

            sk[:,:] = 0.5*(sk3[c] + sk3[c]')
            for spin = 1:nspin
                spins = min(spin, spin_size)
                hk[:,:] = 0.5*(hk3[c]+ hk3[c]')
                
                if !ismissing(h1)
                    hk += 0.5*sk .* (h1+h1')
                end
                if nspin == 2
                    hk += 0.5*sk .* (h1spin[spin,:,:] + h1spin[spin,:,:]')
                end
                
                vals, vects = eigen(Hermitian(hk),Hermitian( sk))
                VALS[c, :,spin] = real(vals)

                if return_more_info
                    VECTS[c,spin, :,:] = vects
                end
            end
            #                    println("tb check ", sum(vects' * sk * vects))
            #                    println("tb check2    ", sum(VECTS[c,:,:]' * sk3[:,:,k1,k2,k3] * VECTS[c,:,:]))

        catch e
            if e isa InterruptException
                println("user interrupt")
                rethrow(InterruptException)
            end

            println("error calc_energy_fft $k1 $k2 $k3 usually due to negative overlap eigenvalue")
            sk[:,:] = 0.5*(sk3[c] + sk3[c]')
            valsS, vectsS = eigen(Hermitian(sk))
            println(valsS)

            rethrow(error("BadOverlap"))

        end
    end
end

"""
         function calc_energy_charge_fft(tbc::tb_crys_sparse; grid=missing, smearing=0.01)

Sparse matrix non-scf energy/eigenvectors.
"""
function calc_energy_charge_fft(tbc::tb_crys_sparse; grid=missing, smearing=0.01)

    #     println("asdf")
    etypes = types_energy(tbc.crys)

    if ismissing(grid)
        grid = get_grid(tbc.crys)
    else
        if length(grid) != 3
            grid = get_grid(tbc.crys, grid) 
        end
    end
    #    println("calc_energy_charge_fft grid $grid")

    nk_red, grid_ind, kpts, kweights = get_kgrid_sym(tbc.crys, grid=grid)
    hk3, sk3 = myfft_R_to_K_sparse(tbc, nk_red, kpts)

    if tbc.scf
        h1 = tbc.tb.h1
        echarge, pot = ewald_energy(tbc)
    else
        h1 = missing
        echarge = 0.0
    end

    if tbc.nspin == 2 || tbc.tb.scfspin
        #         h1up, h1dn = get_spin_h1(tbc)
        #         h1spin = [h1up, h1dn]
        h1spin = get_spin_h1(tbc)
        emag = magnetic_energy(tbc)
    else
        h1spin = missing
        emag = 0.0
    end


    eband, efermi, chargeden, VECTS, VALS, error_flag  =  calc_energy_charge_fft_band2_sym_sparse(hk3, sk3, tbc.nelec, smearing=smearing, h1 = h1, h1spin = h1spin, nk_red = nk_red, kweights=kweights)
    tbc.efermi = efermi
    tbc.eden[:,:] = chargeden[:,:]
    #println("energy comps $eband $etypes $echarge $emag")
    energy = eband + etypes + echarge + emag

    #     println("end asdf")

    return energy, efermi, chargeden, VECTS, VALS, error_flag

end

"""
    function convert_sparse_dense(tbc::tb_crys)

Converts a `tb_crys_dense` to `tb_crys_sparse` or vice versa,
depending on what you give it.  Users usually don't have to worry
about using dense or sparse versions; publically exported commands
perform the same either way, and main functions are overloaded or they
accept the supertype `tb_crys`.  Can be useful for testing speed
difference though.  
"""
function convert_sparse_dense(tbc::tb_crys)

    if typeof(tbc) <: tb_crys_dense
        println("convert to sparse")
        H = []
        S = []
        for r in 1:tbc.tb.nr
            push!(H, sparse( real.(tbc.tb.H[1,:,:,r])))
            push!(S, sparse( real.(tbc.tb.S[:,:,r])))
        end
        tb = make_tb_sparse( H , tbc.tb.ind_arr, tbc.tb.r_dict,  S, h1 = tbc.tb.h1, h1spin = tbc.tb.h1spin)
        tb.scfspin = tbc.tb.scfspin
        tbc_new = make_tb_crys_sparse(tb, tbc.crys, tbc.nelec, tbc.dftenergy, scf=tbc.scf, eden=tbc.eden, gamma=tbc.gamma, background_charge_correction=tbc.background_charge_correction, within_fit=tbc.within_fit, screening=1.0, fermi_energy=tbc.efermi)
        
    elseif typeof(tbc) <: tb_crys_sparse
        println("convert to dense")        
        H = zeros(Complex{Float64}, 1, tbc.tb.nwan, tbc.tb.nwan, tbc.tb.nr)
        S = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan, tbc.tb.nr)
        for r in 1:tbc.tb.nr
            H[1,:,:,r] = tbc.tb.H[r]
            S[:,:,r] = tbc.tb.S[r]
        end
        tb = make_tb( H , tbc.tb.ind_arr, tbc.tb.r_dict,  S, h1 = tbc.tb.h1, h1spin = tbc.tb.h1spin)
        tb.scfspin = tbc.tb.scfspin
        tbc_new = make_tb_crys(tb, tbc.crys, tbc.nelec, tbc.dftenergy, scf=tbc.scf, eden=tbc.eden, gamma=tbc.gamma, background_charge_correction=tbc.background_charge_correction, within_fit=tbc.within_fit, screening=1.0, fermi_energy=tbc.efermi)
        
        
    end

    return tbc_new
    
end

###include("Crystal1563.jl")
###include("DFToutMod.jl")
#using XMLDict



####################### Wannier90 specific 
module BandStruct
"""
BS plotting
"""

using LinearAlgebra
using Plots
using Base.Threads
using ..CrystalMod:get_grid
using ..TB:make_kgrid
using ..TB:calc_energy_fft
using ..TB:tb_crys
using ..TB:tb_crys_kspace
using ..CrystalMod:crystal
using ..CrystalMod:orbital_index
using ..TB:summarize_orb

using ..Atomdata:atoms
using ..BandTools:calc_fermi
using ..BandTools:band_energy


using ..ThreeBodyTB:convert_energy
using ..ThreeBodyTB:convert_dos
using ..ThreeBodyTB:global_energy_units

using ..DFToutMod:dftout

using ..TB:find_vbm_cbm
using ..TB:tb
using ..TB:calc_bands
using ..TB:Hk

export plot_compare_tb
export plot_bandstr
export plot_compare_dft
export band_summary


global no_display=false

function set_no_display(d)
    global no_display = d
end


"""
    function plot_compare_tb(h1::tb_crys, h2::tb_crys; h3=missing)

Plot a comparison between different tight binding objects `h1`, `h2`, and optionally `h3`. Options similar to `plot_bandstr` but more limited.

"""
function plot_compare_tb(h1::tb_crys, h2::tb_crys; h3=missing, kpath=[0.5 0 0 ; 0 0 0; 0.5 0.0 0.5], names = missing, npts=30, efermi = missing, yrange=missing, plot_hk=false,  align="vbm")
    if ismissing(h3)
        plot_compare_tb(h1.tb, h2.tb, h3=missing, kpath=kpath, names = names, npts=npts, efermi = efermi, yrange=yrange, plot_hk=plot_hk, align=align)
    else
        plot_compare_tb(h1.tb, h2.tb, h3=h3.tb, kpath=kpath, names = names, npts=npts, efermi = efermi, yrange=yrange, plot_hk=plot_hk, align=align)
    end
end


"""
    function plot_compare_tb(h1::tb, h2::tb; h3=missing)
"""
function plot_compare_tb(h1::tb, h2::tb; h3=missing, kpath=[0.5 0 0 ; 0 0 0; 0.5 0.0 0.5], names = missing, npts=30, efermi = missing, yrange=missing, plot_hk=false, align="vbm")
    println("plot_compare_tb ")
    plot_bandstr(h1, kpath=kpath, names = names, npts=npts, efermi = efermi, color="green", MarkerSize=4, yrange=yrange, plot_hk=plot_hk, align=align, clear_previous=true)
    plot_bandstr(h2, kpath=kpath, names = names, npts=npts, efermi = efermi, color="yellow", MarkerSize=2, yrange=yrange, plot_hk=plot_hk, align=align, clear_previous=false)    
    if !ismissing(h3)
        plot_bandstr(h3, kpath=kpath, names = names, npts=npts, efermi = efermi, color="magenta", MarkerSize=1, yrange=yrange, plot_hk=plot_hk, align=align, clear_previous=false)
    end

end

"""
    function plot_bandstr(h::tb_crys; kpath, names = missing, proj_types=missing, proj_orbs = missing, proj_nums=missing)

Plot the band structure of a `tb_crys` object. Can also perform a projected band structure if you specify at least one of `proj_types`, `proj_orbs`, `proj_nums`.

k-path specified by a kpath array and names.

Must do scf calculation before plotting.

# Arguments
- `h::tb_crys` - The tight-biding object we want to plot bands from. Only required argument.
- `kpath=[0.5 0 0 ; 0 0 0; 0.5 0.5 0.5; 0 0.5 0.5; 0 0 0 ;0 0 0.5]` - `nk` × 3 array k-point path (high symmetry points).
- `npts=30,` - number of points between high-symmetry k-points.
- `names=missing` - `nk` string array. Names of the high-symmetry k-points 
- `proj_types=missing` - types to project onto. Either `proj_types="H"` or `proj_types=["H", "O"]` are valid.
- `proj_orbs=missing` - orbitals to project onto. either `proj_orbs=:s` or `proj_orbs=[:s, :p]`.
- `proj_nums=missing` - atom numbers to project onto. Either `proj_nums=1` or `proj_nums=[1, 2]`
- `efermi=missing` - allows you to specify fermi energy. Default is to take from `h`
- `color="blue"` - specify line color
- `MarkerSize=missing"` - specify markersize
- `yrange=missing"` - specify y-range. e.g. `yrange=[-0.7, 0.3]`
- `plot_hk=false` - plot things besides the normal band structure. Can be one of `:Seig, :Heig, :Hreal, :Himag, :Sreal, :Simag` to plot H or S eigvals or components. Primarily for debugging.
- `align="vbm"` - default or `"valence"` is to align valence band max to zero energy. Can also be `"min"`, which aligns on the minimum eigenvalue, or `"fermi"` or `"ef"`, which align on the Fermi level, 
- `clear_pervious=true` - clears the plot before adding new stuff.
- `do_display=true` - display the plot. If `false`, can be used with display-less nodes. You can still use `savefig` from `Plots` to produce saved images.
"""
function plot_bandstr(h::tb_crys; kpath=[0.5 0 0 ; 0 0 0; 0.5 0.5 0.5; 0 0.5 0.5; 0 0 0 ;0 0 0.5], names = missing, npts=30, efermi = missing, color="blue", MarkerSize=missing, yrange=missing, plot_hk=false, align = "vbm", proj_types = missing, proj_orbs = missing, proj_nums=missing, clear_previous=true, do_display=true)


    proj_inds = setup_proj(h.crys, h.tb.nwan, proj_types, proj_orbs, proj_nums)

    if h.scf && (h.energy - -999.0) < 1e-5
        println("WARNING - you have to do scf_energy before plotting to get accurate results")
    end
    if ismissing(efermi)
        efermi = h.efermi
    end
    plot_bandstr(h.tb; kpath=kpath, names = names, npts=npts, efermi = efermi, color=color, MarkerSize=MarkerSize, yrange=yrange, plot_hk=plot_hk, align=align, proj_inds=proj_inds, clear_previous=clear_previous, do_display=do_display)
    
end


function plot_bandstr(h::tb_crys_kspace; efermi = missing, color="blue", MarkerSize=missing, yrange=missing, plot_hk=false, align = "vbm", proj_types = missing, proj_orbs = missing, proj_nums=missing, clear_previous=true, do_display=true)

    kpath=h.tb.K

    proj_inds = setup_proj(h.crys, h.tb.nwan, proj_types, proj_orbs, proj_nums)

    if ismissing(efermi)
        VALS = calc_bands(h, kpath)
        energy, efermi, occs = band_energy(VALS, h.tb.kweights, h.nelec, 0.01, returnboth=true)
    end

    plot_bandstr(h.tb; kpath=kpath, npts = 1, efermi = efermi, color=color, MarkerSize=MarkerSize, yrange=yrange, plot_hk=plot_hk, align=align, proj_inds=proj_inds, clear_previous=clear_previous, do_display=do_display)
    
end



"""
    function setup_proj(crys, nwan, proj_types, proj_orbs, proj_nums)

Figure out projection indexes
"""
function setup_proj(crys, nwan, proj_types, proj_orbs, proj_nums)

    if !ismissing(proj_types) || !ismissing(proj_orbs) || !ismissing(proj_nums)
        
        if ismissing(proj_types)
            proj_types = crys.types
        end
        if ismissing(proj_orbs)
            proj_orbs = [:s, :p, :d, :f]
        end
        if ismissing(proj_nums)
            proj_nums = collect(1:crys.nat)
        end
        if typeof(proj_types) == String || typeof(proj_types) == Symbol
            proj_types = [proj_types]
        end
        if typeof(proj_orbs) == Symbol
            proj_orbs = [proj_orbs]
        end
        if typeof(proj_nums) == Int64
            proj_nums = [proj_nums]
        end

        proj_types = Symbol.(proj_types)
        
        ind2orb, orb2ind, etotal, nval = orbital_index(crys)
        
        proj_inds = Int64[]
        for n = 1:nwan
            at,t, orb = ind2orb[n]
            sorb = summarize_orb(orb)
            if (t in proj_types) && (orb in proj_orbs || sorb in proj_orbs) && (at in proj_nums)
                push!(proj_inds, n)
            end
        end
        println("proj_inds ", proj_inds)
    else
        proj_inds = missing
    end
    return proj_inds

end


"""
    function get_kpath(kpath=[0.5 0 0 ; 0 0 0; 0.5 0.5 0.5], names = missing, npts=30)

Construct a k_path for a band structure calculations. Very simple.

- `kpath` high symmetry k-points in fractional BZ coordinates.
- `names` names of kpoints like `["Γ", "X"]`
- `npts` number of points between high-symmetry kpoints
"""
function get_kpath(kpath=[0.5 0 0 ; 0 0 0; 0.5 0.5 0.5], names = missing, npts=30)

    NK = size(kpath)[1]
    K = zeros(npts *(NK-1)+1 , 3)

    nk = 0

    locs = []
    
    for kp = 1:(NK-1)

        push!(locs,nk)
        
        for x = 0:1/npts:(1-1e-5)
            nk += 1
            kt = kpath[kp,:] * (1-x) .+ kpath[kp+1,:] * x
#            println(kt)
            K[nk,:] = kt
        end

    end

#final point
    kp = NK-1
    x = 1.0
    nk += 1
    kt = kpath[kp,:] * (1-x) .+ kpath[kp+1,:] * x
#    println(kt)
    K[nk,:] = kt
    
    push!(locs,nk-1)

#    println("locs: ", locs)
    
    return K, locs
    
end



"""
    function plot_bandstr(h::tb)

Plots using `tb`
"""
function plot_bandstr(h; kpath=[0.5 0 0 ; 0 0 0; 0.5 0.5 0.5; 0 0.5 0.5; 0 0 0 ;0 0 0.5], names = missing, npts=30, efermi = missing, color="blue", MarkerSize=missing, yrange=missing, plot_hk=false, align="vbm", proj_inds=missing, clear_previous=true, do_display=true)
#function plot_bandstr( kpath; names = missing, npts=30, efermi = missing)

                
    println("plot_bandstr ")
    if clear_previous
        #println("clear")
        #        display(plot(legend=false, grid=false, framestyle=:box))
        plot(legend=false, grid=false, framestyle=:box, xtickfontsize=12, ytickfontsize=12, legendfontsize=16)
    end

    if ismissing(MarkerSize)
        if ismissing(proj_inds)
          MarkerSize = 0
        else
            MarkerSize = 5
        end
    end


                

    NK = size(kpath)[1]

    K, locs = get_kpath(kpath, names , npts)
    nk = size(K)[1]
    
    proj = zeros(nk, h.nwan)

    if plot_hk != false

        if plot_hk == :Seig ||  plot_hk == :Heig
            vals = zeros(nk, h.nwan)
        else
            vals = zeros(nk, h.nwan*h.nwan)
        end
            
        for i = 1:nk
            vect, vals_t, hk, sk, vals0 = Hk(h, K[i,:])
            if plot_hk == :Hreal
                vals[i,:] = real(hk[:])
            elseif plot_hk == :Himag
                vals[i,:] = imag(hk[:])
            elseif plot_hk == :Sreal
                vals[i,:] = real(sk[:])
            elseif plot_hk == :Simag
                vals[i,:] = imag(sk[:])
            elseif plot_hk == :Seig
                F=eigen(sk)
                vals[i,:] = real(F.values)
            elseif plot_hk == :Heig
                F=eigen(hk)
                vals[i,:] = real(F.values)
            else
                vals = calc_bands(h, K)
            end
        end

    elseif !ismissing(proj_inds )
        vals = zeros(nk, h.nwan)
        temp = zeros(Complex{Float64}, h.nwan, h.nwan)
#        println("proj inds $proj_inds")
        for i = 1:nk
            vect, vals_t, hk, sk, vals0 = Hk(h, K[i,:])

            for p in proj_inds
                for a = 1:h.nwan
                    for j = 1:h.nwan
                        t = vect[p,a]*conj(vect[j,a])
                        proj[i, a] += 0.5*real( (t*sk[j,p]  + conj(t)* conj(sk[j,p])))
                    end
                end
            end
            
#            sk5 = ( 0.5*(sk+sk'))^0.5
#            sv = sk5 * vect
            
            #if i == 1
            #    println("K ", K[i,:])
            #    println("vect")
            #    println(vect)
            #    println("vals_t")
            #    println(vals_t)
            #    println("sk")
            #    println(sk)
            #end

#            for p in proj_inds
#                proj[i,:] +=  real(conj(sv[p,:]) .* sv[p,:])
#            end
                
                #                for a = 1:h.nwan
 #                   for b = 1:h.nwan
#                        temp[a,b] = vect[a,p]'*sk[a, b]* vect[b,p]

                        #temp[a,b] = vect[p,a]'*sk[a, b]* vect[p, b]
#                    end
#                end
#                temp = temp + conj(temp)
                
#                proj[i,:] += 0.5*sum(real(temp), dims=1)[:]

            #               proj[i,:] += real(vect[:,

#sum(vect[:,2:4] .* ( sk * vect[:,2:4]), dims=2)
#            proj[i,:] = real(sum( conj(vect[:,proj_inds]) .* (sk * vect[:,proj_inds]) , dims=2))
#            end
#            println("$i proj ", proj[i,:])
#            println("check ", sum(vect'*sk*vect))
            
            #            println(vect[proj_inds,:] .* conj(vect[proj_inds,:]))
#            println()
#            proj[i,:] = sum(vect[proj_inds,:]*vect[proj_inds, :]', dims=1)
            #            proj[i,:] = real(diag(vect[proj_inds,:]'*sk[proj_inds, proj_inds]*vect[proj_inds,:]))
            vals[i,:] = vals_t
        end

        
    else
        vals = calc_bands(h, K)
    end


    if global_energy_units == "eV"
        units = "eV"
    else
        units = "Ryd"
    end
        
    alignstr = "Energy ($units)"


    if !ismissing(align)
        println("align ", align)
        align = lowercase(align)
        if align == "min" || align == "minimum"
            vmin = minimum(vals)
            vals = vals .- vmin
            alignstr = "Energy - Emin ($units)"
        elseif align == "fermi" || align == "ef" 
            vals = vals .- efermi
            alignstr = "Energy - E_F ($units)"
        elseif align == "vbm" || align == "valence"
            vbm = maximum(vals[vals .< efermi])
            vals = vals .- vbm
            alignstr = "Energy - VBM ($units)"
        end            
    end


    #    println("yyyyyyyyyyyyyyyyy")
    if ismissing(proj_inds)
        println("color = $color markersize = $MarkerSize")
        plot!(convert_energy(vals), color=color, lw=2.0, marker=(:circle), markersize=MarkerSize, markerstrokecolor=color, legend=false, grid=false)
#        plot!(vals, color="red", marker=(:circle), markersize=MarkerSize, markeredgecolor=color, legend=false, grid=false)
#        plot(vals, color=color, marker=(:circle), markersize=0.1)
    else
        X=repeat(1:(nk), 1,h.nwan)
        
        plot!(X, convert_energy(vals), color="grey", grid=false, legend=false)
        scatter!(X, convert_energy(vals), marker_z = proj,markersize=MarkerSize, markerstrokewidth=0.0, alpha=0.6, legend=false)

#        @save "t.jld" vals proj
    end


    if (ismissing(names))
        names = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"]
        names = [names;names;names;names;names;names;names; names; names; names; names;names;names;names ]
        names = names[1:NK]
    end

    if !(ismissing(yrange))
        if size(yrange) == (1,2)
            yrange=yrange'
        end
        ylims!(convert_energy(yrange[1]), convert_energy(  yrange[2]))

    end

    if ismissing(yrange)
        maxV = min(maximum(vals), 1.0)
        r = maxV - minimum(vals) 
        yrange = [minimum(vals)-0.05*r,  maxV+0.05*r]
    end
    
    ylims!(convert_energy(yrange[1]), convert_energy(yrange[2]))
    xlims!(1, nk)

    if length(locs) < 15
        for l in locs
            lt=l+1
            if lt == 1
                lt == lt + 1e-5
            elseif lt == nk
                lt = lt - 1e-5
            end
            plot!([lt, lt], convert_energy(yrange), linestyle=:dash, color="black")
        end
    end

    if !ismissing(align) 
#        println("nk ", nk)
        if align == "vbm" || align == "valence"
            plot!([0, nk+1], convert_energy([efermi - vbm, efermi - vbm]), color="gray", linestyle=:dot)
        end
        if align == "min" || align == "minimum"
            plot!([0, nk+1], convert_energy([efermi - vmin, efermi - vmin]), color="gray", linestyle=:dot)
        end
        if align == "fermi" || align == "ef"
            plot!([0, nk+1], [0.0, 0.0], color="black", linestyle=:dot)
        end

    end

#    println(locs)
#    println(names)    
    xticks!(Float64.(locs).+1.0, names, xtickfontsize=12)

    if do_display && no_display == false
        println("display")
        display(ylabel!(alignstr, guidefontsize=16))
    else
        println("no display")
        ylabel!(alignstr, guidefontsize=16)
    end
    
    #    println("end plot")
    
end


"""
    function plot_compare_dft(tbc, bs; tbc2=missing)
    
Plots a band structure comparison between a tight-binding crystal object (`tb_crys`) and a
band structure directly from dft (either a `dftout` or `bs` object). 

The k-points are fixed by the `bs` object.

`tbc2` is an optional second `tbc_crys`.
"""
function plot_compare_dft(tbc, bs; tbc2=missing)

    if typeof(bs) == dftout
        bs = bs.bandstruct
    end

    kpts = bs.kpts
    kweights = bs.kweights

    
    vals = calc_bands(tbc.tb, kpts)
    

    nelec = tbc.nelec

    nsemi = 0
    for t in tbc.crys.types
        atom = atoms[t]
        nsemi += atom.nsemicore
    end
    nsemi = Int64(nsemi / 2)

    efermi_dft = calc_fermi(bs.eigs[:,1:end], kweights, nelec+nsemi*2)
    efermi_tbc = calc_fermi(vals, kweights, nelec)

    if !ismissing(tbc2)
        vals2 = calc_bands(tbc2.tb, kpts)
        efermi_tbc2 = calc_fermi(vals2, kweights, nelec)
#        println("efermi_dft $efermi_dft efermi_tbc $efermi_tbc efermi_tbc2 $efermi_tbc2")

        en_dft = band_energy(bs.eigs[:,1:end], kweights, nelec+nsemi)
        en_1 = band_energy(vals, kweights, nelec)
        en_2 = band_energy(vals2, kweights, nelec)
#        println("energy_dft $en_dft en_1 $en_1 en_2 $en_2")

        e_smear_dft  = smearing_energy(bs.eigs[:,nsemi+1:end], kweights, efermi_dft)
        e_smear_1  = smearing_energy(vals, kweights, efermi_tbc)
        e_smear_2  = smearing_energy(vals2, kweights, efermi_tbc2)
        
#        println("smear_dft $e_smear_dft sm_1 $e_smear_1 sm_2 $e_smear_2")

    end


    nk = size(kpts)[1]

#    wan, semicore, nwan, nsemi, wan_atom, atom_wan = tb_indexes(tbc.crys)
    println("nsemi ", nsemi)

    
    
#    if align_min
#        a = minimum(bs.eigs[:,nsemi+1:end])
#        b = minimum(vals)
#        plot(bs.eigs .- a, color="orange",  LineWidth=1.5)    
#        plot!(vals .- b, color="blue", LineWidth=1.0)
#
#        plot!([0, nk], [efermi_tbc - b, efermi_tbc - b], color="red", linesyle=:dash)
#        if !ismissing(tbc2)
#            c = minimum(vals2)
#            plot!(vals2 .- c, color="magenta", linestyle=:dash, markersize=3)
#        end
#
#        display(plot!([0, nk], [efermi_dft - a, efermi_dft - a], color="cyan", linesyle=:dot))



#    else
    
    vbmD, cbmD = find_vbm_cbm(bs.eigs[:,1:end] , efermi_dft)
    vbm, cbm = find_vbm_cbm(vals , efermi_tbc)

#    vbm = 0.0
#    vbmD = 0.0

    plot(legend=true, grid=false, framestyle=:box, xtickfontsize=16, ytickfontsize=16)

    
    plot!(convert_energy( bs.eigs[:,1] .- vbmD) , color="orange", lw=3, label="DFT", grid=false, legend=:topright)    
    if bs.nbnd > 1
        plot!(convert_energy(bs.eigs[:,2:end] .- vbmD) , color="orange", lw=3, label=missing)    
    end

    plot!( convert_energy(vals[:,1] .- vbm), color="blue",  lw=2, label="TB")
    if tbc.tb.nwan > 1
        plot!(convert_energy(vals[:,2:end] .- vbm), color="blue",  lw=2, label=missing)
    end



    if !ismissing(tbc2)
        vbm2, cbm2 = find_vbm_cbm(vals2 , efermi_tbc2)
        plot!(convert_energy(vals2[:,1] .- vbm2),  color="cyan", lw=0.5, label="TB2")
        if tbc.tb.nwan > 1
            plot!(convert_energy(vals2[:,2:end] .- vbm2), lw=0.5, color="cyan", label=missing)
        end
    end


    plot!([0, nk], convert_energy([efermi_dft, efermi_dft] .- vbmD), color="red", linestyle=:dot, label="E_F DFT")
    plot!([0, nk], convert_energy([efermi_tbc, efermi_tbc] .- vbm), color="black", linestyle=:dash, label="E_F TB")

    if global_energy_units == "eV"
        ylabel!("Energy - VBM (eV)", guidefontsize=16)
    else
        ylabel!("Energy - VBM (Ryd.)", guidefontsize=16)
    end

    if ! no_display
        display(ylims!(convert_energy(minimum(vals .- vbm)*1.05), convert_energy(min(maximum(vals .- vbm), 0.8) * 1.05 )))
    else
        ylims!(convert_energy(minimum(vals .- vbm)*1.05), convert_energy(min(maximum(vals .- vbm), 0.8) * 1.05 ))
    end
#    end


end


"""
    function band_summary(tbc, kgrid, fermi=missing)
    
Produces summary of band structure. See below functions for more
specific versions of function that automatically generate the k-points.

Note: gaps are not well-defined for non-magnetic systems with odd
numbers of electrons, as they are required to be metals.

Returns `direct_gap, indirect_gap, gaptype, bandwidth`

-`direct_gap`: minimum gap at one k-point between nominally filled and empty bands. Can be non-zero in metals.  
-`indirect_gap`: LUMO - HOMO. Can be negative if material has a direct gap everywhere, but the conduction band at some k-point is below the valence band at a different k-point. Physically these are indirect gap semimetals.  
-`gaptype` : is `:metal` for all metals, `:direct` or `:indirect` for insulators.
-`bandwidth` : HOMO - minimum_band_energy. Included semicore states if they are in the TB calculation.

"""
function band_summary(tbc, kpts, kweights, fermi=missing)


    vals = calc_bands(tbc, kpts)    

    if ismissing(fermi)
        fermi = calc_fermi(vals, kweights, tbc.nelec)
    end


    
    minband = minimum(vals[:,1])
    
    nelec = Int64(round(tbc.nelec/2.0))

    if abs(nelec - tbc.nelec / 2.0) > 1e-5
        println("WARNING, band_summary gaps don't work well with non-integer or odd electrons: ", tbc.nelec)
        directgap = 0.0
        indirectgap = 0.0
        vbm = fermi

        return directgap, indirectgap, :metal, convert_energy(vbm - minband)
        
    end

    #this if the even number of electrons case
    
    vbm = maximum(vals[:,nelec])
    vbm = min(vbm, fermi)
    
    if nelec+1 <= tbc.tb.nwan
        cbm = minimum(vals[:,nelec+1])
        directgap = minimum(vals[:,nelec+1] - vals[:,nelec])

    else
        cbm = vbm
        directgap = 0.0
    end
    
    indirectgap = cbm - vbm

    if directgap - indirectgap > 1e-5
        gaptype = :indirect
    else
        gaptype = :direct
    end

    if indirectgap < 0.0
        gaptype = :metal
    end

    bandwidth = vbm - minband
    
    return convert_energy(directgap), convert_energy(indirectgap), gaptype, convert_energy(bandwidth)

end

"""
    function band_summary(tbc::tb_crys; kgrid=missing, kpts=missing)

Will automatically generate standard k-grid by default.
"""
function band_summary(tbc::tb_crys; kpts=missing, kgrid=missing)
    if ismissing(kgrid)
        kgrid = get_grid(tbc.crys)
    end
    if ismissing(kpts)
        println("using kgrid $kgrid")
        kpts, kweights = make_kgrid(kgrid)
    end
    
    return band_summary(tbc, kpts, tbc.efermi)
end

"""
    function band_summary(tbc::tb_crys_kspace)

Will use internal k-points by default.
"""
function band_summary(tbc::tb_crys_kspace)
    kpts = tbc.tb.K
    return band_summary(tbc, kpts, missing)
end




end #end module

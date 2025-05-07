###include("Crystal1563.jl")
###include("DFToutMod.jl")
#using XMLDict



#######################
module DOS
"""
DOS
"""

using LinearAlgebra
using Plots
using Base.Threads
using ..CrystalMod:get_grid
using ..TB:calc_energy_fft
using ..TB:tb_crys
using ..TB:tb_crys_sparse
using ..TB:tb_crys_kspace
using ..CrystalMod:crystal
using ..CrystalMod:orbital_index
using ..TB:summarize_orb
using ..TB:get_energy_electron_density_kspace
using ..ThreeBodyTB:convert_energy
using ..ThreeBodyTB:convert_length
using ..ThreeBodyTB:convert_dos
using ..ThreeBodyTB:global_energy_units
using ..ThreeBodyTB:global_length_units
using ..ThreeBodyTB:no_display
using ..Symmetry:get_kgrid_sym
using ..DFToutMod:dftout
using ..DFToutMod:bandstructure
using ..TB:find_vbm_cbm
using ..AtomicProj:proj_dat
using ..TB:tb_indexes

function get_projtype(tbc, ptype=missing)

    ind2orb, orb2ind, etotal, nval = orbital_index(tbc.crys)
    
    if ismissing(ptype)
        if length(Set(tbc.crys.types)) != 1
            ptype=:atomic
        else
            ptype=:orbs
        end
        println("Projection type: ", ptype)
    end

    names = []
    PROJ = []       
    pwan=[]
    
    if ptype == :all || ptype == :All || ptype == "All" || ptype == "all"
        orblist = []
        for n = 1:tbc.tb.nwan
            at,t, orb = ind2orb[n]
            push!(names, "$at $t $orb")
            push!(PROJ, [n])
            push!(pwan, 1)
        end
    elseif ptype == :atomic || ptype == "atomic" || ptype == :atoms || ptype == "atoms" || ptype == :atom || ptype == "atom" || ptype == :Atomic || ptype == "Atomic"
        for ti in Set(tbc.crys.stypes)
            
            proj_inds = Int64[]
            for n = 1:tbc.tb.nwan
                at,t, orb = ind2orb[n]
                sorb = summarize_orb(orb)
                if t == ti
                    push!(proj_inds, n)
                end
            end

            push!(names, String(ti))
            push!(PROJ, proj_inds)
            push!(pwan, length(proj_inds))
        end

        println("PROJ")
        for a in zip(names, PROJ, pwan)
            println(a)
        end
        
    else
        
        for o in [:s, :p, :d, :f]
            proj_inds = Int64[]
            for n = 1:tbc.tb.nwan
                at,t, orb = ind2orb[n]
                sorb = summarize_orb(orb)
                if o == sorb
                    push!(proj_inds, n)
                end
            end
            if length(proj_inds) > 0
                push!(names, String(o))
                push!(PROJ, proj_inds)
                push!(pwan, length(proj_inds))                
            end

        end
    end

    return names, PROJ, pwan

end


"""
    function projection(tbc::tb_crys, vects, sk3, grid; ptype=missing)

Figures out the projections.
`ptype` can be `:atomic` or `:orbs` for atom projection or orbital projection (:s,:p,:d)
Default is to choose `:atomic` except for elemental systems.
"""
function projection(tbc::tb_crys, vects, sk3, grid; ptype=missing, use_sym=false, nk_red=1, grid_ind=missing, kweights = missing)    

    nspin = tbc.nspin
    if tbc.tb.scfspin == true
        nspin = 2
    end

    names, PROJ, pwan = get_projtype(tbc, ptype)

    nk = size(vects)[1]

    temp = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)

#    sk5 = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)
    sk = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)    

    vectsS = zeros(size(vects))

    if use_sym
        proj = zeros(nk_red, tbc.tb.nwan, length(PROJ), nspin)
        
        for spin = 1:nspin
            for c = 1:nk_red
                k1,k2,k3 = grid_ind[c,:]
#                k3 = mod(c-1 , grid[3])+1
#                k2 = 1 + mod((c-1) ÷ grid[3], grid[2])
#                k1 = 1 + (c-1) ÷ (grid[2]*grid[3])
                
                #            skt = 0.5*(sk3[:,:,k1,k2,k3] + sk3[:,:,k1,k2,k3]')
                #            skt5 = skt^0.5
                if length(size(sk3)) == 5
                    v5 = (0.5*( ( sk3[:,:,k1,k2,k3]) .+ ( sk3[:,:,k1,k2,k3])'))^0.5 * ( vects[c,spin,:,:])
                else

                    v5 = (0.5*( collect( sk3[c] .+  sk3[c]')))^0.5 * ( vects[c,spin,:,:])
                end
                          

                
                vectsS[c,spin,:,:] = real(v5.*conj(v5))
            end
        end

    else
        proj = zeros(nk, tbc.tb.nwan, length(PROJ), nspin)
        
        for spin = 1:nspin
            @threads for c = 1:grid[1]*grid[2]*grid[3]
                
                k3 = mod(c-1 , grid[3])+1
                k2 = 1 + mod((c-1) ÷ grid[3], grid[2])
                k1 = 1 + (c-1) ÷ (grid[2]*grid[3])
                
                #            skt = 0.5*(sk3[:,:,k1,k2,k3] + sk3[:,:,k1,k2,k3]')
                #            skt5 = skt^0.5
                if length(size(sk3)) == 5 
                    v5 = (0.5*( ( sk3[:,:,k1,k2,k3]) .+ ( sk3[:,:,k1,k2,k3])'))^0.5 * ( vects[c,spin,:,:])
                else
                    v5 = (0.5*( ( sk3[:,:,k1,k2,k3]) .+ ( sk3[:,:,k1,k2,k3])'))^0.5 * ( vects[c,spin,:,:])                    
                end


                
                vectsS[c,spin,:,:] = real(v5.*conj(v5))
            end
        end
    end
#    println("sum vS  ", sum(vectsS))
    
   #= for spin = 1:nspin
        for (pind, proj_inds) in enumerate(PROJ)
            for p in proj_inds
                for a = 1:tbc.tb.nwan
                    for j = 1:tbc.tb.nwan
                        t = conj(vects[:,spin,p,a]).*vects[:,spin,j,a]
                        @threads for c = 1:grid[1]*grid[2]*grid[3]
                            
                            k3 = mod(c-1 , grid[3])+1
                            k2 = 1 + mod((c-1) ÷ grid[3], grid[2])
                            k1 = 1 + (c-1) ÷ (grid[2]*grid[3])
                            
                            proj[c,a,pind, spin] += real( t[c] * sk3[p,j,k1,k2,k3])
                            #proj[c,a, pind, spin] += 0.5*real( (t[c]*sk3[j,p,k1,k2,k3]  + conj(t[c])* conj(sk3[j,p, k1,k2,k3])))
                            
                            
                        end
                    end
                end
            end
        end
    end
=#

    if use_sym
        for spin = 1:nspin
            for (pind, proj_inds) in enumerate(PROJ)
                println("proj_inds $proj_inds")
                for p1 in proj_inds
                    #                for p2 in proj_inds
                    #                    for a = 1:tbc.tb.nwan
                    @threads for c = 1:nk_red
                        
                        for a = 1:tbc.tb.nwan
                            
                            proj[c,a,pind, spin] += kweights[c]*vectsS[c,spin,p1,a] * nk/2.0

                        end
                        
                    end
                    #                    end
                end
            end
        end

    else

        for spin = 1:nspin
            for (pind, proj_inds) in enumerate(PROJ)
                println("proj_inds $proj_inds")
                for p1 in proj_inds
                    #                for p2 in proj_inds
                    #                    for a = 1:tbc.tb.nwan
                    @threads for c = 1:grid[1]*grid[2]*grid[3]
                        
                        for a = 1:tbc.tb.nwan
                            
                            proj[c,a,pind, spin] += vectsS[c,spin,p1,a]

                        end
                        
                    end
                    #                    end
                end
            end
        end
    end
    return proj, names, pwan
    
    
end
            

"""
    function dos_realspace(tbc; direction=3, grid=missing, smearing=0.005, npts=missing, do_display=true, use_sym=false, energy_grid=100, width=missing, energy_lims=missing, scissors_shift=0.0, scissors_shift_atoms=[])


Makes an atom and spatially resoloved DOS-style plot. Along `direction=1,2,3` (lattice vector 1,2,3) the atoms are plotted, using transparancy to denote DOS amplitude. Try it out to see. Potentially useful for interfaces. Example usage:
`c = makecrys([8.0 0 0; 0 8.0 0; 0 0 10.0], [0 0 0; 0 0 0.5], [:Li, :Cl]) #quasi 1D LiCl chain`
`en, tbc, flag = scf_energy(c*[1,1,10])                                   #ten unit cells`
`dos_realspace(tbc, direction=3)`

Adjust the `energy_grid` and `width` to adjust plotting parameters.
Scissors shift in energy unit to certain atoms will add scissors shift to open band gaps is desired.
"""
function dos_realspace(tbc; direction=3, grid=missing, smearing=0.005, npts=missing, do_display=true, use_sym=false, energy_grid=100, width=missing, energy_lims=missing, scissors_shift=0.0, scissors_shift_atoms=[])

    if ismissing(width)
        width= convert_length(3.0)
    end
    if ismissing(energy_lims)
        energy_lims = convert_energy([-0.3, 0.3])
    end
    
    colors = ["blue", "orange", "green", "magenta", "cyan", "red", "yellow"]
    colors = [colors; colors; colors]
    atom_colors = Dict()
    for c in colors
        for s in Set(tbc.crys.stypes)
            if !(s in keys(atom_colors))
                atom_colors[s] = c
                println("color $s $c ")
                break
            end
        end
    end

    ind2orb, orb2ind, etotal, nval = orbital_index(tbc.crys)
    
    energies,DOS, pdos, names, proj = dos(tbc, grid=missing, smearing=smearing, npts=missing, do_display=false, use_sym=use_sym, proj_type=:all, scissors_shift=scissors_shift, scissors_shift_atoms=scissors_shift_atoms)

    pdos=sum(pdos, dims=[3])
    c = tbc.crys
    coords_cart = convert_length(c.coords[:,direction] * sqrt(sum(c.A[direction,:].^2)))


    

    dos_atom = zeros(length(energies), c.nat)
    for n = 1:tbc.tb.nwan
        at,t, orb = ind2orb[n]
        dos_atom[:,at] += pdos[:,n]
    end
    dos_atoms_grid = zeros(energy_grid, c.nat)
    #new_grid = collect(minimum(energies) : (maximum(energies)-minimum(energies)) / energy_grid: (maximum(energies)+1e-10))

    println("grid width ",  (energy_lims[2]-energy_lims[1]) / energy_grid)
    
    new_grid = collect( energy_lims[1] : (energy_lims[2]-energy_lims[1]) / energy_grid: (energy_lims[2]+1e-10) )

    println("size(new_grid), $(size(new_grid))")
    for en = 1:energy_grid
        ind = energies .> new_grid[en] .&& energies .<= new_grid[en+1]
        dos_atoms_grid[en,:] = sum(dos_atom[ind,:], dims=[1])
    end
    max_den = maximum(dos_atom) * 1.0

    rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

    height = (energy_lims[2]-energy_lims[1]) / energy_grid
    xmin = minimum(coords_cart) - width/2
    xmax = maximum(coords_cart) + width/2

    plot([xmin, xmax], [0.0, 0.0], linestyle=:dash, linewidth=1, color="grey", label="")
    xlims!(xmin, xmax + width*2)
    ylims!(energy_lims[1], energy_lims[2])

    println([ minimum(coords_cart) - width,maximum(coords_cart) + width])
    println([minimum(new_grid) - 0.01, maximum(new_grid) + 0.01])
    
    for at = 1:c.nat
        for ind = 1:energy_grid
            x = coords_cart[at] - width/2
            plot!(rectangle(width, height, x, new_grid[ind] - height/2), opacity= dos_atoms_grid[ind, at] / max_den , color=atom_colors[c.stypes[at]], label="" )
            #            println("plot $width, $height, $x, $(new_grid[ind] - height/2)")
        end
    end

    for at = keys(atom_colors)
        plot!(rectangle(0.0,0.0 , coords_cart[1] - width/2, new_grid[1] - height/2), color=atom_colors[at], label=String(at))
    end

    xlabel!("Atom positions ($global_length_units)",  guidefontsize=12)
    display(ylabel!("Energy - Fermi ($global_energy_units)", guidefontsize=12))
    return new_grid, dos_atoms_grid, coords_cart

end


"""
    function dos_realspace(tbc; direction=3, grid=missing, smearing=0.005, npts=missing, do_display=true, use_sym=false, energy_grid=100, width=missing, energy_lims=missing, scissors_shift=0.0, scissors_shift_atoms=[])


Makes an atom and spatially resoloved DOS-style plot. Along `direction=1,2,3` (lattice vector 1,2,3) the atoms are plotted, using transparancy to denote DOS amplitude. Try it out to see. Potentially useful for interfaces. Example usage:
`c = makecrys([8.0 0 0; 0 8.0 0; 0 0 10.0], [0 0 0; 0 0 0.5], [:Li, :Cl]) #quasi 1D LiCl chain`
`en, tbc, flag = scf_energy(c*[1,1,10])                                   #ten unit cells`
`dos_realspace(tbc, direction=3)`

Adjust the `energy_grid` and `width` to adjust plotting parameters.
Scissors shift in energy unit to certain atoms will add scissors shift to open band gaps is desired.
"""
function dos_realspace(p::proj_dat, dft::dftout; direction=3, grid=missing, smearing=0.005, npts=missing, do_display=true, use_sym=false, energy_grid=100, width=missing, energy_lims=missing, scissors_shift=0.0, scissors_shift_atoms=[])

    if ismissing(width)
        width= convert_length(3.0)
    end
    if ismissing(energy_lims)
        energy_lims = convert_energy([-0.3, 0.3])
    end
    
    colors = ["blue", "orange", "green", "magenta", "cyan", "red", "yellow"]
    colors = [colors; colors; colors]
    atom_colors = Dict()
    for c in colors
        for s in Set(dft.crys.stypes)
            if !(s in keys(atom_colors))
                atom_colors[s] = c
                println("color $s $c ")
                break
            end
        end
    end

    ind2orb, orb2ind, etotal, nval = orbital_index(dft.crys)
    wan, semicore, nwan, nsemi, wan_atom, atom_wan = tb_indexes(dft)
    
    vmin = minimum(p.bs.eigs) - 0.1
    vmax = maximum(p.bs.eigs) + 0.1

    println("v $vmin $vmax ")

#    vmax = min(maximum(vals), 5.0)
    r = vmax - vmin

    if ismissing(npts)
        npts = Int64(round(r * 100 ))
    end
    
    energies = collect(vmin - r*0.02 : r*1.04 / npts    : vmax + r*0.02 + 1e-7)

    energies = energies .- dft.bandstruct.efermi

#    energies,DOS, pdos, names, proj = dos(tbc, grid=missing, smearing=smearing, npts=missing, do_display=false, use_sym=use_sym, proj_type=:all, scissors_shift=scissors_shift, scissors_shift_atoms=scissors_shift_atoms)
#
 #   pdos=sum(pdos, dims=[3])
    c = dft.crys
    coords_cart = convert_length(c.coords[:,direction] * sqrt(sum(c.A[direction,:].^2)))

#    pdos = zeros( length(energies), c.nat, 
#    for k = 1:p.bs.nks
#        for spin = 1:p.bs.nspin
#            for orb = 1:p.natwfc
#                atom = wan_atom[orb]
#                
#            end
#        end
#    end
    
    pdos = zeros(length(energies), nwan, p.bs.nspin)
    
    for spin = 1:p.bs.nspin
#        for (c,e) in enumerate(energies)
#            DOS[c,spin] = sum(KW .*  exp.( -0.5 * (vals[:,:,spin] .- e).^2 / smearing^2 ) ) / 2.0 * nk
#        end
        
#        if do_proj
#            println(size(proj))
        for (i,w) in enumerate(wan)
            for (c,e) in enumerate(energies)
                for k = 1:p.bs.nks
                    for b = 1:p.bs.nbnd
                        pdos[c, i] += real(p.proj[k,w,spin,b] .*  exp.( -0.5 * (p.bs.eigs[k,b,spin] - dft.bandstruct.efermi .- e).^2 / smearing^2 ) ) * p.bs.kweights[k]
                    end
                end
            end
        end
    end
    
    pdos = pdos / smearing / (2.0*pi)^0.5 / p.bs.nks
        
    

    dos_atom = zeros(length(energies), c.nat)
    for n = 1:nwan
        at,t, orb = ind2orb[n]
        dos_atom[:,at] += pdos[:,n]
    end
    dos_atoms_grid = zeros(energy_grid, c.nat)
    #new_grid = collect(minimum(energies) : (maximum(energies)-minimum(energies)) / energy_grid: (maximum(energies)+1e-10))

    println("grid width ",  (energy_lims[2]-energy_lims[1]) / energy_grid)
    
    new_grid = collect( energy_lims[1] : (energy_lims[2]-energy_lims[1]) / energy_grid: (energy_lims[2]+1e-10) )

    println("size(new_grid), $(size(new_grid))")
    for en = 1:energy_grid
        ind = energies .> new_grid[en] .&& energies .<= new_grid[en+1]
        dos_atoms_grid[en,:] = sum(dos_atom[ind,:], dims=[1])
    end
    max_den = maximum(dos_atom) * 1.0

    rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

    height = (energy_lims[2]-energy_lims[1]) / energy_grid
    xmin = minimum(coords_cart) - width/2
    xmax = maximum(coords_cart) + width/2

    plot([xmin, xmax], [0.0, 0.0], linestyle=:dash, linewidth=1, color="grey", label="")
    xlims!(xmin, xmax + width*2)
    ylims!(energy_lims[1], energy_lims[2])

    println([ minimum(coords_cart) - width,maximum(coords_cart) + width])
    println([minimum(new_grid) - 0.01, maximum(new_grid) + 0.01])
    
    for at = 1:c.nat
        for ind = 1:energy_grid
            x = coords_cart[at] - width/2
            plot!(rectangle(width, height, x, new_grid[ind] - height/2), opacity= dos_atoms_grid[ind, at] / max_den , color=atom_colors[c.stypes[at]], label="" )
            #            println("plot $width, $height, $x, $(new_grid[ind] - height/2)")
        end
    end

    for at = keys(atom_colors)
        plot!(rectangle(0.0,0.0 , coords_cart[1] - width/2, new_grid[1] - height/2), color=atom_colors[at], label=String(at))
    end

    xlabel!("Atom positions ($global_length_units)",  guidefontsize=12)
    display(ylabel!("Energy - Fermi ($global_energy_units)", guidefontsize=12))
    return new_grid, dos_atoms_grid, coords_cart

end

function projection(tbcK::tb_crys_kspace, vects, SK; ptype=missing)    

    names, PROJ, pwan = get_projtype(tbcK, ptype)

    nk = size(vects)[1]
    nspin = size(vects)[2]
    temp = zeros(Complex{Float64}, tbcK.tb.nwan, tbcK.tb.nwan)
    proj = zeros(nk, tbcK.tb.nwan, length(PROJ), nspin)

#    sk5 = zeros(Complex{Float64}, tbcK.tb.nwan, tbcK.tb.nwan)
    sk = zeros(Complex{Float64}, tbcK.tb.nwan, tbcK.tb.nwan)    
    c=0

    for k = 1:nk
        c += 1
        if length(size(SK)) == 3
            sk[:,:] = ( 0.5 * (SK[ :, :,k] + SK[ :, :, k]'))
        else
            sk[:,:] = ( 0.5 * (SK[k] + SK[k]'))
        end            
        for spin = 1:nspin
            for (pind, proj_inds) in enumerate(PROJ)
                for p in proj_inds
                    for a = 1:tbcK.tb.nwan
                        for j = 1:tbcK.tb.nwan
                            t = vects[c,spin, p,a]*conj(vects[c,spin, j,a])
                            proj[c,a, pind, spin] += 0.5*real( (t*sk[j,p]  + conj(t)* conj(sk[j,p])))
                        end
                    end
                end
            end
        end
    end
        
    return proj, names, pwan
    
    
end
            
    
"""
    function gaussian_dos(tbc::tb_crys)

Gaussian smearing DOS is now the main DOS. Tetraheral method is still coded, but it is a bit wonky.
"""
function gaussian_dos(tbc::tb_crys; grid=missing, smearing=0.04, npts=missing, proj_type=missing, do_display=true, scissors_shift=0.0, scissors_shift_atoms=[], align=:efermi)

    return dos(tbc, grid=grid, smearing=smearing, npts=npts, proj_type=proj_type, do_display=do_display, scissors_shift=scissors_shift, scissors_shift_atoms=scissors_shift_atoms, align=align)
    
end

"""
    function dos(tbc::tb_crys; grid=missing, smearing=0.04, npts=missing, proj_type=missing, do_display=true)

DOS, mostly for testing.

- `npts` is number of energies
- `proj_type` can be `"none"`, `"atomic"`, or `"orbital"`
- `do_display=false` will suppress the actual plot

The combination of `smearing` and `grid` are important to get converged results.

See also `dos`

`return energies, dos, projected_dos, pdos_names`
"""
function dos(tbc::tb_crys; grid=missing, smearing=0.03, npts=missing, proj_type=missing, do_display=true, use_sym=false, scissors_shift=0.0, scissors_shift_atoms=[], align=:efermi)


    if typeof(tbc) <: tb_crys_sparse
        use_sym=true
    end
        
    println("proj_type $proj_type")
    if ismissing(grid)
        grid = get_grid(tbc.crys)
        grid = Int64.(round.(grid * 1.6))
        println("grid $grid")
    end

    etot, efermi, vals, vects, hk3, sk3 = calc_energy_fft(tbc, grid=grid, smearing=smearing, return_more_info=true, use_sym=use_sym, scissors_shift= scissors_shift, scissors_shift_atoms=scissors_shift_atoms)

    println("dos fermi $efermi xxxxxxxxxxxxx---------------")
    
    nspin = size(vals)[end]
    println("nspin $nspin")
    
#    println("precheck ", sum(vects[1,:,:]' * sk3[:,:,1,1,1] * vects[1,:,:]))
    
    #prelim
    if align == :efermi
        vals = vals .- efermi
    else
        vbm, cbm = find_vbm_cbm(vals, efermi)
        vals = vals .- vbm
    end
    
    nk = size(vals)[1]
    
    vmin = minimum(vals) - 0.05
    vmax = min(maximum(vals)+0.05, 5.0)
    r = vmax - vmin

    if ismissing(npts)
        npts = 500
    end
    
    energies = collect(vmin - r*0.2 : r*1.04 / npts    : vmax + r*0.02 + 1e-7)
    
    DOS = zeros(length(energies))

    if use_sym
        nk_red, grid_ind, kpts, kweights = get_kgrid_sym(tbc.crys, grid=grid)
    else
        nk_red=1
        grid_ind =missing
        kpts = missing
        kweights = ones(prod(grid)) * 2 / nk
    end
        
    
    if ismissing(proj_type)
        do_proj=true
        proj, names, pwan =  projection(tbc, vects, sk3, grid, ptype=proj_type, use_sym=use_sym, nk_red=nk_red, grid_ind=grid_ind, kweights = kweights)
        nproj = size(proj)[3]
#        println("proj $proj names $names pwan $pwan")
        pdos = zeros(length(energies),nproj, nspin)

    elseif proj_type == "none" ||  proj_type == :none
        do_proj=false
        nproj=0
        pdos=missing
        names=missing
        proj = missing
    else
        do_proj=true
        proj, names, pwan =  projection(tbc, vects, sk3, grid, ptype=proj_type, use_sym=use_sym, nk_red=nk_red, grid_ind=grid_ind, kweights = kweights)
        nproj = size(proj)[3]
#        println("proj $proj names $names pwan $pwan")
        pdos = zeros(length(energies),nproj, nspin)
    end


    

    KW = repeat(kweights, 1, tbc.tb.nwan)
                
    DOS = zeros(length(energies), nspin)

#    println("size(KW) $(size(KW))")
#    println("size(vals) $(size(vals))")
#    println("size(DOS) $(size(DOS))")
    
    for spin = 1:nspin
        for (c,e) in enumerate(energies)
            DOS[c,spin] = sum(KW .*  exp.( -0.5 * (vals[:,:,spin] .- e).^2 / smearing^2 ) ) / 2.0 * nk
        end

        if do_proj
            println(size(proj))
            for i = 1:nproj
                for (c,e) in enumerate(energies)
                    pdos[c, i,spin] = sum(proj[:,:,i,spin] .*  exp.( -0.5 * (vals[:,:,spin] .- e).^2 / smearing^2 ) )
                end
            end
        end
    end
    DOS = DOS / smearing / (2.0*pi)^0.5 / nk
    pdos = pdos / smearing / (2.0*pi)^0.5 / nk
    

    println("Int DOS " , sum(DOS) * (energies[2]-energies[1]) )

    ind = energies .< 0

    println("Int DOS occ " , sum(DOS[ind,:]) * (energies[2]-energies[1]) )
    for p in 1:nproj
        println("Int pDOS occ $p : " , sum(pdos[ind, p,:]) * (energies[2]-energies[1]) )
    end

    for p in 1:nproj
        println("Int pDOS all $p : " , sum(pdos[:, p,:]) * (energies[2]-energies[1]) )
    end

    energies = convert_energy(energies)
    DOS = convert_dos(DOS)
    pdos = convert_dos(pdos)
    
    plot_dos(energies, DOS, pdos, names, do_display=do_display)


    
    return energies,DOS, pdos, names, proj

    
end


function gaussian_dos(tbcK::tb_crys_kspace; smearing=0.03, npts=missing, proj_type=missing, do_display=true)

    return dos(tbcK; smearing=smearing, npts=npts, proj_type=proj_type, do_display=do_display)
    
end




function dos(tbcK::tb_crys_kspace; smearing=0.03, npts=missing, proj_type=missing, do_display=true, align=:efermi)


    @time bandenergy, eden, vects, vals, efermi, error_flag = get_energy_electron_density_kspace(tbcK.tb, tbcK.nelec, smearing=0.01)


#    etot, efermi, vals, vects, hk3, sk3 = calc_energy_fft(tbc, grid=grid, smearing=smearing, return_more_info=true)

#    println("precheck ", sum(vects[1,:,:]' * sk3[:,:,1,1,1] * vects[1,:,:]))
    
    #prelim
    if align == :efermi
        vals = vals .- efermi
    else
        vbm, cbm = find_vbm_cbm(vals, efermi)
        vals = vals .- vbm
    end


    nk = size(vals)[1]
    nspin = size(vects)[2]
    
    vmin = minimum(vals) - 0.1
    vmax = maximum(vals) + 0.1

    println("v $vmin $vmax ")

#    vmax = min(maximum(vals), 5.0)
    r = vmax - vmin

    if ismissing(npts)
        npts = Int64(round(r * 100 ))
    end
    
    energies = collect(vmin - r*0.02 : r*1.04 / npts    : vmax + r*0.02 + 1e-7)
    
    dos = zeros(length(energies), nspin)

    if ismissing(proj_type)  
        do_proj=true
        proj, names, pwan =  projection(tbcK, vects, tbcK.tb.Sk, ptype=proj_type)
        nproj = size(proj)[3]
        
        pdos = zeros(length(energies),nproj, nspin)
        
    elseif proj_type == "none" ||  proj_type == :none
        do_proj=false
        nproj=0
        pdos=missing
        names=missing
    else
        #    if ismissing(proj_type) || !(  proj_type != "none"  || proj_type != :none)
        do_proj=true
        proj, names, pwan =  projection(tbcK, vects, tbcK.tb.Sk, ptype=proj_type)
        nproj = size(proj)[3]
        
        pdos = zeros(length(energies),nproj, nspin)

    end
    
    W = repeat(tbcK.tb.kweights, 1, tbcK.tb.nwan)

    
    for spin = 1:nspin
        for (c,e) in enumerate(energies)
            dos[c] = sum(exp.( -0.5 * (vals[:,:, spin] .- e).^2 / smearing^2 ) .* W ) 
        end
    end
    dos = dos / smearing / (2.0*pi)^0.5 / 2

    if do_proj
        for spin = 1:nspin
            for i = 1:nproj
                for (c,e) in enumerate(energies)
                    pdos[c, i, spin] = sum(proj[:,:,i, spin] .*  exp.( -0.5 * (vals[:,:, spin] .- e).^2 / smearing^2 ) .* W)  
                end
            end
        end
        pdos = pdos / smearing / (2.0*pi)^0.5 / 2
    end

    

    println("Int DOS " , sum(dos) * (energies[2]-energies[1]) )

    ind = energies .< 0

    for spin = 1:nspin
        println("$spin Int DOS occ " , sum(dos[ind, spin]) * (energies[2]-energies[1]) )
        for p in 1:nproj
            println("$spin Int pDOS occ $p : " , sum(pdos[ind, p, spin]) * (energies[2]-energies[1]) )
        end
    end    
    energies = convert_energy(energies)
    dos = convert_dos(dos)
    pdos = convert_dos(pdos)
    
    plot_dos(energies, dos, pdos, names, do_display=do_display)


    
    return energies, dos, pdos, names

    
end

"""
    function setup_tetra(grid)

Setup simple tetrahedron method
based on tetra.f90 in QE. Simple tetra method.

"""
function setup_tetra(grid)

    ntetra = prod(grid) * 6

    tetra_map = zeros(Int64, 4, ntetra)

    kpts = zeros(Int64, prod(grid), 3)
    kpts_dict = Dict()
    
    c=0    
    for k1 = 0:grid[1]-1
        for k2 = 0:grid[2]-1
            for k3 = 0:grid[3]-1
                c+=1
                kpts[c,:] = [k1+1, k2+1, k3+1]
                kpts_dict[ [k1+1, k2+1, k3+1]] = c

            end
        end
    end

    c=1
    
    for k1 = 0:grid[1]-1
        for k2 = 0:grid[2]-1
            for k3 = 0:grid[3]-1
                
                k1a = (k1) % grid[1] + 1
                k2a = (k2) % grid[2] + 1
                k3a = (k3) % grid[3] + 1 

                k1b = (k1+1) % grid[1] + 1
                k2b = (k2+1) % grid[2] + 1
                k3b = (k3+1) % grid[3] + 1 

                #corners of cube (eight)
                n1 = kpts_dict[[k1a, k2a, k3a]]
                n2 = kpts_dict[[k1a, k2a, k3b]]
                n3 = kpts_dict[[k1a, k2b, k3a]]
                n4 = kpts_dict[[k1a, k2b, k3b]]
                n5 = kpts_dict[[k1b, k2a, k3a]]
                n6 = kpts_dict[[k1b, k2a, k3b]]
                n7 = kpts_dict[[k1b, k2b, k3a]]
                n8 = kpts_dict[[k1b, k2b, k3b]]

                #six tetra per cube
                tetra_map[1, c] = n1
                tetra_map[2, c] = n2
                tetra_map[3, c] = n3
                tetra_map[4, c] = n6

                tetra_map[1, c+1] = n2
                tetra_map[2, c+1] = n3
                tetra_map[3, c+1] = n4
                tetra_map[4, c+1] = n6
                
                tetra_map[1, c+2] = n1
                tetra_map[2, c+2] = n3
                tetra_map[3, c+2] = n5
                tetra_map[4, c+2] = n6
                
                tetra_map[1, c+3] = n3
                tetra_map[2, c+3] = n4
                tetra_map[3, c+3] = n6
                tetra_map[4, c+3] = n8
                
                tetra_map[1, c+4] = n3
                tetra_map[2, c+4] = n6
                tetra_map[3, c+4] = n7
                tetra_map[4, c+4] = n8
                
                tetra_map[1, c+5] = n3
                tetra_map[2, c+5] = n5
                tetra_map[3, c+5] = n6
                tetra_map[4, c+5] = n7
                
                c+=6

            end
        end
    end

    return kpts, kpts_dict, tetra_map

end


"""
    function dos_tetra(tbc::tb_crys; grid=missing, npts=missing, proj_type=missing, do_display=true)

DOS, using tetrahedral integration

- `grid` is the k-point grid. Defaults to 1.6 times the default grid for energy integration.
- `npts` is number of energies
- `proj_type` can be `"none"`, `"atomic"`, or `"orbital"`. Defaults to `atomic` if more than one atom type.
- `do_display=false` will suppress plotting

`return energies, dos, projected_dos, pdos_names`
"""
function dos_tetra(tbc::tb_crys; grid=missing, npts=missing, proj_type=missing, do_display=true)

    if ismissing(grid)
        grid = get_grid(tbc.crys)
        grid = Int64.(round.(grid * 1.6))
        println("grid $grid")
    end


    @time etot, efermi, vals, vects,hk3, sk3 = calc_energy_fft(tbc, grid=grid, return_more_info=true)

#    println("vects ", size(vects))
#    println(vects[1,:,:])
#    println("vals ", size(vals))
#    println(vals[1, :])
#    println("sk ", size(sk3))
#    println(sk3[:,:,1,1,1])

    nspin = size(vects)[2]
    
    begin
    
    vals = vals .- efermi
    
    nk = size(vals)[1]
    
    vmin = minimum(vals)
    vmax = maximum(vals)
    vmax = min(maximum(vals), 5.0)
    r = vmax - vmin

    if ismissing(npts)
        npts = 500
    end
        
    println("vmin $vmin $vmax $vmax r $r")
        
    energies = collect(vmin - r*0.05 : r*1.02 / npts    : vmax + r*0.05 + 1e-7)

    println("le ", length(energies))
        
    dos = zeros(size(energies)[1], nspin )
    dos_id = zeros( length(energies), nthreads())

    if ismissing(proj_type) ||  (proj_type != "none"   && proj_type != :none)
        do_proj=true
        #println("Projection")
        proj, names, pwan =  projection(tbc, vects, sk3, grid, ptype=proj_type, use_sym=use_sym)
        nproj = size(proj)[3]

        pdos = zeros(length(energies),nproj, nspin)
        pdos_id = zeros(length(energies),nproj, nthreads())

    else
        do_proj=false
        nproj=0
        pdos=missing
        names=missing
    end

        #println("setup tetra")
        kpts, kpts_dict, tetra_map = setup_tetra(grid)

    #    e=zeros(4)

    ntetra = nk*6
    norm = Float64(1/ntetra)
    

    #    ex=zeros(4)
#    px=zeros(4)    
#    ktet = zeros(Int64, 4)

    
    eps = 1e-7

    #println("calc dos tetra")

    end
        
    range = 1:length(energies)

    for spin = 1:nspin
        
        @threads for nt = 1: ntetra
            id = threadid()
            #id = 1
            ex=zeros(4)
            px=zeros(4)    
            e=zeros(4, tbc.tb.nwan)
            p=zeros(4, tbc.tb.nwan, nproj)
            ktet = zeros(Int64, 4)
            wt = zeros(4)


            ktet[:] .= tetra_map[:, nt]
            
            e[:,:] .= vals[ktet, :, spin]        

            if do_proj
                p[:,:, :] .= proj[ktet,:, :, spin]
                f14 = 0.0
                f24 = 0.0
                f35 = 0.0
                G = 0.0
                f13 = 0.0
                f31 = 0.0 
                f14 = 0.0 
                f41 =  0.0
                f23 = 0.0
                f32 = 0.0
                f24 = 0.0
                f42 = 0.0
                
                f12 =  0.0
                f21 = 0.0
                

            end
            
            for nband = 1:tbc.tb.nwan

                ex[:] = e[:,nband]


                if do_proj==false  ####
                    
                    sort!(ex)
                    


                    inds = (ex[1] .< energies .< ex[4] )
                    
                    for c in range[inds]

                        en = energies[c]
                        
                        if en < ex[4] && en >= ex[3]
                            
                            dos_id[c, id] += 3.0 * (ex[4] - en)^2 / (ex[4] - ex[1] + eps) / (ex[4] - ex[2] + eps) / (ex[4] - ex[3] + eps)
                            
                        elseif en < ex[3] && en >= ex[2]
                            
                            dos_id[c, id] += 1.0 / (ex[3] - ex[1] + eps) / (ex[4] - ex[1] + eps) * (3.0 * (ex[2] - ex[1]) + 6.0 * (en - ex[2]) - 3.0 * (ex[3] - ex[1] + ex[4] - ex[2]) / (ex[3]-ex[2] + eps) / (ex[4] - ex[2] + eps) * (en -ex[2])^2)
                            
                        elseif en < ex[2] && en >= ex[1]
                            
                            dos_id[c, id] += 3.0 * (en-ex[1])^2 / (ex[2] - ex[1] + eps) / (ex[3] - ex[1] + eps) / (ex[4] - ex[1] + eps)
                            
                        end
                        
                    end

                elseif do_proj == true ####

                    
                    perm = sortperm(ex)
                    ex[:] = ex[perm]

                    inds = (ex[1] .< energies .< ex[4] )
                    
                    for c in range[inds]
                        #                for c in 1:length(energies)
                        en = energies[c]

                        
                        if en < ex[4] && en >= ex[1]
                            
                            if en < ex[4] && en >= ex[3]
                                
                                f14 = (en-ex[4])/(ex[1]-ex[4] - eps)
                                f24 = (en-ex[4])/(ex[2]-ex[4] - eps)
                                f34 = (en-ex[4])/(ex[3]-ex[4] - eps)
                                
                                G = 3.0 * f14 * f24 * f34 / (ex[4] - en + eps)
                                
                                @inbounds wt[1] = f14 / 3.0
                                @inbounds wt[2] = f24 / 3.0
                                @inbounds wt[3] = f34 / 3.0
                                @inbounds wt[4] = (3.0 - f14 - f24 - f34) / 3.0
                                
                            elseif en < ex[3] && en >= ex[2]
                                
                                f13 = (en-ex[3])/(ex[1]-ex[3] - eps)
                                f31 = 1.0 - f13
                                f14 = (en-ex[4])/(ex[1]-ex[4] - eps)
                                f41 = 1.0 - f14
                                f23 = (en-ex[3])/(ex[2]-ex[3] - eps)
                                f32 = 1.0 - f23
                                f24 = (en-ex[4])/(ex[2]-ex[4] - eps)
                                f42 = 1.0 - f24
                                
                                G   =  3.0 * (f23*f31 + f32*f24)
                                
                                @inbounds wt[1]  =  f14 / 3.0 + f13*f31*f23 / (G + eps)
                                @inbounds wt[2]  =  f23 / 3.0 + f24*f24*f32 / (G + eps)
                                @inbounds wt[3]  =  f32 / 3.0 + f31*f31*f23 / (G + eps)
                                @inbounds wt[4]  =  f41 / 3.0 + f42*f24*f32 / (G + eps)
                                
                                G   =  G / (ex[4]-ex[1] + eps)
                                
                                
                            elseif en < ex[2] && en >= ex[1]
                                
                                
                                f12 = (en-ex[2])/(ex[1]-ex[2] - eps)
                                f21 = 1.0 - f12
                                f13 = (en-ex[3])/(ex[1]-ex[3] - eps)
                                f31 = 1.0 - f13
                                f14 = (en-ex[4])/(ex[1]-ex[4] - eps)
                                f41 = 1.0 - f14
                                
                                G  =  3.0 * f21 * f31 * f41 / (en-ex[1] + eps)
                                
                                @inbounds wt[1] =  (f12 + f13 + f14) / 3.0
                                @inbounds wt[2] =  f21  / 3.0
                                @inbounds wt[3] =  f31  / 3.0
                                @inbounds wt[4] =  f41  / 3.0
                                
                                
                            end
                            
                            
                            
                            @inbounds dos_id[c, id] += sum(wt) * G
                            
                            for projind = 1:nproj
                                px[:] = @view p[perm,nband, projind]
                                @inbounds pdos_id[c,projind, id] += ( px' * wt )  * G
                            end
                            
                        end
                        
                    end

                end #### if do_proj

            end

        end

    
        dos[:, spin] = sum(dos_id, dims=2) * norm
        dos[:,spin] = max.(dos[:,spin], 0.0)
        pdos[:,:,spin] = sum(pdos_id, dims=3) * norm
        
    end
        
    #println("Int DOS " , sum(dos) * (energies[2]-energies[1]) )
    correct =  tbc.tb.nwan / ( sum(dos) * (energies[2]-energies[1]) ) 

    dos = dos * correct
        
    if nproj > 0
        pdos = max.(pdos, 0.0)
        
#        pdos[:,:] = sum(pdos_id, dims=3) * norm

        pdos = pdos * correct

        #println("Int PDOS " , sum(dos, dims=1) * (energies[2]-energies[1]) )

        t = sum(pdos, dims=1) * (energies[2]-energies[1])
        fix = pwan' ./ t
        pdos = pdos .* fix
        pdos2 = zeros(size(pdos)[1],size(pdos)[2], nspin)
        
    else
        pdos2=missing
    end

    
    #simple smoothing
    dos2 = zeros(size(dos)[1], nspin)


    for spin = 1:tbc.nspin
        for i = 3:(length(dos[:,1])-2)
            dos2[i, spin] =  0.125*dos[i-2, spin] + 0.25*dos[i-1, spin] + 0.25*dos[i, spin] + 0.25*dos[i+1, spin] + 0.125*dos[i+2, spin]
            if nproj>0
                #            pdos2[i,:] = 0.25*pdos[i-1,:] + 0.5*pdos[i,:] + 0.25*pdos[i+1,:]
                pdos2[i,:, spin] = 0.125*pdos[i-2,:, spin] + 0.25*pdos[i-1,:, spin] + 0.25*pdos[i,:, spin] + 0.25*pdos[i+1,:, spin] + 0.125*pdos[i+2,:, spin]            
            end
        end
    end
    
    if nproj > 0
        pdos2 = max.(pdos2, 0.0)
        pdos2 = convert_dos(pdos2)
    end
    
    energies = convert_energy(energies)
    dos2 = convert_dos(dos2)

    plot_dos(energies, dos2, pdos2, names, do_display=do_display)

    
    return energies, dos2, pdos2, names
    
end

"""
    function plot_dos(energies, dos, pdos, names; filename=missing, do_display=true)

Does the actual DOS plotting, called by `dos` or `gaussian_dos`
"""
function plot_dos(energies, dos, pdos, names; filename=missing, do_display=true)

    if ismissing(names)
        plot(legend=false, grid=false, framestyle=:box)
    else 
        plot(legend=true, grid=false, framestyle=:box)
    end

    nspin = size(dos)[2]

    if nspin == 1
    
        plot!(energies, dos[:,1], color="black", lw=4, label="Total", legend=:topleft, xtickfontsize=12,ytickfontsize=12, legendfontsize=10)

        if !ismissing(names)
            colors = ["blue", "orange", "green", "magenta", "cyan", "red", "yellow"]
            if !ismissing(pdos)
                for i in 1:size(pdos)[2]
                    color = colors[i%7+1]
                    #            println("i $i $color", names[i])
                    plot!(energies, pdos[:,i,1], color=color, lw=3, label=names[i])
                end
            end
        end    

    elseif nspin == 2

        plot!(energies, dos[:,1], color="black", lw=4, label="Total", legend=:topleft, xtickfontsize=12,ytickfontsize=12, legendfontsize=10)

        if !ismissing(names)
            colors = ["blue", "orange", "green", "magenta", "cyan", "red", "yellow"]
            if !ismissing(pdos)
                for i in 1:size(pdos)[2]
                    color = colors[i%7+1]
                    #            println("i $i $color", names[i])
                    plot!(energies, pdos[:,i,1], color=color, lw=3, label=names[i])
                end
            end
        end    

        plot!(energies, -dos[:,2], color="black", lw=4, label=false)

        if !ismissing(names)
            colors = ["blue", "orange", "green", "magenta", "cyan", "red", "yellow"]
            if !ismissing(pdos)
                for i in 1:size(pdos)[2]
                    color = colors[i%7+1]
                    #            println("i $i $color", names[i])
                    plot!(energies, -pdos[:,i,2], color=color, lw=3, label=false)
                end
            end
        end    
    end

        
    
    if global_energy_units == "eV"
        ylabel!("DOS  ( 1 / eV )", guidefontsize=16)
        xlabel!("Energy - E_F ( eV )", guidefontsize=16)
    else
        ylabel!("DOS  ( 1 / Ryd. )", guidefontsize=16)
        xlabel!("Energy - E_F ( Ryd. )", guidefontsize=16)
    end

    xl1 = minimum(energies)
    xl2 = min(10.0, maximum(energies))
#    xl2 = ( maximum(energies))

    xlims!(xl1, xl2)

    if nspin == 1
        ylims!(0.0, maximum(dos) * 1.1)
    elseif nspin == 2
        ylims!(-maximum(dos) * 1.1, maximum(dos) * 1.1)
    end
    
    if do_display && ! no_display
        display(plot!([0,0], [-maximum(dos) * 1.1, maximum(dos) * 1.1], color="black", linestyle=:dash, label=""))
    else

        plot!([0,0], [-maximum(dos) * 1.1, maximum(dos) * 1.1], color="black", linestyle=:dash, label="")

    end        

    if !ismissing(filename)
        savefig(filename)
    end
    
end


"""
    function plot_dos_flip(energies, dos, pdos, names; filename=missing, do_display=true, yrange=[-6, 4.0])

For plotting with x and y axes reversed, useful for bandstucture plots
"""
function plot_dos_flip(energies, dos, pdos, names; filename=missing, do_display=true, yrange=[-6, 4.0])

    
    if ismissing(names)
        plot(legend=false, grid=false, framestyle=:box)
    else 
        plot(legend=true, grid=false, framestyle=:box)
    end

    nspin = size(dos)[2]

    if nspin == 1
    
        plot!( dos[:,1], energies, color="black", lw=4, label="Total", legend=:topright, xtickfontsize=12,ytickfontsize=12, legendfontsize=10)

        if !ismissing(names)
            colors = ["blue", "orange", "green", "magenta", "cyan", "red", "yellow"]
            if !ismissing(pdos)
                for i in 1:size(pdos)[2]
                    color = colors[i%7+1]
                    #            println("i $i $color", names[i])
                    plot!( pdos[:,i,1], energies, color=color, lw=3, label=names[i])
                end
            end
        end    

    elseif nspin == 2

        plot!( dos[:,1], energies, color="black", lw=4, label="Total", legend=:topright, xtickfontsize=12,ytickfontsize=12, legendfontsize=10)

        if !ismissing(names)
            #colors = ["blue", "orange", "green", "magenta", "cyan", "red", "yellow"]
            colors = ["blue", "magenta", "cyan", "orange", "green",  "red", "yellow"]
            if !ismissing(pdos)
                for i in 1:size(pdos)[2]
                    color = colors[i%7+1]
                    #            println("i $i $color", names[i])
                    plot!(pdos[:,i,1], energies,  color=color, lw=3, label=names[i])
                    p = ylims!(yrange[1], yrange[2])
                    
                end
            end
        end    

        plot!( -dos[:,2], energies, color="black", lw=4, label=false)

        if !ismissing(names)
            colors = ["blue", "magenta", "cyan", "orange", "green",  "red", "yellow"]
            if !ismissing(pdos)
                for i in 1:size(pdos)[2]
                    color = colors[i%7+1]
                    #            println("i $i $color", names[i])
                    plot!( -pdos[:,i,2], energies, color=color, lw=3, label=false)
                    p = ylims!(yrange[1], yrange[2])
                    
                end
            end
        end    
    end

        
    
    if global_energy_units == "eV"
        xlabel!("DOS  ( 1 / eV )", guidefontsize=10)
        #ylabel!("Energy - E_F ( eV )", guidefontsize=16)
    else
        xlabel!("DOS  ( 1 / Ryd. )", guidefontsize=10)
        #ylabel!("Energy - E_F ( Ryd. )", guidefontsize=16)
    end

    xl1 = minimum(energies)
    xl2 = min(10.0, maximum(energies))
#    xl2 = ( maximum(energies))

    p = ylims!(yrange[1], yrange[2])

    println("DOS flit range ", yrange)
    
    if nspin == 1
        xlims!(0.0, maximum(dos) * 1.1)
    elseif nspin == 2
        xlims!(-maximum(dos) * 1.1, maximum(dos) * 1.1)
    end

    yticks!(Float64[],String[])
    
    if do_display && ! no_display
        display(plot!( [-maximum(dos) * 1.1, maximum(dos) * 1.1], [0,0], color="black", linestyle=:dot, label=""))
    else

        plot!( [-maximum(dos) * 1.1, maximum(dos) * 1.1], [0,0], color="black", linestyle=:dot, label="")

    end        

    return p
    
#    if !ismissing(filename)
#        savefig(filename)
    #    end
    
end


function dos(dft::dftout; smearing=0.03, npts=missing, proj_type=missing, do_display=true, align=:efermi)

    return dos(dft.bandstruct, smearing=smearing, do_display=do_display, align=align)
end

function dos(bs::bandstructure; smearing=0.03, npts=missing, proj_type=missing, do_display=true, align=:efermi)
    
    
    #prelim
    if align == :efermi
        vals = bs.eigs .- bs.efermi
    else
        vbm, cbm = find_vbm_cbm(bs.eigs, bs.efermi)
        vals = bs.eigs .- vbm
    end
    
    

    nk = size(vals)[1]
    nspin = size(vals)[3]
    
    vmin = minimum(vals) - 0.1
    vmax = maximum(vals) + 0.1

    r = vmax - vmin

    if ismissing(npts)
        npts = Int64(round(r * 100 ))
    end
    
    energies = collect(vmin - r*0.02 : r*1.04 / npts    : vmax + r*0.02 + 1e-7)
    
    dos = zeros(length(energies), nspin)

    W = repeat(bs.kweights, 1, bs.nbnd)

    for spin = 1:nspin
        for (c,e) in enumerate(energies)
            dos[c, spin] = sum(exp.( -0.5 * (vals[:,:, spin] .- e).^2 / smearing^2 ) .* W ) 
        end
    end
    dos = dos / smearing / (2.0*pi)^0.5 / 2

    println("Int DOS " , sum(dos) * (energies[2]-energies[1]) )

    ind = energies .< 0

    for spin = 1:nspin
        println("$spin Int DOS occ " , sum(dos[ind, spin]) * (energies[2]-energies[1]) )
        
    end    
    energies = convert_energy(energies)
    dos = convert_dos(dos)
    
    plot_dos(energies, dos, missing, [], do_display=do_display)


    
    return energies, dos * bs.nspin

    
end

function compare_dos_dft(tbc::tb_crys, dft::dftout; smearing = 0.03, do_display=true, align=:vbm)

    grid = dft.bandstruct.kgrid
    energies_dft, dos_dft = dos(dft, do_display=false, align=align)
    energies_tb , dos_tb  = dos(tbc, grid = grid, do_display=false, align=align)
    
    if ismissing(names)
        plot(legend=false, grid=false, framestyle=:box)
    else 
        plot(legend=true, grid=false, framestyle=:box)
    end

    if dft.nspin == 1
        plot!(energies_dft, dos_dft[:,1], color="black", lw=4, label="DFT", legend=:topleft, xtickfontsize=12,ytickfontsize=12, legendfontsize=10)
        plot!(energies_tb, dos_tb[:,1], color="red", lw=4, linestyle=:dash, label="TB3", legend=:topleft, xtickfontsize=12,ytickfontsize=12, legendfontsize=10)
    else

        plot!(energies_dft, dos_dft[:,1], color="black", lw=4, label="DFT", legend=:topleft, xtickfontsize=12,ytickfontsize=12, legendfontsize=10)
        plot!(energies_tb, dos_tb[:,1], color="red", lw=4, linestyle=:dash, label="TB3", legend=:topleft, xtickfontsize=12,ytickfontsize=12, legendfontsize=10)
        plot!(energies_dft, -dos_dft[:,2], color="black", lw=4, label="DFT", legend=:topleft, xtickfontsize=12,ytickfontsize=12, legendfontsize=10)
        plot!(energies_tb, -dos_tb[:,2], color="red", lw=4, linestyle=:dash, label="TB3", legend=:topleft, xtickfontsize=12,ytickfontsize=12, legendfontsize=10)

    end

    if global_energy_units == "eV"
        ylabel!("DOS  ( 1 / eV )", guidefontsize=16)
        if align == :efermi
            xlabel!("Energy - E_F ( eV )", guidefontsize=16)
        else
            xlabel!("Energy - VBM ( eV )", guidefontsize=16)            
        end
    else
        ylabel!("DOS  ( 1 / Ryd. )", guidefontsize=16)
        if align == :efermi
            xlabel!("Energy - E_F ( Ryd. )", guidefontsize=16)
        else
            xlabel!("Energy - VBM ( Ryd. )", guidefontsize=16)            
        end

    end

        
    if dft.nspin == 1
        ylims!(0.0, max(maximum(dos_dft), maximum(dos_tb)) * 1.1)
    elseif dft.nspin == 2
        d = max(maximum(abs.(dos_dft)), maximum(abs.(dos_tb)))
        ylims!(-d*1.1, d*1.1)
    end

    d = max(maximum(abs.(dos_dft)),maximum(abs.(dos_tb)))
    if do_display && ! no_display
        display(plot!([0,0], [-d * 1.1, d * 1.1], color="black", linestyle=:dash, label=""))
    else

        plot!([0,0], [- d * 1.1, d * 1.1], color="black", linestyle=:dash, label="")

    end        
                    
        
    
end


end #end module

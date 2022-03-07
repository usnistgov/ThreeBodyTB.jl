###include("Crystal1563.jl")
###include("DFToutMod.jl")
#using XMLDict



####################### Wannier90 specific 
module TB
"""
Scripts to run tight-binding from wannier90
"""

using LinearAlgebra
using SpecialFunctions
using EzXML
using XMLDict
using GZip
using Printf
using Plots
using FFTW
#using JLD
using Base.Threads

#include("Atomdata.jl")
using ..Atomdata:atoms
#include("Commands.jl")

using ..DFToutMod:bandstructure
using ..DFToutMod:dftout
using ..AtomicMod:atom
using ..CrystalMod:crystal
using ..CrystalMod:makecrys
using ..Utility:arr2str
using ..Utility:str_w_spaces
using ..Utility:parse_str_ARR_float
using ..BandTools:calc_fermi
using ..BandTools:band_energy
using ..BandTools:gaussian
using ..BandTools:smearing_energy

using ..CrystalMod:get_grid
using ..Ewald:electrostatics_getgamma
using ..CrystalMod:orbital_index

export tb
export tb_crys
export tb_k
export tb_crys_kspace
export read_tb_crys
export read_tb_crys_kspace
export write_tb_crys
export write_tb_crys_kspace
export make_tb_crys
export make_tb_crys_kspace
export make_tb
export make_tb_k
export load_hr_dat
export Hk

export calc_bands



symbol_dict = Dict()
#        1   2    3    4    5     6    7     8      9   
list = [:s, :pz, :px, :py, :dz2,:dxz,:dyz,:dx2_y2,:dxy,:fz3,:fxz2,:fyz2,:fz_x2y2,:fxyz,:fx_x2_3y2,:fy_3x2_y2]

for i = 1:16
    symbol_dict[list[i]] = i
end


    
"""
    mutable struct tb{T}

Holds key tight-binding information in real-space. Like `_hr.dat` file from Wannier90. Also part of the `tb_crys` object.

# Holds
- `H::Array{Complex{T},4}` Hamiltonian. `nwan`×`nwan`×`nr`×`nspin`
- `ind_array::Array{Int64,3}` `nr`×3 , holds the r-space supercells of the TB object.
- `r_dict::Dict` keys are three Ints like [0,0,0], returns the corresponding `ind_array` index.
- `nwan::Int` Number of orbitals (generalized wannier functions).
- `nspin::Int` Number of spins (2=magnetic)
- `nr::Int64` number of R-space supercells.
- `nonorth::Bool` :  `true` if non-orthogonal. Almost always `true` in this code.
- `S::Array{Complex{T},3}` : Overlap matrix, organized like `H`
- `scf::Bool` equal to `true` if requires self-consistency (usually `true` for fit `tb`, `false` for direct from DFT)
- `h1::Array{T,3}` Has the term determined by scf calculations, if calculated already.
"""
mutable struct tb{T}

    #    H::Array{Complex{Float64},3}
    H::Array{Complex{T},4}    
    ind_arr::Array{Int64,2}
    r_dict::Dict
    nwan::Int
    nr::Int
    nspin::Int              
    nonorth::Bool
    #    S::Array{Complex{Float64},3}
    S::Array{Complex{T},3}    
    scf::Bool
    h1::Array{T,2} #scf term

end

Base.show(io::IO, h::tb) = begin
    println(io)
    nwan=h.nwan
    nr=h.nr
    nonorth = h.nonorth
    scf = h.scf
    nspin=h.nspin
    println(io, "tight binding real space object; nwan = $nwan, nr = $nr, nonorth = $nonorth, scf = $scf, nspin = $nspin" )
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
mutable struct tb_crys{T}

    tb::tb
    crys::crystal
    nelec::Float64
    dftenergy::Float64
    scf::Bool
    gamma::Array{T, 2}
    eden::Array{Float64,2}
    within_fit::Bool
    energy::Float64
    efermi::Float64
    nspin::Int64
end



Base.show(io::IO, x::tb_crys) = begin
    println(io)
    println(io, "tb_crys object ")
    println(io)    
    println(io, x.crys)
    println(io)
    println(io, "nelec: ", x.nelec, "; nspin: ", x.nspin)
    println(io, "within_fit: ", x.within_fit,"  ; scf: ", x.scf)
    println("calculated energy: ", round(x.energy*1000)/1000)
    println("efermi: ", round(x.efermi*1000)/1000)
    dq = get_dq(x)
    println(io, "charges: ", round.(dq * 100)/100)
    println(io)
    println(io, x.tb)    
    println(io)
    
end   


"""
    mutable struct tb_k{T}

Tight binding object in k-space. Can be from direct import of DFT band
structure using atomic proj, or from fft'ed tb object. Similar to real-space version.

# Holds
- `Hk::Array{Complex{T},4}` Hamiltonian in k-space
- `K::Array{Float64,2}` K point array. In fractional coordinates of BZ.
- `kweights::Array{Float64,1}` K point weights.
- `k_dict::Dict` Dictionary from k-point like [0,0,0] to index.
- `nwan::Int64` Number of orbitals / generalized wannier functions.
- `nk::Int64` Number of k-points.
- `nspin::Int64` Number of spins (2=magnetic)
- `nonorth::Bool` 
- `Sk::Array{Complex{T},3}` Overlap matrix.
- `scf::Bool` needs self-consistency?
- `h1::Array{T,3} #scf term` holds scf term if present
- `grid::Array{Int64,1}` dimensions of k-point grid, from regular grid like `[8,8,8]`


""" 
mutable struct tb_k{T}


    Hk::Array{Complex{T},4}    
    K::Array{Float64,2}
    kweights::Array{Float64,1}
    k_dict::Dict
    nwan::Int64
    nk::Int64
    nspin::Int64
    nonorth::Bool
    Sk::Array{Complex{T},3}    
    scf::Bool
    h1::Array{T,2} #scf term
    grid::Array{Int64,1}

end


Base.show(io::IO, h::tb_k) = begin
    println(io)
    nwan=h.nwan
    nk=h.nk
    nonorth = h.nonorth
    scf = h.scf
    grid = h.grid
    nspin = h.nspin
    println(io, "tight binding k-space object; nwan = $nwan, nk = $nk, nonorth = $nonorth, scf = $scf, grid = $grid, nspin = $nspin" )
    println(io)
    
end   


"""
   mutable struct tb_crys_kspace{T}

Hold k-point tight binding and crystal structure. Similar to `tb_crys`

# Holds
- `tb::tb_k`
- `crys::crystal`
- `nelec::Float64`
- `nspin::Int64`
- `dftenergy::Float64`
- `scf::Bool`
- `gamma::Array{T, 2}`
- `eden::Array{Float64,2}`
"""
mutable struct tb_crys_kspace{T}

    tb::tb_k
    crys::crystal
    nelec::Float64
    nspin::Int64
    dftenergy::Float64
    scf::Bool
    gamma::Array{T, 2}
    eden::Array{Float64,2}
    energy::Float64

end



Base.show(io::IO, x::tb_crys_kspace) = begin
    println(io)
    println(io, x.crys)
    println(io)
    println(io, "nelec: ", x.nelec, "; nspin: ", x.nspin)
    println(io)
    println(io, x.tb)    
    println(io)
    
end   


# old utility function
#function get_el(tbc::tb_crys, n1,n2,na,nb,r)
#
#    ind = 1
#    try
#        ind=tbc.tb.r_dict[r]
#        
#    catch
#        try
#            ind=tbc.tb.r_dict[r']
#        catch
#            println("missing")
#            return missing, missing, missing, missing
#        end
#    end
#    
#    R = tbc.tb.ind_arr[ind,:]
#    dR = tbc.crys.A'*(-tbc.crys.coords[na,:] + tbc.crys.coords[nb,:] + R) 
#    dist = sum(dR.^2)^0.5
#
#    
#    return tbc.tb.H[n1,n2,ind], tbc.tb.S[n1,n2,ind], dist, dR
#    
#end

#function get_el(tbc::tb_crys, n1,n2,r)
#
#    return get_el(tbc.tb, n1,n2,r)
#    
#end

"""
    function read_tb_crys(filename, tbc::tb_crys)

Reads and returns from `filename` a `tb_crys` object. See `write_tb_crys`

If cannot find `"filename"`, will look for `"filename.xml"`, `"filename.gz"`, `"filename.xml.gz"`

Can read gzipped files directly.

"""
function read_tb_crys(filename; directory=missing)
    """
    get tbc object from xml file, written by write_tb_crys (see below)
    """
    if ismissing(directory)
        println("read ", filename)
    else
        println("read ", directory*"/"*filename)
    end

    #    try
    f = missing
    if !ismissing(directory)
        filename=directory*"/"*filename
    end
    if !isfile(filename)
        if isfile(filename*".xml")
            filename=filename*".xml"
        elseif isfile(filename*".gz")
            filename=filename*".gz"
        elseif isfile(filename*".xml.gz")
            filename=filename*".xml.gz"
        end
    end
    if !isfile(filename)
        println("warning error read_tb_crys $filename $directory not found")
    else
        println("found $filename")
    end

    try
        f = gzopen(filename, "r")
    catch
        println("error opening $filename")
    end
    
    fs = read(f, String)
    close(f)

    d = xml_dict(fs)["root"]

    ##crys
    nat = parse(Int64, d["crystal"]["nat"])


    A = parse_str_ARR_float(d["crystal"]["A"])
    A = reshape(A, 3,3)'
    
    coords = parse_str_ARR_float(d["crystal"]["coords"])
    
    coords = reshape(coords, 3,nat)'
    
    types = split(d["crystal"]["types"])
    
    crys = makecrys(A,coords,types, units="Bohr")
    ##
    ##nelec
    nelec = parse(Float64,d["nelec"])
    nspin = 1
    if "nspin" in keys(d)
        nspin = parse(Int64,d["nspin"])
    end

    dftenergy = parse(Float64,d["dftenergy"])        

    if "energy" in keys(d)
        energy = parse(Float64,d["energy"])
    else
        energy = -999.0
    end

    if "efermi" in keys(d)
        efermi = parse(Float64,d["efermi"])
    else
        efermi = 0.0
    end

    

    #SCF stuff is optional
    scf = false
    gamma = missing
    eden = missing

    if "scf" in keys(d) scf = parse(Bool,d["scf"])  end
    
    if "gamma" in keys(d)
        gamma = parse_str_ARR_float(d["gamma"])
        s = Int64(round(sqrt(size(gamma)[1])))
        gamma = reshape(gamma, s, s)
        
    end
    

    ##tb
    
    nr = parse(Int64,d["tightbinding"]["nr"])
    nwan = parse(Int64,d["tightbinding"]["nwan"])
    
    nonorth = parse(Bool,d["tightbinding"]["nonorth"])
    
    ind_arr = parse_str_ARR_float(d["tightbinding"]["ind_arr"])
    ind_arr = reshape(ind_arr,3, nr)'
    
    r_dict = Dict()
    for i in 1:nr
        r_dict[copy(ind_arr[i,:])] = i
    end
    #    r_dict = Dict(d["tightbinding"]["r_dict"])

    if "scf" in keys(d["tightbinding"])
        tb_scf = parse(Bool,d["tightbinding"]["scf"])
    else
        tb_scf = false
    end
    if "h1" in keys(d["tightbinding"])
        h1 = parse_str_ARR_float(d["tightbinding"]["h1"])
        h1 = reshape(h1, nwan, nwan)
    else
        h1 = missing
    end
    if tb_scf == false
        h1 =  missing
    end


    if "eden" in keys(d)
        eden  = parse_str_ARR_float(d["eden"])
        eden = reshape(eden, nspin, nwan)
    end
    
    function readstr(st)

        H = zeros(Complex{Float64}, nspin, nwan,nwan,nr)
        S = zeros(Complex{Float64}, nwan,nwan,nr)        
        
        lines = split(st, "\n")
        for line in lines
            sp = split(line)
            if length(sp) == 7 || length(sp) == 9
                m = parse(Int64,sp[2])
                n = parse(Int64,sp[3])
                r = parse(Int64,sp[1])
                if nspin == 1
                    H[1, m,n,r] = parse(Float64,sp[4]) + im*parse(Float64,sp[5])
                    S[ m,n,r] = parse(Float64,sp[6]) + im*parse(Float64,sp[7])
                elseif nspin == 2
                    H[1, m,n,r] = parse(Float64,sp[4]) + im*parse(Float64,sp[5])
                    H[2, m,n,r] = parse(Float64,sp[6]) + im*parse(Float64,sp[7])
                    S[ m,n,r] = parse(Float64,sp[8]) + im*parse(Float64,sp[9])
                end
            end
        end
        return H, S
    end
        
    H,S = readstr(d["tightbinding"]["H"])

        #    println("read ")
        #    println("h1")
        #    println(h1)

    if nonorth
        tb = make_tb(H, ind_arr, r_dict, S, h1=h1)
    else
        tb = make_tb(H, ind_arr, r_dict, h1=h1)
    end    
    
    tbc = make_tb_crys(tb, crys, nelec, dftenergy, scf=scf, eden=eden, gamma=gamma, tb_energy=energy, fermi_energy=efermi)

    return tbc
        
    #    catch
        #        println("ERROR: Failed to read ", filename)
        #        return -1
        #    end
        
end



"""
    function read_tb_crys_kspace(filename; directory=missing)

Reads and returns from `filename` a `tb_crys_kspace` object. See `write_tb_crys_kspace`

If cannot find `"filename"`, will look for `"filename.xml"`, `"filename.gz"`, `"filename.xml.gz"`

Can read gzipped files directly.


"""
function read_tb_crys_kspace(filename; directory=missing)
"""
get tbc object from xml file, written by write_tb_crys (see below)
"""
    if ismissing(directory)
        println("read ", filename)
    else
        println("read ", directory*"/"*filename)
    end

    f = missing
    if !ismissing(directory)
        filename=directory*"/"*filename
    end
    if !isfile(filename)
        if isfile(filename*".xml")
            filename=filename*".xml"
        elseif isfile(filename*".gz")
            filename=filename*".gz"
        elseif isfile(filename*".xml.gz")
            filename=filename*".xml.gz"
        end
    end
    if !isfile(filename)
        println("warning error read_tb_crys $filename $directory not found")
    else
        println("read_tb_crys_kspace found $filename")
    end

    try
        f = gzopen(filename, "r")
    catch
        println("error opening $filename")
    end


    #    try
#    f = missing
#    try
#        if ismissing(directory)
#            f = gzopen(filename, "r")
#        else
#            f = gzopen(directory*"/"*filename, "r")
#        end
#    catch
#        filename=filename*".xml"
#            println("trying ", filename)
#        if ismissing(directory)
#            f = gzopen(filename, "r")
#        else
#            f = gzopen(directory*"/"*filename, "r")
#       end
#    end
    
    fs = read(f, String)
    close(f)

    d = xml_dict(fs)["root"]

    ##crys
    nat = parse(Int64, d["crystal"]["nat"])

    A = parse_str_ARR_float(d["crystal"]["A"])
    A = reshape(A, 3,3)'
    
    coords = parse_str_ARR_float(d["crystal"]["coords"])
    
    coords = reshape(coords, 3,nat)'
    
    types = split(d["crystal"]["types"])
    
    crys = makecrys(A,coords,types, units="Bohr")
    ##
    ##nelec
    nelec = parse(Float64,d["nelec"])
    dftenergy = parse(Float64,d["dftenergy"])        


    nspin = 1
    if "nspin" in keys(d)
        nspin = parse(Int64,d["nspin"])
    end


    #SCF stuff is optional
    scf = false
    gamma = missing
    eden = missing

    if "scf" in keys(d) scf = parse(Bool,d["scf"])  end
    
    if "gamma" in keys(d)
        gamma = parse_str_ARR_float(d["gamma"])
        s = Int64(round(sqrt(size(gamma)[1])))
        gamma = reshape(gamma, s, s)
        
    end
    if "eden" in keys(d)
        eden  = parse_str_ARR_float(d["eden"])
    end
    

    ##tb
    
    nk = parse(Int64,d["tightbinding"]["nk"])
    nwan = parse(Int64,d["tightbinding"]["nwan"])
    
    nonorth = parse(Bool,d["tightbinding"]["nonorth"])
        
    kind_arr = parse_str_ARR_float(d["tightbinding"]["kind_arr"])
    kind_arr = reshape(kind_arr,3, nk)'

    kweights = parse_str_ARR_float(d["tightbinding"]["kweights"])

    if "grid" in keys(d["tightbinding"])
#        grid = parse(Int64,d["tightbinding"]["grid"])
#        grid = parse_str_ARR_float(d["tightbinding"]["kweights"])
        grid = Int64.(parse_str_ARR_float(d["tightbinding"]["grid"]))

    else
        grid = [0,0,0]
    end
    
    k_dict = Dict()
    for i in 1:nk
        k_dict[copy(kind_arr[i,:])] = i
    end
    #    r_dict = Dict(d["tightbinding"]["r_dict"])

    if "scf" in keys(d["tightbinding"])
        tb_scf = parse(Bool,d["tightbinding"]["scf"])
    else
        tb_scf = false
    end
    if "h1" in keys(d["tightbinding"])
        h1 = parse_str_ARR_float(d["tightbinding"]["h1"])
        h1 = reshape(h1, nwan, nwan)
    else
        h1 = missing
    end
    if tb_scf == false
        h1 =  missing
    end

    
    function readstr(st)

        H = zeros(Complex{Float64}, nwan,nwan,nk, nspin)
        S = zeros(Complex{Float64}, nwan,nwan,nk)        
        
        lines = split(st, "\n")
        for line in lines
            sp = split(line)
            if length(sp) == 7 || length(sp) == 9
                m = parse(Int64,sp[2])
                n = parse(Int64,sp[3])
                r = parse(Int64,sp[1])
                if nspin == 1
                    H[m,n,r,1] = parse(Float64,sp[4]) + im*parse(Float64,sp[5])
                    S[m,n,r] = parse(Float64,sp[6]) + im*parse(Float64,sp[7])
                else
                    H[m,n,r,1] = parse(Float64,sp[4]) + im*parse(Float64,sp[5])
                    H[m,n,r,2] = parse(Float64,sp[6]) + im*parse(Float64,sp[7])
                    S[m,n,r] = parse(Float64,sp[8]) + im*parse(Float64,sp[9])
                end
            end
        end

        return H, S
    end
        
    Hk,Sk = readstr(d["tightbinding"]["Hk"])

    println("hk ", sum(abs.(Hk)))
    
#    println("read ")
#    println("h1")
#    println(h1)

    tb = make_tb_k(Hk, kind_arr, kweights, Sk, h1=h1, grid=grid, nonorth=nonorth)

    tbck = make_tb_crys_kspace(tb, crys, nelec, dftenergy, scf=scf, eden=eden, gamma=gamma)

    return tbck
        
#    catch
#        println("ERROR: Failed to read ", filename)
#        return -1
#    end
    
end


"""
    function write_tb_crys(filename, tbc::tb_crys)

Writes to `filename` a `tb_crys` object, using xml formatting. See `read_tb_crys`

"""
function write_tb_crys(filename, tbc::tb_crys)
    """
    write xml tb_crys object
    """

    doc = XMLDocument()
    root = ElementNode("root")
    setroot!(doc, root)

    crystal = ElementNode("crystal")
    link!(root, crystal)
    
    addelement!(crystal, "A", arr2str(tbc.crys.A))
    addelement!(crystal, "nat", string(tbc.crys.nat))
    addelement!(crystal, "coords", arr2str(tbc.crys.coords))
    addelement!(crystal, "types", str_w_spaces(tbc.crys.types))

    addelement!(root, "nelec", string(tbc.nelec))
    addelement!(root, "nspin", string(tbc.nspin))
    addelement!(root, "dftenergy", string(tbc.dftenergy))
    addelement!(root, "scf", string(tbc.scf))
    addelement!(root, "gamma", arr2str(tbc.gamma))
    addelement!(root, "eden", arr2str(tbc.eden))

    addelement!(root, "efermi", string(tbc.energy))
    addelement!(root, "energy", string(tbc.efermi))

    
    tightbinding = ElementNode("tightbinding")
    link!(root, tightbinding)
    
    addelement!(tightbinding, "nr", string(tbc.tb.nr))
    addelement!(tightbinding, "nwan", string(tbc.tb.nwan))
    addelement!(tightbinding, "nspin", string(tbc.tb.nspin))
    addelement!(tightbinding, "nonorth", string(tbc.tb.nonorth))
    addelement!(tightbinding, "ind_arr",arr2str(tbc.tb.ind_arr))
    addelement!(tightbinding, "r_dict",string(tbc.tb.r_dict))

    addelement!(tightbinding, "scf",string(tbc.tb.scf))
    addelement!(tightbinding, "h1",arr2str(tbc.tb.h1))
    
    function makestr(H,S, nonorth, ind_arr)

        nw = size(H)[2]
        nr = size(H)[4]
        nspin = size(H)[1]
        #        st = ""
        io_tmp = IOBuffer() #writing to iostream is much faster than making a giant string 
        #strings in julia are immutable, so the entire thing had to be copied over and over

        for r = 1:nr
            for m = 1:nw
                for n = 1:nw            
                    if nonorth 
                        Sr = real(S[m,n,r])
                        Si = imag(S[m,n,r])
                    else
                        if (m == n) && ind_arr[r,1] == 0 && ind_arr[r,2] == 0 && ind_arr[r,3] == 0 
                            Sr = 1.0
                        else
                            Sr = 0.0
                        end

                        Si = 0.0
                    end 

                    if nspin == 2
                        @printf(io_tmp, "% 3s % 3s % 3s % 2.10f % 2.10f % 2.10f % 2.10f % 2.10f % 2.10f \n", r, m,n,real(H[1, m,n,r]), imag(H[1, m,n,r]),real(H[2, m,n,r]), imag(H[2, m,n,r]), Sr, Si)

                    else
                        @printf(io_tmp, "% 3s % 3s % 3s % 2.10f % 2.10f % 2.10f % 2.10f \n", r, m,n,real(H[1, m,n,r]), imag(H[1, m,n,r]),Sr, Si)                   
                    end
                end
            end
        end
        st= String(take!(io_tmp))
        close(io_tmp)
        
        return st

    end
    println("a")
    st = makestr(tbc.tb.H, tbc.tb.S, tbc.tb.nonorth, tbc.tb.ind_arr)
    
    addelement!(tightbinding, "H", st)
    
#    println(doc)
#
#    println()

    #prettyprint
    io=open(filename, "w")
    prettyprint(io, doc);
    close(io)

    #  write(io, doc);
#    write(filename, doc);    

    return Nothing
    
#    return doc
    
end

"""
    function write_tb_crys_kspace(filename, tbc::tb_crys_kspace)

Save a `tb_crys_kspace` object to xml format. See `read_tb_crys_kspace`
"""
function write_tb_crys_kspace(filename, tbc::tb_crys_kspace)
    """
    write xml tb_crys kspace object
    """

    doc = XMLDocument()
    root = ElementNode("root")
    setroot!(doc, root)

    crystal = ElementNode("crystal")
    link!(root, crystal)
    
    addelement!(crystal, "A", arr2str(tbc.crys.A))
    addelement!(crystal, "nat", string(tbc.crys.nat))
    addelement!(crystal, "coords", arr2str(tbc.crys.coords))
    addelement!(crystal, "types", str_w_spaces(tbc.crys.types))

    addelement!(root, "nelec", string(tbc.nelec))
    addelement!(root, "nspin", string(tbc.nspin))
    addelement!(root, "dftenergy", string(tbc.dftenergy))
    addelement!(root, "scf", string(tbc.scf))
    addelement!(root, "gamma", arr2str(tbc.gamma))
    addelement!(root, "eden", arr2str(tbc.eden))

    tightbinding = ElementNode("tightbinding")
    link!(root, tightbinding)
    
    addelement!(tightbinding, "nk", string(tbc.tb.nk))
    addelement!(tightbinding, "nwan", string(tbc.tb.nwan))
    addelement!(tightbinding, "nonorth", string(tbc.tb.nonorth))
    addelement!(tightbinding, "kind_arr",arr2str(tbc.tb.K))
    addelement!(tightbinding, "kweights",arr2str(tbc.tb.kweights))
    addelement!(tightbinding, "k_dict",string(tbc.tb.k_dict))

    addelement!(tightbinding, "scf",string(tbc.tb.scf))
    addelement!(tightbinding, "h1",arr2str(tbc.tb.h1))

    addelement!(tightbinding, "grid",arr2str(tbc.tb.grid))
    
    function makestr(H,S, nonorth)

        nspin = size(H)[4]        
        nw = size(H)[1]
        nr = size(H)[3]

#        st = ""
        io_tmp = IOBuffer() #writing to iostream is much faster than making a giant string 
                            #strings in julia are immutable, so the entire thing had to be copied over and over

        for r = 1:nr
            for m = 1:nw
                for n = 1:nw            
                    if nonorth
                        Sr = real(S[m,n,r])
                        Si = imag(S[m,n,r])
                    else
                        if (m == n) 
                            Sr = 1.0
                        else
                            Sr = 0.0
                        end
                        Si = 0.0
                    end
                    
                    if  nspin == 1
                        @printf(io_tmp, "% 3s % 3s % 3s % 2.10f % 2.10f % 2.10f % 2.10f \n", r, m,n,real(H[m,n,r,1]), imag(H[m,n,r,1]),Sr, Si)
                    else
                        @printf(io_tmp, "% 3s % 3s % 3s % 2.10f % 2.10f % 2.10f % 2.10f  % 2.10f % 2.10f \n", r, m,n,real(H[m,n,r,1]), imag(H[m,n,r,1]),real(H[ m,n,r,2]), imag(H[m,n,r,2]), Sr, Si)
                    end
#                    st = st*t

                end
            end
        end
        st= String(take!(io_tmp))
        close(io_tmp)

        return st

    end

    st = makestr(tbc.tb.Hk, tbc.tb.Sk, tbc.tb.nonorth)
    
    addelement!(tightbinding, "Hk", st)
    
#    println(doc)
#
#    println()

    #prettyprint
    io=open(filename, "w")
    prettyprint(io, doc);
    close(io)

    #  write(io, doc);
#    write(filename, doc);    

    return Nothing
    
#    return doc
    
end
    
"""
    function make_tb_crys(ham::tb,crys::crystal, nelec::Float64, dftenergy::Float64; scf=false, eden = missing, gamma=missing, within_fit=true, screening=1.0, tb_energy=-999, fermi_energy=0.0 )

Constructor function for `tb_crys` object
"""
function make_tb_crys(ham::tb,crys::crystal, nelec::Float64, dftenergy::Float64; scf=false, eden = missing, gamma=missing, within_fit=true, screening=1.0, tb_energy=-999, fermi_energy=0.0 )

    T = typeof(crys.coords[1,1])
    nspin = ham.nspin
    if ismissing(eden)
        if scf == false
            eden = zeros(nspin,ham.nwan)
        else
            eden = get_neutral_eden(crys, ham.nwan)
            t = get_neutral_eden(crys, hamk.nwan)
            eden = zeros(nspin, hamk.nwan)
            if nspin == 1
                eden[1,:] = t
            elseif nspin == 2
                eden[1,:] = t
                eden[2,:] = t
            end
            
        end
    end

    if ismissing(gamma) 
#        println("gamma")
        gamma = electrostatics_getgamma(crys, screening=screening) #do this once and for all
    end

    nspin = ham.nspin
    
#    println("type scf " , typeof(scf))
#    println("type gamma " , typeof(gamma))
#    println("type eden " , typeof(eden))
    
#    return tb_crys{T}(ham,crys,nelec, dftenergy, scf, gamma, eden)
    return tb_crys{T}(ham,crys,nelec, dftenergy, scf, gamma, eden, within_fit, tb_energy, fermi_energy, nspin)
end

"""
    function make_tb_crys_kspace(hamk::tb_k,crys::crystal, nelec::Float64, dftenergy::Float64; scf=false, eden = missing, gamma=missing, screening=1.0)

Constructor function for `tb_crys_kspace` object
"""
function make_tb_crys_kspace(hamk::tb_k,crys::crystal, nelec::Float64, dftenergy::Float64; scf=false, eden = missing, gamma=missing, screening=1.0)

    nspin = hamk.nspin
    
    T = typeof(crys.coords[1,1])

    if ismissing(eden)
        if scf == false
            eden = zeros(nspin, hamk.nwan)
        else
            t = get_neutral_eden(crys, hamk.nwan)
            eden = zeros(nspin, hamk.nwan)
            if nspin == 1
                eden[1,:] = t
            elseif nspin == 2
                eden[1,:] = t
                eden[2,:] = t
            end

        end
    end
    
    if ismissing(gamma) 
        gamma = electrostatics_getgamma(crys, screening=screening) #do this once and for all
    end

    return tb_crys_kspace{T}(hamk,crys,nelec,nspin, dftenergy, scf, gamma, eden, -999.0)
end


"""
    function make_tb(H, ind_arr, r_dict::Dict; h1=missing)

Constructor function for `tb`
"""
function make_tb(H, ind_arr, r_dict::Dict; h1=missing)
    nw=size(H,2)
    if nw != size(H)[3]
        s=size(H)
        error("error in make_tb $s")
    end
    nr=size(H,4)
    nspin = size(H, 1)

    
    if size(ind_arr) != (nr,3)
        error("make_tb ind_arr ", size(ind_arr))
    end

    T=typeof(real(H[1,1,1,1]))
    
    S =zeros(Complex{T}, size(H))

    if ismissing(h1) 
        h1 =zeros(T, nw, nw)
        scf = false
    else
        scf = true
    end

    return tb{T}(H, ind_arr, r_dict,nw, nr, nspin, false, S, scf, h1)
end

"""
    function make_tb(H, ind_arr, r_dict::Dict, S; h1=missing)

Constructor function for `tb` with overlaps
"""
function make_tb(H, ind_arr, r_dict::Dict, S; h1=missing)
    nw=size(H,2)
    if nw != size(H,3) 
        exit("error in make_tb ", size(H))
    end
    nr=size(H,4)
    nspin=size(H,1)

    if size(ind_arr) != (nr,3) 
        error("make_tb ind_arr size ", size(ind_arr), size(H))
    end

    if size(S) != (nw,nw,nr)
        error("make_tb size S doesn't match size H", size(S), size(H))
    end
    
    T=typeof(real(H[1,1,1,1]))

    if ismissing(h1) 
        h1 =zeros(T, nw, nw)
        scf = false
    else
        scf = true
    end

#    println("size H ", size(H), " ", typeof(H), " npin $nspin ")
#    println("T ", T)
#    println("S ", size(S), " " , typeof(S))
#    println("size h1 ", size(h1), "  ", typeof(h1))
    return tb{T}(H, ind_arr,  r_dict,nw, nr, nspin, true, S, scf, h1)
end





"""
    function make_tb_k(Hk, K, kweights, Sk; h1=missing, grid=[0,0,0], nonorth=true)

Constructor for `tb_kspace`
"""
function make_tb_k(Hk, K, kweights, Sk; h1=missing, grid=[0,0,0], nonorth=true)
    nw=size(Hk,1)
    nk=size(Hk,3)
    nspin=size(Hk,4)

#    println("make_tb_k $nw $nk xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx") 
#    for k in 1:nk
#        println("K $k ", K[k,:])
#    end

    if size(K) != (nk,3) 
        error("make_tb_k kind_arr size ", size(K), size(Hk))
    end

    if size(Sk) != (nw,nw,nk)
        error("make_tb_k size S doesn't match size H", size(Sk), size(Hk))
    end
    
    T=typeof(real(Hk[1,1,1,1]))

    K2 = zeros(Float64,size(K))
    k_dict = Dict()  
    for i in 1:nk
        K2[i,:] = Int64.(round.(K[i,:] * 100000000))/100000000
        k_dict[K2[i,:]] = i
    end

    if ismissing(h1) 
        h1 =zeros(T, nw, nw)
        scf = false
    else
        scf = true
    end
#=
    println("make_tb_k")
    println(typeof(Hk))
    println(typeof(K2))
    println(typeof(kweights))
    println(typeof(k_dict))
    println(typeof(nw))
    println(typeof(nk))
    println(typeof(true))
    println(typeof(Sk))
    println(typeof(scf))
    println(typeof(h1))
    println(typeof(grid))
=#
    return tb_k{T}(Hk, K2, kweights, k_dict,nw, nk,nspin,  nonorth, Sk, scf, h1, grid)
end


function make_rdict(ind_arr)
    r_dict = Dict()
    for c = 1:size(ind_arr)[1]
        r_dict[ind_arr[c,:]] = c
    end
    return r_dict
end

        


"""
    function make_tb(H, ind_arr, S; h1=missing)

Constructor function for `tb`, better programming.
"""
function make_tb(H, ind_arr, S; h1=missing)

    r_dict = make_rdict(ind_arr)
    
    nw=size(H,2)
    if nw != size(H,3) 
        exit("error in make_tb ", size(H))
    end
    nr=size(H,4)
    nspin=size(H,1)

    if size(ind_arr) != (nr,3) 
        error("make_tb ind_arr size ", size(ind_arr), size(H))
    end

    if size(S) != (nw,nw,nr)
        error("make_tb size S doesn't match size H", size(S), size(H))
    end
    
    vartype=typeof(H[1,1,1, 1])

    if ismissing(h1)
        h1 =zeros(Float64, nw, nw)
        scf = false
    else
        scf = true
    end

    
    return tb{vartype}(H, ind_arr, r_dict,nw, nr, nspin, true, S, scf, h1)
end



"""
    function load_hr_dat(filename, directory="")
Load a wannier90 hr.dat file
Not currently a major part of program, but you can use if you want.
"""
function load_hr_dat(filename, directory="")
"""
Load hr.dat file from wannier90
"""

    convert_ev_ryd = 1.0/13.605693122
    if directory != ""
        f=gzopen("$directory/$filename", "r")
    else
        f=gzopen("$filename", "r")
    end
    
    hr_file = readlines(f)
    close(f)

    println(hr_file[1]) #header
    nwan = parse(Int,hr_file[2])          
    nr = parse(Int,hr_file[3])          

    println("nwan:", nwan)
    println("nr:", nr)    


    #parse weights
    if nr % 15 == 0
        wlines = convert(Int, floor(nr / 15))
    else
        wlines = convert(Int, floor(nr / 15)) + 1
    end
    println("wlines ", wlines)
    weights = zeros(nr)
    for i = 4:(3+wlines)
        line = hr_file[i]

        ind_start = (i-4)*15 + 1
        ind_end = (i-3)*15

#        println(line)
#        println(i, " ", ind_start," ", ind_end)

        if ind_end > nr
            ind_end = nr
        end
        weights[ind_start:ind_end] = map(x->parse(Float64,x),split(line))

    end

#    println(weights)

    #now parse main Ham

    H = zeros(Complex{Float64},1, nwan,nwan,nr)

    r_dict = Dict()
    ind_arr = zeros(nr,3)
    
    nind = 0

    nonorth=false
    sp = split(hr_file[4+wlines])
    if length(sp) == 9
        nonorth = true
        S = zeros(Complex{Float64}, nwan,nwan, nr)
        
    end


    
    for i = (4+wlines):length(hr_file)
        sp = split(hr_file[i])

        rind = [parse(Int, sp[1]), parse(Int, sp[2]),parse(Int, sp[3])]
        
        if rind in keys(r_dict)
            ind = r_dict[rind]
        else
            nind += 1
            r_dict[rind] = nind
            ind = nind
            ind_arr[ind,:] = rind
        end

        
        
        nw1 = parse(Int64,sp[4])
        nw2 = parse(Int64,sp[5])
        
        h = (parse(Float64, sp[6]) + parse(Float64, sp[7])*im) / weights[ind] * convert_ev_ryd
        H[1,nw1,nw2,ind] = h

        if nonorth
            s = (parse(Float64, sp[8]) + parse(Float64, sp[9])*im) / weights[ind]
            S[nw1,nw2,ind] = s

        end
        
    end

#    println(H[1,1,1])
#    println(H[2,2,2])

    if nonorth
        return make_tb(H, ind_arr, r_dict, S)
    else
        return make_tb(H, ind_arr, r_dict)
    end
end


"""
    function Hk(h::tb_crys_kspace, kpoint)

 Calculate band structure at a k-point from a `tb_crys_kspace` object.
 Note, can only return precalculated k-points.
 Need real-space version to get arbitrary k-points.

 #Returns
 - `vect` - Eigenvectors num_wan × num_wan complex matrix at kpoint
 - `vals` - Eigenvalues (num_wan)
 - `hk` - Hamiltonian at kpoint
 - `sk` - Overlap matrix at kpoint
 - `vals0` - <vect | Hk0 | vect> where Hk0 is the non-scf part of the Hamiltonian.

 #Arguments
 - `h::tb_crys_kspace` - tb_crys_kspace object
 - `kpoint` - e.g. `[0.0,0.0,0.0]`
 - `scf=missing` - default is to take SCF from h.
 """
 function Hk(h::tb_crys_kspace, kpoint; scf=missing, spin=1)

     return Hk(h.tb, kpoint, scf=scf, spin=spin)

 end

 """
     function Hk(h::tb_k, kpoint; scf=missing, spin=1)

 Calculate band structure at a k-point from `tb_k`. Must be pre-calculated k-point.
 """
 function Hk(h::tb_k, kpoint; scf=missing, spin=1) 

     if ismissing(scf)
         scf = h.scf
     end

     if spin == 2 && h.nspin == 1
         println("ERROR, asking for spin 2 from nsp ham")
     end

     
     kpoint = vec(kpoint)

     hk = zeros(Complex{Float64}, h.nwan, h.nwan)
     sk = zeros(Complex{Float64}, h.nwan, h.nwan)
     hk0 = zeros(Complex{Float64}, h.nwan, h.nwan)

     if kpoint in keys(h.k_dict)
         ind = h.k_dict[kpoint]
     else
         ind = 0
         for i = 1:h.nk
 #            println(kpoint')
 #            println(size(kpoint'))
 #            println(h.K[i,:])
 #            println(size(h.K[i,:]))
             if sum(abs.(kpoint - h.K[i,:])) < 5e-4
                 ind = i
                 break
             end
         end
         if ind == 0
             println("error, k not found!!! , $kpoint tb_k")
             println(tb)
             for i = 1:h.nk
                 println(h.K[i,:])
             end

             return
         end
     end

     hk0[:,:] = h.Hk[:,:,ind, spin]
     hk0 = 0.5*(hk0 + hk0')

     sk[:,:] = h.Sk[:,:,ind]
     sk = 0.5*(sk + sk')

     if scf
 #        println("add SCF")
         hk = hk0 + 0.5 * sk .* (h.h1 + h.h1')
     else
         hk[:,:] = hk0[:,:]
     end
     hk = 0.5*(hk[:,:] + hk[:,:]')

     nw=size(hk)[1]
     vects = zeros(nw,nw)
     vals = zeros(nw)
     vals0 = zeros(nw)

     try
         if h.nonorth
             sk = 0.5*(sk[:,:] + sk[:,:]')
             F=eigen(hk[:,:], sk[:,:])
         else
             F=eigen(Hermitian(hk)) #orthogonal
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
     return vects, vals, hk, sk, vals0

 end


 """
     function Hk(hk,sk, h::tb, kpoint; spin=1)

 Hk function with pre-allocated memory hk, sk
 """
 function Hk(hk,sk, h::tb, kpoint; spin=1)

     kpoint = vec(kpoint)

 #    println("Hk kpoint : $kpoint")

 #    println("Hk nonorth: " , h.nonorth)

 #    println("hk", typeof(hk), size(hk))

 #    println("zero")

     hk0 = zeros(Complex{Float64}, size(hk))

     fill!(hk, zero(Float64))

     if h.nonorth  
         fill!(sk, zero(Float64))
     end

 ##    Hk = zeros(Complex{Float64}, h.nwan, h.nwan)

     twopi_i = -1.0im*2.0*pi

 #    println("repeat")
     kmat = repeat(kpoint', h.nr,1)

 #    println("exp")
     exp_ikr = exp.(twopi_i * sum(kmat .* h.ind_arr,dims=2))

     for m in 1:h.nwan
         for n in 1:h.nwan        
             hk0[m,n] = h.H[spin,m,n,:]'*exp_ikr[:]
             if h.nonorth
                 sk[m,n] = h.S[m,n,:]'*exp_ikr[:]
             end
         end
     end
     hk0 = 0.5*(hk0 + hk0')

     if h.scf
 #        println("add SCF")
         hk = hk0 + sk .* h.h1
     else
         hk[:,:] = hk0[:,:]
     end
     hk = 0.5*(hk[:,:] + hk[:,:]')

 #    ex = 0.0 + im*0.0
 #    for n = 1:h.nr
  #       ex = exp(-twopi_i*(transpose(h.ind_arr[n,:])*kpoint))
  #       hk .+= ex .* h.H[:,:,n]
  #       if h.nonorth
  #           sk .+= ex .* h.S[:,:,n]
  #       end
  #   end
     #end




     nw=size(hk)[1]
     vects = zeros(nw,nw)
     vals = zeros(nw)
     vals0 = zeros(nw)

     try
         if h.nonorth
             sk = 0.5*(sk[:,:] + sk[:,:]')
             F=eigen(hk[:,:], sk[:,:])
         else
             #        println("orth")
             #        println(typeof(hk))
             hk = 0.5*(hk[:,:] + hk[:,:]')            
             F=eigen(hk[:,:]) #orthogonal
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

     #    F=eigen(hk, sk)

 #    println("hk")
 #    println(hk)

 #    println("hk vals")
 #    println(vals)

     return vects, sort!(vals), hk, sk, vals0

 end



 """
     function calc_bands(tbc::tb_crys, kpoints::Array{Float64,2}; spin=1)

 Calculate bandstructure for k-points from k-point array. Returns eigenvalues.

 # Arguments
 - `tbc::tb_crys` - The tight binding object
 - `kpoints::Array{Float64,2}` - k-point array. e.g. `[0.0 0.0 0.0; 0.0 0.0 0.1]`
 - `spin::Int ` - spin component (1 or 2) if magnetic, only 1 is non-sp.
 """
 function calc_bands(tbc::tb_crys, kpoints::Array{Float64,2})
 #    if tbc.scf == true
 #        h1, dq = get_h1(tbc, tbc.eden)
 #    else
 #        h1 == missing
 #    end
     return calc_bands(tbc.tb, kpoints)
 end

 function calc_bands(tbc::tb_crys_kspace, kpoints::Array{Float64,2})

     return calc_bands(tbc.tb, kpoints)

 end

 """
     function calc_bands(h, kpoints::Array{Float64,2}; spin=1)

 Calculate bandstructure for k-points from k-point array. h is a `tb` or `tb_k` object.
 """
 function calc_bands(h, kpoints::Array{Float64,2})
 """
 calculate band structure at group of points
 """
 #    println("calc_bands")
 #    println("h.scf ", h.scf)
 #    println("temp arrays")

 #    T=Float64
 #    
 #    hktemp= zeros(Complex{T}, h.nwan, h.nwan)
 #    sktemp= zeros(Complex{T}, h.nwan, h.nwan)


     nk = size(kpoints,1)

 #    println("vals time")

     Vals = zeros(Float64,  nk,h.nwan,h.nspin )

     for spin = 1:h.nspin
         for i = 1:nk
             
             vect, vals, hk, sk, vals0 = Hk(h, kpoints[i,:], spin=spin)
             Vals[i,:, spin] = vals
             #        if sum(abs.(kpoints[i,:])) == 0
             #            println("k0")
             #            println(vals)
             #        end
             
         end
     end
     return Vals

 end    


 """
     function summarize_orb(orb::Symbol)

 Convert exact orbital (:px) to type of orbital :p .
 """
 function summarize_orb(orb::Symbol)
     if orb == :s
         return :s
     elseif orb == :p || orb == :px || orb == :py || orb == :pz
         return :p
     elseif orb == :d || orb == :dz2 || orb == :dxz || orb == :dyz || orb == :dx2_y2 || orb == :dxy
         return :d
     else
         println("not coded yet $orb summarize_orb")
     end
 end

 function summarize_orb_num(orb::Symbol)
     if orb == :s
         return 0
     elseif orb == :p || orb == :px || orb == :py || orb == :pz
         return 1
     elseif orb == :d || orb == :dz2 || orb == :dxz || orb == :dyz || orb == :dx2_y2 || orb == :dxy
         return 2
     else
         println("not coded yet $orb summarize_orb")
     end
 end

 function orb_num(orb::Symbol)
     if orb == :s
         return 1
     elseif orb == :pz
         return 2
     elseif orb == :px
         return 3
     elseif orb == :py
         return 4
     elseif orb == :dz2
         return 5
     elseif orb == :dxz
         return 6
     elseif orb == :dyz
         return 7
     elseif orb == :dx2_y2
         return 8
     elseif orb == :dxy
         return 9
     elseif orb == :fz3
         return 10
     elseif orb == :fxz2
         return 11
     elseif orb == :fyz2
         return 12
     elseif orb == :fz_x2y2
         return 13
     elseif orb == :fxyz
         return 14
     elseif orb == :fx_x2_3y2
         return 15
     elseif orb == :fy_3x2_y2
         return 16
     else
         println("not coded yet $orb orb_num")
     end
 end



 #function calc_bands(h::tb, kpoints::Array{Float64,2})
 #"""
 #calculate band structure at group of points
 #"""
 #
 ##    println("temp arrays")
 #    hktemp= zeros(Complex{Float64}, h.nwan, h.nwan)
 #    sktemp= zeros(Complex{Float64}, h.nwan, h.nwan)
 #
 #    nk = size(kpoints,1)
 #
 ##    println("vals time")
 #    
 #    Vals = zeros(Float64, nk,h.nwan)
 #    
 #    for i = 1:nk
 #
 #        vect, vals, hk, sk = Hk(hktemp, sktemp, h, kpoints[i,:])
 #        Vals[i,:] = vals
 #
 #    end
 #
 #    return Vals
 #
 #end    
 #


 """
     function Hk(h::tb, kpoint; spin=1)

 Calculate band structure at a k-point from `tb`
 """
 function Hk(h::tb, kpoint; spin=1)

     T=typeof(real(h.H[1,1,1,1]))
     #if we don't have a temp array
     hktemp= zeros(Complex{T}, h.nwan, h.nwan)
     sktemp= zeros(Complex{T}, h.nwan, h.nwan)    
     return Hk(hktemp, sktemp, h, kpoint, spin=spin)
 end

 """
     function Hk(h::tb_crys, kpoint; spin=1)

 Calculate band structure at a k-point from `tb_crys`
 """
 function Hk(tbc::tb_crys, kpoint; spin=1 )


     return Hk(tbc.tb, kpoint, spin=spin)

 end


 """
     function calc_energy_fft(tbc::tb_crys; grid=missing, smearing=0.01, return_more_info=false)

 Get energy using fft.

 returns energy

 if `return_more_info==true` then returns
 `etot, efermi, vals, vects`

 """
 function calc_energy_fft(tbc::tb_crys; grid=missing, smearing=0.01, return_more_info=false)

     etypes = types_energy(tbc.crys)


     hk3, sk3 = myfft_R_to_K(tbc, grid)

     ret =  calc_energy_fft_band(hk3, sk3, tbc.nelec, smearing=smearing, return_more_info=return_more_info, h1=tbc.tb.h1)


     if return_more_info
         energy, efermi, vals, vects  = ret

 #        println("PRE-precheck ", sum(vects[1,:,:]' * sk3[:,:,1,1,1] * vects[1,:,:]))

         etot = etypes + energy
         return etot, efermi, vals, vects, hk3, sk3

     end
     return ret + etypes

 end

 """
     function calc_energy_fft_band(hk3, sk3, nelec; smearing=0.01, return_more_info=false, h1 = missing)

 Return energy from hamiltonian `hk3`, overlap `sk3`, `nelec`, etc.
 Primarly for internal calling after fft.
 """
 function calc_energy_fft_band(hk3, sk3, nelec; smearing=0.01, return_more_info=false, h1 = missing)
 #h1 is the scf contribution

     grid = size(sk3)[3:5]
     nk = prod(grid)
     nwan = size(hk3)[1]
     nspin = size(hk3)[3]

     VALS = zeros(Float64, nk,nwan,nspin)

     if return_more_info == true
         VECTS = zeros(Complex{Float64}, nk,nspin, nwan, nwan)
     end

     c=0

     thetype=typeof(real(sk3[1,1,1,1,1]))
     sk = zeros(Complex{thetype}, nwan, nwan)
     hk = zeros(Complex{thetype}, nwan, nwan)


     for k1 = 1:grid[1]
         for k2 = 1:grid[2]
             for k3 = 1:grid[3]
                 c += 1
                 try

                     sk[:,:] = 0.5*(sk3[:,:,k1,k2,k3] + sk3[:,:,k1,k2,k3]')
                     for spin = 1:nspin
                         hk[:,:] = 0.5*(hk3[:,:,spin, k1,k2,k3] + hk3[:,:,spin, k1,k2,k3]')
                         
                         if !ismissing(h1)
                             hk += sk .* h1
                         end
                         
                         vals, vects = eigen(hk, sk)
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
                     sk[:,:] = 0.5*(sk3[:,:,k1,k2,k3] + sk3[:,:,k1,k2,k3]')
                     valsS, vectsS = eigen(sk)
                     println(valsS)

                     rethrow(error("BadOverlap"))

                 end
             end
         end
     end

     band_en, efermi = band_energy(VALS, ones(nk), nelec, smearing, returnef=true)
     energy_smear = smearing_energy(VALS, ones(nk), efermi, smearing)

     if return_more_info
         return  band_en + energy_smear, efermi, VALS, VECTS
     else
         return  band_en + energy_smear
     end

 end 


 """
     function calc_energy_charge_fft_band(hk3, sk3, nelec; smearing=0.01, h1 = missing)

 Calculate energy and charge density. For internal use.

 `return energy0, efermi, chargeden[:], VECTS, VALS, error_flag`
 """
 function calc_energy_charge_fft_band(hk3, sk3, nelec; smearing=0.01, h1 = missing)

#     println("memory")
     if true

         
         grid = size(sk3)[3:5]
         #    print("calc_energy_charge_fft_band grid $grid")
         nk = prod(grid)
         nwan = size(hk3)[1]
         nspin = size(hk3)[3]
         VALS = zeros(Float64, nk, nwan, nspin)
         VALS0 = zeros(Float64, nk,nwan, nspin)
#         c=0

         thetype=typeof(real(sk3[1,1,1,1,1]))
#         sk = zeros(Complex{thetype}, nwan, nwan)
#         hk = zeros(Complex{thetype}, nwan, nwan)
#         hk0 = zeros(Complex{thetype}, nwan, nwan)

         VECTS = zeros(Complex{thetype}, nk, nspin, nwan, nwan)
         cVECTS = zeros(Complex{thetype}, nk, nspin,  nwan, nwan)
         SK = zeros(Complex{thetype}, nk, nwan, nwan)

         error_flag = false

         if ismissing(h1)
             h1 = zeros(nwan,nwan)
         end

         
     end
#
#     println("grid", grid)

     @threads for c = 1:grid[1]*grid[2]*grid[3]

         k3 = mod(c-1 , grid[3])+1
         k2 = 1 + mod((c-1) ÷ grid[3], grid[2])
         k1 = 1 + (c-1) ÷ (grid[2]*grid[3])

#         k3 = mod(kind-1 , grid[1]*grid[2])+1
#         k1 = 1 + (kind-1) ÷ grid[2]*grid[[3]
#                                          
#     @time for k1 = 1:grid[1]
#         @threads for k2 = 1:grid[2]
#             for k3 = 1:grid[3]
#                 c = k3 + (k2-1) * grid[3] + (k1-1) * grid[2]*grid[3] 
#                 c += 1
#                 println("c $c cx $cx")
                 try
                     sk = 0.5*( (@view sk3[:,:,k1,k2,k3]) + (@view sk3[:,:,k1,k2,k3])')
                     SK[c,:,:] = sk
                     for spin = 1:nspin
                         hk0 = 0.5*( (@view hk3[:,:,spin, k1,k2,k3]) + (@view hk3[:,:,spin, k1,k2,k3])')
                         hk = hk0 + sk .* h1
                         vals, vects = eigen(hk, sk)
                         
                         VALS[c,:, spin] = vals
                         VALS0[c,:, spin] = real.(diag(vects'*hk0*vects))
                         VECTS[c,spin, :,:] = vects
                     end
                 catch e
                     if e isa InterruptException
                         println("user interrupt")
                         rethrow(InterruptException)
                     end

                     if error_flag == false
                         println("error calc_energy_fft $k1 $k2 $k3 usually negative overlap eig")
                         sk = 0.5*(sk3[:,:,k1,k2,k3] + sk3[:,:,k1,k2,k3]')
                         valsS, vectsS = eigen(sk)
                         println(valsS)
                         error_flag=true
                         rethrow(error("BadOverlap"))
                     end
                     error_flag=true

#                 end
#             end
         end
     end

     #println("stuff")
     begin
         maxSK = maximum(abs.(SK), dims=1)[1,:,:]

         
         energy, efermi = band_energy(VALS, ones(nk), nelec, smearing, returnef=true)
         
         occ = gaussian.(VALS.-efermi, smearing)

         max_occ = findlast(sum(occ, dims=[1,3]) .> 1e-10)[2]
         
#         println("occ ", size(occ))
#         println(occ)
         
         energy_smear = smearing_energy(VALS, ones(nk), efermi, smearing)
         #    println("energy smear $energy_smear")
         
         energy0 = sum(occ .* VALS0) / nk * 2.0
         
         #    println("ENERGY 0 : $energy0")
         
         energy0 += energy_smear
         
         #    println("sum occ ", sum(occ), "  ", sum(occ) / (grid[1]*grid[2]*grid[3]))#


         denmat = zeros(Float64, nspin, nwan, nwan)

     end
         
     if true
         TEMP = zeros(Complex{Float64}, nwan, nwan) 

         for spin = 1:nspin
             TEMP[:,:] .= 0.0
             pVECTS = permutedims(VECTS[:,spin,:,:][:,:,:], [1,3,2])
             pVECTS_C = conj(pVECTS)
             
             #        cVECTS = conj(VECTS)
             grid3 = prod(grid)
             t_temp = zeros(Complex{Float64},grid3 )
             @threads for i = 1:nwan
                 for j = 1:nwan
                     
                     if maxSK[i,j] > 1e-7
                         
                         TEMP[i,j] += sum( sum(  (@view occ[:,1:max_occ,spin]).* (@view pVECTS_C[:,1:max_occ,i])  .* (@view pVECTS[:,1:max_occ,j]), dims=2) .* (@view SK[:,i,j]))
                         
                     end
                     
                 end
             end
             TEMP =   (TEMP + conj(TEMP))
             denmat[spin,:,:] += 0.5* real.( TEMP)
             denmat[spin,:,:] = denmat[spin,:,:] / (grid[1]*grid[2]*grid[3])
         end
         chargeden = zeros(nspin, nwan)
         for spin = 1:nspin
             chargeden[spin,:] = sum(denmat[spin,:,:][:,:], dims=1)
         end
         
     end

     if nspin == 2
         energy0 = energy0 / 2.0
     end
         
     
     return energy0, efermi, chargeden, VECTS, VALS, error_flag


 end 

 """
     function calc_energy_charge_fft(tbc::tb_crys; grid=missing, smearing=0.01)

 Do fft, then calculate energy and charge.
 """
 function calc_energy_charge_fft(tbc::tb_crys; grid=missing, smearing=0.01)

     etypes = types_energy(tbc.crys)

     if ismissing(grid)
         grid = get_grid(tbc.crys)
     else
         if length(grid) != 3
             grid = get_grid(tbc.crys, grid) 
         end
     end
 #    println("calc_energy_charge_fft grid $grid")
     hk3, sk3 = myfft_R_to_K(tbc, grid)

     if tbc.scf
         h1 = tbc.tb.h1
         echarge, pot = ewald_energy(tbc)
     else
         h1 = missing
         echarge = 0.0
     end
     eband, efermi, chargeden, VECTS, VALS, error_flag  =  calc_energy_charge_fft_band(hk3, sk3, tbc.nelec, smearing=smearing, h1 = h1)
     tbc.efermi = efermi
     tbc.eden = chargeden
     #println("energy comps $eband $etypes $echarge")
     energy = eband + etypes + echarge

     return energy, efermi, chargeden, VECTS, VALS, error_flag

 end

"""
    function dumb_cd(tbc; grid=missing, smearing=0.01)

I don't remember what this function is for. It seems to have something to do with charge density, but it doesn't pay attention to the occupations
"""
 function dumb_cd(tbc; grid=missing, smearing=0.01)

     if ismissing(grid)
         grid = get_grid(tbc.crys)
     end

     hk3, sk3 = myfft_R_to_K(tbc, grid)

     h1 = tbc.tb.h1

     eband, efermi, chargeden, VECTS, VALS, error_flag  =  calc_energy_charge_fft_band(hk3, sk3, tbc.nelec, smearing=smearing, h1 = h1)


     K1 = size(sk3)[3]
     K2 = size(sk3)[3]
     K3 = size(sk3)[3]

     nw = size(sk3)[1]

     nspin = tbc.nspin
     
     charge = zeros(nspin, nw, nw)
     chargeX = zeros(nspin, nw)

     for spin = 1:nspin

         c=0
         for k1 in 1:K1
             for k2 in 1:K2
                 for k3 in 1:K3
                     c += 1
                     sk = 0.5*(sk3[:,:,k1,k2,k3] + sk3[:,:,k1,k2,k3]')
                     
                     for a = 1:nw
                         for ind1 = 1:nw
                             for ind2 = 1:nw
                                 charge[spin, ind1,ind2] += real( 0.5* ( VECTS[c,spin, ind1,a]' * VECTS[c,spin,ind2,a] * sk[ind1,ind2]  + VECTS[c,spin,ind1,a] * VECTS[c,spin,ind2,a]' * conj(sk[ind1,ind2])))
                             end
                         end
                     end

                     #                sk = 0.5*(sk3[:,:,k1,k2,k3] + sk3[:,:,k1,k2,k3]')
                     #                st = (sk ^ -0.5)
                     #                ht =  st *   (sk.*h1 + hk3[:,:,k1,k2,k3]) * st
                     #                val, vect = eigen(0.5*(ht + ht'))
                     #                println(val - VALS[c,:])
                     #                for a = 1:nband
                     #                    charge2 += 0.5*real(conj(vect[:,a]) .* vect[:,a] + (vect[:,a]) .* conj(vect[:,a]) )
                     #                end
                     

                 end
             end
         end
         charge[spin,:,:] = charge[spin,:,:] / K1/K2/K3
         chargeX[spin,:] =  sum(charge[spin,:,:][:,:], dims=1)[:]
     end
     return chargeX, charge
     
 end

 """
     function calc_energy(tbc::tb_crys; smearing=0.01, returnk=false)

 Calculate energy without fft.
 """
 function calc_energy(tbc::tb_crys; smearing=0.01, returnk=false)
     """
     calculate the energy from a kgrid
     """
     kgrid = get_grid(tbc.crys)

     return calc_energy(tbc, kgrid, smearing=smearing,returnk=returnk)

 end 

 """
     function types_energy(tbc::tb_crys)

 Calculate the reference energy of each atom in crystal. This results in the 
 total energy being indexed to seperated non-spin-polarized atoms. This is arbirary.
 """
 function types_energy(tbc::tb_crys)

     return types_energy(tbc.crys)

 end

 """
     function types_energy(tbc::tb_crys)
 """
 function types_energy(tbc::tb_crys_kspace)

     return types_energy(tbc.crys)

 end

 """
     function types_energy(c::crystal)
 """
 function types_energy(c::crystal)

     return types_energy(c.types)
 end

 """
     function types_energy(types)
 """
 function types_energy(types)

     et = 0.0
     for t in types
         et += atoms[t].energy_offset
     end
     return et

 end

 """
     function calc_energy(h::tb_crys, kgrid; smearing=0.01, returnk=false)

 Calculate energy no fft
 """
 function calc_energy(h::tb_crys, kgrid; smearing=0.01, returnk=false)
     """
     calculate the energy from a kgrid
     """
     if ismissing(kgrid)
         kgrid = get_grid(tbc.crys)
     end

     etypes = types_energy(h.crys)
     if h.scf
         energy_charge, pot = ewald_energy(h)
     else
         energy_charge = 0.0
     end

     ret = calc_energy_band(h.tb, h.nelec, kgrid, smearing=smearing,returnk=returnk)
 #    if returnk
 #        eband = ret[1]
 #        ek = ret[2]
 #        return eband+etypes, ek
 #    end

     return ret + etypes + energy_charge

 end 

 """
     function calc_energy_band(h::tb, nelec, kgrid; smearing=0.01, returnk=false)

 calculate energy no fft
 """
 function calc_energy_band(h::tb, nelec, kgrid; smearing=0.01, returnk=false)

     kpts, kweights = make_kgrid(kgrid)

 #    println(size(kpts), " ksize " , sum(kweights))
     eigs = calc_bands(h, kpts)
     en, efermi =  band_energy(eigs, kweights, nelec, smearing, returnef=true)
     e_smear  = smearing_energy(eigs, kweights, efermi)

     return en+e_smear

 end 

 """
     function make_kgrid(kgrid)

 -`kgrid` is an array of 3 integers like `[8,8,8]` 

 returns regular MP Gamma-centered k-point grid and (equal) k-weights.
 """
function make_kgrid(kgrid)
    kpts = zeros(Float64, prod(kgrid),3)
    kweights = ones(Float64, prod(kgrid))/prod(kgrid) * 2.0
    c=0
    for k1 = 0:kgrid[1]-1
        for k2 = 0:kgrid[2]-1
            for k3 = 0:kgrid[3]-1
                c+=1
                kpts[c,1] = Float64(k1)/Float64(kgrid[1])
                kpts[c,2] = Float64(k2)/Float64(kgrid[2])
                kpts[c,3] = Float64(k3)/Float64(kgrid[3])
            end
        end
    end
    return kpts, kweights
    
end

function make_kgrid(tbc::tb_crys)

    return make_krid(tbc.crys)
end

function  make_kgrid(c::crystal)

    return make_kgrid( get_grid(c))
    
end


 #=
 function trim_shift(tbc::tb_crys, tol=0.0002)

     energy_orig = tbc.dftenergy
     trim(tbc.tb, tol)
     energy_new = calc_energy(tbc)
     shift = (energy_orig - energy_new)/tbc.nelec
     c_zero = tbc.tb.r_dict[[0,0,0]]
     for i = 1:tbc.tb.nwan
         tbc.tb.H[i,i,c_zero] += shift
     end 
     energy_final = calc_energy(tbc)
     println("trim shift energy_orig $energy_orig energy_trimmed $energy_new energy_final $energy_final")
 end
 =#


 """
     function trim(h::tb, tol=0.0002)

 Remove terms in `tb` with abs value smaller than `tol`. (Ryd)
 Will speed calculations at cost of accuracy. USE WITH CARE.
 Usually not necessary except for very large or very detailed calculations.
 """
 function trim(h::tb, tol=0.0002)


     c=0
     keep = []
     npin = h.nspin
     for spin = 1:nspin
         for r = 1:h.nr
             if sum(abs.(h.H[spin,:,:,r]) .> tol) > 0 || h.ind_arr[r,:] == [0, 0, 0] 
                 if  ! (r in keep)
                     push!(keep, r)
                     c += 1
                 end
             end
         end
     end
     println("trim, out of ", h.nr, " keep $c")

     h.H = h.H[:,:,:,keep]
     if h.nonorth
         h.S = h.S[:,:,keep]
     end
     h.nr = c
     h.ind_arr = h.ind_arr[keep,:]

     h.r_dict = Dict()
     for i in 1:h.nr
         h.r_dict[copy(h.ind_arr[i,:])] = i
     end




 end    

 """
     function renormalize_tb(d::dftout, h::tb) 

 Shift eigenvalues from DFT calculation so that the band energy matches
 the DFT total energy. This doesn't do anything useful right now.
 """
 function renormalize_tb(d::dftout, h::tb)
 """
 Changes the tight binding matrix elements so that the total band energy 
 equals the total atomization energy . 


 """

     ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(d.crys)

     #eigs = calc_bands(h, d.bandstruct.kpts)

     band_en = band_energy(eigs, d.bandstruct.kweights, nval)

     println("original band energy ", band_en)

     etot_dft = d.energy

     atomization_energy = etot_dft - etotal_atoms

     println("atomization_enegy: ", atomization_energy)

     shift = (atomization_energy - band_en)/nval

     println("shift $shift nval $nval")

     ind = h.r_dict[[0,0,0]]

     println("ind $ind ", h.ind_arr[ind,:])

     for spin = 1:h.nspin
         for i = 1:h.nwan
             h.H[spin, i,i,ind] = h.H[spin, i,i,ind] + shift
         end
     end
     
     eigs = calc_bands(h, d.bandstruct.kpts)
     band_en_new = band_energy(eigs, d.bandstruct.kweights, nval)

     println("new band_energy ", band_en_new)

 end


 """
    function organizedata(tbc::tb_crys)

 Rearrange the data in the tbc as a function of distances for plotting purposes.

 `    return data_onsite, data_arr`

 Returns two arrays. The first has data on the onsite elements.

 - column 1 and 3 have atom indexes
 - column 2 and 4 have orbital numbers
 - columns 5 and 6 have real and imaginary parts of H
 - columns 11 and 12 have real and imaginary parts of S
 - column 7 has the closest inter-atomic distance
 - column 8 has the index of the closest atom.

 The second has intersite data

 - column 1 and 3 have atom indexes for the atom pairs
 - column 2 and 4 have orbital numbers
 - columns 5 and 6 have real and imaginary parts of H
 - columns 11 and 12 have real and imaginary parts of S
 - column 7 has the inter-atomic distance
 - column 8,9,10 have the the direction cosines lmn for the atom pair.

 """
 function organizedata(tbc::tb_crys; spin=1)

     return organizedata(tbc.crys, tbc.tb,spin=spin)

 end

 function plot_organizedata_offsite(data_arr, at1,o1,at2,o2; data = :H, color="blue", marker=:x)

     rows = size(data_arr)[1]
     inds = []
     for i = 1:rows
         (AT1,O1,AT2,O2) = data_arr[i,1:4]
         if at1 == AT1 && o1 == O1 && at2 == AT2 && o2 == O2 
             push!(inds, i)
         end
     end

     if data == :H
         println("scatter")
         display(scatter!(data_arr[inds,7], data_arr[inds,5], color=color, marker=marker))
     end

     return inds
     
 end

 """
     function organizedata(crys::crystal, h::tb)
 """
 function organizedata(crys::crystal, h::tb; spin=1)

     ind2orb, orb2ind, etotal, nval = orbital_index(crys)

     data = []

     data_arr = zeros(h.nwan*h.nwan*h.nr,12)

     At = transpose(crys.A)

     counter = 0
     ind_dir = Dict()
     for n = 1:h.nwan
         (na, nta, norb) = ind2orb[n]
         for m = 1:h.nwan        
             (ma, mta, morb) = ind2orb[m]            
             if na == ma
                 counter += 1
                 ind_dir[(n,m)] = counter
             end
         end
     end

     data_onsite = zeros(counter,12)
     data_onsite[:,7] .= 10000000000.0

     lmn = zeros(3)
     c=0
     c2 = 0

     c_zero = 0

     for n = 1:h.nwan
         (na, nta, norb) = ind2orb[n]
         for m = 1:h.nwan        
             (ma, mta, morb) = ind2orb[m]            
             for i in 1:h.nr

                 c+=1

                 #                println("n m i $n $m $i")

                 R = h.ind_arr[i,:]

                 if R[1] == 0 && R[2] == 0 && R[3] == 0
                     c_zero = i
                 end

 #                dR = At*(-crys.coords[na,:] .+ crys.coords[ma,:] .+ R) 
                 dR = At*(crys.coords[na,:] .- crys.coords[ma,:] .+ R) 

                 dist = sum(dR.^2)^0.5

                 if dist > 1e-7
                     c2+=1
                     if dist > 1e-7
                         lmn[:] = dR/dist
                     else
                         lmn .= 0.0
                     end
                     push!(data, [na, norb, ma, morb, h.H[spin, n,m,i], dist, lmn])

                     data_arr[c2, 1] = na
                     data_arr[c2, 2] = symbol_dict[norb]              
                     data_arr[c2, 3] = ma
                     data_arr[c2, 4] = symbol_dict[morb]
                     data_arr[c2, 5] = real(h.H[spin,  n,m, i])
                     data_arr[c2, 6] = imag(h.H[spin, n,m, i])
                     data_arr[c2, 7] = dist
                     data_arr[c2, 8] = lmn[1]
                     data_arr[c2, 9] = lmn[2]
                     data_arr[c2, 10] = lmn[3]
                     data_arr[c2, 11] = real(h.S[ n,m, i])
                     data_arr[c2, 12] = imag(h.S[ n,m, i])
                 end
             end
         end
     end
     for (n,m) in keys(ind_dir)
         c2 = ind_dir[(n,m)]
         (na, nta, norb) = ind2orb[n]
         (ma, mta, morb) = ind2orb[m]            

         data_onsite[c2,1] = na
         data_onsite[c2,3] = ma
         data_onsite[c2,2] = symbol_dict[norb]
         data_onsite[c2,4] = symbol_dict[morb]

         data_onsite[c2, 5] = real(h.H[spin, n,m, c_zero])
         data_onsite[c2, 6] = imag(h.H[spin, n,m, c_zero])

         data_onsite[c2, 11] = real(h.S[ n,m, c_zero])
         data_onsite[c2, 12] = imag(h.S[ n,m, c_zero])

         for x in 1:crys.nat
             for rx = [-1,0,1]
                 for ry = [-1,0,1]
                     for rz = [-1,0,1]
                         R = [rx,ry,rz]
                         dR = At*(-crys.coords[na,:] .+ crys.coords[x,:] .+ R) 
                         dist = sum(dR.^2)^0.5
                         if dist > 1e-7 && dist <= data_onsite[c2,7]
                             data_onsite[c2,7] = dist
                             data_onsite[c2,8] = x
                         end
                     end
                 end
             end
         end
     end
     return data_onsite, data_arr

 end


 """
     function myfft_R_to_K(tbc, grid=missing) 

 Does Fourier Transform R->k (fft) using FFTW for tb_crys
 """
 function myfft_R_to_K(tbc, grid=missing)


 #    wan, semicore, nwan, nsemi, wan_atom, atom_wan = tb_indexes(tbc.crys)

     if ismissing(grid)
         grid = get_grid(tbc.crys)
     else
         if length(grid) != 3
             grid = get_grid(tbc.crys, grid)
         end
     end

 #    println("grid ", grid)

 #    kgrid = make_kgrid(grid)

     nwan = tbc.tb.nwan
     nspin = tbc.nspin

     nr = prod(grid)

 #    if tbc.tb.nonorth

     ham_R = zeros(Complex{Float64},  nwan, nwan, nspin, grid[1], grid[2], grid[3])
     S_R = zeros(Complex{Float64}, nwan, nwan,  grid[1], grid[2], grid[3])

     ind = zeros(Int64, 3)
     new_ind = zeros(Int64, 3)

     for c in 1:size(tbc.tb.ind_arr)[1]

         ind[:] = tbc.tb.ind_arr[c,:]

         new_ind[1] = mod(ind[1], grid[1])+1
         new_ind[2] = mod(ind[2], grid[2])+1
         new_ind[3] = mod(ind[3], grid[3])+1

 #        println(ind, " new ", new_ind)

         for spin = 1:nspin
             ham_R[:,:,spin, new_ind[1], new_ind[2], new_ind[3]] += tbc.tb.H[spin, :,:,c]
         end
         S_R[:,:,new_ind[1], new_ind[2], new_ind[3]] += tbc.tb.S[:,:,c]

     end

     hamK3 = fft(ham_R, [4,5,6])
     SK3 = fft(S_R, [3,4,5])

     return hamK3, SK3

 end    

 """
     function myfft(crys, nonorth, grid, kpts,ham_kS, Sk=missing)

 Does Fourier Transform K->R (ifft) using FFTW.

 Arguments
 - `crys` crystal
 - `nonorth` nonorogonal bool
 - `grid` k-point grid size
 - `kpts` the k-points `nkpts`×`3` in the original order, to be rearranged into grid
 - `ham_kS` hamiltonian in k space (`nw`×`nw`×`nkpts`)
 - `Sk` overlaps in k space

 """
 function myfft(crys, nonorth, grid, kpts,ham_kS, Sk=missing)

     wan, semicore, nwan, nsemi, wan_atom, atom_wan = tb_indexes(crys)

     nspin = size(ham_kS)[4]
     nwan = size(ham_kS)[1]
     nks  = size(ham_kS)[3]

     println("myfft $nspin $nwan $nks")
     
     if length(size(ham_kS)) == 5 #we don't have to calculate, already in fft form
         println("don't calc")
         ham_k_fftw = ham_kS
         S_k_fftw = Sk

     else

         ham_k_fftw = zeros(Complex{Float64}, nspin, nwan, nwan, grid[1], grid[2], grid[3])
         S_k_fftw = zeros(Complex{Float64}, nwan, nwan, grid[1], grid[2], grid[3])

         for k = 1:nks
             #        println("k $k ", p.bs.kpts[k,:], " " , grid)
             kn = Int.(round.(kpts[k,:] .* grid)).+1
             for i in 1:3
                 if kn[i] <= 0
                     kn[i] = kn[i] + grid[i]
                 end
             end
             if nonorth
                 S_k_fftw[:,:,kn[1], kn[2], kn[3]] = Sk[:,:,k]            
             end
             for spin = 1:nspin
                 if nonorth
                     ham_k_fftw[spin, :,:,kn[1], kn[2], kn[3]] = ham_kS[:,:,k,spin]
                 else
                     ham_k_fftw[spin, :,:,kn[1], kn[2], kn[3]] = ham_kS[:,:,k,spin]
                 end
             end
         end

         println("ham_k_fftw ", size(ham_k_fftw))
     end    


     #does the actual ifft
     hamR3 = ifft(ham_k_fftw, [4,5,6])
     if nonorth
         SR3 = ifft(S_k_fftw, [3,4,5])    
     end


 #    println("myfft SR3 0 0 0 ")
 #    println(SR3[:,:,1,1,1])
 #    println()

 #    twopi_i = 1.0im*2.0*pi


     r_dict = Dict()

     R_grid, R_int_grid, sym_R = get_sym_R(crys, grid)
     ngrid2 = size(R_grid)[1]

 #    print("R_int_grid")
 #    print(R_int_grid[1:20, :])
 #    println()

     ind_arr = R_grid
     ham_r = zeros(Complex{Float64},  nspin, nwan, nwan,ngrid2)

     if nonorth
         S_r = zeros(Complex{Float64},  nwan, nwan,ngrid2)
     end

     sym = zero(Float64)
     exp_val = zero(Complex{Float64})
     c=0

     sym_mat = zeros(Float64, ngrid2)

 #    println("R_int_grid")

     for c = 1:ngrid2
         r_dict[ind_arr[c,:]] = c
         rint = R_int_grid[c,:]
         #println("$c ", R_int_grid[c,:] )
         #put ifft'd arrays into new arrays
         for nw_1 = 1:nwan
             a1 = wan_atom[nw_1]
             for nw_2 = 1:nwan
                 a2 = wan_atom[nw_2]
                 if nonorth
                     S_r[nw_1,nw_2,c] += SR3[nw_1,nw_2,rint[1], rint[2], rint[3]] * sym_R[a1,a2,c]
                 end
                 for spin = 1:nspin
                     ham_r[spin, nw_1,nw_2,c] += hamR3[spin, nw_1,nw_2,rint[1], rint[2], rint[3]] * sym_R[a1,a2,c]
                 end
                 
                 #                println("adding $a1 $a2 $c ", hamR3[a1,a2,rint[1], rint[2], rint[3]], " " , sym_R[a1,a2,c]
 #                        )
             end
         end
     end

     if nonorth
         return ham_r, S_r, r_dict, ind_arr
     else
         return ham_r, r_dict, ind_arr
     end             


 end

 """
     function get_sym_R(crys, grid, sss = 1.0)

 Figures out the r-space grid using Wigner-Seitz like construction to figure out the
 best arrangement of r-grid points to keep periodic copies closest to the original atom 
 and take into account symmetry.

 `return R_grid, R_int_grid, sym_R`

 returns the R_grid, the integer version, and the symmetry factor of each point.

 """
 function get_sym_R(crys, grid, sss = 1.0)


     grid2 = [0,0,0]
     for i = 1:3
         if grid[i] == 1
             grid2[i] = 3
         elseif grid[i]%2 == 0
             grid2[i] = Int(grid[i]/2)+4
         else        
             grid2[i] = Int((grid[i]+1)/2)+4
         end
     end

  #   println("grid2: " , grid2)

     A = crys.A
     At = transpose(crys.A)
 #    At = crys.A
     c=0


     r = zeros(Float64, 1,3)
     R = zeros(Float64, 1, 3)
     rint = zeros(Int64, 1, 3)

     dR = zeros(Float64, 1, 3)

     dist_arr = zeros(5^3)

     ngrid2 = (grid2[1]*2+1)*(grid2[2]*2+1)*(grid2[3]*2+1)
 #    println("ngrid2: $ngrid2")
     sym_R = zeros(Float64,  crys.nat, crys.nat, ngrid2)

     R_grid = zeros(Int64, ngrid2, 3)
     R_int_grid = zeros(Int64, ngrid2, 3)


     for a = 1:crys.nat
         for b = 1:crys.nat        
             c=0
             for r1 = -grid2[1]:(grid2[1])
                 r[1] = r1
 #                rint[1]= Int.(-r[1] .+ 1)
                 for r2 = -grid2[2]:(grid2[2])
                     r[2] = r2
 #                    rint[2]= Int.(-r[2] .+ 1)
                     for r3 = -grid2[3]:(grid2[3])

                         c+=1

                         r[3] = r3
 #                        rint[3]= Int.(-r[3] .+ 1)

 #####                   #prepare Rgrids
                         if (a == 1 && b == 1)

                             R_grid[c,:] = [r1,r2,r3]

                             rint[:] = Int.(r[:] .+ 1)
                             for i = 1:3
                                 if rint[i] <= 0
                                     rint[i] = rint[i] + grid[i]
                                 elseif rint[i] <= 0
                                     rint[i] = rint[i] + grid[i]
                                 elseif rint[i] > grid[i]
                                     rint[i] = rint[i] - grid[i]
                                 end
                             end
                             R_int_grid[c,:] = rint[:]
                         end
 ###############################

 #                        dR[:] = At*(-crys.coords[a,:] + crys.coords[b,:] + r) A
                         dR[:] = (crys.coords[[a],:] - crys.coords[[b],:] + r) * A
                         dist = sum(dR.^2)^0.5
                         cd = 0
 #                        if dist < 7.0
 #                            println("AP $a $b [$r1 $r2 $r3] $dist $dR $r")
 #                        end
                         #                        for R1 = [0]
 #                            for R2 = [0]                            
 #                        for R1 = [-2*grid[1], -grid[1], 0, grid[1], 2*grid[1]]
 #                            for R2 = [-2*grid[1], -grid[2], 0, grid[2], 2*grid[2]]
 #                                for R3 = [-2*grid[3], -grid[3], 0, grid[3], 2*grid[3]]
                         for R1 = -2:2
                             for R2 = -2:2
                                 for R3 = -2:2
                                     R[1,:] = [R1*grid[1] R2*grid[2] R3*grid[3]]
                                     cd += 1
                                     dR[:] = (crys.coords[[a],:] - crys.coords[[b],:] + r + R) * A
                                     dist_arr[cd] = sum(dR.^2)^0.5
                                 end
                             end
                         end
                         if isapprox(dist, minimum(dist_arr))
                             cd = 0
                             for i = 1:125
                                 if isapprox(dist, dist_arr[i])
                                     cd += 1
                                 end
                             end
                             sym_R[a,b,c] = 1.0/Float64(cd)
                         end


                     end
                 end
             end
         end
     end


     for a = 1:crys.nat
         for b = 1:crys.nat        
 #            for i in 1:ngrid2
 #                println("sym_R[$a,$b,$i] = ", sym_R[a,b,i], "  R= ", R_grid[i,:])                
 #            end
             println("SUM $a $b : ", sum(sym_R[a,b,:]))
         end
     end

     keep = []
     nkeep = 0
     for i = 1:ngrid2 
         if maximum(sym_R[:,:,i]) > 1e-5
             push!(keep, i)
             nkeep += 1
         end
     end
     println("nkeep : $nkeep")
     R_grid = R_grid[keep,:]
     R_int_grid = R_int_grid[keep,:]
     sym_R = sym_R[:,:,keep]

     return R_grid, R_int_grid, sym_R

 end    

 """
     function tb_indexes(d::dftout)

 Figures out mapping between DFT projected hamiltonian orbitals and crystal and the wannier orbitals we want.

 `return wan, semicore, nwan, nsemi, wan_atom, atom_wan`

 - `wan` has the indexes of the wannier orbitals
 - `semicore` has the indexes of semicore states.
 - `nwan` number of wannier orbs
 - `nsemi` number of semicore states
 - `wan_atom` dictionary wannier to atom numbers
 - `atom_wan` dictionary atom numbers to wannier orbitals
 """
 function tb_indexes(d::dftout)
     crys = d.crys
     return tb_indexes(crys)
 end

 function tb_indexes(crys::crystal)

     semicore = []
     wan = []

     wan_atom = Dict()
     atom_wan = Dict()

     n = 1
     n2 = 1
     c=0
     for (at_num, t) in enumerate(crys.types)
 #        println(t, " ",atoms[t].nsemicore/2," ",atoms[t].nwan/2)

         n2 += Int(atoms[t].nsemicore/2)
         if n2 > n
             append!(semicore, collect(n:n2-1))
         end
         n = n2
         n2 += Int(atoms[t].nwan/2)
         if n2 > n
             append!(wan, collect(n:n2-1))
         end
         n = n2

         atom_wan[at_num] = []
         for nw = 1:Int(atoms[t].nwan/2)
             c += 1
             wan_atom[c] = at_num
             append!(atom_wan[at_num], c)
         end

     end

     nwan = length(wan)
     nsemi = length(semicore)



     return wan, semicore, nwan, nsemi, wan_atom, atom_wan
 end

 """
     function symm_by_orbitals(crys::crystal, mat)

 Helper function to re-symmeterize the overlap matrix properly, starting from incomplete k-points.
 """
 function symm_by_orbitals(crys::crystal, mat)
     c=0
     for (at_num, t) in enumerate(crys.types)
         for o in atoms[t].orbitals
             if o == :s
                 c += 1
             elseif o == :p
                 t = sum(mat[c+1:c+3])/3.0
                 mat[c+1:c+3] .= t
                 c += 3
             elseif o == :d
                 t = sum(mat[c+1:c+5])/5.0
                 mat[c+1:c+5] .= t
                 c += 5
             elseif o == :f
                 t = sum(mat[c+1:c+7])/7.0
                 mat[c+1:c+7] .= t
                 c += 7
             else
                 println("WARNING symm_by_orbitals expecting orbital got : $o")
                 c += 1
             end
         end
     end
     return mat
 end



 """
     function find_vbm_cbm(eigs, fermi)

 Find the valence band max and conduction band minimum from eigs, relative to Fermi level.
 """
 function find_vbm_cbm(eigs, fermi)

     vbm = -100000.0
     cbm = 100000.0
     for eig in eigs[:]
         if eig < fermi && eig  > vbm
             vbm = eig
         end

         if eig > fermi && eig  < cbm
             cbm = eig
         end
     end
     return vbm, cbm
 end


 """
     function ewald_energy(tbc::tb_crys, delta_q=missing)

 Return ewald energy term from tbc. If `delta_q`, the atomic charge density, is missing,
 loads from `tbc`.
 """
 function ewald_energy(tbc::tb_crys, delta_q=missing)

     gamma = tbc.gamma 
     crys = tbc.crys

     if ismissing(delta_q)
         delta_q =  get_dq(crys , tbc.eden)
     end

     return ewald_energy(crys, gamma, delta_q)

 end

 """
     function ewald_energy(tbc::tb_crys_kspace, delta_q=missing)
 """
 function ewald_energy(tbc::tb_crys_kspace, delta_q=missing)

     gamma = tbc.gamma 
     crys = tbc.crys

     if ismissing(delta_q)
         delta_q =  get_dq(crys , tbc.eden)
     end

     return ewald_energy(crys, gamma, delta_q)

 end

 """
     function ewald_energy(crys::crystal, gamma, delta_q::Array{Float64,1})

 Does the actual calculation.
 """
 function ewald_energy(crys::crystal, gamma, delta_q::Array{Float64,1})

     T = typeof(crys.coords[1,1])
     pot = zeros(T, crys.nat, crys.nat)

     for i = 1:crys.nat
         for j = 1:crys.nat
             pot[i,j] = gamma[i,j] * delta_q[i] * delta_q[j] 
         end
     end


     energy = 0.5*sum(pot)

 #    println("ewald_energy ", energy, " " , delta_q, " ", gamma[1,1], " ", gamma[1,2], " ", gamma[2,1], " ", gamma[2,2])

     return energy, pot

 end

 """
     function get_neutral_eden(tbc::tb_crys)

 Gets a neutral charge density (no charge transfer) to start SCF calculation.
 """
 function get_neutral_eden(tbc::tb_crys)

     return get_neutral_eden(tbc.crys, tbc.tb.nwan)

 end

 """
     function get_neutral_eden(crys::crystal, nwan=missing)
 """
 function get_neutral_eden(crys::crystal, nwan=missing)

     if ismissing(nwan)
         ind2orb, orb2ind, etotal, nval = orbital_index(crys)
         nwan = length(keys(ind2orb))
     end

     eden = zeros(nwan)
     counter = 0
     for (i, t) in enumerate(crys.types)
         at = atoms[t]
         z_ion = at.nval
         still_need = z_ion
         for o in at.orbitals
             if o == :s
                 counter += 1
                 if still_need <= 2.0 && still_need > 1e-5
                     eden[counter] = still_need/2.0
                     still_need = 0.0
                 elseif still_need >= 2.0 && still_need > 1e-5
                     eden[counter] = 1.0
                     still_need = still_need - 2.0
                 else
                     counter += 1
                 end
             elseif o == :p
 #                println("p")
                 if still_need <= 6.0 && still_need > 1e-5
                     eden[counter+1] = (still_need/2.0)/3.0
                     eden[counter+2] = (still_need/2.0)/3.0
                     eden[counter+3] = (still_need/2.0)/3.0
                     still_need = 0.0
                     counter += 3
                 elseif still_need >= 6.0 && still_need > 1e-5
                     eden[counter+1] = 1.0
                     eden[counter+2] = 1.0
                     eden[counter+3] = 1.0
                     still_need = still_need - 6.0
                     counter += 3
                 else
                     counter += 3
                 end
             elseif  o == :d
                 if still_need <= 10.0 && still_need > 1e-5
                     eden[counter+1] = (still_need/2.0)/5.0
                     eden[counter+2] = (still_need/2.0)/5.0
                     eden[counter+3] = (still_need/2.0)/5.0
                     eden[counter+4] = (still_need/2.0)/5.0
                     eden[counter+5] = (still_need/2.0)/5.0
                     still_need = 0.0
                     counter += 5
                 elseif still_need >= 10.0 && still_need > 1e-5
                     eden[counter+1] = 1.0
                     eden[counter+2] = 1.0
                     eden[counter+3] = 1.0
                     eden[counter+4] = 1.0
                     eden[counter+5] = 1.0
                     still_need = still_need - 10.0
                     counter += 5
                 else
                     counter += 5
                 end
             else
                 println("bad orbital get_neutral_eden $o")
             end
         end
     end
     return eden

 end

 """
     function get_dq(tbc::tb_crys_kspace)

 Get atomic charge density from `tb_crys` or `tb_crys_kspace` or `crys + eden`
 """
 function get_dq(tbc::tb_crys_kspace)
     return get_dq(tbc.crys, tbc.eden)
 end

 """
     function get_dq(tbc::tb_crys)
 """
 function get_dq(tbc::tb_crys)
     return get_dq(tbc.crys, tbc.eden)
 end

 """
     function get_dq(crys::crystal, chargeden::Array{Float64,1})
 """
 function get_dq(crys::crystal, chargeden::Array{Float64,2})

     nspin = size(chargeden)[1]

     e_den = zeros(Float64, crys.nat)
     z_ion = zeros(Float64, crys.nat)

     for spin = 1:nspin
         counter = 0
         for (i, t) in enumerate(crys.types)
             at = atoms[t]
             z_ion[i] = at.nval
             for o = 1:Int64(at.nwan/2)
                 counter += 1
                 e_den[i] += chargeden[spin,counter] 
             end
         end
     end
     if nspin == 1
         e_den = e_den * 2.0
     end

     
     dq = -z_ion + e_den

     dq = dq .- sum(dq)/crys.nat #charge sum rule

 #    println("e_den ", e_den)
 #    println("z_ion ", z_ion)
 #    println("dq ", dq)


     return dq

 end

 """
     function get_h1(tbc::tb_crys)

 Get H1, the potential term added to tight binding in SCF calculation.
 """
 function get_h1(tbc::tb_crys)
     return get_h1(tbc, tbc.eden)
 end


 """
     function get_h1(tbc::tb_crys, chargeden::Array{Float64,1})
 """
 function get_h1(tbc::tb_crys, chargeden::Array{Float64,1})

     dq = get_dq(tbc.crys, chargeden)

     gamma = tbc.gamma

     epsilon = gamma * dq

     h1 = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)
     o1 = 1
     for i = 1:tbc.crys.nat
         at1 = atoms[tbc.crys.types[i]  ]
         nw1 = Int64(at1.nwan/2)
         o2 = 1
         for j = 1:tbc.crys.nat
             at2 = atoms[tbc.crys.types[j]]
             nw2 = Int64(at2.nwan/2)
             for c1 = o1:o1+nw1-1
                 for c2 = o2:o2+nw2-1
                     h1[c1,c2] = 0.5 * (epsilon[i] + epsilon[j])
                 end
             end
             o2 += nw2

         end
         o1 += nw1
     end

     return 0.5*(h1 + h1'), dq

 end


 """
     function get_h1(tbc::tb_crys_kspace)
 """
 function get_h1(tbc::tb_crys_kspace)
     return get_h1(tbc, tbc.eden)
 end

 """
     function get_h1(tbc::tb_crys_kspace, chargeden::Array{Float64,1})
 """
 function get_h1(tbc::tb_crys_kspace, chargeden::Array{Float64,1})

     dq = get_dq(tbc.crys, chargeden)

     gamma = tbc.gamma

     epsilon = gamma * dq

     h1 = zeros(Complex{Float64}, tbc.tb.nwan, tbc.tb.nwan)
     o1 = 1
     for i = 1:tbc.crys.nat
         at1 = atoms[tbc.crys.types[i]  ]
         nw1 = Int64(at1.nwan/2)
         o2 = 1
         for j = 1:tbc.crys.nat
             at2 = atoms[tbc.crys.types[j]]
             nw2 = Int64(at2.nwan/2)
             for c1 = o1:o1+nw1-1
                 for c2 = o2:o2+nw2-1
                     h1[c1,c2] = 0.5 * (epsilon[i] + epsilon[j])
                 end
             end
             o2 += nw2

         end
         o1 += nw1
     end

     return 0.5*(h1 + h1'), dq

 end

 """
     function get_energy_electron_density_kspace(tbcK::tb_crys_kspace; smearing = 0.01)

 Get energy / charge density from k-space tight binding object.

 `return bandenergy + etypes + echarge + energy_smear, eden, VECTS, VALS, error_flag`

 """
 function get_energy_electron_density_kspace(tbcK::tb_crys_kspace; smearing = 0.01)

     bandenergy, eden, VECTS, VALS, efermi, error_flag = get_energy_electron_density_kspace(tbcK.tb, tbcK.nelec, smearing=smearing)

     if tbcK.nspin == 1
         tbcK.eden[1,:] = eden[1,:]
     elseif tbcK.nspin == 2
         tbcK.eden[1,:] = eden[1,:]
         tbcK.eden[2,:] = eden[2,:]
     end
     
     if tbcK.scf
 #        h1 = tbc.tb.h1
         echarge, pot = ewald_energy(tbcK)
     else
 #        h1 = missing
         echarge = 0.0
     end

     etypes = types_energy(tbcK.crys)

     println("efermi $efermi")
     
     energy_smear = smearing_energy(VALS, tbcK.tb.kweights, efermi, smearing)
     println("CALC ENERGIES $etypes $echarge $bandenergy $energy_smear = ", bandenergy + etypes + echarge + energy_smear)

     return bandenergy + etypes + echarge + energy_smear, eden, VECTS, VALS, error_flag


 end

 """
     function get_energy_electron_density_kspace(tb_k::tb_k, nelec; smearing = 0.01)

 K-space get energy and electron density from `tb_k`
 """
 function get_energy_electron_density_kspace(tb_k::tb_k, nelec; smearing = 0.01)

     temp = zeros(Complex{Float64}, tb_k.nwan, tb_k.nwan)
     denmat = zeros(Float64, tb_k.nspin, tb_k.nwan, tb_k.nwan)

     VALS = zeros(Float64, tb_k.nk,  tb_k.nwan,tb_k.nspin)
     VALS0 = zeros(Float64, tb_k.nk,  tb_k.nwan,tb_k.nspin)
     VECTS = zeros(Complex{Float64}, tb_k.nk, tb_k.nspin, tb_k.nwan, tb_k.nwan)
     SK = zeros(Complex{Float64}, tb_k.nk,  tb_k.nwan, tb_k.nwan)

     error_flag = false

     for spin = 1:tb_k.nspin
         for k in 1:tb_k.nk
#             try
                 vects, vals, hk, sk, vals0 = Hk(tb_k, tb_k.K[k,:], spin=spin)
                 VALS[k, :, spin] = vals
                 VECTS[k,spin, :,:] = vects
                 SK[k,:,:] = sk
                 VALS0[k, :,spin] = vals0
#             catch e
#                 if e isa InterruptException
#                     println("user interrupt")
#                     rethrow(InterruptException)
#                 end
#                 println("warning, get_energy_electron_density_kspace error")
#                 error_flag = true
#             end
         end
     end
     
     energy, efermi, occs = band_energy(VALS, tb_k.kweights, nelec, smearing, returnboth=true)

     energy0 = sum(occs .* VALS0 .* tb_k.kweights) / sum(tb_k.kweights)
     if tb_k.nspin == 1
         energy0 = energy0 * 2
     end
         

         
     for spin = 1:tb_k.nspin
         for k in 1:tb_k.nk
             for a = 1:tb_k.nwan
                 for i = 1:tb_k.nwan
                     for j = 1:tb_k.nwan
                         temp[i,j] = VECTS[k,spin, i,a]' * SK[k,i,j] * VECTS[k,spin, j,a]
                     end
                 end
                 temp = temp + conj(temp)
                 denmat[spin,:,:] += 0.5 * occs[k,a, spin] * real.(temp) * tb_k.kweights[k]

             end
         end
     end

     
#     println("denmat")
#     println(denmat)
     
     
     electron_den = sum(denmat, dims=2) / sum(tb_k.kweights)
     electron_den = electron_den[:,1,:]
     if tb_k.nspin == 2
         electron_den = electron_den*2.0
         energy0 = energy0
     end


     toobig = electron_den .> 1.0
     normal = electron_den .< 1.0
     electron_den[toobig] .= 1.0
     
     if tb_k.nspin == 2
         tot = sum(electron_den)
         correction = tot - nelec
     elseif tb_k.nspin == 1
         tot = sum(electron_den) * 2.0
     end
     correction = tot - nelec
     electron_den[normal] .-= correction / sum(normal)

     
     return energy0, electron_den, VECTS, VALS, efermi, error_flag

 end






end #end module

##include("Crystal.jl")
###include("DFToutMod.jl")
#using XMLDict



####################### QE specific 
"""
    module AtomicProj

Analyze projections from QE projwfc.x
"""
module AtomicProj
"""
Scripts to analyze atomic projections
"""


#export load_xml
using LinearAlgebra
using DelimitedFiles

#include("Atomdata.jl")
using ..Atomdata:atoms
#using JLD
#using ExXML
using XMLDict
using GZip

using FFTW
using ..DFToutMod
#using ..CrystalMod:crystal
#using ..Utility:arr2str
#using ..Utility:str_w_spaces
using ..Utility:parse_str_ARR_float
using ..Utility:parse_str_ARR_complex
using ..Utility:write_to_file
using ..DFToutMod:bandstructure
using ..TB:make_tb
using ..TB:make_tb_crys
using ..BandTools:band_energy
#using ..TB:orbital_index
using ..CrystalMod:orbital_index

using ..TB:write_tb_crys
#using ..TB:orbital_index
using ..TB:trim
using ..QE:loadXML_bs
using ..QE:loadXML
using ..TB:calc_energy
using ..TB:types_energy
using ..CrystalMod:crystal
using ..DFT:runSCF
using ..CrystalMod:get_grid
using ..TB:myfft
using ..TB:tb_indexes
using ..TB:symm_by_orbitals

using ..TB:tb_crys_kspace
using ..TB:make_tb_crys_kspace
using ..TB:make_tb_k
using ..TB:write_tb_crys_kspace


include("Commands.jl")
using ..ThreeBodyTB:TEMPLATES




"""
    mutable struct proj_dat

Holds data from projwfc.x

- `bs::bandstructure`
- `natwfc::Int64` number of atomic wavefunctions, from dft projection
- `proj::Array{Complex{Float64}, 3}`  atomc projections `nk × natwfc × nbnd` , where `nbnd` is the number of bands in DFT
- `overlaps::Array{Complex{Float64}, 3} overlap matrix from dft `nk × natwfc × natwfc`
- `nspin::Int64` number of psins
This is created from `loadXML_proj`, which calls `make_proj`
"""
mutable struct proj_dat

    bs::bandstructure
    natwfc::Int64
    nspin::Int64
    proj::Array{Complex{Float64}, 4}
    overlaps::Array{Complex{Float64}, 3}    
end

#print proj dat
Base.show(io::IO, d::proj_dat) = begin
    println(io,"projection data: nbnd = ", d.bs.nbnd, "; nkpts = ", d.bs.nks, "; natwfc = ", d.natwfc, "; nspin = ", d.nspin)
end   



"""
    function makedict_proj(savedir)

Return xml proj data from QE savedir
"""
function makedict_proj(savedir)

    if isfile(savedir*"/atomic_proj.xml")
        filename=savedir*"/atomic_proj.xml"
    elseif isfile(savedir*"/atomic_proj.xml.gz")
        filename=savedir*"/atomic_proj.xml.gz"
    else
        println("error warning no atomic_proj.xml.gz or atomic_proj.xml")
        filename = missing
    end
    
    f = gzopen(filename, "r")
    fs = read(f, String)
    close(f)
    
    d = xml_dict(fs)

    return d
end

"""
    function make_projwfcx(prefix, tmpdir)

Create the inputfile for running projwfc.x from QE
"""
function make_projwfcx(prefix, tmpdir)
"""
makes the proj input file
"""
    c_dict = make_commands(1)


    template_file=open("$TEMPLATES/template.proj")
    temp = read(template_file, String)
    close(template_file)
    
    temp = replace(temp, "TMPDIR" => tmpdir)
    temp = replace(temp, "PREFIX" => prefix)
    
    return temp
    
end


"""
    function run_projwfcx(projfile="proj.in"; directory="./", nprocs=1)

Run projwfc.x from QE code
"""
function run_projwfcx(projfile="proj.in"; directory="./", nprocs=1)
"""
run projwfc.x QE command
"""

#    nprocs = 1
#    println("set nprocs $nprocs")
    
    c_dict = make_commands(nprocs)
    proj = c_dict["proj"]
    command = `$proj $directory/$projfile `
    println("projwfc.x command 1")
    println(command)
    println()
    flush(stdout)
    try
        println("Running projwfc.x")
        s = read(command, String)
        f = open(directory*"/"*projfile*".out", "w")
        write(f, s)
        close(f)
        
        println("Ran projwfc.x")
        println()
    catch

        try
                proj = c_dict["proj_serial"]
                command = `$proj $directory/$projfile `
                println("projwfc.x command 4")
                println(command)
                flush(stdout)
                s = read(command, String)
                f = open(directory*"/"*projfile*".out", "w")
                write(f, s)
                close(f)
                flush(stdout)
                command = `$proj $directory/$projfile `
        catch
            println("failed proj twice")
            try
                proj = c_dict["proj_serial_backup"]
                command = `$proj $directory/$projfile `
                println("projwfc.x command 4")
                println(command)
                flush(stdout)
                s = read(command, String)
                f = open(directory*"/"*projfile*".out", "w")
                write(f, s)
                close(f)
                flush(stdout)
                command = `$proj $directory/$projfile `

            catch
                try
                    proj = c_dict["proj_backup"]
                    command = `$proj $directory/$projfile `
                    println("projwfc.x command 4")
                    println(command)
                    flush(stdout)
                    s = read(command, String)
                    f = open(directory*"/"*projfile*".out", "w")
                    write(f, s)
                    close(f)
                    flush(stdout)
                    command = `$proj $directory/$projfile `

                catch
                    println("all proj failed")
                end
            end
        end
        

    end
    
    
    return -1
end


"""
    function makeOG(prefix, tmpdir )

open_grid.x inputfile. My workflow doesn't currently use this, as open_grid.x was unreliable, and
I've rewritten the code to avoid it and work directly with independent k-points only.
"""
function makeOG(prefix, tmpdir )
"""
makes the OG file
"""
    template_file=open("$TEMPLATES/template_og.in")
    temp = read(template_file, String)
    close(template_file)
    
    temp = replace(temp, "TMPDIR" => tmpdir)
    temp = replace(temp, "PREFIX" => prefix)
    
    return temp
    
end

"""
    function run_og(filename="og.in";  directory="./", nprocs=1)

runs open_grid.x
"""
function run_og(filename="og.in";  directory="./", nprocs=1)
"""
run open_grid.x command
"""
    c_dict = make_commands(nprocs)
    og = c_dict["og"]
    command = `$og $directory/$filename`
    println("opengrid command")
    println(command)
    try
        println("Running open_grid.x")

        s = read(command, String)
        f = open(directory*"/"*filename*".out", "w")
        write(f, s)
        close(f)
        
        println("Ran open_grid.x")
        return 0
    catch
        println("Failed to run open_grid.x")
        return -1
    end
        
end


function run_nscf(dft, directory; tmpdir="./", nprocs=1, prefix="qe", min_nscf=false, only_kspace=false, klines=missing, gamma_only=false)

    olddir = "$directory/$prefix.save"
    println("doing nscf")
    nscfdir = "$directory/$prefix.nscf.save"
    if isdir(nscfdir)
        println("removing original nscf directory $nscfdir")
        rm(nscfdir,recursive=true)
    end
    mkpath(nscfdir)
    #        println("after mkdpath")
    
    files = readdir(olddir)
    for f in files
        if length(f) > 3 && f[end-2] == 'U' && f[end-1] == 'P' && f[end] == 'F'
            #                println("found $f , copy")
            cp("$olddir/$f", "$nscfdir/$f")            
        end
    end
    
    try
        if isfile("$olddir/charge-density.dat")
            cp("$olddir/charge-density.dat",  "$nscfdir/charge-density.dat")
        end
        if isfile("$olddir/data-file-schema.xml")
            cp("$olddir/data-file-schema.xml", "$nscfdir/data-file-schema.xml" )
        end
        if isfile("$olddir/data-file-schema.xml.gz")
            cp("$olddir/data-file-schema.xml.gz", "$nscfdir/data-file-schema.xml.gz" )
            
            tounzip = "$nscfdir/data-file-schema.xml.gz"
            command = `gunzip $tounzip`
            s = read(command, String)
            println("gunzipped $tounzip")
        end
        if isfile("$olddir/charge-density.dat.gz")
            cp("$olddir/charge-density.dat.gz", "$nscfdir/charge-density.dat.gz" )
            
            tounzip = "$nscfdir/charge-density.dat.gz"
            command = `gunzip $tounzip`
            s = read(command, String)
            println("gunzipped $tounzip")
            
        end
        if isfile("$olddir/charge-density.hdf5")
            cp("$olddir/charge-density.hdf5", "$nscfdir/charge-density.hdf5" )
        end
        if isfile("$olddir/charge-density.hdf5.gz")
            cp("$olddir/charge-density.hdf5.gz", "$nscfdir/charge-density.hdf5.gz" )
            
            
            tounzip = "$nscfdir/charge-density.hdf5.gz"
            command = `gunzip $tounzip`
            s = read(command, String)
            println("gunzipped $tounzip")
            
            
        end

        
        
    catch
        println("missing charge density or xml file, cannot run nscf")
    end
    
    crys = dft.crys
    tot_charge = dft.tot_charge
    grid=dft.bandstruct.kgrid
    
    if gamma_only
        grid = [1,1,1]
    elseif min_nscf
        println("minimize kgrid")
        grid = max.(grid .- 2, 2)
    end
                

    if only_kspace
        calc = "nscf-sym"
    else
        calc = "nscf"
    end

    dft_nscf = missing

    magnetic = false
    if dft.nspin == 2
        magnetic = true
    end
    
    try
        dft_nscf = runSCF(crys, prefix="$prefix.nscf", directory=directory,tmpdir=directory, wannier=2, nprocs=nprocs, skip=false, calculation=calc, tot_charge=tot_charge, grid=grid, klines=klines, magnetic=magnetic)
    catch
        println()
        println("first nscf failed, trying backup nscf with fewer extra bands")
        println()
        try
            dft_nscf = runSCF(crys, prefix="$prefix.nscf", directory=directory,tmpdir=directory, wannier=1, nprocs=nprocs, skip=false, calculation=calc, tot_charge=tot_charge, use_backup=true, grid=grid, klines=klines, magnetic=magnetic)
        catch
            println("try 2")
            dft_nscf = runSCF(crys, prefix="$prefix.nscf", directory=directory,tmpdir=directory, wannier=-1, nprocs=nprocs, skip=false, calculation=calc, tot_charge=tot_charge, use_backup=true, grid=grid, klines=klines, magnetic=magnetic)
        end
    end
    
    dft_nscf.energy = dft.energy
    dft_nscf.energy_smear = dft.energy_smear
    dft_nscf.forces = dft.forces
    dft_nscf.stress = dft.stress        
    
    prefix = prefix*".nscf"

    return dft_nscf, prefix
    
end




"""
    function projwfc_workf(dft::dftout)

This is the main workflow for the creation of TB matrix elements from QE DFT calculations.

Starting from a converged QE scf calculation...

1) Run NSCF calculation with extra empty bands. If you want a real space TB object, these need to be a full k-point grid not using symmtery. k-space only can use irreducible k-points

2) Run projwfc.x

3) Construct the TB hamiltonian in k-space

4) (Optional) FT to real space

5) (Optional) cleanup wavefunctions.

# Arguments

- `dft::dftout` The starting scf calculation
- `directory="./"`
- `nprocs=1` number of processors
- `freeze=true` Keep occupied eigenvalues fixed to exact DFT values
- `writefile="projham.xml"` output file for real-space TB
- `writefilek="projham_K.xml"` output file for k-space TB
- `skip_og=true`  Not used anymore
- `skip_proj=true` If projections are already run, don't run them again, load from file.
- `shift_energy=true` Shift energy of eigenvalues s.t. total energy equal band energy.
- `cleanup=true` Remove the large wavefunction files from disk, keeping nscf/projection files.
- `skip_nscf=true` If nscf calculation is already done, load results from file.
- `localized_factor = 0.15` Adjust extent of overlap matrix. 0.0 uses full atomic wavefunctions, which can be overly delocalized. 1.0 is fully localized.
- `only_kspace=false` Do not create real-space tb. Usually true in current code, as I can fit directly from k-space tb only.
- `screening = 1.0` If use a screening factor to reduce value of U in Ewald calculation. Usually leave at 1.0.
"""
function projwfc_workf(dft::dftout; directory="./", nprocs=1, freeze=true, writefile="projham.xml",writefilek="projham_K.xml", skip_og=true, skip_proj=true, shift_energy=true, cleanup=true, skip_nscf=true, localized_factor = 0.15, only_kspace=false, screening = 1.0, min_nscf=false, gamma_only=false)
"""

Steps:

2) run projwfc.x
2a) load output
3) make projected ham in k-space
4) do fourier transform
5) optionally write results to file
"""
    prefix=dft.prefix
    outdir=dft.outdir
    println()
    println("projwfc_workflow---------------------------------------------------------------------------")
    println()
    println("Step #0 nscf---------------------------------------------------------------------------------")
    s=""
    prefix_orig = deepcopy(prefix)

#    try

    olddir = "$directory/$prefix.save"
    nscfdir = "$directory/$prefix.nscf.save"

    
    if !(skip_nscf) || !(isdir(nscfdir)) ||  ( !isfile(nscfdir*"/atomic_proj.xml") && !isfile(nscfdir*"/atomic_proj.xml.gz"))  #change

        dft_nscf, prefix = run_nscf(dft, directory; tmpdir=directory, nprocs=nprocs, prefix=prefix, min_nscf=min_nscf, only_kspace=only_kspace, gamma_only=gamma_only)

#        println("DO COPY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! kfg")
#        cp("/home/kfg/codes/TB_run/optim/testing2/t.UPF" , "$nscfdir/h.pbesol.UPF", force=true)
        
    else
        if (isdir(nscfdir))
            crys = dft.crys
            #            dft_nscf = runSCF(crys, prefix="$prefix.nscf", directory=directory,tmpdir=directory, wannier=true, nprocs=nprocs, skip=false, calculation="nscf")
            dft_nscf = loadXML("$directory/$prefix.nscf.save")            
            dft_nscf.energy = dft.energy
            dft_nscf.forces = dft.forces
            dft_nscf.stress = dft.stress        
            prefix = prefix*".nscf"
            println("loaded nscf aaaaaaaaaaaaaaaaaaaaaa")
        else
            crys = dft.crys
            dft_nscf = dft
        end



    end
        
#    catch
#        println("failed nscf")
#        println(s)
#    end

#    return

    if maximum(abs.(dft_nscf.bandstruct.eigs[:,1,:] - dft.bandstruct.eigs[:,1,:])) > 1e-2
        println("dft_nscf and dft eigenvalues do not match. catastrophic error. we refuse to continue")
        throw("dft_nscf and dft eigenvalues do not match. catastrophic error. we refuse to continue")
    end

    newprefix = prefix
    
    useog = false
    if useog
        println("Step #1 og---------------------------------------------------------------------------------")
        println()

        newprefix=prefix*"_open"
        
        function og()
            og= makeOG(prefix, outdir)
            write_to_file(og, "og.in", directory)
            run_og("og.in", directory=directory, nprocs=nprocs)
        end
        
        if skip_og
            try
                dftopen_bs = loadXML_bs("$directory/$newprefix.save") #get unfolded k-points and bs
                println("using precomputed opengrid dir")
            catch
                println("precomputinga opengrid dir")
                og()
                dftopen_bs = loadXML_bs("$directory/$newprefix.save") #get unfolded k-points and bs
                
            end
        else
            og()
            dftopen_bs = loadXML_bs("$directory/$newprefix.save") #get unfolded k-points and bs
        
        end

#        kgrid = copy(dft.bandstruct.kgrid)
#        dft_nscf.bandstruct = dftopen_bs
#        dft_nscf.bandstruct.kgrid = kgrid

    else
        println("DON'T USE OG, since it works poorly")
        dftopen_bs = loadXML_bs("$directory/$newprefix.save")
    end

    

        
    println()
    println("Step #2 run projwfc.x----------------------------------------------------------------------")
    println()

    println("A dft")
    println(dft.crys.A)
    println("A dft nscf")
    println(dft_nscf.crys.A)

    A = dft.crys.A             #needed to convert kpoints to crystal units
    a1=sum(A[1,:].^2)^0.5
    B = inv(A./a1)'
    
    function proj()
        projstr = make_projwfcx(newprefix, directory)
        write_to_file(projstr, "proj.in", directory)
        run_projwfcx("proj.in", directory=directory, nprocs=nprocs)
    end
    
    if skip_proj
        try
            println("using precomputed proj file $directory/$newprefix.save $B")
            p = loadXML_proj("$directory/$newprefix.save", B)
            println("using precomputed proj file")
        catch
            println("computing proj file")
            proj()
            p = loadXML_proj("$directory/$newprefix.save", B)
        end
    else
        proj()
        p = loadXML_proj("$directory/$newprefix.save", B)
        
    end
    
    if sum(abs.(dft_nscf.crys.A - dft.crys.A)) > 1e-4
        println("Warning, dft and dft_nscf may not match")
        println(dft.crys)
        println(dft_nscf.crys)
    end

    p.bs.kpts[:,:] = dft_nscf.bandstruct.kpts #fix kpoint rounding issue

    println("P OVERLAPS ", size(p.overlaps))
    
    #check if p matches dft_nscf
    if sum(abs.(dft_nscf.bandstruct.eigs[1,:,1] - p.bs.eigs[1,:,1])) > 1e-5
        println("warning, projected eigs and dft_nscf eigs do not match ", sum(abs.(dft_nscf.bandstruct.eigs[1,:,1] - p.bs.eigs[1,:,1])))
        println("rerun nscf")
        prefix = deepcopy(prefix_orig)
        dft_nscf = run_nscf(dft, directory; tmpdir=directory, nprocs=nprocs, prefix=prefix, min_nscf=min_nscf, only_kspace=only_kspace)
#        dft_nscf = run_nscf()
        println("rerun proj")
        proj()
        p = loadXML_proj("$directory/$newprefix.save", B)

    end


    println()
    println("Step #3 ham_k------------------------------------------------------------------------------")
    println()
    projection_warning = false
    if freeze

        
        band_froz = Int64(round(dft_nscf.bandstruct.nelec/2))+1


        en_froz = minimum(dft_nscf.bandstruct.eigs[:,band_froz,:])
        en_froz = max(en_froz, dft.bandstruct.efermi + 0.05)
        println("en_froz: ", en_froz, " band froz $band_froz")
        println("efermi energy ", dft.bandstruct.efermi)
        println("nk ", size(dft_nscf.bandstruct.kpts), "     ", size(dft_nscf.bandstruct.kweights))
        ham_k, EIG, Pmat, Nmat, VAL, projection_warning = AtomicProj.create_tb(p, dft_nscf, energy_froz=en_froz+.05, shift_energy=shift_energy);
        #A,B,C = AtomicProj.create_tb(p, dft, energy_froz=en_froz+.01); 
        #return A,B,C        
    else
        ham_k, EIG, Pmat, Nmat, VAL, projection_warning = AtomicProj.create_tb(p, dft_nscf);
    end

#    println("done")
#    return 0

    println()
    println("Step #4 ham_k/ham_r------------------------------------------------------------------------------")
    println()

    println("after create P OVERLAPS ", size(p.overlaps))


    #get tight binding
    tbck = prepare_ham_k(p, dft_nscf, dft_nscf.bandstruct.kgrid ,ham_k, nonorth=true, localized_factor = localized_factor, screening=screening)

#    return tbck
    
    println(tbck.nspin, " after prepare P OVERLAPS ", size(p.overlaps))


    if !ismissing(writefilek)
        println("Step #5-0 write_tb_k-------")
        write_tb_crys_kspace("$directory/$writefilek", tbck);
    end

    if !only_kspace
        tbc = AtomicProj.get_ham_r(tbck, dft.tot_charge)
        if !ismissing(writefile)
            println()
            println("Step #5 write_tb---------------------------------------------------------------------------")
            println()
            write_tb_crys("$directory/$writefile", tbc);
        end
    else
        tbc = missing
    end



#    tbc = AtomicProj.get_ham_r(p, dft_nscf, dft_nscf.bandstruct.kgrid ,ham_k, nonorth=true, localized_factor = localized_factor);

    #tbc = AtomicProj.get_ham_r( p, dft_nscf, dft_nscf.bandstruct.kgrid ,ham_k, nonorth=false);    


    #tbc = AtomicProj.get_ham_r_slow(p, dft_nscf, dft.bandstruct.kgrid ,ham_k, nonorth=true);    


#    trim(tbc.tb, 0.00005)
    
    println()
    println("Done proj----------------------------------------------------------------------------------")
    println()

    if cleanup
        println("clean: warning, removing wavefunctions to save space")
        try
            newdir = "$directory/$newprefix.save"
            olddir = "$directory/$prefix.save"
            parentdir = "$directory"

            if (prefix_orig != prefix) && (newprefix != prefix_orig)
                origdir = "$directory/$prefix_orig.save"
                toclean = [newdir,olddir,origdir]
            else
                toclean = [newdir,olddir]
            end

            if (isdir(nscfdir))
                toclean = [toclean; nscfdir]
            end
            if (isdir(parentdir))
                toclean = [toclean; parentdir]
            end


            tot = 0
            for d in toclean
                files=read(`ls $d/`, String)
                for f in split(files, "\n")
                    if occursin("wfc", f)
                        if tot <= 2
                            println("removing $d/$f")
                        elseif tot==3
                            println("removing $d/$f")
                            println("...")
                        end
                        try
                            cmd=`\rm $d/$f`
                            s = read(cmd, String)
                            tot += 1
                        catch
                            println("failed to delete $d/$f")
                        end
                    end
                end
            end
            println("clean-up done: $tot wfc files deleted")
        catch
            println("catch, something wrong with cleanup")
        end
    else
        println("don't clean wfcs")
    end        
    
    return tbc, tbck, projection_warning
    
end
    


"""
    function loadXML_proj(savedir, B=missing)

Load proj from QE output file. Need the QE save dir like "qe.save". `B` are reciprocal lattice vectors.
"""
function loadXML_proj(savedir, B=missing)
    
    d= makedict_proj(savedir)

#    println("KEYS : ", keys(d))

    da = missing
    if "ATOMIC_PROJECTIONS" in keys(d)
        da = d["ATOMIC_PROJECTIONS"] #has all the real data in the xml
    else
        da = d["PROJECTIONS"] #has all the real data in the xml
    end

    if "NUMBER_OF_BANDS" in keys(da["HEADER"])
    
        nbnd = parse(Int, da["HEADER"]["NUMBER_OF_BANDS"][""])
        nk = parse(Int, da["HEADER"]["NUMBER_OF_K-POINTS"][""])    
        nspin = parse(Int, da["HEADER"]["NUMBER_OF_SPIN_COMPONENTS"][""])
        natwfc = parse(Int, da["HEADER"]["NUMBER_OF_ATOMIC_WFC"][""])
        nelec = parse(Float64, da["HEADER"]["NUMBER_OF_ELECTRONS"][""])
        efermi = parse(Float64, da["HEADER"]["FERMI_ENERGY"][""])                
        
#        units_energy = da["HEADER"]["UNITS_FOR_ENERGY"][:UNITS]
#        units_kpt = da["HEADER"]["UNITS_FOR_K-POINTS"][:UNITS]    

    else
        
        nbnd = parse(Int, da["HEADER"][:NUMBER_OF_BANDS])
        nk = parse(Int, da["HEADER"][Symbol("NUMBER_OF_K-POINTS")])    
        nspin = parse(Int, da["HEADER"][:NUMBER_OF_SPIN_COMPONENTS])
        natwfc = parse(Int, da["HEADER"][:NUMBER_OF_ATOMIC_WFC])
        nelec = parse(Float64, da["HEADER"][:NUMBER_OF_ELECTRONS])
        efermi = parse(Float64, da["HEADER"][:FERMI_ENERGY])                
        
#        units_energy = da["HEADER"]["UNITS_FOR_ENERGY"][:UNITS]
#        units_kpt = da["HEADER"]["UNITS_FOR_K-POINTS"][:UNITS]    

    end
    
#    println(units_energy, " " , units_kpt)
#    println("nbands: ", nbnd," nk: ", nk)

    proj = zeros(Complex{Float64}, nk, natwfc, nspin, nbnd)
    eigs = zeros(nk,nbnd, nspin)

    kpts = zeros(Float64, nk,3)
    weights = zeros(Float64, nk)

    overlaps = zeros(Complex{Float64}, nk, natwfc, natwfc)

    
    if "K-POINTS" in keys(da) #old format

        kpts = parse_str_ARR_float(da["K-POINTS"][""])
        weights = parse_str_ARR_float(da["WEIGHT_OF_K-POINTS"][""])

        dap = da["PROJECTIONS"]
        dae = da["EIGENVALUES"]
        dao = da["OVERLAPS"]
        kpts = reshape(kpts, 3,nk)'
        

        if nspin == 1
            for k in 1:nk
                t = dap["K-POINT."*string(k)]
                for a in 1:natwfc
                    proj[k, a,1, :] =  parse_str_ARR_complex(t["ATMWFC."*string(a)][""])
                end
            end
            
            for k in 1:nk
                eigs[k, :] = parse_str_ARR_float(dae["K-POINT."*string(k)]["EIG"][""])
                
            end
            
            for k in 1:nk
                t = dao["K-POINT."*string(k)]
                overlaps[k, :,:] = reshape(parse_str_ARR_complex(t["OVERLAP.1"][""]), natwfc, natwfc)'
            end
        elseif nspin == 2
            for k in 1:nk
                t = dap["K-POINT."*string(k)]
                for a in 1:natwfc
                    proj[k, a,1, :] =  parse_str_ARR_complex(t["SPIN.1"]["ATMWFC."*string(a)][""])
                end
                for a in 1:natwfc
                    proj[k, a,2, :] =  parse_str_ARR_complex(t["SPIN.2"]["ATMWFC."*string(a)][""])
                end
            end
            for k in 1:nk
                eigs[k, :, 1] = parse_str_ARR_float(dae["K-POINT."*string(k)]["EIG.1"][""])
            end
            for k in 1:nk
                eigs[k, :, 2] = parse_str_ARR_float(dae["K-POINT."*string(k)]["EIG.2"][""])
            end
            
            for k in 1:nk
                t = dao["K-POINT."*string(k)]
                overlaps[k, :,:] = reshape(parse_str_ARR_complex(t["OVERLAP.1"][""]), natwfc, natwfc)'
            end
            
            
        end

        
    else #other format that QE switched to in order to break my code :<(

#        println(da["EIGENSTATES"])
        
#        println("asdf")
        if "" in keys(da["EIGENSTATES"])
            d_eigstates = da["EIGENSTATES"][""]


            n = length(d_eigstates)
            c = 0
            for n = 2:6:n

                c += 1


                
                weights[c] = parse(Float64, d_eigstates[n]["K-POINT"][:Weight])
                kpts[c,:] = parse_str_ARR_float(d_eigstates[n]["K-POINT"][""])

                eigs[c,:] = parse_str_ARR_float(d_eigstates[n+2]["E"])

                t =  d_eigstates[n+4]["PROJS"]["ATOMIC_WFC"]

                #            println("t ", t)
                ##            println("t[1]")
                #           println(parse_str_ARR_complex(t[1][""]))
                
                try
                    for a in 1:natwfc
                        proj[c,a,1,:] =  parse_str_ARR_complex(t[a][""])
                    end
                catch
                    for a in 1:natwfc
                        proj[c,a,1,:] =  parse_str_ARR_complex(t[""])
                    end
                end                
            end

            d_over = da["OVERLAPS"]["OVPS"]

            for k = 1:nk
                overlaps[k, :,:] =  reshape(parse_str_ARR_complex(d_over[k][""]), natwfc, natwfc)'
            end
            
        else
            d_eigstates = da["EIGENSTATES"]

            #------
            c = 1

            
            weights[c] = parse(Float64, d_eigstates["K-POINT"][:Weight])
            kpts[c,:] = parse_str_ARR_float(d_eigstates["K-POINT"][""])
            
            eigs[c,:] = parse_str_ARR_float(d_eigstates["E"])

            t =  d_eigstates["PROJS"]["ATOMIC_WFC"]

            try
                for a in 1:natwfc
                    proj[c,a,1,:] =  parse_str_ARR_complex(t[a][""])
                end
            catch
                for a in 1:natwfc
                    proj[c,a,1,:] =  parse_str_ARR_complex(t[""])
                end
            end                
            
            d_over = da["OVERLAPS"]["OVPS"]
            
            for k = 1:nk
                overlaps[k, :,:] =  reshape(parse_str_ARR_complex(d_over[""]), natwfc, natwfc)'
            end
            #0-------
        end        
    end
    
    if !(ismissing(B))
        
        kpts = kpts * inv(B)
    end
    
    


    #    println(overlaps)
    
    bs = DFToutMod.makebs(nelec, efermi, kpts, weights, [0,0,0], eigs, nspin=nspin)
    
    return make_proj(bs, proj, overlaps)
    
end

"""
    function make_proj(bs, proj, overlaps)

Constructor for proj_dat.
"""
function make_proj(bs, proj, overlaps)

    if bs.nks != size(proj)[1] || bs.nks != size(overlaps)[1]
        error("make_proj something wrong nks ", bs.nks," ",size(proj)[1]," ",size(overlaps)[1])
    end
    if bs.nbnd != size(proj)[4]
        error("make_proj something wrong nbnd ", bs.nbnd," ",size(proj)[4])
    end
    if size(proj)[2] != size(overlaps)[2] ||  size(proj)[2] !=	size(overlaps)[3]
        error("make_proj something wrong in natwfc ",  size(proj)[2]," ",size(overlaps)[2]," " ,size(overlaps)[3])
    end

    return proj_dat(bs, size(proj)[2], bs.nspin, proj, overlaps)
    
end


function shift_eigenvalues(d::dftout)

    println("sdf")
    wan, semicore, nwan, nsemi, wan_atom, atom_wan = tb_indexes(d)

    efermi_dft = d.bandstruct.efermi

#=    #decide semicore
    p2 = zeros(p.bs.nbnd)
    INDSEMI = zeros(Int64, p.bs.nks, nsemi)
    for k = 1:p.bs.nks
        p2[:] = real.(sum(p.proj[k,semicore,:] .* conj.(p.proj[k,semicore,:]) , dims=1))
        INDSEMI[k, :] = sortperm(p2, rev=true)[1:nsemi]
        if sum(p2) < nsemi - 0.5
            println("warning, identify semicore")
        end
#        if k < 10
#            println("p2 ", p2)
#            println(k, " indsemi " , INDSEMI[k, :])
#        end

    end
    println("INDSEMI ", INDSEMI)
=#
    NBND = d.bandstruct.nbnd - nsemi
    EIGS = zeros(d.bandstruct.nks, NBND)
#    PROJ = zeros(Complex{Float64}, p.bs.nks, nwan, NBND)

    # setup data
#    println("nsemi $nsemi")
    for k = 1:d.bandstruct.nks
        counter = 0
        for n = (1+nsemi):d.bandstruct.nbnd

            #if !(n in INDSEMI[k, :])
            if true
                counter += 1
                #println("k $k counter $counter k $k n  $n  EIGS  $(size(EIGS)) d $(size(d.bandstruct.eigs))  EIGS $(EIGS[k,counter]) = eigs $(d.bandstruct.eigs[k,n])")
                EIGS[k,counter] = d.bandstruct.eigs[k,n]
#                PROJ[k,:,counter] = p.proj[k,wan,n]
                
            end
        end
#        if k < 20
#            println("$k EIGS ", EIGS[k,1:5])
#        end
        

#        if k < 2
#            println("ALLEIGS k $k ", p.bs.eigs[k,1:4])
#            println("EIGS    k $k ", EIGS[k,1:4])#
#
#        end
    end

#    println("EIGS")
#    println(EIGS)
#    println("shifting eigenvalues to match dft atomization energy")
    ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(d.crys)
    #        band_en = band_energy(d.bandstruct.eigs[:,nsemi+1:end], d.bandstruct.kweights, nval)
    band_en = band_energy(EIGS, d.bandstruct.kweights, nval - d.tot_charge )

    #        println("d.crys")
    #        println(d.crys)
    
    etypes = types_energy(d.crys)
    
    
    etot_dft = d.energy
    e_smear = d.energy_smear
    
    atomization_energy = etot_dft - etotal_atoms - etypes  - e_smear
    
    
#    println("atomization_energy $atomization_energy")
    
    band_en = band_en 
    shift = (atomization_energy - band_en  )/ (nval - d.tot_charge)
#    println("shift $shift")
    
    EIGS = EIGS .+ shift
    
    return EIGS, shift

end


#function create_temp(p::proj_dat, d::dftout; energy_froz=missing, nfroz=0, shift_energy=true)
##
#
#    return
#end

"""
    function create_tb(p::proj_dat, d::dftout; energy_froz=missing, nfroz=0, shift_energy=true)

Does the main creation of TB hamiltonian from DFT projection data in k-space.

# Arguments
- `p::proj_dat` Projection data
- `d::dftout` DFT scf data
- `energy_froz=missing` Energy to start freezing eigenvalues
- `nfroz=0` number of frozen bands.
- `shift_energy=true` if `true` shift eigenvalues so band energy `==` total energy
"""
function create_tb(p::proj_dat, d::dftout; energy_froz=missing, nfroz=0, shift_energy=true)




    wan, semicore, nwan, nsemi, wan_atom, atom_wan = tb_indexes(d)

    efermi_dft = d.bandstruct.efermi

    #decide semicore
    p2 = zeros(p.bs.nbnd)
    INDSEMI = zeros(Int64, p.bs.nks, p.nspin, nsemi)
    for k = 1:p.bs.nks
        for spin = 1:p.nspin
            p2[:] = real.(sum(p.proj[k,semicore,spin, :] .* conj.(p.proj[k, semicore,spin, :]) , dims=1))
            INDSEMI[k,spin,  :] = sortperm(p2, rev=true)[1:nsemi]
            if sum(p2) < nsemi - 0.5
                println("warning, identify semicore")
            end
#            if k < 10
#                println("p2 ", p2)
#                println(k, " indsemi " , INDSEMI[k, :])
#            end
        end

    end
    #    println("INDSEMI ", INDSEMI)

    NBND = p.bs.nbnd - nsemi
    EIGS = zeros(p.bs.nks,  NBND,p.nspin)
    PROJ = zeros(Complex{Float64}, p.bs.nks, nwan, p.nspin, NBND)

    # setup data
    for k = 1:p.bs.nks
        for spin = 1:p.nspin
            counter = 0
            for n = 1:p.bs.nbnd
                if !(n in INDSEMI[k,spin, :])
                    counter += 1
                    EIGS[k,counter,spin] = p.bs.eigs[k,n, spin]
                    PROJ[k,:,spin, counter] = p.proj[k,wan,spin, n]
                    #if k == 1
                    #    println(PROJ[1,1,1,1], " $k $spin $counter  add p ", p.proj[k,wan,spin, n])
                    #end
                end
            end
        end
#        if k < 20
#            println("$k EIGS ", EIGS[k,1:5])
#        end
        

#        if k < 2
#            println("ALLEIGS k $k ", p.bs.eigs[k,1:4])
#            println("EIGS    k $k ", EIGS[k,1:4])#
#
#        end
    end

    #    println(p.nspin, " PROJ check ", PROJ[1,1,1,1])
    
    if shift_energy
        
        println("shifting eigenvalues to match dft atomization energy")
        ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(d.crys)
#        band_en = band_energy(d.bandstruct.eigs[:,nsemi+1:end], d.bandstruct.kweights, nval)
        band_en = band_energy(EIGS, d.bandstruct.kweights, nval)

        println("band_energy $band_en")
        println("EIGS  ", EIGS[1,:,1])
        if p.nspin == 2
            println("EIGS2 ", EIGS[1,:,2])
        end            
        
#        println("d.crys")
#        println(d.crys)

        etypes = types_energy(d.crys)

        
        etot_dft = d.energy
        e_smear = d.energy_smear

        atomization_energy = etot_dft - etotal_atoms - etypes  - e_smear

#        println("$etot_dft $etotal_atoms $etypes $e_smear")
#        println("atomization_energy $atomization_energy")

        band_en = band_en 
        shift = (atomization_energy - band_en  )/(nval - d.tot_charge)

#        println("shift $shift")
#        return 0

        EIGS = EIGS .+ shift
        

        if !(ismissing(energy_froz))
            println("energy_froz before $energy_froz")
#            energy_froz = energy_froz + shift
            energy_froz = efermi_dft + shift
            

        end
        println("energy_froz after $energy_froz")
        #println("sum shifted EIGS ", sum(EIGS))
    else
        println("no shift: match dft eigenvals")
    end
    

    P = zeros(Complex{Float64}, NBND, NBND)  #p.bs.nbnd-nsemi,p.bs.nbnd-nsemi)

    Pmat = zeros(Float64, p.bs.nks, NBND)
    Nmat = zeros(Float64, p.bs.nks, nwan)

    Btilde = zeros(Complex{Float64}, nwan, NBND)
    B = zeros(Complex{Float64}, nwan, NBND)

    
    ham_dft = zeros(Float64, NBND)

    ham_k = zeros(Complex{Float64},  p.natwfc-nsemi, p.natwfc-nsemi,p.bs.nks, p.nspin)


    htemp = zeros(Complex{Float64}, nwan, nwan)

    max_c = minimum(EIGS[:,end,:])
    min_c = max_c - 0.1

    
    
    min_c = max(maximum(EIGS[:,1:nwan,:])+.001, min_c) #ensure we don't cut off needed states


    #simple cutoff function
    function cutoff(num, min_c, max_c)
        if num < min_c
            return 1.0
        elseif num > max_c
            return 0.0
        else
            t = (num - min_c)/(max_c-min_c)
            return 1.0 - 10.0 * t^3 + 15.0 *  t^4  - 6.0 * t^5
        end
    end


    Palt_min = 1.0
    
    badk=0
    PROJECTABILITY = zeros(p.nspin, p.bs.nks, NBND)

    
    for spin = 1:p.nspin
        for k = 1:p.bs.nks
            #    for k = 1:nks
            
            #this avoids breaking symmetry by cutting off a symmetrically equivalent pair/triplet, etc because we include a finite number of bands

            energy_cutoff = EIGS[k,end,spin]
            max_ind = 0

            for i=max(1,NBND-5):NBND
                #            if p.bs.eigs[k,i] < energy_cutoff - 1e-3
                if EIGS[k,i,spin] < energy_cutoff - 1e-3
                    max_ind = i
                end
            end

            
            #        B[:,1:max_ind-nsemi] = p.proj[k,wan, 1+nsemi:max_ind]
            B[:,1:max_ind] = PROJ[k,:,spin, 1:max_ind]

            #smoothing
            for i in 1:max_ind
                #            B[:,i] *= cutoff(p.bs.eigs[k,i+nsemi], min_c, max_c)
                B[:,i] *= cutoff(EIGS[k,i,spin], min_c, max_c)
            end
            
            P[:,:] .= 0.0

            P[1:max_ind,1:max_ind] = B[:,1:max_ind]'*B[:,1:max_ind]
            PROJECTABILITY[spin,k,:] = real.(diag(P))

            println("PROJECTABILITY spin $spin k $k ", PROJECTABILITY[spin,k,:])
            

            Palt = B[:,1:max_ind]*B[:,1:max_ind]'

#            if k == 1
#                println("PALT $spin")
#                println(Palt)
#                println("stuff")
#                println(PROJ[k,:,spin, 1:max_ind])
#            end
            
            for i in 1:nwan


                if real(Palt[i,i]) < Palt_min
                    Palt_min = real(Palt[i,i])
                end

                #warn user if projection is bad, i.e. projectors aren't projecting onto anything. 
                if real(Palt[i,i]) < 0.80
                    println("Warning, difficult to project atomic wfc ", i," ",  real(Palt[i,i]), " k ", k, " tr(P) = " , tr(real(P)))
                end
                if real(Palt[i,i]) < 0.45
                    badk += 1
                    println("kinda bad k $k spin $spin :", badk)
                end
                if real(Palt[i,i]) < 0.35
                    badk += 1
                    println("very  bad k $k spin $spin :", badk)
                end
                
            end

            #eigenvalues of projection matrix are key to this method
            val, vect = eigen(Hermitian( 0.5 * (P[1:max_ind,1:max_ind] + P[1:max_ind,1:max_ind]')    ) ) # 

            good_proj = (max_ind - nwan + 1 ) : max_ind


            Btilde[:,1:max_ind] = vect[:,good_proj]'

            #approximate hamiltonian using highest eigenvalues of projection matrix.
            
            htemp[:,:] = Btilde[:,1:max_ind] * Diagonal(EIGS[k,1:max_ind, spin])  * Btilde[:,1:max_ind]'

            neweigs, newvect = eigen(Hermitian((htemp+htemp')/2.0))

            Pmat[k, :] = real(diag(P))
            Nmat[k, :] = neweigs
            

            ham_k[:,:,k,spin] = B[:,1:max_ind] * Btilde[:,1:max_ind]' * htemp * Btilde[:,1:max_ind] *  B[:,1:max_ind]'        
            
            
            ham_k[:,:,k,spin] = (ham_k[:,:,k,spin]  + ham_k[:,:,k,spin]')/2.0


        end
    end

    #decide if need to send serious warning to user.
    if badk / p.bs.nks > 0.15
        println("warning lots of bad kpoints ", badk, " of ", p.bs.nks)
        warn_badk=true
    else
        warn_badk=false
    end

    println("Min Palt: $Palt_min")

    println("EIG TEST ", eigvals(ham_k[:,:,1,1]))
    

    #here we shift the eigenvalues around to match DFT eigenvalues below a cutoff.
    #this requires identifying which bands are supposed to match which eigenvlues, which
    #can be tricky. I use a projection heuriestic, but it can fail at band crossings that mix bands.
    if !(ismissing(energy_froz))
        energy_froz2 = energy_froz+2.0
        println("energy_froz: $energy_froz , $energy_froz2")
        for spin = 1:p.nspin
            for k = 1:p.bs.nks
                val_tbt, vect = eigen(Hermitian(ham_k[:,:,k, spin] ))
                val_tb = real(val_tbt)
                val_tb_new = deepcopy(val_tb)
                #            val_pw = p.bs.eigs[k,nsemi+1:nsemi+nwan]

                nxxx = size(PROJECTABILITY)[3]
                val_pw = EIGS[k,1:nxxx, spin]
            
                order = Dict()
                
                score_mat = zeros(nxxx, nwan)
                

                for n2 = 1:nwan
                    cd_vect = real(vect[:,n2].*conj(vect[:,n2]))
                    
                    dist_min = 1000.0
                    nmin = 0
                    
                    for n1 = 1:nxxx
                        
                        #                    t = p.proj[k,wan, nsemi + n1] 
                        t = PROJ[k,:,spin, n1]
                        x = PROJECTABILITY[spin,k,n1]
                        
                        cd_dft = real(t .* conj(t))
                        
#                        dist_en = (val_tb[n2] - val_pw[n1]).^2 * 1.0
                        dist_cd = sum((cd_dft - cd_vect).^2)
                        #                        dist = dist_en + dist_cd + 0.1*(n1 - n2)^2
                        dist = (val_tb[n2] - val_pw[n1]).^2 + dist_cd + 0.05*(n1 - n2)^2 + 1 / (x + 1e-3)

                        if k == 1
                            println("order k $k n1 $n1 n2 $n2 dist $dist val_tb $(val_tb[n2]) val_pw $(val_pw[n1])  x $(PROJECTABILITY[spin,k,n1])  cd_dft $(cd_dft)")
                        end
                        score_mat[n1,n2] = dist 
                        #                    if dist < dist_min
                        #                        dist_min = dist
                        #                        nmin = n1
                        #                    end
                    end
                end
                for n = 1:nwan
                    s = sortperm(score_mat[:,n], rev=false)
                    for ss in s
                        if !(ss in keys(order))
                            order[n] = ss
                            break
                        end
                    end
                end

#acutally do the change
                for n in 1:nwan
                    if val_pw[order[n]] < energy_froz
                        val_tb_new[n] = val_pw[order[n]]
                        #elseif val_pw[order[n]] < energy_froz2
                    else
                        #x=cutoff(val_tb[n], energy_froz, energy_froz2)
                        x = PROJECTABILITY[spin,k,order[n]]
                        #      if x > 0.1
                        #x = x^2
                        val_tb_new[n] = val_pw[order[n]] * x + val_tb[n]*(1.0-x)
                        println("n $n k $k x $x xsqrt $(sqrt(x)) val_pw $(val_pw[order[n]]) val_tb $(val_tb[n])")
                  #      end
                    end
                    
                end

#resym
                for n = 1:nwan-1
                    c= [n]
                    for n2 = n+1:nwan
                        if abs(val_tb_new[n] - val_tb_new[n2]) < 1e-4
                            push!(c, n2)
                        end
                    end
                    if length(c) > 1
                        t = sum(val_tb_new[c]) / length(c)
                        val_tb_new[c] .= t
                    end
                end
                
                ham_k[:,:,k, spin] = vect*Diagonal(val_tb_new)*vect'
                ham_k[:,:,k, spin] = (ham_k[:,:,k, spin]  + ham_k[:,:,k, spin]')/2.0
                val_tb_new, vect = eigen(Hermitian(ham_k[:,:,k, spin] ))
            
            end
        end
    end
    
#alternate eigenvalue fix method.
    if nfroz >= 1 && (ismissing(energy_froz))
        println("nfroz: ", nfroz)
        for spin = 1:p.nspin
            for k = 1:p.bs.nks
                val, vect = eigen(Hermitian(ham_k[:,:,k,spin] ))
                val[1:nfroz] = EIGS[k,1:nfroz,spin]
                if (abs(EIGS[k, nfroz+1,spin] - EIGS[k, nfroz,spin])< 1e-5) && nwan >= nfroz+1
                    val[nfroz+1] = EIGS[k,nfroz+1,spin]
                    if (abs(EIGS[k, nfroz+2,spin] - EIGS[k, nfroz,spin])< 1e-5) && nwan >= nfroz+2
                        val[nfroz+2] = EIGS[k,nfroz+2,spin]
                    end
                end
                
                


                ham_k[:,:,k,spin] = vect*Diagonal(val)*vect'
                ham_k[:,:,k,spin] = (ham_k[:,:,k,spin]  + ham_k[:,:,k,spin]')/2.0

            end
        end
    end

   
    #reshift to match dft total energy
    
    VAL = zeros(Float64, p.bs.nks, nwan, p.nspin)
    VECT = zeros(Complex{Float64}, p.bs.nks, p.nspin, nwan, nwan)
    for spin = 1:p.nspin
        for k = 1:p.bs.nks
            val, vect = eigen(Hermitian(ham_k[:,:,k, spin]))
            VAL[k,:, spin] = val
            VECT[k,spin,:,:] = vect
        end
    end
        
    if shift_energy
        println("reshift")
#        println("size VAL ", size(VAL), " size kweights ", size(d.bandstruct.kweights))

        println("VAL  ", VAL[1,:,1])
        if p.nspin == 2
            println("VAL  ", VAL[1,:,2])
        end
            
        
        band_en = band_energy(VAL, d.bandstruct.kweights, nval - d.tot_charge)

        println("band_en_old " , band_en)

        shift = (atomization_energy - band_en)/(nval - d.tot_charge)
        VAL = VAL .+ shift
        for spin = 1:p.nspin
            for k = 1:p.bs.nks
                
                val = VAL[k,:, spin]
                vect = VECT[k,spin,:,:]
                ham_k[:,:,k,spin] = vect*Diagonal(val)*vect'
                
            end
        end
        println("done reshift")

        for spin = 1:p.nspin
            for k = 1:p.bs.nks
                val, vect = eigen(Hermitian((ham_k[:,:,k, spin] + ham_k[:,:,k,spin]')/2.0))
                VAL[k,:,spin] = val
                VECT[k,spin,:,:] = vect
            end
        end
            
        band_en = band_energy(VAL, d.bandstruct.kweights, nval - d.tot_charge)
        println("band_en_new " , band_en, " " , band_en+etypes)


    end

    println("SIZE ham_k ", size(ham_k))
    
    return ham_k, EIGS, Pmat, Nmat, VAL, warn_badk




end
    

"""
    function get_ham_r_slow(p::proj_dat, d::dftout, grid, ham_k::Array{Complex{Float64}, 3}; nonorth=true, K=missing)

Slow method for doing fourier transform. Uses standard ft formula, not fft.

Use fast version instead.
"""
function get_ham_r_slow(p::proj_dat, d::dftout, grid, ham_k::Array{Complex{Float64}, 3}; nonorth=true, K=missing)
"""
DOESN'T USE FFTW. MOSTLY FOR DEBUGGING / missing fft libraries

Do Fourier transform k->R. grid is the size of the k-space grid. 
Optionally you can include the k space grid explictly, but normally 
K is already inside p. K must match ham_k
"""

    
    wan, semicore, nwan, nsemi, wan_atom, atom_wan = tb_indexes(d)
    
#    grid = d.bandstruct.kgrid

    println("grid ", grid)
    println("ham_k size " , size(ham_k))
    println("nonorth: ", nonorth)
    nwan = size(ham_k)[2]

    if ismissing(K)
        K = p.bs.kpts
        nks = size(K)[1]
    else
        nks = size(K)[1]
    end
    
    if nonorth 
        #nonorthogonal!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ham_kS = zeros(Complex{Float64},  nwan, nwan,nks)
        S = zeros(Complex{Float64}, nwan, nwan)
        Sk = zeros(Complex{Float64}, nwan, nwan,nks)    
        Sk2 = zeros(Complex{Float64},  nwan, nwan,nks)
        for k in 1:nks
            S[:,:] = p.overlaps[k, wan,wan]
            S = (S+S')/2.0
            Sk[:,:,k] =  S
        end
        for i=1:nwan #normalize overlaps
            sumi = sum(Sk[i,i,:])
            Sk[i,i,:] = Sk[i,i,:] / real(sumi) * Float64(nks)
        end
        for k in 1:nks
            S[:,:] = Sk[:,:,k]
            S = (S + S')/2.0
            Sk2[:,:,k] = sqrt(S)
            
            ham_kS[:,:,k] = Sk2[:,:,k]*ham_k[:,:,k]*Sk2[:,:,k]
            ham_kS[:,:,k] = (ham_kS[:,:,k]  + ham_kS[:,:,k]')/2.0
        
        end
        #end nonorthogonal!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end
        
    
    r=zeros(Float64,3)
    twopi_i = 1.0im*2.0*pi

    grid2 = [0,0,0]
    for i = 1:3
        if grid[i]%2 == 0
            grid2[i] = Int(grid[i]/2)
        else        
            grid2[i] = Int((grid[i]-1)/2)
        end
    end
    
    r_dict = Dict()


    ngrid2 = (grid2[1]*2+1)*(grid2[2]*2+1)*(grid2[3]*2+1)
    println("new real space grid: " , -grid2, " to " , grid2, " " , ngrid2)

    ind_arr = zeros(ngrid2,3)
    ham_r = zeros(Complex{Float64},  nwan, nwan,ngrid2)

    if nonorth
        S_r = zeros(Complex{Float64},  nwan, nwan,ngrid2)
    end
    
#    return 0

    sym = zero(Float64)
    exp_val = zero(Complex{Float64})
    c=0

    sym_mat = zeros(Float64, ngrid2)

    for r1 = -grid2[1]:(grid2[1])

        if (r1 + grid[1] == grid2[1]) ||  (r1 - grid[1] == -grid2[1])
            sym1 = 2.0
        else
            sym1 = 1.0
        end
            
        for r2 = -grid2[2]:(grid2[2])
            if (r2 + grid[2] == grid2[2]) ||  (r2 - grid[2] == -grid2[2])
                sym2 = 2.0
            else
                sym2 = 1.0
            end

            
            for r3 = -grid2[3]:(grid2[3])

                if (r3 + grid[3] == grid2[3]) ||  (r3 - grid[3] == -grid2[3])
                    sym3 = 2.0
                else
                    sym3 = 1.0
                end


#                println("rs $r1 $r2 $r3 $sym1 $sym2 $sym3 $sym")
                
                c+=1 
                sym_mat[c] = sym1 * sym2 * sym3
                r[:] = [r1 r2 r3]
                ind_arr[c,:] = r[:]
                r_dict[ind_arr[c,:]] = c
            end
        end
    end


    if nonorth
        for k in 1:nks
            for c = 1:ngrid2            
                exp_val = exp(twopi_i * (ind_arr[c,:]'*K[k,:])) / sym_mat[c]
                ham_r[:,:, c ] .+= ham_kS[:,:,k] .* exp_val
                S_r[:,:, c] .+=  Sk[:,:,k] .* exp_val

            end
        end
            
    else
        for k in 1:nks
            for c = 1:ngrid2            
                exp_val = exp(twopi_i * (ind_arr[c,:]'*K[k,:])) / sym_mat[c]
                ham_r[ :,:, c] .+= ham_k[:,:,k] .* exp_val
            end
        end
    end



    #turn ham_r into tight-binding object
    if nonorth
        tb = make_tb(ham_r/ Float64(nks), ind_arr, r_dict, S_r / Float64(nks))
    else
        tb = make_tb(ham_r/ Float64(nks), ind_arr, r_dict)
    end

#    renormalize_tb(d, tb)

    nelec=d.bandstruct.nelec - nsemi * 2.0

    ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(d.crys)
    dftenergy=d.energy - etotal_atoms

    model = make_tb_crys(tb, d.crys, nelec, dftenergy)

    
    return model
    
end    
    
"""
    function prepare_ham_k(p::proj_dat, d::dftout, grid, ham_k::Array{Complex{Float64}, 4}; nonorth=true, K=missing, localized_factor = 0.05, screening=1.0)

Constructs the actual k-space hamiltonain, which involves lots of putting matricies in the correct form and de-orthogonalizing if desired.

#arguments
- `p::proj_dat` projection data
- `d::dftout` dft data
- `grid` kpoint grid spacing 
- `ham_k::Array{Complex{Float64}, 3}` Hamiltonian
- `nonorth=true` make a non-orthogonal TB 
- `K=missing` k-point array, usually get from `p`
- `localized_factor = 0.15` increase localization of overlaps.
- `screening=1.0` mulitply U by this factor. usually not used.
"""
function prepare_ham_k(p::proj_dat, d::dftout, grid, ham_k::Array{Complex{Float64}, 4}; nonorth=true, K=missing, localized_factor = 0.15, screening=1.0)

    wan, semicore, nwan, nsemi, wan_atom, atom_wan = tb_indexes(d)
    
    println("grid ", grid)
    println("ham_k size " , size(ham_k))
    println("nonorth: ", nonorth)
    nwan = size(ham_k)[2]
    
    if ismissing(K)
        K = p.bs.kpts
        nks = size(K)[1]
    else
        nks = size(K)[1]
    end
    
    #fix mysterious sign error K-space. QE overlaps have weird signs in x and y direction
    # I don't know why, you have to ask them.
    
    ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(d.crys)

    nval = nval - d.tot_charge
    
    OVERLAPS = deepcopy(p.overlaps)
    println("wan ", wan)
    println("so ", size(OVERLAPS))
    OVERLAPS = OVERLAPS[:,wan, wan]
    println("so2 ", size(OVERLAPS))
    for n = 1:nwan
        a1,t1,orb = ind2orb[n]
        if orb == :px || orb == :py
            println("fix mysterious sign issue K-space $n, $a1, $t1, $orb")
            ham_k[n,:,:,:] = -1.0 * ham_k[n,:,:,:] #spin
            ham_k[:,n,:,:] = -1.0 * ham_k[:,n,:,:]

            if nonorth
                OVERLAPS[:,n,:] = -1.0*OVERLAPS[:,n,:]
                OVERLAPS[:,:,n] = -1.0*OVERLAPS[:,:,n]
                
            end
        elseif orb == :dxz || orb == :dyz
            println("fix mysterious sign issue K-space D ORBITALS $n, $a1, $t1, $orb")
            ham_k[n,:,:,:] = -1.0 * ham_k[n,:,:,:] #spin
            ham_k[:,n,:,:] = -1.0 * ham_k[:,n,:,:]

            if nonorth
                OVERLAPS[:,n,:] = -1.0*OVERLAPS[:,n,:]
                OVERLAPS[:,:,n] = -1.0*OVERLAPS[:,:,n]
        
            end

        end
    end




    if nonorth 
        #nonorthogonal!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ham_kS = zeros(Complex{Float64},  nwan, nwan,nks, p.nspin)
        S = zeros(Complex{Float64}, nwan, nwan)
        Sk = zeros(Complex{Float64}, nwan, nwan,nks)    
        Sk2 = zeros(Complex{Float64},  nwan, nwan,nks)
            

        
        for k in 1:nks
#            S[:,:] = p.overlaps[k, wan,wan]
            S[:,:] = OVERLAPS[k, :,:]

            S = (S+S')/2.0           #+ I(nwan) * localized_factor
            Sk[:,:,k] =  S
        end

#        println("fourteen")
#        println(Sk[:,:,14])

        wtot = sum(p.bs.kweights)
        SUMI = zeros(nwan)
        for i=1:nwan #normalize overlaps
            SUMI[i] = sum(Sk[i,i,:] .* p.bs.kweights) / wtot
        end
        println("old SUMI ", SUMI)
        SUMI = symm_by_orbitals(d.crys, SUMI) #this handles the fact the k-grid can break orbital symmetry
        println("new SUMI ", SUMI)
        for i=1:nwan #normalize overlaps
            Sk[i,i,:] = Sk[i,i,:] / SUMI[i]
        end


        for k in 1:nks
            S = Sk[:,:,k]
            S = (S + S')/2.0
            Sk[:,:,k] = S
        end            

        if localized_factor > 1e-5
            for k = 1:nks
                Sk[:,:,k] = Sk[:,:,k] * (1.0 - localized_factor) + I(nwan) * localized_factor 
            end
        end

        #does the de-orthogonalization
        for k in 1:nks
            S[:,:] = Sk[:,:,k]
            S = (S+S')/2.0
            Sk2[:,:,k] = sqrt(S)

            for spin = 1:p.nspin
                ham_kS[:,:,k, spin] = Sk2[:,:,k]*ham_k[:,:,k, spin]*Sk2[:,:,k]
                ham_kS[:,:,k, spin] = (ham_kS[:,:,k, spin]  + ham_kS[:,:,k, spin]')/2.0
            end
        end

        tbk = make_tb_k(ham_kS, K, d.bandstruct.kweights, Sk, grid=grid, nonorth=true)
        tbck = make_tb_crys_kspace(tbk, d.crys, nval, d.energy - etotal_atoms, scf=false, screening=screening)

#        println("nspin k ", tbck.nspin)
        
        return tbck

        #end nonorthogonal!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
        
        Sk = zeros(Complex{Float64}, nwan, nwan,nks)    
        for i = 1:nks
            Sk[i,:,:] = I(nwan)
        end

        tbk = make_tb_k(ham_k, K, d.bandstruct.kweights, Sk, grid=grid, nonorth=false)
        tbck = make_tb_crys_kspace(tbk, d.crys, nval, d.energy - etotal_atoms, scf=false)
        
        return tbck

    end

    ###########


end



"""
    function get_ham_r(tbck::tb_crys_kspace)

Do fft to get the real-space ham. Requires tbck to be on a regular grid centered at Gamma with no symmetry.
"""
function get_ham_r(tbck::tb_crys_kspace, tot_charge)
"""
Do Fourier transform k->R. 
"""

    nonorth = tbck.tb.nonorth
    ham_k = tbck.tb.Hk
    Sk = tbck.tb.Sk
    K = tbck.tb.K
    grid = tbck.tb.grid
    nwan = tbck.tb.nwan
    nspin = tbck.tb.nspin
    
    println("$nspin get_ham_r size ham_k ", size(ham_k))

    if nonorth
        ham_r, S_r, r_dict, ind_arr = myfft(tbck.crys, nonorth, grid, K,ham_k, Sk)
    else
        ham_r, r_dict, ind_arr =      myfft(tbck.crys, nonorth, grid, K,ham_k, missing)
    end             

#    println("size ham_r ", size(ham_r))
    
    ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(tbck.crys)

    #turn ham_r into tight-binding object
    if nonorth
        tb = make_tb(ham_r, ind_arr, r_dict, S_r )        
    else
        S_r = zeros(size(ham_r))

        println(size(S_r))

        r = r_dict[[0,0,0]]
        for i =  1:nwan
            S_r[i,i,r] = 1.0  #trivial overlaps
        end
        tb = make_tb(ham_r, ind_arr, r_dict, S_r)        
    end

#    renormalize_tb(d, tb)
#    nelec=.bandstruct.nelec ##- nsemi * 2.0
    nelec = nval - tot_charge

    dftenergy=tbck.dftenergy
    ind2orb, orb2ind, etotal_atoms, nval =  orbital_index(tbck.crys)


    println("make_tb_crys tb nspin ", tb.nspin)
    model = make_tb_crys(tb, tbck.crys, nelec, dftenergy )
    
    return model
    
end    
    








##function run_proj(inputstr, outputstr, nprocs=1, directory="./")
##"""
##run command
##"""
##    directory=rstrip(directory, '/')
##    
##    #get commandline
##    c_dict = make_commands(nprocs)
##    qe = c_dict["qe"]
##    
##
##    command = `$qe $directory/$inputstr`
##    println("actual command")
##    println(command)
##    s=""
##    try
##        println("Running DFT")
##
##        s = read(command, String)
##        
##        println("Writing output")
##        f = open(directory*"/"*outputstr, "w")
##        write(f, s)
##        close(f)
##        
##        println("Ran DFT")
##        return 0
##    catch
##        println("Failed to run DFT")
##        println(s)
##        return -1
##    end
##        
##end
##


end

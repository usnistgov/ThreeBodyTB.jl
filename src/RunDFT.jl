###include("Crystal.jl")
###include("DFToutMod.jl")
#using XMLDict



####################### QE specific 
"""
    module QE

Module for running Quantum Espresso. Generic DFT version below.
"""
module QE
"""
Scripts to run Quantum Espresso
"""
#export load_xml
using LinearAlgebra

#using ExXML
using XMLDict
using GZip
using ..DFToutMod
using ..CrystalMod:crystal
using ..Utility:arr2str
using ..Utility:str_w_spaces
using ..Utility:parse_str_ARR_float
#using ..TB:get_grid
using ..CrystalMod:get_grid
using ..Atomdata:atoms

include("Commands.jl")

using ..ThreeBodyTB:TEMPLATES
using ..ThreeBodyTB:PSEUDOS


"""
    function run_pwscf(inputstr, outputstr, nprocs=1, directory="./", use_backup=false)

Run the pw.x code from QE on `inputstr`
"""
function run_pwscf(inputstr, outputstr, nprocs=1, directory="./", use_backup=false)
"""
run command
"""
    directory=rstrip(directory, '/')
    
    #get commandline
    c_dict = make_commands(nprocs)

    if use_backup
        qe = c_dict["qe_backup"]
    else
        qe = c_dict["qe"]
    end    

    command = `$qe $directory/$inputstr`
    println("actual command")
    println(command)
    s=""
    try
        println("Running DFT")

        s = read(command, String)
        
        println("Writing output")
        f = open(directory*"/"*outputstr, "w")
        write(f, s)
        close(f)
        
        println("Ran DFT")
        return 0
    catch err
        println(err)
        println("Failed to run DFT")
        println(s)
        return -1
    end
        
end


"""
    function runSCF(crys::crystal, inputstr=missing, prefix=missing, tmpdir="./", directory="./", functional="PBESOL", wannier=0, nprocs=1, skip=false, calculation="scf", dofree="all", tot_charge = 0.0, smearing = 0.01, magnetic=false, cleanup=false, use_backup=false)

Workflow for doing SCF DFT calculation on `crys`

Return `dftout`
"""
function runSCF(crys::crystal, inputstr=missing, prefix="qe", tmpdir="./", directory="./", functional="PBESOL", wannier=0, nprocs=1, skip=false, calculation="scf", dofree="all", tot_charge = 0.0, smearing = 0.01, magnetic=false, cleanup=false, use_backup=false, grid=missing, klines=missing, nstep=30, startingpot=missing, hybrid=false)
"""
Run SCF calculation using QE
"""
    println("runSCF hybrid ", hybrid)
    qeout = []
#    println("runSCF 1")
    if skip == true
        try
            savedir="$tmpdir/$prefix.save"
            qeout = loadXML(savedir)
            println("skip: we loaded SCF data instead of rerunning")
            return qeout
        catch
            println("skipping failed, continue scf calculation")
        end
    end

    if ismissing(inputstr)
        if calculation == "nscf" || calculation == "nscf-sym"
            inputstr="qe.nscf.in"
        elseif calculation == "relax"
            inputstr="qe.relax.in"
        else
            inputstr="qe.in"
        end
    end
    outputstr = inputstr*".out"

    if !isdir(directory)
        println("warning, making directory path $directory")
        mkpath(directory)
    end

    savedir="$tmpdir/$prefix.save"


#    println("runSCF 2")
    
    tmpdir, prefix, inputfile =  makeSCF(crys, directory=directory, prefix=prefix, tmpdir=tmpdir, functional=functional, wannier=wannier, calculation=calculation, dofree=dofree, tot_charge=tot_charge, smearing=smearing, magnetic=magnetic, grid=grid, klines=klines, nstep=nstep, startingpot=startingpot, hybrid=hybrid)
    
    f = open(directory*"/"*inputstr, "w")
    write(f, inputfile)
    close(f)

#    println("runSCF 3")

    ret = run_pwscf(inputstr, outputstr, nprocs, directory, use_backup)    
#    println("DFT complete with error")
#    println("runSCF 4")

    updated = false
    if calculation == "relax" || calculation == "vc-relax"
        try
            qeout = loadXML(savedir)
            crys = qeout.crys
            println("updated crys")
            updated = true
        catch
            println("couldn't open results")
        end
    end

    if ret != 0 && updated == false
        println("failed DFT, trying with different mixing")
        #tmpdir, prefix, inputfile =  makeSCF(crys, directory, prefix, tmpdir, functional, wannier, calculation, dofree, tot_charge, smearing, magnetic, mixing="TF", grid=grid, klines=klines)
        tmpdir, prefix, inputfile =  makeSCF(crys, directory=directory, prefix=prefix, tmpdir=tmpdir, functional=functional, wannier=wannier, calculation=calculation, dofree=dofree, tot_charge=tot_charge, smearing=smearing, magnetic=magnetic, mixing="TF",  grid=grid, klines=klines, nstep=nstep, hybrid=hybrid)
        

        f = open(directory*"/"*inputstr, "w")
        write(f, inputfile)
        close(f)

        ret = run_pwscf(inputstr, outputstr, nprocs, directory, use_backup)    

        if ret != 0
            println("warning, run_pwscf threw an error again: $ret")
            error("runSCF")
        end
       
    end
    
#    println("runSCF 5")


#    savedir="$tmpdir/$prefix.save"

    println("directory: ", directory)
    println("outputfile: ", outputstr)
    println("savedir: ", savedir)

    try
        qeout = loadXML(savedir)
    catch
        error("runSCF couldn't open output final try")
    end

    if cleanup
        println("call cleanup")
        doclean(savedir)
    end

#    println("runSCF 6")


    return qeout
    
end

"""
    function doclean(d)

Clean up wavefunctions in directory `d`
"""
function doclean(d)
    tot = 0

    files=read(`ls $d/`, String)
    for f in split(files, "\n")
        if occursin("wfc", f)
            if tot <= 2
                println("removing $d/$f")
            elseif tot==3
                println("removing $d/$f")
                println("etc...")
                println()
            end
            cmd=`\rm $d/$f`
            s = read(cmd, String)
            tot += 1
        end
    end
    println("clean-up done: $tot wfc files deleted")
end

"""
    function makeSCF(crys::crystal, directory="./", prefix=missing, tmpdir=missing, functional="PBESOL", wannier=0, calculation="scf", dofree="all", tot_charge = 0.0, smearing = 0.01, magnetic=false; mixing="local-TF")

Make QE inputfile for SCF DFT calculation.
"""

function makeSCF(crys::crystal; directory="./", prefix=missing, tmpdir=missing, functional="PBESOL", wannier=0, calculation="scf", dofree="all", tot_charge = 0.0, smearing = 0.01, magnetic=false, mixing="local-TF", grid=missing, klines = missing, nstep=30, startingpot=missing, hybrid=false)
"""
Make inputfile for SCF calculation
"""
    println("makeSCF hybrid ", hybrid)

    c_dict = make_commands(1)

    template_file=open("$TEMPLATES/template_qe.in")
    temp = read(template_file, String)
    close(template_file)

    temp = replace(temp, "PSEUDODIR" => PSEUDOS)
    
    temp = replace(temp, "JULIANAT" => crys.nat)


    c = calculation
    if calculation == "nscf-sym"
        c = "nscf"
    end
    temp = replace(temp, "SCF" => c)

    temp = replace(temp, "MIXING" => mixing)

    if !ismissing(startingpot) && startingpot == "file"
        temp = replace(temp, "STARTATOMIC" => startingpot)
    else
        temp = replace(temp, "STARTATOMIC" => "atomic")
    end

    if calculation == "nscf" 
        temp = replace(temp, "nosym = false" => "nosym = true")
        temp = replace(temp, "noinv = false" => "noinv = true")
    end
        
    
    settypes = []
    for t in crys.types
        if !(t in settypes)
            push!(settypes, t)
        end
    end

#    println("settypes", settypes)
    
    ntypes = length(settypes)
    temp = replace(temp, "JULIANTYPE" => ntypes)

    
    
    brav="ibrav = 0"

#    if abs(tot_charge) > 1e-5
#        a = crys.A[1,1]
#        brav="ibrav = 1, celldm(1) = $a "
#    end

    temp = replace(temp, "JULIABRAV" => brav)

    function convert_to_qe_type(t)
        if occursin("_",t)
            return t[1:end-2]*"D"
        else
            return t
        end
        return t
    end
    
    t=""
    for i = 1:ntypes
        if functional=="PBESOL"
            psp=lowercase(settypes[i])*".pbesol.UPF"
        elseif functional=="LDA"
            psp=lowercase(settypes[i])*".lda.UPF"
        else
            psp=lowercase(settypes[i])*".pbe.UPF"
        end
        
        mass = atoms[settypes[i]].mass
        t = t * str_w_spaces([convert_to_qe_type(settypes[i]), mass, psp, "\n"])
    end
    temp = replace(temp, "JULIAPSP\n" => t) ###strip(t,"\n"))
    

#    if abs(tot_charge) < 1e-5

    temp = replace(temp, "JULIACELL\n" => arr2str(crys.A))#strip(arr2str(crys.A), "\n"))

#    else
 
    #temp = replace(temp, "JULIACELL\n" => "")#strip(arr2str(crys.A), "\n"))        
    #temp = replace(temp, "CELL_PARAMETERS\n"  => "")
    
    #end

    t=""

    for i = 1:crys.nat
        t=t*convert_to_qe_type(crys.types[i])*" "
        t=t*arr2str(crys.coords[i,:])
        if i != crys.nat
            t=t*"\n"
        end

    end
    temp = replace(temp, "JULIACOORDS" => t)

    
    temp = replace(temp, "JULIASMEARING" => smearing)


#=    B = transpose(inv(crys.A))
    kden = 55.0 
    b1 = norm(B[1,:])
    b2 = norm(B[2,:])
    b3 = norm(B[3,:])

    k1 = convert(Int, round(kden * b1))
    k2 = convert(Int, round(kden * b2))
    k3 = convert(Int, round(kden * b3))

    if k1%2 == 1
        k1 = k1 + 1
    end
    if k2%2 == 1
        k2 = k2 + 1
    end
    if k3%2 == 1
        k3 = k3 + 1
    end

    
    kpoints = [k1, k2, k3]
=#
    qpoints = [1,1,1]
    if ismissing(klines)
        if ismissing(grid)
            kpoints = get_grid(crys)
        else
            if length(grid) == 3
                kpoints =  grid
            else
                kpoints = get_grid(crys, grid)
            end
        end
        
        if hybrid
            qpoints = [1,1,1] #use fewer qpoints than kpoints
            for i = 1:3
                if kpoints[i] == 1 || kpoints[i] == 2
                    qpoints[i] = 1
                elseif kpoints[i] == 3
                    qpoints[i] = 2
                else
                    qpoints[i] = Int64(ceil(kpoints[i] / 3))
                end
            end
            #qpoints must divide kpoints evenly
            for i = 1:3
                for j = 1:5
                    if abs(mod(kpoints[i]/qpoints[i], 1.0) ) > 1e-5
                        kpoints[i] += 1
                    else
                        break
                    end
                end
            end
        end

        
#    if calculation == "nscf"
        #        kpoints = get_grid(c, kden=50.0)
        #    end
        temp = replace(temp, "JULIAKTYPE" => "K_POINTS automatic")
        temp = replace(temp, "JULIAKPOINTS" => arr2str(kpoints)*" 0 0 0 ")

    else
        
#        println("KLINES xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
#        println(klines)

        if calculation != "nscf" 
            println("warning, you are specifying the klines without using nscf, i hope you know what you are doing")
        end

        nk = size(klines)[1]
        if size(klines)[2] == 3
            klines2 = ones(nk,4)
            klines2[:,1:3] = klines
            klines = deepcopy(klines2)
        end

        kstr = "$nk\n"
        kstr *= arr2str(klines)

#        println("before ", temp)

        temp = replace(temp, "JULIAKTYPE" => "K_POINTS crystal")
        temp = replace(temp, "JULIAKPOINTS" => kstr)

        if hybrid
            println("warning, klines and hybrid don't currently work")
        end
        
#        println("after", temp)

    end
        
    if ismissing(tmpdir)
        tmpdir=string(rand(1)[1])
    end
    
    if ismissing(prefix)
        prefix=string(rand(1)[1])    
    end
    
    temp = replace(temp, "PREFIX" => prefix)
    temp = replace(temp, "TMPDIR" => tmpdir)    

    other=""
    println("wannier $wannier ")
    if (typeof(wannier) == Bool && wannier==true) || wannier == 2
        nbandsemi = 0
        nbandval = 0        
        
        for t in crys.types
            nbandsemi += atoms[t].nsemicore
            nbandval += atoms[t].nwan
        end
        nbandtot = 4+1+convert(Int64, round(nbandsemi/2.0 + nbandval * 3.5 / 2.0))   #no spin yet
        println("nbandtot ", [4,1,convert(Int64, round(nbandsemi/2.0)), convert(Int64, nbandval * 3.5 / 2.0)])
        if "La" in crys.types
            nbandtot += 7*crys.nat
        end

        println("nbandtot $nbandtot")
        
        other *="nbnd = "*string(nbandtot)*"\n"
    elseif wannier == 1 #fewer bands can be more stable
        nbandsemi = 0
        nbandval = 0        
        
        for t in crys.types
            nbandsemi += atoms[t].nsemicore
            nbandval += atoms[t].nwan
        end
        nbandtot = 0+convert(Int64, round(nbandsemi/2.0 + nbandval * 2.0 / 2.0))   #no spin yet
        other *="nbnd = "*string(nbandtot)*"\n"
    elseif wannier == -1
        nbandsemi = 0
        nbandval = 0        
        
        for t in crys.types
            nbandsemi += atoms[t].nsemicore
            nbandval += atoms[t].nwan
        end
        nbandtot = 0+convert(Int64, round(nbandsemi/2.0 + nbandval * 1.5 / 2.0))   #no spin yet
        other *="nbnd = "*string(nbandtot)*"\n"
    end

    if abs(tot_charge) > 1e-5
        #other *= "   tot_charge = $tot_charge , assume_isolated = 'mp' \n"
        other *= "   tot_charge = $tot_charge \n"
    end

    if magnetic
        other *= "   nspin = 2 \n"
        for c in 1:ntypes
            other *= "  starting_magnetization( $c ) = 0.7 \n"
        end
    end

    if hybrid
        other *= "input_dft='gaupbesol' \n" #requires modified QE, ask kfg for help.
        println("warning - hybrid functional pbesol requires modified QE")
        
        other *= "nqx1 = $(qpoints[1]) , nqx2 = $(qpoints[1]) , nqx3 = $(qpoints[1]) ,\n"
        other *= "x_gamma_extrapolation = .false.\n"
        other *= "exxdiv_treatment = 'none' \n"

        temp = replace(temp, "conv_thr = 1d-9" => "conv_thr = 1d-4")
    end
    
    temp = replace(temp, "JULIAOTHER" => other)

    st_free = "cell_dofree = 'all'"
    st_free2 = replace(st_free, "all"=> dofree)
    
    temp = replace(temp, "JULIADOFREE" => st_free2)
    
    if calculation == "vc-relax"
        temp = replace(temp, "JULIACELLFACTOR" => "cell_factor = 3.0")
    else
        temp = replace(temp, "JULIACELLFACTOR" => "cell_factor = 1.0")
    end

    temp = replace(temp, "NSTEP" => nstep)
                             
    return tmpdir, prefix, temp
    
end

"""
   function makedict(savedir)

Load a data-file-schema.xml as julia dictionary.
"""  
function makedict(savedir)
    
    if isfile(savedir*"/data-file-schema.xml")
        filename=savedir*"/data-file-schema.xml"
    elseif isfile(savedir*"/data-file-schema.xml.gz")
        filename=savedir*"/data-file-schema.xml.gz"
    else
        println("ERROR warning missing data-file-schema.xml or data-file-schema.xml.gz")
        println()
        filename=missing
    end


    f = gzopen(filename, "r")
    fs = read(f, String)
    close(f)
    
    d = xml_dict(fs)

    return d
end


"""
    function loadXML(savedir)

Load a QE SCF DFT calculation into a `dftout` object.
"""
function loadXML(savedir)

#    println("start loadXML")
    convert_ha_ryd = 2.0
    
    d= makedict(savedir)

    calctype = d["espresso"]["input"]["control_variables"]["calculation"]
    
    if "convergence_achieved" in keys(d["espresso"]["output"]["convergence_info"]["scf_conv"])
        convergence = d["espresso"]["output"]["convergence_info"]["scf_conv"]["convergence_achieved"]
        if convergence == "false" && calctype != "nscf"
            println("error dft $savedir convergence false")
            return missing
        end
    elseif "scf_error" in keys(d["espresso"]["output"]["convergence_info"]["scf_conv"])
        scf_error = parse(Float64, d["espresso"]["output"]["convergence_info"]["scf_conv"]["scf_error"])
        if abs(scf_error) > 1e-6
            println("error dft $savedir convergence false $scf_error")
            return missing
        end
    end
    

    stin = d["espresso"]["output"]["atomic_structure"]
    stout = d["espresso"]["output"]["atomic_structure"]
    a1 = stout["cell"]["a1"]
    a2 = stout["cell"]["a2"]
    a3 = stout["cell"]["a3"]

#    println("a1 ", a1)
#    println("a2 ", a2)
#    println("a3 ", a3)

    A = zeros((3,3))
    A[1,:] = parse_str_ARR_float(a1)
    A[2,:] = parse_str_ARR_float(a2)
    A[3,:] = parse_str_ARR_float(a3)
    
#    println(A)

    nat = parse(Int, stin[:nat])
    types = []
    pos = zeros(nat,3)

#    println(d["espresso"]["input"]["k_points_IBZ"])

    if "monkhorst_pack" in keys(d["espresso"]["input"]["k_points_IBZ"])
        nk1 = parse(Int, d["espresso"]["input"]["k_points_IBZ"]["monkhorst_pack"][:nk1])
        nk2 = parse(Int, d["espresso"]["input"]["k_points_IBZ"]["monkhorst_pack"][:nk2])
        nk3 = parse(Int, d["espresso"]["input"]["k_points_IBZ"]["monkhorst_pack"][:nk3])
    else
        nk1 = 0
        nk2 = 0
        nk3 = 0
    end
    
    prefix = d["espresso"]["input"]["control_variables"]["prefix"]
    outdir = d["espresso"]["input"]["control_variables"]["outdir"]    

    tot_charge = parse(Float64, d["espresso"]["input"]["bands"]["tot_charge"])
    #println("tot_charge $tot_charge")
    
    kgrid = [nk1, nk2,nk3]
#    println("nat ", nat)
#    println(st["atomic_positions"]["atom"])
    
    if nat == 1
        atom =stout["atomic_positions"]["atom"]
        push!(types, atom[:name])
        pos[1,:] = parse_str_ARR_float(atom[""])
        
    else
        for at = 1:nat
#            println(at)
            atom =stout["atomic_positions"]["atom"][at]
            push!(types, atom[:name])
            pos[at,:] = parse_str_ARR_float(atom[""])
        end
    end
    
    coords_crys = pos * inv(A)

    if "forces" in keys(d["espresso"]["output"])
        out = d["espresso"]["output"]["forces"]
        f_lin = parse_str_ARR_float(out[""])
        forces = transpose(reshape(f_lin , 3,nat)) * convert_ha_ryd
    else
        forces = -99.0*ones(nat, 3)
    end

    if "stress" in keys(d["espresso"]["output"])
    
        out = d["espresso"]["output"]["stress"]
        s_lin = parse_str_ARR_float(out[""])
        stress = transpose(reshape(s_lin , 3,3)) * convert_ha_ryd
    else
        stress = -99.0*ones(3,3)
    end
        
    out = d["espresso"]["output"]["total_energy"]
    energy = parse(Float64,out["etot"]) * convert_ha_ryd
    energy_smear =  0.0
    if "demet" in keys(out)
        energy_smear = parse(Float64,out["demet"]) * convert_ha_ryd
    end

    nspin = 1
    mag_tot = 0.0
    mag_abs = 0.0
    if "spin" in keys(d["espresso"]["input"])
        lsda = d["espresso"]["input"]["spin"]["lsda"]
        if lsda == "true"
            nspin = 2
        end
        if "magnetization" in keys(d["espresso"]["output"])
            if "total" in keys(d["espresso"]["output"]["magnetization"])
                mag_tot = parse(Float64, d["espresso"]["output"]["magnetization"]["total"])
            else
                mat_tot = 0.0
            end
            mag_abs = parse(Float64, d["espresso"]["output"]["magnetization"]["absolute"])
        end
    end

    bs = loadXML_bs(d)
    bs.kgrid=kgrid
    
    d = DFToutMod.makedftout(A, coords_crys, types, energy, energy_smear, forces, stress, bs, prefix=prefix, outdir=outdir, tot_charge=tot_charge, nspin=nspin, mag_tot=mag_tot, mag_abs=mag_abs)

#    println("end loadXML")

    return d
    
end

"""
    function loadXML_bs(savedir::String)

Load `bandstructure` from QE xml file.
"""
function loadXML_bs(savedir::String)

    d=makedict(savedir)

    return loadXML_bs(d)
    
end
    

"""
    function loadXML_bs(savedir::String)

Load `bandstructure` from QE xml file that was converted to a dict already
"""
function loadXML_bs(d)

    convert_ha_ryd = 2.0
    
    out = d["espresso"]["output"]["band_structure"]
    nspin = 1
    if "nbnd" in keys(out)
        nbnd = parse(Int,out["nbnd"])
        nspin = 1
    elseif "nbnd_up" in keys(out)
        nspin = 2
        nbnd_up = parse(Int,out["nbnd_up"])
        nbnd_dn = parse(Int,out["nbnd_dw"])

        if nbnd_up != nbnd_up
            println("band numbers don't match for some reason !!!!! $nbnd_up $nbnd_dn")
        end
        nbnd = nbnd_up
    end



    nelec = parse(Float64,out["nelec"])
    efermi = parse(Float64,out["fermi_energy"]) * convert_ha_ryd


    st = d["espresso"]["output"]
    
    b1 = st["basis_set"]["reciprocal_lattice"]["b1"]
    b2 = st["basis_set"]["reciprocal_lattice"]["b2"]
    b3 = st["basis_set"]["reciprocal_lattice"]["b3"]

    B = zeros((3,3))
    B[1,:] = parse_str_ARR_float(b1)
    B[2,:] = parse_str_ARR_float(b2)
    B[3,:] = parse_str_ARR_float(b3)

    
    nks = parse(Int, out["nks"])
    bandstruct = zeros(nks, nbnd, nspin)
    kpts = zeros(nks, 3)
    weights = zeros(nks)    
    for b = 1:nks
#        println("out[ks_energies]")
#        println(out["ks_energies"])
        if b in keys(out["ks_energies"])
            ks = out["ks_energies"][b]
        else
            ks = out["ks_energies"]
        end

        weights[b] = parse(Float64, ks["k_point"][:weight] )
        kpts[b,:] = parse_str_ARR_float(ks["k_point"][""])
        if nspin == 1
            bandstruct[b,:,1] = parse_str_ARR_float(ks["eigenvalues"][""])
        else
            eig = parse_str_ARR_float(ks["eigenvalues"][""])
            bandstruct[b,:,1] = eig[1:nbnd] #up
            bandstruct[b,:,2] = eig[1+nbnd:2*nbnd] #down
        end
    end

    kpts = kpts * inv(B)
    
    bandstruct *= convert_ha_ryd

    bs = DFToutMod.makebs(nelec,efermi, kpts, weights, [0,0,0], bandstruct; nspin=nspin)

    return bs
    
end



end #module


############################## General DFT functions
"""
    module DFT

This is the generic DFT interface. Only QE is currently implemented however.
"""
module DFT
"""
Scripts to run DFT codes 
"""
#using ..CrystalMod
using ..QE
using ..CrystalMod:crystal

#include("Atomdata.jl")
#using ..Atomdata:atoms


"""
    function runSCF(crys::crystal; inputstr=missing, prefix=missing, tmpdir="./", directory="./", functional="PBESOL", wannier=0, nprocs=1, code="QE", skip=false, calculation="scf", dofree="all", tot_charge = 0.0, smearing = missing, magnetic=false, cleanup=false, use_backup=false)

Workflow for generic DFT SCF calculation. `code` can only by "QE"
"""
function runSCF(crys::crystal; inputstr=missing, prefix=missing, tmpdir="./", directory="./", functional="PBESOL", wannier=0, nprocs=1,procs=1, code="QE", skip=false, calculation="scf", dofree="all", tot_charge = 0.0, smearing = missing, magnetic=false, cleanup=false, use_backup=false, grid=missing, klines=missing, nstep=30, startingpot=missing, hybrid=false)
    
    nprocs = max(nprocs, procs)

    if ismissing(smearing)
        if calculation=="scf"
            smearing = 0.01
        else
            smearing = 0.02
        end
    end

    println(functional)
    if code == "QE"
        if ismissing(prefix)
            prefix="qe"
        end
        qeout = missing
        try
            qeout = QE.runSCF(crys, inputstr, prefix, tmpdir, directory, functional, wannier, nprocs, skip, calculation, dofree, tot_charge, smearing, magnetic, cleanup, use_backup, grid, klines, nstep, startingpot, hybrid)
            return qeout
        catch
            println("WARNING failure, restart qe with higher smearing, hope that helps!!")
            qeout = QE.runSCF(crys, inputstr, prefix, tmpdir, directory, functional, wannier, nprocs, skip, calculation, dofree, tot_charge, smearing*5, magnetic, cleanup, true, grid, klines, nstep, startingpot, hybrid)
            return qeout

        end

    else
        println("code variable not recognized", code)
        println("currently supported: QE")
    end        


    
end
    
end #end RunDFT

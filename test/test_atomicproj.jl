using ThreeBodyTB
using Test
using Suppressor

TESTDIR=ThreeBodyTB.TESTDIR


@testset "test atomicproj" begin

    
    
    @suppress begin

        units_old = ThreeBodyTB.set_units()
        ThreeBodyTB.set_units(both="atomic")


        filname = "$TESTDIR/data_forces/znse.in_vnscf_vol_2/projham_K.xml.gz"
        tbc_ref = ThreeBodyTB.TB.read_tb_crys_kspace(filname)

        fil = "$TESTDIR/data_forces/znse.in_vnscf_vol_2/"
        dft = ThreeBodyTB.QE.loadXML(fil*"/qe.save/")

        tbc, tbck, proj_warn = ThreeBodyTB.AtomicProj.projwfc_workf(dft, directory=fil, nprocs=1, writefile=missing, writefilek=missing, skip_og=true, skip_proj=true, freeze=true, localized_factor = 0.15, only_kspace = false, screening = 1.0 )


        
        energy_tot, tbc_, conv_flag = scf_energy(tbc)
        
        @test (energy_tot - dft.atomize_energy) < 1e-4

#        energy_tot_ref, tbc_, conv_flag = scf_energy(tbc_ref)
        energy_tot_ref, X = ThreeBodyTB.TB.get_energy_electron_density_kspace(tbc_ref)

        @test (energy_tot_ref - dft.atomize_energy) < 1e-4

        ThreeBodyTB.set_units(energy = units_old[1], length=units_old[2])
        
    end

    
end

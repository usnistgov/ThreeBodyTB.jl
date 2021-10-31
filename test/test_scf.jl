using Test
using ThreeBodyTB
using Suppressor

#include("../includes_laguerre.jl")
#include("../Ewald.jl")

TESTDIR=ThreeBodyTB.TESTDIR

function test_basics()

    @testset "test basic scf Hk TB manipulation" begin

         @suppress begin
#        if true

            units_old = ThreeBodyTB.set_units()
            ThreeBodyTB.set_units(both="atomic")

            
            types=["Al"];
            
            #positions, crystal units
            pos=zeros((1,3));
            
            #lattice vectors, in Bohr units currently
            A=[ [3.7962751610 3.7962751610 0]; [3.7962751610 0 3.7962751610 ]; [ 0 3.7962751610 3.7962751610]];
            
            #makes the crystal
            c=makecrys(A, pos, types, units="Bohr")

            energy, tbc, flag = scf_energy(c)
            energyS, tbcS, flagS = scf_energy(c, mixing_mode=:simple)

            
#            @test abs(energy - -0.28986347423986025) < 1e-2
            @test abs(energy - -0.28707042) < 1e-2
            @test abs(energyS - energy) < 1e-5

            energy_fft = ThreeBodyTB.TB.calc_energy_fft(tbc)

            directgap, indirectgap, gaptype, bandwidth = ThreeBodyTB.BandStruct.band_summary(tbc)
            @test gaptype == :metal
            @test abs(indirectgap) < 1e-5
            @test abs(directgap) < 1e-5
            @test abs(bandwidth - 0.8334073432114558) < 1e-5

            
            @test abs(energy - energy_fft) < 1e-5

            energy_calc = ThreeBodyTB.TB.calc_energy(tbc)

            @test abs(energy - energy_calc) < 1e-5

            energy2, force, stress, tbc2 = scf_energy_force_stress(tbc)
            
            @test abs(energy - energy2) < 1e-4  
            @test sum(abs.(force)) < 1e-4
            @test sum(abs.(stress))  < 1e-2

            vects, vals, ham, S,e =  Hk(tbc, [0 0 0 ])

#            @test abs( (vals[2] - vals[1]) - 1.9352469029181867) < 1e-2   #basic eigen value testing
            @test abs( (vals[2] - vals[1]) - 1.9585164167567397) < 1e-2   #basic eigen value testing
            @test abs( (vals[3] - vals[2]) ) < 1e-4
            @test abs( (vals[4] - vals[2]) ) < 1e-4

            vects2, vals2, ham2, S2,e2 =  Hk(tbc.tb, [0 0 0 ])
            
            @test sum(abs.(ham2 - ham)) < 1e-7

            tbck = ThreeBodyTB.TB.read_tb_crys_kspace("$TESTDIR/data_forces/al_fcc_projham_K.xml.gz")

            directgap, indirectgap, gaptype, bandwidth = ThreeBodyTB.BandStruct.band_summary(tbck)
            @test gaptype == :metal
            @test abs(indirectgap) < 1e-5
            @test abs(directgap) < 1e-5
            @test abs(bandwidth - 0.8418566664739338) < 1e-2
            
            vects3, vals3, ham3, S3,e3 =  Hk(tbck, [0 0 0 ])

            @test sum(abs.(vals3[1] - vals[1])) < 0.05

            tb = deepcopy(tbc.tb)
            ThreeBodyTB.TB.trim(tb)

            vectsT, valsT, hamT, ST,eT =  Hk(tb, [0 0 0 ])

            @test sum(abs.(valsT - vals2)) < 2e-2
#            println("VALST $valsT")
#            println("VALS2 $vals2")

            data_onsite, data_arr = ThreeBodyTB.TB.organizedata(tbc.crys, tbc.tb)

#            @test abs(data_onsite[1,5] - -0.677103) < 1e-2  #onsite s orbital
#            @test abs(data_onsite[1,7] - 5.36874) < 1e-2 #n.n. distance


            @test abs(data_onsite[1,5] - -0.696210333954235) < 1e-2  #onsite s orbital
            @test abs(data_onsite[1,7] - 5.368743819186305) < 1e-2 #n.n. distance
            @test abs(data_onsite[1,11] - 1.0) < 1e-3 # overlap onsite

            h1, dq = ThreeBodyTB.TB.get_h1(tbc)
            
            @test dq[1] < 1e-5

            en, eden, VECTS, VALS, err = ThreeBodyTB.TB.get_energy_electron_density_kspace(tbck)

            @test abs(en - -0.28962628722369677) < 1e-3

            tbck2 = ThreeBodyTB.SCF.remove_scf_from_tbc(tbck)
            tbc2 = ThreeBodyTB.SCF.remove_scf_from_tbc(tbc)

            en2, eden2, VECTS2, VALS2, err2 = ThreeBodyTB.TB.get_energy_electron_density_kspace(tbck2)

            @test abs(en - en2) < 1e-3

            H = ThreeBodyTB.Atomdata.atoms["H"]

            @test H.Z == 1 

            energies, dos, pdos, names  =  ThreeBodyTB.DOS.gaussian_dos(tbc, do_display=false)

            @test length(energies) == length(dos)

            energies, dos, pdos, names  =  ThreeBodyTB.DOS.gaussian_dos(tbck, do_display=false)

            @test length(energies) == length(dos)
            @test abs(sum(dos) * (energies[2] - energies[1]) - 4.0) < 1e-2
            ThreeBodyTB.set_units(energy = units_old[1], length=units_old[2])

            
        end
    end
end


test_basics()

using Test
using ThreeBodyTB
using Suppressor

#include("../includes_laguerre.jl")
#include("../Ewald.jl")

function testU()

    @testset "testing U" begin

#        Atomdata.atoms["X"].energy_offset = 0.0
#        Atomdata.atoms["Xa"].energy_offset = 0.0
        @suppress begin

            units_old = ThreeBodyTB.set_units()
            ThreeBodyTB.set_units(both="atomic")


            cscl = makecrys([1.0 0 0; 0 1.0 0; 0 0 1.0]*10.0, [0 0 0; 0.5 0.5 0.5], ["X", "Xa"], units="Bohr");
            coefXX = ThreeBodyTB.CalcTB.make_coefs(Set((:X, :X)),2, fillzeros=true)
            coefXXa = ThreeBodyTB.CalcTB.make_coefs(Set((:X, :Xa)),2, fillzeros=true)
            coefXaXa = ThreeBodyTB.CalcTB.make_coefs(Set((:Xa, :Xa)),2, fillzeros=true)

            database = Dict()
            database[(:X, :X)] = coefXX
            database[(:X, :Xa)] = coefXXa
            database[(:Xa, :X)] = coefXXa
            database[(:Xa, :Xa)] = coefXaXa
            database["scf"] = true
            database["SCF"] = true

            tbc = ThreeBodyTB.CalcTB.calc_tb_fast(cscl, database, use_threebody=false, use_threebody_onsite=false)
            tbc.gamma = ThreeBodyTB.Ewald.electrostatics_getgamma(cscl, onlyU=true) #only use onsite U, no electrostatics

            @test size(tbc.gamma) == (2,2)
            @test tbc.gamma[1,1] == 0.5
            @test tbc.gamma[2,2] == 0.5
            @test tbc.gamma[1,2] == 0.0

            tbc.scf = true

            energy, efermi, e_den, dq, error_flag = ThreeBodyTB.SCF.scf_energy(tbc, mix=0.1, grid=[1 1 1], conv_thr = 1e-8, iters=1000)
            println("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

            #        energy, efermi, e_den, dq = SCF.scf_energy(tbc, mix=0.5, grid=[1 1 1], conv_thr = 1e-8, mixing_mode=:pulay, iters=1000)
            
            #        println("energy $energy")
            #        println("efermi $efermi")
            #        println("dq $dq")

            @test abs(energy - 0.5) ≤ 1e-6
            @test abs(efermi - 1.5) ≤ 1e-5
            @test abs(dq[1] - 1.0) ≤ 1e-5
            @test abs(dq[2] + 1.0) ≤ 1e-5
            @test abs(e_den[1] - 1.0) ≤ 1e-5
            @test abs(e_den[2] - 0.0) ≤ 1e-5


            cscl2 = makecrys([1.0 0 0; 0 1.0 0; 0 0 1.0]*10.0, [0 0 0; 0.5 0.5 0.5], ["X", "X"], units="Bohr");
            tbc2 = ThreeBodyTB.CalcTB.calc_tb_fast(cscl2, database, use_threebody=false, use_threebody_onsite=false)
            tbc2.gamma = ThreeBodyTB.Ewald.electrostatics_getgamma(cscl, onlyU=true) #only use onsite U, no electrostatics
            tbc2.scf = true
            energy, efermi, e_den, dq, error_flag = ThreeBodyTB.SCF.scf_energy(tbc2, mix=0.75, grid=[1 1 1], conv_thr = 1e-8)
            @test abs(energy - -0.011283791670955126) ≤ 1e-6

            ThreeBodyTB.set_units(energy = units_old[1], length=units_old[2])

        end
    end
end

function test_nacl()

    @testset "testing nacl" begin
        @suppress begin
            units_old = ThreeBodyTB.set_units()
            ThreeBodyTB.set_units(both="atomic")

            cscl = makecrys([1.0 0 0; 0 1.0 0; 0 0 1.0]*19.5509319833, [0 0 0; 0.5 0.5 0.5], ["Na", "Cl"], units="Bohr");

            #no tb params
            coefXX = ThreeBodyTB.CalcTB.make_coefs(Set(("Na", "Na")),2, fillzeros=true)
            coefXXa = ThreeBodyTB.CalcTB.make_coefs(Set(("Na", "Cl")),2, fillzeros=true)
            coefXaXa = ThreeBodyTB.CalcTB.make_coefs(Set(("Cl", "Cl")),2, fillzeros=true)

            database = Dict()
            database[(:Na, :Na)] = coefXX
            database[(:Na, :Cl)] = coefXXa
            database[(:Cl, :Na)] = coefXXa
            database[(:Cl, :Cl)] = coefXaXa
            database["scf"] = true
            database["SCF"] = true

            tbc = ThreeBodyTB.CalcTB.calc_tb_fast(cscl, database, use_threebody=false, use_threebody_onsite=false)

            tbc.scf = true

            energy, efermi, e_den, dq, error_flag = ThreeBodyTB.SCF.scf_energy(tbc, mix=0.005, grid=[1 1 1], conv_thr = 1e-8, smearing = 0.005, iters = 1000)

            #        println("energy $energy")
            #        println("efermi $efermi")
            #        println("e_den $e_den")
            #        println("dq $dq")


            tbc_noscf = ThreeBodyTB.CalcTB.calc_tb_fast(cscl, database, use_threebody=false, use_threebody_onsite=false)
            tbc_noscf.scf = false
            energy_no, efermi_no, e_den_no, dq_no, error_flag = ThreeBodyTB.SCF.scf_energy(tbc_noscf, mix=0.005, grid=[1 1 1], conv_thr = 1e-8, smearing = 0.005, iters = 1000)

            #        println("energy_no $energy_no")
            #        println("efermi_no $efermi_no")
            #        println("e_den_no $e_den_no")
            #        println("dq_no $dq_no")


            #        println("the actual DFT atomization  energy of this configuration is -0.13722807000000614")
            #        println("energy_scf is $energy")
            #        println("energy_noscf is $energy_no")

            @test abs(-0.13722807000000614 - energy) < 5e-2

            ThreeBodyTB.set_units(energy = units_old[1], length=units_old[2])

        end
    end
end

function test_nacl2()

    @testset "testing nacl2" begin
        @suppress begin
            units_old = ThreeBodyTB.set_units()
            ThreeBodyTB.set_units(both="atomic")

            rs = makecrys([1.0 1.0 0; 1.0 0 1.0; 0 1.0 1.0]*10.6062839325, [0 0 0; 0.5 0.5 0.5], ["Na", "Cl"], units="Bohr");

            #no tb params
            coefXX = ThreeBodyTB.CalcTB.make_coefs(Set(("Na", "Na")),2, fillzeros=true)
            coefXXa = ThreeBodyTB.CalcTB.make_coefs(Set(("Na", "Cl")),2, fillzeros=true)
            coefXaXa = ThreeBodyTB.CalcTB.make_coefs(Set(("Cl", "Cl")),2, fillzeros=true)

            database = Dict()
            database[(:Na, :Na)] = coefXX
            database[(:Na, :Cl)] = coefXXa
            database[(:Cl, :Na)] = coefXXa
            database[(:Cl, :Cl)] = coefXaXa
            database["scf"] = true
            database["SCF"] = true

            tbc = ThreeBodyTB.CalcTB.calc_tb_fast(rs, database, use_threebody=false, use_threebody_onsite=false)

            tbc.scf = true

            energy, efermi, e_den, dq, error_flag = ThreeBodyTB.SCF.scf_energy(tbc, mix=0.005, grid=[1 1 1], conv_thr = 1e-8, smearing = 0.01, iters = 1000)

            #        println("energy $energy")
            #        println("efermi $efermi")
            #        println("e_den $e_den")
            #        println("dq $dq")

            tbc_noscf = ThreeBodyTB.CalcTB.calc_tb_fast(rs, database, use_threebody=false, use_threebody_onsite=false)
            tbc_noscf.scf = false
            energy_no, efermi_no, e_den_no, dq_no, error_flag = ThreeBodyTB.SCF.scf_energy(tbc_noscf, mix=0.5, grid=[1 1 1], conv_thr = 1e-8, smearing = 0.01, iters = 1000) #, mixing_mode=:simple)

            #        println("energy_no $energy_no")
            #        println("efermi_no $efermi_no")
            #        println("e_den_no $e_den_no")
            #        println("dq_no $dq_no")


            #        println("the actual DFT atomization energy of this configuration is -0.23582715999999948")
            #        println("energy_scf is $energy")
            #        println("energy_noscf is $energy_no")

            @test abs(-0.23582715999999948 - energy) < 6e-2

            ThreeBodyTB.set_units(energy = units_old[1], length=units_old[2])


        end
    end
end

testU()
#println("sleep ")
#sleep(3)
test_nacl()
#println("sleep ")
#sleep(3)
test_nacl2()

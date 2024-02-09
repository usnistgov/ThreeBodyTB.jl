using Test
using ThreeBodyTB
using Suppressor


TESTDIR=ThreeBodyTB.TESTDIR

@suppress begin
    ThreeBodyTB.set_units(both="atomic")
end


function test_force_sym()

    @testset "testing force sym" begin
        @suppress begin

            #hexagonal symmetry is tricky
            c = makecrys([4 0 0; -0.5*4 4*sqrt(3)/2 0;0 0 4.2]*2.0, [0.5 0.5 0.5; 0.7 0.5 0.5; 0.5 0.7 0.5; 0.3 0.3 0.5], ["N", "H", "H", "H"]);

            enS,  fcartS, stressS,  tbcS = scf_energy_force_stress(c, use_sym=true)
            en,  fcart, stress,  tbc = scf_energy_force_stress(c, use_sym=false)

            @test abs(enS - en) <  1e-7
            @test sum(abs.(fcart - fcartS)) < 1e-7
            @test sum(abs.(stressS - stress)) < 1e-5

            @test sum(abs.(tbc.eden - tbcS.eden)) < 1e-5
            
        end
    end
end

function test_force_sym_single()

    @testset "testing force single" begin
        @suppress begin
            
            #hexagonal symmetry is tricky
            c = makecrys([4 0 0; -0.5*4 4*sqrt(3)/2 0;0 0 4.2]*2.0, [0.5 0.5 0.5; 0.7 0.5 0.5; 0.5 0.7 0.5; 0.3 0.3 0.5], ["N", "H", "H", "H"]);

            en, tbc, flag = scf_energy(c)
            
            en, force, stress = ThreeBodyTB.Force_Stress.get_energy_force_stress_fft_LV_sym(deepcopy(tbc), ThreeBodyTB.ManageDatabase.database_cached)
            enS, forceS, stressS = ThreeBodyTB.Force_Stress.get_energy_force_stress_fft_LV_sym_SINGLE(deepcopy(tbc), ThreeBodyTB.ManageDatabase.database_cached)
            
            @test abs(enS - en) <  1e-6
            @test sum(abs.(force - forceS)) < 1e-5
            @test sum(abs.(stressS - stress)) < 1e-5

            # println("res")
            #println([enS , en])
            #println(force - forceS)
            #println(stressS - stress)
            
        end
    end
end


test_force_sym()
test_force_sym_single()

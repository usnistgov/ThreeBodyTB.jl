using Test
using ThreeBodyTB
using Suppressor

function test1()
    tbc_list = []
    @testset "testing force fast SPARSE" begin

        @suppress begin
            units_old = ThreeBodyTB.set_units()
            for SCF in [true, false]
                ThreeBodyTB.set_units(both="atomic")
                
                c = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.45 0 0 ], ["Hx", "Hx"], units="Bohr");

                database = Dict()
                database[(:Hx, :Hx)] = ThreeBodyTB.CalcTB.make_coefs(Set(["Hx", "Hx"]), 2)
                
                database["scf"] = SCF
                database["SCF"] = SCF

                en, f_cart, stress, tbc = scf_energy_force_stress(c, database=database, repel=false)
                
                enFD, fFD = ThreeBodyTB.Force_Stress.finite_diff(c, database,1, 1, stress_mode=false, repel=false);
                
                @test abs(fFD - f_cart[1,1]) < 1e-5

                enFD, sFD = ThreeBodyTB.Force_Stress.finite_diff(c, database,1, 1, stress_mode=true, repel=false);
                
                @test abs(sFD - stress[1,1]) < 1e-5

                tbc_sparse = ThreeBodyTB.CalcTB.calc_tb_LV_sparse(c, database)
                en_sparse, force_sparse, stress_sparse = ThreeBodyTB.Force_Stress.get_energy_force_stress_fft_LV_sym_SINGLE(tbc_sparse, database)

                @test sum(abs.(f_cart - force_sparse)) < 1e-8
                @test sum(abs.(stress - stress_sparse)) < 1e-8
                @test abs(en - en_sparse) < 1e-8
                
                
            end
            ThreeBodyTB.set_units(energy = units_old[1], length=units_old[2])
        
        end
    end
    return tbc_list
end


tbc_list = test1();


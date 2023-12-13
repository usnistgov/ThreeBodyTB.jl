using Test
using ThreeBodyTB
using Suppressor
#=
include("../Utility.jl")
include("../BandTools.jl")
include("../Atomic.jl")
include("../Atomdata.jl")
include("../Crystal.jl")
include("../DFToutMod.jl")
include("../TB.jl")
include("../CalcTB.jl")
include("../FitTB.jl")
=#

#basic loading

#include("../includes_laguerre.jl")

function test1()
    tbc_list = []
    @testset "testing force fast" begin

        @suppress begin
            for SCF in [true, false]
                units_old = ThreeBodyTB.set_units()
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
                
            end
            ThreeBodyTB.set_units(energy = units_old[1], length=units_old[2])
        
        end
    end
    return tbc_list
end


tbc_list = test1();


function test1_charge()
    tbc_list = []
    @testset "testing force fast charge" begin

        @suppress begin
            units_old = ThreeBodyTB.set_units()
            ThreeBodyTB.set_units(both="atomic")
        
            c = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.45 0 0 ], ["Hx", "Hx"], units="Bohr");

            database = Dict()
            database[(:Hx, :Hx)] = ThreeBodyTB.CalcTB.make_coefs(Set(["Hx", "Hx"]), 2)
            
            database["scf"] = true
            database["SCF"] = true

            en, f_cart, stress, tbc = scf_energy_force_stress(c, database=database, repel=false, tot_charge=0.1)
            
            enFD, fFD = ThreeBodyTB.Force_Stress.finite_diff(c, database,1, 1, stress_mode=false, repel=false, tot_charge=0.1);
            
            @test abs(fFD - f_cart[1,1]) < 1e-5

            enFD, sFD = ThreeBodyTB.Force_Stress.finite_diff(c, database,1, 1, stress_mode=true, repel=false, tot_charge=0.1);
            
            @test abs(sFD - stress[1,1]) < 1e-5
            
            ThreeBodyTB.set_units(energy = units_old[1], length=units_old[2])
        
        end
    end
    return tbc_list
end

test1_charge();

Nothing;

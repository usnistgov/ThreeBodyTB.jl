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

#            database[(:Hx, :Hx)].U = 0.2
#            database[(:Hx, :Hx)].U3 = 0.1
            
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

function test2_charge()
    tbc_list = []
    @testset "testing force fast charge" begin

        @suppress begin
            units_old = ThreeBodyTB.set_units()
            ThreeBodyTB.set_units(both="atomic")
        
            c = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.45 0 0 ], ["Hx", "H"], units="Bohr");

            #make a second type of "fictitious" hydrogen with different eigenvalue
            ThreeBodyTB.Atomdata.atoms["Hx"] = ThreeBodyTB.AtomicMod.makeatom("Hx" , 100, 1, 1, 1.0079,   1.0,   0,   [:s ],     -0.90531367, [-1.5], 0.0, 0.9,0.1)
            ThreeBodyTB.Atomdata.atoms[:Hx ] = ThreeBodyTB.Atomdata.atoms["Hx"]

            
            database = Dict()
            database[(:Hx, :Hx)] = ThreeBodyTB.CalcTB.make_coefs(Set(["Hx", "Hx"]), 2)
            database[(:Hx, :H)] = ThreeBodyTB.CalcTB.make_coefs(Set(["Hx", "H"]), 2, datU = 10.0 * ones(ThreeBodyTB.CalcTB.n_ufit * 2))
            database[(:H, :Hx)] = ThreeBodyTB.CalcTB.make_coefs(Set(["Hx", "H"]), 2, datU = 10.0 * ones(ThreeBodyTB.CalcTB.n_ufit * 2))
            database[(:H, :H)] = ThreeBodyTB.CalcTB.make_coefs(Set(["H", "H"]), 2)

#            database[(:Hx, :Hx)].U = 0.2
#            database[(:Hx, :Hx)].U3 = 0.1
            
            database["scf"] = true
            database["SCF"] = true

            en, f_cart, stress, tbc = scf_energy_force_stress(c, database=database, repel=false, tot_charge=0.0)
            
            enFD, fFD = ThreeBodyTB.Force_Stress.finite_diff(c, database,1, 1, stress_mode=false, repel=false, tot_charge=0.0);
            
            @test abs(fFD - f_cart[1,1]) < 1e-5

            enFD, sFD = ThreeBodyTB.Force_Stress.finite_diff(c, database,1, 1, stress_mode=true, repel=false, tot_charge=0.0);
            
            @test abs(sFD - stress[1,1]) < 1e-5
            
            ThreeBodyTB.set_units(energy = units_old[1], length=units_old[2])

            #put things back
            ThreeBodyTB.Atomdata.atoms["Hx"] = ThreeBodyTB.AtomicMod.makeatom("Hx" , 100, 1, 1, 1.0079,   1.0,   0,   [:s ],     -0.90531367, [-5.5], 0.0, 0.9194820210526328, 0.0)
            ThreeBodyTB.Atomdata.atoms[:Hx ] = ThreeBodyTB.Atomdata.atoms["Hx"]
            
            
        end
    end
    return tbc_list
end

test2_charge();

Nothing;

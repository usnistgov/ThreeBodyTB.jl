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

    @testset "testing classical elemental" begin
        @suppress begin

            database = Dict()
            coef2 = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Al], 2)
            coef2.datH[:] = [1.0, 0.8, 0.6, 0.4, 0.2, 0.1]
            database[(:Al, :Al)] = coef2
            
            c2 = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.15], ["Al", "Al"])
            v = 1.0;
            C = ThreeBodyTB.CrystalMod.crystal[];
            EN = Float64[];
            for x in 1:12;
                v = v * 1.12;
                for c in [c2];
                    push!(C, c*v);
                    en, _ = ThreeBodyTB.Classical.calc_energy_cl(c*v, database=database, use_threebody=false, use_em=false);
                    push!(EN, en / c.nat ); #normalize per atom
                end;
            end
            V, _ = ThreeBodyTB.ClassicalFit.prepare_fit_cl(C, use_threebody=false, get_force=false, use_em=false);
            x = V \ EN
            @test sum(abs.(x - [1.0, 0.8, 0.6, 0.4, 0.2, 0.1])) < 1e-8
            
            
        end
    end
end

function test3()
    @testset "testing classical elemental 3body" begin
        @suppress begin

            database = Dict()
            coef2 = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Al], 2)
            coef2.datH[:] = [1.0, -0.8, 0.6, -0.4, 0.2, -0.1]
            database[(:Al, :Al)] = coef2

            coef3 = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Al, :Al], 3)
            coef3.datH[:] = 1 ./ [1,2,2,2,3,3,3,4,5,5,5,6,6,6,-7,-7,-7,-7,-7,-7, 8, 8, 8]
            database[(:Al, :Al, :Al)] = coef3
            

            
            c2 = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10], ["Al", "Al"])
            c3a = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10; 0.1 0 0.05 ], ["Al", "Al", "Al"])
            c3 = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10; 0.1 0 0 ], ["Al", "Al", "Al"])

            v = 1.0;
            C = ThreeBodyTB.CrystalMod.crystal[];
            EN = Float64[];
            for x in 1:12;
                v = v * 1.12;
                for c in [c2, c3 ,c3a];
                    push!(C, c*v);
                    en, _ = ThreeBodyTB.Classical.calc_energy_cl(c*v, database=database, use_threebody=true, use_em = false)
                    push!(EN, en / c.nat);
                end;
            end


            V, _ = ThreeBodyTB.ClassicalFit.prepare_fit_cl(C, use_threebody=true, get_force=false, use_em = false);
            x = V \ EN

            @test sum(abs.(x[1:6] - [1.0, -0.8, 0.6, -0.4, 0.2, -0.1])) < 1e-8
            @test sum(abs.(x[7:13] - 1 ./ [1,2,3,4,5,6,-7] )) < 1e-8
            
            
        end
    end
end


function test3_2atoms()

    @testset "testing classical elemental 3body 2atoms" begin
        @suppress begin

            database = Dict()
            coef2a = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Al], 2)
            coef2b = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Ga], 2)
            coef2c = ThreeBodyTB.Classical.make_coefs_cl([:Ga, :Ga], 2)            
            coef2a.datH[:] = [1.0, -0.8, 0.6, -0.4, 0.2, -0.1]
            coef2b.datH[:] = [1.0, 0.8, 0.6, -0.4, 0.2, -0.1]
            coef2c.datH[:] = [1.0, 0.8, 0.6, -0.4, -0.2, 0.1]

            database[(:Al, :Al)] = coef2a

            database[(:Al, :Ga)] = coef2b
            database[(:Ga, :Al)] = coef2b

            database[(:Ga, :Ga)] = coef2c
            

            coef3 = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Al, :Al], 3)
            coef3.datH[:] = 1 ./ [1,2,2,2,3,3,3,4,5,5,5,6,6,6,-7,-7,-7,-7,-7,-7, 8, 8, 8] 
            database[(:Al, :Al, :Al)] = coef3

            coef3b = ThreeBodyTB.Classical.make_coefs_cl([:Ga, :Ga, :Ga], 3)
            coef3b.datH[:] = 1 ./ [1,2,2,2,3,3,3,4,5,5,5,6,6,6,-7,-7,-7,-7,-7,-7, 8, 8, 8] * 2.0 
            database[(:Ga, :Ga, :Ga)] = coef3b

            coef3c = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Ga, :Ga], 3)
            coef3c.datH[:] = ones(23)*0.37
            database[(:Al, :Ga, :Ga)] = coef3c
            database[(:Ga, :Al, :Ga)] = coef3c
            database[(:Ga, :Ga, :Al)] = coef3c

            database[(:Al, :Al, :Ga)] = coef3c
            database[(:Al, :Ga, :Al)] = coef3c
            database[(:Ga, :Al, :Al)] = coef3c
            

            
            c2 = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10], ["Al", "Al"])
            c2a = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10], ["Al", "Ga"])
            c2b = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10], ["Ga", "Ga"])
            
            c3a = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10; 0.1 0 0.05 ], ["Al", "Al", "Al"])
            c3 = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10; 0.1 0 0 ], ["Al", "Al", "Al"])

            c3aB = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10; 0.1 0 0.05 ], ["Ga", "Ga", "Ga"])
            c3B = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10; 0.1 0 0 ], ["Ga", "Ga", "Ga"])

            c3aC = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10; 0.1 0 0.05 ], ["Al", "Ga", "Ga"])
            c3C = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10; 0.1 0 0 ], ["Ga", "Al", "Ga"])

            c3aD = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10; 0.1 0 0.05 ], ["Al", "Al", "Ga"])
            c3D = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10; 0.1 0 0 ], ["Ga", "Al", "Al"])

            c3aE = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.11; 0.06 0 0.06 ], ["Al", "Al", "Ga"])
            c3E = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10; 0.11 0 0 ], ["Ga", "Al", "Al"])

            c3aF = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.11; 0.06 0 0.06 ], ["Al", "Ga", "Ga"])
            c3F = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10; 0.11 0 0 ], ["Ga", "Ga", "Al"])
            
            v = 1.0;C = ThreeBodyTB.CrystalMod.crystal[];  EN = Float64[];
            for x in 1:11; v = v * 1.12;
                for c in [c2, c2a, c2b, c3 ,c3a, c3aB, c3B, c3aC, c3C,c3aD, c3D, c3aE, c3E, c3aF, c3F];
                    push!(C, c*v);
                    en, _ = ThreeBodyTB.Classical.calc_energy_cl(c*v, database=database, use_threebody=true, use_em =false)                    
                    push!(EN, en / c.nat);
                end
            end

            database_NEW = ThreeBodyTB.ClassicalFit.do_fit_cl(C,ENERGIES=EN, use_threebody=true, use_em=false, lambda = 0.0);

            @test sum(abs.(database[(:Al, :Ga)].datH - database_NEW[(:Al, :Ga)].datH)) < 5e-3
            @test sum(abs.(database[(:Al, :Al, :Al)].datH - database_NEW[(:Al, :Al, :Al)].datH)) < 5e-3
            @test sum(abs.(database[(:Ga, :Al, :Al)].datH - database_NEW[(:Ga, :Al, :Al)].datH)) < 5e-3
            @test sum(abs.(database[(:Ga, :Ga, :Al)].datH - database_NEW[(:Ga, :Ga, :Al)].datH)) < 5e-3
            
            return database_NEW
            
        end
    end
end

function test3_3atoms()

    @testset "testing classical elemental 3body 3atoms" begin
        @suppress begin

            database = Dict()
            coef2a = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Al], 2)
            coef2b = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Ga], 2)
            coef2c = ThreeBodyTB.Classical.make_coefs_cl([:Ga, :Ga], 2)            
            coef2d = ThreeBodyTB.Classical.make_coefs_cl([:In, :In], 2)
            coef2e = ThreeBodyTB.Classical.make_coefs_cl([:In, :Al], 2)
            coef2f = ThreeBodyTB.Classical.make_coefs_cl([:In, :Ga], 2)            

            coef2a.datH[:] .= 0.0
            coef2b.datH[:] .= 0.0
            coef2c.datH[:] .= 0.0
            coef2d.datH[:] .= 0.0
            coef2e.datH[:] .= 0.0
            coef2f.datH[:] .= 0.0


            database[(:Al, :Al)] = coef2a
            database[(:Al, :Ga)] = coef2b
            database[(:Ga, :Al)] = coef2b
            database[(:Ga, :Ga)] = coef2c
            database[(:In, :In)] = coef2d
            database[(:In, :Al)] = coef2e
            database[(:Al, :In)] = coef2e
            database[(:In, :Ga)] = coef2f
            database[(:Ga, :In)] = coef2f
            

            coef3 = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Ga, :In], 3)
            coef3.datH[:] = 1 ./ [1,2,2,2,3,3,3,4,5,5,5,6,6,6,-7,-7,-7,-7,-7,-7, 8, 8, 8] 
            database[(:Al, :Ga, :In)] = coef3
            database[(:Al, :In, :Ga)] = coef3
            database[(:Ga, :Al, :In)] = coef3
            database[(:In, :Al, :Ga)] = coef3
            database[(:In, :Ga, :Al)] = coef3
            database[(:Ga, :In, :Al)] = coef3

            c3aC = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10; 0.1 0 0.05 ], ["Al", "Ga", "In"])
            c3C = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10; 0.1 0 0 ], ["Al", "Ga", "In"])

            c3aE = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.11; 0.06 0 0.07 ], ["Ga", "In", "Al"])
            c3E = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10; 0.13 0 0 ], ["Ga", "In", "Al"])

            c3aF = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.11; 0.06 0 0.08 ], ["Al", "In", "Ga"])
            c3F = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.10; 0.12 0 0 ], ["Al", "In", "Ga"])
            
            v = 1.0;C = ThreeBodyTB.CrystalMod.crystal[];  EN = Float64[];
            for x in 1:11; v = v * 1.12;
                #                for c in [c2, c2a, c2b, c3 ,c3a, c3aB, c3B, c3aC, c3C,c3aD, c3D, c3aE, c3E, c3aF, c3F];
                for c in [ c3aC, c3C, c3aE, c3E, c3aF, c3F]
                    push!(C, c*v);
                    en, _ = ThreeBodyTB.Classical.calc_energy_cl(c*v, database=database, use_threebody=true, use_em = false)
                    push!(EN, en / c.nat);
                end
            end

            database_NEW = ThreeBodyTB.ClassicalFit.do_fit_cl(C,ENERGIES=EN, use_threebody=true, lambda = 0.0);

            @test sum(abs.(database[(:Al, :Ga)].datH - database_NEW[(:Al, :Ga)].datH)) < 1e-3
            @test sum(abs.(database[(:In, :Ga)].datH - database_NEW[(:In, :Ga)].datH)) < 1e-3
            @test sum(abs.(database[(:Al, :Ga, :In)].datH - database_NEW[(:Al, :Ga, :In)].datH)) < 1e-3
            @test sum(abs.(database[(:Al, :In, :Ga)].datH - database_NEW[(:Al, :In, :Ga)].datH)) < 1e-3
            
            return database_NEW
            
        end
    end
end


test1();
test3();
database_NEW2 = test3_2atoms();
database_NEW3 = test3_3atoms();
Nothing;

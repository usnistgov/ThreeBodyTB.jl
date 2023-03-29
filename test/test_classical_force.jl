using Test
using ThreeBodyTB
using Suppressor
using LinearAlgebra

function test1()

    @testset "testing classical force" begin
        @suppress begin

            database = Dict()
            coef2 = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Al], 2)
            coef2.datH[:] = [1.0, 0.8, 0.6, 0.4, 0.2, 0.1]
            database[(:Al, :Al)] = coef2
            
            c2 = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.1], ["Al", "Al"])
            c2a = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.100001], ["Al", "Al"])

            energy2, force2, stress2 = ThreeBodyTB.Classical.energy_force_stress(c2, database=database);
            energy2a, force2a, stress2a = ThreeBodyTB.Classical.energy_force_stress(c2a, database=database);
            
            finitediff = (energy2a - energy2) / (0.100001 - 0.1) / 30.0

            @test abs(finitediff + force2[2,3]) < 1e-5
            @test abs(finitediff + stress2[3,3] * abs(det(c2.A))/3) < 1e-5

            println("finitediff ", finitediff)
            println("force2")
            println(force2)
            println("stress")
            println(stress)
        end
    end
end


function test3()

    @testset "testing classical force 3body" begin
        @suppress begin

            database = Dict()
            coef2 = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Al], 2)
            coef2.datH[:] = [1.0, 0.8, 0.6, 0.4, 0.2, 0.1]
            database[(:Al, :Al)] = coef2

            coef3 = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Al, :Al], 3)
            coef3.datH[:] = 1 ./ [1,2,2,2,3,3,3,4,5,5,5,6,6,6,-7,-7,-7,-7,-7,-7]
            database[(:Al, :Al, :Al)] = coef3
            
            c2 = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.1; 0.1 0 0], ["Al", "Al", "Al"])
            c2a = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.100001; 0.1 0 0], ["Al", "Al", "Al"])
            c2b = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0.000001 0 0.10000; 0.1 0 0], ["Al", "Al", "Al"])

            energy2, force2, stress2 = ThreeBodyTB.Classical.energy_force_stress(c2, database=database);
            energy2a, force2a, stress2a = ThreeBodyTB.Classical.energy_force_stress(c2a, database=database);
            energy2b, force2b, stress2b = ThreeBodyTB.Classical.energy_force_stress(c2b, database=database);
            
            finitediff = (energy2a - energy2) / (0.100001 - 0.1) / 30.0

            @test abs(finitediff + force2[2,3]) < 1e-4

            finitediff = (energy2b - energy2) / (0.000001 - 0.0) / 30.0

            @test abs(finitediff + force2[2,1]) < 1e-4
            
            println("finitediff ", finitediff)
            println("force2")
            println(force2)
            println("stress")
            println(stress)
        end
    end
end

test1()
test3()
Nothing;

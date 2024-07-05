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

    @testset "testing classical elemental fit to force stress 2body" begin
        @suppress begin

            database = Dict()
            coef2 = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Al], 2)
            coef2.datH[:] = [1.0, 0.8, 0.6, 0.4, 0.2, 0.1]
            database[(:Al, :Al)] = coef2
            
            c2 = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.15], ["Al", "Al"])
            v = 1.0;
            C = []; EN = Float64[]; FORCES= Float64[]; STRESSES = Float64[]; EN2 = Float64[];
            for x in 1:12;
                v = v * 1.12;
                for c in [c2];
                    push!(C, c*v);
                    energy2,_ = ThreeBodyTB.Classical.calc_energy_cl(c*v, database=database, use_threebody=false, use_em = false, use_charges = false)
                    energy, force, stress = ThreeBodyTB.Classical.energy_force_stress_cl(c*v, database=database, use_threebody=false, use_em = false, use_charges = false)
                    push!(EN, energy)
                    push!(EN2, energy2)
                    FORCES=vcat(FORCES, force[:])
                    STRESSES = vcat(STRESSES, [stress[1,1], stress[1,2],stress[1,3],stress[2,2],stress[2,3],stress[3,3]])
                end;
            end
            V,Vf,Vs, _ = ThreeBodyTB.ClassicalFit.prepare_fit_cl(C, use_threebody=false, get_force=true, use_em = false, use_charges = false);
            x = V \ (EN / c2.nat)
            @test sum(abs.(x - [1.0, 0.8, 0.6, 0.4, 0.2, 0.1])) < 1e-8

            xf = Vf \ FORCES
            @test sum(abs.(xf - [1.0, 0.8, 0.6, 0.4, 0.2, 0.1])) < 1e-8

            xs = Vs \ STRESSES
            @test sum(abs.(xs - [1.0, 0.8, 0.6, 0.4, 0.2, 0.1])) < 1e-8
            
            xtot = [V;Vf;Vs] \ [EN/c2.nat; FORCES; STRESSES]
            @test sum(abs.(xtot - [1.0, 0.8, 0.6, 0.4, 0.2, 0.1])) < 1e-8
            
        end
    end
end

test1();

function test1a()

    @testset "testing classical elemental fit to force stress 2body charges" begin
        @suppress begin

#            database = Dict()
#            coef2 = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Al], 2)
#            coef2.datH[:] = [1.0, 0.8, 0.6, 0.4, 0.2, 0.1]
#            database[(:Al, :Al)] = coef2

            databaseC = Dict()
            coef2z = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Al], 2)
            coef2z.datH[:] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            databaseC[(:Al, :Al)] = coef2z
            coef2_charges = ThreeBodyTB.Classical.make_coefs_cl([:Al, :Al], 2, charges=true)
            coef2_charges.datH[:] = -1.0*[0.0, 1.0, 0.8, 0.6, 0.4, 0.2, 0.1]
            databaseC[(:Al, :Al, :charges)] = coef2_charges
            
            c2 = makecrys([30 0 0; 0 30 0; 0 0 30], [0 0 0; 0 0 0.15], ["Al", "Al"])
            v = 1.0;
            C = []; EN = Float64[]; FORCES= Float64[]; STRESSES = Float64[];
            for x in 1:12;
                v = v * 1.12;
                for c in [c2];
                    push!(C, c*v);
                    energy, force, stress = ThreeBodyTB.Classical.energy_force_stress_cl(c*v, database=databaseC, use_threebody=false, use_em=false, use_charges = true, fixed_charges = [1.0, -1.0])
                    push!(EN, energy / c.nat)
                    FORCES=vcat(FORCES, force[:])
                    STRESSES = vcat(STRESSES, [stress[1,1], stress[1,2],stress[1,3],stress[2,2],stress[2,3],stress[3,3]])
                end;
            end
            V,Vf,Vs, _ = ThreeBodyTB.ClassicalFit.prepare_fit_cl(C, use_threebody=false, get_force=true, use_em = false, use_charges= true, fixed_charges = [1.0, -1.0]);

            ref = [ 0.5,   0.4,   0.3,   0.2,   0.1,   0.05,   0.0,   -0.5,   -0.4,   -0.3,   -0.2,   -0.1,   -0.05]
            
            x = V \ EN
            @test sum(abs.(x - ref)) < 1e-8

            xf = Vf \ FORCES
            @test sum(abs.(xf - ref)) < 1e-8

            xs = Vs \ STRESSES
            @test sum(abs.(xs - ref)) < 1e-8
            
            xtot = [V;Vf;Vs] \ [EN; FORCES; STRESSES]
            @test sum(abs.(xtot - ref)) < 1e-8
            
        end
    end
end

test1a();

#test3();
#database_NEW2 = test3_2atoms();
#database_NEW3 = test3_3atoms();
Nothing;

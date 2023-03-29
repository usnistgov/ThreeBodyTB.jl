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
            C = []; EN = Float64[]; FORCES= Float64[]; STRESSES = Float64[];
            for x in 1:12;
                v = v * 1.12;
                for c in [c2];
                    push!(C, c*v);
                    energy, force, stress = ThreeBodyTB.Classical.energy_force_stress(c*v, database=database, use_threebody=false)
                    push!(EN, energy)
                    FORCES=vcat(FORCES, force[:])
                    STRESSES = vcat(STRESSES, [stress[1,1], stress[1,2],stress[1,3],stress[2,2],stress[2,3],stress[3,3]])
                end;
            end
            V,Vf,Vs, _ = ThreeBodyTB.Classical.prepare_fit_cl(C, use_threebody=false, get_force=true);
            x = V \ EN
            @test sum(abs.(x - [1.0, 0.8, 0.6, 0.4, 0.2, 0.1])) < 1e-8

            xf = Vf \ FORCES
            @test sum(abs.(xf - [1.0, 0.8, 0.6, 0.4, 0.2, 0.1])) < 1e-8

            xs = Vs \ STRESSES
            @test sum(abs.(xs - [1.0, 0.8, 0.6, 0.4, 0.2, 0.1])) < 1e-8
            
            xtot = [V;Vf;Vs] \ [EN; FORCES; STRESSES]
            @test sum(abs.(xtot - [1.0, 0.8, 0.6, 0.4, 0.2, 0.1])) < 1e-8
            
        end
    end
end

test1();

#test3();
#database_NEW2 = test3_2atoms();
#database_NEW3 = test3_3atoms();
Nothing;

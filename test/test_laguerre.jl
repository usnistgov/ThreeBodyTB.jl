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
    @testset "testing laguerre fitting fake example" begin

        @suppress begin
        
        c1 = makecrys([7.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0], ["Li"], units="Bohr");
        c2 = makecrys([8.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0], ["Li"], units="Bohr");
        c3 = makecrys([8.5 0 0; 0 8.5 0; 0 0 7.0], [0 0 0], ["Li"], units="Bohr");
        c4 = makecrys([6.0 0 0; 0 6.0 0; 0 0 6.0], [0 0 0], ["Li"], units="Bohr");

        c5 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.49 0 0 ], ["Li", "Li"], units="Bohr");
        c6 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.48 0 0 ], ["Li", "Li"], units="Bohr");
        c7 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.47 0 0 ], ["Li", "Li"], units="Bohr");
        c8 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.46 0 0 ], ["Li", "Li"], units="Bohr");
        c9 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.45 0 0 ], ["Li", "Li"], units="Bohr");

            c4a = makecrys([5.0 0 0; 0 5.0 0; 0 0 5.0], [0 0 0], ["Li"], units="Bohr");
            c4b = makecrys([4.0 0 0; 0 4.0 0; 0 0 4.0], [0 0 0], ["Li"], units="Bohr");
            c4c = makecrys([3.0 0 0; 0 3.0 0; 0 0 3.0], [0 0 0], ["Li"], units="Bohr");

            
        begin
            #        if true

            database = Dict()
            database[(:Li, :Li)] = ThreeBodyTB.CalcTB.make_coefs(Set(["Li", "Li"]), 2)
            
            tbc1 = ThreeBodyTB.CalcTB.calc_tb_fast(c1, database, use_threebody=false);
            tbc2 = ThreeBodyTB.CalcTB.calc_tb_fast(c2, database, use_threebody=false);
            tbc3 = ThreeBodyTB.CalcTB.calc_tb_fast(c3, database, use_threebody=false);
            tbc4 = ThreeBodyTB.CalcTB.calc_tb_fast(c4, database, use_threebody=false);
            tbc4a = ThreeBodyTB.CalcTB.calc_tb_fast(c4a, database, use_threebody=false);
            tbc4b = ThreeBodyTB.CalcTB.calc_tb_fast(c4b, database, use_threebody=false);
            tbc4c = ThreeBodyTB.CalcTB.calc_tb_fast(c4c, database, use_threebody=false);
            tbc5 = ThreeBodyTB.CalcTB.calc_tb_fast(c5, database, use_threebody=false);
            tbc6 = ThreeBodyTB.CalcTB.calc_tb_fast(c6, database, use_threebody=false);
            tbc7 = ThreeBodyTB.CalcTB.calc_tb_fast(c7, database, use_threebody=false);
            tbc8 = ThreeBodyTB.CalcTB.calc_tb_fast(c8, database, use_threebody=false);
            tbc9 = ThreeBodyTB.CalcTB.calc_tb_fast(c9, database, use_threebody=false);
            
            tbc_list = [tbc1, tbc2, tbc3, tbc4, tbc5, tbc6, tbc7, tbc8, tbc9, tbc4a, tbc4b, tbc4c]

            newdatabase = ThreeBodyTB.FitTB.do_fitting(tbc_list, fit_threebody=false, do_plot=false)
#            newdatabase = ThreeBodyTB.FitTB.do_fitting_recursive(tbc_list, fit_threebody=false, do_plot=false)


            #        println("newdartabase")
            #        println(newdatabase[("Li", "Li")])
            @test sum(abs.(newdatabase[(:Li, :Li)].datH .- database[(:Li, :Li)].datH)) ≤ 1.5e-3
        end
        
        end
    end
    return tbc_list
end


tbc_list = test1();
Nothing

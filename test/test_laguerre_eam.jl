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
    @testset "testing laguerre fitting fake example eam" begin

      @suppress  begin
        
        c1 = makecrys([7.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0], ["Li"], units="Bohr");
        c2 = makecrys([8.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0], ["Li"], units="Bohr");
        c3 = makecrys([8.5 0 0; 0 8.5 0; 0 0 7.0], [0 0 0], ["Li"], units="Bohr");
        c4 = makecrys([6.0 0 0; 0 6.0 0; 0 0 6.0], [0 0 0], ["Li"], units="Bohr");

        c5 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.49 0 0 ], ["Li", "Li"], units="Bohr");
        c6 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.45 0 0 ], ["Li", "Li"], units="Bohr");
        c7 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.40 0 0 ], ["Li", "Li"], units="Bohr");
        c8 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.35 0 0 ], ["Li", "Li"], units="Bohr");
        c9 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.47 0 0 ], ["Li", "Li"], units="Bohr");

            c4a = makecrys([5.0 0 0; 0 5.0 0; 0 0 5.0], [0 0 0], ["Li"], units="Bohr");
            c4b = makecrys([4.0 0 0; 0 4.0 0; 0 0 4.0], [0 0 0], ["Li"], units="Bohr");
            c4c = makecrys([3.0 0 0; 0 3.0 0; 0 0 3.0], [0 0 0], ["Li"], units="Bohr");

          cX1 = makecrys([100.0 0 0; 0 100.0 0; 0 0 100.0], [0 0 0; 0 0 5.0/100; 5.0/100 0 0], [:Li, :Li, :Li], units="Bohr");
          cX2 = makecrys([100.0 0 0; 0 100.0 0; 0 0 100.0], [0 0 0; 0 0 5.5/100; 5.0/100 0 0], [:Li, :Li, :Li], units="Bohr");
          cX3 = makecrys([100.0 0 0; 0 100.0 0; 0 0 100.0], [0 0 0; 0 0 5.5/100; 5.5/100 0 0], [:Li, :Li, :Li], units="Bohr");

          cX4 = makecrys([100.0 0 0; 0 100.0 0; 0 0 100.0], [0 0 0; 0 0 5.0/100], [:Li, :Li], units="Bohr");
          cX5 = makecrys([100.0 0 0; 0 100.0 0; 0 0 100.0], [0 0 0; 0 0 5.5/100], [:Li, :Li], units="Bohr");
          

            
        begin
            #        if true

            database = Dict()
            database[(:Li, :Li)] = ThreeBodyTB.CalcTB.make_coefs(Set(["Li", "Li"]), 2, use_eam=true)
            database[(:Li, :Li)].datS[1] = 0.0
            database[(:Li, :Li)].datS[2:end] .= 0.0
            database[(:Li, :Li)].datH .= 0.0
            database[(:Li, :Li)].datH[end-5:end] .= [1.0, 0.9, 0.8, 0.7, 0.6 ,0.5]
            
            tbc1 = ThreeBodyTB.CalcTB.calc_tb_LV(c1, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);
            tbc2 = ThreeBodyTB.CalcTB.calc_tb_LV(c2, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);
            tbc3 = ThreeBodyTB.CalcTB.calc_tb_LV(c3, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);
            tbc4 = ThreeBodyTB.CalcTB.calc_tb_LV(c4, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);
            tbc4a = ThreeBodyTB.CalcTB.calc_tb_LV(c4a, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);
            tbc4b = ThreeBodyTB.CalcTB.calc_tb_LV(c4b, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);
            tbc4c = ThreeBodyTB.CalcTB.calc_tb_LV(c4c, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);
            tbc5 = ThreeBodyTB.CalcTB.calc_tb_LV(c5, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);
            tbc6 = ThreeBodyTB.CalcTB.calc_tb_LV(c6, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);
            tbc7 = ThreeBodyTB.CalcTB.calc_tb_LV(c7, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);
            tbc8 = ThreeBodyTB.CalcTB.calc_tb_LV(c8, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);
            tbc9 = ThreeBodyTB.CalcTB.calc_tb_LV(c9, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);

            tbcX1 = ThreeBodyTB.CalcTB.calc_tb_LV(cX1, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);
            tbcX2 = ThreeBodyTB.CalcTB.calc_tb_LV(cX2, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);
            tbcX3 = ThreeBodyTB.CalcTB.calc_tb_LV(cX3, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);
            tbcX4 = ThreeBodyTB.CalcTB.calc_tb_LV(cX4, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);
            tbcX5 = ThreeBodyTB.CalcTB.calc_tb_LV(cX5, database, use_threebody=false, use_threebody_onsite=false, use_eam=true);
            
            
            tbc_list = [tbc1, tbc2, tbc3, tbc4, tbc5, tbc6, tbc7, tbc8, tbc9, tbc4a, tbc4b, tbc4c, tbcX1, tbcX2, tbcX3, tbcX4, tbcX5]

            newdatabase = ThreeBodyTB.FitTB.do_fitting(tbc_list, fit_threebody=false, fit_threebody_onsite=false, do_plot=false, fit_eam=true)

            #println(newdatabase[(:Li, :Li)].datH)
            #println(database[(:Li, :Li)].datH)
            println(sum(abs.(newdatabase[(:Li, :Li)].datH .- database[(:Li, :Li)].datH)))
            @test sum(abs.(newdatabase[(:Li, :Li)].datH .- database[(:Li, :Li)].datH)) ≤ 1.5e-3

        end
        
        end
    end
    return tbc_list
end


tbc_list = test1();
Nothing

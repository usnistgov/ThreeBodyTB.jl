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

function testEAM()
    tbc_list = []
    @testset "testing laguerre fitting fake example" begin

        atom = :H
        
      @suppress  begin
        
        c1 = makecrys([7.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0], [atom], units="Bohr");
        c2 = makecrys([8.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0], [atom], units="Bohr");
        c3 = makecrys([8.5 0 0; 0 8.5 0; 0 0 7.0], [0 0 0], [atom], units="Bohr");
        c4 = makecrys([6.0 0 0; 0 6.0 0; 0 0 6.0], [0 0 0], [atom], units="Bohr");

        c5 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.49 0 0 ], [atom, atom], units="Bohr");
        c6 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.45 0 0 ], [atom, atom], units="Bohr");
        c7 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.40 0 0 ], [atom, atom], units="Bohr");
        c8 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.35 0 0 ], [atom, atom], units="Bohr");
        c9 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.47 0 0 ], [atom, atom], units="Bohr");

          c4a = makecrys([5.0 0 0; 0 5.0 0; 0 0 5.0], [0 0 0], [atom], units="Bohr");
          c4b = makecrys([4.0 0 0; 0 4.0 0; 0 0 4.0], [0 0 0], [atom], units="Bohr");
          c4c = makecrys([3.0 0 0; 0 3.0 0; 0 0 3.0], [0 0 0], [atom], units="Bohr");

            
        begin
            #        if true

            database = Dict()
            database[(atom, atom)] = ThreeBodyTB.CalcTB.make_coefs(Set([atom, atom]), 2)
            database[(:eam, atom)] = ThreeBodyTB.CalcTB.make_coefs(Set([atom, atom]), 0)
#            database[(atom, atom)].datH .= 0.0
            
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
            
            tbc_list = [tbc1, tbc2, tbc3, tbc4, tbc5, tbc6, tbc7, tbc8, tbc9, tbc4a, tbc4b, tbc4c]

            newdatabase = ThreeBodyTB.FitTB.do_fitting(tbc_list, fit_threebody=false, fit_threebody_onsite=false,do_plot=false, fit_eam=true)
#            newdatabase = ThreeBodyTB.FitTB.do_fitting_recursive(tbc_list, fit_threebody=false, do_plot=false)


            #        println("newdartabase")
            #        println(newdatabase[("Li", "Li")])
#            println(newdatabase[(atom, atom)].datH)
#            println(database[(atom, atom)].datH)
#            println(sum(abs.(newdatabase[(atom, atom)].datH .- database[(atom, atom)].datH)))
            @test sum(abs.(newdatabase[(atom, atom)].datH .- database[(atom, atom)].datH)) ≤ 1.5e-3

            @test sum(abs.(newdatabase[(:eam, atom)].datH .- database[(:eam, atom)].datH)) ≤ 1.5e-3

        end
        
        end
    end
    return tbc_list
end


tbc_list = testEAM();
Nothing

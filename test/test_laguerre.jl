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

      @suppress  begin
          database = Dict()
          function prepare_data()
              c1 = makecrys([7.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0], ["Li"], units="Bohr");
              c2 = makecrys([8.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0], ["Li"], units="Bohr");
              c3 = makecrys([8.5 0 0; 0 8.5 0; 0 0 7.0], [0 0 0], ["Li"], units="Bohr");
              c4 = makecrys([6.0 0 0; 0 6.0 0; 0 0 6.0], [0 0 0], ["Li"], units="Bohr");

              c5 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.49 0 0 ], ["Li", "Li"], units="Bohr");
              c6 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.45 0 0 ], ["Li", "Li"], units="Bohr");
              c7 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.40 0 0 ], ["Li", "Li"], units="Bohr");
              c8 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.35 0 0 ], ["Li", "Li"], units="Bohr");
              c9 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.47 0 0 ], ["Li", "Li"], units="Bohr");

              c10 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.30 0 0 ], ["Li", "Li"], units="Bohr");
              c11 = makecrys([20.0 0 0; 0 8.0 0; 0 0 8.0], [0 0 0; 0.49 0 0 ], ["Li", "Li"], units="Bohr");

              
              c4a = makecrys([5.0 0 0; 0 5.0 0; 0 0 5.0], [0 0 0], ["Li"], units="Bohr");
              c4b = makecrys([5.5 0 0; 0 5.5 0; 0 0 5.5], [0 0 0], ["Li"], units="Bohr");
              #            c4c = makecrys([3.0 0 0; 0 3.0 0; 0 0 3.0], [0 0 0], ["Li"], units="Bohr");
              #            c4b = makecrys([8.0 0 0; 0 8.0 0; 0 0 8.0], [0 0 0], ["Li"], units="Bohr");
              c4c = makecrys([9.0 0 0; 0 9.0 0; 0 0 9.0], [0 0 0], ["Li"], units="Bohr");

              
              #        begin
              #        if true

#              database = Dict()
              database[(:Li, :Li)] = ThreeBodyTB.CalcTB.make_coefs(Set(["Li", "Li"]), 2)
              
              tbc1 = ThreeBodyTB.CalcTB.calc_tb_LV(c1, database, use_threebody=false, use_threebody_onsite=false);
              tbc2 = ThreeBodyTB.CalcTB.calc_tb_LV(c2, database, use_threebody=false, use_threebody_onsite=false);
              tbc3 = ThreeBodyTB.CalcTB.calc_tb_LV(c3, database, use_threebody=false, use_threebody_onsite=false);
              tbc4 = ThreeBodyTB.CalcTB.calc_tb_LV(c4, database, use_threebody=false, use_threebody_onsite=false);
              tbc4a = ThreeBodyTB.CalcTB.calc_tb_LV(c4a, database, use_threebody=false, use_threebody_onsite=false);
              tbc4b = ThreeBodyTB.CalcTB.calc_tb_LV(c4b, database, use_threebody=false, use_threebody_onsite=false);
              tbc4c = ThreeBodyTB.CalcTB.calc_tb_LV(c4c, database, use_threebody=false, use_threebody_onsite=false);
              tbc5 = ThreeBodyTB.CalcTB.calc_tb_LV(c5, database, use_threebody=false, use_threebody_onsite=false);
              tbc6 = ThreeBodyTB.CalcTB.calc_tb_LV(c6, database, use_threebody=false, use_threebody_onsite=false);
              tbc7 = ThreeBodyTB.CalcTB.calc_tb_LV(c7, database, use_threebody=false, use_threebody_onsite=false);
              tbc8 = ThreeBodyTB.CalcTB.calc_tb_LV(c8, database, use_threebody=false, use_threebody_onsite=false);
              tbc9 = ThreeBodyTB.CalcTB.calc_tb_LV(c9, database, use_threebody=false, use_threebody_onsite=false);

              tbc10 = ThreeBodyTB.CalcTB.calc_tb_LV(c10, database, use_threebody=false, use_threebody_onsite=false);
              tbc11 = ThreeBodyTB.CalcTB.calc_tb_LV(c11, database, use_threebody=false, use_threebody_onsite=false);

              tbc_list = [tbc1, tbc2, tbc3, tbc4, tbc5, tbc6, tbc7, tbc8, tbc9, tbc4a, tbc4b, tbc4c, tbc10, tbc11]

              return tbc_list
          end
          tbc_list = prepare_data()
            newdatabase = ThreeBodyTB.FitTB.do_fitting(tbc_list, fit_threebody=false, fit_threebody_onsite=false, do_plot=false)
#            newdatabase = ThreeBodyTB.FitTB.do_fitting_recursive(tbc_list, fit_threebody=false, do_plot=false)


            #        println("newdartabase")
            #        println(newdatabase[("Li", "Li")])
            println(newdatabase[(:Li, :Li)].datH)
            println(database[(:Li, :Li)].datH)
            println(sum(abs.(newdatabase[(:Li, :Li)].datH .- database[(:Li, :Li)].datH)))
            @test sum(abs.(newdatabase[(:Li, :Li)].datH .- database[(:Li, :Li)].datH)) ≤ 1.5e-3

          if false
              for t in tbc_list
                  en = 0.0
                  en1 = 0.0
                  @suppress begin
                      en, tbc, flag = scf_energy(t)
                      en1, tbc, flag1 = scf_energy(t.crys, database=newdatabase, repel=false)
                  end
                  println("en $en   en1 $en1  diff  $(en1 - en)")
              end
          end
          
#        end
#        
        end
    end
    return tbc_list
end


tbc_list = test1();

Nothing

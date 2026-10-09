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
        
        begin
            database = Dict()
            database["SCF"] = false
            database["scf"] = false
            
            c_temp = ThreeBodyTB.CalcTB.make_coefs(Set(["Li", "Li"]), 2)
            for k in keys(c_temp.inds)
                if k[end] == :O || k[end] == :A
                    c_temp.datH[ c_temp.inds[k]] .= 0.2
                end
            end
            
            datH_arr = [c_temp.datH c_temp.datH * 0.5 c_temp.datH * 0.1]
            
            database[(:Li, :Li)] = ThreeBodyTB.CalcTB.make_coefs(Set(["Li", "Li"]), 2, neighbor_dist = 6.0, n_bound = 2, bound_arr = [2.0 3.0; 100.0  101.0], datH_arr = datH_arr , use_arr=true)
            #database[(:Li, :Li)] = ThreeBodyTB.CalcTB.make_coefs(Set(["Li", "Li"]), 2, neighbor_dist = 6.0, n_bound = 1, bound_arr = [2.0 3.0], datH_arr = datH_arr[:,:] , use_arr=true)

            #database2 = Dict()
            #database2["SCF"] = false
            #database2["scf"] = false
            
            #c_temp2 = ThreeBodyTB.CalcTB.make_coefs(Set(["Li", "Li"]), 2)
            #for k in keys(c_temp2.inds)
            #    if k[end] == :O || k[end] == :A
            #        c_temp2.datH[ c_temp2.inds[k]] .= 0.2
            #    end
            #end
            
            #datH_arr = [c_temp2.datH c_temp2.datH * 0.5 c_temp2.datH * 0.1]
            
            #database2[(:Li, :Li)] = ThreeBodyTB.CalcTB.make_coefs(Set(["Li", "Li"]), 2, neighbor_dist = 6.0, n_bound = 2, bound_arr = [2.0 3.0; 100.0  101.0], datH_arr = datH_arr[:,:] , use_arr=true)
            #database[(:Li, :Li)] = ThreeBodyTB.CalcTB.make_coefs(Set(["Li", "Li"]), 2, neighbor_dist = 6.0, n_bound = 1, bound_arr = [2.0 3.0], datH_arr = datH_arr , use_arr=true)

            #function f()
            #database3 = Dict()
            #database3["SCF"] = false
            #database3["scf"] = false
           # 
            #c_temp3 = ThreeBodyTB.CalcTB.make_coefs(Set(["Li", "Li"]), 2)
            #for k in keys(c_temp3.inds)
            #    if k[end] == :O || k[end] == :A
            #        c_temp3.datH[ c_temp3.inds[k]] .= 0.2
            #    end
            #end
           # 
           # datH_arr = [c_temp3.datH c_temp3.datH * 0.5 c_temp3.datH * 0.1]
           # 
           # database3[(:Li, :Li)] = ThreeBodyTB.CalcTB.make_coefs(Set(["Li", "Li"]), 2, neighbor_dist = 6.0, n_bound = 2, bound_arr = [2.0 3.0; 5.5  6.5], datH_arr = datH_arr[:,:] , use_arr=true)
           # end
            
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
              
              tbc1 = ThreeBodyTB.CalcTB.calc_tb_LV(c1, database, use_threebody=false, use_threebody_onsite=false, use_arr= true);
              tbc2 = ThreeBodyTB.CalcTB.calc_tb_LV(c2, database, use_threebody=false, use_threebody_onsite=false, use_arr= true);
              tbc3 = ThreeBodyTB.CalcTB.calc_tb_LV(c3, database, use_threebody=false, use_threebody_onsite=false, use_arr= true);
              tbc4 = ThreeBodyTB.CalcTB.calc_tb_LV(c4, database, use_threebody=false, use_threebody_onsite=false, use_arr= true);
              tbc4a = ThreeBodyTB.CalcTB.calc_tb_LV(c4a, database, use_threebody=false, use_threebody_onsite=false, use_arr= true);
              tbc4b = ThreeBodyTB.CalcTB.calc_tb_LV(c4b, database, use_threebody=false, use_threebody_onsite=false, use_arr= true);
              tbc4c = ThreeBodyTB.CalcTB.calc_tb_LV(c4c, database, use_threebody=false, use_threebody_onsite=false, use_arr= true);
              tbc5 = ThreeBodyTB.CalcTB.calc_tb_LV(c5, database, use_threebody=false, use_threebody_onsite=false, use_arr= true);
              tbc6 = ThreeBodyTB.CalcTB.calc_tb_LV(c6, database, use_threebody=false, use_threebody_onsite=false, use_arr= true);
              tbc7 = ThreeBodyTB.CalcTB.calc_tb_LV(c7, database, use_threebody=false, use_threebody_onsite=false, use_arr= true);
              tbc8 = ThreeBodyTB.CalcTB.calc_tb_LV(c8, database, use_threebody=false, use_threebody_onsite=false, use_arr= true);
              tbc9 = ThreeBodyTB.CalcTB.calc_tb_LV(c9, database, use_threebody=false, use_threebody_onsite=false, use_arr= true);

              tbc10 = ThreeBodyTB.CalcTB.calc_tb_LV(c10, database, use_threebody=false, use_threebody_onsite=false, use_arr= true);
              tbc11 = ThreeBodyTB.CalcTB.calc_tb_LV(c11, database, use_threebody=false, use_threebody_onsite=false, use_arr= true);

              tbc_list = [tbc1, tbc2, tbc3, tbc4, tbc5, tbc6, tbc7, tbc8, tbc9, tbc4a, tbc4b, tbc4c, tbc10, tbc11]

              for v = 6.25:0.02:6.4
                  cX = makecrys([v 0 0; 0 v 0; 0 0 v], [0 0 0], ["Li"], units="Bohr");
                  tbcX = ThreeBodyTB.CalcTB.calc_tb_LV(cX, database, use_threebody=false, use_threebody_onsite=false, use_arr= true);
                  push!(tbc_list, tbcX)
              end
                  
              
              return tbc_list
          end
          tbc_list = prepare_data()
          #newdatabase = ThreeBodyTB.FitTB.do_fitting(tbc_list, fit_threebody=false, fit_threebody_onsite=false, do_plot=false)

            #ch = database[(:Li, :Li)].datH[:]
            #ch_arr = [ch ch]
            ch_arr =  database[(:Li, :Li)].datH_arr[:,:]
            #database_arr = Dict()

          
          
            newdatabase_direct, err = ThreeBodyTB.FitTB.do_fitting_direct(tbc_list, fit_threebody=false, fit_threebody_onsite=false, do_plot=false, database_arr = database, arr_ind = 2, ch_arr = ch_arr, use_arr=true )

          
            #        println("newdartabase")
            #        println(newdatabase[("Li", "Li")])
          #  println(newdatabase[(:Li, :Li)].datH)
          #  println(database[(:Li, :Li)].datH)
          #  println(sum(abs.(newdatabase[(:Li, :Li)].datH .- database[(:Li, :Li)].datH)))
          #  @test sum(abs.(newdatabase[(:Li, :Li)].datH .- database[(:Li, :Li)].datH)) ≤ 1.5e-3

            @test sum(abs.(newdatabase_direct[(:Li, :Li)].datH .- database[(:Li, :Li)].datH_arr[:,2])) ≤ 1.5e-3
          
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

using Test
using ThreeBodyTB
using Suppressor


function test1()
    tbc_list = []
    @testset "testing laguerre fitting fake example" begin

      @suppress  begin

          coef = ThreeBodyTB.CalcTB.make_coefs(Set([:H, :H]), 2, fillzeros=true)
          coefHLi = ThreeBodyTB.CalcTB.make_coefs(Set([:Li, :H]), 2, fillzeros=true)
          coefLi = ThreeBodyTB.CalcTB.make_coefs(Set([:Li, :Li]), 2, fillzeros=true)
          

          database_start = Dict()
          database_start[(:H, :Li)] = coefHLi
          database_start[(:Li, :H)] = coefHLi
          database_start[(:Li, :Li)] = coefLi
          
          
          database_trivial = Dict(); database_trivial["scf"] = true; database_trivial["SCF"] = true;
          c0 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5], [:H])
          en, tbc0, flag = scf_energy(c0, database=database_trivial, repel=false)
          
          x = 0.02;  c1 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:H, :H])
          x = 0.03; c2 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:H, :H])
          x = 0.04; c3 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:H, :H])
          x = 0.05; c4 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:H, :H])
          x = 0.06; c5 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:H, :H])
          x = 0.10; c6 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:H, :H])
          x = 0.12; c7 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:H, :H])
          x = 0.035; c8 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:H, :H])
          
          tbc_list = [tbc0]
#          kpoints = []

          for c in [c1,c2,c3,c4,c5,c6,c7,c8]
              en, tbc, flag = scf_energy(c, database=database_trivial, repel=false)
              push!(tbc_list, tbc)
#              push!(kpoints, zeros(1,3))
          end

          #          newdatabase = ThreeBodyTB.FitTB.do_fitting_recursive(tbc_list, fit_threebody=false, fit_threebody_onsite=false, do_plot=false, fit_umat = false)
          newdatabase = ThreeBodyTB.FitTB.do_fitting_recursive(tbc_list, fit_threebody=false, fit_threebody_onsite=false, do_plot=false, fit_umat = false, lambda = 1e-3, starting_database = database_start);   
          @test sum(abs.(newdatabase[(:H, :H)].datU )) < 1e-7

          begin
              x = 0.03; c1a = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x; 0.0 0.0 0.0], [:H, :H, :Li])
              x = 0.04; c3a = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x; 0.0 0.0 0.0], [:H, :H, :Li])
              x = 0.05; c5a = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x; 0.0 0.0 0.0], [:H, :H, :Li])
              x = 0.06; c7a = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x; 0.0 0.0 0.0], [:H, :H, :Li])
              x = 0.035; c8a = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x; 0.0 0.0 0.0], [:H, :H, :Li])
              x = 0.045; c9a = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x; 0.0 0.0 0.0], [:H, :H, :Li])          

              database_trivialU = Dict(); database_trivialU["scf"] = true; database_trivialU["SCF"] = true;

              coef = ThreeBodyTB.CalcTB.make_coefs(Set([:H, :H]), 2, fillzeros=true)
              coefHLi = ThreeBodyTB.CalcTB.make_coefs(Set([:Li, :H]), 2, fillzeros=true)
              coefLi = ThreeBodyTB.CalcTB.make_coefs(Set([:Li, :Li]), 2, fillzeros=true)

              coef.datU = [10.0, 10.0]
              #coef.datU = [0.0, 0.0]
              coefHLi.datU = [0.0, 0.0, 0.0, 0.0]
              coefLi.datU = [0.0, 0.0]

              database_trivialU[(:H, :H)] = coef
              database_trivialU[(:H, :Li)] = coefHLi
              database_trivialU[(:Li, :H)] = coefHLi
              database_trivialU[(:Li, :Li)] = coefLi

              database_start = Dict()
              database_start[(:H, :Li)] = coefHLi
              database_start[(:Li, :H)] = coefHLi
              database_start[(:Li, :Li)] = coefLi
              
              EN = []
              EN_OLD = []
              tbc_list_new = deepcopy(tbc_list)
              for c in [c1,c2,c3,c4,c5,c6,c7,c8,c1a,c3a,c5a,c7a,c8a,c9a]
                  tbc = ThreeBodyTB.CalcTB.calc_tb_LV(c, database_trivialU, use_threebody=false, use_threebody_onsite=false, repel = false, use_umat = true) 
                  en, tbc, flag = scf_energy(tbc , mix = 0.05, iters = 300, repel=false,   conv_thr = 1e-8)
                  push!(EN, en)
                  push!(tbc_list_new, tbc)
                  tbc = ThreeBodyTB.CalcTB.calc_tb_LV(c, newdatabase, use_threebody=false, use_threebody_onsite=false, repel = false, use_umat = false) 
                  en, tbc, flag = scf_energy(tbc , mix = 0.05, iters = 300, repel = false, conv_thr = 1e-8)
                  push!(EN_OLD, en)
              end
          end
          newdatabaseU = ThreeBodyTB.FitTB.do_fitting_recursive(deepcopy(tbc_list_new), fit_threebody=false, fit_threebody_onsite=false, do_plot=false, fit_umat = true, lambda = 1e-10, starting_database=database_start);          
          
          EN_TEST = []

          #newdatabaseU[(:H, :H)].datU = [10.0, 10.0, 10.0]
          for c in [c1,c2,c3,c4,c5,c6,c7,c8,c1a,c3a,c5a,c7a,c8a,c9a]
              tbc = ThreeBodyTB.CalcTB.calc_tb_LV(c, newdatabaseU, use_threebody=false, use_threebody_onsite=false, repel = false, use_umat = true) 
              en, tbc, flag = scf_energy(tbc , mix = 0.05, iters = 300, conv_thr = 1e-8)
              push!(EN_TEST, en)
          end
          [EN_OLD - EN EN_TEST - EN]
          
        
        end
    end
    return tbc_list
end


tbc_list = test1();
Nothing


function test2()
    tbc_list = []
    @testset "testing laguerre fitting fake example" begin

      @suppress  begin

          coef = ThreeBodyTB.CalcTB.make_coefs(Set([:H, :H]), 2, fillzeros=true)
          coefHLi = ThreeBodyTB.CalcTB.make_coefs(Set([:Li, :H]), 2, fillzeros=true)
          coefLi = ThreeBodyTB.CalcTB.make_coefs(Set([:Li, :Li]), 2, fillzeros=true)
          

          database_start = Dict()
          database_start[(:H, :Li)] = coefHLi
          database_start[(:Li, :H)] = coefHLi
          database_start[(:Li, :Li)] = coefLi
          
          
          database_trivial = Dict(); database_trivial["scf"] = true; database_trivial["SCF"] = true;
          c0 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5], [:H])
          en, tbc0, flag = scf_energy(c0, database=database_trivial, repel=false)
          
          x = 0.02;  c1 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:H, :H])
          x = 0.03; c2 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:H, :H])
          x = 0.04; c3 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:H, :H])
          x = 0.05; c4 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:H, :H])
          x = 0.06; c5 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:H, :H])
          x = 0.10; c6 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:H, :H])
          x = 0.12; c7 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:H, :H])
          x = 0.035; c8 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:H, :H])
          x = 0.035; c9 = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x], [:Li, :Li])

          
          tbc_list = [tbc0]
#          kpoints = []

          for c in [c1,c2,c3,c4,c5,c6,c7,c8,c9]
              en, tbc, flag = scf_energy(c, database=database_trivial, repel=false)
              push!(tbc_list, tbc)
#              push!(kpoints, zeros(1,3))
          end

          #          newdatabase = ThreeBodyTB.FitTB.do_fitting_recursive(tbc_list, fit_threebody=false, fit_threebody_onsite=false, do_plot=false, fit_umat = false)
          newdatabase = ThreeBodyTB.FitTB.do_fitting_recursive(tbc_list, fit_threebody=false, fit_threebody_onsite=false, do_plot=false, fit_umat = false, lambda = 1e-3, starting_database = database_start);   
          @test sum(abs.(newdatabase[(:H, :H)].datU )) < 1e-7

          begin
              x = 0.03; c1a = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x; 0.0 0.0 0.0], [:H, :H, :Li])
              x = 0.04; c3a = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x; 0.0 0.0 0.0], [:H, :H, :Li])
              x = 0.05; c5a = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x; 0.0 0.0 0.0], [:H, :H, :Li])
              x = 0.06; c7a = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x; 0.0 0.0 0.0], [:H, :H, :Li])
              x = 0.035; c8a = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x; 0.0 0.0 0.0], [:H, :H, :Li])
              x = 0.045; c9a = makecrys([100 0 0; 0 100 0; 0 0 100], [0.5 0.5 0.5-x; 0.5 0.5 0.5+x; 0.0 0.0 0.0], [:H, :H, :Li])          

              database_trivialU = Dict(); database_trivialU["scf"] = true; database_trivialU["SCF"] = true;

              coef = ThreeBodyTB.CalcTB.make_coefs(Set([:H, :H]), 2, fillzeros=true)
              coefHLi = ThreeBodyTB.CalcTB.make_coefs(Set([:Li, :H]), 2, fillzeros=true)
              coefLi = ThreeBodyTB.CalcTB.make_coefs(Set([:Li, :Li]), 2, fillzeros=true)

              coef.datU = [10.0, 10.0]
              #coef.datU = [0.0, 0.0]
              coefHLi.datU = [0.0, 0.0, 0.0, 0.0]
              coefLi.datU = [0.0, 0.0]

              database_trivialU[(:H, :H)] = coef
              database_trivialU[(:H, :Li)] = coefHLi
              database_trivialU[(:Li, :H)] = coefHLi
              database_trivialU[(:Li, :Li)] = coefLi

              database_start = Dict()
              database_start[(:H, :Li)] = coefHLi
              database_start[(:Li, :H)] = coefHLi
              database_start[(:H, :H)] = coef
              

              
              EN = []
              EN_OLD = []
              tbc_list_new = deepcopy(tbc_list)
              for c in [c1,c2,c3,c4,c5,c6,c7,c8,c9,c1a,c3a,c5a,c7a,c8a,c9a]
                  tbc = ThreeBodyTB.CalcTB.calc_tb_LV(c, database_trivialU, use_threebody=false, use_threebody_onsite=false, repel = false, use_umat = true) 
                  en, tbc, flag = scf_energy(tbc , mix = 0.05, iters = 300, repel=false,   conv_thr = 1e-8)
                  push!(EN, en)
                  push!(tbc_list_new, tbc)
                  tbc = ThreeBodyTB.CalcTB.calc_tb_LV(c, newdatabase, use_threebody=false, use_threebody_onsite=false, repel = false, use_umat = false) 
                  en, tbc, flag = scf_energy(tbc , mix = 0.05, iters = 300, repel = false, conv_thr = 1e-8)
                  push!(EN_OLD, en)
              end
          end
          newdatabaseU = ThreeBodyTB.FitTB.do_fitting_recursive(deepcopy(tbc_list_new), fit_threebody=false, fit_threebody_onsite=false, do_plot=false, fit_umat = true, lambda = 1e-7, starting_database=database_start);          
          
          EN_TEST = []

          #newdatabaseU[(:H, :H)].datU = [10.0, 10.0, 10.0]
          for c in [c1,c2,c3,c4,c5,c6,c7,c8,c9, c1a,c3a,c5a,c7a,c8a,c9a]
              tbc = ThreeBodyTB.CalcTB.calc_tb_LV(c, newdatabaseU, use_threebody=false, use_threebody_onsite=false, repel = false, use_umat = true) 
              en, tbc, flag = scf_energy(tbc , mix = 0.05, iters = 300, conv_thr = 1e-8)
              push!(EN_TEST, en)
          end
          [EN_OLD - EN EN_TEST - EN]
          
        
        end
    end
    return tbc_list
end


tbc_list = test2();
Nothing

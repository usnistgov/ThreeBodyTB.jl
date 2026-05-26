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

            database = Dict()
            database[(:H, :H)] = ThreeBodyTB.CalcTB.make_coefs(Set(["H", "H"]), 2)
            database[(:H, :H)].datH[:] .= 0.0            
            database[(:H, :H, :H)] = ThreeBodyTB.CalcTB.make_coefs(Set(["H", "H", "H"]), 3)

            database[(:H, :H, :H)].datH[:] .= 0.1
            #database[(:H, :H, :H)].datH[1] = 0.1
            function prepare_data()
                c1 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0], ["H"], units="Bohr");
                c2 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.49], ["H", "H"], units="Bohr");
                c3 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.45], ["H", "H"], units="Bohr");
                c4 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.4], ["H", "H"], units="Bohr");
                c5 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.35], ["H", "H"], units="Bohr");
                c6 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.3], ["H", "H"], units="Bohr");
                c7 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.25], ["H", "H"], units="Bohr");
                c8 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.2], ["H", "H"], units="Bohr");
                c9 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.15], ["H", "H"], units="Bohr");
                c10 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0]*0.8, [0 0 0], ["H"], units="Bohr");

                c20 = makecrys([5.0 0 0; 0 5.0 0; 0 0 5.0], [0 0 0], ["H"], units="Bohr");
                c21 = makecrys([5.0 0 0; 0 5.0 0; 0 0 5.0]*1.1, [0 0 0], ["H"], units="Bohr");
                c22 = makecrys([5.0 5.0 0; 5.0 0 5.0 ; 0 5.0 5.0]*0.9, [0 0 0], ["H"], units="Bohr");
                c23 = makecrys([5.0 5.0 0; 5.0 0 5.0 ; 0 5.0 5.0]*0.95, [0 0 0], ["H"], units="Bohr");
                c24 = makecrys([5.0 5.0 0; 5.0 0 5.0 ; 0 5.0 5.0]*0.85, [0 0 0], ["H"], units="Bohr");

                c25 = makecrys([5.0 5.0 -5; 5.0 -5 5.0 ; -5 5.0 5.0]*0.9, [0 0 0], ["H"], units="Bohr");
                c26 = makecrys([5.0 5.0 -5; 5.0 -5 5.0 ; -5 5.0 5.0]*0.95, [0 0 0], ["H"], units="Bohr");
                c27 = makecrys([5.0 5.0 -5; 5.0 -5 5.0 ; -5 5.0 5.0]*0.85, [0 0 0], ["H"], units="Bohr");

                c6a = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.3; 0 0 0.6], ["H", "H", "H"], units="Bohr");
                c6b = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.3; 0 0 0.65], ["H", "H", "H"], units="Bohr");
                c6c = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.3; 0 0 0.55], ["H", "H", "H"], units="Bohr");


#            begin
                #        if true

                
                tbc1 = ThreeBodyTB.CalcTB.calc_tb_LV(c1, database, use_threebody=true, use_threebody_onsite=true);
                tbc2 = ThreeBodyTB.CalcTB.calc_tb_LV(c2, database, use_threebody=true, use_threebody_onsite=true);
                tbc3 = ThreeBodyTB.CalcTB.calc_tb_LV(c3, database, use_threebody=true, use_threebody_onsite=true);
                tbc4 = ThreeBodyTB.CalcTB.calc_tb_LV(c4, database, use_threebody=true, use_threebody_onsite=true);
                tbc5 = ThreeBodyTB.CalcTB.calc_tb_LV(c5, database, use_threebody=true, use_threebody_onsite=true);
                tbc6 = ThreeBodyTB.CalcTB.calc_tb_LV(c6, database, use_threebody=true, use_threebody_onsite=true);
                tbc7 = ThreeBodyTB.CalcTB.calc_tb_LV(c7, database, use_threebody=true, use_threebody_onsite=true);
                tbc8 = ThreeBodyTB.CalcTB.calc_tb_LV(c8, database, use_threebody=true, use_threebody_onsite=true);
                tbc9 = ThreeBodyTB.CalcTB.calc_tb_LV(c9, database, use_threebody=true, use_threebody_onsite=true);
                tbc10 = ThreeBodyTB.CalcTB.calc_tb_LV(c10, database, use_threebody=true, use_threebody_onsite=true);

                tbc20 = ThreeBodyTB.CalcTB.calc_tb_LV(c20, database, use_threebody=true, use_threebody_onsite=true);
                tbc21 = ThreeBodyTB.CalcTB.calc_tb_LV(c21, database, use_threebody=true, use_threebody_onsite=true);
                tbc22 = ThreeBodyTB.CalcTB.calc_tb_LV(c22, database, use_threebody=true, use_threebody_onsite=true);
                tbc23 = ThreeBodyTB.CalcTB.calc_tb_LV(c23, database, use_threebody=true, use_threebody_onsite=true);
                tbc24 = ThreeBodyTB.CalcTB.calc_tb_LV(c24, database, use_threebody=true, use_threebody_onsite=true);
                tbc25 = ThreeBodyTB.CalcTB.calc_tb_LV(c25, database, use_threebody=true, use_threebody_onsite=true);
                tbc26 = ThreeBodyTB.CalcTB.calc_tb_LV(c26, database, use_threebody=true, use_threebody_onsite=true);
                tbc27 = ThreeBodyTB.CalcTB.calc_tb_LV(c27, database, use_threebody=true, use_threebody_onsite=true);


                tbc6a = ThreeBodyTB.CalcTB.calc_tb_LV(c6a, database, use_threebody=true, use_threebody_onsite=true);
                tbc6b = ThreeBodyTB.CalcTB.calc_tb_LV(c6b, database, use_threebody=true, use_threebody_onsite=true);
                tbc6c = ThreeBodyTB.CalcTB.calc_tb_LV(c6c, database, use_threebody=true, use_threebody_onsite=true);

                
                tbc_list = [tbc1, tbc2, tbc3, tbc4, tbc5, tbc6, tbc7, tbc8, tbc9, tbc10, tbc20, tbc21, tbc22, tbc23, tbc24, tbc25, tbc26, tbc27, tbc6a, tbc6b, tbc6c]
                return tbc_list
            end
            tbc_list = prepare_data()
                
#            kpts, kweights = ThreeBodyTB.TB.make_kgrid([4,4, 4])
            #KPOINTS, KWEIGHTS, nk_max = ThreeBodyTB.FitTB.get_k_simple(kpts, tbc_list)

            #newdatabase, ch, cs, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, HON, ind_BIG, KEYS, HIND, SIND, DMIN\
            #_TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, YS_new, cs , ch_refit, SPIN = ThreeBodyTB.FitTB.do_fitting_linear(tbc_list, fit_threebody=true, do_plot=true, mode=:kspace, kpoints=KPOINTS);

            newdatabase = ThreeBodyTB.FitTB.do_fitting(tbc_list, fit_threebody=true, fit_threebody_onsite=true, fit_eam=false, do_plot=false)


            #newdatabase = ThreeBodyTB.FitTB.do_fitting_recursive(tbc_list, fit_threebody=true, fit_threebody_onsite=true, do_plot=false, fit_eam=false)


                #        println("newdartabase")
            #        println(newdatabase[("Li", "Li")])
            @test sum(abs.(newdatabase[(:H, :H)].datH .- database[(:H, :H)].datH)) ≤ 1e-3
            @test sum(abs.(newdatabase[(:H, :H)].datS .- database[(:H, :H)].datS)) ≤ 1e-3
            @test sum(abs.(newdatabase[(:H, :H, :H)].datH .- database[(:H, :H, :H)].datH)) ≤ 1e-2
        end
        
        #        end
    end

end


tbc_list = test1();
Nothing;

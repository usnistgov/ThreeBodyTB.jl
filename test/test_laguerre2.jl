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

        begin
            #        if true

            database = Dict()
            database[(:H, :H)] = ThreeBodyTB.CalcTB.make_coefs(Set(["H", "H"]), 2)
            
            tbc1 = ThreeBodyTB.CalcTB.calc_tb_fast(c1, database, use_threebody=false);
            tbc2 = ThreeBodyTB.CalcTB.calc_tb_fast(c2, database, use_threebody=false);
            tbc3 = ThreeBodyTB.CalcTB.calc_tb_fast(c3, database, use_threebody=false);
            tbc4 = ThreeBodyTB.CalcTB.calc_tb_fast(c4, database, use_threebody=false);
            tbc5 = ThreeBodyTB.CalcTB.calc_tb_fast(c5, database, use_threebody=false);
            tbc6 = ThreeBodyTB.CalcTB.calc_tb_fast(c6, database, use_threebody=false);
            tbc7 = ThreeBodyTB.CalcTB.calc_tb_fast(c7, database, use_threebody=false);
            tbc8 = ThreeBodyTB.CalcTB.calc_tb_fast(c8, database, use_threebody=false);
            tbc9 = ThreeBodyTB.CalcTB.calc_tb_fast(c9, database, use_threebody=false);
            tbc10 = ThreeBodyTB.CalcTB.calc_tb_fast(c10, database, use_threebody=false);
            
            tbc_list = [tbc1, tbc2, tbc3, tbc4, tbc5, tbc6, tbc7, tbc8, tbc9, tbc10]

            kpts, kweights = ThreeBodyTB.TB.make_kgrid([2,2,2])
            KPOINTS, KWEIGHTS, nk_max = ThreeBodyTB.FitTB.get_k_simple(kpts, tbc_list)

            newdatabase, ch, cs, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, HON, ind_BIG, KEYS, HIND, SIND, DMIN\
            _TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, YS_new, cs , ch_refit, SPIN = ThreeBodyTB.FitTB.do_fitting_linear(tbc_list, fit_threebody=false, do_plot=true, mode=:kspace, kpoints=KPOINTS);
            #            newdatabase = ThreeBodyTB.FitTB.do_fitting(tbc_list, fit_threebody=false, do_plot=false)
#            newdatabase = ThreeBodyTB.FitTB.do_fitting_recursive(tbc_list, fit_threebody=false, do_plot=false)


            #        println("newdartabase")
            #        println(newdatabase[("Li", "Li")])
            @test sum(abs.(newdatabase[(:H, :H)].datH .- database[(:H, :H)].datH)) ≤ 1e-3
        end
        
        end
    end
    return tbc_list
end


tbc_list = test1();
Nothing;

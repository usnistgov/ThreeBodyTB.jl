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
    
    

    @testset "testing laguerre fitting fake example recursive with spin" begin
    @suppress     begin
        units_old = ThreeBodyTB.set_units()

        begin 
            ThreeBodyTB.set_units(both="atomic")

            c1 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0]*1.5, [0 0 0], ["Hx"], units="Bohr");
            c2 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.49], ["Hx", "Hx"], units="Bohr");
            c3 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.45], ["Hx", "Hx"], units="Bohr");
            c4 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.4], ["Hx", "Hx"], units="Bohr");
            c5 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.35], ["Hx", "Hx"], units="Bohr");
            c6 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.3], ["Hx", "Hx"], units="Bohr");
            c7 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.25], ["Hx", "Hx"], units="Bohr");
            c8 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.2], ["Hx", "Hx"], units="Bohr");
            c9 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0;0 0 0.15], ["Hx", "Hx"], units="Bohr");
            c10 = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0]*1.6, [0 0 0], ["Hx"], units="Bohr");
        end

        
        begin
            #        if true

            database = Dict()
            database["SCF"] = true
            database["scf"] = true
            database[(:Hx, :Hx)] = ThreeBodyTB.CalcTB.make_coefs(Set(["Hx", "Hx"]), 2)
            
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

        end

        tbc_list = []
        energy_list = []
        for tbc in [tbc1, tbc2, tbc3, tbc4, tbc5, tbc6, tbc7, tbc8, tbc9, tbc10]
            en, tbcX, flag = scf_energy(tbc.crys, database=database, nspin=1, conv_thr=1e-7, grid = [2,2,2], mix = 0.3)
            if flag == false
                en, tbcX, flag = scf_energy(tbc.crys, database=database, nspin=1, conv_thr=1e-7, grid = [2,2,2], mix = 0.7)
            end
            push!(energy_list, en)

            #
            #                    
            if true
                H = tbcX.tb.H
                H2 = zeros(2, size(H)[2], size(H)[3], size(H)[4])
                H2[1,:,:,:] = H
                H2[2,:,:,:] = H
                #            tbcX.nspin = 2
                #            tbcX.tb.nspin = 2
                #            tbcX.tb.H = H2
                #            tbcX.tb.scfspin = true
                tbt = ThreeBodyTB.TB.make_tb(H2, tbcX.tb.ind_arr, tbcX.tb.S)
                tbcX = ThreeBodyTB.TB.make_tb_crys(tbt, tbcX.crys, tbcX.nelec, tbcX.dftenergy, scf=tbcX.scf, eden = tbcX.eden, tb_energy=tbcX.energy[1])
            end
            push!(tbc_list, tbcX)
        end
        

            
        kpts, kweights = ThreeBodyTB.TB.make_kgrid([2,2,2])
        KPOINTS, KWEIGHTS, nk_max = ThreeBodyTB.FitTB.get_k_simple(kpts, tbc_list)
        
        newdatabase, ch, cs, X_Hnew_BIG, Xc_Hnew_BIG, Xc_Snew_BIG, X_H, X_Snew_BIG, Y_H, Y_S, HON, ind_BIG, KEYS, HIND, SIND, DMIN\
        _TYPES, DMIN_TYPES3, keepind, keepdata, Y_Hnew_BIG, Y_Snew_BIG, YS_new, cs , ch_refit, SPIN = ThreeBodyTB.FitTB.do_fitting_linear(tbc_list, fit_threebody=false, do_plot=false, mode=:kspace, kpoints=KPOINTS);
        
        @test sum(abs.(newdatabase[(:Hx, :Hx)].datH .- database[(:Hx, :Hx)].datH)) ≤ 2e-3

        newdatabase_rec = ThreeBodyTB.FitTB.do_fitting_recursive(tbc_list, fit_threebody=false, do_plot=false, kpoints=kpts);
        
        @test sum(abs.(newdatabase_rec[(:Hx, :Hx)].datH .- database[(:Hx, :Hx)].datH)) ≤ 0.1

        

        
        EDIFF = []


        #                for c in [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10, c11]
        for c in [c1,c3,c5,c7,c9]
            en1, tbc1, flag1 = scf_energy(c, database=database, nspin=2, conv_thr=1e-7, mix = 0.7, repel=false)
            if flag1 == false
                en1, tbc1, flag1 = scf_energy(c, database=database, nspin=2, conv_thr=1e-7, mix = 0.3, repel=false)
            end
            en2, tbc2, flag2 = scf_energy(c, database=newdatabase, nspin=2, conv_thr=1e-7, mix = 0.7, repel=false)
            if flag2 == false
                en2, tbc2, flag2 = scf_energy(c, database=newdatabase, nspin=2, conv_thr=1e-7, mix = 0.3, repel=false)
            end                    
            @test abs(en1 - en2 ) < 0.01
            push!(EDIFF, abs(en1 - en2 ))
        end
        @test sum(abs.(newdatabase[(:Hx, :Hx)].datH - database[(:Hx, :Hx)].datH)) < 1e-2
        println("EDIFF ", EDIFF)
        ThreeBodyTB.set_units(energy = units_old[1], length=units_old[2])
        
    end
    end    
    
    
end

#=
#            c1 = makecrys([5.0 0 0; 0 5.0 0; 0 0 5.0], [0 0 0], ["Hx"], units="Bohr");
#            c2 = makecrys([8.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0], ["Hx"], units="Bohr");
#            c3 = makecrys([8.5 0 0; 0 8.5 0; 0 0 7.0], [0 0 0], ["Hx"], units="Bohr");
#            c4 = makecrys([3.0 0 0; 0 3.0 0; 0 0 3.0], [0 0 0], ["Hx"], units="Bohr");
#
#
#            c5 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.49 0 0 ], ["Hx", "Hx"], units="Bohr");
#            c6 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.4 0 0 ], ["Hx", "Hx"], units="Bohr");
#            c7 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.35 0 0 ], ["Hx", "Hx"], units="Bohr");
#            c8 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.30 0 0 ], ["Hx", "Hx"], units="Bohr");
#            c9 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.25 0 0 ], ["Hx", "Hx"], units="Bohr");
#
#            c10 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.20 0 0 ], ["Hx", "Hx"], units="Bohr");
#            c11 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.15 0 0 ], ["Hx", "Hx"], units="Bohr");
#            
            #        c10 = c1 * 1.3
            #        c11 = c1 * 1.4

            
#            begin
#                #        if true
#
#                database = Dict()
#
#                database["SCF"] = true
#                database["scf"] = true
#                
#                database[(:Hx, :Hx)] = ThreeBodyTB.CalcTB.make_coefs(Set(["Hx", "Hx"]), 2)
#                
##                tbc1 = ThreeBodyTB.CalcTB.calc_tb_fast(c1, database, use_threebody=false);
##                tbc2 = ThreeBodyTB.CalcTB.calc_tb_fast(c2, database, use_threebody=false);
##                tbc3 = ThreeBodyTB.CalcTB.calc_tb_fast(c3, database, use_threebody=false);
##                tbc4 = ThreeBodyTB.CalcTB.calc_tb_fast(c4, database, use_threebody=false);
##                tbc5 = ThreeBodyTB.CalcTB.calc_tb_fast(c5, database, use_threebody=false);
##                tbc6 = ThreeBodyTB.CalcTB.calc_tb_fast(c6, database, use_threebody=false);
##                tbc7 = ThreeBodyTB.CalcTB.calc_tb_fast(c7, database, use_threebody=false);
##                tbc8 = ThreeBodyTB.CalcTB.calc_tb_fast(c8, database, use_threebody=false);
##                tbc9 = ThreeBodyTB.CalcTB.calc_tb_fast(c9, database, use_threebody=false);
##                tbc10 = ThreeBodyTB.CalcTB.calc_tb_fast(c10, database, use_threebody=false);
##                tbc11 = ThreeBodyTB.CalcTB.calc_tb_fast(c11, database, use_threebody=false);
#
#
#                
#                tbc_list = []
#
#                energy_list = []
#                for tbc in [tbc1, tbc2, tbc3, tbc4, tbc5, tbc6, tbc7, tbc8, tbc9, tbc10, tbc11]
#                    #            for tbc in [tbc1, tbc2, tbc3, tbc10]
#                    en, tbcX, flag = scf_energy(tbc.crys, database=database, nspin=1, conv_thr=1e-7, grid = [2,2,2], mix = 0.3)
#                    if flag == false
#                        en, tbcX, flag = scf_energy(tbc.crys, database=database, nspin=1, conv_thr=1e-7, grid = [2,2,2], mix = 0.7)
#                    end
#                    push!(energy_list, en)
#                    #                println("tbcX ", tbcX.nspin, " " , size(tbcX.eden))
#
#                    
#                    H = tbcX.tb.H
#                    H2 = zeros(2, size(H)[2], size(H)[3], size(H)[4])
#                    H2[1,:,:,:] = H
#                    H2[2,:,:,:] = H
#                    tbcX.nspin = 2
#                    tbcX.tb.nspin = 2
#                    tbcX.tb.H = H2
#                    tbcX.tb.scfspin = true
#
#                    push!(tbc_list, tbcX)
#                end
#
#                
#                #            tbc_list = [tbc1, tbc2, tbc3, tbc4, tbc5, tbc6, tbc7, tbc8, tbc9, tbc10, tbc11]
#
#                kpts, kweights = ThreeBodyTB.TB.make_kgrid([2,2,2])
#                
#                #            newdatabase = ThreeBodyTB.FitTB.do_fitting(tbc_list, fit_threebody=false, do_plot=false)
#                newdatabase = ThreeBodyTB.FitTB.do_fitting_recursive(tbc_list, fit_threebody=false, do_plot=false, niters=20, kpoints=kpts)
#
#
                #        println("newdartabase")
                #        println(newdatabase[("Li", "Li")])
                #            @test sum(abs.(newdatabase[(:Hx, :Hx)].datH .- database[(:Hx, :Hx)].datH)) ≤ 1e-5




    return tbc_list
end


=#

tbc_list = test1();
Nothing;

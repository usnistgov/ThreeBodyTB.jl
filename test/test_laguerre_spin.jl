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
        @suppress begin 
            units_old = ThreeBodyTB.set_units()
            ThreeBodyTB.set_units(both="atomic")
            
            c1 = makecrys([7.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0], ["Hx"], units="Bohr");
            c2 = makecrys([8.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0], ["Hx"], units="Bohr");
            c3 = makecrys([8.5 0 0; 0 8.5 0; 0 0 7.0], [0 0 0], ["Hx"], units="Bohr");
            c4 = makecrys([6.0 0 0; 0 6.0 0; 0 0 6.0], [0 0 0], ["Hx"], units="Bohr");

            c5 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.49 0 0 ], ["Hx", "Hx"], units="Bohr");
            c6 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.4 0 0 ], ["Hx", "Hx"], units="Bohr");
            c7 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.35 0 0 ], ["Hx", "Hx"], units="Bohr");
            c8 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.30 0 0 ], ["Hx", "Hx"], units="Bohr");
            c9 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.25 0 0 ], ["Hx", "Hx"], units="Bohr");

            c10 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.20 0 0 ], ["Hx", "Hx"], units="Bohr");
            c11 = makecrys([14.0 0 0; 0 7.0 0; 0 0 7.0], [0 0 0; 0.15 0 0 ], ["Hx", "Hx"], units="Bohr");
            
            #        c10 = c1 * 1.3
            #        c11 = c1 * 1.4

            
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
                tbc11 = ThreeBodyTB.CalcTB.calc_tb_fast(c11, database, use_threebody=false);


                
                tbc_list = []

                energy_list = []
                for tbc in [tbc1, tbc2, tbc3, tbc4, tbc5, tbc6, tbc7, tbc8, tbc9, tbc10, tbc11]
                    #            for tbc in [tbc1, tbc2, tbc3, tbc10]
                    en, tbcX, flag = scf_energy(tbc.crys, database=database, nspin=2, conv_thr=1e-7, grid = [4,4,4], mix = 0.3)
                    if flag == false
                        en, tbc, flag = scf_energy(tbc.crys, database=database, nspin=2, conv_thr=1e-7, mix = 0.7)
                    end
                    push!(energy_list, en)
                    #                println("tbcX ", tbcX.nspin, " " , size(tbcX.eden))
                    H = tbcX.tb.H
                    H2 = zeros(2, size(H)[2], size(H)[3], size(H)[4])
                    H2[1,:,:,:] = H
                    H2[2,:,:,:] = H
                    tbcX.nspin = 2
                    tbcX.tb.nspin = 2
                    tbcX.tb.H = H2
                    tbcX.tb.scfspin = true
                    push!(tbc_list, tbcX)
                end

                
                #            tbc_list = [tbc1, tbc2, tbc3, tbc4, tbc5, tbc6, tbc7, tbc8, tbc9, tbc10, tbc11]

                kpts, kweights = ThreeBodyTB.TB.make_kgrid([4,4,4])
                
                #            newdatabase = ThreeBodyTB.FitTB.do_fitting(tbc_list, fit_threebody=false, do_plot=false)
                newdatabase = ThreeBodyTB.FitTB.do_fitting_recursive(tbc_list, fit_threebody=false, do_plot=false, niters=20, kpoints=kpts)


                #        println("newdartabase")
                #        println(newdatabase[("Li", "Li")])
                #            @test sum(abs.(newdatabase[(:Hx, :Hx)].datH .- database[(:Hx, :Hx)].datH)) ≤ 1e-5
                EDIFF = []
                for c in [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10, c11]
                    en1, tbc1, flag1 = scf_energy(c, database=database, nspin=2, conv_thr=1e-7, mix = 0.7)
                    if flag1 == false
                        en1, tbc1, flag1 = scf_energy(c, database=database, nspin=2, conv_thr=1e-7, mix = 0.3)
                    end
                    en2, tbc2, flag2 = scf_energy(c, database=newdatabase, nspin=2, conv_thr=1e-7, mix = 0.7)
                    if flag2 == false
                        en2, tbc2, flag2 = scf_energy(c, database=newdatabase, nspin=2, conv_thr=1e-7, mix = 0.3)
                    end                    
                    @test abs(en1 - en2 ) < 0.01
                    push!(EDIFF, abs(en1 - en2 ))
                end
                @test sum(abs.(newdatabase[(:Hx, :Hx)].datH - database[(:Hx, :Hx)].datH)) < 1e-2
                println("EDIFF ", EDIFF)
            end
            
            ThreeBodyTB.set_units(energy = units_old[1], length=units_old[2])
        end
    end
    return tbc_list
end


tbc_list = test1();

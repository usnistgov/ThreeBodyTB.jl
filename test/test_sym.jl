

using Test
using ThreeBodyTB
using Suppressor


TESTDIR=ThreeBodyTB.TESTDIR


function test_Symmetry()

    @testset "test Symmetry" begin
        @suppress begin 

            c_hex = makecrys([1.0 0 0; -0.5 sqrt(3)/2 0 ; 0 0 1.4] * 3, [0 0 0], ["Al"])
            c_hex2 = makecrys([0.0 1.4 0; -0.5 0 sqrt(3)/2  ; 1.0 0 0 ] * 3, [0 0 0], ["Al"])

            c_std = ThreeBodyTB.Symmetry.get_standard_crys(c_hex, sym_prec = 1e-4);
            c_std2 = ThreeBodyTB.Symmetry.get_standard_crys(c_hex2, sym_prec = 1e-4);

            @test sum(abs.(c_std.A - c_std2.A)) < 1e-5

            s1, dat1 = ThreeBodyTB.Symmetry.get_symmetry(c_hex)
            s2, dat2 = ThreeBodyTB.Symmetry.get_symmetry(c_hex2)

            @test s1 == 191
            @test s2 == 191

            kpath, names, c_stdX = ThreeBodyTB.Symmetry.get_kpath_sym(c_hex)

            @test names[1] == "Γ"

            c = makecrys([6 0 0; 0 6 0;0 0 6], [0.0 0 0.0], ["H"])

            nk, grid_ind, kpts, kweights = ThreeBodyTB.Symmetry.get_kgrid_sym(c, grid=[6,6,6])

            @test nk == 20
            @test isapprox(sum(kweights) , 2.0, atol=1e-10)

            v,t = ThreeBodyTB.Symmetry.symmetrize_vector_tensor(rand(1,3), [1.0 0 0; 0 1.0 0; 0 0 1.0], c; sym_prec = 5e-4)

            @test isapprox(sum(abs.(v)) , 0.0, atol=1e-10)
            @test isapprox(sum(t) , 3.0, atol=1e-10)

            en, tbc, flag = scf_energy(c, use_sym = false)
            enS, tbcS, flagS = scf_energy(c, use_sym = true)

            @test isapprox(enS - en, 0.0, atol=1e-6)
            
        end
    end
end


test_Symmetry()

function test_Symmetry_Plot()

    @testset "test Symmetry Plot" begin
        @suppress begin 

            c_hex = makecrys([1.0 0 0; -0.5 sqrt(3)/2 0 ; 0 0 1.4] * 5, [0 0 0], ["Al"], units=:Bohr)
            plot_bandstr_sym(c_hex, do_display=false);
            plot_bandstr_dos(c_hex, do_display=false);

            @test 3 == 3
            

        end
    end
end

test_Symmetry_Plot()  

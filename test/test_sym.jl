

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

            s1 = ThreeBodyTB.Symmetry.get_symmetry(c_hex)
            s2 = ThreeBodyTB.Symmetry.get_symmetry(c_hex2)

            @test s1 == 191
            @test s2 == 191

            kpath, names, c_stdX = ThreeBodyTB.Symmetry.get_kpath_sym(c_hex)

            @test names[1] == "Γ"
            

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

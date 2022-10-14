using Test
using ThreeBodyTB
using Suppressor


TESTDIR=ThreeBodyTB.TESTDIR


function test_wan()

    @testset "test wan read write" begin
        @suppress begin 
            types=["Li"]
            pos=zeros((1,3))
            A=[ [10.0 0 0]; [0 10.0 0]; [ 0 0 10.0]]
            c=makecrys(A, pos, types, units="Bohr")
            en, tbc, flag = scf_energy(c)

            ThreeBodyTB.TB.write_hr_dat(tbc, directory=TESTDIR)
            
            tb = ThreeBodyTB.TB.load_hr_dat("wannier90_hr.dat", directory=TESTDIR)

            @test rm("$TESTDIR/wannier90_hr.dat") == nothing
        end
    end
end


test_wan()


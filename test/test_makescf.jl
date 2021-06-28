using ThreeBodyTB
using Test
using Suppressor

@testset "testing makescf" begin
    
    @suppress begin 
        types=["Li"]
        pos=zeros((1,3))
        A=[ [10.0 0 0]; [0 10.0 0]; [ 0 0 10.0]]
        c=makecrys(A, pos, types, units="Bohr")

        tmpdir, prefix, inputfile_str = ThreeBodyTB.QE.makeSCF(c)

        @test occursin("Li", inputfile_str)

    end
end

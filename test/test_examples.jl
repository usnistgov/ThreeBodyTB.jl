using ThreeBodyTB
using Test
using Suppressor

EXAMPLESDIR=ThreeBodyTB.EXAMPLESDIR

disp_old = deepcopy(ThreeBodyTB.BandStruct.no_display)
ThreeBodyTB.set_no_display(true)

for f in readdir("$EXAMPLESDIR")
    if occursin(".jl", f) && !occursin("~", f)
        @testset "example $f" begin
            @suppress begin 
                include("$EXAMPLESDIR/$f")
                @test 1 == 1
            end
        end        
    end
end

ThreeBodyTB.set_no_display(disp_old)


@testset "rm example pdfs generated" begin

    rm("mgs_compare_dft.pdf")
    rm("mgs_compare_self.pdf")
    rm("mgs_proj.pdf")
    rm("ScP_rocksalt.pdf")
    rm("al_fcc.pdf")
    rm("ScP_rocksalt_DOS.pdf")
##    rm("ScP_rocksalt_DOS_orbital.pdf")
    @test 1 == 1
end

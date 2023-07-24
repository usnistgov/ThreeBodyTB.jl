using ThreeBodyTB
using Test
using Suppressor


SRCDIR=ThreeBodyTB.SRCDIR

@testset "testing prototypes" begin
    @suppress begin 
        include("$SRCDIR/../code_for_dataset_gen/Prototypes.jl")
        t =  oxidation_guess("Na", "Cl")
        tref = [["Na", "Cl", :core_binary],["Na", "Cl", :A1B1], ["Na", "Cl", "pd3te"]]

        @test t == tref

        t =  oxidation_guess("Na", "Mg")
        #tref = [["Na", "Mg", :core_binary],["Na", "Mg", :metals], ["Na", "Mg", "hh_oxygen"], ["Na", "Mg", "hh_oxygen2"], ["Na", "Mg", "hex_oxygen"]]
        tref = [["Na", "Mg", :core_binary],["Na", "Mg", :metals]]
        
        @test t == tref

        #ThreeBodyTB.set_units(both="eVAng")

    end
end



using Test
using ThreeBodyTB
using Suppressor
using LinearAlgebra

TESTDIR=ThreeBodyTB.TESTDIR

function test_sparse()

    @testset "test plotting SPARSE" begin
        
      @suppress  begin
#        if true

          units_old = ThreeBodyTB.set_units()
          ThreeBodyTB.set_units(both="atomic")

          
          types=["Al"];
          
          #positions, crystal units
          pos=zeros((1,3));
          
          #lattice vectors, in Bohr units currently
          A=[ [3.7962751610 3.7962751610 0]; [3.7962751610 0 3.7962751610 ]; [ 0 3.7962751610 3.7962751610]];
          
          #makes the crystal
          c=makecrys(A, pos, types, units="Bohr")
          c = ThreeBodyTB.Symmetry.get_standard_crys(c)
          
          en_s, tbc_s, flag_s = scf_energy(c, sparse=true)
          @test typeof(tbc_s) <: ThreeBodyTB.TB.tb_crys_sparse
          @test typeof(tbc_s) <: ThreeBodyTB.TB.tb_crys
          @test flag_s == true
          
          plot_bandstr(tbc_s)
          plot_bandstr_dos(tbc_s)
          plot_bandstr_sym(tbc_s)
          dos(tbc_s)

          @test 1 == 1 #we at least check to see if there are errors running the codes.
          
        end
    end
end


test_sparse();
Nothing;

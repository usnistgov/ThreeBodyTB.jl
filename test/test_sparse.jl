using Test
using ThreeBodyTB
using Suppressor


TESTDIR=ThreeBodyTB.TESTDIR

function test_sparse()

    @testset "test SCF SPARSE" begin
        
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

          energy, tbc, flag = scf_energy(c)

          tbc_sparse = ThreeBodyTB.CalcTB.calc_tb_LV_sparse(c, ThreeBodyTB.ManageDatabase.database_cached)
          en_sparse, ret = ThreeBodyTB.SCF.scf_energy(tbc_sparse)
          
          @test abs(energy - -0.28707042) < 1e-2
          @test abs(en_sparse - energy) < 1e-5

          c2 = makecrys([5 0 0; 0 5 0; 0 0 10.0], [0 0 0; 0 0 0.51], [:Al, :P])

          energy2, tbc, flag = scf_energy(c2)

          tbc_sparse2 = ThreeBodyTB.CalcTB.calc_tb_LV_sparse(c2, ThreeBodyTB.ManageDatabase.database_cached)
          en_sparse2, ret2 = ThreeBodyTB.SCF.scf_energy(tbc_sparse2)

          @test abs(en_sparse2 - energy2) < 1e-5
          
        end
    end
end


test_sparse();
Nothing;

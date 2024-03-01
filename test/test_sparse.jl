using Test
using ThreeBodyTB
using Suppressor


TESTDIR=ThreeBodyTB.TESTDIR

function test_sparse()

    #sparsity should be basically seamless to users of code, as long as they don't look at the struct internal data representation.
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
          @test typeof(tbc) <: ThreeBodyTB.TB.tb_crys_dense

          en_sparse, tbc_sparse, flag_sparse = scf_energy(c, sparse=true)          
#          tbc_sparse = ThreeBodyTB.CalcTB.calc_tb_LV_sparse(c, ThreeBodyTB.ManageDatabase.database_cached)
#          en_sparse, ret = ThreeBodyTB.SCF.scf_energy(tbc_sparse)
          
          @test abs(energy - -0.28707042) < 1e-2
          @test abs(en_sparse - energy) < 1e-5

          a,b,c = Hk(tbc_sparse, [0.0,0.0,0.5])
          aS, bS, cS = Hk(tbc, [0.0,0.0,0.5])

          @test sum(abs.(b - bS)) < 1e-3

          #---------------
          c2 = makecrys([5 0 0; 0 5 0; 0 0 10.0], [0 0 0; 0 0 0.51], [:Al, :P])

          energy2, tbc, flag = scf_energy(c2)

          tbc_sparse2 = ThreeBodyTB.CalcTB.calc_tb_LV_sparse(c2, ThreeBodyTB.ManageDatabase.database_cached)
          en_sparse2, ret2 = ThreeBodyTB.SCF.scf_energy(tbc_sparse2)

          @test abs(en_sparse2 - energy2) < 1e-5


          #------------- check if auto set to sparse for nat > 100
          
          c = makecrys([10.0 0 0; 0 10.0 0; 0 0 10.0], [0 0 0 ], [:H], units="Bohr");
          en, tbc_auto, flag = scf_energy(c*[5,5,5])

          @test typeof(tbc_auto) <: ThreeBodyTB.TB.tb_crys_sparse
          
        end
    end
end


test_sparse();
Nothing;

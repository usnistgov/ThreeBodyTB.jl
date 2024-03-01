using Test
using ThreeBodyTB
using Suppressor
using LinearAlgebra

TESTDIR=ThreeBodyTB.TESTDIR

function test_sparse()

    @testset "test relax SPARSE" begin
        
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
          en, force, stress, tbc = scf_energy_force_stress(c)
          en_s, force_s, stress_s, tbc_s = scf_energy_force_stress(c, sparse=true)

          @test typeof(tbc_s) <: ThreeBodyTB.TB.tb_crys_sparse
          @test typeof(tbc) <: ThreeBodyTB.TB.tb_crys_dense
          
          @test abs(en_s - en) < 1e-5
          @test sum(abs.(force - force_s)) < 1e-5
          @test sum(abs.(stress - stress_s)) < 1e-5
          
          ret = relax_structure(c)
          ret_s = relax_structure(c, sparse = true)

          @test abs(abs(det(ret[1].A)) - abs(det(ret_s[1].A))) < 1e-3

          @test typeof(ret_s[2]) <: ThreeBodyTB.TB.tb_crys_sparse
          @test typeof(ret[2]) <: ThreeBodyTB.TB.tb_crys_dense

          @test typeof(ret[2]) <: ThreeBodyTB.TB.tb_crys   #both subtypes of the abstract supertype tb_crys
          @test typeof(ret_s[2]) <: ThreeBodyTB.TB.tb_crys
          
        end
    end
end


test_sparse();
Nothing;

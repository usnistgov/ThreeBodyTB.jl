using Test
using ThreeBodyTB
using Suppressor

#basic loading

#name of test directory
TESTDIR=ThreeBodyTB.TESTDIR

@testset "basic crystal dftout tests" begin
    types=["Li"]
    pos=zeros((1,3))
    A=[ [10.0 0 0]; [0 10.0 0]; [ 0 0 10.0]]

    c=makecrys(A, pos, types, units="Bohr")

    @test c.A == A
    @test c.coords == pos
    @test c.types == types

    c2=makecrys("$TESTDIR/POSCAR_testing")
    @test c2.types == ["Cl", "N","N","Sr", "Sr"]

    c3=makecrys("$TESTDIR/qe_testing.in")
    @test c3.types == ["Si", "C"]

    
    forces=ones((1,3))
    energy=0.0
    energy_smear = 0.0
    d=makedftout(c, energy, energy_smear, forces, zeros(3,3))
    @test energy == d.energy
    @test forces == d.forces

    d2=makedftout(A,pos,types, energy, energy_smear, forces, zeros(3,3))
    @test d.forces == d2.forces
    
    #supposed to throw an error
    @suppress begin
        forces2=zeros((2,3))
        @test_throws ErrorException("Error forces") makedftout(c, energy,energy_smear, forces2, zeros(3,3))
    end

    dft_out = ThreeBodyTB.QE.loadXML("$TESTDIR/dimer.save/")

    convert_ryd_ha = 2.0

    @test convert_ryd_ha*-1.428445472333963e1 ≈ dft_out.energy
    @test convert_ryd_ha*-6.425031283667987e-2 ≈ dft_out.bandstruct.efermi
    @test convert_ryd_ha*-1.832485909355843e0 ≈ dft_out.bandstruct.eigs[1,1]
    
end

@testset "advanced crystal tests" begin
    @suppress begin
        types=["Li"]
        pos=zeros((1,3))
        A=[ [10.0 0 0]; [0 10.0 0]; [ 0 0 10.0]]
        
        c=makecrys(A, pos, types, units="Bohr")
        
#        println(c)
        c2 = c * [2,2,2]
     
        ThreeBodyTB.CrystalMod.generate_random_distortion(c2,0.01,0.01)
        ThreeBodyTB.CrystalMod.write_poscar(c2,"$TESTDIR/POSCAR_tmp")
        rm("$TESTDIR/POSCAR_tmp")

        ThreeBodyTB.CrystalMod.write_efs(c, 0.0, zeros(1,3), zeros(3,3), "$TESTDIR/qe.out")
        rm("$TESTDIR/qe.out")
        

        
        @test 1 == 1

        
    end
end



if false

    #run QE . You actually need quantum espresso working to run this test, not appropriate for CI.
    @testset "run QE" begin
        types=["Li"]
        pos=zeros((1,3))
        A=[ [6.0 0 0]; [0 6.0 0]; [ 0 0 6.0]]
        
        c=makecrys(A, pos, types, units="Bohr")
        
        d = ThreeBodyTB.DFT.runSCF(c, nprocs=4)
        
        @test -14.35143737 ≈ d.energy

    end
end

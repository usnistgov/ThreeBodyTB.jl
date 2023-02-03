using ThreeBodyTB
using Test
using Suppressor

function test1()

    @testset "testing ewald" begin
        
        @suppress begin 
            cscl = makecrys([1.0 0 0; 0 1.0 0; 0 0 1.0]*2/sqrt(3), [0 0 0; 0.5 0.5 0.5], ["Na", "Cl"], units="Bohr");
            
            gamma_ij_tot_k1, bc = ThreeBodyTB.Ewald.electrostatics_getgamma(cscl, kappa = 1.0, noU=true)
            gamma_ij_tot_k2, bc = ThreeBodyTB.Ewald.electrostatics_getgamma(cscl, kappa = 2.0, noU=true)
            
            rs = makecrys([1.0 1.0 0; 1.0 0.0 1.0; 0 1.0 1.0], [0 0 0; 0.5 0.5 0.5], ["Na", "Cl"], units="Bohr");
            
            kappa = 3.0
            gamma_ij_tot_rs, bc = ThreeBodyTB.Ewald.electrostatics_getgamma(rs, kappa = kappa, noU=true)


            @test  isapprox(gamma_ij_tot_k1[1,2] - gamma_ij_tot_k1[1,1], gamma_ij_tot_k2[1,2] - gamma_ij_tot_k2[1,1], rtol=1e-5) #do different kappa values agree?
            @test  isapprox(sum( (gamma_ij_tot_k1[1,2] - gamma_ij_tot_k1[1,1])), 2*1.76268, rtol = 1e-4) #madelung constant cscl
            @test  isapprox(sum( (gamma_ij_tot_rs[1,2] - gamma_ij_tot_rs[1,1])), 2*1.74757, rtol = 1e-4) #madelung constant nacl

        end
        

#        rs_real = CrystalMod.makecrys([10.0 10.0 0; 10.0 0.0 10.0; 0 10.0 10.0], [0 0 0; 0.5 0.5 0.5], ["Na", "Cl"], units="Bohr");
#        gamma_ij_tot_rs = Ewald.electrostatics_getgamma(rs_real)
#        rs_real = CrystalMod.makecrys([3.0 3.0 0; 3.0 0.0 3.0; 0 3.0 3.0], [0 0 0; 0.5 0.5 0.5], ["Na", "Cl"], units="Bohr");
#        gamma_ij_tot_rs = Ewald.electrostatics_getgamma(rs_real)

#        for kappa = exp10.(range(log10(0.05), stop=log10(10), length=30))
#            println("KAPPA $kappa--------------------------------")
#            gamma_ij_tot_rs = Ewald.electrostatics_getgamma(rs_real, kappa = kappa)
#        end

    end
end

function test2()

    @testset "testing ewald units - dipole" begin

        @suppress begin 
        
            c = makecrys([1.0 0 0; 0 1.0 0; 0 0 1.0]*100, [0 0 0; 0.0 0.0 0.01], ["Na", "Cl"], units="Bohr");
            
            gamma, bc = ThreeBodyTB.Ewald.electrostatics_getgamma(c, noU=true)
            
            energy, pot = ThreeBodyTB.TB.ewald_energy(c, gamma, bc, [1.0, -1.0])

#        println("energy $energy")

            @test  isapprox(energy, -2.0, rtol=1e-5) #energy in rydberg e2=2; r=1;q=1; e2*q*-q/r  == -2, other atoms far away
        end
        


    end
end


test1()
test2()

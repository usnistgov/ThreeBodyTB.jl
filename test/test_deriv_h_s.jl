using Test
using ThreeBodyTB
using Suppressor

#include("../includes_laguerre.jl")
#include("../Ewald.jl")

TESTDIR=ThreeBodyTB.TESTDIR

function test_deriv()

    @testset "test derivative of ham" begin
        
      @suppress  begin
#        if true

          ThreeBodyTB.set_units(both="atomic")

          
          types=["Al"];
          
          #positions, crystal units
          pos=zeros((1,3));
          
          #lattice vectors, in Bohr units currently
          A=[ [3.7962751610 3.7962751610 0]; [3.7962751610 0 3.7962751610 ]; [ 0 3.7962751610 3.7962751610]];
          
          #makes the crystal
          c=makecrys(A, pos, types, units="Bohr")

          energy, tbc, flag = scf_energy(c)
          
          dk = 0.002
          K = -dk:dk:0.1
          nk =length(K)
          
          ind1 = 2
          ind2 = 2
          HK = zeros(Complex{Float64}, nk)
          SK = zeros(Complex{Float64}, nk)
          dHK = zeros(Complex{Float64}, nk)
          dSK = zeros(Complex{Float64}, nk)
          
          for (counter, k) in enumerate(K)
#              println([0 0 k])
              vects, vals, hk, sk, vals0,dhk, dsk = ThreeBodyTB.TB.Hk_derivative(tbc, [0.1 0.1 k])
              HK[counter] = hk[ind1, ind2]
              SK[counter] = sk[ind1, ind2]
              dHK[counter] = dhk[ind1, ind2, 3]
              dSK[counter] = dsk[ind1, ind2, 3]
          end

          dHK_finited = zeros(Complex{Float64}, nk-2)
          dSK_finited = zeros(Complex{Float64}, nk-2)
          
          for counter in 2:(nk-1)
              println("counter $counter")
              dHK_finited[counter-1] = (HK[counter+1] - HK[counter-1]) / (2*dk)
              dSK_finited[counter-1] = (SK[counter+1] - SK[counter-1]) / (2*dk)
          end

          plot(K[2:end-1], real( dHK_finited), label="dHk finite diff $ind1 $ind2", linewidth=2)
          plot!(K[2:end-1], real( dHK[2:end-1]), label="dHk formula $ind1 $ind2", linestyle=:dash, linewidth=2)

          plot!(K[2:end-1], real( dSK_finited), label="dSk finite diff $ind1 $ind2", linewidth=2)
          plot!(K[2:end-1], real( dSK[2:end-1]), label="dSk formula $ind1 $ind2", linestyle=:dash, linewidth=2)
          xlabel!("K (fractional)")
          ylabel!("d (Matrix Element) / dk ")
              
          
    end
end


test_deriv()

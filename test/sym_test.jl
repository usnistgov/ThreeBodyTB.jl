using Test
using ThreeBodyTB


function test1()

    @testset "testing sym1" begin
        abc = rand(3)
        lmn = rand(3)
        lmn = lmn / sqrt(sum(lmn.^2))

        lmnx = [1, 0, 0]
        lmny = [0, 1, 0]

        tot1 = 0.0
        tot2 = 0.0
        for orb1 = [:s, :pz, :px, :py, :dz2,:dxz,:dyz,:dx2_y2,:dxy]
            for orb2 = [:s, :pz, :px, :py, :dz2,:dxz,:dyz,:dx2_y2,:dxy]

                t = ThreeBodyTB.CalcTB.symmetry_factor(orb1, orb2, lmn, abc)
                t2 = ThreeBodyTB.CalcTB.symmetry_factor_fit(orb1, orb2, lmn)


#                println(orb1," ", orb2, "   ", t, "   " , sum(t2.*abc) )
#                if isapprox(t, sum(t2.*abc), rtol=1e-9) == false
#                    println("ERROR $orb1 $orb2")
#                end
                @test isapprox(t, sum(t2.*abc), rtol=1e-9)
                   
            end
        end             
        for orb1 = [:s, :pz, :px, :py, :dxy, :dz2, :dxz, :dyz]
            for orb2 = [:s, :pz, :px, :py, :dxy, :dz2, :dxz, :dyz]
                tot1 += ThreeBodyTB.CalcTB.symmetry_factor(orb1, orb2, lmnx, abc)
                tot2 += ThreeBodyTB.CalcTB.symmetry_factor(orb1, orb2, lmny, abc)
            end
        end

#        println("tot1 $tot1 tot2 $tot2 ")
        @test isapprox(tot1, tot2, rtol=1e-5)
        
    end
end


test1()
#test2()

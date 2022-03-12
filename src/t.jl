using LinearAlgebra


function go(H, S, W, X, nel, start, mix=0.5)

    charge = deepcopy(start)
    nat = Int64(size(H)[1]/2)

    del = [1,-1]
    energy = 0.0

    for iter = 1:20
        
        h1 = zeros(2*nat,2*nat)
        for spin = 1:2
            

            for i = 1:nat
                h1[i+(spin-1)*nat,i+(spin-1)*nat] += del[spin] * (charge[i] - charge[i+nat]) * W +  abs(charge[i] - charge[i+nat]) * X
            end

        end
        println(h1)
        Hm = H + h1
        vals, vects = eigen(Hm)
        
        eband = sum(real.(diag(vects[:,1:nel]' * H * vects[:,1:nel])))
#        println("nat $nat size charge ", size(charge))
        m = charge[1:nat] - charge[nat .+ (1:nat)]
        emag = sum(0.5 * (m.^2) * W+ 0.5 * (abs.(m).^2) * X)
        
        energy = emag + eband
        println("energy $energy emag $emag eband $eband m $m charge $charge")

        charge = charge * (1-mix) + mix * sum(vects[:,1:nel].*conj(vects[:,1:nel]), dims=2)

    end

    return energy, charge

end

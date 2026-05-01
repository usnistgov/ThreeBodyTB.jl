using Polynomials
using SpecialPolynomials
using StaticUnivariatePolynomials

const chebmax = 20

cheb_basis = basis.(Chebyshev{Float64}, 0:chebmax)
cheb_poly = Polynomials.Polynomial[]
cheb_spoly = StaticUnivariatePolynomials.Polynomial[]

for b in cheb_basis
    push!(cheb_poly, Polynomials.Polynomial(b))
    push!(cheb_spoly, StaticUnivariatePolynomials.Polynomial(Int64.(collect(Polynomials.Polynomial(b)))  ...,))
end


cheb_poly_sub = Polynomials.Polynomial[]
cheb_spoly_sub = StaticUnivariatePolynomials.Polynomial[]

for n in 3:length(cheb_poly)
    if mod(n, 2)  == 0
        push!(cheb_poly_sub, cheb_poly[n] + 1)
        push!(cheb_spoly_sub, cheb_spoly[n] + 1)
    elseif mod(n, 2)  == 1
        push!(cheb_poly_sub, -cheb_poly[n] + 1)
        push!(cheb_spoly_sub, -cheb_spoly[n] + 1)
    end
end



function cheb_energy_fn(coefs, rho, nmax, rho_max)
    en = 0.0
#    println("cheb_energy_fn")
    for i = 1:nmax
#        println("i $i $(coefs[i]) rho $rho x $(-1.0 + rho/rho_max * 2) cheb $(cheb_spoly[i](-1.0 + rho/rho_max * 2 )) tot $(coefs[i] * cheb_spoly[i](-1.0 + rho/rho_max * 2 ))")
        en += coefs[i] * cheb_spoly_sub[i](-1.0 + rho/rho_max * 2 )
    end
    return en #* (1 - cutoff_fn_fast(rho, 0, 0.1))
end

function cheb_energy_fn_prepare(rho, N_cheb, rho_max)
    en = 0.0

    ret = zeros(N_cheb)
    
    for i = 1:N_cheb
        ret[i] = cheb_spoly_sub[i](-1.0 + rho/rho_max * 2 )
    end
    return ret
end


function get_g(dist, r_loc)

#    println("dist $dist r_loc $r_loc m $m n_cheb $n_cheb")
#    println("dist convert $(-1 + 2*dist/r_loc)")
    ##gg = g(dist, m, r_loc, norm=1.0)

    gg = exp(-r_loc * dist) 
    
    return gg
    #println("gg $gg")
    #rho = zeros(n_cheb)
    #gg_norm = -1 + 2*dist/r_loc

    #for n = 1:n_cheb
    #    rho[n] += cheb_spoly[n](gg)
    #    println("rho n $n $(cheb_spoly[n](gg))")
    #end        
    return rho
end


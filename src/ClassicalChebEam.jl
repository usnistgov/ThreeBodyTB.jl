
using Polynomials
using SpecialPolynomials
using StaticUnivariatePolynomials


cheb_basis = basis.(Chebyshev{Float64}, 0:chebmax)
cheb_poly = Polynomials.Polynomial[]
cheb_spoly = StaticUnivariatePolynomials.Polynomial[]

for b in cheb_basis
    push!(cheb_poly, Polynomials.Polynomial(b))
    push!(cheb_spoly, StaticUnivariatePolynomials.Polynomial(Int64.(collect(Polynomials.Polynomial(b)))  ...,))
end


function cheb_energy_fn(coefs, rho, nmax, rho_max)
    en = 0.0
#    println("cheb_energy_fn")
    for i = 1:nmax
#        println("i $i $(coefs[i]) rho $rho x $(-1.0 + rho/rho_max * 2) cheb $(cheb_spoly[i](-1.0 + rho/rho_max * 2 )) tot $(coefs[i] * cheb_spoly[i](-1.0 + rho/rho_max * 2 ))")
        en += coefs[i] * cheb_spoly[i](-1.0 + rho/rho_max * 2 )
    end
    return en #* (1 - cutoff_fn_fast(rho, 0, 0.1))
end


function get_g(dist, r_loc, m, n_cheb)

#    println("dist $dist r_loc $r_loc m $m n_cheb $n_cheb")
#    println("dist convert $(-1 + 2*dist/r_loc)")
    ##gg = g(dist, m, r_loc, norm=1.0)

    gg = exp(-r_loc * dist) * cutoff_fn_fast(dist, 14.0, 16.0)
    
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

function g(dist, m, rmax; norm = 1.0)
#    println("inside g dist $dist m $m rmax $rmax")
    if dist < rmax
        return norm * (rmax - dist)^m
    else
        return 0.0
    end
end


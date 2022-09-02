module MyOptim

using LinearAlgebra
using Suppressor

function conjgrad(fn, grad, x0; maxstep=2.0, niters=50, conv_thr = 1e-2, fn_conv = 1e-4, verbosity="low")

    println()
    println("Conj Grad START")
    x = deepcopy(x0)
    storage = zeros(size(x))

    xold = deepcopy(x)

#    println("x0 ", x)

#    println("CG RUN GRAD1")
    #iter 1
    f, g = grad(storage, x)

    fold1 = 1e8
    fold2 = 1e8

#    println("f ", f, " x   ",  x)

    dx = deepcopy(-g)

    dx1 = deepcopy(dx)

    s = deepcopy(dx)
    
    step_size = maxstep/10.0

#    println("CG RUN LS1")

    old_step=step_size
    x, step_size = linesearch(x, dx, fn, f, step_size, verbosity)
    step_size = min(step_size, maxstep)



    for i = 1:niters
        

        if sum(abs.(xold - x)) > 1e-7  #only update grad if x changes

            fold2 = fold1
            fold1 = f

            f, g = grad(storage, x)
            xold = deepcopy(x)
#        else
#            println("nochange to grad")
        end

        if (sqrt(sum(g.^2)) < conv_thr && f < 0.1) || ( abs(fold1 - fold2) < fn_conv &&  abs(fold1 - f) < fn_conv) 
#            println("yes conv sum_grad ", sqrt(sum(g.^2)), " fn_diff_1 ", abs(fold1 - fold2), " fn_diff_2 ", abs(fold1 - f))
            return x, f, g
#        else
#            println("no conv ", sqrt(sum(g.^2)), " ", conv_thr, " f ", f)
        end

        if step_size < 0.5e-5
            println("WARNING, step size too small: $step_size , we give up")
            return x, f, g
        end
#        println("MY CONJGRAD $i ", f, " sg ", sqrt(sum(g.^2)), " step_size $step_size --------------------")


        dx1[:] = dx[:]
        dx[:] = -g[:]

#        println("dx ", dx)
#        println("dx1 ", dx1)

        beta = ( dx' * (dx - dx1)) / (dx1' * dx1)  #PR
#        println("beta $beta")
        beta = max(0, beta)

        s = dx + beta * s

        old_step=step_size

        x, step_size, good = linesearch(x, s, fn, f, step_size, verbosity)
        if good
            step_size =min(step_size, maxstep)
        else
            beta = 0.0
        end
        println("MY CG quadratic linesearch iter $i fn val: $f | rms grad: ", sqrt(sum(g.^2)), " within_stepsize: $good , old step: $old_step , new step_size: $step_size ")
        

    end

    println("WARNING, conv not achieved")
    return x, f, g

end

function linesearch(x, dx, fn, f0, step_size, verbosity)
#    println("CG FN0 $f0")
#    println("CG dx $dx")
#    println("ls $x $dx $f0 $step_size")
#    println("CG FN1 ss $step_size")

    if verbosity=="low"
        @suppress f0r, flag0r = fn( x + dx * step_size * 0.0)
    else
        f0r, flag0r = fn( x + dx * step_size * 0.0)
    end

#    println("CG real $f0r $flag0r")

    if verbosity=="low"
        @suppress f1, flag1 = fn( x + dx * step_size * 0.5)
    else
        f1, flag1 = fn( x + dx * step_size * 0.5)
    end

#    println("CG FN1 $f1 $flag1")
    if f1 > f0 || flag1 == false
#        println("MY LS uphill $f1 > $f0 flag $flag1, reduce step")
        return x, step_size/3.1, false
    end
#    println("CG FN2")    

    if verbosity == "low"
        @suppress f2, flag2 = fn( x + dx * step_size)
    else
        f2, flag2 = fn( x + dx * step_size)
    end

#    println("CG FN2 $f2 $flag2")

    if flag2 == false
#        println("MY LS flag2 $flag2")
        return x + dx * step_size * 0.5, step_size/1.75, false
    end

    c = f0
    a = 2 * (f2-2*f1+f0)
    b = f2 - a - c

#    println("a b c $a $b $c ", -b/(2*a))

    if a < 0
#        println("MY LS a negative $a, take max step increase stepsize")
        return x + dx * step_size * 1.0, step_size * 2.0, true
    else
        step = min(-b/(2*a), 1.0)
#        println("MY LS NORMAL $step")
        return x + dx * step * step_size, step * step_size * 1.5, true
    end
    

end


end #end module

#=
function myf(x)

    return -0.5*x'*A*x + 0.1*sum(x.^4)

end

using ForwardDiff



function grad(storage, x)
    sleep(0.1)

    return myf(x), ForwardDiff.gradient(myf, x)

end

function grad_only(storage, x)
    storage[:] = ForwardDiff.gradient(myf, x)
    sleep(0.1)
    return storage
end
=#

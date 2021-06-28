using ThreeBodyTB
using JLD2

fillist = readdir()

for f in fillist
    if occursin("jld2",f)
        if occursin("2body" ,f )
            dim = 2
        elseif occursin("3body", f)
            dim = 3
        else
            println("err $f")
        end

        t = @load f
        savename=f[1:end-5]
        c = ThreeBodyTB.CalcTB.make_coefs(dat.names, dim, datH=dat.datH, datS=dat.datS, min_dist = dat.min_dist, maxmin_val_train=dat.maxmin_val_train, dist_frontier=dat.dist_frontier)
        ThreeBodyTB.CalcTB.write_coefs("$savename.xml", c)

    end
end

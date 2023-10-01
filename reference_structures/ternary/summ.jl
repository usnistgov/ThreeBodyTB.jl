using ThreeBodyTB

set_units(both="atomic")
d = readdir("/home/kfg/codes/jarvis_test/tern_7_12")

for line in d
#    println(line)
#    try
    #sp = split(line)
    #    sp2 = split(line, "/")
     #   println(sp)
    c = makecrys("/home/kfg/codes/jarvis_test/tern_7_12/$line")
        dict = Dict()
        for t in c.types
            if !(t in keys(dict))
                dict[t] = 0
            end
            dict[t] += 1
        end
        st = ""
        st2 = ""
        for k  in sort(collect(keys(dict)))
            st *= k*"$(dict[k])"
            st2 *= "$k "
        end
        
        println("/home/kfg/codes/jarvis_test/tern_7_12/"*line ," ",  line, " " , st , " $(c.nat) " , st2)

#    cs = Set(c.stypes)
#    for cc in cs
#        ind = c.stypes .!= cc
#        c2 = makecrys(c.A, c.coords[ind,:], c.stypes[ind])
#        ThreeBodyTB.CrystalMod.write_poscar(c2, "/home/kfg/codes/jarvis_test/tern_7_12_remove1/rm.$cc.$line")
#    end

#    catch
#        continue
#    end
end

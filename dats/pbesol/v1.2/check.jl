using ThreeBodyTB

function go()

    res = []
    data = open("mindist2.csv", "w")

    ThreeBodyTB.ManageDatabase.add_to_database(Set([:H,:H]), directory="../v1.3", verbose=false)
    m = ThreeBodyTB.ManageDatabase.database_cached[(:H, :H)].min_dist


    for x in  [:H,  :He,  :Li,  :Be,  :B,  :C,  :N,  :O,  :F,  :Ne,  :Na,  :Mg,  :Al,  :Si,  :P,  :S,  :Cl,  :Ar,  :K,  :Ca,  :Sc,  :Ti,  :V,  :Cr,  :Mn,  :Fe,  :Co,  :Ni,  :Cu,  :Zn,  :Ga,  :Ge,  :As,  :Se,  :Br,  :Kr,  :Rb,  :Sr,  :Y,  :Zr,  :Nb,  :Mo,  :Tc,  :Ru,  :Rh,  :Pd,  :Ag,  :Cd,  :In,  :Sn,  :Sb,  :Te,  :I,  :Xe,  :Cs,  :Ba,  :La,  :Ce,  :Pr,  :Nd,  :Pm,  :Sm,  :Eu,  :Gd,  :Tb,  :Dy,  :Ho,  :Er,  :Tm,  :Yb,  :Lu,  :Hf,  :Ta,  :W,  :Re,  :Os,  :Ir,  :Pt,  :Au,  :Hg,  :Tl,  :Pb,  :Bi]
        if x in keys(ThreeBodyTB.Atomdata.atoms)
            for y in  [:H,  :He,  :Li,  :Be,  :B,  :C,  :N,  :O,  :F,  :Ne,  :Na,  :Mg,  :Al,  :Si,  :P,  :S,  :Cl,  :Ar,  :K,  :Ca,  :Sc,  :Ti,  :V,  :Cr,  :Mn,  :Fe,  :Co,  :Ni,  :Cu,  :Zn,  :Ga,  :Ge,  :As,  :Se,  :Br,  :Kr,  :Rb,  :Sr,  :Y,  :Zr,  :Nb,  :Mo,  :Tc,  :Ru,  :Rh,  :Pd,  :Ag,  :Cd,  :In,  :Sn,  :Sb,  :Te,  :I,  :Xe,  :Cs,  :Ba,  :La,  :Ce,  :Pr,  :Nd,  :Pm,  :Sm,  :Eu,  :Gd,  :Tb,  :Dy,  :Ho,  :Er,  :Tm,  :Yb,  :Lu,  :Hf,  :Ta,  :W,  :Re,  :Os,  :Ir,  :Pt,  :Au,  :Hg,  :Tl,  :Pb,  :Bi]
                if y in keys(ThreeBodyTB.Atomdata.atoms)
                    try
                        ThreeBodyTB.ManageDatabase.add_to_database(Set([x,y]), directory="../v1.3", verbose=false)
                        m = ThreeBodyTB.ManageDatabase.database_cached[(x, y)].min_dist
                    catch
                        ThreeBodyTB.ManageDatabase.add_to_database(Set([x,y]), directory="./", verbose=false)
                        m = ThreeBodyTB.ManageDatabase.database_cached[(x, y)].min_dist
                    end                    
                    push!(res, [x, y, m])
                    println("$x $y $m")
                    println(data, "$x $y $m")
                end
            end
        end
    end
    close(data)
end
go()

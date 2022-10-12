
module Symmetry

using Spglib
using ..CrystalMod:crystal
using ..CrystalMod:make_crys


function get_symmetry(c::crystal; verbose=true, sym_prec = 1e-5, magmoms=missing)

    coords = [c.coords[i,:] for i in 1:c.nat]
    
    if ismissing(magmoms)
        cell = Spglib.Cell(c.A, coords, c.types)
    else
        cell = Spglib.Cell(c.A, c.coords, c.types, magmoms=magmoms)
    end
    
    dat = Spglib.get_dataset(cell, sypprec = sym_prec)
    
    if verbose
        println("Symmetry info")
        println("Space group # $(dat.spacegroup_number)     $(dat.international_symbol)      $(dat.hall_symbol)")
        println("Point group   $(dat.pointgroup_symbol)")
    end

    return dat.spacegroup_number

end


function get_standard_crys(c::crystal; sym_prec = 1e-5, magmoms=missing)

    coords = [c.coords[i,:] for i in 1:c.nat]
    if ismissing(magmoms)
        cell = Spglib.Cell(c.A, coords, c.types)
    else
        cell = Spglib.Cell(c.A, coords, c.types, magmoms=magmoms)
    end
    
    cell_std = Spglib.standardize_cell(cell, to_primitive=false, symprec=symprec)

    #cell_std = Spglib.standardize_cell(cell, to_primitive=false, symprec=symprec)
    #cell_std = Spglib.find_primitive(cell_std)

    coords = zeros(c.nat, 3)
    t = cell_std.positions
    for i = 1:length(t)
        coords[i,:] = t[i][:]
    end

    c_std = make_crys(collect(cell_std.lattice), coords, cell.types)

    return c_std
    
end

function get_kpath(c::crystal; sym_prec = 1e-5, magmoms = missing)

    c_std = get_standard_crys(c, sym_prec = 1e-5, magmoms=magmoms)
    sym = get_symmetry(c_std, verbose=false, sym_prec = sym_prec, magmoms=magmoms)
    
    spg_sym = Dict()
    spg_sym[1] = "P1"
    spg_sym[2] = "P-1"
    spg_sym[3] = "P2"
    spg_sym[4] = "P21"
    spg_sym[5] = "C2"
    spg_sym[6] = "Pm"
    spg_sym[7] = "Pc"
    spg_sym[8] = "Cm"
    spg_sym[9] = "Cc"
    spg_sym[10] = "P2/m"
    spg_sym[11] = "P21/m"
    spg_sym[12] = "C2/m"
    spg_sym[13] = "P2/c"
    spg_sym[14] = "P21/c"
    spg_sym[15] = "C2/c"
    spg_sym[16] = "P222"
    spg_sym[17] = "P2221"
    spg_sym[18] = "P21212"
    spg_sym[19] = "P212121"
    spg_sym[20] = "C2221"
    spg_sym[21] = "C222"
    spg_sym[22] = "F222"
    spg_sym[23] = "I222"
    spg_sym[24] = "I212121"
    spg_sym[25] = "Pmm2"
    spg_sym[26] = "Pmc21"
    spg_sym[27] = "Pcc2"
    spg_sym[28] = "Pma2"
    spg_sym[29] = "Pca21"
    spg_sym[30] = "Pnc2"
    spg_sym[31] = "Pmn21"
    spg_sym[32] = "Pba2"
    spg_sym[33] = "Pna21"
    spg_sym[34] = "Pnn2"
    spg_sym[35] = "Cmm2"
    spg_sym[36] = "Cmc21"
    spg_sym[37] = "Ccc2"
    spg_sym[38] = "Amm2"
    spg_sym[39] = "Aem2"
    spg_sym[40] = "Ama2"
    spg_sym[41] = "Aea2"
    spg_sym[42] = "Fmm2"
    spg_sym[43] = "Fdd2"
    spg_sym[44] = "Imm2"
    spg_sym[45] = "Iba2"
    spg_sym[46] = "Ima2"
    spg_sym[47] = "Pmmm"
    spg_sym[48] = "Pnnn"
    spg_sym[49] = "Pccm"
    spg_sym[50] = "Pban"
    spg_sym[51] = "Pmma"
    spg_sym[52] = "Pnna"
    spg_sym[53] = "Pmna"
    spg_sym[54] = "Pcca"
    spg_sym[55] = "Pbam"
    spg_sym[56] = "Pccn"
    spg_sym[57] = "Pbcm"
    spg_sym[58] = "Pnnm"
    spg_sym[59] = "Pmmn"
    spg_sym[60] = "Pbcn"
    spg_sym[61] = "Pbca"
    spg_sym[62] = "Pnma"
    spg_sym[63] = "Cmcm"
    spg_sym[64] = "Cmce"
    spg_sym[65] = "Cmmm"
    spg_sym[66] = "Cccm"
    spg_sym[67] = "Cmme"
    spg_sym[68] = "Ccce"
    spg_sym[69] = "Fmmm"
    spg_sym[70] = "Fddd"
    spg_sym[71] = "Immm"
    spg_sym[72] = "Ibam"
    spg_sym[73] = "Ibca"
    spg_sym[74] = "Imma"
    spg_sym[75] = "P4"
    spg_sym[76] = "P41"
    spg_sym[77] = "P42"
    spg_sym[78] = "P43"
    spg_sym[79] = "I4"
    spg_sym[80] = "I41"
    spg_sym[81] = "P-4"
    spg_sym[82] = "I-4"
    spg_sym[83] = "P4/m"
    spg_sym[84] = "P42/m"
    spg_sym[85] = "P4/n"
    spg_sym[86] = "P42/n"
    spg_sym[87] = "I4/m"
    spg_sym[88] = "I41/a"
    spg_sym[89] = "P422"
    spg_sym[90] = "P4212"
    spg_sym[91] = "P4122"
    spg_sym[92] = "P41212"
    spg_sym[93] = "P4222"
    spg_sym[94] = "P42212"
    spg_sym[95] = "P4322"
    spg_sym[96] = "P43212"
    spg_sym[97] = "I422"
    spg_sym[98] = "I4122"
    spg_sym[99] = "P4mm"
    spg_sym[100] = "P4bm"
    spg_sym[101] = "P42cm"
    spg_sym[102] = "P42nm"
    spg_sym[103] = "P4cc"
    spg_sym[104] = "P4nc"
    spg_sym[105] = "P42mc"
    spg_sym[106] = "P42bc"
    spg_sym[107] = "I4mm"
    spg_sym[108] = "I4cm"
    spg_sym[109] = "I41md"
    spg_sym[110] = "I41cd"
    spg_sym[111] = "P-42m"
    spg_sym[112] = "P-42c"
    spg_sym[113] = "P-421m"
    spg_sym[114] = "P-421c"
    spg_sym[115] = "P-4m2"
    spg_sym[116] = "P-4c2"
    spg_sym[117] = "P-4b2"
    spg_sym[118] = "P-4n2"
    spg_sym[119] = "I-4m2"
    spg_sym[120] = "I-4c2"
    spg_sym[121] = "I-42m"
    spg_sym[122] = "I-42d"
    spg_sym[123] = "P4/mmm"
    spg_sym[124] = "P4/mcc"
    spg_sym[125] = "P4/nbm"
    spg_sym[126] = "P4/nnc"
    spg_sym[127] = "P4/mbm"
    spg_sym[128] = "P4/mnc"
    spg_sym[129] = "P4/nmm"
    spg_sym[130] = "P4/ncc"
    spg_sym[131] = "P42/mmc"
    spg_sym[132] = "P42/mcm"
    spg_sym[133] = "P42/nbc"
    spg_sym[134] = "P42/nnm"
    spg_sym[135] = "P42/mbc"
    spg_sym[136] = "P42/mnm"
    spg_sym[137] = "P42/nmc"
    spg_sym[138] = "P42/ncm"
    spg_sym[139] = "I4/mmm"
    spg_sym[140] = "I4/mcm"
    spg_sym[141] = "I41/amd"
    spg_sym[142] = "I41/acd"
    spg_sym[143] = "P3"
    spg_sym[144] = "P31"
    spg_sym[145] = "P32"
    spg_sym[146] = "R3"
    spg_sym[147] = "P-3"
    spg_sym[148] = "R-3"
    spg_sym[149] = "P312"
    spg_sym[150] = "P321"
    spg_sym[151] = "P3112"
    spg_sym[152] = "P3121"
    spg_sym[153] = "P3212"
    spg_sym[154] = "P3221"
    spg_sym[155] = "R32"
    spg_sym[156] = "P3m1"
    spg_sym[157] = "P31m"
    spg_sym[158] = "P3c1"
    spg_sym[159] = "P31c"
    spg_sym[160] = "R3m"
    spg_sym[161] = "R3c"
    spg_sym[162] = "P-31m"
    spg_sym[163] = "P-31c"
    spg_sym[164] = "P-3m1"
    spg_sym[165] = "P-3c1"
    spg_sym[166] = "R-3m"
    spg_sym[167] = "R-3c"
    spg_sym[168] = "P6"
    spg_sym[169] = "P61"
    spg_sym[170] = "P65"
    spg_sym[171] = "P62"
    spg_sym[172] = "P64"
    spg_sym[173] = "P63"
    spg_sym[174] = "P-6"
    spg_sym[175] = "P6/m"
    spg_sym[176] = "P63/m"
    spg_sym[177] = "P622"
    spg_sym[178] = "P6122"
    spg_sym[179] = "P6522"
    spg_sym[180] = "P6222"
    spg_sym[181] = "P6422"
    spg_sym[182] = "P6322"
    spg_sym[183] = "P6mm"
    spg_sym[184] = "P6cc"
    spg_sym[185] = "P63cm"
    spg_sym[186] = "P63mc"
    spg_sym[187] = "P-6m2"
    spg_sym[188] = "P-6c2"
    spg_sym[189] = "P-62m"
    spg_sym[190] = "P-62c"
    spg_sym[191] = "P6/mmm"
    spg_sym[192] = "P6/mcc"
    spg_sym[193] = "P63/mcm"
    spg_sym[194] = "P63/mmc"
    spg_sym[195] = "P23"
    spg_sym[196] = "F23"
    spg_sym[197] = "I23"
    spg_sym[198] = "P213"
    spg_sym[199] = "I213"
    spg_sym[200] = "Pm-3"
    spg_sym[201] = "Pn-3"
    spg_sym[202] = "Fm-3"
    spg_sym[203] = "Fd-3"
    spg_sym[204] = "Im-3"
    spg_sym[205] = "Pa-3"
    spg_sym[206] = "Ia-3"
    spg_sym[207] = "P432"
    spg_sym[208] = "P4232"
    spg_sym[209] = "F432"
    spg_sym[210] = "F4132"
    spg_sym[211] = "I432"
    spg_sym[212] = "P4332"
    spg_sym[213] = "P4132"
    spg_sym[214] = "I4132"
    spg_sym[215] = "P-43m"
    spg_sym[216] = "F-43m"
    spg_sym[217] = "I-43m"
    spg_sym[218] = "P-43n"
    spg_sym[219] = "F-43c"
    spg_sym[220] = "I-43d"
    spg_sym[221] = "Pm-3m"
    spg_sym[222] = "Pn-3n"
    spg_sym[223] = "Pm-3n"
    spg_sym[224] = "Pn-3m"
    spg_sym[225] = "Fm-3m"
    spg_sym[226] = "Fm-3c"
    spg_sym[227] = "Fd-3m"
    spg_sym[228] = "Fd-3c"
    spg_sym[229] = "Im-3m"
    spg_sym[230] = "Ia-3d"

    centering = spg_sym[sym][1]
    
    

    a = sqrt(sum(c_std.A[1,:].^2))
    b = sqrt(sum(c_std.A[2,:].^2))
    c = sqrt(sum(c_std.A[3,:].^2))

    alpha = acos(c_std.A[2,:] .* c_std.A[3,:] / b / c) * 180/pi
    beta =  acos(c_std.A[1,:] .* c_std.A[3,:] / a / c) * 180/pi
    gamma = acos(c_std.A[1,:] .* c_std.A[2,:] / a / b) * 180/pi


    if sym in 1:2
        lattice="triclinic"
    elseif sym in 3:15
        lattice="monoclinic"
    elseif sym in 16:74
        lattice="orthorhombic"
    elseif sym in 75:142
        lattice="tetragonal"
    elseif sym in 143:167
        if sym in [146, 148, 155, 160, 161, 166, 167]
            lattice = "rhombohedral"
        else
            lattice="hexagonal"
        end
    elseif sym in 168:194
        lattice="hexagonal"
    elseif sym in 195:230
        lattice="cubic"
    else
        println("sym error $sym")
        lattice="triclinic"
    end
    
    if lattice == "cubic"
        if centering == "F"
            cat = "FCC"
        elseif centering == "I"
            cat = "BCC"
        else
            cat = "SC"
        end
    
    elseif lattice == "tetragonal"
        if centering == "P"
            cat = "tet"
        else
            if c < a
                cat = "bctet1"
            else
                cat = "bctet2"
            end
        end
        
    elseif lattice == "orthorhombic"
        if centering == "P"
            cat = "orc"
        elseif centering == "I"
            cat = "orci"
        elseif centering == "C"
            cat = "orcc"
        else
            if 1/a^2 ≈ 1/b^2 + 1/c^2
                cat = "orcf3"
            elseif 1/a^2 > 1/b^2 + 1/c^2
                cat = "orcf1"
            elseif 1/a^2 < 1/b^2 + 1/c^2
                cat = "orcf2"
            end
        end

    elseif lattice == "hexagonal"
        cat = "hex"

    elseif lattice == "rhombohedral"
        if alpha < 90.0
            cat = "rhl1"
        else
            cat = "rhl2"
        end

    elseif lattice == "monoclinic"
        if centering == "P"
            cat = "mcl"
        else
            if gamma ≈ 90.0
                cat = "mclc2"
            elseif gamma > 90
                cat = "mclc1"
            else
                if     b * cos(alpha * pi / 180) / c + b ^ 2 * sin(alpha * pi / 180) ^2 2 / a ^2 2  ≈ 1
                    cat = "mclc4"
                
                elseif b * cos(alpha * pi / 180) / c + b ^ 2 * sin(alpha * pi / 180) ^ 2 / a ^ 2  < 1
                    cat = "mclc3"
                else
                    cat = "mclc5"
                end
            end
        end

    elseif lattice == "triclinic"
        
        if alpha > 90 && beta > 90 && gamma > 90
            cat = "tria"
        elseif alpha < 90 && beta < 90 && gamma < 90
            cat = "trib"
        elseif alpha > 90 && beta > 90 && gamma ≈ 90
            cat = "tria"
        elseif alpha < 90 && beta < 90 && gamma ≈ 90
            cat = "trib"
        else
            cat = "tria"
        end
        
    end

    #########
    
    kdict = Dict()
    if cat == "SC"
        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["X"] = [0.0, 0.5, 0.0]
        kdict["R"] = [0.5, 0.5, 0.5]
        kdict["M"] = [0.5, 0.5, 0.0]

        names = ["Γ", "X", "M", "Γ", "R", "X"]
    
    elseif cat == "FCC"
        
        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["K"] = [3.0 / 8.0, 3.0 / 8.0, 3.0 / 4.0]
        kdict["L"] = [0.5, 0.5, 0.5]
        kdict["U"] = [5.0 / 8.0, 1.0 / 4.0, 5.0 / 8.0]
        kdict["W"] = [0.5, 1.0 / 4.0, 3.0 / 4.0]
        kdict["X"] = [0.5, 0.0, 0.5]

        names =  ["Γ", "X", "W", "K", "Γ", "L", "U", "W", "L", "K"]

    elseif cat == "BCC"

        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["H"] = [0.5, -0.5, 0.5]
        kdict["P"] = [0.25, 0.25, 0.25]
        kdict["N"] = [0.0, 0.0, 0.5]
        
        names = ["Γ", "H", "N", "Γ", "P", "H"]

    elseif cat == "tet"

        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["A"] = [0.5, 0.5, 0.5]
        kdict["M"] = [0.5, 0.5, 0.0]
        kdict["R"] = [0.0, 0.5, 0.5]
        kdict["X"] = [0.0, 0.5, 0.0]
        kdict["Z"] = [0.0, 0.0, 0.5]

        names = ["Γ", "X", "M", "Γ", "Z", "R", "A", "Z"]

    elseif cat == "bctet1"

        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["M"] = [-0.5, 0.5, 0.5]
        kdict["N"] = [0.0, 0.5, 0.0]
        kdict["P"] = [0.25, 0.25, 0.25]
        kdict["X"] = [0.0, 0.0, 0.5]
        kdict["Z"] = [eta, eta, -eta]
        kdict["Z_1"] = [-eta, 1 - eta, eta]
        
        names = ["Γ", "X", "M", "Γ", "Z", "P", "N", "Z_1", "M"]

    elseif cat == "bctet2"

        eta = (1 + a^2/c^2) / 4
        zeta = a^2 / (2 * c^2)
        
        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["N"] = [0.0, 0.5, 0.0]
        kdict["P"] = [0.25, 0.25, 0.25]
        kdict["Σ"] = [-eta, eta, eta]
        kdict["Σ_1"] = [eta, 1 - eta, -eta]
        kdict["X"] = [0.0, 0.0, 0.5]
        kdict["Y"] = [-zeta, zeta, 0.5]
        kdict["Y_1"] = [0.5, 0.5, -zeta]
        kdict["Z"] = [0.5, 0.5, -0.5]
        
        names =    [ "Γ", "X",  "Y", "Σ", "Γ","Z", "Σ_1", "N", "P", "Y_1",   "Z",  ]

    elseif cat == "orc"
        
        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["R"] = [0.5, 0.5, 0.5]
        kdict["S"] = [0.5, 0.5, 0.0]
        kdict["T"] = [0.0, 0.5, 0.5]
        kdict["U"] = [0.5, 0.0, 0.5]
        kdict["X"] = [0.5, 0.0, 0.0]
        kdict["Y"] = [0.0, 0.5, 0.0]
        kdict["Z"] = [0.0, 0.0, 0.5]

        names =   ["Γ", "X", "S", "Y", "Γ", "Z", "U", "R", "T", "Z"]

    elseif cat = "orcf1"

        zeta = (1 + a ^ 2 / b ^ 2 - a ^ 2 / c ^ 2) / 4
        eta = (1 + a ^ 2 / b ^ 2 + a ^ 2 / c ^ 2) / 4

        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["A"] = [0.5, 0.5 + zeta, zeta]
        kdict["A_1"] = [0.5, 0.5 - zeta, 1 - zeta]
        kdict["L"] = [0.5, 0.5, 0.5]
        kdict["T"] = [1, 0.5, 0.5]
        kdict["X"] = [0.0, eta, eta]
        kdict["X_1"] = [1, 1 - eta, 1 - eta]
        kdict["Y"] = [0.5, 0.0, 0.5]
        kdict["Z"] = [0.5, 0.5, 0.0]

        names =     ["Γ", "Y", "T", "Z", "Γ", "X", "A_1", "Y"]

    elseif cat == "orcf2"

        phi = (1 + c ^ 2 / b ^ 2 - c ^ 2 / a ^ 2) / 4
        eta = (1 + a ^ 2 / b ^ 2 - a ^ 2 / c ^ 2) / 4
        delta = (1 + b ^ 2 / a ^ 2 - b ^ 2 / c ^ 2) / 4
        
        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["C"] = [0.5, 0.5 - eta, 1 - eta]
        kdict["C_1"] = [0.5, 0.5 + eta, eta]
        kdict["D"] = [0.5 - delta, 0.5, 1 - delta]
        kdict["D_1"] = [0.5 + delta, 0.5, delta]
        kdict["L"] = [0.5, 0.5, 0.5]
        kdict["H"] = [1 - phi, 0.5 - phi, 0.5]
        kdict["H_1"] = [phi, 0.5 + phi, 0.5]
        kdict["X"] = [0.0, 0.5, 0.5]
        kdict["Y"] = [0.5, 0.0, 0.5]
        kdict["Z"] = [0.5, 0.5, 0.0]

        names =    ["Γ", "Y", "C", "D", "X", "Γ", "Z", "D_1", "H", "C"]

    elseif cat == "orcf3"

        zeta = (1 + a ^ 2 / b ^ 2 - a ^ 2 / c ^ 2) / 4
        eta = (1 + a ^ 2 / b ^ 2 + a ^ 2 / c ^ 2) / 4

        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["A"] = [0.5, 0.5 + zeta, zeta]
        kdict["A_1"] = [0.5, 0.5 - zeta, 1 - zeta]
        kdict["L"] = [0.5, 0.5, 0.5]
        kdict["T"] = [1, 0.5, 0.5]
        kdict["X"] = [0.0, eta, eta]
        kdict["X_1"] = [1, 1 - eta, 1 - eta]
        kdict["Y"] = [0.5, 0.0, 0.5]
        kdict["Z"] = [0.5, 0.5, 0.0]

        names = ["Γ", "Y", "T", "Z", "Γ", "X", "A_1", "Y"]

    elseif cat == "orci"

        zeta = (1 + a ^ 2 / c ^ 2) / 4
        eta = (1 + b ^ 2 / c ^ 2) / 4
        delta = (b ^ 2 - a ^ 2) / (4 * c ^ 2)
        mu = (a ^ 2 + b ^ 2) / (4 * c ^ 2)

        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["L"] = [-mu, mu, 0.5 - delta]
        kdict["L_1"] = [mu, -mu, 0.5 + delta]
        kdict["L_2"] = [0.5 - delta, 0.5 + delta, -mu]
        kdict["R"] = [0.0, 0.5, 0.0]
        kdict["S"] = [0.5, 0.0, 0.0]
        kdict["T"] = [0.0, 0.0, 0.5]
        kdict["W"] = [0.25, 0.25, 0.25]
        kdict["X"] = [-zeta, zeta, zeta]
        kdict["X_1"] = [zeta, 1 - zeta, -zeta]
        kdict["Y"] = [eta, -eta, eta]
        kdict["Y_1"] = [1 - eta, eta, -eta]
        kdict["Z"] = [0.5, 0.5, -0.5]

        names =   [ "Γ", "X","L", "T", "W", "R", "X_1", "Z", "Γ", "Y", "S", "W" ]

    elseif cat == "orcc"

        zeta = (1 + a ^ 2 / b ^ 2) / 4

        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["A"] = [zeta, zeta, 0.5]
        kdict["A_1"] = [-zeta, 1 - zeta, 0.5]
        kdict["R"] = [0.0, 0.5, 0.5]
        kdict["S"] = [0.0, 0.5, 0.0]
        kdict["T"] = [-0.5, 0.5, 0.5]
        kdict["X"] = [zeta, zeta, 0.0]
        kdict["X_1"] = [-zeta, 1 - zeta, 0.0]
        kdict["Y"] = [-0.5, 0.5, 0]
        kdict["Z"] = [0.0, 0.0, 0.5]

        names = ["Γ",  "X",   "S",  "R",   "A", "Z",  "Γ", "Y",   "X_1",  "A_1", "T",  "Y"  ]
        

    elseif cat == "hex"

        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["A"] = [0.0, 0.0, 0.5]
        kdict["H"] = [1.0 / 3.0, 1.0 / 3.0, 0.5]
        kdict["K"] = [1.0 / 3.0, 1.0 / 3.0, 0.0]
        kdict["L"] = [0.5, 0.0, 0.5]
        kdict["M"] = [0.5, 0.0, 0.0]

        names = ["Γ", "M", "K", "Γ", "A", "L", "H", "A"]
        
    elseif cat == "rhl1"

        eta = (1 + 4 * cos(alpha)) / (2 + 4 * cos(alpha))
        nu = 3.0 / 4.0 - eta / 2.0

        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["B"] = [eta, 0.5, 1.0 - eta]
        kdict["B_1"] = [1.0 / 2.0, 1.0 - eta, eta - 1.0]
        kdict["F"] = [0.5, 0.5, 0.0]
        kdict["L"] = [0.5, 0.0, 0.0]
        kdict["L_1"] = [0.0, 0.0, -0.5]
        kdict["P"] = [eta, nu, nu]
        kdict["P_1"] = [1.0 - nu, 1.0 - nu, 1.0 - eta]
        kdict["P_2"] = [nu, nu, eta - 1.0]
        kdict["Q"] = [1.0 - nu, nu, 0.0]
        kdict["X"] = [nu, 0.0, -nu]
        kdict["Z"] = [0.5, 0.5, 0.5]

        names = ["Γ", "L", "B_1", "Γ", "B", "Z", "Γ", "X", "Q", "F", "L", "P"]

    elseif cat == "rhl2"
        
        eta = 1 / (2 * tan(alpha / 2.0) ^ 2)
        nu = 3.0 / 4.0 - eta / 2.0

        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["F"] = [0.5, -0.5, 0.0]
        kdict["L"] = [0.5, 0.0, 0.0]
        kdict["P"] = [1 - nu, -nu, 1 - nu]
        kdict["P_1"] = [nu, nu - 1.0, nu - 1.0]
        kdict["Q"] = [eta, eta, eta]
        kdict["Q_1"] = [1.0 - eta, -eta, -eta]
        kdict["Z"] = [0.5, -0.5, 0.5]

        names = ["Γ", "P", "Z", "Q", "Γ", "F", "P_1", "Q_1", "L", "Z"]

    elseif cat == "mcl"

        eta = (1 - b * cos(beta) / c) / (2 * sin(beta) ^ 2)
        nu = 0.5 - eta * c * cos(beta) / b
        
        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["A"] = [0.5, 0.5, 0.0]
        kdict["C"] = [0.0, 0.5, 0.5]
        kdict["D"] = [0.5, 0.0, 0.5]
        kdict["D_1"] = [0.5, 0.5, -0.5]
        kdict["E"] = [0.5, 0.5, 0.5]
        kdict["H"] = [0.0, eta, 1.0 - nu]
        kdict["H_1"] = [0.0, 1.0 - eta, nu]
        kdict["H_2"] = [0.0, eta, -nu]
        kdict["M"] = [0.5, eta, 1.0 - nu]
        kdict["M_1"] = [0.5, 1 - eta, nu]
        kdict["M_2"] = [0.5, 1 - eta, nu]
        kdict["X"] = [0.0, 0.5, 0.0]
        kdict["Y"] = [0.0, 0.0, 0.5]
        kdict["Y_1"] = [0.0, 0.0, -0.5]
        kdict["Z"] = [0.5, 0.0, 0.0]

        names = ["Γ", "Y", "H", "C", "E", "M_1", "A", "X", "H_1", "Γ", "M", "D", "Z", "Y"]

    elseif cat == "mclc1"

        zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha) ^ 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        psi = 0.75 - a ^ 2 / (4 * b ^ 2 * sin(alpha) ^ 2)
        phi = psi + (0.75 - psi) * b * cos(alpha) / c

        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["N"] = [0.5, 0.0, 0.0]
        kdict["N_1"] = [0.0, -0.5, 0.0]
        kdict["F"] = [1 - zeta, 1 - zeta, 1 - eta]
        kdict["F_1"] = [zeta, zeta, eta]
        kdict["F_2"] = [-zeta, -zeta, 1 - eta]
        kdict["I"] = [phi, 1 - phi, 0.5]
        kdict["I_1"] = [1 - phi, phi - 1, 0.5]
        kdict["L"] = [0.5, 0.5, 0.5]
        kdict["M"] = [0.5, 0.0, 0.5]
        kdict["X"] = [1 - psi, psi - 1, 0.0]
        kdict["X_1"] = [psi, 1 - psi, 0.0]
        kdict["X_2"] = [psi - 1, -psi, 0.0]
        kdict["Y"] = [0.5, 0.5, 0.0]
        kdict["Y_1"] = [-0.5, -0.5, 0.0]
        kdict["Z"] = [0.0, 0.0, 0.5]

        names = ["Γ", "Y", "F", "L", "I", "Γ", "N", "M", "L"]

    elseif cat == "mclc2"
        
        zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha) ^ 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        psi = 0.75 - a ^ 2 / (4 * b ^ 2 * sin(alpha) ^ 2)
        phi = psi + (0.75 - psi) * b * cos(alpha) / c
        
        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["N"] = [0.5, 0.0, 0.0]
        kdict["N_1"] = [0.0, -0.5, 0.0]
        kdict["F"] = [1 - zeta, 1 - zeta, 1 - eta]
        kdict["F_1"] = [zeta, zeta, eta]
        kdict["F_2"] = [-zeta, -zeta, 1 - eta]
        kdict["F_3"] = [1 - zeta, -zeta, 1 - eta]
        kdict["I"] = [phi, 1 - phi, 0.5]
        kdict["I_1"] = [1 - phi, phi - 1, 0.5]
        kdict["L"] = [0.5, 0.5, 0.5]
        kdict["M"] = [0.5, 0.0, 0.5]
        kdict["X"] = [1 - psi, psi - 1, 0.0]
        kdict["X_1"] = [psi, 1 - psi, 0.0]
        kdict["X_2"] = [psi - 1, -psi, 0.0]
        kdict["Y"] = [0.5, 0.5, 0.0]
        kdict["Y_1"] = [-0.5, -0.5, 0.0]
        kdict["Z"] = [0.0, 0.0, 0.5]
        
        names =  ["Γ", "Y", "F", "L", "I", "N", "Γ", "M"]

    elseif cat == "mclc3"

        mu = (1 + b ^ 2 / a ^ 2) / 4.0
        delta = b * c * cos(alpha) / (2 * a ^ 2)
        zeta = mu - 0.25 + (1 - b * cos(alpha) / c) / (4 * sin(alpha) ^ 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        phi = 1 + zeta - 2 * mu
        psi = eta - 2 * delta

        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["F"] = [1 - phi, 1 - phi, 1 - psi]
        kdict["F_1"] = [phi, phi - 1, psi]
        kdict["F_2"] = [1 - phi, -phi, 1 - psi]
        kdict["H"] = [zeta, zeta, eta]
        kdict["H_1"] = [1 - zeta, -zeta, 1 - eta]
        kdict["H_2"] = [-zeta, -zeta, 1 - eta]
        kdict["I"] = [0.5, -0.5, 0.5]
        kdict["M"] = [0.5, 0.0, 0.5]
        kdict["N"] = [0.5, 0.0, 0.0]
        kdict["N_1"] = [0.0, -0.5, 0.0]
        kdict["X"] = [0.5, -0.5, 0.0]
        kdict["Y"] = [mu, mu, delta]
        kdict["Y_1"] = [1 - mu, -mu, -delta]
        kdict["Y_2"] = [-mu, -mu, -delta]
        kdict["Y_3"] = [mu, mu - 1, delta]
        kdict["Z"] = [0.0, 0.0, 0.5]

        names = ["Γ", "Y", "F", "H", "Z", "I", "F_1", "H_1", "Y_1", "X", "Γ", "N"]

    elseif cat == "mclc4"

        mu = (1 + b ^ 2 / a ^ 2) / 4.0
        delta = b * c * cos(alpha) / (2 * a ^ 2)
        zeta = mu - 0.25 + (1 - b * cos(alpha) / c) / (4 * sin(alpha) ^ 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        phi = 1 + zeta - 2 * mu
        psi = eta - 2 * delta

        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["F"] = [1 - phi, 1 - phi, 1 - psi]
        kdict["F_1"] = [phi, phi - 1, psi]
        kdict["F_2"] = [1 - phi, -phi, 1 - psi]
        kdict["H"] = [zeta, zeta, eta]
        kdict["H_1"] = [1 - zeta, -zeta, 1 - eta]
        kdict["H_2"] = [-zeta, -zeta, 1 - eta]
        kdict["I"] = [0.5, -0.5, 0.5]
        kdict["M"] = [0.5, 0.0, 0.5]
        kdict["N"] = [0.5, 0.0, 0.0]
        kdict["N_1"] = [0.0, -0.5, 0.0]
        kdict["X"] = [0.5, -0.5, 0.0]
        kdict["Y"] = [mu, mu, delta]
        kdict["Y_1"] = [1 - mu, -mu, -delta]
        kdict["Y_2"] = [-mu, -mu, -delta]
        kdict["Y_3"] = [mu, mu - 1, delta]
        kdict["Z"] = [0.0, 0.0, 0.5]

        names = ["Γ", "Y", "F", "H", "Z", "I", "H_1", "Y_1", "X", "Γ", "N"]

    elseif cat == "mclc5"

        zeta = (
            b ^ 2 / a ^ 2 + (1 - b * cos(alpha) / c) / sin(alpha) ^ 2
        ) / 4
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        mu = (
            eta / 2 + b ^ 2 / (4 * a ^ 2) - b * c * cos(alpha) / (2 * a ^ 2)
        )
        nu = 2 * mu - zeta
        rho = 1 - zeta * a ^ 2 / b ^ 2
        omega = (
            (4 * nu - 1 - b ^ 2 * sin(alpha) ^ 2 / a ^ 2)
            * c
            / (2 * b * cos(alpha))
        )

        delta = zeta * c * cos(alpha) / b + omega / 2 - 0.25

        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["F"] = [nu, nu, omega]
        kdict["F_1"] = [1 - nu, 1 - nu, 1 - omega]
        kdict["F_2"] = [nu, nu - 1, omega]
        kdict["H"] = [zeta, zeta, eta]
        kdict["H_1"] = [1 - zeta, -zeta, 1 - eta]
        kdict["H_2"] = [-zeta, -zeta, 1 - eta]
        kdict["I"] = [rho, 1 - rho, 0.5]
        kdict["I_1"] = [1 - rho, rho - 1, 0.5]
        kdict["L"] = [0.5, 0.5, 0.5]
        kdict["M"] = [0.5, 0.0, 0.5]
        kdict["N"] = [0.5, 0.0, 0.0]
        kdict["N_1"] = [0.0, -0.5, 0.0]
        kdict["X"] = [0.5, -0.5, 0.0]
        kdict["Y"] = [mu, mu, delta]
        kdict["Y_1"] = [1 - mu, -mu, -delta]
        kdict["Y_2"] = [-mu, -mu, -delta]
        kdict["Y_3"] = [mu, mu - 1, delta]
        kdict["Z"] = [0.0, 0.0, 0.5]

        names = ["Γ", "Y", "F", "L", "I", "I_1", "Z", "H", "F_1", "H_1", "Y_1", "X", "Γ", "N", "M"]

    elseif cat == "tria"
        
        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["L"] = [0.5, 0.5, 0.0]
        kdict["M"] = [0.0, 0.5, 0.5]
        kdict["N"] = [0.5, 0.0, 0.5]
        kdict["R"] = [0.5, 0.5, 0.5]
        kdict["X"] = [0.5, 0.0, 0.0]
        kdict["Y"] = [0.0, 0.5, 0.0]
        kdict["Z"] = [0.0, 0.0, 0.5]

        names = ["X", "Γ", "Y", "L", "Γ", "Z", "N", "Γ", "M", "R", "Γ"]

    else ## cat == "trib"

        kdict["Γ"] = [0.0, 0.0, 0.0]
        kdict["L"] = [0.5, -0.5, 0.0]
        kdict["M"] = [0.0, 0.0, 0.5]
        kdict["N"] = [-0.5, -0.5, 0.5]
        kdict["R"] = [0.0, -0.5, 0.5]
        kdict["X"] = [0.0, -0.5, 0.0]
        kdict["Y"] = [0.5, 0.0, 0.0]
        kdict["Z"] = [-0.5, 0.0, 0.5]
        
        names = ["X", "Γ", "Y", "L", "Γ", "Z", "N", "Γ", "M", "R", "Γ"],

    end

    kpts = kassemble(kdict, names)

    return kpts, names, c_std


end




function kassemble(kdict, names)

    kpts = []
    for n in names
        push!(kpts, kdict[n])
    end

    return kpts

end




end #end module

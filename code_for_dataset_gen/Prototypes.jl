using ThreeBodyTB
using LinearAlgebra
using DelimitedFiles
using ..ThreeBodyTB.Atomdata:atom_prefered_oxidation
using ..ThreeBodyTB.Atomdata:min_dimer_dist_dict
using ..ThreeBodyTB.Atomdata:sub_list
using ..ThreeBodyTB.Atomdata:electronegativity
using ..ThreeBodyTB.QE:loadXML
using ..ThreeBodyTB.CalcTB:distances_etc_3bdy_parallel
using ..ThreeBodyTB.CalcTB:calc_frontier
using ..ThreeBodyTB.ManageDatabase:prepare_database

ThreeBodyTB.set_units(both="atomic")

struct proto_data

    CalcD::Dict
    core_mono::Array{String}
    core_binary::Array{String}

    A0::Array{String}
    A1B1::Array{String}
    A1B2::Array{String}
    A1B3::Array{String}
    A1B4::Array{String}
    A1B5::Array{String}
    A1B6::Array{String}
    A2B3::Array{String}
    A2B5::Array{String}
    A3B5::Array{String}
    metals::Array{String}
    short_bonds::Array{String}
    core_ternary::Array{String}
    core_mono_mag::Array{String}
    core_mono_mag2::Array{String}

end


function check_frontier(crys)

    prepare_database(crys)
    database = ThreeBodyTB.ManageDatabase.database_cached
    violation_list, vio_bool = calc_frontier(crys, database, test_frontier=true, verbose=false)
    
    return vio_bool
end


function get_twobody_dist(A,B)

    ab = 4.7
    #prepare_database([A,B])
    #database = ThreeBodyTB.ManageDatabase.database_cached
    
    #ab = database[(A,B)].min_dist * 1.0001
    #println("min dist $A $B $ab")
    return ab
end
#=    println("get_twobody_dist $A $B")

dir1 = "/wrk/kfg/julia_data/$A"*"_$B"*"_kspace"
dir2 = "/wrk/kfg/julia_data/$B"*"_$A"*"_kspace"

if isdir(dir1)
dir = dir1
else
dir = dir2
end

println(dir)
ab = 1.0
try
name_t="dimer.in"
ncalc_t = 1
newst_t = "coords"
ab_dir="$dir/$name_t"*"_vnscf_"*"$newst_t"*"_"*"1"        
println("try $ab_dir")
println("dft ",ab_dir*"/qe.save")

dft = loadXML(ab_dir*"/qe.save")
ab = -dft.crys.coords[1,3] * dft.crys.A[3,3] * 2.0
ab2 = (min_dimer_dist_dict[A] + min_dimer_dist_dict[B]) / 2.0 
ab = min(ab,ab2)
println("try $ab")
catch
ab = (min_dimer_dist_dict[A] + min_dimer_dist_dict[B]) / 2.0 
println("catch $ab")
end
println("ab $ab")

return ab
end
=#

function check_twobody_dist(crys)
    R_keep, R_keep_ab, array_ind3, array_floats3, dist_arr, c_zero, dmin_types, dmin_types3 = distances_etc_3bdy_parallel(crys,10.0, 0.0)
    ret = true
    D = Dict()
    for t1 in crys.types
        for t2 in crys.types
            D[(t1,t2)] = 100.0
        end
    end
    for t1 in crys.types
        for t2 in crys.types
            ab = get_twobody_dist(t1,t2)
            D[(t1,t2)] = min(D[(t1,t2)], ab)
        end
    end
    return D
end
#            if dmin_types[Set([t1,t2])] < ab
#                ret = false
#                break
#            end
#        end
#    end
#    return ret
#end



function setup_proto_data()

    STRUCTDIR = ThreeBodyTB.STRUCTDIR

    CalcD = Dict()

    CalcD["big1"] = ["$STRUCTDIR/binary/big1", "scf", "all", "scf", "nscf", false]
    CalcD["big2"] = ["$STRUCTDIR/binary/big2", "scf", "all", "scf", "nscf", false]
    CalcD["big3"] = ["$STRUCTDIR/binary/big3", "scf", "all", "scf", "nscf", false]
    CalcD["big4"] = ["$STRUCTDIR/binary/big4", "scf", "all", "scf", "nscf", false]


    CalcD["atom2_mag"] = ["$STRUCTDIR/atom_small.in", "scf", "all", "scf", "nscf", true]
    CalcD["atom_mag"] = ["$STRUCTDIR/atom.in",        "scf", "all", "scf", "nscf", true]
    CalcD["sc_mag"] = ["$STRUCTDIR/sc.in.up", "vc-relax", "all", "vol-mag", "nscf", true]
    CalcD["bcc_mag"] = ["$STRUCTDIR/bcc.in.up", "vc-relax", "all", "vol-mag", "nscf", true]
    CalcD["fcc_mag"] = ["$STRUCTDIR/fcc.in.up", "vc-relax", "all", "vol-mag", "nscf", true]
    CalcD["line_mag"] = ["$STRUCTDIR/line.in.up", "vc-relax", "z", "vol-mag", "nscf", true]
    CalcD["dimer_mag"] =       ["$STRUCTDIR/dimer.in", "relax", "2Dxy", "coords", "nscf", true]

    CalcD["sc_mag2"] = ["$STRUCTDIR/sc.in.up", "vc-relax", "all", "vol-mag-more", "nscf", true]
    CalcD["bcc_mag2"] = ["$STRUCTDIR/bcc.in.up", "vc-relax", "all", "vol-mag-more", "nscf", true]
    CalcD["fcc_mag2"] = ["$STRUCTDIR/fcc.in.up", "vc-relax", "all", "vol-mag-more", "nscf", true]


    CalcD["rsO_mag"] = ["$STRUCTDIR/binary/rocksaltO.in", "vc-relax", "all", "vol-mag", "nscf", true]
    CalcD["rsF_mag"] = ["$STRUCTDIR/binary/rocksaltF.in", "vc-relax", "all", "vol-mag", "nscf", true]

    CalcD["dimer_fe_mag"] =       ["$STRUCTDIR/dimer.in.fe", "relax", "2Dxy", "coords-small2", "nscf", true]
    CalcD["dimer_mn_mag"] =       ["$STRUCTDIR/dimer.in.mn", "relax", "2Dxy", "coords-small2", "nscf", true]
    CalcD["dimer_f_mag"] =       ["$STRUCTDIR/dimer.in.f", "relax", "2Dxy", "coords-small2", "nscf", true]
    CalcD["dimer_h_mag"] =       ["$STRUCTDIR/dimer.in.h", "relax", "2Dxy", "coords-small2", "nscf", true]

    #    CalcD["dimer_fe_nomag"] =       ["$STRUCTDIR/dimer.in.fe", "relax", "2Dxy", "coords-small2", "nscf", false]
    #    CalcD["dimer_mn_nomag"] =       ["$STRUCTDIR/dimer.in.mn", "relax", "2Dxy", "coords-small2", "nscf", false]
    #    CalcD["dimer_f_nomag"] =       ["$STRUCTDIR/dimer.in.f", "relax", "2Dxy", "coords-small2", "nscf", false]
    #    CalcD["dimer_h_nomag"] =       ["$STRUCTDIR/dimer.in.h", "relax", "2Dxy", "coords-small2", "nscf", false]


    CalcD["b6"] = ["$STRUCTDIR/POSCAR_b6", "vc-relax", "all", "vol", "nscf", false]

    CalcD["atom"] = ["$STRUCTDIR/atom.in", "scf", "all", "scf", "nscf", false]
    CalcD["sc"] = ["$STRUCTDIR/sc.in.up", "vc-relax", "all", "vol-big", "nscf", false]
    CalcD["sc_inv"] = ["$STRUCTDIR/sc.in.up", "vc-relax", "all", "break_inv", "nscf", false]
    CalcD["bcc"] = ["$STRUCTDIR/bcc.in.up", "vc-relax", "all", "vol-mid", "nscf", false]
    CalcD["bcc_inv"] = ["$STRUCTDIR/bcc_atom2.in", "vc-relax", "all", "break_inv", "nscf", false]
    CalcD["fcc"] = ["$STRUCTDIR/fcc.in.up", "vc-relax", "all", "vol-mid", "nscf", false]

    CalcD["line_bin"] = ["$STRUCTDIR/binary/line.in", "vc-relax", "z", "scf", "nscf", false]

    CalcD["line"] = ["$STRUCTDIR/line.in.up", "vc-relax", "z", "vol2", "nscf", false]
    CalcD["line_rumple"] = ["$STRUCTDIR/line.in.rumple.up", "vc-relax", "z", "scf", "nscf", false]


    CalcD["sc_verydense"] = ["$STRUCTDIR/sc.in.up", "vc-relax", "all", "vol-verydense", "nscf", false]

    CalcD["fcc_verydense"] = ["$STRUCTDIR/fcc.in.up", "vc-relax", "all", "vol-verydense", "nscf", false]
    CalcD["bcc_verydense"] = ["$STRUCTDIR/bcc.in.up", "vc-relax", "all", "vol-verydense", "nscf", false]
    CalcD["diamond_verydense"] = ["$STRUCTDIR/diamond.in.up", "vc-relax", "all", "vol-verydense", "nscf", false]


    CalcD["fcc_dense"] = ["$STRUCTDIR/fcc.in.up", "vc-relax", "all", "vol-dense", "nscf", false]
    CalcD["bcc_dense"] = ["$STRUCTDIR/bcc.in.up", "vc-relax", "all", "vol-dense", "nscf", false]

    CalcD["hcp"] = ["$STRUCTDIR/hcp.in.up", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["diamond"] = ["$STRUCTDIR/diamond.in.up", "vc-relax", "all", "vol-mid", "nscf", false]
    CalcD["graphene"] = ["$STRUCTDIR/fake_graphene.in", "vc-relax", "2Dxy", "2D-mid", "nscf", false]
    CalcD["hex"] = ["$STRUCTDIR/hex.in.up", "vc-relax", "2Dxy", "2D-mid", "nscf", false]
    CalcD["hex_short"] = ["$STRUCTDIR/hex.in.up", "vc-relax", "2Dxy", "2D-short", "nscf", false]
    CalcD["square"] = ["$STRUCTDIR/square.in.up", "vc-relax", "2Dxy", "scf", "nscf", false]


    CalcD["el_party"] =       ["$STRUCTDIR/dimer.in.small", "relax", "2Dxy", "el_party", "nscf", false]
    CalcD["bin_party"] =       ["$STRUCTDIR/binary/dimer.in.small", "relax", "2Dxy", "bin_party", "nscf", false]
    CalcD["tern_party"] =       ["$STRUCTDIR/ternary/POSCAR_rand", "relax", "2Dxy", "tern_party", "nscf", false]

    CalcD["fourAgraph"] =       ["$STRUCTDIR/fourAgraph.in", "relax", "2Dxy", "vol-mid", "scf", false]
    CalcD["fourAsquare"] =       ["$STRUCTDIR/fourAsquare.in", "relax", "2Dxy", "vol-mid", "scf", false]
    CalcD["fourAtet"] =       ["$STRUCTDIR/fourAtet.in", "relax", "2Dxy", "vol-mid", "scf", false]


    CalcD["dimer"] =       ["$STRUCTDIR/dimer.in", "relax", "2Dxy", "coords", "nscf", false]
    CalcD["dimer_short"] = ["$STRUCTDIR/dimer.in", "relax", "2Dxy", "coords-short", "nscf", false]
    CalcD["dimer_small"] = ["$STRUCTDIR/dimer_small.in", "relax", "2Dxy", "scf_small", "nscf", false]
    CalcD["dimer_super"] = ["$STRUCTDIR/dimer.in", "relax", "2Dxy", "coords_super", "nscf", false]
    CalcD["dimer2_super"] = ["$STRUCTDIR/binary/dimer.in", "scf", "2Dxy", "coords_super", "nscf", false]
    CalcD["hcp_shape"] = ["$STRUCTDIR/hcp.in.up", "vc-relax", "all", "shape", "nscf", false]

    CalcD["hex_2lay"] = ["$STRUCTDIR/hex_2layers.in.up", "vc-relax", "2Dxy", "2D", "nscf", false]
    CalcD["bcc_2lay"] = ["$STRUCTDIR/bcc_2layers.in.up", "vc-relax", "2Dxy", "2D", "nscf", false]

    CalcD["fcc_huge"] = ["$STRUCTDIR/fcc.in.up", "vc-relax", "all", "vol-huge", "nscf", false]

    CalcD["bcc_tet"] = ["$STRUCTDIR/bcc_tet.in", "vc-relax", "all", "scf", "nscf", false]

    CalcD["trimer"] =       ["$STRUCTDIR/atom_small.in", "none", "2Dxy", "coords_trimer", "nscf", false]
    CalcD["trimerX"] =       ["$STRUCTDIR/dimer.in", "relax", "2Dxy", "coords_trimerX", "scf", false]
    CalcD["trimerY"] =       ["$STRUCTDIR/dimer.in", "relax", "2Dxy", "coords_trimerY", "scf", false]

#    CalcD["trimer"] =       ["$STRUCTDIR/trimer.in", "none", "2Dxy", "coords_trimer", "nscf", false]
    CalcD["trimer_dense"] =       ["$STRUCTDIR/trimer.in", "none", "2Dxy", "coords_trimer_dense", "nscf", false]
    CalcD["trimer2"] =       ["$STRUCTDIR/trimer.in2", "none", "2Dxy", "coords_trimer2", "nscf", false]
    CalcD["trimer3"] =       ["$STRUCTDIR/trimer.in3", "none", "2Dxy", "coords_trimer3", "nscf", false]

    CalcD["trimer_ab2"] =       ["$STRUCTDIR/binary/trimer.in.ab2", "none", "2Dxy", "coords_trimer_ab", "nscf", false]
    CalcD["trimer2_ab2"] =       ["$STRUCTDIR/binary/trimer.in2.ab2", "none", "2Dxy", "coords_trimer_ab", "nscf", false]

    CalcD["trimer_ba2"] =       ["$STRUCTDIR/binary/trimer.in.ba2", "none", "2Dxy", "coords_trimer_ab", "nscf", false]
    CalcD["trimer2_ba2"] =       ["$STRUCTDIR/binary/trimer.in2.ba2", "none", "2Dxy", "coords_trimer_ab", "nscf", false]


    CalcD["trimer_ab_new"] =       ["$STRUCTDIR/binary/trimer.in.ab2", "none", "2Dxy", "coords_trimer_ab_new", "nscf", false]
    CalcD["trimer_a_new"] =       ["$STRUCTDIR/binary/trimer.in.ab2", "none", "2Dxy", "coords_trimer_a_new", "nscf", false]
    

    CalcD["trimer_ab2_dense"] =       ["$STRUCTDIR/binary/trimer.in.ab2", "none", "2Dxy", "coords_trimer_ab_dense", "nscf", false]
    CalcD["trimer2_ab2_dense"] =       ["$STRUCTDIR/binary/trimer.in2.ab2", "none", "2Dxy", "coords_trimer_ab_dense", "nscf", false]

    CalcD["trimer_ba2_dense"] =       ["$STRUCTDIR/binary/trimer.in.ba2", "none", "2Dxy", "coords_trimer_ab_dense", "nscf", false]
    CalcD["trimer2_ba2_dense"] =       ["$STRUCTDIR/binary/trimer.in2.ba2", "none", "2Dxy", "coords_trimer_ab_dense", "nscf", false]


    CalcD["trimer_ba2_big"] =       ["$STRUCTDIR/binary/trimer.in.ba2.big", "none", "2Dxy", "coords_trimer_ab_big", "nscf", false]


    CalcD["as_221"] = ["$STRUCTDIR/POSCAR_As_221", "vc-relax", "all", "scf", "nscf", false]
    CalcD["as_orth"] = ["$STRUCTDIR/POSCAR_As_ortho", "vc-relax", "all", "scf", "nscf", false]
    CalcD["ga_tet"] = ["$STRUCTDIR/POSCAR_ga_tet", "vc-relax", "all", "scf", "nscf", false]
    CalcD["ge_wurtz"] = ["$STRUCTDIR/POSCAR_ge_wurtz", "vc-relax", "all", "scf", "nscf", false]
    CalcD["pb_r3m"] = ["$STRUCTDIR/POSCAR_pb_r3m_2atom", "vc-relax", "all", "scf", "nscf", false]
    CalcD["beta_sn"] = ["$STRUCTDIR/beta_sn.in", "vc-relax", "all", "scf", "nscf", false]
    CalcD["n"] = ["$STRUCTDIR/n.in", "vc-relax", "all", "scf", "nscf", false]
    CalcD["f"] = ["$STRUCTDIR/f.in", "vc-relax", "all", "scf", "nscf", false]
    CalcD["ga"] = ["$STRUCTDIR/POSCAR_ga", "vc-relax", "all", "scf", "nscf", false]
    CalcD["bi"] = ["$STRUCTDIR/bi.in", "vc-relax", "all", "scf", "nscf", false]
    CalcD["te"] = ["$STRUCTDIR/te.in", "vc-relax", "all", "scf", "nscf", false]
    CalcD["in"] = ["$STRUCTDIR/in.in", "vc-relax", "all", "scf", "nscf", false]
    CalcD["i2"] = ["$STRUCTDIR/i2.in", "vc-relax", "all", "scf", "nscf", false]
    CalcD["li_p6mmm"] = ["$STRUCTDIR/POSCAR_li_p6mmm", "vc-relax", "all", "scf", "nscf", false]

    CalcD["sc_shape"] = ["$STRUCTDIR/sc.in.up", "vc-relax", "all", "shape", "nscf", false]
    CalcD["diamond_shear"] = ["$STRUCTDIR/diamond.in.up", "vc-relax", "all", "shear", "nscf", false]

    CalcD["hh_mono"] = ["$STRUCTDIR/POSCAR_hh_mono", "vc-relax", "all", "vol-mid", "nscf", false]


    #    CalcD["cscl"] = ["$STRUCTDIR/binary/cscl.in", "vc-relax", "all", "vol-big", "nscf", false]
    CalcD["cscl"] = ["$STRUCTDIR/binary/cscl.in", "vc-relax", "all", "scf", "nscf", false]
    CalcD["cscl_layers"] = ["$STRUCTDIR/binary/cscl_layers.in", "vc-relax", "all", "scf", "nscf", false]
    CalcD["cscl_inv"] = ["$STRUCTDIR/binary/cscl.in", "vc-relax", "all", "break_inv", "nscf", false]
    CalcD["rocksalt"] = ["$STRUCTDIR/binary/rocksalt.in", "vc-relax", "all", "vol-mid", "nscf", false]
    CalcD["rocksalt_inv"] = ["$STRUCTDIR/binary/rocksalt.in", "vc-relax", "all", "break_inv", "nscf", false]

    CalcD["rocksalt_dense"] = ["$STRUCTDIR/binary/rocksalt.in", "vc-relax", "all", "vol-dense", "nscf", false]
    CalcD["hcp_v2_dense"] = ["$STRUCTDIR/binary/hcp.in2", "vc-relax", "all", "vol-dense", "nscf", false]
    CalcD["znse_dense"] = ["$STRUCTDIR/binary/znse.in", "vc-relax", "all", "vol-dense", "nscf", false]

    CalcD["nias"] = ["$STRUCTDIR/binary/POSCAR_nias_hex", "vc-relax", "all", "scf", "nscf", false]
    CalcD["227"] = ["$STRUCTDIR/binary/POSCAR_227", "vc-relax", "all", "scf", "nscf", false]
    CalcD["laga"] = ["$STRUCTDIR/binary/POSCAR_laga_63", "vc-relax", "all", "scf", "nscf", false]

    CalcD["distort"] = ["$STRUCTDIR/binary/POSCAR_distort", "vc-relax", "all", "scf", "nscf", false]
    CalcD["distort_ab2"] = ["$STRUCTDIR/binary/POSCAR_distort_ab2", "vc-relax", "all", "scf", "nscf", false]
    CalcD["distort_ba2"] = ["$STRUCTDIR/binary/POSCAR_distort_ba2", "vc-relax", "all", "scf", "nscf", false]

    CalcD["h2_atom"] = ["$STRUCTDIR/binary/h2_atom.in", "vc-relax", "all", "scf", "nscf", false]
    CalcD["h1_atom"] = ["$STRUCTDIR/binary/h1_atom.in", "vc-relax", "all", "scf", "nscf", false]

    ##    CalcD["227_BA"] = ["$STRUCTDIR/binary/POSCAR_227_BA", "vc-relax", "all", "vol-mid", "nscf", false]


    CalcD["znse_shear"] = ["$STRUCTDIR/binary/znse.in", "vc-relax", "all", "shear", "nscf", false]
    CalcD["hbn"] =  ["$STRUCTDIR/binary/hbn.in", "vc-relax", "2Dxy", "scf", "nscf", false]
    CalcD["znse"] = ["$STRUCTDIR/binary/znse.in", "vc-relax", "all", "vol-mid", "nscf", false]

    CalcD["dimer2"] = ["$STRUCTDIR/binary/dimer2.in", "relax", "all", "coords", "nscf", false]
#    CalcD["dimer2"] = ["$STRUCTDIR/binary/dimer2.in", "relax", "all", "scf", "nscf", false]

    CalcD["dimer2_min"] = ["$STRUCTDIR/binary/dimer.in", "none", "all", "coords_min", "nscf", false]
    CalcD["tri2_min"] = ["$STRUCTDIR/binary/dimer.in", "none", "all", "coords_min_tri", "nscf", false]
    CalcD["dimer_min"] = ["$STRUCTDIR/dimer.in", "none", "all", "coords_min", "nscf", false]
    CalcD["tri_min"] = ["$STRUCTDIR/binary/dimer.in", "none", "all", "coords_min_tri", "nscf", false]


    CalcD["dimer2_rev"] = ["$STRUCTDIR/binary/dimer_rev.in", "relax", "all", "coords", "nscf", false]
    CalcD["square2"] = ["$STRUCTDIR/binary/square.in", "vc-relax", "2Dxy", "scf", "nscf", false]
    CalcD["caf2"] = ["$STRUCTDIR/binary/POSCAR_caf2", "vc-relax", "all", "vol2", "nscf", false]

    CalcD["co2"] = ["$STRUCTDIR/binary/co2.in", "relax", "all", "coords-small", "nscf", false]
    CalcD["co2_v2"] = ["$STRUCTDIR/binary/co2_v2.in", "relax", "all", "coords-small", "nscf", false]

    CalcD["square_ab2"] = ["$STRUCTDIR/binary/square_ab2.in", "vc-relax", "2Dxy", "2D", "nscf", false]
    CalcD["mgf2"] = ["$STRUCTDIR/binary/POSCAR_mgf2", "vc-relax", "all", "scf", "nscf", false]
    CalcD["hcp_v2"] = ["$STRUCTDIR/binary/hcp.in2", "vc-relax", "all", "vol", "nscf", false]

    CalcD["p2ca3"] = ["$STRUCTDIR/binary/POSCAR_p2ca3", "vc-relax", "all", "scf", "nscf", false]

    CalcD["triangle"] = ["$STRUCTDIR/binary/triangle_r.in", "vc-relax", "all", "scf", "nscf", false]
    CalcD["triangle2"] = ["$STRUCTDIR/binary/triangle_r2.in", "vc-relax", "all", "scf", "nscf", false]


    CalcD["beta_sn2"] = ["$STRUCTDIR/binary/beta_sn.in.up", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["hbn_real"] = ["$STRUCTDIR/binary/hbn_real.in", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["znseAAAB"] = ["$STRUCTDIR/binary/znse.in.super.AAAB", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["rocksaltAAAB"] = ["$STRUCTDIR/binary/rocksalt.in.super.AAAB", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["znseABBB"] = ["$STRUCTDIR/binary/znse.in.super.ABBB", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["rocksaltABBB"] = ["$STRUCTDIR/binary/rocksalt.in.super.ABBB", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["al2o3"] = ["$STRUCTDIR/binary/POSCAR_al2o3", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["sis2"] = ["$STRUCTDIR/binary/POSCAR_sis2", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["tio2_rutile"] = ["$STRUCTDIR/binary/POSCAR_tio2_rutile", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["bi2se3"] = ["$STRUCTDIR/binary/POSCAR_bi2se3", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["sis_2d"] = ["$STRUCTDIR/binary/sis_2d.in", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["sio2_224"] = ["$STRUCTDIR/binary/sio2.224.in", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["ges"] = ["$STRUCTDIR/binary/POSCAR_ges", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["alf3"] = ["$STRUCTDIR/binary/POSCAR_alf3", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["lipd3"] = ["$STRUCTDIR/binary/POSCAR_lipd3", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["mgf2_v2"] = ["$STRUCTDIR/binary/POSCAR_mgf2_v2", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["sns"] = ["$STRUCTDIR/binary/POSCAR_sns", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["sns2"] = ["$STRUCTDIR/binary/POSCAR_sns2", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["y2o3"] = ["$STRUCTDIR/binary/POSCAR_y2o3", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["wurtz"] = ["$STRUCTDIR/binary/POSCAR_wurtz", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["mgcl2"] = ["$STRUCTDIR/binary/POSCAR_mgcl2", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["mgcl2_tet"] = ["$STRUCTDIR/binary/POSCAR_mgcl2_tet", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["asna3_2d"] = ["$STRUCTDIR/binary/POSCAR_asna3", "vc-relax", "2Dxy", "vol2", "nscf", false]
    CalcD["gain3"] = ["$STRUCTDIR/binary/POSCAR_gain3", "vc-relax", "all", "vol2", "nscf", false]

    CalcD["li3n_hex"] = ["$STRUCTDIR/binary/POSCAR_li3n", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["nan3"] = ["$STRUCTDIR/binary/POSCAR_nan3", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["rbo2"] = ["$STRUCTDIR/binary/POSCAR_rbo2", "vc-relax", "all", "vol2", "nscf", false]



    CalcD["nb2o5"] = ["$STRUCTDIR/binary/POSCAR_nb2o5", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["sif4"] = ["$STRUCTDIR/binary/POSCAR_sif4", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["ticl2"] = ["$STRUCTDIR/binary/POSCAR_ticl2", "vc-relax", "all", "vol2", "nscf", false]

    CalcD["snf4"] = ["$STRUCTDIR/binary/POSCAR_snf4", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["ga2s3"] = ["$STRUCTDIR/binary/POSCAR_ga2s3", "vc-relax", "all", "vol2", "nscf", false]

    CalcD["bcc_13"] = ["$STRUCTDIR/binary/POSCAR_bcc_13", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["bcc_31"] = ["$STRUCTDIR/binary/POSCAR_bcc_31", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["mg2si"] = ["$STRUCTDIR/binary/POSCAR_mg2si", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["simg2"] = ["$STRUCTDIR/binary/POSCAR_simg2", "vc-relax", "all", "vol2", "nscf", false]

    CalcD["fcc_12"] = ["$STRUCTDIR/binary/fcc_12.in", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["fcc_21"] = ["$STRUCTDIR/binary/fcc_21.in", "vc-relax", "all", "vol2", "nscf", false]

    CalcD["fcc_conv_ABBB"] = ["$STRUCTDIR/binary/fcc_conv_ABBB.in", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["fcc_conv_BAAA"] = ["$STRUCTDIR/binary/fcc_conv_BAAA.in", "vc-relax", "all", "vol2", "nscf", false]

    CalcD["mgb2_12"] = ["$STRUCTDIR/binary/POSCAR_mgb2", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["mgb2_21"] = ["$STRUCTDIR/binary/POSCAR_mgb2_21", "vc-relax", "all", "vol2", "nscf", false]

    CalcD["mgb2_AB"] = ["$STRUCTDIR/binary/POSCAR_mgb2_AB", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["mgb2_BA"] = ["$STRUCTDIR/binary/POSCAR_mgb2_BA", "vc-relax", "all", "vol2", "nscf", false]

    CalcD["ab2_71"] = ["$STRUCTDIR/binary/POSCAR_ab2_71", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["ba2_71"] = ["$STRUCTDIR/binary/POSCAR_ba2_71", "vc-relax", "all", "vol2", "nscf", false]

    CalcD["irn2_38"] = ["$STRUCTDIR/binary/POSCAR_irn2_38", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["cao2_12"] = ["$STRUCTDIR/binary/POSCAR_cao2_12", "vc-relax", "all", "vol2", "nscf", false]

    CalcD["p2o5"] = ["$STRUCTDIR/binary/POSCAR_p2o5", "vc-relax", "all", "vol2", "nscf", false]


    CalcD["znseAABB"] = ["$STRUCTDIR/binary/znse.AABB.in", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["squareAABB"] = ["$STRUCTDIR/binary/square.in.AABB", "vc-relax", "2Dxy", "2D", "nscf", false]


    CalcD["dimer_pair"] = ["$STRUCTDIR/binary/POSCAR_dimer_pair", "vc-relax", "all", "vol2", "nscf", false]


    CalcD["mgb2_AB-mid"] = ["$STRUCTDIR/binary/POSCAR_mgb2_AB", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["mgb2_BA-mid"] = ["$STRUCTDIR/binary/POSCAR_mgb2_BA", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["mgb2_12-mid"] = ["$STRUCTDIR/binary/POSCAR_mgb2", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["mgb2_21-mid"] = ["$STRUCTDIR/binary/POSCAR_mgb2_21", "vc-relax", "all", "vol2", "nscf", false]


    #    CalcD["i4mmm_tet_AB"] = ["$STRUCTDIR/binary/POSCAR_i4mmm_4atom_AB", "vc-relax", "all", "vol2", "nscf", false]
    #    CalcD["i4mmm_tet_BA"] = ["$STRUCTDIR/binary/POSCAR_i4mmm_4atom_BA", "vc-relax", "all", "vol2", "nscf", false]



    CalcD["gei2"] = ["$STRUCTDIR/binary/POSCAR_GeI2", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["gei4_mol"] = ["$STRUCTDIR/binary/gei4_molecule.in", "relax", "all", "vol2", "nscf", false]

    CalcD["ab3_mol"] = ["$STRUCTDIR/binary/ab3_molecule.in", "relax", "all", "vol2", "nscf", false]


    CalcD["tet"] = ["$STRUCTDIR/binary/tet.in", "vc-relax", "all", "vol2", "nscf", false]


    CalcD["mof6"] = ["$STRUCTDIR/binary/POSCAR_mof6", "vc-relax",  "all", "vol2", "nscf", false]


    CalcD["ascl5"] = ["$STRUCTDIR/binary/POSCAR_ascl5", "vc-relax",  "2Dxy", "2D", "nscf", false]
    CalcD["rocksalt_shape"] = ["$STRUCTDIR/binary/rocksalt.in", "vc-relax", "all", "shape", "nscf", false]
    CalcD["rocksalt_2lay"] = ["$STRUCTDIR/binary/rocksalt.in.2lay", "vc-relax", "2Dxy", "2D", "nscf", false]

    CalcD["gas"] = ["$STRUCTDIR/binary/POSCAR_gas", "vc-relax", "all", "vol2", "nscf", false]

    CalcD["hex12"] = ["$STRUCTDIR/binary/hex_trim_12.in", "vc-relax",  "2Dxy", "2D", "nscf", false]
    CalcD["hex21"] = ["$STRUCTDIR/binary/hex_trim_21.in", "vc-relax",  "2Dxy", "2D", "nscf", false]

    CalcD["hex12a"] = ["$STRUCTDIR/binary/hex_trim_12.in", "vc-relax",  "2Dxy", "2D-mid", "nscf", false]
    CalcD["hex21a"] = ["$STRUCTDIR/binary/hex_trim_21.in", "vc-relax",  "2Dxy", "2D-mid", "nscf", false]

    #
    CalcD["sn2o2"] = ["$STRUCTDIR/binary/POSCAR_sn2o2", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["te2o6"] = ["$STRUCTDIR/binary/POSCAR_te2o6", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["a1f6"] = ["$STRUCTDIR/binary/POSCAR_a1f6", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["a2f6"] = ["$STRUCTDIR/binary/POSCAR_a2f6", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["ca1pd5"] = ["$STRUCTDIR/binary/POSCAR_ca1pd5", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["anatase"] = ["$STRUCTDIR/binary/POSCAR_anatase", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["ga4sr"] = ["$STRUCTDIR/binary/POSCAR_ga4sr", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["bi1f5"] = ["$STRUCTDIR/binary/POSCAR_bi1f5", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["cu2o1"] = ["$STRUCTDIR/binary/POSCAR_cu2o1", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["nb1f4"] = ["$STRUCTDIR/binary/POSCAR_nb1f4", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["rb1in4"] = ["$STRUCTDIR/binary/POSCAR_rb1in4", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["ni4la2"] = ["$STRUCTDIR/binary/POSCAR_ni4la2", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["re3n3"] = ["$STRUCTDIR/binary/POSCAR_re3n3", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["re3n"] = ["$STRUCTDIR/binary/POSCAR_re3n", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["cab6"] = ["$STRUCTDIR/binary/POSCAR_cab6", "vc-relax", "all", "vol2", "nscf", false]

    CalcD["c2ca"] = ["$STRUCTDIR/binary/POSCAR_c2ca", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["osn2"] = ["$STRUCTDIR/binary/POSCAR_osn2", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["pt2o2"] = ["$STRUCTDIR/binary/POSCAR_pt2o2", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["p3n5"] = ["$STRUCTDIR/binary/POSCAR_p3n5", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["pd3te"] = ["$STRUCTDIR/binary/POSCAR_pd3te", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["bpt"] = ["$STRUCTDIR/binary/POSCAR_bpt", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["pd3s"] = ["$STRUCTDIR/binary/POSCAR_JVASP-pd3s", "vc-relax", "all", "vol2", "nscf", false]


    CalcD["bcc_5lay"] = ["$STRUCTDIR/binary/POSCAR_bcc_5lay", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["fcc_5lay"] = ["$STRUCTDIR/binary/POSCAR_fcc_5lay", "vc-relax", "all", "vol2", "nscf", false]

    CalcD["bcc_4lay"] = ["$STRUCTDIR/binary/POSCAR_bcc_4lay", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["fcc_4lay"] = ["$STRUCTDIR/binary/POSCAR_fcc_4lay", "vc-relax", "all", "vol2", "nscf", false]

#    CalcD["bcc_3lay"] = ["$STRUCTDIR/binary/POSCAR_bcc_3lay", "vc-relax", "all", "vol2", "nscf", false]
#    CalcD["fcc_3lay"] = ["$STRUCTDIR/binary/POSCAR_fcc_3lay", "vc-relax", "all", "vol2", "nscf", false]


    CalcD["k2n6"] = ["$STRUCTDIR/binary/POSCAR_k2n6", "vc-relax", "all", "vol2", "nscf", false]

    CalcD["mol_h2o"] =    ["$STRUCTDIR/binary/POSCAR_mol_h2o", "relax", "all", "coords-small2", "nscf", false]
    CalcD["mol_h2o_v1"] = ["$STRUCTDIR/binary/POSCAR_mol_h2o_v1", "relax", "all", "coords-small2", "nscf", false]
    CalcD["mol_h2o_v2"] = ["$STRUCTDIR/binary/POSCAR_mol_h2o_v2", "relax", "all", "coords-small2", "nscf", false]

    CalcD["abb_line"] = ["$STRUCTDIR/binary/abb_line.in", "relax", "all", "coords-small2", "nscf", false]
    CalcD["baa_line"] = ["$STRUCTDIR/binary/baa_line.in", "relax", "all", "coords-small2", "nscf", false]

    CalcD["abb_tri"] = ["$STRUCTDIR/binary/abb_tri.in", "relax", "all", "coords-small2", "nscf", false]
    CalcD["baa_tri"] = ["$STRUCTDIR/binary/baa_tri.in", "relax", "all", "coords-small2", "nscf", false]

    CalcD["quad"] = ["$STRUCTDIR/binary/quad.in", "relax", "all", "vol2", "nscf", false]


    #ternary
    CalcD["abc_line"] = ["$STRUCTDIR/ternary/abc_line.in", "relax", "all", "coords-small2", "nscf", false]
    CalcD["bac_line"] = ["$STRUCTDIR/ternary/bac_line.in", "relax", "all", "coords-small2", "nscf", false]
    CalcD["cab_line"] = ["$STRUCTDIR/ternary/cab_line.in", "relax", "all", "coords-small2", "nscf", false]

    CalcD["fcc_tern"] = ["$STRUCTDIR/ternary/fcc_tern.in", "vc-relax", "all", "vol-mid", "nscf", false]
    CalcD["hex_trim"] = ["$STRUCTDIR/ternary/hex_trim_3.in", "vc-relax",  "2Dxy", "2D", "nscf", false]

    CalcD["hh1"] = ["$STRUCTDIR/ternary/POSCAR_hh1", "vc-relax", "all", "vol-mid", "nscf", false]
    CalcD["hh2"] = ["$STRUCTDIR/ternary/POSCAR_hh2", "vc-relax", "all", "vol-mid", "nscf", false]
    CalcD["hh3"] = ["$STRUCTDIR/ternary/POSCAR_hh3", "vc-relax", "all", "vol-mid", "nscf", false]

    CalcD["stuffhex_1"] = ["$STRUCTDIR/ternary/POSCAR_mgb2_1", "vc-relax", "all", "vol-mid", "nscf", false]
    CalcD["stuffhex_2"] = ["$STRUCTDIR/ternary/POSCAR_mgb2_2", "vc-relax", "all", "vol-mid", "nscf", false]
    CalcD["stuffhex_3"] = ["$STRUCTDIR/ternary/POSCAR_mgb2_3", "vc-relax", "all", "vol-mid", "nscf", false]

    CalcD["stuffhex_z_1"] = ["$STRUCTDIR/ternary/POSCAR_z_mgb2_1", "vc-relax", "all", "vol-mid", "nscf", false]
    CalcD["stuffhex_z_2"] = ["$STRUCTDIR/ternary/POSCAR_z_mgb2_2", "vc-relax", "all", "vol-mid", "nscf", false]
    CalcD["stuffhex_z_3"] = ["$STRUCTDIR/ternary/POSCAR_z_mgb2_3", "vc-relax", "all", "vol-mid", "nscf", false]

    CalcD["rocksalt_2lay_abo2"] = ["$STRUCTDIR/ternary/rocksalt.in.2lay_abo2", "vc-relax", "2Dxy", "2D_tern", "nscf", false]
    CalcD["p4mmm"] = ["$STRUCTDIR/ternary/POSCAR_p4mmm", "vc-relax", "all", "3D_tern", "nscf", false]
    CalcD["caf2_abc"] = ["$STRUCTDIR/ternary/POSCAR_caf2_ABC", "vc-relax", "all", "vol2", "nscf", false]

    CalcD["perov"] = ["$STRUCTDIR/ternary/POSCAR_abo3", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["perov2"] = ["$STRUCTDIR/ternary/POSCAR_abo3_2", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["perov3"] = ["$STRUCTDIR/ternary/POSCAR_abo3_3", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["perov4"] = ["$STRUCTDIR/ternary/POSCAR_abo3_4", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["perov5"] = ["$STRUCTDIR/ternary/POSCAR_abo3_5", "vc-relax", "all", "vol2", "nscf", false]
    CalcD["perov6"] = ["$STRUCTDIR/ternary/POSCAR_abo3_6", "vc-relax", "all", "vol2", "nscf", false]


    CalcD["simple_hex"] = ["$STRUCTDIR/simple_hex.in", "vc-relax", "all", "flyaway", "nscf", false]

    CalcD["trimer_rand"] =       ["$STRUCTDIR/ternary/POSCAR_rand", "relax", "all", "vol2", "nscf", false]


    CalcD["trimer_tern"] =       ["$STRUCTDIR/ternary/POSCAR_trimer_tern", "none", "2Dxy", "trimer_tern", "nscf", false]
    CalcD["trimer_tern_right"] =       ["$STRUCTDIR/ternary/POSCAR_trimer_tern_right", "none", "2Dxy", "trimer_tern_right", "nscf", false]
    CalcD["trimer_tern_angle"] =       ["$STRUCTDIR/ternary/POSCAR_trimer_tern_angle", "none", "2Dxy", "trimer_tern_angle", "nscf", false]

    CalcD["trimer_tern_line"] =       ["$STRUCTDIR/ternary/POSCAR_trimer_tern_line", "none", "2Dxy", "trimer_tern_line", "nscf", false]


    CalcD["hh_oxygen"] = ["$STRUCTDIR/ternary/POSCAR_hh_oxygen", "vc-relax", "all", "vol-oxygen", "nscf", false]
    CalcD["hh_oxygen2"] = ["$STRUCTDIR/ternary/POSCAR_hh_oxygen2", "vc-relax", "all", "vol-oxygen", "nscf", false]
    CalcD["hex_oxygen"] = ["$STRUCTDIR/ternary/hex_trim_3.in_oxygen", "vc-relax",  "2Dxy", "2D-oxygen", "nscf", false]

    CalcD["sipau"] = ["$STRUCTDIR/ternary/POSCAR_SiPAu", "vc-relax",  "all", "vol2", "nscf", false]

    CalcD[""] = ["$STRUCTDIR/ternary/hex_trim_3.in_oxygen", "vc-relax",  "2Dxy", "2D-oxygen", "nscf", false]







    core_mono_mag = [   "atom_mag", "atom2_mag", "sc_mag" ,"bcc_mag","fcc_mag","line_mag" ,"dimer_mag"]
    core_mono_mag2 = [   "bcc_mag2","fcc_mag2", "sc_mag2"]



    core_mono = [     "sc", "atom",     "bcc",     "bcc_inv",     "fcc",     "hcp",  "diamond",     "graphene",     "hex",     "square",     "dimer" ,"tri_min", "bcc_5lay", "fcc_5lay", "bcc_4lay", "fcc_4lay", "hex_2lay", "bcc_2lay", "fcc_dense", "bcc_dense", "znse_dense", "line"]




    #core_binary = [   "dimer2",   "cscl",          "hbn",     "rocksalt",  "znse",      "dimer2_min", "tri2_min", "square2", "hcp_v2", "rocksalt_2lay", "cscl_layers", "fcc_12", "fcc_21", "bcc_13", "bcc_31", "rocksalt_dense", "znse_dense", "znseAABB"]

    core_binary = [ "rocksalt", "cscl", "znse",  "hbn", "square2", "rocksalt_2lay", "hcp_v2", "dimer2", "line_bin"]

    A0 = [   "as_orth",    "ga_tet",    "ge_wurtz",    "pb_r3m",    "beta_sn",   "ga",    "bi",    "te",    "in",    "i2",    "li_p6mmm",    "hcp_shape",   "n" ] #"bcc_tet.in",  same as POSCAR_ga_tet   #"as_221",  is simple cubic
    
    #A1B1 = ["hbn_real", "sis_2d", "ges", "sns"] #nias #
    A1B1 = ["ges", "sns"] #nias #
    A1B2 = ["caf2", "tio2_rutile", "anatase"]   # "sns2" duplictes ticl2 #"mgf2_v2" is rutile
#    A1B3 = ["alf3", "asna3_2d", "ab3_mol", "gain3"]
    A1B3 = ["alf3"]
    #A1B4 = ["sif4", "snf4", "gei4_mol"]
    A1B4 = ["sif4"]
#    A1B5 = ["ascl5", "bi1f5"]
    A1B5 = [ "bi1f5"]
#    A1B6 = ["mof6", "a1f6"]
    A1B6 = [ "a1f6"]
#    A2B3 = ["p2ca3", "al2o3", "bi2se3", "ga2s3"]
    A2B3 = ["al2o3"]
#    A2B5 = ["nb2o5", "p2o5"]
    A2B5 = ["nb2o5"]
    A3B5 = [] #"p3n5"

    short_bonds = ["dimer_pair", "cao2_12", "irn2_38"]
    

    #    metals = [ "mg2si", "simg2",  "mgb2_12", "mgb2_21", "fcc_conv_ABBB","fcc_conv_BAAA", "ab2_71", "ba2_71"]
    #    metals = [ "mg2si", "simg2",  "mgb2_12", "mgb2_21",  "fcc_conv_ABBB","fcc_conv_BAAA",  "ab2_71", "ba2_71"]
    #metals = [ "mg2si", "simg2",  "mgb2_12", "mgb2_21",    "ab2_71", "ba2_71"]
    metals = [   "ab2_71", "ba2_71"]

    # all_ternary = ["abc_line", "bac_line", "cab_line", "fcc_tern", "hex_trim", "hh1", "hh2", "hh3", "stuffhex_1", "stuffhex_2", "stuffhex_3","stuffhex_z_1", "stuffhex_z_2", "stuffhex_z_3", "rocksalt_2lay_abo2", "caf2_abc", "perov", "perov2", "perov3",  "perov4",  "perov5",  "perov6"  ]
#    core_ternary = ["abc_line", "bac_line", "cab_line", "fcc_tern", "hex_trim", "hh1", "hh2", "hh3", "stuffhex_1", "stuffhex_2", "stuffhex_3", "trimer_tern","trimer_tern_right","trimer_tern_line", "trimer_tern_angle", "p4mmm" ]
    core_ternary = ["abc_line", "bac_line", "cab_line", "fcc_tern", "hex_trim", "trimer_rand", "sipau"]

    pd = proto_data(CalcD, core_mono, core_binary, A0, A1B1, A1B2, A1B3, A1B4, A1B5, A1B6, A2B3, A2B5,A3B5, metals, short_bonds, core_ternary, core_mono_mag, core_mono_mag2)

    return pd

end

function  do_run(pd, T1, T2, T3, tmpname, dir, procs, torun; nscf_only = false, only_kspace = false, check_only = false, min_nscf = false, dft_only=false)

    TORUN = []
    for t in torun
        if t == :core_mono
            TORUN = [TORUN;pd.core_mono]
        elseif t == :core_binary
            TORUN = [TORUN;pd.core_binary]
        elseif t == :A0
            TORUN = [TORUN;pd.A0]
        elseif t == :A1B1
            TORUN = [TORUN;pd.A1B1]
        elseif t == :A1B2
            TORUN = [TORUN;pd.A1B2]
        elseif t == :A1B3
            TORUN = [TORUN;pd.A1B3]
        elseif t == :A1B4
            TORUN = [TORUN;pd.A1B4]
        elseif t == :A1B5
            TORUN = [TORUN;pd.A1B5]
        elseif t == :A1B6
            TORUN = [TORUN;pd.A1B6]
        elseif t == :A2B3
            TORUN = [TORUN;pd.A2B3]
        elseif t == :A2B5
            TORUN = [TORUN;pd.A2B5]
        elseif t == :A3B5
            TORUN = [TORUN;pd.A3B5]
        elseif t == :metals
            TORUN = [TORUN;pd.metals]
        elseif t == :short_bonds
            TORUN = [TORUN;pd.short_bonds]
        elseif t == :core_ternary
            TORUN = [TORUN;pd.core_ternary]
        elseif t == :core_mono_mag
            TORUN = [TORUN;pd.core_mono_mag]
        elseif t == :core_mono_mag2
            TORUN = [TORUN;pd.core_mono_mag2]
        else
            TORUN = [TORUN; t]
        end
    end

    println("TORUN")
    for t in TORUN
        println(t)
    end
    println("--")
    #    sleep(1)

    already_done = []
    not_done = []

    for st in TORUN

        println("start ", st)

        file, scf, free, newst, calc_mode, magnetic = pd.CalcD[st]

        if nscf_only 
            if calc_mode != "nscf"
                println("nscf_only $nscf_only skip $st")
                continue
            end
        end

        arr = split(file, '/')
        name = arr[end]
        #        println("st $st $scf")

        #############
        #check if already done

        if newst == "vol"
            ncalc = length([ 0.95 1.0 1.05])
        elseif newst == "vol2"
            ncalc = length([ 0.94 1.0 ])
        elseif newst == "vol-mid"
            ncalc = length( [ 0.9 0.95 1.0 1.05 1.1 ])
        elseif newst == "vol-mag"
            ncalc = length( [ 1.0 ])
        elseif newst == "vol-mag-more"
            ncalc = length( [0.9 0.95  1.05 1.1])
        elseif newst == "vol-oxygen"
            ncalc = 4
        elseif newst == "2D-oxygen"
            ncalc = 4
        elseif newst == "vol-dense"
            ncalc = length( [ 0.87 ])
        elseif newst == "vol-verydense"
            ncalc = length( [0.77 0.82 0.75 0.70 0.65 0.60 0.55 0.50])
        elseif newst == "vol-big"
            ncalc = length( [0.80 0.85 0.9 0.95 1.0 1.05 1.1 1.2 1.3 1.5 ])
        elseif newst == "vol-huge"
            ncalc = length( [0.9 0.95 1.0 1.05 1.1 1.2 1.3 1.5 2.0 2.5 3.0 3.5 4.0 5.0])
        elseif newst == "2D_tern"
            ncalc = length( [0.90 0.95 1.0 1.05 ])
        elseif newst == "3D_tern"
            ncalc = length( [ 0.95 1.0 1.05 ])
        elseif newst == "2D"
            #ncalc = length( [0.90 0.95 1.0 1.05 1.10])
            ncalc = length( [0.95 1.0 1.05])
        elseif newst == "2D-mid"
            ncalc = length( [0.86 0.88 0.91 0.96 1.0 1.05 1.10])
        elseif newst == "2D-short"
            ncalc = length( [0.80 0.83])
        elseif newst == "shape"
            ncalc = length( [-0.06 -0.03 0.03 0.06])
        elseif newst == "coords"
#            ncalc = length( [-0.20 -0.17 -0.14 -0.10 -0.07 -0.03 0.0 0.03 0.07 0.10 0.15 0.2 0.25 0.35 0.5])
            ncalc = length([-0.10 -0.07 -0.03 0.0 0.03 0.07 0.10 0.15 0.2 0.25 0.35 0.5])
        elseif newst == "el_party"
            ncalc = length( [-0.20 -0.17 -0.14 -0.10 -0.07 -0.03 0.0 0.03 0.07 0.10 0.15 0.2 0.25 0.35 0.5]) * 2
        elseif newst == "bin_party"
            ncalc = length( [-0.20 -0.17 -0.14 -0.10 -0.07 -0.03 0.0 0.03 0.07 0.10 0.15 0.2 0.25 0.35 0.5]) * 2
        elseif newst == "tern_party"
            ncalc = length( [-0.20 -0.17 -0.14 -0.10 -0.07 -0.03 0.0 0.03 0.07 0.10 0.15 0.2 0.25 0.35 0.5]) * 2
        elseif newst == "coords-scf"
            ncalc = 1
        elseif newst == "coords_min"
            ncalc = 2
        elseif newst == "coords_min_tri"
            if T1 != T2
                ncalc = 6
            else
                ncalc = 1
            end
        elseif newst == "coords-short"
            ncalc = length( [-0.3 -0.27 -0.25])
        elseif newst == "coords-small"
            ncalc = length( [ -0.15  -0.10 -0.05  0.0 0.05 0.10  0.15  ])
        elseif newst == "coords-small2"
            ncalc = length( [ -0.08, 0.0   ])
        elseif newst == "coords_super"
            ncalc = length( 0.08:.03:0.2)
        elseif newst == "break_inv"
            ncalc = length( [0.01 0.02 0.05 0.07 ])
        elseif newst == "shear"
            ncalc = length([0.01 0.02 ])
        elseif newst == "flyaway"
            ncalc = length( [-0.05, 0.0, 0.01, 0.05, 0.2, 0.5, 1.0, 2.0, 3.0, 4.0])
        elseif newst == "scf"
            ncalc = 1
        elseif newst == "scf_small"
            ncalc = 1
        elseif (newst == "coords_trimer"  ||  newst == "coords_trimer2" || newst == "coords_trimer_ab" || newst == "coords_trimer_ab_big" )
            #            ncalc = length([1.05, 1.1, 1.15, 1.2, 1.25, 1.3])
            ncalc = length([1.05, 1.1,  1.2,  1.3])
        elseif newst == "coords_trimer_ab_new"  
            #            ncalc = length([1.05, 1.1, 1.15, 1.2, 1.25, 1.3])
            ncalc = 4*4
        elseif newst == "coords_trimer_a_new"  
            #            ncalc = length([1.05, 1.1, 1.15, 1.2, 1.25, 1.3])
            ncalc = 4*2
        elseif newst == "coords_trimer3"
            ncalc = 5*4*2
        elseif (newst == "coords_trimer_dense" || newst == "coords_trimer_ab_dense" )
            ncalc = length([1.0, 0.95])
        elseif newst == "trimer_tern"
            ncalc = 2
        elseif newst == "trimer_tern_right"
            ncalc = 9
        elseif newst == "trimer_tern_line"
            ncalc = 3
        elseif newst == "trimer_tern_angle"
            ncalc = 3
        else
            println("error newst: ", newst)
            ncalc = 1
        end

        for n in 1:ncalc
            if magnetic
                d="$dir/mag_$name"*"_vnscf_"*"$newst"*"_"*"$ncalc"        
                #                d2="$dir/no_mag_$name"*"_vnscf_"*"$newst"*"_"*"$ncalc"        
                if dft_only
                    if isdir(d*"/qe.save")
                        push!(already_done, d)
                    else
                        push!(not_done, d)
                    end
                else
                    if (isfile(d*"/projham_K.xml") ||  isfile(d*"/projham_K.xml.gz")) # && ( isfile(d2*"/projham_K.xml") ||  isfile(d2*"/projham_K.xml.gz"))
                        push!(already_done, d)
                    else
                        push!(not_done, d)
                    end
                end
            else
                d="$dir/$name"*"_vnscf_"*"$newst"*"_"*"$n"        
                if dft_only
                    if isdir(d*"/qe.save")
                        push!(already_done, d)
                    else
                        push!(not_done, d)
                    end
                else
                    if isfile(d*"/projham_K.xml") ||  isfile(d*"/projham_K.xml.gz")
                        push!(already_done, d)
                    else
                        push!(not_done, d)
                    end
                end
            end
        end
        if check_only
            continue
        end

        if magnetic
            d="$dir/mag_$name"*"_vnscf_"*"$newst"*"_"*"$ncalc"        
            #            d2="$dir/no_mag_$name"*"_vnscf_"*"$newst"*"_"*"$ncalc"        
        else
            d="$dir/$name"*"_vnscf_"*"$newst"*"_"*"$ncalc"        
        end

        
        if isfile(d*"/projham_K.xml") ||  isfile(d*"/projham_K.xml.gz")
            println("everything is already done! $d")
            println("we can move on")
            continue
        else
            println("not done yet")
        end
        
        if dft_only
            if isdir(d*"/qe.save")
                continue
            end
        end

        #######################


        randi = Int64(round(rand()*1000000))
        try
            println("read crys")
            c = ThreeBodyTB.CrystalMod.makecrys(file)

            if newst == "2D_tern" || newst == "3D_tern"
                #arrange in electronegativity order
                en = [electronegativity[T1], electronegativity[T2], electronegativity[T3]]
                ind = sortperm(en)
                T = [deepcopy(T1), deepcopy(T2), deepcopy(T3)]
                T1 = T[ind[1]]
                T2 = T[ind[2]]
                T3 = T[ind[3]]
                println("electronegativity order $T1 $T2 $T3")
            end

            
            for i in 1:c.nat
                if c.types[i] == "A"
                    c.types[i] = T1
                elseif c.types[i] == "B"
                    c.types[i] = T2
                elseif c.types[i] == "C"
                    c.types[i] = T3
                elseif c.types[i] == "O"
                    c.types[i] = "O"
                elseif c.types[i] == "F"
                    c.types[i] = "F"
                elseif c.types[i] == "Fe"
                    c.types[i] = "Fe"
                elseif c.types[i] == "Mn"
                    c.types[i] = "Mn"
                elseif c.types[i] == "H"
                    c.types[i] = "H"
                else
                    c.types[i] = T2
                end
            end

            println(c)
            println("did read crys, run dft")

            if scf != "none"
                

                #preadjust vol
                avg_rad = 0.0
                for t in c.types
                    avg_rad += ThreeBodyTB.Atomdata.atom_radius[t] / 100.0 / 0.529177
                end
                avg_rad = avg_rad / c.nat
                
                vol_peratom = abs(det(c.A)) / c.nat
                
                ratio = (avg_rad^3 * (4 * pi / 3)) / vol_peratom
                println("ratio $ratio")
                if ratio > 1.1
                    c.A = c.A * (ratio^(1.0/3.0))
                    println("new starting c.A")
                    println(c.A)
                end
                

                println("START DFT.runSCF")
                
                if scf != "scf"
                    dft_ref = ThreeBodyTB.DFT.runSCF(c, inputstr=name, nprocs=procs, prefix="$name.qe.relax", directory="$dir", tmpdir="/$tmpname/$name.$randi", wannier=false, code="QE", skip=true, calculation=scf, dofree=free, cleanup=true, magnetic=magnetic)
                    
                    println("did dft, get new struct")
                    
                    if (scf == "relax" || scf == "vc-relax") && maximum(abs.(dft_ref.stress)) > 2e-5
                        println()
                        println("structure convergence not reached, relax again")
                        c2 = dft_ref.crys
                        println(c2)
                        println()
                        
                        dft_ref = ThreeBodyTB.DFT.runSCF(c2, inputstr=name, nprocs=procs, prefix="$name.qe.relax", directory="$dir", tmpdir="/$tmpname/$name.$randi", wannier=false, code="QE", skip=true, calculation=scf, dofree=free, cleanup=true, magnetic=magnetic)
                    end
                    println("END DFT.runSCF")
                    println(dft_ref)

                    cnew = deepcopy(dft_ref.crys)

                else
                    cnew = deepcopy(c)
                end
                    

                #            cnew = deepcopy(c)

            else
                cnew = deepcopy(c)
            end

            torun = []
            if newst == "vol"
                for x in [ 0.95 1.0 1.05]
                    c = deepcopy(cnew)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "vol2"
                for x in [ 0.91 1.0 ]
                    c = deepcopy(cnew)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "vol-mag"
                for x in [ 1.0  ]
                    c = deepcopy(cnew)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "vol-mag-more"
                for x in [0.9 0.95 1.05 1.1  ]
                    c = deepcopy(cnew)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "vol-mid"
                for x in [ 0.9 0.95 1.0 1.05 1.1 ]
                    c = deepcopy(cnew)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "vol-oxygen"
                cx = makecrys(cnew.A, cnew.coords[1:2,:], cnew.types[1:2])

                for x in [ 0.90, 0.96, 1.0  ]
                    c = deepcopy(cx)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end
                cx.coords[1,3] += 0.02
                cx.A = cx.A * 0.93
                push!(torun, deepcopy(cx))


            elseif newst == "2D-oxygen"
                cx = makecrys(cnew.A, cnew.coords[2:3,:], cnew.types[2:3])

                for x in [ 0.91, 0.95, 1.0, 0.85  ]
                    c = deepcopy(cx)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end

            elseif newst == "vol-dense"
                for x in [ 0.87 ]
                    c = deepcopy(cnew)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "vol-verydense"
                for x in [0.77 0.82 0.75 0.70 0.65 0.60 0.55 0.50]
                    c = deepcopy(cnew)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "vol-big"
                for x in [0.80 0.85 0.9 0.95 1.0 1.05 1.1 1.2 1.3 1.5 ]
                    c = deepcopy(cnew)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "vol-huge"
                for x in [0.9 0.95 1.0 1.05 1.1 1.2 1.3 1.5 2.0 2.5 3.0 3.5 4.0 5.0]
                    c = deepcopy(cnew)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "2D_tern"
                for x in [0.90 0.95 1.0 1.05 ]
                    c = deepcopy(cnew)
                    c.A[1:2,:] = c.A[1:2,:] * x
                    push!(torun, deepcopy(c))
                end 
            elseif newst == "3D_tern"
                for x in [ 0.95 1.0 1.05 ]
                    c = deepcopy(cnew)
                    c.A = c.A * x
                    push!(torun, deepcopy(c))
                end 
            elseif newst == "2D"
#                for x in [0.90 0.95 1.0 1.05 1.10]
                for x in [ 0.95 1.0 1.05 ]
                    c = deepcopy(cnew)
                    c.A[1:2,:] = c.A[1:2,:] * x
                    push!(torun, deepcopy(c))
                end 
            elseif newst == "2D-mid"
                for x in [0.86 0.88 0.91 0.96 1.0 1.05 1.10]
                    c = deepcopy(cnew)
                    c.A[1:2,:] = c.A[1:2,:] * x
                    push!(torun, deepcopy(c))
                end 
            elseif newst == "2D-short"
                for x in [0.80 0.83]
                    c = deepcopy(cnew)
                    c.A[1:2,:] = c.A[1:2,:] * x
                    push!(torun, deepcopy(c))
                end 
            elseif newst == "shape"
                for x in [-0.06 -0.03 0.03 0.06]
                    c = deepcopy(cnew)
                    c.A[1:2,:] = c.A[1:2,:] * (1+x)
                    c.A[3,:] = c.A[3,:] * (1- 2.0 * x)
                    push!(torun, deepcopy(c))
                end 
            elseif newst == "coords_trimer_ab_new"

                aa = min_dimer_dist_dict[ (T1,T1)]
                ab = min_dimer_dist_dict[ (T1,T2)]
                bb = min_dimer_dist_dict[ (T2,T2)]

                
                A = diagm([16,8,16.0])
                coords = zeros(3,3)
                t = [T1, T1, T2]
                coords[1,1] = aa / A[1,1] / 2.0
                coords[2,1] = 1.0 - aa / A[1,1] / 2.0
                coords[3,3] = sqrt(ab^2 -  (aa/2 )^2) / A[3,3]
                coords[3,1] = 0.00005

                cnew = makecrys(A, coords, t)
                
                for x in [1.02, 1.07, 1.15 ,  1.3]
                    c = deepcopy(cnew) * x
                    push!(torun, deepcopy(c))
                end

                A = diagm([16,8,16.0])
                coords = zeros(3,3)
                t = [T2, T2, T1]
                coords[1,1] = bb / A[1,1] / 2.0
                coords[2,1] = 1.0 - bb / A[1,1] / 2.0
                coords[3,3] = sqrt(ab^2 -  (bb/2 )^2) / A[3,3]
                coords[3,1] = 0.00005
                
                cnew = makecrys(A, coords, t)
                
                for x in [1.02, 1.07, 1.15 ,  1.3]
                    c = deepcopy(cnew) * x
                    push!(torun, deepcopy(c))
                end
                
                
                A = diagm([16,8,16.0])
                coords = zeros(3,3)
                t = [T1, T1, T2]
                coords[1,1] = aa / A[1,1] / 2.0
                coords[2,1] = 1.0 - aa / A[1,1] / 2.0
                coords[3,3] = ab / A[3,3]
                coords[3,1] = aa / A[1,1] / 2.0

                cnew = makecrys(A, coords, t)
                
                for x in [1.02, 1.07, 1.15, 1.3]
                    c = deepcopy(cnew) * x
                    push!(torun, deepcopy(c))
                end

                A = diagm([16,8,16.0])
                coords = zeros(3,3)
                t = [T2, T2, T1]
                coords[1,1] = bb / A[1,1] / 2.0
                coords[2,1] = 1.0 - bb / A[1,1] / 2.0
                coords[3,3] = ab / A[3,3]
                coords[3,1] = bb / A[1,1] / 2.0

                cnew = makecrys(A, coords, t)
                
                for x in [1.02, 1.07, 1.15, 1.3]
                    c = deepcopy(cnew) * x
                    push!(torun, deepcopy(c))
                end

            elseif newst == "coords_trimer_a_new"

                aa = min_dimer_dist_dict[ (T1,T1)]

                
                A = diagm([15,7,15.0])
                coords = zeros(3,3)
                t = [T1, T1, T2]
                coords[1,1] = aa / A[1,1] / 2.0
                coords[2,1] = 1.0 - aa / A[1,1] / 2.0
                coords[3,3] = sqrt(aa^2 -  (aa/2 )^2) / A[3,3]
                coords[3,1] = 0.00005

                cnew = makecrys(A, coords, t)
                
                for x in [1.02, 1.07, 1.15 ,  1.3]
                    c = deepcopy(cnew) * x
                    push!(torun, deepcopy(c))
                end

                
                A = diagm([15,7,15.0])
                coords = zeros(3,3)
                t = [T1, T1, T2]
                coords[1,1] = aa / A[1,1] / 2.0
                coords[2,1] = 1.0 - aa / A[1,1] / 2.0
                coords[3,3] = aa / A[3,3]
                coords[3,1] = aa / A[1,1] / 2.0

                cnew = makecrys(A, coords, t)
                
                for x in [1.02, 1.07, 1.15, 1.3]
                    c = deepcopy(cnew) * x
                    push!(torun, deepcopy(c))
                end


                
            elseif newst == "coords_trimer_ab"

                a = min_dimer_dist_dict[ cnew.types[1]]
                println("coords_trimer_ab $a ", cnew.types[1], " " , cnew.types[3] )
                ab = 1.0
                try
                    name_t="dimer.in"
                    ncalc_t = 1
                    newst_t = "coords"
                    ab_dir="$dir/$name_t"*"_vnscf_"*"$newst_t"*"_"*"$ncalc_t"        
                    println("try $ab_dir")
                    dft = loadXML(ab_dir*"/qe.save")

                    ab = -dft.crys.coords[1,3] * dft.crys.A[3,3] * 2.0
                    println("ab $ab loaded")

                catch
                    ab = (min_dimer_dist_dict[T1] + min_dimer_dist_dict[T2]) / 2.0 
                    println("ab $ab estimated")

                end
                
                #                for x in [1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
                for x in [1.01, 1.05, 1.1,  1.2,  1.3]

                    c = deepcopy(cnew)
                    if c.coords[1,3] == 0.8
                        c.A[1,1] = ab / (0.4^2 + 0.2^2)^0.5 * x
                    else
                        c.A[1,1] = ab / 0.4 * x
                    end

                    c.A[2,2] = ab * 2.0 * x
                    c.A[3,3] = a / 0.4 * x
                    push!(torun, deepcopy(c))
                end


            elseif newst == "coords_trimer_ab_dense"

                a = min_dimer_dist_dict[ cnew.types[1]]
                println("coords_trimer_ab $a ", cnew.types[1], " " , cnew.types[3] )
                ab = 1.0
                try
                    name_t="dimer.in"
                    ncalc_t = 1
                    newst_t = "coords"
                    ab_dir="$dir/$name_t"*"_vnscf_"*"$newst_t"*"_"*"$ncalc_t"        
                    println("try $ab_dir")
                    dft = ThreeBodyTB.QE.loadXML(ab_dir*"/qe.save")

                    ab = -dft.crys.coords[1,3] * dft.crys.A[3,3] * 2.0
                    println("ab $ab loaded")

                catch
                    ab = (min_dimer_dist_dict[T1] + min_dimer_dist_dict[T2]) / 2.0 
                    println("ab $ab estimated")

                end
                
                #                for x in [1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
                for x in [1.0, 0.95]

                    c = deepcopy(cnew)
                    if c.coords[1,3] == 0.8
                        c.A[1,1] = ab / (0.4^2 + 0.2^2)^0.5 * x
                    else
                        c.A[1,1] = ab / 0.4 * x
                    end

                    c.A[2,2] = ab * 2.0 * x
                    c.A[3,3] = a / 0.4 * x
                    push!(torun, deepcopy(c))
                end


            elseif newst == "coords_trimer_ab_big"

                a = min_dimer_dist_dict[ cnew.types[1]]
                println("coords_trimer_ab big $a ", cnew.types[1], " " , cnew.types[3] )
                ab = 1.0
                try
                    name_t="dimer.in"
                    ncalc_t = 1
                    newst_t = "coords"
                    ab_dir="$dir/$name_t"*"_vnscf_"*"$newst_t"*"_"*"$ncalc_t"        
                    println("try $ab_dir")
                    dft = ThreeBodyTB.QE.loadXML(ab_dir*"/qe.save")

                    ab = -dft.crys.coords[1,3] * dft.crys.A[3,3] * 2.0
                    println("ab $ab loaded")

                catch
                    ab = (min_dimer_dist_dict[T1] + min_dimer_dist_dict[T2]) / 2.0 
                    println("ab $ab estimated")

                end
                
                for x in [1.05, 1.15, 1.25, 1.35, 1.45, 1.55]
                    #                for x in [1.05, 1.1,  1.2,  1.3]

                    c = deepcopy(cnew)
                    c.A[1,1] = ab / (0.4) * x

                    c.A[2,2] = ab * 3.0 * x
                    c.A[3,3] = a / 0.2 * x
                    push!(torun, deepcopy(c))
                end


            elseif newst == "coords_trimer"
                a = min_dimer_dist_dict[(T1,T1)]
                c = makecrys([12 0 0; 0 12 0; 0 0 12]*1.0, [0 0 0; 0 0 a/12; a/12 0 0], [T1, T1, T1])
                for x in [1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
                    push!(torun, deepcopy(c*x))
                end
                c = makecrys([12 0 0; 0 12 0; 0 0 12]*1.0, [0 0 0; a/24 0 a/12; a/12 0 0], [T1, T1, T1])
                for x in [1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
                    push!(torun, deepcopy(c*x))
                end
            elseif newst == "coords_trimerY"
#                a = min_dimer_dist_dict[(T1,T1)]
                a = -1.0*cnew.A[3,3] * (cnew.coords[1,3] - cnew.coords[2,3] - 1)  / 20.0
                for y = [0.83, 0.85, 0.90]
                    b = a *y
                    c = makecrys([20 0 0; 0 14 0; 0 0 20]*1.0, [0 0 0; 0 0 b], [T1, T1])
                    push!(torun, deepcopy(c))
                end
                for y = [0.85, 0.90]
                    b = a *y
                    push!(torun, deepcopy(c))
                    for x in [0.85, 0.90]
                        c = makecrys([20 0 0; 0 14 0; 0 0 20]*1.0, [0 0 0; 0 0 b; a*x 0 0], [T1, T1, T1])
                        push!(torun, deepcopy(c))
                    end

                    for x in [0.85, 0.90]
                        c = makecrys([20 0 0; 0 14 0; 0 0 20]*1.0, [0 0 0; 0 0 b; sqrt(3)/2*a*x 0 a/2*x], [T1, T1, T1])
                        push!(torun, deepcopy(c))
                    end

                    for x in [0.85, 0.90]
                        c = makecrys([20 0 0; 0 14 0; 0 0 20]*1.0, [0 0 0; 0 0 b; 0 0 b+a*x], [T1, T1, T1])
                        push!(torun, deepcopy(c))
                    end
                end
            elseif newst == "coords_trimerX"
#                a = min_dimer_dist_dict[(T1,T1)]
                a = -1.0*cnew.A[3,3] * (cnew.coords[1,3] - cnew.coords[2,3] - 1)  / 20.0
                for y = [0.90, 0.95, 1.0, 1.1, 1.2, 1.4, 1.6, 2.0, 2.5]
                    b = a *y
                    c = makecrys([20 0 0; 0 14 0; 0 0 20]*1.0, [0 0 0; 0 0 b], [T1, T1])
                    push!(torun, deepcopy(c))
                end
                for y = [0.95, 1.0, 1.1]
                    b = a *y
                    push!(torun, deepcopy(c))
                    for x in [0.95, 1.0, 1.1, 1.3]
                        c = makecrys([20 0 0; 0 14 0; 0 0 20]*1.0, [0 0 0; 0 0 b; a*x 0 0], [T1, T1, T1])
                        push!(torun, deepcopy(c))
                    end

                    for x in [0.95, 1.0, 1.1, 1.3]
                        c = makecrys([20 0 0; 0 14 0; 0 0 20]*1.0, [0 0 0; 0 0 b; sqrt(3)/2*a*x 0 a/2*x], [T1, T1, T1])
                        push!(torun, deepcopy(c))
                    end

                    for x in [0.95, 1.0, 1.1, 1.3]
                        c = makecrys([20 0 0; 0 14 0; 0 0 20]*1.0, [0 0 0; 0 0 b; 0 0 b+a*x], [T1, T1, T1])
                        push!(torun, deepcopy(c))
                    end
                end

            elseif newst == "coords_trimer3"
                a = min_dimer_dist_dict[T1]
                for x in [1.15, 1.2, 1.25, 1.3, 1.35]
                    for y in [1.15, 1.2, 1.25, 1.3]
                        for z in [0, 0.5]
                            c = deepcopy(cnew)
                            c.A[1,1] = 15
                            c.A[2,2] = 10
                            c.A[3,3] = 15
                            c.coords[:,:] = zeros(3,3)
                            c.coords[2,1] = a*x/15
                            c.coords[3,3] = a*y/15
                            c.coords[3,1] = a*z/15
                            push!(torun, deepcopy(c))
                        end
                    end
                end

            elseif newst == "coords_trimer_dense"
                a = min_dimer_dist_dict[T1]
                for x in [1.0, 0.95]
                    c = deepcopy(cnew)
                    c.A[1,1] = a / (0.4^2 + 0.2^2)^0.5 * x
                    c.A[2,2] = a * 1.5 * x
                    c.A[3,3] = a / 0.4 * x
                    push!(torun, deepcopy(c))
                end

            elseif newst == "coords_trimer2"
                a = min_dimer_dist_dict[T1]
                for x in [1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
                    c = deepcopy(cnew)
                    c.A[1,1] = a / 0.4 * x
                    c.A[2,2] = a * 1.5 * x
                    c.A[3,3] = a / 0.4 * x
                    push!(torun, deepcopy(c))
                end

            elseif newst == "coords"
#                for x in [-0.20 -0.17 -0.14 -0.10 -0.07 -0.03 0.0 0.03 0.07 0.10 0.15 0.2 0.25 0.35 0.5]
                for x in [-0.10 -0.07 -0.03 0.0 0.03 0.07 0.10 0.15 0.2 0.25 0.35 0.5]
                    c = deepcopy(cnew)
                    c.coords = c.coords * (1+x)
                    push!(torun, deepcopy(c))
                end

            elseif newst == "el_party"
                c = deepcopy(cnew)
                push!(torun, deepcopy(c))
                ThreeBodyTB.ManageDatabase.prepare_database(c)
                COEF = ThreeBodyTB.ManageDatabase.database_cached[(Symbol(c.types[1]), Symbol(c.types[1]))]
                min_dist = COEF.min_dist / 12 /2
                println("min_dist ", min_dist)
                for x in [0, 0.03, 0.05, 0.07, 0.09, 0.12, 0.15, 0.20, 0.25, 0.3, 0.4]
                    
                    c = deepcopy(cnew)
                    c.coords[1,3] = min_dist * (1+x)
                    c.coords[2,3] = -min_dist * (1+x)
                    push!(torun, deepcopy(c))
                end
                c0 = min_dist * 2
                t = cnew.types[1]
                min_dist = c0 * 2 * 0.20
                
                for x in [0.03, 0.05, 0.07]
                    for theta in [  pi/4, pi/2, 3*pi/4, pi]
                        for z2 in [0.01, 0.08, 0.12, 0.15]
#                for x in [ 0.0, 0.25, 0.5, 0.75, 1.0, 1.25]
#                    for z1 in [-1.0, -0.6, -0.3, 0.0, 0.25, 0.5]
#                        for z2 in [-0.15, -0.07, 0.0, 0.07]
                            
                            X = c0  * (1+x) * sin(theta)
                            Z1 = c0 * (1+x) * cos(theta)
                            Z2 = c0 * (1+z2)
 
                            c = makecrys([12.0 0 0 ; 0 7 0; 0 0 12.0], [0 0 0; 0 0 Z2; X 0 Z1] , [t,t,t])
                            D = check_twobody_dist(c)
                            if D[(t,t)] > min_dist
                                push!(torun, deepcopy(c))
                            end

                        end
                    end
                end

            elseif newst == "tern_party"
                for x in [-0.05,  0.0, 0.05]
                    c = deepcopy(cnew)
                    c = c*(1+x)
                    push!(torun, deepcopy(c))
                end
                t1 = cnew.types[1]
                t2 = cnew.types[2]
                t3 = cnew.types[3]

                ab = min_dimer_dist_dict[(t1,t2)]/12.0
                ac = min_dimer_dist_dict[(t1,t3)]/12.0
                bc = min_dimer_dist_dict[(t2,t3)]/12.0

                for x in [0.02, 0.1]
                    for theta in [  pi/4, pi/2, 3*pi/4, pi]
                        for z2 in [ 0.05,  0.2]

                            X = ac  * (1+x) * sin(theta)
                            Z1 = ac * (1+x) * cos(theta)
                            Z2 = ab * (1+z2)
                            c = makecrys([12.0 0 0 ; 0 7 0; 0 0 12.0], [0 0 0; 0 0 Z2; X 0 Z1] , [t1,t2,t3])
                            
                            D = check_twobody_dist(c)
                            good = true
                            for key in keys(D)
                                if D[(t1,t2)] < min_dimer_dist_dict[key]
                                    good = false
                                end
                            end
                            if good
                                push!(torun, deepcopy(c))
                            end
###############


                            X = ab  * (1+x) * sin(theta)
                            Z1 = ab * (1+x) * cos(theta)
                            Z2 = bc * (1+z2)

                            c = makecrys([12.0 0 0 ; 0 7 0; 0 0 12.0], [0 0 0; 0 0 Z2; X 0 Z1] , [t2,t3,t1])
                            
                            D = check_twobody_dist(c)
                            good = true
                            for key in keys(D)
                                if D[(t1,t2)] < min_dimer_dist_dict[key]
                                    good = false
                                end
                            end
                            if good
                                push!(torun, deepcopy(c))
                            end


##
                            X = ab  * (1+x) * sin(theta)
                            Z1 = ab * (1+x) * cos(theta)
                            Z2 = bc * (1+z2)

                            c = makecrys([12.0 0 0 ; 0 7 0; 0 0 12.0], [0 0 0; 0 0 Z2; X 0 Z1] , [t3,t1,t2])
                            
                            D = check_twobody_dist(c)
                            good = true
                            for key in keys(D)
                                if D[(t1,t2)] < min_dimer_dist_dict[key]
                                    good = false
                                end
                            end
                            if good
                                push!(torun, deepcopy(c))
                            end

                            
                        end
                    end
                end


            elseif newst == "bin_party"
                for x in [-0.22 -0.20 -0.17 -0.14 -0.10 -0.07 -0.03 0.0 0.03 0.07 0.10 0.2 0.35 0.5]
                    c = deepcopy(cnew)
                    c.coords = c.coords * (1+x)
                    push!(torun, deepcopy(c))
                end
                c0 = abs(cnew.coords[2,3])
                t1 = cnew.types[1]
                t2 = cnew.types[2]

                

                for types in [[t1,t2,t1], [t2,t1,t1], [t1,t2,t2], [t2,t1,t2] ]

                    ab = min_dimer_dist_dict[(types[1],types[2])]/12.0
                    ac = min_dimer_dist_dict[(types[1],types[3])]/12.0
                    bc = min_dimer_dist_dict[(types[2],types[3])]/12.0
                    
                    for x in [0.01, 0.15]
                        for theta in [  pi/4, pi/2, 3*pi/4, pi]
                            for z2 in [ 0.05,  0.2]

                                X = ac  * (1+x) * sin(theta)
                                Z1 = ac * (1+x) * cos(theta)
                                Z2 = ab * (1+z2)
                                c = makecrys([12.0 0 0 ; 0 7 0; 0 0 12.0], [0 0 0; 0 0 Z2; X 0 Z1] , [types[1], types[2], types[3]])
                                
                                D = check_twobody_dist(c)
                                good = true
                                for key in keys(D)
                                    if D[(t1,t2)] < min_dimer_dist_dict[key]
                                        good = false
                                    end
                                end
                                if good
                                    push!(torun, deepcopy(c))
                                end
                                ###############

                                
                            end
                        end
                    end
                end                
            elseif newst == "coords-short"
                for x in [-0.3 -0.27 -0.25]
                    c = deepcopy(cnew)
                    c.coords = c.coords * (1+x)
                    push!(torun, deepcopy(c))
                end
            elseif newst == "coords-small"
                for x in [ -0.15  -0.10 -0.05  0.0 0.05 0.10  0.15  ]
                    c = deepcopy(cnew)
                    c.coords = c.coords * (1+x)
                    push!(torun, deepcopy(c))
                end
            elseif newst == "coords-small2"
                for x in [ -0.08, 0.0  ]
                    c = deepcopy(cnew)
                    c.coords = c.coords * (1+x)
                    push!(torun, deepcopy(c))
                end
            elseif newst == "coords_super"
                for x in 0.08:.03:0.2
                    c = deepcopy(cnew)
                    c.coords[1,3] = -1.0*x
                    c.coords[2,3] = x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "break_inv"
                for x in [0.01 0.02 0.05 0.07 ]
                    c = deepcopy(cnew)
		    if c.nat == 1
		        c = c * [2 1 1]
                    end
                    c.A = c.A * 1.02
                    c.coords[1,1] = c.coords[1,1] + x
                    push!(torun, deepcopy(c))
                end
            elseif newst == "shear"
                for x in [0.01 0.02 ]
                    c = deepcopy(cnew)
                    c = c 
                    c.A = c.A * (I +  [0 x 0; x 0 0; 0 0 0])

                    push!(torun, deepcopy(c))
                end 

            elseif newst == "flyaway"
                for x in [-0.05, 0.0, 0.01, 0.05, 0.2, 0.5, 1.0, 2.0, 3.0, 4.0]
                    c = deepcopy(cnew)
                    c = c 
                    c.A = c.A * (I +  [0 0 0; 0 0 0; 0 0 x])

                    push!(torun, deepcopy(c))
                end 
                
            elseif newst == "scf"
                push!(torun, deepcopy(cnew))
            elseif newst == "scf_small"
                push!(torun, deepcopy(cnew*0.99))
            elseif newst == "trimer_tern"
                println("trimer_tern")
                counter = 0
                for x in [0.5, 0.7, 0.9, 1.1, 1.3,  1.5,  1.7, 1.9, 2.1]
                    c = deepcopy(cnew)
                    c.coords = c.coords * x
                    if check_twobody_dist(c)
                        counter += 1
                        push!(torun, deepcopy(c))
                        if counter >= 2
                            break
                        end
                    end
                end
            elseif newst == "coords_min"
                c = deepcopy(cnew)
                t = c.types
                a12 = get_twobody_dist(t[1], t[2])
                c.A[:,:] = [10.0 0 0; 0 10.0 0; 0 0 16.0]
                c.coords[1,:] .= 0.0
                c.coords[2,:] .= 0.0
                c.coords[2,3] = a12 / 16.0  * 0.99
                push!(torun, deepcopy(c))
                
                c2 = deepcopy(c)
                c2.coords[2,3] = a12 / 16.0 * 1.035
                push!(torun, deepcopy(c2))

            elseif newst == "coords_min_tri"

                c = deepcopy(cnew)
                t = c.types
                a12 = get_twobody_dist(t[1], t[2])

                A = [16.0 0 0; 0 8.0 0; 0 0 16.0]
                coords1= zeros(3,3)
                coords2= zeros(3,3)
                tt1 = [t[1], t[2], t[1]]
                tt2 = [t[2], t[1], t[2]]

                a11 = get_twobody_dist(t[1], t[1])
                a22 = get_twobody_dist(t[2], t[2])

                coords1[2,3] = a12 / 16.0 * 1.05
                coords2[2,3] = a12 / 16.0 * 1.05

                coords1[3,1] = a11 / 16.0 * 1.05
                coords2[3,1] = a22 / 16.0 * 1.05
                
                c1 = makecrys(A, coords1, tt1)
                c2 = makecrys(A, coords2, tt2)

                push!(torun, deepcopy(c1))

                if t[1] != t[2]
                    push!(torun, deepcopy(c2))
                end

                if t[1] != t[2]
                    A = [8.0 0 0; 0 8.0 0; 0 0 22.0]

                    tt1 = [t[1], t[2], t[2]]
                    tt2 = [t[2], t[1], t[1]]

                    coords1= zeros(3,3)
                    coords2= zeros(3,3)

                    coords1[2,3] = a12 / 22.0 * 1.05
                    coords1[3,3] = -a12 / 22.0 * 1.05
                    
                    coords2[2,3] = a12 / 22.0 * 1.05
                    coords2[3,3] = -a12 / 22.0 * 1.05
                    
                    c1 = makecrys(A, coords1, tt1)
                    c2 = makecrys(A, coords2, tt2)

                    push!(torun, deepcopy(c1))
                    push!(torun, deepcopy(c2))

                    A = [16.0 0 0; 0 8.0 0; 0 0 16.0]

                    tt1 = [t[1], t[2], t[2]]
                    tt2 = [t[2], t[1], t[1]]

                    coords1= zeros(3,3)
                    coords2= zeros(3,3)

                    coords1[2,1] = a12 / 16.0 * 1.05
                    coords1[3,3] = -a12 / 16.0 * 1.05
                    
                    coords2[2,1] = a12 / 16.0 * 1.05
                    coords2[3,3] = -a12 / 16.0 * 1.05
                    
                    c1 = makecrys(A, coords1, tt1)
                    c2 = makecrys(A, coords2, tt2)

                    push!(torun, deepcopy(c1))
                    push!(torun, deepcopy(c2))


                end
                
                
            elseif newst == "trimer_tern_right"
                println("trimer_tern right")
                counter = 0
                c = deepcopy(cnew)
                
                orders = [[1,2,3], [2,1,3], [3,1,2]]
                for o in orders
                    t = c.types[o]
                    a12 = get_twobody_dist(t[1], t[2])
                    a13 = get_twobody_dist(t[1], t[3])
                    coords = zeros(3,3)
                    coords[2,1] = a12/12.0 * 1.02
                    coords[3,2] = a13/12.0 * 1.02
                    c2 = makecrys(c.A, coords, t)
                    counterX = 0
                    for x in 1.0:0.02:5.0
                        c3 = deepcopy(c2)
                        c3.A[1,:] = c3.A[1,:] * x
                        c3.A[2,:] = c3.A[2,:] * x
                        if check_frontier(c3)
                            push!(torun, deepcopy(c3))
                            counterX += 1
                            if counterX >= 3
                                break
                            end
                        end
                    end
                end

            elseif newst == "trimer_tern_line"
                println("trimer_tern line")
                counter = 0
                c = deepcopy(cnew)
                
                orders = [[1,2,3], [2,1,3], [3,1,2]]
                for o in orders
                    t = c.types[o]
                    a12 = get_twobody_dist(t[1], t[2])
                    a13 = get_twobody_dist(t[1], t[3])
                    coords = zeros(3,3)
                    coords[2,3] = a12/16.0 * 1.04
                    coords[3,3] = -a13/16.0 * 1.04
                    c2 = makecrys(c.A, coords, t)
                    counterX = 0
                    for x in [1.0]
                        c3 = deepcopy(c2)
                        c3.A[1,:] = c3.A[1,:] * x
                        c3.A[2,:] = c3.A[2,:] * x
                        if check_frontier(c3)
                            push!(torun, deepcopy(c3))
                            counterX += 1
                            if counterX >= 3
                                break
                            end
                        end
                    end
                end

            elseif newst == "trimer_tern_angle"
                println("trimer_tern angle")
                counter = 0
                c = deepcopy(cnew)
                
                orders = [[1,2,3], [2,1,3], [3,1,2]]
                for o in orders
                    t = c.types[o]
                    a12 = get_twobody_dist(t[1], t[2])
                    a13 = get_twobody_dist(t[1], t[3])
                    coords = zeros(3,3)
                    coords[2,1] = a12/12.0 * 1.02
                    coords[3,2] = a13/12.0 * 1.02 * cos(30/180)
                    coords[3,1] = a13/12.0 * 1.02 * sin(30/180)
                    c2 = makecrys(c.A, coords, t)
                    counterX = 0
                    for x in 1.0:0.02:5.0
                        c3 = deepcopy(c2)
                        c3.A[1,:] = c3.A[1,:] * x
                        c3.A[2,:] = c3.A[2,:] * x
                        if check_frontier(c3)
                            push!(torun, deepcopy(c3))
                            counterX += 1
                            if counterX >= 1
                                break
                            end
                        end
                    end
                end
                
            else
                println("error newst: ", newst)
            end
            

            for (i,c) in enumerate(torun)
                println("start torun $name  $i")
                try

                    if magnetic
                        d="$dir/mag_$name"*"_vnscf_"*"$newst"*"_"*"$i"        
                        #                        d2="$dir/no_mag_$name"*"_vnscf_"*"$newst"*"_"*"$i"        
                    else
                        d="$dir/$name"*"_vnscf_"*"$newst"*"_"*"$i"        
                    end
                    println(d)
                    println(c)
                    if isfile(d*"/projham_K.xml") ||  isfile(d*"/projham_K.xml.gz")
                        println("continue")
                        continue
                    end

                    dft = ThreeBodyTB.DFT.runSCF(c, nprocs=procs, prefix="qe", directory="$d", tmpdir="$d", wannier=false, code="QE", skip=true, cleanup=true, magnetic=magnetic)
                    if calc_mode == "nscf"

                        try
                            if !dft_only 
                                tbc, tbck = ThreeBodyTB.AtomicProj.projwfc_workf(dft, nprocs=procs, directory=d, skip_og=true, skip_proj=true, freeze=true, localized_factor = 0.15, cleanup=true, only_kspace=only_kspace, min_nscf = min_nscf)
                            end
                        catch err3
                            println("err3")
                            println(err3)
                            println("skip nscf $i ")##
                        end
                    end

                    if magnetic && false
                        dft = ThreeBodyTB.DFT.runSCF(c, nprocs=procs, prefix="qe", directory="$d2", tmpdir="$d2", wannier=false, code="QE", skip=true, cleanup=true, magnetic=false)
                        if calc_mode == "nscf"
                            try
                                if !dft_only 
                                    tbc, tbck = ThreeBodyTB.AtomicProj.projwfc_workf(dft, nprocs=procs, directory=d2, skip_og=true, skip_proj=true, freeze=true, localized_factor = 0.15, cleanup=true, only_kspace=only_kspace, min_nscf = min_nscf)
                                end
                                catch err3
                                println("err3")
                                println(err3)
                                println("skip nscf $i ")##
                            end
                            
                        end
                    end

                catch err2
                    println("err2")
                    println(err2)
                    println("skip $i dft $i ")##
                    
                end
            end
        catch err
            println(err)
            println("skip everything  $st ")
        end
    end
    
    if check_only
        return already_done, not_done
    end

end


pd = setup_proto_data()

function structure_substitute(atom1, atom2, atom3)

    transmetals = ["Sc", "Y", "La", "Ti", "Zr", "Hf", "V", "Nb", "Ta", "Cr", "Mo", "W", "Mn", "Tc", "Re", "Fe", "Ru", "Os", "Co", "Rh", "Ir", "Ni", "Pd", "Pt", "Cu", "Ag", "Au", "Zn", "Cd", "Hg", "B", "Al", "Ga", "In", "Tl"]

    STRUCTDIR = ThreeBodyTB.STRUCTDIR
    summ = readdlm("$STRUCTDIR/ternary/summ.csv")

    satom = sort([atom1, atom2, atom3])

    function arediff(st)
        return (st[1] != st[2]) || (st[1] != st[3]) || (st[2] != st[3])
    end

    orders = [[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]]


    sl = []
    for i in 1:3
        push!(sl, sub_list[satom[i]])
    end
    slalt = []
    for i in 1:3
        if satom[i] in transmetals
            push!(slalt,  transmetals)
        else
            push!(slalt,  sub_list[satom[i]])
        end
    end

    CR0 = []
    SUB0 = []

    #find exact matches
    for i in 1:size(summ)[1]
        at = summ[i,5:7]
        if satom == sort(at)
            c = makecrys(String(summ[i,1]))
            push!(CR0, c)
            push!(SUB0, [String(summ[i,2]), c.types[1], c.types[2], c.types[3]])

        end
    end

    #generate all single subs
    SUB1 = []
    SUB2 = []
    SUBTRANS = []
    println("generate")
    
    CR1 = []
    CR2 = []
    CRTRANS = []

    for i in 1:3
        if i == 1
            others = [2,3]
        elseif i == 2
            others = [1,3]
        elseif i == 3
            others = [1,2]
        end

        for s in sl[i]
            temp = [s, satom[others[1]], satom[others[2]]]
            orig = [satom[i], satom[others[1]], satom[others[2]]]
            if arediff(temp)
                for o in orders
                    to = temp[o]
                    for i in 1:size(summ)[1]
#                        nat = parse(Int64, summ[i,4])
                        nat = summ[i,4]
                        if nat > 6
                            continue
                        end
                        at = summ[i,5:7]
                        if to == at
                            #                            println(at, " ", orig[o])
                            c = makecrys(String(summ[i,1]))
                            for (ind,t) in enumerate(c.types)
                                for ii = 1:3
                                    if t == at[ii]
                                        c.types[ind] = orig[o[ii]]
                                    end
                                end
                            end
                            push!(CR1, c)

                            push!(SUB1, [String(summ[i,2]), orig[o[1]], orig[o[2]], orig[o[3]]])
                        end
                    end
                end
            end
        end

        for s in slalt[i]
            temp = [s, satom[others[1]], satom[others[2]]]
            orig = [satom[i], satom[others[1]], satom[others[2]]]
            if arediff(temp)
                for o in orders
                    to = temp[o]
                    for i in 1:size(summ)[1]
                        at = summ[i,5:7]
                        nat = summ[i,4]
#                        nat = parse(summ[i,4], Int64)
                        if nat > 6
                            continue
                        end

                        if to == at
                            c = makecrys(String(summ[i,1]))
                            for (ind,t) in enumerate(c.types)
                                for ii = 1:3
                                    if t == at[ii]
                                        c.types[ind] = orig[o[ii]]
                                    end
                                end
                            end

                            push!(CRTRANS, c)
                            #                            println(at, " t ", orig[o])
                            push!(SUBTRANS, [String(summ[i,2]), orig[o[1]], orig[o[2]], orig[o[3]]])
                        end
                    end
                end
            end
        end



    end

    #generate all double subs

    println("generate")
    for i in 1:3
        if i == 1
            others = [2,3]
        elseif i == 2
            others = [1,3]
        elseif i == 3
            others = [1,2]
        end
        for s1 in sl[others[1]]
            for s2 in sl[others[2]]
                temp = [s1,s2,satom[i]]
                #            temp = [s, satom[others[1]], satom[others[2]]]
                orig = [satom[others[1]], satom[others[2]], satom[i]]
                if arediff(temp)
                    for o in orders
                        to = temp[o]
                        for i in 1:size(summ)[1]
                            at = summ[i,5:7]
                            #nat = parse(summ[i,4], Int64)
                            nat = summ[i,4]
                            if nat > 6
                                continue
                            end
                            if to == at
                                c = makecrys(String(summ[i,1]))
                                for (ind,t) in enumerate(c.types)
                                    for ii = 1:3
                                        if t == at[ii]
                                            c.types[ind] = orig[o[ii]]
                                        end
                                    end
                                end

                                push!(CR2, c)

                                #                                println(at, " d ", orig[o])
                                push!(SUB2, [String(summ[i,2]), orig[o[1]], orig[o[2]], orig[o[3]]])
                            end
                        end
                    end
                end
            end
        end
    end
    return CR0, CR1, CR2, CRTRANS, SUB0, SUB1, SUB2, SUBTRANS


end

#    println("EXACT")
#    for i in exact
#        println(summ[i,:])
#    end
#    println()

#        if sat == satom
#            exact = [exact; i]
#        else
#            at1 == atom1

function do_run_ternary_sub(at1, at2, at3, dir,procs, n1=6, n2 = 12; min_nscf = false, makecrys_only=false)

    CR0, CR1, CR2, CRTRANS, SUB0, SUB1, SUB2, SUBTRANS = structure_substitute(at1, at2, at3)

    c0 = length(CR0)
    c1 = length(CR1)
    c2 = length(CR2)
    ct = length(CRTRANS)

    CBIG = [CR1 ;  CR2 ; CRTRANS]
    SBIG = [SUB1 ; SUB2 ; SUBTRANS]
    
    satom = sort([at1,at2, at3])
    
    DONE = []
    DONE_TYPES = []

    DO_RELAX = []


    torun = []
    n = 0

    #add exact atom matches at 95 percent, 105 %,  assume that 100 percent is already done somewhere
    for (i,c) in enumerate(CR0)
        name = SUB0[i][1]
        if name in DONE
            continue
        end
        n += 1

        if n <= n1
            push!(torun, deepcopy(c) * 0.94)
            push!(DONE, name)
            push!(DO_RELAX, false)

            push!(torun, deepcopy(c) * 0.99)
            push!(DONE, name)
            push!(DO_RELAX, false)

            push!(torun, deepcopy(c) * 1.04)
            push!(DONE, name)
            push!(DO_RELAX, false)

            push!(torun, deepcopy( ThreeBodyTB.CrystalMod.generate_random_distortion(c, 0.2, 0.001)))
            push!(DONE, name)
            push!(DO_RELAX, false)

            push!(torun, deepcopy(c) * 0.90)
            push!(DONE, name)
            push!(DO_RELAX, false)

        end
        n += 1
        #if n >= n2
        #    break
        #end
        push!(DONE_TYPES, sort(c.types))
    end

    #add others with differing stoichiometry
    for (i,c) in enumerate(CBIG)
        name = SBIG[i][1]
        if name in DONE
            continue
        end
        if sort(c.types) in DONE_TYPES
            continue
        end

        n += 1
        push!(DO_RELAX, true)
        push!(torun, deepcopy(c) )
        push!(DONE, name)
        if n >= n2
            break
        end
        push!(DONE_TYPES, sort(c.types))
    end

    #add anything else if we still need more structures to reach n2
    if n < n2
        for (i,c) in enumerate(CBIG)
            name = SBIG[i][1]
            if name in DONE #&& n >= n1
                continue
            end
            n += 1
            push!(torun, deepcopy(c) )
            push!(DONE, name)
            push!(DO_RELAX, true)
            if n >= n2
                break
            end

        end            
    end
    println("DONE")
    for d in DONE
        println(d)
    end

    if makecrys_only
        return torun
    end

    for (i,c) in enumerate(torun)
        name = DONE[i]

        do_relax = DO_RELAX[i]

        println("running dft ternary")
        println(c)
        println()
        d="$dir/$name"*"_vnscf_"*"$i"        
        try
            if !isdir("$d/qe.save")
                if do_relax
                    dftR = ThreeBodyTB.DFT.runSCF(c, calculation = "vc-relax", nstep=3, nprocs=procs, prefix="qe", directory="$d", tmpdir="$d", wannier=false, code="QE", skip=true, cleanup=true)
                    cr = deepcopy(dftR.crys)
                else
                    cr = c
                end
                dft = ThreeBodyTB.DFT.runSCF(cr, nprocs=procs, prefix="qe", directory="$d", tmpdir="$d", wannier=false, code="QE", skip=false, cleanup=true)
            else
                dft = ThreeBodyTB.DFT.runSCF(cr, nprocs=procs, prefix="qe", directory="$d", tmpdir="$d", wannier=false, code="QE", skip=true, cleanup=true)
            end
            tbc, tbck = ThreeBodyTB.AtomicProj.projwfc_workf(dft, nprocs=procs, directory=d, skip_og=true, skip_proj=true, freeze=true, localized_factor = 0.15, cleanup=true, only_kspace=true, min_nscf = min_nscf)
        catch err
            println("err dft $d")
            println(err)
        end
        println("done run")
        flush(stdout)

    end

    return torun

end



function oxidation_guess(atom1, atom2, atom3)

    possible_configs = zeros(Int64, 0, 7)

    for (c1, o1) = enumerate(atom_prefered_oxidation[atom1])
        for (c2,o2) = enumerate(atom_prefered_oxidation[atom2])
            for (c3,o3) = enumerate(atom_prefered_oxidation[atom3])


                if c1 == 1
                    score1=1
                elseif o1 == 0
                    score1=4
                else
                    score1=2
                end

                if c2 == 1
                    score2=1
                elseif o2 == 0
                    score2=4
                else
                    score2=2
                end

                if c3 == 1
                    score3=1
                elseif o3 == 0
                    score3=4
                else
                    score3=2
                end


                for n1 = 1:6
                    for n2 = 1:6
                        for n3 = 1:6
                            
                            if (n1+n2+n3) >= 6
                                continue
                            end

                            if n1 == 1 && n2 == 1 && n3 == 1
                                continue
                            end

                            if gcd(n1,n2) > 1 && gcd(n1,n3) > 1 && gcd(n2,n3) > 1
                                continue
                            end
                            if o1 * n1 + o2 * n2 + o3 * n3 == 0
                                possible_configs = [possible_configs; n1 n2 n3 (score1+score2+score3)*2+n1+n2+n3   score1  score2 score3]
                            elseif  abs(o1 * n1 + o2 * n2 + o3 * n3) == 1
                                possible_configs = [possible_configs; n1 n2 n3 (score1+score2+score3)*2+n1+n2+n3+5  score1  score2 score3] 
                            end
                        end
                    end
                end
            end
        end
    end
    possible_configs = possible_configs[sortperm(possible_configs[:, 4]), :]
    

end



function oxidation_guess(atom1, atom2)

    if atom1 == atom2
        keep = [[atom1, atom1, :core_mono]]
        keep = push!(keep , [atom1, atom1, :A0])

        bigmetals = ["Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", "Sc", "Y", "La", "Hg", "H"]

        if atom1 in bigmetals
            push!(keep, [atom1, atom1, "sc_verydense"])
            push!(keep, [atom1, atom1, "fcc_verydense"])
            push!(keep, [atom1, atom1, "bcc_verydense"])
            push!(keep, [atom1, atom1, "diamond_verydense"])
        end


        return keep
    end


    possible_configs = zeros(Int64, 0, 5)

    for (c1, o1) = enumerate(atom_prefered_oxidation[atom1])
        for (c2,o2) = enumerate(atom_prefered_oxidation[atom2])

            if c1 == 1
                score1=1
            elseif o1 == 0
                score1=4
            else
                score1=2
            end

            if c2 == 1
                score2=1
            elseif o2 == 0
                score2=4
            else
                score2=2
            end


            for n1 = 1:6
                for n2 = 1:6
                    
                    if gcd(n1,n2) > 1
                        continue
                    end

                    if o1 * n1 + o2 * n2 == 0
                        possible_configs = [possible_configs; n1 n2 score1+score2   score1  score2]
                    elseif  abs(o1 * n1 + o2 * n2) == 1
                        possible_configs = [possible_configs; n1 n2 score1+score2+10  score1  score2] 
                    end
                end
            end
        end
    end

    possible_configs = possible_configs[sortperm(possible_configs[:, 3]), :]

    #    println("possible configs $atom1 $atom2 score score1 score2")
    #    for p in 1:min(size(possible_configs,1), 6)
    #        println(possible_configs[p,:])
    #    end

    n_config = 0
#    keep = []
    
    keep = [[atom1, atom2, :core_binary]]

    for p in 1:size(possible_configs,1)
        use = false
        if possible_configs[p,3] == 2 || possible_configs[p,3] == 3
            #            println("a ",  possible_configs[p,:])
            use = true
            n_config += 1
        elseif possible_configs[p,3] == 4 && n_config <= 2
            use = true
            #            println("b ",  possible_configs[p,:])
            n_config += 1
            #        else
            #            println("R ",  possible_configs[p,:])

        end
        if use
            if possible_configs[p,1] == 1 && possible_configs[p,2] == 1
                push!(keep, [atom1, atom2, :A1B1])
            elseif possible_configs[p,1] == 1 && possible_configs[p,2] == 2
                push!(keep, [atom1, atom2, :A1B2])
            elseif possible_configs[p,1] == 1 && possible_configs[p,2] ==3
                push!(keep, [atom1, atom2, :A1B3])
            elseif possible_configs[p,1] == 1 && possible_configs[p,2] == 4
                push!(keep, [atom1, atom2, :A1B4])
            elseif possible_configs[p,1] == 1 && possible_configs[p,2] == 5
                push!(keep, [atom1, atom2, :A1B5])
            elseif possible_configs[p,1] == 1 && possible_configs[p,2] == 6
                push!(keep, [atom1, atom2, :A1B6])
            elseif possible_configs[p,1] == 2 && possible_configs[p,2] == 3
                push!(keep, [atom1, atom2, :A2B3])
            elseif possible_configs[p,1] == 2 && possible_configs[p,2] == 5
                push!(keep, [atom1, atom2, :A2B5])
            elseif possible_configs[p,1] == 3 && possible_configs[p,2] == 5
                push!(keep, [atom1, atom2, :A3B5])


            elseif possible_configs[p,2] == 1 && possible_configs[p,1] == 2
                push!(keep, [atom2, atom1, :A1B2])
            elseif possible_configs[p,2] == 1 && possible_configs[p,1] == 3
                push!(keep, [atom2, atom1, :A1B3])
            elseif possible_configs[p,2] == 1 && possible_configs[p,1] == 4
                push!(keep, [atom2, atom1, :A1B4])
            elseif possible_configs[p,2] == 1 && possible_configs[p,1] == 5
                push!(keep, [atom2, atom1, :A1B5])
            elseif possible_configs[p,2] == 1 && possible_configs[p,1] == 6
                push!(keep, [atom2, atom1, :A1B6])
            elseif possible_configs[p,2] == 2 && possible_configs[p,1] == 3
                push!(keep, [atom2, atom1, :A2B3])
            elseif possible_configs[p,2] == 2 && possible_configs[p,1] == 5
                push!(keep, [atom2, atom1, :A2B5])
            elseif possible_configs[p,2] == 3 && possible_configs[p,1] == 5
                push!(keep, [atom2, atom1, :A3B5])
            end


        end
    end
    if length(keep) <= 1
#        if maximum(atom_prefered_oxidation[atom1]) > 0 || maximum(atom_prefered_oxidation[atom2]) > 0
        push!(keep, [atom1, atom2, :metals])
#        end
    end


    #special configurations


    if atom1 in ["B"] #&&  atom2 in ["H", "Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", "Y", "La", "Sc", "Tl", "In", "Ga", "Si", "Ge", "C"]
        push!(keep, [atom2, atom1, "cab6"])
    end
    if atom2 in ["B"] #&&  atom1 in ["H", "Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", "Y", "La", "Sc", "Tl", "In", "Ga", "Si", "Ge", "C"]
        push!(keep, [atom1, atom2, "cab6"])
    end 

    if atom1 in ["Cr", "Mo", "W", "Mn", "Tc", "Re", "Fe", "Ru", "Os"] &&  atom2 in ["N"]
        push!(keep, [atom1, atom2, "re3n"])
        push!(keep, [atom1, atom2, "re3n3"])
    end 

    if atom2 in ["Cr", "Mo", "W", "Mn", "Tc", "Re", "Fe", "Ru", "Os"] &&  atom1 in ["N"]
        push!(keep, [atom2, atom1, "re3n"])
        push!(keep, [atom2, atom1, "re3n3"])
    end 

    if atom1 in ["O", "H", "F", "S", "Cl", "C", "N", "P"] &&  atom2 in ["O", "H", "F", "S", "Cl", "C", "N", "P"] && atom1 != atom2
        push!(keep, [atom1, atom2, "mol_h2o_v1"])
        push!(keep, [atom2, atom1, "mol_h2o_v2"])
    end 

    if atom1 in ["La", "Y", "Sc"] &&  atom2 in ["Ni", "Pd", "Pt"] 
        push!(keep, [atom2, atom1, "ni4la2"])
    end 
    if atom2 in ["La", "Y", "Sc"] &&  atom1 in ["Ni", "Pd", "Pt"] 
        push!(keep, [atom1, atom2, "ni4la2"])
    end 

    if atom1 in ["Rb", "K", "Cs", "Na"] &&  atom2 in ["Al", "Ga", "In", "Tl"] 
        push!(keep, [atom1, atom2, "rb1in4"])
    end 
    if atom2 in ["Rb", "K", "Cs", "Na"] &&  atom1 in ["Al", "Ga", "In", "Tl"] 
        push!(keep, [atom2, atom1, "rb1in4"])
    end 

    if atom1 in ["Nb", "V", "Ta"] &&  atom2 in ["F"]
        push!(keep, [atom1, atom2, "nb1f4"])
    end 
    if atom2 in ["Nb", "V", "Ta"] &&  atom1 in ["F"]
        push!(keep, [atom2, atom1, "nb1f4"])
    end 

    if atom1 in ["Li", "Na", "Cu", "Ag", "Au", "Hg"] &&  atom2 in ["O", "S", "Se"]
        push!(keep, [atom1, atom2, "cu2o1"])
    end 
    if atom2 in ["Li", "Na", "Cu", "Ag", "Au", "Hg"] &&  atom1 in ["O", "S", "Se"]
        push!(keep, [atom2, atom1, "cu2o1"])
    end 

    if atom1 in ["Al", "Ga", "In"] &&  atom2 in ["Ca", "Sr", "Ba"]
        push!(keep, [atom1, atom2, "ga4sr"])
    end 
    if atom2 in ["Al", "Ga", "In"] &&  atom1 in ["Ca", "Sr", "Ba"]
        push!(keep, [atom2, atom1, "ga4sr"])
    end 

    if atom1 in ["H", "Ca", "Sr", "Ba"] &&  atom2 in ["Pd", "Ni", "Pt"]
        push!(keep, [atom1, atom2, "ca1pd5"])
    end 
    if atom2 in ["H", "Ca", "Sr", "Ba"] &&  atom1 in ["Pd", "Ni", "Pt"]
        push!(keep, [atom2, atom1, "ca1pd5"])
    end 

    if atom1 in ["B", "Al", "Ga", "In", "Tl", "Sc", "Y", "La", "Sb", "Bi", "Co", "Fe", "Ni", "Mn", "Cr", "Ti", "V", "Co"] &&  atom2 in ["F"]
        push!(keep, [atom1, atom2, "a2f6"])
    end 
    if atom2 in ["B", "Al", "Ga", "In", "Tl", "Sc", "Y", "La", "Sb", "Bi", "Co", "Fe", "Ni", "Mn", "Cr", "Ti", "V", "Co" ] &&  atom1 in ["F"]
        push!(keep, [atom2, atom1, "a2f6"])
    end 

    if atom1 in ["Fe", "Mn", "Co", "V", "Cr", "Mo", "W", "Mn", "Tc", "Re", "Ru", "Os", "Se", "Te"] &&  atom2 in ["F"]
        push!(keep, [atom1, atom2, "a1f6"])
        push!(keep, [atom1, atom2, "mof6"])
    end 
    if atom2 in ["Fe", "Mn", "Co", "V", "Cr", "Mo", "W", "Mn", "Tc", "Re", "Ru", "Os", "Se", "Te"] &&  atom1 in ["F"]
        push!(keep, [atom2, atom1, "a1f6"])
        push!(keep, [atom2, atom1, "mof6"])
    end 


    if atom1 in ["Cr", "Mo", "W","Tc", "Re","Mn", "Se", "Te"] &&  atom2 in ["O"]
        push!(keep, [atom1, atom2, "te2o6"])
    end 
    if atom2 in ["Cr", "Mo", "W","Tc", "Re","Mn","Se", "Te"] &&  atom1 in ["O"]
        push!(keep, [atom2, atom1, "te2o6"])
    end 

    if atom1 in ["Ge", "Sn", "Pb"] &&  atom2 in ["O"]
        push!(keep, [atom1, atom2, "sn2o2"])
    end 
    if atom2 in ["Ge", "Sn", "Pb"] &&  atom1 in ["O"]
        push!(keep, [atom2, atom1, "sn2o2"])
    end 

    if atom1 in ["Li", "Na", "K", "Rb", "Cs"] &&  atom2 in ["N"]
        push!(keep, [atom1, atom2, "k2n6"])
    end 
    if atom2 in ["Li", "Na", "K", "Rb", "Cs"] &&  atom1 in ["N"]
        push!(keep, [atom2, atom1, "k2n6"])
    end 
    
    transmetals = ["Sc", "Y", "La", "Ti", "Zr", "Hf", "V", "Nb", "Ta", "Cr", "Mo", "W", "Mn", "Tc", "Re", "Fe", "Ru", "Os", "Co", "Rh", "Ir", "Ni", "Pd", "Pt", "Cu", "Ag", "Au", "Zn", "Cd", "Hg", "B", "Al", "Ga", "In", "Tl"]
    anions  = [ "C", "Si", "Ge", "Sn", "Pb", "N", "P", "As", "Sb", "Bi", "O", "S", "Se", "Te", "F", "Cl", "Br", "I", "H"]
    stronganions  = [  "N", "P", "As", "O", "S", "Se", "F", "Cl", "Br", "I"]
    othermetals = ["Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba"]
    metals = [transmetals; othermetals]

    if atom1 in transmetals && atom2 == "H" && atom1 != atom2
        push!(keep, [atom2, atom1, "lipd3"])
    end
    if atom2 in transmetals && atom1 == "H"  && atom1 != atom2
        push!(keep, [atom1, atom2, "lipd3"])
    end


    if atom1 in ["Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", "Y", "La", "Sc"] &&  atom2 in ["C"]
        push!(keep, [atom2, atom1, "c2ca"])
    end 
    if atom2 in ["Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", "Y", "La", "Sc"] &&  atom1 in ["C"]
        push!(keep, [atom1, atom2, "c2ca"])
    end 

    if atom1 in transmetals &&  atom2 in ["N"]
        push!(keep, [atom1, atom2, "osn2"])
    end 
    if atom2 in transmetals &&  atom1 in ["N"]
        push!(keep, [atom2, atom1, "osn2"])
    end 

    if atom1 in ["Pt", "Pd", "Ni", "Co", "Rh", "Ir", "Cu", "Ag", "Au"] &&  atom2 in ["O", "S"]
        push!(keep, [atom1, atom2, "pt2o2"])
    end 
    if atom2 in ["Pt", "Pd", "Ni", "Co", "Rh", "Ir", "Cu", "Ag", "Au"] &&  atom1 in ["O", "S"]
        push!(keep, [atom2, atom1, "pt2o2"])
    end 


    if atom1 in metals &&  atom2 in anions
        push!(keep, [atom1, atom2, "pd3te"])
    end 
    if atom2 in metals &&  atom1 in anions
        push!(keep, [atom2, atom1, "pd3te"])
    end 

    if atom1 in metals &&  atom2 in ["B", "C", "N"]
        push!(keep, [atom2, atom1, "bpt"])
    end 
    if atom2 in metals &&  atom1 in ["B", "C", "N"]
        push!(keep, [atom1, atom2, "bpt"])
    end 

    if atom1 in ["Ru", "Os"] &&  atom2 in ["O"]
        push!(keep, [atom1, atom2, "sif4"])
        push!(keep, [atom1, atom2, "gei4_mol"])
    end 
    if atom2 in ["Ru", "Os"] &&  atom1 in ["O"]
        push!(keep, [atom2, atom1, "sif4"])
        push!(keep, [atom2, atom1, "gei4_mol"])
    end 
#
#    if atom1 in anions && !(atom2 in anions)
#        push!(keep, [atom2, atom1, "bcc_5lay"])
#        push!(keep, [atom2, atom1, "fcc_5lay"])
#    end
#    if atom2 in anions && !(atom1 in anions)
#        push!(keep, [atom1, atom2, "bcc_5lay"])
#        push!(keep, [atom1, atom2, "fcc_5lay"])
#    end

    if atom1 in ["Re", "Tc", "Ru", "Os", "Rh", "Ir", "Pd", "Pt", "Ag", "Au", "Ni"] && atom2 in ["Te", "Se", "S", "Be", "B"]
        push!(keep, [atom1, atom2, "pd3s"])
    end
    if atom2 in ["Re", "Tc", "Ru", "Os", "Rh", "Ir", "Pd", "Pt", "Ag", "Au", "Ni"] && atom1 in ["Te", "Se", "S", "Be", "B"]
        push!(keep, [atom2, atom1, "pd3s"])
    end

    bigmetals = ["Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", "Sc", "Y", "La", "Hg", "H"]

#    if (atom1 in bigmetals && !(atom2 in stronganions)) || (atom2 in bigmetals && !(atom1 in stronganions)) 
#        push!(keep, [atom1, atom2, "hh_oxygen"])
#        push!(keep, [atom1, atom2, "hh_oxygen2"])
#        push!(keep, [atom1, atom2, "hex_oxygen"])
#    end

    if atom1 == atom2 && atom1 in bigmetals
        push!(keep, [atom1, atom1, "sc_verydense"])
        push!(keep, [atom1, atom1, "fcc_verydense"])
        push!(keep, [atom1, atom1, "bcc_verydense"])
        push!(keep, [atom1, atom1, "diamond_verydense"])
    end

    if atom1 in metals && atom2 in ["H", "N", "O"]
        push!(keep, [atom2, atom1, "h2_atom"])
    end
    if atom2 in metals && atom1 in ["H", "N", "O"]
        push!(keep, [atom1, atom2, "h2_atom"])
    end
    

    #    for k in keep
    #        println(k)
    #    end
    return keep 
end

#    A1B1 = ["hbn_real", "sis_2d", "ges", "sns", "wurtz"]
#    A1B2 = ["mgcl2", "mgcl2_tet", "caf2", "sis2", "tio2_rutile", "co2", "square_ab2", "mgf2",  "ticl2", "gei2"]   # "sns2" duplictes ticl2 #"mgf2_v2" is rutile
#    A1B3 = ["alf3", "asna3_2d", "ab3_mol"]
#    A1B4 = ["sif4", "snf4", "gei4_mol"]
#    A1B5 = ["ascl5"]
#    A2B3 = ["y2o3", "p2ca3", "al2o3", "bi2se3", "ga2s3", "gas"]
#    A2B5 = ["nb2o5"]




#=

end
end

#            println(i, " s1 " , SUB1[end])
end
end
#generate all double subs
SUB2 = []
println("generate2")
for i in 1:3
if i == 1
others = [2,3]
elseif i == 2
others = [1,3]
elseif i == 3
others = [1,2]
end


for s1 in sl[others[1]]
for s2 in sl[others[2]]
temp = [s1,s2,satom[i]]
if arediff(temp)
for o in orders
push!(SUB2,deepcopy(temp[o]))
end
end
#                println(i, " s2 " , SUB2[end])
end
end
end

TRANSMETAL = []
println("generate all trans metal subs")
for i in 1:3
if satom[i] in transmetals
sl[i] = transmetals
end
end
for i in 1:3
if i == 1
others = [2,3]
elseif i == 2
others = [1,3]
elseif i == 3
others = [1,2]
end

for s in sl[i]
temp = [s, satom[others[1]], satom[others[2]]]
if arediff(temp)
for o in orders
push!(TRANSMETAL,deepcopy(temp[o]))
end
end
end

for s in sl[i]
if arediff([s, satom[others[1]], satom[others[2]]])
push!(TRANSMETAL, sort([s, satom[others[1]], satom[others[2]]]))
end
#            println(i, " st " , TRANSMETAL[end])
end


#        for s1 in sl[others[1]]
#            for s2 in sl[others[2]]
#                push!(TRANSMETAL, sort([s1,s2,satom[i]]))
#            end
#        end
end



exact = []
sub1 = []
sub2 = []
subtrans = []
for i in 1:size(summ)[1]
at1 = summ[i,5]
at2 = summ[i,6]
at3 = summ[i,7]
at = summ[i,5:7]

sat = sort([at1, at2, at3])
ex = 0
if at1 in satom
ex += 1
end
if at2 in satom
ex += 1
end
if at3 in satom
ex += 1
end

if ex == 3
exact = [exact; i]
elseif ex == 2
if sat in SUB1
sub1 = [sub1; i]
bestorder(summ[i,5:7], satom)

end
if sat in TRANSMETAL
subtrans = [subtrans; i]
end

elseif ex == 1
if sat in SUB2
sub2 = [sub2; i]
end

end
end

println("EXACT ", length(exact))
for i in exact
println(summ[i,:])
end
println()

println("SUB1  ", length(sub1))
for i in sub1
println(summ[i,:])
end
println()

println("SUB2  ", length(sub2))
for i in sub2
println(summ[i,:])
end
println()

println("SUBTRNS  ", length(subtrans))
for i in subtrans
println(summ[i,:])
end
println()

return summ[exact, :], summ[sub1, :], summ[sub2, :], summ[subtrans, :]
=#

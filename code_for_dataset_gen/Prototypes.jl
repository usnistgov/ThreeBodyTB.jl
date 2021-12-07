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

end


function check_frontier(crys)

    prepare_database(crys)
    database = ThreeBodyTB.ManageDatabase.database_cached
    violation_list, vio_bool = calc_frontier(crys, database, test_frontier=true, verbose=false)
    
    return vio_bool
end


function get_twobody_dist(A,B)

    prepare_database([A,B])
    database = ThreeBodyTB.ManageDatabase.database_cached
    
    ab = database[(A,B)].min_dist * 1.0001
    println("min dist $A $B $ab")
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
    for t1 in crys.types
        for t2 in crys.types
            ab = get_twobody_dist(t1,t2)
            if dmin_types[Set([t1,t2])] < ab
                ret = false
                break
            end
        end
    end
    return ret
end



function setup_proto_data()

    STRUCTDIR = ThreeBodyTB.STRUCTDIR

    CalcD = Dict()


    CalcD["atom"] = ["$STRUCTDIR/atom.in", "scf", "all", "scf", "nscf"]
    CalcD["sc"] = ["$STRUCTDIR/sc.in.up", "vc-relax", "all", "vol-big", "nscf"]
    CalcD["sc_inv"] = ["$STRUCTDIR/sc.in.up", "vc-relax", "all", "break_inv", "nscf"]
    CalcD["bcc"] = ["$STRUCTDIR/bcc.in.up", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["bcc_inv"] = ["$STRUCTDIR/bcc_atom2.in", "vc-relax", "all", "break_inv", "nscf"]
    CalcD["fcc"] = ["$STRUCTDIR/fcc.in.up", "vc-relax", "all", "vol-mid", "nscf"]

    CalcD["line"] = ["$STRUCTDIR/line.in.up", "vc-relax", "z", "vol", "nscf"]
    CalcD["line_rumple"] = ["$STRUCTDIR/line.in.rumple.up", "vc-relax", "z", "vol", "nscf"]


    CalcD["sc_verydense"] = ["$STRUCTDIR/sc.in.up", "vc-relax", "all", "vol-verydense", "nscf"]

    CalcD["fcc_verydense"] = ["$STRUCTDIR/fcc.in.up", "vc-relax", "all", "vol-verydense", "nscf"]
    CalcD["bcc_verydense"] = ["$STRUCTDIR/bcc.in.up", "vc-relax", "all", "vol-verydense", "nscf"]
    CalcD["diamond_verydense"] = ["$STRUCTDIR/diamond.in.up", "vc-relax", "all", "vol-verydense", "nscf"]


    CalcD["fcc_dense"] = ["$STRUCTDIR/fcc.in.up", "vc-relax", "all", "vol-dense", "nscf"]
    CalcD["bcc_dense"] = ["$STRUCTDIR/bcc.in.up", "vc-relax", "all", "vol-dense", "nscf"]

    CalcD["hcp"] = ["$STRUCTDIR/hcp.in.up", "vc-relax", "all", "vol", "nscf"]
    CalcD["diamond"] = ["$STRUCTDIR/diamond.in.up", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["graphene"] = ["$STRUCTDIR/fake_graphene.in", "vc-relax", "2Dxy", "2D-mid", "nscf"]
    CalcD["hex"] = ["$STRUCTDIR/hex.in.up", "vc-relax", "2Dxy", "2D-mid", "nscf"]
    CalcD["hex_short"] = ["$STRUCTDIR/hex.in.up", "vc-relax", "2Dxy", "2D-short", "nscf"]
    CalcD["square"] = ["$STRUCTDIR/square.in.up", "vc-relax", "2Dxy", "2D", "nscf"]


    CalcD["dimer"] =       ["$STRUCTDIR/dimer.in", "relax", "2Dxy", "coords", "nscf"]
    CalcD["dimer_short"] = ["$STRUCTDIR/dimer.in", "relax", "2Dxy", "coords-short", "nscf"]
    CalcD["dimer_small"] = ["$STRUCTDIR/dimer_small.in", "relax", "2Dxy", "scf_small", "nscf"]
    CalcD["dimer_super"] = ["$STRUCTDIR/dimer.in", "relax", "2Dxy", "coords_super", "nscf"]
    CalcD["dimer2_super"] = ["$STRUCTDIR/binary/dimer.in", "scf", "2Dxy", "coords_super", "nscf"]
    CalcD["hcp_shape"] = ["$STRUCTDIR/hcp.in.up", "vc-relax", "all", "shape", "nscf"]

    CalcD["hex_2lay"] = ["$STRUCTDIR/hex_2layers.in.up", "vc-relax", "2Dxy", "2D", "nscf"]
    CalcD["bcc_2lay"] = ["$STRUCTDIR/bcc_2layers.in.up", "vc-relax", "2Dxy", "2D", "nscf"]

    CalcD["fcc_huge"] = ["$STRUCTDIR/fcc.in.up", "vc-relax", "all", "vol-huge", "nscf"]

    CalcD["bcc_tet"] = ["$STRUCTDIR/bcc_tet.in", "vc-relax", "all", "vol", "nscf"]


    CalcD["trimer"] =       ["$STRUCTDIR/trimer.in", "none", "2Dxy", "coords_trimer", "nscf"]
    CalcD["trimer_dense"] =       ["$STRUCTDIR/trimer.in", "none", "2Dxy", "coords_trimer_dense", "nscf"]
    CalcD["trimer2"] =       ["$STRUCTDIR/trimer.in2", "none", "2Dxy", "coords_trimer2", "nscf"]
    CalcD["trimer3"] =       ["$STRUCTDIR/trimer.in3", "none", "2Dxy", "coords_trimer3", "nscf"]

    CalcD["trimer_ab2"] =       ["$STRUCTDIR/binary/trimer.in.ab2", "none", "2Dxy", "coords_trimer_ab", "nscf"]
    CalcD["trimer2_ab2"] =       ["$STRUCTDIR/binary/trimer.in2.ab2", "none", "2Dxy", "coords_trimer_ab", "nscf"]

    CalcD["trimer_ba2"] =       ["$STRUCTDIR/binary/trimer.in.ba2", "none", "2Dxy", "coords_trimer_ab", "nscf"]
    CalcD["trimer2_ba2"] =       ["$STRUCTDIR/binary/trimer.in2.ba2", "none", "2Dxy", "coords_trimer_ab", "nscf"]


    CalcD["trimer_ab2_dense"] =       ["$STRUCTDIR/binary/trimer.in.ab2", "none", "2Dxy", "coords_trimer_ab_dense", "nscf"]
    CalcD["trimer2_ab2_dense"] =       ["$STRUCTDIR/binary/trimer.in2.ab2", "none", "2Dxy", "coords_trimer_ab_dense", "nscf"]

    CalcD["trimer_ba2_dense"] =       ["$STRUCTDIR/binary/trimer.in.ba2", "none", "2Dxy", "coords_trimer_ab_dense", "nscf"]
    CalcD["trimer2_ba2_dense"] =       ["$STRUCTDIR/binary/trimer.in2.ba2", "none", "2Dxy", "coords_trimer_ab_dense", "nscf"]


    CalcD["trimer_ba2_big"] =       ["$STRUCTDIR/binary/trimer.in.ba2.big", "none", "2Dxy", "coords_trimer_ab_big", "nscf"]


    CalcD["as_221"] = ["$STRUCTDIR/POSCAR_As_221", "vc-relax", "all", "vol", "nscf"]
    CalcD["as_orth"] = ["$STRUCTDIR/POSCAR_As_ortho", "vc-relax", "all", "vol", "nscf"]
    CalcD["ga_tet"] = ["$STRUCTDIR/POSCAR_ga_tet", "vc-relax", "all", "vol", "nscf"]
    CalcD["ge_wurtz"] = ["$STRUCTDIR/POSCAR_ge_wurtz", "vc-relax", "all", "vol", "nscf"]
    CalcD["pb_r3m"] = ["$STRUCTDIR/POSCAR_pb_r3m_2atom", "vc-relax", "all", "vol", "nscf"]
    CalcD["beta_sn"] = ["$STRUCTDIR/beta_sn.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["n"] = ["$STRUCTDIR/n.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["f"] = ["$STRUCTDIR/f.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["ga"] = ["$STRUCTDIR/POSCAR_ga", "vc-relax", "all", "vol", "nscf"]
    CalcD["bi"] = ["$STRUCTDIR/bi.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["te"] = ["$STRUCTDIR/te.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["in"] = ["$STRUCTDIR/in.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["i2"] = ["$STRUCTDIR/i2.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["li_p6mmm"] = ["$STRUCTDIR/POSCAR_li_p6mmm", "vc-relax", "all", "vol", "nscf"]

    CalcD["sc_shape"] = ["$STRUCTDIR/sc.in.up", "vc-relax", "all", "shape", "nscf"]
    CalcD["diamond_shear"] = ["$STRUCTDIR/diamond.in.up", "vc-relax", "all", "shear", "nscf"]

    CalcD["hh_mono"] = ["$STRUCTDIR/POSCAR_hh_mono", "vc-relax", "all", "vol-mid", "nscf"]


    CalcD["cscl"] = ["$STRUCTDIR/binary/cscl.in", "vc-relax", "all", "vol-big", "nscf"]
    CalcD["cscl_layers"] = ["$STRUCTDIR/binary/cscl_layers.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["cscl_inv"] = ["$STRUCTDIR/binary/cscl.in", "vc-relax", "all", "break_inv", "nscf"]
    CalcD["rocksalt"] = ["$STRUCTDIR/binary/rocksalt.in", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["rocksalt_inv"] = ["$STRUCTDIR/binary/rocksalt.in", "vc-relax", "all", "break_inv", "nscf"]

    CalcD["rocksalt_dense"] = ["$STRUCTDIR/binary/rocksalt.in", "vc-relax", "all", "vol-dense", "nscf"]
    CalcD["hcp_v2_dense"] = ["$STRUCTDIR/binary/hcp.in2", "vc-relax", "all", "vol-dense", "nscf"]
    CalcD["znse_dense"] = ["$STRUCTDIR/binary/znse.in", "vc-relax", "all", "vol-dense", "nscf"]

    CalcD["nias"] = ["$STRUCTDIR/binary/POSCAR_nias_hex", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["227"] = ["$STRUCTDIR/binary/POSCAR_227", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["laga"] = ["$STRUCTDIR/binary/POSCAR_laga_63", "vc-relax", "all", "vol-mid", "nscf"]

    CalcD["distort"] = ["$STRUCTDIR/binary/POSCAR_distort", "vc-relax", "all", "vol", "nscf"]
    CalcD["distort_ab2"] = ["$STRUCTDIR/binary/POSCAR_distort_ab2", "vc-relax", "all", "vol", "nscf"]
    CalcD["distort_ba2"] = ["$STRUCTDIR/binary/POSCAR_distort_ba2", "vc-relax", "all", "vol", "nscf"]

    CalcD["h2_atom"] = ["$STRUCTDIR/binary/h2_atom.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["h1_atom"] = ["$STRUCTDIR/binary/h1_atom.in", "vc-relax", "all", "vol", "nscf"]

##    CalcD["227_BA"] = ["$STRUCTDIR/binary/POSCAR_227_BA", "vc-relax", "all", "vol-mid", "nscf"]


    CalcD["znse_shear"] = ["$STRUCTDIR/binary/znse.in", "vc-relax", "all", "shear", "nscf"]
    CalcD["hbn"] = ["$STRUCTDIR/binary/hbn.in", "vc-relax", "2Dxy", "2D-mid", "nscf"]
    CalcD["znse"] = ["$STRUCTDIR/binary/znse.in", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["dimer2"] = ["$STRUCTDIR/binary/dimer2.in", "relax", "all", "coords", "nscf"]

    CalcD["dimer2_min"] = ["$STRUCTDIR/binary/dimer.in", "none", "all", "coords_min", "nscf"]
    CalcD["tri2_min"] = ["$STRUCTDIR/binary/dimer.in", "none", "all", "coords_min_tri", "nscf"]
    CalcD["dimer_min"] = ["$STRUCTDIR/dimer.in", "none", "all", "coords_min", "nscf"]
    CalcD["tri_min"] = ["$STRUCTDIR/binary/dimer.in", "none", "all", "coords_min_tri", "nscf"]


    CalcD["dimer2_rev"] = ["$STRUCTDIR/binary/dimer_rev.in", "relax", "all", "coords", "nscf"]
    CalcD["square2"] = ["$STRUCTDIR/binary/square.in", "vc-relax", "2Dxy", "2D", "nscf"]
    CalcD["caf2"] = ["$STRUCTDIR/binary/POSCAR_caf2", "vc-relax", "all", "vol-mid", "nscf"]

    CalcD["co2"] = ["$STRUCTDIR/binary/co2.in", "relax", "all", "coords-small", "nscf"]
    CalcD["co2_v2"] = ["$STRUCTDIR/binary/co2_v2.in", "relax", "all", "coords-small", "nscf"]

    CalcD["square_ab2"] = ["$STRUCTDIR/binary/square_ab2.in", "vc-relax", "2Dxy", "2D", "nscf"]
    CalcD["mgf2"] = ["$STRUCTDIR/binary/POSCAR_mgf2", "vc-relax", "all", "vol", "nscf"]
    CalcD["hcp_v2"] = ["$STRUCTDIR/binary/hcp.in2", "vc-relax", "all", "vol-mid", "nscf"]

    CalcD["p2ca3"] = ["$STRUCTDIR/binary/POSCAR_p2ca3", "vc-relax", "all", "vol-mid", "nscf"]

    CalcD["triangle"] = ["$STRUCTDIR/binary/triangle.in", "relax", "all", "coords-small", "nscf"]
    CalcD["triangle2"] = ["$STRUCTDIR/binary/triangle2.in", "relax", "all", "coords-small", "nscf"]


    CalcD["beta_sn2"] = ["$STRUCTDIR/binary/beta_sn.in.up", "vc-relax", "all", "vol", "nscf"]
    CalcD["hbn_real"] = ["$STRUCTDIR/binary/hbn_real.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["znseAAAB"] = ["$STRUCTDIR/binary/znse.in.super.AAAB", "vc-relax", "all", "vol", "nscf"]
    CalcD["rocksaltAAAB"] = ["$STRUCTDIR/binary/rocksalt.in.super.AAAB", "vc-relax", "all", "vol", "nscf"]
    CalcD["znseABBB"] = ["$STRUCTDIR/binary/znse.in.super.ABBB", "vc-relax", "all", "vol", "nscf"]
    CalcD["rocksaltABBB"] = ["$STRUCTDIR/binary/rocksalt.in.super.ABBB", "vc-relax", "all", "vol", "nscf"]
    CalcD["al2o3"] = ["$STRUCTDIR/binary/POSCAR_al2o3", "vc-relax", "all", "vol", "nscf"]
    CalcD["sis2"] = ["$STRUCTDIR/binary/POSCAR_sis2", "vc-relax", "all", "vol", "nscf"]
    CalcD["tio2_rutile"] = ["$STRUCTDIR/binary/POSCAR_tio2_rutile", "vc-relax", "all", "vol", "nscf"]
    CalcD["bi2se3"] = ["$STRUCTDIR/binary/POSCAR_bi2se3", "vc-relax", "all", "vol", "nscf"]
    CalcD["sis_2d"] = ["$STRUCTDIR/binary/sis_2d.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["sio2_224"] = ["$STRUCTDIR/binary/sio2.224.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["ges"] = ["$STRUCTDIR/binary/POSCAR_ges", "vc-relax", "all", "vol", "nscf"]
    CalcD["alf3"] = ["$STRUCTDIR/binary/POSCAR_alf3", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["lipd3"] = ["$STRUCTDIR/binary/POSCAR_lipd3", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["mgf2_v2"] = ["$STRUCTDIR/binary/POSCAR_mgf2_v2", "vc-relax", "all", "vol", "nscf"]
    CalcD["sns"] = ["$STRUCTDIR/binary/POSCAR_sns", "vc-relax", "all", "vol", "nscf"]
    CalcD["sns2"] = ["$STRUCTDIR/binary/POSCAR_sns2", "vc-relax", "all", "vol", "nscf"]
    CalcD["y2o3"] = ["$STRUCTDIR/binary/POSCAR_y2o3", "vc-relax", "all", "vol", "nscf"]
    CalcD["wurtz"] = ["$STRUCTDIR/binary/POSCAR_wurtz", "vc-relax", "all", "vol", "nscf"]
    CalcD["mgcl2"] = ["$STRUCTDIR/binary/POSCAR_mgcl2", "vc-relax", "all", "vol", "nscf"]
    CalcD["mgcl2_tet"] = ["$STRUCTDIR/binary/POSCAR_mgcl2_tet", "vc-relax", "all", "vol", "nscf"]
    CalcD["asna3_2d"] = ["$STRUCTDIR/binary/POSCAR_asna3", "vc-relax", "2Dxy", "2D", "nscf"]
    CalcD["gain3"] = ["$STRUCTDIR/binary/POSCAR_gain3", "vc-relax", "all", "vol", "nscf"]

    CalcD["li3n_hex"] = ["$STRUCTDIR/binary/POSCAR_li3n", "vc-relax", "all", "vol", "nscf"]
    CalcD["nan3"] = ["$STRUCTDIR/binary/POSCAR_nan3", "vc-relax", "all", "vol", "nscf"]
    CalcD["rbo2"] = ["$STRUCTDIR/binary/POSCAR_rbo2", "vc-relax", "all", "vol", "nscf"]



    CalcD["nb2o5"] = ["$STRUCTDIR/binary/POSCAR_nb2o5", "vc-relax", "all", "vol", "nscf"]
    CalcD["sif4"] = ["$STRUCTDIR/binary/POSCAR_sif4", "vc-relax", "all", "vol", "nscf"]
    CalcD["ticl2"] = ["$STRUCTDIR/binary/POSCAR_ticl2", "vc-relax", "all", "vol", "nscf"]

    CalcD["snf4"] = ["$STRUCTDIR/binary/POSCAR_snf4", "vc-relax", "all", "vol", "nscf"]
    CalcD["ga2s3"] = ["$STRUCTDIR/binary/POSCAR_ga2s3", "vc-relax", "all", "vol", "nscf"]

    CalcD["bcc_13"] = ["$STRUCTDIR/binary/POSCAR_bcc_13", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["bcc_31"] = ["$STRUCTDIR/binary/POSCAR_bcc_31", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["mg2si"] = ["$STRUCTDIR/binary/POSCAR_mg2si", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["simg2"] = ["$STRUCTDIR/binary/POSCAR_simg2", "vc-relax", "all", "vol-mid", "nscf"]

    CalcD["fcc_12"] = ["$STRUCTDIR/binary/fcc_12.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["fcc_21"] = ["$STRUCTDIR/binary/fcc_21.in", "vc-relax", "all", "vol", "nscf"]

    CalcD["fcc_conv_ABBB"] = ["$STRUCTDIR/binary/fcc_conv_ABBB.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["fcc_conv_BAAA"] = ["$STRUCTDIR/binary/fcc_conv_BAAA.in", "vc-relax", "all", "vol", "nscf"]

    CalcD["mgb2_12"] = ["$STRUCTDIR/binary/POSCAR_mgb2", "vc-relax", "all", "vol", "nscf"]
    CalcD["mgb2_21"] = ["$STRUCTDIR/binary/POSCAR_mgb2_21", "vc-relax", "all", "vol", "nscf"]

    CalcD["mgb2_AB"] = ["$STRUCTDIR/binary/POSCAR_mgb2_AB", "vc-relax", "all", "vol", "nscf"]
    CalcD["mgb2_BA"] = ["$STRUCTDIR/binary/POSCAR_mgb2_BA", "vc-relax", "all", "vol", "nscf"]

    CalcD["ab2_71"] = ["$STRUCTDIR/binary/POSCAR_ab2_71", "vc-relax", "all", "vol", "nscf"]
    CalcD["ba2_71"] = ["$STRUCTDIR/binary/POSCAR_ba2_71", "vc-relax", "all", "vol", "nscf"]

    CalcD["irn2_38"] = ["$STRUCTDIR/binary/POSCAR_irn2_38", "vc-relax", "all", "vol", "nscf"]
    CalcD["cao2_12"] = ["$STRUCTDIR/binary/POSCAR_cao2_12", "vc-relax", "all", "vol", "nscf"]

    CalcD["p2o5"] = ["$STRUCTDIR/binary/POSCAR_p2o5", "vc-relax", "all", "vol", "nscf"]


    CalcD["znseAABB"] = ["$STRUCTDIR/binary/znse.AABB.in", "vc-relax", "all", "vol", "nscf"]
    CalcD["squareAABB"] = ["$STRUCTDIR/binary/square.in.AABB", "vc-relax", "2Dxy", "2D", "nscf"]


    CalcD["dimer_pair"] = ["$STRUCTDIR/binary/POSCAR_dimer_pair", "vc-relax", "all", "vol", "nscf"]


    CalcD["mgb2_AB-mid"] = ["$STRUCTDIR/binary/POSCAR_mgb2_AB", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["mgb2_BA-mid"] = ["$STRUCTDIR/binary/POSCAR_mgb2_BA", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["mgb2_12-mid"] = ["$STRUCTDIR/binary/POSCAR_mgb2", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["mgb2_21-mid"] = ["$STRUCTDIR/binary/POSCAR_mgb2_21", "vc-relax", "all", "vol-mid", "nscf"]


#    CalcD["i4mmm_tet_AB"] = ["$STRUCTDIR/binary/POSCAR_i4mmm_4atom_AB", "vc-relax", "all", "vol", "nscf"]
#    CalcD["i4mmm_tet_BA"] = ["$STRUCTDIR/binary/POSCAR_i4mmm_4atom_BA", "vc-relax", "all", "vol", "nscf"]



    CalcD["gei2"] = ["$STRUCTDIR/binary/POSCAR_GeI2", "vc-relax", "all", "vol", "nscf"]
    CalcD["gei4_mol"] = ["$STRUCTDIR/binary/gei4_molecule.in", "relax", "all", "coords-small", "nscf"]

    CalcD["ab3_mol"] = ["$STRUCTDIR/binary/ab3_molecule.in", "relax", "all", "coords-small", "nscf"]


    CalcD["tet"] = ["$STRUCTDIR/binary/tet.in", "vc-relax", "all", "vol", "nscf"]


    CalcD["mof6"] = ["$STRUCTDIR/binary/POSCAR_mof6", "vc-relax",  "all", "vol-mid", "nscf"]


    CalcD["ascl5"] = ["$STRUCTDIR/binary/POSCAR_ascl5", "vc-relax",  "2Dxy", "2D", "nscf"]
    CalcD["rocksalt_shape"] = ["$STRUCTDIR/binary/rocksalt.in", "vc-relax", "all", "shape", "nscf"]
    CalcD["rocksalt_2lay"] = ["$STRUCTDIR/binary/rocksalt.in.2lay", "vc-relax", "2Dxy", "2D", "nscf"]

    CalcD["gas"] = ["$STRUCTDIR/binary/POSCAR_gas", "vc-relax", "all", "vol", "nscf"]

    CalcD["hex12"] = ["$STRUCTDIR/binary/hex_trim_12.in", "vc-relax",  "2Dxy", "2D", "nscf"]
    CalcD["hex21"] = ["$STRUCTDIR/binary/hex_trim_21.in", "vc-relax",  "2Dxy", "2D", "nscf"]

    CalcD["hex12a"] = ["$STRUCTDIR/binary/hex_trim_12.in", "vc-relax",  "2Dxy", "2D-mid", "nscf"]
    CalcD["hex21a"] = ["$STRUCTDIR/binary/hex_trim_21.in", "vc-relax",  "2Dxy", "2D-mid", "nscf"]

#
    CalcD["sn2o2"] = ["$STRUCTDIR/binary/POSCAR_sn2o2", "vc-relax", "all", "vol", "nscf"]
    CalcD["te2o6"] = ["$STRUCTDIR/binary/POSCAR_te2o6", "vc-relax", "all", "vol", "nscf"]
    CalcD["a1f6"] = ["$STRUCTDIR/binary/POSCAR_a1f6", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["a2f6"] = ["$STRUCTDIR/binary/POSCAR_a2f6", "vc-relax", "all", "vol", "nscf"]
    CalcD["ca1pd5"] = ["$STRUCTDIR/binary/POSCAR_ca1pd5", "vc-relax", "all", "vol", "nscf"]
    CalcD["anatase"] = ["$STRUCTDIR/binary/POSCAR_anatase", "vc-relax", "all", "vol", "nscf"]
    CalcD["ga4sr"] = ["$STRUCTDIR/binary/POSCAR_ga4sr", "vc-relax", "all", "vol", "nscf"]
    CalcD["bi1f5"] = ["$STRUCTDIR/binary/POSCAR_bi1f5", "vc-relax", "all", "vol", "nscf"]
    CalcD["cu2o1"] = ["$STRUCTDIR/binary/POSCAR_cu2o1", "vc-relax", "all", "vol", "nscf"]
    CalcD["nb1f4"] = ["$STRUCTDIR/binary/POSCAR_nb1f4", "vc-relax", "all", "vol", "nscf"]
    CalcD["rb1in4"] = ["$STRUCTDIR/binary/POSCAR_rb1in4", "vc-relax", "all", "vol", "nscf"]
    CalcD["ni4la2"] = ["$STRUCTDIR/binary/POSCAR_ni4la2", "vc-relax", "all", "vol", "nscf"]
    CalcD["re3n3"] = ["$STRUCTDIR/binary/POSCAR_re3n3", "vc-relax", "all", "vol", "nscf"]
    CalcD["re3n"] = ["$STRUCTDIR/binary/POSCAR_re3n", "vc-relax", "all", "vol", "nscf"]
    CalcD["cab6"] = ["$STRUCTDIR/binary/POSCAR_cab6", "vc-relax", "all", "vol", "nscf"]

    CalcD["c2ca"] = ["$STRUCTDIR/binary/POSCAR_c2ca", "vc-relax", "all", "vol", "nscf"]
    CalcD["osn2"] = ["$STRUCTDIR/binary/POSCAR_osn2", "vc-relax", "all", "vol", "nscf"]
    CalcD["pt2o2"] = ["$STRUCTDIR/binary/POSCAR_pt2o2", "vc-relax", "all", "vol", "nscf"]
    CalcD["p3n5"] = ["$STRUCTDIR/binary/POSCAR_p3n5", "vc-relax", "all", "vol", "nscf"]
    CalcD["pd3te"] = ["$STRUCTDIR/binary/POSCAR_pd3te", "vc-relax", "all", "vol", "nscf"]
    CalcD["bpt"] = ["$STRUCTDIR/binary/POSCAR_bpt", "vc-relax", "all", "vol", "nscf"]
    CalcD["pd3s"] = ["$STRUCTDIR/binary/POSCAR_JVASP-pd3s", "vc-relax", "all", "vol", "nscf"]


    CalcD["bcc_5lay"] = ["$STRUCTDIR/binary/POSCAR_bcc_5lay", "vc-relax", "all", "vol", "nscf"]
    CalcD["fcc_5lay"] = ["$STRUCTDIR/binary/POSCAR_fcc_5lay", "vc-relax", "all", "vol", "nscf"]


    CalcD["k2n6"] = ["$STRUCTDIR/binary/POSCAR_k2n6", "vc-relax", "all", "vol", "nscf"]

    CalcD["mol_h2o"] =    ["$STRUCTDIR/binary/POSCAR_mol_h2o", "relax", "all", "coords-small", "nscf"]
    CalcD["mol_h2o_v1"] = ["$STRUCTDIR/binary/POSCAR_mol_h2o_v1", "relax", "all", "coords-small", "nscf"]
    CalcD["mol_h2o_v2"] = ["$STRUCTDIR/binary/POSCAR_mol_h2o_v2", "relax", "all", "coords-small", "nscf"]


    #ternary
    CalcD["abc_line"] = ["$STRUCTDIR/ternary/abc_line.in", "relax", "all", "coords-small2", "nscf"]
    CalcD["bac_line"] = ["$STRUCTDIR/ternary/bac_line.in", "relax", "all", "coords-small2", "nscf"]
    CalcD["cab_line"] = ["$STRUCTDIR/ternary/cab_line.in", "relax", "all", "coords-small2", "nscf"]

    CalcD["fcc_tern"] = ["$STRUCTDIR/ternary/fcc_tern.in", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["hex_trim"] = ["$STRUCTDIR/ternary/hex_trim_3.in", "vc-relax",  "2Dxy", "2D-mid", "nscf"]

    CalcD["hh1"] = ["$STRUCTDIR/ternary/POSCAR_hh1", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["hh2"] = ["$STRUCTDIR/ternary/POSCAR_hh2", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["hh3"] = ["$STRUCTDIR/ternary/POSCAR_hh3", "vc-relax", "all", "vol-mid", "nscf"]

    CalcD["stuffhex_1"] = ["$STRUCTDIR/ternary/POSCAR_mgb2_1", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["stuffhex_2"] = ["$STRUCTDIR/ternary/POSCAR_mgb2_2", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["stuffhex_3"] = ["$STRUCTDIR/ternary/POSCAR_mgb2_3", "vc-relax", "all", "vol-mid", "nscf"]

    CalcD["stuffhex_z_1"] = ["$STRUCTDIR/ternary/POSCAR_z_mgb2_1", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["stuffhex_z_2"] = ["$STRUCTDIR/ternary/POSCAR_z_mgb2_2", "vc-relax", "all", "vol-mid", "nscf"]
    CalcD["stuffhex_z_3"] = ["$STRUCTDIR/ternary/POSCAR_z_mgb2_3", "vc-relax", "all", "vol-mid", "nscf"]

    CalcD["rocksalt_2lay_abo2"] = ["$STRUCTDIR/ternary/rocksalt.in.2lay_abo2", "vc-relax", "2Dxy", "2D_tern", "nscf"]
    CalcD["p4mmm"] = ["$STRUCTDIR/ternary/POSCAR_p4mmm", "vc-relax", "all", "3D_tern", "nscf"]
    CalcD["caf2_abc"] = ["$STRUCTDIR/ternary/POSCAR_caf2_ABC", "vc-relax", "all", "vol", "nscf"]

    CalcD["perov"] = ["$STRUCTDIR/ternary/POSCAR_abo3", "vc-relax", "all", "vol", "nscf"]
    CalcD["perov2"] = ["$STRUCTDIR/ternary/POSCAR_abo3_2", "vc-relax", "all", "vol", "nscf"]
    CalcD["perov3"] = ["$STRUCTDIR/ternary/POSCAR_abo3_3", "vc-relax", "all", "vol", "nscf"]
    CalcD["perov4"] = ["$STRUCTDIR/ternary/POSCAR_abo3_4", "vc-relax", "all", "vol", "nscf"]
    CalcD["perov5"] = ["$STRUCTDIR/ternary/POSCAR_abo3_5", "vc-relax", "all", "vol", "nscf"]
    CalcD["perov6"] = ["$STRUCTDIR/ternary/POSCAR_abo3_6", "vc-relax", "all", "vol", "nscf"]


    CalcD["simple_hex"] = ["$STRUCTDIR/simple_hex.in", "vc-relax", "all", "flyaway", "nscf"]

    CalcD["trimer_tern"] =       ["$STRUCTDIR/ternary/POSCAR_trimer_tern", "none", "2Dxy", "trimer_tern", "nscf"]
    CalcD["trimer_tern_right"] =       ["$STRUCTDIR/ternary/POSCAR_trimer_tern_right", "none", "2Dxy", "trimer_tern_right", "nscf"]
    CalcD["trimer_tern_angle"] =       ["$STRUCTDIR/ternary/POSCAR_trimer_tern_angle", "none", "2Dxy", "trimer_tern_angle", "nscf"]

    CalcD["trimer_tern_line"] =       ["$STRUCTDIR/ternary/POSCAR_trimer_tern_line", "none", "2Dxy", "trimer_tern_line", "nscf"]


    CalcD["hh_oxygen"] = ["$STRUCTDIR/ternary/POSCAR_hh_oxygen", "vc-relax", "all", "vol-oxygen", "nscf"]
    CalcD["hh_oxygen2"] = ["$STRUCTDIR/ternary/POSCAR_hh_oxygen2", "vc-relax", "all", "vol-oxygen", "nscf"]
    CalcD["hex_oxygen"] = ["$STRUCTDIR/ternary/hex_trim_3.in_oxygen", "vc-relax",  "2Dxy", "2D-oxygen", "nscf"]

    CalcD[""] = ["$STRUCTDIR/ternary/hex_trim_3.in_oxygen", "vc-relax",  "2Dxy", "2D-oxygen", "nscf"]






    core_mono = [     "sc", "atom",     "sc_inv",     "bcc",     "bcc_inv",     "fcc",     "hcp", "hcp_shape",      "diamond",     "graphene",     "hex",     "square",     "dimer" ,"tri_min", "dimer_min", "trimer", "trimer2", "bcc_5lay", "fcc_5lay",  "hex_2lay", "bcc_2lay", "fcc_dense", "bcc_dense", "znse_dense"]




    core_binary = [    "cscl",          "hbn",     "rocksalt",  "rocksalt_inv",   "znse",     "dimer2", "dimer2_min", "tri2_min", "square2", "hcp_v2", "rocksalt_shape", "rocksalt_2lay", "cscl_layers", "fcc_12", "fcc_21", "bcc_13", "bcc_31", "rocksalt_dense", "hcp_v2_dense", "znse_dense",  "distort", "227", "znseAABB"]

    A0 = [   "as_orth",    "ga_tet",    "ge_wurtz",    "pb_r3m",    "beta_sn",   "ga",    "bi",    "te",    "in",    "i2",    "li_p6mmm",    "hcp_shape",   "n" ] #"bcc_tet.in",  same as POSCAR_ga_tet   #"as_221",  is simple cubic
    
    A1B1 = ["hbn_real", "sis_2d", "ges", "sns", "wurtz"] #nias #
    A1B2 = ["mgcl2", "mgcl2_tet", "caf2", "sis2", "tio2_rutile", "co2",  "mgf2",  "ticl2", "gei2", "anatase"]   # "sns2" duplictes ticl2 #"mgf2_v2" is rutile
    A1B3 = ["alf3", "asna3_2d", "ab3_mol", "gain3", "li3n_hex"]
    A1B4 = ["sif4", "snf4", "gei4_mol"]
    A1B5 = ["ascl5", "bi1f5"]
    A1B6 = ["mof6", "a1f6"]
    A2B3 = ["y2o3", "p2ca3", "al2o3", "bi2se3", "ga2s3", "gas"]
    A2B5 = ["nb2o5", "p2o5"]
    A3B5 = [] #"p3n5"

    short_bonds = ["dimer_pair", "cao2_12", "irn2_38"]
    

#    metals = [ "mg2si", "simg2",  "mgb2_12", "mgb2_21", "fcc_conv_ABBB","fcc_conv_BAAA", "ab2_71", "ba2_71"]
#    metals = [ "mg2si", "simg2",  "mgb2_12", "mgb2_21",  "fcc_conv_ABBB","fcc_conv_BAAA",  "ab2_71", "ba2_71"]
    metals = [ "mg2si", "simg2",  "mgb2_12", "mgb2_21",    "ab2_71", "ba2_71"]

   # all_ternary = ["abc_line", "bac_line", "cab_line", "fcc_tern", "hex_trim", "hh1", "hh2", "hh3", "stuffhex_1", "stuffhex_2", "stuffhex_3","stuffhex_z_1", "stuffhex_z_2", "stuffhex_z_3", "rocksalt_2lay_abo2", "caf2_abc", "perov", "perov2", "perov3",  "perov4",  "perov5",  "perov6"  ]
    core_ternary = ["abc_line", "bac_line", "cab_line", "fcc_tern", "hex_trim", "hh1", "hh2", "hh3", "stuffhex_1", "stuffhex_2", "stuffhex_3", "trimer_tern","trimer_tern_right","trimer_tern_line", "trimer_tern_angle", "p4mmm" ]

    pd = proto_data(CalcD, core_mono, core_binary, A0, A1B1, A1B2, A1B3, A1B4, A1B5, A1B6, A2B3, A2B5,A3B5, metals, short_bonds, core_ternary)

    return pd

end

function  do_run(pd, T1, T2, T3, tmpname, dir, procs, torun; nscf_only = false, only_kspace = false, check_only = false)

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

        file, scf, free, newst, calc_mode = pd.CalcD[st]

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
        elseif newst == "vol-mid"
            ncalc = length( [ 0.9 0.95 1.0 1.05 1.1 ])
        elseif newst == "vol-oxygen"
            ncalc = 4
        elseif newst == "2D-oxygen"
            ncalc = 4
        elseif newst == "vol-dense"
            ncalc = length( [0.8 0.83 0.87 ])
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
            ncalc = length( [0.90 0.95 1.0 1.05 1.10])
        elseif newst == "2D-mid"
            ncalc = length( [0.86 0.88 0.91 0.96 1.0 1.05 1.10])
        elseif newst == "2D-short"
            ncalc = length( [0.80 0.83])
        elseif newst == "shape"
            ncalc = length( [-0.06 -0.03 0.03 0.06])
        elseif newst == "coords"
            ncalc = length( [-0.20 -0.17 -0.14 -0.10 -0.07 -0.03 0.0 0.03 0.07 0.10 0.15 0.2 0.25 0.35 0.5])
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
            ncalc = length( [ -0.12  -0.06 0.0  0.06   ])
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
        elseif newst == "coords_trimer3"
            ncalc = 5*4*2
        elseif (newst == "coords_trimer_dense" || newst == "coords_trimer_ab_dense" )
            ncalc = length([1.0, 0.95])
        elseif newst == "trimer_tern"
            ncalc = 4
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
            d="$dir/$name"*"_vnscf_"*"$newst"*"_"*"$n"        
            if isfile(d*"/projham_K.xml") ||  isfile(d*"/projham_K.xml.gz")
                push!(already_done, d)
            else
                push!(not_done, d)
            end
        end
        if check_only
            continue
        end

        d="$dir/$name"*"_vnscf_"*"$newst"*"_"*"$ncalc"        
        
        if isfile(d*"/projham_K.xml") ||  isfile(d*"/projham_K.xml.gz")
            println("everything is already done! $d")
            println("we can move on")
            continue
        else
            println("not done yet")
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
                
                dft_ref = ThreeBodyTB.DFT.runSCF(c, inputstr=name, nprocs=procs, prefix="$name.qe.relax", directory="$dir", tmpdir="/$tmpname/$name.$randi", wannier=false, code="QE", skip=true, calculation=scf, dofree=free, cleanup=true)

                println("did dft, get new struct")

                if (scf == "relax" || scf == "vc-relax") && maximum(abs.(dft_ref.stress)) > 2e-5
                    println()
                    println("structure convergence not reached, relax again")
                    c2 = dft_ref.crys
                    println(c2)
                    println()

                    dft_ref = ThreeBodyTB.DFT.runSCF(c2, inputstr=name, nprocs=procs, prefix="$name.qe.relax", directory="$dir", tmpdir="/$tmpname/$name.$randi", wannier=false, code="QE", skip=true, calculation=scf, dofree=free, cleanup=true)
                end

                println("END DFT.runSCF")
                println(dft_ref)

                cnew = deepcopy(dft_ref.crys)
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
                for x in [0.8 0.83 0.87 ]
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
                for x in [0.90 0.95 1.0 1.05 1.10]
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
                for x in [1.05, 1.1,  1.2,  1.3]

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
                a = min_dimer_dist_dict[T1]
                for x in [1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
                    c = deepcopy(cnew)
                    c.A[1,1] = a / (0.4^2 + 0.2^2)^0.5 * x
                    c.A[2,2] = a * 1.5 * x
                    c.A[3,3] = a / 0.4 * x
                    push!(torun, deepcopy(c))
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
                for x in [-0.20 -0.17 -0.14 -0.10 -0.07 -0.03 0.0 0.03 0.07 0.10 0.15 0.2 0.25 0.35 0.5]
                    c = deepcopy(cnew)
                    c.coords = c.coords * (1+x)
                    push!(torun, deepcopy(c))
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
                for x in [ -0.12,  -0.06, 0.0, 0.06  ]
                    c = deepcopy(cnew)
                    c.coords[:,3] = c.coords[:,3] .- c.coords[1,3]
#                    c.coords = c.coords * (1+x)
#                    if c.coords[2,3] > 0.35
#                        c.coords[2,3] = 0.35
#                    end
#                    if c.coords[3,3] < -0.35
#                        c.coords[3,3] = -0.35
#                    end
                    c.A[3,3] = c.A[3,3] * (1+x)

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
                        if counter >= 4
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

                    d="$dir/$name"*"_vnscf_"*"$newst"*"_"*"$i"        
                    println(d)
                    println(c)
                    if isfile(d*"/projham_K.xml") ||  isfile(d*"/projham_K.xml.gz")
                        println("continue")
                        continue
                    end

                    dft = ThreeBodyTB.DFT.runSCF(c, nprocs=procs, prefix="qe", directory="$d", tmpdir="$d", wannier=false, code="QE", skip=true, cleanup=true)
                    if calc_mode == "nscf"

                        try
                            tbc, tbck = ThreeBodyTB.AtomicProj.projwfc_workf(dft, nprocs=procs, directory=d, skip_og=true, skip_proj=true, freeze=true, localized_factor = 0.15, cleanup=true, only_kspace=true)
                        catch err3
                            println("err3")
                            println(err3)
                            println("skip nscf $i ")##
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

function do_run_ternary_sub(at1, at2, at3, dir,procs, n1=6, n2 = 12)

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
            push!(torun, deepcopy(c) * 0.95)
            push!(DONE, name)
        end
        n += 1
        if n <= n1
            push!(torun, deepcopy(c) * 1.05)
            push!(DONE, name)
        end
        if n >= n2
            break
        end
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
        if n <= n1
            push!(torun, deepcopy(c) * 0.95)
            push!(DONE, name)
        end
        n += 1
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
            if name in DONE
                continue
            end
            n += 1
            push!(torun, deepcopy(c) )
            push!(DONE, name)
            if n >= n2
                break
            end

        end            
    end
    println("DONE")
    for d in DONE
        println(d)
    end

    for (i,c) in enumerate(torun)
        name = DONE[i]

        println("running dft ternary")
        println(c)
        println()
        d="$dir/$name"*"_vnscf_"*"$i"        
        try
            dft = ThreeBodyTB.DFT.runSCF(c, nprocs=procs, prefix="qe", directory="$d", tmpdir="$d", wannier=false, code="QE", skip=true, cleanup=true)
            tbc, tbck = ThreeBodyTB.AtomicProj.projwfx_workf(dft, nprocs=procs, directory=d, skip_og=true, skip_proj=true, freeze=true, localized_factor = 0.15, cleanup=true, only_kspace=true)
        catch
            println("err dft $d")
        end
        println("done run")


    end

    return torun

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
    keep = []
    
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
    if length(keep) == 1
        if maximum(atom_prefered_oxidation[atom1]) > 0 || maximum(atom_prefered_oxidation[atom2]) > 0
            push!(keep, [atom1, atom2, :metals])
        end
    end


    #special configurations


    if atom1 in ["B"] &&  atom2 in ["H", "Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", "Y", "La", "Sc", "Tl", "In", "Ga"]
        push!(keep, [atom2, atom1, "cab6"])
    end
    if atom2 in ["B"] &&  atom1 in ["H", "Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", "Y", "La", "Sc", "Tl", "In", "Ga"]
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

    if atom1 in ["B", "Al", "Ga", "In", "Tl", "Sc", "Y", "La", "Sb", "Bi", "Co", "Fe", "Ni", "Mn", "Cr", "Ti", "V"] &&  atom2 in ["F"]
        push!(keep, [atom1, atom2, "a2f6"])
    end 
    if atom2 in ["B", "Al", "Ga", "In", "Tl", "Sc", "Y", "La", "Sb", "Bi", "Co", "Fe", "Ni", "Mn", "Cr", "Ti", "V"] &&  atom1 in ["F"]
        push!(keep, [atom2, atom1, "a2f6"])
    end 

    if atom1 in ["Cr", "Mo", "W", "Mn", "Tc", "Re", "Ru", "S", "Se", "Te"] &&  atom2 in ["F"]
        push!(keep, [atom1, atom2, "a1f6"])
        push!(keep, [atom1, atom2, "mof6"])
    end 
    if atom2 in ["Cr", "Mo", "W", "Mn", "Tc", "Re", "Ru", "Os", "S", "Se", "Te"] &&  atom1 in ["F"]
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

    if atom1 in anions && !(atom2 in anions)
        push!(keep, [atom2, atom1, "bcc_5lay"])
        push!(keep, [atom2, atom1, "fcc_5lay"])
    end
    if atom2 in anions && !(atom1 in anions)
        push!(keep, [atom1, atom2, "bcc_5lay"])
        push!(keep, [atom1, atom2, "fcc_5lay"])
    end

    if atom1 in ["Re", "Tc", "Ru", "Os", "Rh", "Ir", "Pd", "Pt", "Ag", "Au", "Ni"] && atom2 in ["Te", "Se", "S", "Be", "B"]
        push!(keep, [atom1, atom2, "pd3s"])
    end
    if atom2 in ["Re", "Tc", "Ru", "Os", "Rh", "Ir", "Pd", "Pt", "Ag", "Au", "Ni"] && atom1 in ["Te", "Se", "S", "Be", "B"]
        push!(keep, [atom2, atom1, "pd3s"])
    end

    bigmetals = ["Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba", "Sc", "Y", "La", "Hg", "H"]

    if (atom1 in bigmetals && !(atom2 in stronganions)) || (atom2 in bigmetals && !(atom1 in stronganions)) 
        push!(keep, [atom1, atom2, "hh_oxygen"])
        push!(keep, [atom1, atom2, "hh_oxygen2"])
        push!(keep, [atom1, atom2, "hex_oxygen"])
    end

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

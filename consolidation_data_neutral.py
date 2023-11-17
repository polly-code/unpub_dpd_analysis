import numpy as np
import pandas as pd
import glob
import sys


def msd_slope(path):
    msd = np.loadtxt(path)
    slope1, d1 = np.polyfit(np.log10(msd[:5, 0]), np.log10(msd[:5, 1]), 1)
    slope2, d2 = np.polyfit(np.log10(msd[300:500, 0]), np.log10(msd[300:500, 1]), 1)
    return slope1, 10 ** d1, slope2, 10 ** d2


def pc_slope(path):
    pc = pd.read_csv(path)
    temp_lim = pc.count()["in"]
    slope_in, _ = np.polyfit(
        np.log10(pc["in_x"].values[:temp_lim]), np.log10(pc["in"].values[:temp_lim]), 1
    )
    temp_lim = pc.count()["out"]
    slope_out, _ = np.polyfit(
        np.log10(pc["out_x"].values[:temp_lim]),
        np.log10(pc["out"].values[:temp_lim]),
        1,
    )
    temp_lim = pc.count()["across"]
    slope_across, _ = np.polyfit(
        np.log10(pc["across_x"].values[:temp_lim]),
        np.log10(pc["across"].values[:temp_lim]),
        1,
    )
    return slope_in, slope_out, slope_across


def aver_duration(path):
    ad = pd.read_csv(path)

    # define tail as last quartile
    lim = int(float(len(ad["x"].values)) / 4)
    # calculate integral for the last quartile
    d_in1_sum = np.sum(ad["inside_1"].values[-lim:] * ad["x"].values[-lim:])
    d_ou1_sum = np.sum(ad["outside_1"].values[-lim:] * ad["x"].values[-lim:])
    d_across1_sum = np.sum(ad["across_1"].values[-lim:] * ad["x"].values[-lim:])
    # calculate duration times
    d_in1 = np.average(ad["x"].values, weights=ad["inside_1"].values)
    d_out1 = np.average(ad["x"].values, weights=ad["outside_1"].values)
    d_across1 = np.average(ad["x"].values, weights=ad["across_1"].values)
    d_in2 = np.average(ad["x"].values, weights=ad["inside_2"].values)
    d_out2 = np.average(ad["x"].values, weights=ad["outside_2"].values)
    d_across2 = np.average(ad["x"].values, weights=ad["across_2"].values)
    d_in3 = np.average(ad["x"].values, weights=ad["inside_3"].values)
    d_out3 = np.average(ad["x"].values, weights=ad["outside_3"].values)
    d_across3 = np.average(ad["x"].values, weights=ad["across_3"].values)
    d_ctcf1 = np.average(ad["x_ctcf"].values, weights=ad["ctcf_1"].values)
    d_ctcf2 = np.average(ad["x_ctcf"].values, weights=ad["ctcf_2"].values)
    d_ctcf3 = np.average(ad["x_ctcf"].values, weights=ad["ctcf_3"].values)
    # calculate passage times
    pt_in1 = np.average(ad["x_pt"].values, weights=ad["pt_inside_1"].values)
    pt_out1 = np.average(ad["x_pt"].values, weights=ad["pt_outside_1"].values)
    pt_across1 = np.average(ad["x_pt"].values, weights=ad["pt_across_1"].values)
    pt_in2 = np.average(ad["x_pt"].values, weights=ad["pt_inside_2"].values)
    pt_out2 = np.average(ad["x_pt"].values, weights=ad["pt_outside_2"].values)
    pt_across2 = np.average(ad["x_pt"].values, weights=ad["pt_across_2"].values)
    pt_in3 = np.average(ad["x_pt"].values, weights=ad["pt_inside_3"].values)
    pt_out3 = np.average(ad["x_pt"].values, weights=ad["pt_outside_3"].values)
    pt_across3 = np.average(ad["x_pt"].values, weights=ad["pt_across_3"].values)
    pt_ctcf1 = np.average(ad["x_pt_ctcf"].values, weights=ad["pt_ctcf_1"].values)
    pt_ctcf2 = np.average(ad["x_pt_ctcf"].values, weights=ad["pt_ctcf_2"].values)
    pt_ctcf3 = np.average(ad["x_pt_ctcf"].values, weights=ad["pt_ctcf_3"].values)
    return (
        d_in1,
        d_out1,
        d_across1,
        d_in2,
        d_out2,
        d_across2,
        d_in3,
        d_out3,
        d_across3,
        d_in1_sum,
        d_ou1_sum,
        d_across1_sum,
        pt_in1,
        pt_out1,
        pt_across1,
        pt_in2,
        pt_out2,
        pt_across2,
        pt_in3,
        pt_out3,
        pt_across3,
        d_ctcf1,
        d_ctcf2,
        d_ctcf3,
        pt_ctcf1,
        pt_ctcf2,
        pt_ctcf3,
    )


def aver_dist(path):
    ad = pd.read_csv(path)
    d_in_w = np.average(ad["bins"].values, weights=ad["inside"].values)
    d_out_w = np.average(ad["bins"].values, weights=ad["outside"].values)
    d_across_w = np.average(ad["bins"].values, weights=ad["across"].values)
    d_ctcf_w = np.average(ad["bins"].values, weights=ad["ctcf"].values)
    return d_in_w, d_out_w, d_across_w, d_ctcf_w


def aver_rgyr(path):
    rgyr = pd.read_csv(path)
    rg_in_w = np.average(rgyr["bins"].values, weights=rgyr["inside"].values)
    rg_out_w = np.average(rgyr["bins"].values, weights=rgyr["outside"].values)
    rg_across_w = np.average(rgyr["bins"].values, weights=rgyr["across"].values)
    return rg_in_w, rg_out_w, rg_across_w


joint_data = pd.DataFrame(
    columns=[
        "ctcf",
        "extruder_speed",
        "loading_rate",
        "unloading_rate",
        "msd_alpha1",
        "msd_D1",
        "msd_alpha2",
        "msd_D2",
        "pc_slope_in",
        "pc_slope_out",
        "pc_slope_across",
        "cont_duration_in1",
        "cont_duration_out1",
        "cont_duration_across1",
        "cont_duration_in2",
        "cont_duration_out2",
        "cont_duration_across2",
        "cont_duration_in3",
        "cont_duration_out3",
        "cont_duration_across3",
        "cont_duration_in1_sum",
        "cont_duration_out1_sum",
        "cont_duration_across1_sum",
        "rgyr_in",
        "rgyr_out",
        "rgyr_across",
        "dsts_in",
        "dsts_out",
        "dsts_across",
        "dsts_ctcf",
        # "num_loops",
        "passage_time_in1",
        "passage_time_out1",
        "passage_time_across1",
        "passage_time_in2",
        "passage_time_out2",
        "passage_time_across2",
        "passage_time_in3",
        "passage_time_out3",
        "passage_time_across3",
        "cont_duration_ctcf1",
        "cont_duration_ctcf2",
        "cont_duration_ctcf3",
        "passage_time_ctcf1",
        "passage_time_ctcf2",
        "passage_time_ctcf3",
    ]
)
p1 = sys.argv[1]
paths = glob.glob(p1)
for path in paths:
    path_msd = path + "/mytraj_c.lammpstrj_msd_a.dat"
    path_pc = path + "/pc_in_out_across.zip"
    path_rgyr = path + "/hist_rgyr_in_out_across.zip"
    path_duration = path + "/hist_duration_in_out_across.zip"
    path_dsts = path + "/hist_dists_in_out_across.zip"
    if not (
        glob.glob(path_msd)
        and glob.glob(path_pc)
        and glob.glob(path_rgyr)
        and glob.glob(path_duration)
        and glob.glob(path_dsts)
    ):
        continue
    msd_a1, msd_d1, msd_a2, msd_d2 = msd_slope(path_msd)
    pc_in, pc_out, pc_across = pc_slope(path_pc)
    (
        cd_in1,
        cd_out1,
        cd_across1,
        cd_in2,
        cd_out2,
        cd_across2,
        cd_in3,
        cd_out3,
        cd_across3,
        cd_in1_sum,
        cd_out1_sum,
        cd_across1_sum,
        pt_in1,
        pt_out1,
        pt_across1,
        pt_in2,
        pt_out2,
        pt_across2,
        pt_in3,
        pt_out3,
        pt_across3,
        cont_duration_ctcf1,
        cont_duration_ctcf2,
        cont_duration_ctcf3,
        passage_time_ctcf1,
        passage_time_ctcf2,
        passage_time_ctcf3,
    ) = aver_duration(path_duration)
    rg_in, rg_out, rg_across = aver_rgyr(path_rgyr)
    d_in, d_out, d_across, d_ctcf = aver_dist(path_dsts)

    s2 = {
        "ctcf": "off",
        "extruder_speed": 0,
        "loading_rate": 0,
        "unloading_rate": 0,
        "msd_alpha1": msd_a1,
        "msd_D1": msd_d1,
        "msd_alpha2": msd_a2,
        "msd_D2": msd_d2,
        "pc_slope_in": pc_in,
        "pc_slope_out": pc_out,
        "pc_slope_across": pc_across,
        "cont_duration_in1": cd_in1,
        "cont_duration_out1": cd_out1,
        "cont_duration_across1": cd_across1,
        "cont_duration_in2": cd_in2,
        "cont_duration_out2": cd_out2,
        "cont_duration_across2": cd_across2,
        "cont_duration_in3": cd_in3,
        "cont_duration_out3": cd_out3,
        "cont_duration_across3": cd_across3,
        "cont_duration_in1_sum": cd_in1_sum,
        "cont_duration_out1_sum": cd_out1_sum,
        "cont_duration_across1_sum": cd_across1_sum,
        "rgyr_in": rg_in,
        "rgyr_out": rg_out,
        "rgyr_across": rg_across,
        "dsts_in": d_in,
        "dsts_out": d_out,
        "dsts_across": d_across,
        "dsts_ctcf": d_ctcf,
        "num_loops": 0,
        "passage_time_in1": pt_in1,
        "passage_time_out1": pt_out1,
        "passage_time_across1": pt_across1,
        "passage_time_in2": pt_in2,
        "passage_time_out2": pt_out2,
        "passage_time_across2": pt_across2,
        "passage_time_in3": pt_in3,
        "passage_time_out3": pt_out3,
        "passage_time_across3": pt_across3,
        "cont_duration_ctcf1": cont_duration_ctcf1,
        "cont_duration_ctcf2": cont_duration_ctcf2,
        "cont_duration_ctcf3": cont_duration_ctcf3,
        "passage_time_ctcf1": passage_time_ctcf1,
        "passage_time_ctcf2": passage_time_ctcf2,
        "passage_time_ctcf3": passage_time_ctcf3,
    }
    joint_data = joint_data.append(s2, ignore_index=True)
    print(s2)
joint_data.to_csv(
    sys.argv[1] + "consolidation_neutral.zip", index=False, compression="zip"
)

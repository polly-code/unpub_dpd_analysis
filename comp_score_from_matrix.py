import numpy as np
import pandas as pd
import argparse
import sys
import time
import copy


def hp1_version(m, c):
    m_dist_norm = copy.deepcopy(m)
    m_dist_norm = iterative_matr_corr(m_dist_norm, n=50)
    mydiag = []
    for i in range(len(m_dist_norm)):
        mydiag.append(np.nanmean(np.diagonal(m_dist_norm, i)))
    for i in range(len(m_dist_norm)):
        if np.isnan(c[i]):
            m_dist_norm[i, :] = np.nan
            continue
        for j in range(len(m_dist_norm)):
            if i == j:
                m_dist_norm[i, j] = np.nan
                continue
            if np.isnan(c[j]):
                m_dist_norm[i, j] = np.nan
                continue
            if mydiag[abs(j - i)] <= 0:
                m_dist_norm[i, j] = np.nan
            else:
                m_dist_norm[i, j] /= mydiag[abs(j - i)]
    aa = []
    ab = []
    bb = []
    for i in range(len(m_dist_norm)):
        if np.isnan(c[i]):
            continue
        if c[i] >= 0:
            temp_aa = []
            temp_ab = []
        elif c[i] < 0:
            temp_bb = []
            temp_ab = []
        for j in range(len(m_dist_norm)):
            if c[i] >= 0 and c[j] >= 0:
                temp_aa.append(m_dist_norm[i, j])
            elif c[i] < 0 and c[j] < 0:
                temp_bb.append(m_dist_norm[i, j])
            elif (c[i] < 0 and c[j] >= 0) or (c[i] >= 0 and c[j] < 0):
                temp_ab.append(m_dist_norm[i, j])
        if c[i] >= 0:
            aa.append(np.nanmean(temp_aa) / np.nanmean(m_dist_norm[i, :]))
            ab.append(np.nanmean(temp_ab) / np.nanmean(m_dist_norm[i, :]))
        elif c[i] < 0:
            bb.append(np.nanmean(temp_bb) / np.nanmean(m_dist_norm[i, :]))
            ab.append(np.nanmean(temp_ab) / np.nanmean(m_dist_norm[i, :]))
    #    print("Compartmental score:")
    #    print(f"AA - {np.nanmean(aa)}")
    #    print(f"AB - {np.nanmean(ab)}")
    #    print(f"BB - {np.nanmean(bb)}")
    return np.nanmean(aa), np.nanmean(ab), np.nanmean(bb)


def iter_loadtxt(filename, delimiter=",", skiprows=0, dtype=float):
    def iter_func():
        with open(filename, "r") as infile:
            for _ in range(skiprows):
                next(infile)
            for line in infile:
                line = line.rstrip().split(delimiter)
                for item in line:
                    yield dtype(item)
        iter_loadtxt.rowlength = len(line)

    data = np.fromiter(iter_func(), dtype=dtype)
    data = data.reshape((-1, iter_loadtxt.rowlength))
    return data


def make_masks(init_c):
    c = copy.deepcopy(init_c)
    c[c >= 0] = 2
    c[c < 0] = 1
    base_mask = np.outer(c, c)
    for i in range(len(base_mask)):
        base_mask[i, i] = np.nan
    aa = copy.deepcopy(base_mask)
    aa[aa != 4] = np.nan
    aa[aa == 4] = 1
    bb = copy.deepcopy(base_mask)
    bb[bb != 1] = np.nan
    ab = copy.deepcopy(base_mask)
    ab[ab != 2] = np.nan
    ab[ab == 2] = 1
    base_mask[base_mask > 0] = 1
    return (aa, bb, ab, base_mask)


def load_compartment_track():
    p2 = "/tungstenfs/scratch/ggiorget/pavel/code_projects/systematic_chrodyn_simulations/experimental_data/out_8k_merged_3reps_2times/compartments.f23.cis.vecs.tsv"
    swap_compartments = 1
    com = pd.read_csv(p2, sep="\t")
    com = com.rename({"chrom": "chr"}, axis=1)
    sub_com = com[com.chr == "chr1"]
    sub_com = sub_com.assign(
        mid=((sub_com["start"].values + sub_com["end"].values) / 2 / 8000).astype(int)
    )
    p3 = "/tungstenfs/scratch/ggiorget/pavel/code_projects/systematic_chrodyn_simulations/experimental_data/ctcf/ctcf_orientation.csv"
    ctcf = pd.read_csv(p3, sep=",")
    sub_ctcf = ctcf[ctcf.chr == "chr1"]
    sub_ctcf = sub_ctcf.assign(
        mid=((sub_ctcf["start"].values + sub_ctcf["end"].values) / 2 / 8000).astype(int)
    )
    for site in sub_ctcf.mid.values:
        sub_com.loc[sub_com.mid == site, "E1"] = np.nan
    return swap_compartments * sub_com["E1"].values


def iterative_matr_corr(inp_m, n=50):
    """
    Iterative matrix correction
    Inputs:
        inp_m: 2d numpy array - contact matrix
        n: int, number of iterations, default 50
    Outputs:
        out_m: 2d numpy array - corrected contact matrix
    """
    out_m = copy.deepcopy(inp_m)
    for _ in range(n):
        average2 = np.nanmean(out_m, axis=0)
        out_m1 = out_m / average2
        average1 = np.nanmean(out_m, axis=1)
        out_m = out_m1 / average1[:, None]
    return out_m


def normalize_map(m, c):
    m_norm = copy.deepcopy(m)
    mydiag = []
    temp_c = copy.deepcopy(c)
    temp_c = temp_c[~np.isnan(temp_c)]
    n_aa = np.count_nonzero(temp_c >= 0) * (np.count_nonzero(temp_c >= 0) - 1)
    n_bb = np.count_nonzero(temp_c < 0) * (np.count_nonzero(temp_c < 0) - 1)
    n_ab = 2 * np.count_nonzero(temp_c < 0) * (np.count_nonzero(temp_c >= 0))
    for i in range(len(m_norm)):
        mydiag.append(np.nanmean(np.diagonal(m_norm, i)))
    for i in range(len(m_norm)):
        if np.isnan(c[i]):
            m_norm[i, :] = np.nan
            continue
        for j in range(len(m_norm)):
            if i == j:
                m_norm[i, j] = np.nan
                continue
            if mydiag[abs(j - i)] < 0.01:
                m_norm[i, j] = np.nan
                continue
            if mydiag[abs(j - i)] >= 0.01:
                if c[i] >= 0 and c[j] >= 0:
                    m_norm[i, j] = m_norm[i, j] / mydiag[abs(j - i)] / n_aa
                elif c[i] < 0 and c[j] < 0:
                    m_norm[i, j] = m_norm[i, j] / mydiag[abs(j - i)] / n_bb
                elif (c[i] < 0 and c[j] >= 0) or (c[i] >= 0 and c[j] < 0):
                    m_norm[i, j] = m_norm[i, j] / mydiag[abs(j - i)] / n_ab
                else:
                    m_norm[i, j] = np.nan
            else:
                m_norm[i, j] = np.nan
    return m_norm


parser = argparse.ArgumentParser(description="Create a single chain in mol2 format")
parser.add_argument(
    "-p",
    "--path",
    type=str,
    default="./mytraj.lammpstrj_cm.dat",
    help="path to contact matrix from lammps trajectory file",
)
path_to_matrix = parser.parse_args().path

c = load_compartment_track()
# mask_aa, mask_bb, mask_ab, base_mask = make_masks(c)
m = iter_loadtxt(path_to_matrix, delimiter="\t")
r_aa, r_ab, r_bb = hp1_version(m, c)

print(f"{path_to_matrix}\naa {r_aa}\nab {r_ab}\nbb {r_bb}")

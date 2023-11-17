import numpy as np
import sys
import mdtraj as md
import pandas as pd
from scipy.spatial import distance
import matplotlib.pyplot as plt
import glob


def comp_score_calculations(arr, c, rc: float, hlim: int):
    raw_mat = np.zeros((len(arr[0]), len(arr[0])))
    for chromosomes in arr:
        m = distance.cdist(chromosomes, chromosomes, "euclidean")
        raw_mat += np.where(m < rc, 1, 0)
    mat = np.zeros((len(c), len(c)))
    for i in range(len(arr[0])):
        for j in range(max(0, i - hlim), min(len(arr[0]), i + hlim)):
            if i >= len(c):
                ni = i - len(c)
            else:
                ni = i
            if j >= len(c):
                nj = j - len(c)
            else:
                nj = j
            mat[ni, nj] += raw_mat[i, j]

    plt.imshow(np.log(mat), cmap="magma_r")
    plt.colorbar()
    plt.savefig(
        f"cm_{rc}.png",
        dpi=300,
    )
    plt.colorbar()
    plt.close()
    # average = mat.mean(axis=0)
    average = np.nanmean(m, axis=0)
    aa = np.zeros((len(mat), 2))
    ab = np.zeros((len(mat), 2))
    bb = np.zeros((len(mat), 2))
    c_aa = []
    c_ab = []
    c_bb = []

    for i in range(len(mat)):
        if np.isnan(c[i]):
            continue
        for j in range(max(0, i - hlim), min(len(mat), i + hlim)):
            if np.isnan(c[j]):
                continue
            if c[i] != c[j]:
                ab[i, 0] += mat[i, j]
                ab[i, 1] += 1
            elif c[i] >= 0:
                aa[i, 0] += mat[i, j]
                aa[i, 1] += 1
            elif c[i] < 0:
                bb[i, 0] += mat[i, j]
                bb[i, 1] += 1
        if c[i] >= 0 and average[i] > 0:
            c_aa.append(aa[i, 0] / aa[i, 1] / average[i])
            c_ab.append(ab[i, 0] / ab[i, 1] / average[i])
        elif average[i] > 0:
            c_bb.append(bb[i, 0] / bb[i, 1] / average[i])
            c_ab.append(ab[i, 0] / ab[i, 1] / average[i])
    return (np.mean(c_aa), np.mean(c_ab), np.mean(c_bb), mat)


p1 = sys.argv[1]
p2 = "/tungstenfs/scratch/ggiorget/pavel/code_projects/systematic_chrodyn_simulations/experimental_data/out_8k/compartments.f2.cis.vecs.tsv"
swap_compartments = -1
hlim = 1000

com = pd.read_csv(p2, sep="\t")
com = com.rename({"chrom": "chr"}, axis=1)
sub_com = com[com.chr == "chr1"]
atom_types = swap_compartments * sub_com["E1"].values

rcs = np.array([1.5])
df = pd.DataFrame(columns=["aa", "bb", "ab", "rc", "var"])
tr = md.load(p1 + "/mytraj.lammpstrj", top=sys.argv[2])
tr.unitcell_lengths = tr.unitcell_lengths * 10

for i in [1, 2, 5, 10]:
    sub_tr = tr[-10 * i :: 10]
    print(f"length {len(sub_tr)}")
    data_scaled = sub_tr.xyz
    if p1[-1] == "/":
        variable_lammps = p1.split("/")[-2]
    else:
        variable_lammps = p1.split("/")[-1]

    for rc in rcs:
        comp_score_aa, comp_score_ab, comp_score_bb, mat = comp_score_calculations(
            data_scaled, atom_types, rc, hlim
        )
        # save pc to dat file
        #    np.savetxt(
        #        f"{sys.argv[1]}/pc_var_{np.around(variable_lammps,2)}_rc_{np.around(rc, 2)}.dat",
        #        calculate_pc(mat),
        #    )

        df = df.append(
            {
                "aa": comp_score_aa,
                "ab": comp_score_ab,
                "bb": comp_score_bb,
                "rc": rc,
                "var": variable_lammps,
            },
            ignore_index=True,
        )
        print(f"Compartmental score:")
        print(f"\tA-A\t\t{comp_score_aa}")
        print(f"\tB-B\t\t{comp_score_bb}")
        print(f"\tA-B\t\t{comp_score_ab}")
    df.to_csv(sys.argv[1] + "entire_dataset.zip", index=False, compression="zip")

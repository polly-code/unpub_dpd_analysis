import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd


def logbin(data_x, data_y, numdots=20, ptf="Deafault path"):
    """
    Replot data from linear space to logspace with evenly distributed dots.
    Inputs:
        data_x, data_y - arrays of x,y
        numdots - number of values on the future graph
    """

    newdata_y = []
    newdata_x = []
    distr = np.logspace(
        np.log10(min(data_x)), np.log10(max(data_x)), num=numdots * 2 + 1
    )
    for i in range(1, len(distr), 2):
        val = 0
        counter = 0
        for j in range(len(data_x)):
            if data_x[j] > distr[i - 1] and data_x[j] < distr[i + 1]:
                val += data_y[j]
                counter += 1
        if counter > 0:
            newdata_y.append(val / counter)
            newdata_x.append(distr[i])
    return newdata_x, newdata_y


def calc_pc(cm, myrange):
    pc = np.zeros((len(cm), 2))
    for i in range(myrange[0][0], myrange[0][1]):
        for j in range(myrange[1][0], myrange[1][1]):
            if j > i:
                pc[j - i, 1] += cm[i, j]
                pc[j - i, 0] += 1
    res = np.divide(
        pc[:, 1], pc[:, 0], out=np.zeros_like(pc[:, 1]), where=pc[:, 0] != 0
    )
    return res / max(res)


cm = np.loadtxt(sys.argv[1] + "mytraj_c.lammpstrj_cm.dat_s")

in_range = [[450, 551], [450, 551]]
in_tad = calc_pc(cm, in_range)
in_tad[in_tad == 0] = np.nan
in_tad_x, in_tad_y = logbin(np.arange(1, len(in_tad)), in_tad[1:], numdots=30)

out_range1 = [[0, 450], [0, 450]]
out_tad1 = calc_pc(cm, out_range1)
out_range2 = [[551, 1000], [551, 1000]]
out_tad2 = calc_pc(cm, out_range2)
out_tad = np.mean([out_tad1, out_tad2], axis=0)
out_tad[out_tad == 0] = np.nan
out_tad_x, out_tad_y = logbin(np.arange(1, len(out_tad)), out_tad[1:], numdots=30)

across_range1 = [[0, 450], [450, 551]]
across_tad1 = calc_pc(cm, across_range1)
across_range2 = [[450, 551], [551, 1000]]
across_tad2 = calc_pc(cm, across_range1)
across_tad = np.mean([across_tad1, across_tad2], axis=0)
across_tad[across_tad == 0] = np.nan
across_tad_x, across_tad_y = logbin(
    np.arange(1, len(across_tad)), across_tad[1:], numdots=30
)

aver_range = [[0, 1000], [0, 1000]]
aver = calc_pc(cm, aver_range)
aver[aver == 0] = np.nan
aver_x, aver_y = logbin(np.arange(1, len(aver)), aver[1:], numdots=30)

print("Pc path", sys.argv[1])

hic_pc = pd.read_csv(
    "/tungstenfs/scratch/ggiorget/pavel/code_projects/systematic_chrodyn_simulations/bonev_track/Bonev.scaling.1kb.csv"
)
hic_x, hic_y = logbin(
    hic_pc["s_bp"].values / 1000,
    hic_pc["balanced.avg"].values / np.max(hic_pc["balanced.avg"].values),
    numdots=30,
)

plt.figure(figsize=(10, 8))
plt.plot(in_tad_x, in_tad_y / max(in_tad_y), label="in")
plt.plot(out_tad_x, out_tad_y / max(out_tad_y), label="out")
plt.plot(across_tad_x, across_tad_y / max(across_tad_y), label="across")
plt.plot(aver_x, aver_y / max(aver_y), label="average")
plt.plot(hic_x, hic_y / max(hic_y), label="Bonev et.al. 2017")
plt.xlabel("Genome distance, a.u. (8 Kb)")
plt.ylabel("Contact probability")
plt.legend()
plt.xscale("log")
plt.yscale("log")
plt.savefig(sys.argv[1] + "pc_contact_prob.png", dpi=300)
plt.close()

data2out = pd.DataFrame()
data2out["in_x"] = in_tad_x
data2out["in"] = in_tad_y
data2out["out_x"] = out_tad_x
data2out["out"] = out_tad_y
data2out["across_x"] = across_tad_x
data2out["across"] = across_tad_y
data2out.to_csv(sys.argv[1] + "pc_in_out_across.zip", index=False, compression="zip")

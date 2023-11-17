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


cm = np.loadtxt(sys.argv[1])
if len(sys.argv) == 3:
    mat_size = int(sys.argv[2])
else:
    mat_size = 1000
aver_range = [[0, mat_size], [0, mat_size]]
aver = calc_pc(cm, aver_range)
aver[aver == 0] = np.nan
aver_x, aver_y = logbin(np.arange(1, len(aver)), aver[1:], numdots=30)


# hic_pc = pd.read_csv(
#    "/tungstenfs/scratch/ggiorget/pavel/code_projects/systematic_chrodyn_simulations/bonev_track/Bonev.scaling.1kb.csv"
# )
# hic_x, hic_y = logbin(
#    hic_pc["s_bp"].values / 1000,
#    hic_pc["balanced.avg"].values / np.max(hic_pc["balanced.avg"].values),
#    numdots=30,
# )
hic_pc = np.loadtxt(
        "/tungstenfs/scratch/ggiorget/pavel/NGS/hic/data/pia_merged_1_2/output_pavel/hic_results/data/higlass/f1_lin_pc.dat")
#        "/tungstenfs/scratch/ggiorget/pavel/code_projects/systematic_chrodyn_simulations/experimental_data/exp_pc.dat"
#    "/tungstenfs/scratch/ggiorget/pavel/hic/practise/pia_merged_1_2/output_pavel/hic_results/data/higlass/f2_lin_pc.dat"
)
hic_x, hic_y = logbin(
    hic_pc[:, 0] / 16,
    hic_pc[:, 1] / np.max(hic_pc[:, 1]),
    numdots=30,
)

plt.figure(figsize=(10, 8))
plt.plot(aver_x, aver_y / max(aver_y), label="average")
plt.plot(hic_x, hic_y / max(hic_y), label="hic +aux 90min")
plt.xlabel("Genome distance, a.u. (8 Kb)")
plt.ylabel("Contact probability")
plt.legend()
plt.xscale("log")
plt.yscale("log")
plt.savefig(sys.argv[1] + "pc_contact_prob.png", dpi=300)
plt.close()

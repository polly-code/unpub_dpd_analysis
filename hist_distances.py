import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import sys
import pandas as pd
import scipy.stats as ss
import itertools


def distr_bursts(arr, rc_min: int = 1, rc_max: int = 4):
    """
    Calculate distribution of bursts considering contact within rc
    Inputs:
        arr - 1D array of floats: distances, order=time
        rc - float: cutting radius
    Outputs:
        distr - 1D arrays of ints: every subarr corresponds to a specific rc and
        every value corresponds to the contact duration
    """
    start_flag_2 = False
    start_flag_1 = False
    start_flag_0 = False
    rcs = np.arange(rc_min, rc_max)
    distr = [[] for _ in range(len(rcs))]
    list_streaks = np.zeros(len(rcs))
    for i in arr:
        if not start_flag_2 and i > rcs[2]:
            start_flag_2 = True
            start_flag_1 = True
            start_flag_0 = True
            list_streaks[2] = 0
            list_streaks[1] = 0
            list_streaks[0] = 0
        elif not start_flag_1 and i > rcs[1]:
            start_flag_1 = True
            start_flag_0 = True
            list_streaks[1] = 0
            list_streaks[0] = 0
        elif not start_flag_0 and i > rcs[0]:
            start_flag_0 = True
            list_streaks[0] = 0

        if i < rcs[0]:
            list_streaks += 1
        elif i < rcs[1]:
            list_streaks[1:] += 1
            if list_streaks[0] > 0 and start_flag_0:
                distr[0].append(list_streaks[0])
                list_streaks[0] = 0
        elif i < rcs[2]:
            list_streaks[2] += 1
            if list_streaks[0] > 0 and start_flag_0:
                distr[0].append(list_streaks[0])
                list_streaks[0] = 0
            if list_streaks[1] > 0 and start_flag_1:
                distr[1].append(list_streaks[1])
                list_streaks[1] = 0
        else:
            if list_streaks[0] > 0 and start_flag_0:
                distr[0].append(list_streaks[0])
                list_streaks[0] = 0
            if list_streaks[1] > 0 and start_flag_1:
                distr[1].append(list_streaks[1])
                list_streaks[1] = 0
            if list_streaks[2] > 0 and start_flag_2:
                distr[2].append(list_streaks[2])
                list_streaks[2] = 0

    return distr


def distr_passage_time(arr, rc_min: int = 1, rc_max: int = 4):
    """
    Calculate distribution of bursts considering contact within rc
    Inputs:
        arr - 1D array of floats: distances, order=time
        rc - float: cutting radius
    Outputs:
        distr - 1D arrays of ints: every subarr corresponds to a specific rc and
        every value corresponds to the contact duration
    """
    start_flag_2 = False
    start_flag_1 = False
    start_flag_0 = False
    rcs = np.arange(rc_min, rc_max)
    distr = [[] for _ in range(len(rcs))]
    list_streaks = np.zeros(len(rcs))
    for i in arr:
        if not start_flag_0 and i < rcs[0]:
            start_flag_2 = True
            start_flag_1 = True
            start_flag_0 = True
            list_streaks[2] = 0
            list_streaks[1] = 0
            list_streaks[0] = 0
        elif not start_flag_1 and i < rcs[1]:
            start_flag_1 = True
            start_flag_2 = True
            list_streaks[1] = 0
            list_streaks[2] = 0
        elif not start_flag_2 and i < rcs[2]:
            start_flag_2 = True
            list_streaks[2] = 0

        if i > rcs[2]:
            list_streaks += 1
        elif i > rcs[1]:
            list_streaks[:2] += 1
            if list_streaks[2] > 0 and start_flag_2:
                distr[2].append(list_streaks[2])
                list_streaks[2] = 0
        elif i > rcs[0]:
            list_streaks[0] += 1
            if list_streaks[2] > 0 and start_flag_2:
                distr[2].append(list_streaks[2])
                list_streaks[2] = 0
            if list_streaks[1] > 0 and start_flag_1:
                distr[1].append(list_streaks[1])
                list_streaks[1] = 0
        else:
            if list_streaks[0] > 0 and start_flag_0:
                distr[0].append(list_streaks[0])
                list_streaks[0] = 0
            if list_streaks[1] > 0 and start_flag_1:
                distr[1].append(list_streaks[1])
                list_streaks[1] = 0
            if list_streaks[2] > 0 and start_flag_2:
                distr[2].append(list_streaks[2])
                list_streaks[2] = 0

    return distr


def ks2_all(arr: list, num: int):
    st_pv = []
    for i in range(int):
        for j in range(i, 16):
            st, pv = stats.ks_2samp(outside[i], outside[j])
            st_pv.append([st, pv])
            print("KS test between ", i, j, " stat=", st, " p-value=", pv)
    return st_pv


def three_tests(a, b):
    stat, p_Student = ss.ttest_ind(a, b)
    stat, p_MW = ss.mannwhitneyu(a, b)
    stat, p_KS = ss.ks_2samp(a, b)
    return p_Student, p_MW, p_KS


def msd_1d(arr):
    res = np.zeros((len(arr) - 1, 2))
    for i in range(len(arr)):
        for j in range(i + 1, len(arr)):
            res[j - i - 1, 0] += (arr[j] - arr[i]) ** 2
            res[j - i - 1, 1] += 1
    return res[:, 0] / res[:, 1]


f = open(sys.argv[1] + "dsts.dat", "r")
for i in f:
    header = i.split()
    f.close()
    break
d = np.loadtxt(sys.argv[1] + "dsts.dat", skiprows=100000)  # , max_rows=1000)
d = d.T
# print("length of dsts ", len(d))


outside = []
# print("outside")
for i in range(0, 64, 8):
    outside.append(d[i])
for i in range(7, 64, 8):
    outside.append(d[i])
plt.figure(figsize=(8, 6))
for i in range(16):
    h, b = np.histogram(outside[i], bins=50, range=(0, 12))
    plt.plot((b[:-1] + b[1:]) / 2, h, label=str(i))
plt.title("Outside TAD")
plt.xlabel("Distance, a.u.")
plt.ylabel("Count")
plt.legend()
plt.savefig(sys.argv[1] + "hist_out.png", dpi=300)
plt.clf()
# print("outside TAD ", np.mean(outside, axis=1))


inside = []
# print("inside")
for i in range(3, 64, 8):
    inside.append(d[i])
for i in range(4, 64, 8):
    inside.append(d[i])
plt.figure(figsize=(8, 6))
for i in range(16):
    h, b = np.histogram(inside[i], bins=50, range=(0, 12))
    plt.plot((b[:-1] + b[1:]) / 2, h, label=str(i))
plt.title("Inside TAD")
plt.xlabel("Distance, a.u.")
plt.ylabel("Count")
plt.legend()
plt.savefig(sys.argv[1] + "hist_in.png", dpi=300)
plt.clf()
# print("inside TAD ", np.mean(inside, axis=1))


across = []
for i in range(64, 80):
    across.append(d[i])
    # print(header[i])
plt.figure(figsize=(8, 6))
for i in range(16):
    h, b = np.histogram(across[i], bins=50, range=(0, 12))
    plt.plot((b[:-1] + b[1:]) / 2, h, label=str(i))
plt.title("Across TAD")
plt.xlabel("Distance, a.u.")
plt.ylabel("Count")
plt.legend()
plt.savefig(sys.argv[1] + "hist_across.png", dpi=300)
plt.clf()
# print("across CTCF ", np.mean(across, axis=1))

ctcf = []
for i in range(80, 88):
    ctcf.append(d[i])
    # print(header[i])
plt.figure(figsize=(8, 6))
for i in range(8):
    h, b = np.histogram(ctcf[i], bins=50, range=(0, 12))
    if i == 0:
        ctcf2out = h
    else:
        ctcf2out += h
    plt.plot((b[:-1] + b[1:]) / 2, h, label=str(i))
plt.title("Between CTCFs")
plt.xlabel("Distance, a.u.")
plt.ylabel("Count")
plt.legend()
plt.savefig(sys.argv[1] + "hist_ctcf.png", dpi=300)
plt.clf()
# print("CTCF ", np.mean(ctcf, axis=1))


if len(sys.argv) == 3 and sys.argv[2] == "msd":
    plt.figure(figsize=(8, 6))
    msds = []
    for i in range(8):
        msds.append(msd_1d(ctcf[i]))
    msd = msds.mean(axis=0)
    plt.plot((np.arange(len(msd)) + 1) * 10000, msd)
    plt.xlabel("Simulation steps")
    plt.ylabel(r"MSD, a.u.$^2$")
    plt.legend()
    plt.savefig(sys.argv[1] + "msd_ctcf.png", dpi=300)
    plt.clf()

data2out = pd.DataFrame()
data2bp = pd.DataFrame()
plt.figure(figsize=(8, 6))
ini = True
merged = list(itertools.chain(*outside))
data2bp["outside"] = merged
for i in range(16):
    h, b = np.histogram(outside[i], bins=50, range=(0, 12))
    if ini:
        th = h
        ini = False
    else:
        th += h
plt.plot((b[:-1] + b[1:]) / 2, th, label="outside")
data2out = data2out.assign(bins=(b[:-1] + b[1:]) / 2)
data2out = data2out.assign(outside=th)

ini = True
merged = list(itertools.chain(*inside))
data2bp["inside"] = merged
for i in range(16):
    h, b = np.histogram(inside[i], bins=50, range=(0, 12))
    if ini:
        th = h
        ini = False
    else:
        th += h
plt.plot((b[:-1] + b[1:]) / 2, th, label="inside")
data2out = data2out.assign(inside=th)

ini = True
merged = list(itertools.chain(*across))
data2bp["across"] = merged
for i in range(16):
    h, b = np.histogram(across[i], bins=50, range=(0, 12))
    if ini:
        th = h
        ini = False
    else:
        th += h
plt.plot((b[:-1] + b[1:]) / 2, th, label="across")
data2out = data2out.assign(across=th)
plt.xlabel("Distance, a.u.")
plt.ylabel("Count")
plt.legend()
plt.savefig(sys.argv[1] + "hist_in_out_across.png", dpi=300)
plt.clf()

# plt.figure(figsize=(8, 6))
# boxplot = data2bp.boxplot(
#    column=["outside", "inside", "across"], grid=False, rot=45, fontsize=15
# )
# p_Student, p_MW, p_KS = three_tests(data2bp.outside.values, data2bp.inside.values)
# plt.text(
#    1.05,
#    0.8,
#    "in / out p-values:\nSt. t-test %e\nMW test %e\nKS test %e"
#    % (p_Student, p_MW, p_KS),
#    transform=plt.gca().transAxes,
# )
# p_Student, p_MW, p_KS = three_tests(data2bp.outside.values, data2bp.across.values)
# plt.text(
#    1.05,
#    0.6,
#    "out / across p-values:\nSt. t-test %e\nMW test %e\nKS test %e"
#    % (p_Student, p_MW, p_KS),
#    transform=plt.gca().transAxes,
# )
# p_Student, p_MW, p_KS = three_tests(data2bp.across.values, data2bp.inside.values)
# plt.text(
#    1.05,
#    0.4,
#    "in / across p-values:\nSt. t-test %e\nMW test %e\nKS test %e"
#    % (p_Student, p_MW, p_KS),
#    transform=plt.gca().transAxes,
# )
# plt.tight_layout()
# plt.savefig(sys.argv[1] + "bp_hist_in_out_across.png", dpi=300)
# plt.clf()
data2out = data2out.assign(ctcf=ctcf2out)
data2out.to_csv(
    sys.argv[1] + "hist_dists_in_out_across.zip", index=False, compression="zip"
)

data2out = pd.DataFrame()
bot = 1
top = 4
# contact durations
db_outside = [[] for _ in range(int(top - bot))]
for i in range(len(outside)):
    data = distr_bursts(outside[i], rc_min=bot, rc_max=top)
    for j in range(top - bot):
        db_outside[j] += data[j]
db_inside = [[] for _ in range(int(top - bot))]
for i in range(len(inside)):
    data = distr_bursts(inside[i], rc_min=bot, rc_max=top)
    for j in range(top - bot):
        db_inside[j] += data[j]
db_across = [[] for _ in range(int(top - bot))]
for i in range(len(across)):
    data = distr_bursts(across[i], rc_min=bot, rc_max=top)
    for j in range(top - bot):
        db_across[j] += data[j]
db_ctcf = [[] for _ in range(int(top - bot))]
for i in range(len(ctcf)):
    data = distr_bursts(ctcf[i], rc_min=bot, rc_max=top)
    for j in range(top - bot):
        db_ctcf[j] += data[j]
# passage times
pt_outside = [[] for _ in range(int(top - bot))]
for i in range(len(outside)):
    data = distr_passage_time(outside[i], rc_min=bot, rc_max=top)
    for j in range(top - bot):
        pt_outside[j] += data[j]
pt_inside = [[] for _ in range(int(top - bot))]
for i in range(len(inside)):
    data = distr_passage_time(inside[i], rc_min=bot, rc_max=top)
    for j in range(top - bot):
        pt_inside[j] += data[j]
pt_across = [[] for _ in range(int(top - bot))]
for i in range(len(across)):
    data = distr_passage_time(across[i], rc_min=bot, rc_max=top)
    for j in range(top - bot):
        pt_across[j] += data[j]
pt_ctcf = [[] for _ in range(int(top - bot))]
for i in range(len(ctcf)):
    data = distr_passage_time(ctcf[i], rc_min=bot, rc_max=top)
    for j in range(top - bot):
        pt_ctcf[j] += data[j]

for rcut in range(top - bot):
    plt.figure(figsize=(8, 6))
    mybins = np.logspace(0, np.log10(200000), 101)
    h, b = np.histogram(db_outside[rcut], bins=mybins)
    data2out = data2out.assign(x=10 * (b[:-1] + b[1:]) / 2)
    data2out = data2out.assign(**{f"outside_{rcut+bot}": h})
    h = h.astype(float)
    h[h == 0] = np.nan
    plt.plot(10 * (b[:-1] + b[1:]) / 2, h, label="outside")

    h, b = np.histogram(db_inside[rcut], bins=mybins)
    data2out = data2out.assign(**{f"inside_{rcut+bot}": h})
    h = h.astype(float)
    h[h == 0] = np.nan
    plt.plot(10 * (b[:-1] + b[1:]) / 2, h, label="inside")

    h, b = np.histogram(db_across[rcut], bins=mybins)
    data2out = data2out.assign(**{f"across_{rcut+bot}": h})
    h = h.astype(float)
    h[h == 0] = np.nan
    plt.plot(10 * (b[:-1] + b[1:]) / 2, h, label="across")
    plt.xlabel("Duration, simulation steps")
    plt.ylabel("Count")
    plt.legend()
    plt.yscale("log")
    plt.xscale("log")
    plt.savefig(
        sys.argv[1]
        + "diff_ranges_in_out_across_contact_length_rc="
        + str(rcut)
        + ".png",
        dpi=300,
    )
    plt.clf()

    # passage times
    plt.figure(figsize=(8, 6))
    mybins = np.logspace(0, np.log10(200000), 101)
    h, b = np.histogram(pt_outside[rcut], bins=mybins)
    data2out = data2out.assign(x_pt=10 * (b[:-1] + b[1:]) / 2)
    data2out = data2out.assign(**{f"pt_outside_{rcut+bot}": h})
    h = h.astype(float)
    h[h == 0] = np.nan
    plt.plot(10 * (b[:-1] + b[1:]) / 2, h, label="pt_outside")

    h, b = np.histogram(pt_inside[rcut], bins=mybins)
    data2out = data2out.assign(**{f"pt_inside_{rcut+bot}": h})
    h = h.astype(float)
    h[h == 0] = np.nan
    plt.plot(10 * (b[:-1] + b[1:]) / 2, h, label="pt_inside")

    h, b = np.histogram(pt_across[rcut], bins=mybins)
    data2out = data2out.assign(**{f"pt_across_{rcut+bot}": h})
    h = h.astype(float)
    h[h == 0] = np.nan
    plt.plot(10 * (b[:-1] + b[1:]) / 2, h, label="pt_across")
    plt.xlabel("Duration, simulation steps")
    plt.ylabel("Count")
    plt.legend()
    plt.yscale("log")
    plt.xscale("log")
    plt.savefig(
        sys.argv[1] + "diff_ranges_in_out_across_passage_time_rc=" + str(rcut) + ".png",
        dpi=300,
    )
    plt.clf()

# ctcf
for rcut in range(top - bot):
    # contact duration
    plt.figure(figsize=(8, 6))
    mybins = np.logspace(0, np.log10(200000), 101)
    h, b = np.histogram(db_ctcf[rcut], bins=mybins)
    data2out = data2out.assign(x_ctcf=10 * (b[:-1] + b[1:]) / 2)
    data2out = data2out.assign(**{f"ctcf_{rcut+bot}": h})
    h = h.astype(float)
    h[h == 0] = np.nan
    plt.plot(10 * (b[:-1] + b[1:]) / 2, h, label="ctcf")
    plt.xlabel("Duration, simulation steps")
    plt.ylabel("Count")
    plt.legend()
    plt.yscale("log")
    plt.xscale("log")
    plt.savefig(
        sys.argv[1] + "diff_ranges_ctcf_contact_length_rc=" + str(rcut) + ".png",
        dpi=300,
    )
    plt.clf()
    # passage times
    plt.figure(figsize=(8, 6))
    mybins = np.logspace(0, np.log10(200000), 101)
    h, b = np.histogram(pt_ctcf[rcut], bins=mybins)
    data2out = data2out.assign(x_pt_ctcf=10 * (b[:-1] + b[1:]) / 2)
    data2out = data2out.assign(**{f"pt_ctcf_{rcut+bot}": h})
    h = h.astype(float)
    h[h == 0] = np.nan
    plt.plot(10 * (b[:-1] + b[1:]) / 2, h, label="pt_ctcf")
    plt.xlabel("Duration, simulation steps")
    plt.ylabel("Count")
    plt.legend()
    plt.yscale("log")
    plt.xscale("log")
    plt.savefig(
        sys.argv[1] + "diff_ranges_ctcf_passage_time_rc=" + str(rcut) + ".png",
        dpi=300,
    )
    plt.clf()


data2out.to_csv(
    sys.argv[1] + "hist_duration_in_out_across.zip", index=False, compression="zip"
)

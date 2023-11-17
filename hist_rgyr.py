import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import sys
import uuid
import pandas as pd
import itertools
import scipy.stats as ss


def three_tests(a, b):
    stat, p_Student = ss.ttest_ind(a, b)
    stat, p_MW = ss.mannwhitneyu(a, b)
    stat, p_KS = ss.ks_2samp(a, b)
    return p_Student, p_MW, p_KS


def plot_distr(arr, label: str = str(uuid.uuid1()), path: str = "./"):
    plt.figure(figsize=(8, 6))
    for i in range(8):
        h, b = np.histogram(arr[i], bins=50, range=(0.5, 5))
        plt.plot((b[:-1] + b[1:]) / 2, h, label=str(i))
    plt.legend()
    plt.title(label)
    plt.xlabel("Radius of gyration")
    plt.ylabel("Count")
    plt.savefig(path + "hist_rgyr_" + label + ".png", dpi=300)
    plt.clf()


def massiv_ks_2samp_test(arr, n):
    for i in range(n):
        for j in range(i, n):
            st, pv = stats.ks_2samp(arr[i], arr[j])
            print("KS test between ", i, j, " stat=", st, " p-value=", pv)


d = np.loadtxt(sys.argv[1] + "rgyr.dat", skiprows=50000)
d = d.T
print("length of rgyr ", len(d))

outside = []
print("outside")
for i in range(8, 16):
    outside.append(d[i])
plot_distr(outside, label="outside", path=sys.argv[1])
print("outside TAD ", np.mean(outside, axis=1))

inside = []
print("inside")
for i in range(0, 8):
    inside.append(d[i])
plot_distr(inside, label="inside", path=sys.argv[1])
print("inside TAD ", np.mean(inside, axis=1))

across = []
print("across")
for i in range(16, 24):
    across.append(d[i])
plot_distr(across, label="across", path=sys.argv[1])
print("across CTCF ", np.mean(across, axis=1))


def trash():
    plt.figure(figsize=(18, 6))
    for i in range(8):
        plt.plot(
            np.mean(inside[i][:-23].reshape(-1, 1000), axis=1),
            c="red",
            alpha=0.3,
        )
        plt.plot(
            np.mean(outside[i][:-23].reshape(-1, 1000), axis=1),
            c="blue",
            alpha=0.3,
        )
        plt.plot(
            np.mean(across[i][:-23].reshape(-1, 1000), axis=1),
            c="green",
            alpha=0.3,
        )
    plt.xlabel(r"Time, $10^4$ steps")
    plt.ylabel("Radius of gyration")
    plt.plot([], color="red", label="inside", alpha=0.3)
    plt.plot([], color="blue", label="outside", alpha=0.3)
    plt.plot([], color="green", label="across", alpha=0.3)
    plt.legend()
    plt.savefig(sys.argv[1] + "rgyr_t_all.png", dpi=300)
    plt.clf()


data2out = pd.DataFrame()
data2bp = pd.DataFrame()
plt.figure(figsize=(8, 6))
ini = True
merged = list(itertools.chain(*outside))
data2bp["outside"] = merged
for i in range(8):
    h, b = np.histogram(outside[i], bins=50, range=(0.5, 5))
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
for i in range(8):
    h, b = np.histogram(inside[i], bins=50, range=(0.5, 5))
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
for i in range(8):
    h, b = np.histogram(across[i], bins=50, range=(0.5, 5))
    if ini:
        th = h
        ini = False
    else:
        th += h
plt.plot((b[:-1] + b[1:]) / 2, th, label="across")
data2out = data2out.assign(across=th)
plt.legend()
plt.xlabel("Radius of gyration")
plt.ylabel("Count")
plt.savefig(sys.argv[1] + "hist_rgyr_in_out_across.png", dpi=300)
plt.clf()

plt.figure(figsize=(8, 6))
boxplot = data2bp.boxplot(
    column=["outside", "inside", "across"], grid=False, rot=45, fontsize=15
)
# p_Student, p_MW, p_KS = three_tests(data2bp.outside.values, data2bp.inside.values)
# plt.text(
#    1.05,
#    0.8,
#    "out / in p-values:\nSt. t-test %e\nMW test %e\nKS test %e"
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
plt.tight_layout()
plt.savefig(sys.argv[1] + "bp_rgyr_in_out_across.png", dpi=300)
data2out.to_csv(
    sys.argv[1] + "hist_rgyr_in_out_across.zip", index=False, compression="zip"
)

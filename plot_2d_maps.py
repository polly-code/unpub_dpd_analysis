import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

df=pd.read_csv("//tungsten-nas.fmi.ch/tungsten/scratch/ggiorget/pavel/code_projects/systematic_chrodyn_simulations/data/consolidation.zip")
subdf = df[(df["ctcf"]=='on') & (df["loading_rate"] == 0.00001)]
sns.lineplot(data=subdf, x="extruder_speed", y="pc_slope_out", hue="unloading_rate", palette="flare", hue_norm=mpl.colors.LogNorm())

for i in [0.01, 0.001, 0.0001, 0.00001]:
    subdf = df[(df["ctcf"]=='on') & (df["loading_rate"] == i)]
    plt.figure()
    splot = sns.lineplot(data=subdf, x="unloading_rate", y="msd_alpha2", hue="extruder_speed", palette="flare", hue_norm=mpl.colors.LogNorm())
    splot.set(xscale="log")
    splot.set(title=f"loading_rate {i}")

for i in [40, 200, 500, 1000]:
    subdf = df[(df["ctcf"]=='on') & (df["extruder_speed"] == i)]
    toplot = subdf.pivot("loading_rate", "unloading_rate", "msd_alpha1")
    plt.figure()
    ax = sns.heatmap(toplot, vmin=0.35, vmax=0.55)
    ax.set(title=f"extruder speed {i}")
    
for i in [40, 200, 500, 1000]:
    subdf = df[(df["ctcf"]=='on') & (df["extruder_speed"] == i)]
    toplot = subdf.pivot("loading_rate", "unloading_rate", "pc_slope_out")
    plt.figure()
    ax = sns.heatmap(toplot, vmin=-2.3, vmax=0.1, cbar_kws={'label': 'pc_out'})
    ax.set(title=f"extruder speed {i}")
    
for i in [40, 200, 500, 1000]:
    subdf = df[(df["ctcf"]=='on') & (df["extruder_speed"] == i)]
    toplot = subdf.pivot("loading_rate", "unloading_rate", "rgyr_in")
    plt.figure()
    ax = sns.heatmap(toplot, vmin=440000, vmax=675000, cmap='gnuplot2', cbar_kws={'label': 'rgyr_in'})
    ax.set(title=f"extruder speed {i}")

# contact duration in out across rc=1
for i in [40, 200, 500, 1000]:
    subdf = df[(df["ctcf"]=='on') & (df["extruder_speed"] == i)]
    toplot = subdf.pivot("loading_rate", "unloading_rate", "cont_duration_in1")
    plt.figure()
    ax = sns.heatmap(toplot, vmin=10, vmax=45, cbar_kws={'label': 'cont_duration_in1'})
    ax.set(title=f"extruder speed {i}")
    
for i in [40, 200, 500, 1000]:
    subdf = df[(df["ctcf"]=='on') & (df["extruder_speed"] == i)]
    toplot = subdf.pivot("loading_rate", "unloading_rate", "cont_duration_out1")
    plt.figure()
    ax = sns.heatmap(toplot, vmin=10, vmax=45, cbar_kws={'label': 'cont_duration_out1'})
    ax.set(title=f"extruder speed {i}")
    
for i in [40, 200, 500, 1000]:
    subdf = df[(df["ctcf"]=='on') & (df["extruder_speed"] == i)]
    toplot = subdf.pivot("loading_rate", "unloading_rate", "cont_duration_across1")
    plt.figure()
    ax = sns.heatmap(toplot, vmin=10, vmax=45, cbar_kws={'label': 'cont_duration_across1'})
    ax.set(title=f"extruder speed {i}")
    
# contact duration in out across rc=2
for i in [40, 200, 500, 1000]:
    subdf = df[(df["ctcf"]=='on') & (df["extruder_speed"] == i)]
    toplot = subdf.pivot("loading_rate", "unloading_rate", "cont_duration_in2")
    plt.figure()
    ax = sns.heatmap(toplot, vmin=10, vmax=45, cbar_kws={'label': 'cont_duration_in2'})
    ax.set(title=f"extruder speed {i}")
    
for i in [40, 200, 500, 1000]:
    subdf = df[(df["ctcf"]=='on') & (df["extruder_speed"] == i)]
    toplot = subdf.pivot("loading_rate", "unloading_rate", "cont_duration_out2")
    plt.figure()
    ax = sns.heatmap(toplot, vmin=10, vmax=45, cbar_kws={'label': 'cont_duration_out2'})
    ax.set(title=f"extruder speed {i}")
    
for i in [40, 200, 500, 1000]:
    subdf = df[(df["ctcf"]=='on') & (df["extruder_speed"] == i)]
    toplot = subdf.pivot("loading_rate", "unloading_rate", "cont_duration_across2")
    plt.figure()
    ax = sns.heatmap(toplot, vmin=10, vmax=45, cbar_kws={'label': 'cont_duration_across2'})
    ax.set(title=f"extruder speed {i}")
# contact duration in out across rc=3
for i in [40, 200, 500, 1000]:
    subdf = df[(df["ctcf"]=='on') & (df["extruder_speed"] == i)]
    toplot = subdf.pivot("loading_rate", "unloading_rate", "cont_duration_in3")
    plt.figure()
    ax = sns.heatmap(toplot, vmin=10, vmax=45, cbar_kws={'label': 'cont_duration_in3'})
    ax.set(title=f"extruder speed {i}")
    
for i in [40, 200, 500, 1000]:
    subdf = df[(df["ctcf"]=='on') & (df["extruder_speed"] == i)]
    toplot = subdf.pivot("loading_rate", "unloading_rate", "cont_duration_out3")
    plt.figure()
    ax = sns.heatmap(toplot, vmin=10, vmax=45, cbar_kws={'label': 'cont_duration_out3'})
    ax.set(title=f"extruder speed {i}")
    
for i in [40, 200, 500, 1000]:
    subdf = df[(df["ctcf"]=='on') & (df["extruder_speed"] == i)]
    toplot = subdf.pivot("loading_rate", "unloading_rate", "cont_duration_across3")
    plt.figure()
    ax = sns.heatmap(toplot, vmin=10, vmax=45, cbar_kws={'label': 'cont_duration_across3'})
    ax.set(title=f"extruder speed {i}")

variables = df.columns.values[4:]
for var in variables:
    temp_min = min(df[var].values)
    temp_max = max(df[var].values)
    for ctcf in ['on', 'off']:
        for i in [40, 200, 500, 1000]:
            subdf = df[(df["ctcf"]==ctcf) & (df["extruder_speed"] == i)]
            if subdf.empty:
                continue
            toplot = subdf.pivot("loading_rate", "unloading_rate", var)
            
            plt.figure()
            ax = sns.heatmap(toplot, vmin=temp_min, vmax=temp_max, cmap='gnuplot2', cbar_kws={'label': var})
            ax.set(title=f"CTCFs {ctcf}, extruder speed {i}")
            plt.savefig(f"//tungsten-nas.fmi.ch/tungsten/scratch/ggiorget/pavel/code_projects/systematic_chrodyn_simulations/data/consolidation_img/{ctcf}_{i}_{var}.png", dpi=300)
            plt.close()
    
    
    

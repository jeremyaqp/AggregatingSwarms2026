import numpy as np
import matplotlib.pyplot as plt
from glob import glob

data_files = glob("../../Data/Aggregation/*.npy")

data = dict(
    {
        "A115": [],
        "F115": []
    }
)
for f_count in data_files:
    f_lit = f_count.replace("count", "lit")
    arr_count = np.load(f_count)
    arr_lit = np.load(f_lit)[:43199]

    dataKey = f_count.split("/")[-1].split("_")[0]
    if dataKey not in data.keys():
        continue

    data[dataKey].append(arr_lit / 64)


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica",
    "font.size": 14
})

colors = {
    "A115": "#6ca6e9",
    "F115": "#6dc268",
    "MEAN_A115": "#12486e",
    "MEAN_F115": "#1b5c1b"
}

labels = {
    "A115": r"Aligner",
    "F115": r"Fronter"
}

subplots = {
    "A115": 211,
    "F115": 212
}

smooth = 91
MAX_LEN = 42834
plt.rcParams['font.family'] = 'Times New Roman'
plt.figure(figsize=(5, 4.5))
for key in data.keys():
    for arr in data[key]:
        arr_smooth = np.convolve(arr, np.ones(smooth) / smooth, mode="valid")[:MAX_LEN]
        plt.plot(np.arange(len(arr_smooth)) / 6.0, arr_smooth, color=colors[key], linewidth=1.0, alpha=0.6)
        plt.xlabel("Time (s.)", fontsize=22)
        plt.ylabel(r"$N_{\circ} / N$", fontsize=22)
    # MEAN
    mean_arr = np.mean(data[key], axis=0)
    mean_arr = np.convolve(mean_arr, np.ones(smooth) / smooth, mode="valid")[:MAX_LEN]
    plt.plot(np.arange(len(mean_arr)) / 6.0, mean_arr, color=colors["MEAN_" + key], linewidth=1.0, alpha=1.0)
    plt.plot(0, 0, color=colors[key], label=labels[key])
handles, labels = plt.gca().get_legend_handles_labels()
order = [1,0]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.grid()
plt.yticks(np.linspace(0, 0.9, num=10), fontsize=10)
plt.ylim(0, 0.6)
plt.xlim(0, 8000)
#plt.axhline(0.84, color='k', linestyle='--', alpha=0.5, linewidth=2.0)
#plt.annotate("0.84 \n (max)", xy=(8000, 0.84), color='k', fontsize=16, ha='left', va='center')
print("FULL MEAN A115 : ", np.mean(data["A115"]))
print("FULL MEAN F115 : ", np.mean(data["F115"]))
plt.savefig(f"../out/AggregationResults_smooth_merged.png", dpi=300, bbox_inches='tight')
plt.show()

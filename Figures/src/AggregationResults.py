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
    "F0": "r",
    "A0": "b"
}

labels = {
    "A115": r"Aligner \textit{move}",
    "F115": r"Fronter \textit{move}",
    "F0": r"Fronter \textit{still}",
    "A0": r"Aligner \textit{still}"
}

subplots = {
    "A115": 211,
    "F115": 212
}

smooth = 7
plt.rcParams['font.family'] = 'Times New Roman'
plt.figure(figsize=(4.5, 4.5))
for key in data.keys():
    plt.subplot(subplots[key])
    for arr in data[key]:
        arr_smooth = np.convolve(arr, np.ones(smooth) / smooth, mode="valid")[:42834]
        plt.plot(np.arange(len(arr_smooth)) / 6.0, arr_smooth, color=colors[key], linewidth=0.8, alpha=0.8)
        plt.xlabel("Time (s.)", fontsize=22)
        plt.ylabel(r"$N_{\circ} / N_{tot}$", fontsize=22)
    plt.plot(0, 0, color=colors[key], label=labels[key])
    #plt.legend()
    plt.grid()
    plt.yticks(np.linspace(0, 0.9, num=10), fontsize=10)
    plt.ylim(0, 1)
    plt.xlim(0, 8000)
    plt.axhline(0.84, color='k', linestyle='--', alpha=0.5)
    plt.annotate("0.84 \n (max)", xy=(8000, 0.84), color='k', fontsize=16, ha='left', va='center')
    print("FULL MEAN A115 : ", np.mean(data["A115"]))
    print("FULL MEAN F115 : ", np.mean(data["F115"]))
plt.savefig("../out/AggregationResults_smooth7.png", dpi=300, bbox_inches='tight')
plt.show()


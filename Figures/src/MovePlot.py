import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Times New Roman",
        "font.size": 14
    })

    final_df = pd.read_csv("../../Data/simu_data_aggregated.csv")
    
    # drop ratio= 1.0 rows if any
    final_df = final_df[final_df["ratio"] != 1.0]
    np.set_printoptions(threshold=np.inf, linewidth=np.inf)
    print(np.unique(np.array(final_df["ratio"])))
    print(len(final_df["ratio"] ))
    print(len(final_df[final_df["ratio"] == 0.]))

    df_mean = final_df.groupby("ratio", as_index=False)["lit"].median()
    df_mean_alignment = final_df.groupby("ratio", as_index=False)["alignment"].median()

    fig, ax = plt.subplots(figsize=(16, 6))
    ax2 = ax.twinx()  # instantiate a second Axes that shares the same x-axis

    plt.scatter(final_df["ratio"], final_df["lit"] / 64, color="darkred", alpha=0.015, s=40, marker="o")
    plt.scatter(final_df["ratio"], final_df["alignment"], color="darkblue", alpha=0.015, s=40, marker="o")
    plt.plot(df_mean["ratio"], df_mean["lit"] / 64, color="red", linewidth=2, alpha=1.0)
    plt.plot(df_mean["ratio"], df_mean_alignment["alignment"], color="blue", linewidth=2, alpha=1.0)

    plt.axhline(y=0.18, color='black', linestyle='--')
    plt.xlim(-5.0, 5.0)
    plt.ylim(0, 1.0)
    ax.set_xlabel(r"$\epsilon / \tau_n$", fontsize=30)
    ax.set_ylabel(r"$N_\circ / N_\bullet$", fontsize=30, color="darkred")
    ax2.set_ylabel(r'$<\Psi>$', fontsize=30, color="darkblue")  # we already handled the x-label with ax1
    ax.annotate("0.18 (avg)", xy=(-4.9, 0.21), color='k', fontsize=16, ha='left', va='center')
    ax.set_yticks(np.linspace(0, 1.0, num=11))
    ax.set_xticks(np.linspace(-5, 5.0, num=11))
    ax.tick_params(axis='y', labelsize=22, labelcolor="darkred")
    ax.tick_params(axis='x', labelsize=22)
    ax2.set_yticks(np.linspace(0, 1.0, num=11))
    ax2.tick_params(axis='y', labelsize=22, labelcolor="darkblue")
    plt.tight_layout()
    plt.savefig(f"../out/SimuAggregation.png", dpi=300)
    plt.show()

# ======================
# ENTRY POINT
# ======================
if __name__ == "__main__":
    main()


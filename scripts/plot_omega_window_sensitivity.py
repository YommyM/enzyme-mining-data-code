#!/usr/bin/env python3
"""Supporting Figure S2: family detail and active-site-window sensitivity."""

from pathlib import Path
import collections
import csv

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "data" / "omega_reference" / "original_workflow_multilabel_profile_benchmark.csv"
RANGE = ROOT / "data" / "omega_pocket_range_current" / "pocket_range_summary.csv"
SIZES = ROOT / "data" / "omega_pocket_range_current" / "fragment_sizes.csv"
OUT = ROOT / "figures" / "omega_window_sensitivity"
GLOBAL = "#777777"
LOCAL = "#0072B2"

SHORT = {
    "N-acyl diamino-acid aminotransferase": "N-acyl diamino acid",
    "GABA/5-aminovalerate--2-oxoglutarate aminotransferase": "GABA/5AV–2OG",
    "GABA--pyruvate aminotransferase": "GABA–pyruvate",
    "terminal-amino/diamino--2-oxoglutarate aminotransferase": "terminal/diamine–2OG",
    "beta-amino-acid aminotransferase": "β-amino acid",
    "noncarboxylate acidic omega-amino aminotransferase": "phosphonate/sulfonate",
    "S-selective aromatic/aliphatic amine transaminase (Fold I)": "S-selective amine",
    "broad-spectrum amine transaminase (Fold I)": "broad-spectrum amine",
}


def family_values(rows, representation, key):
    values = collections.defaultdict(list)
    for row in rows:
        if row["representation"] == representation:
            values[row["family"]].append(float(row[key]))
    return {family: float(np.mean(items)) for family, items in values.items()}


def bootstrap_ci(values, draws=20000):
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(20260923)
    means = np.mean(rng.choice(values, size=(draws, len(values)), replace=True), axis=1)
    return tuple(np.quantile(means, [0.025, 0.975]))


def main():
    rows = list(csv.DictReader(TABLE.open()))
    range_rows = list(csv.DictReader(RANGE.open()))
    size_rows = list(csv.DictReader(SIZES.open()))
    windows = [3, 5, 7, 10, 15, 20]
    by_window = {int(row["representation"].replace("window_", "")): row for row in range_rows}

    global_ap = family_values(rows, "global", "average_precision")

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7.0, "axes.labelsize": 7.0,
        "axes.titlesize": 8.0, "xtick.labelsize": 6.2, "ytick.labelsize": 5.8,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig, (ax_a, ax_b, ax_c) = plt.subplots(
        1, 3, figsize=(7.15, 3.05), gridspec_kw={"width_ratios": [1.0, 1.05, 0.82]}
    )

    global_mean = float(np.mean(list(global_ap.values())))
    global_low, global_high = bootstrap_ci(list(global_ap.values()))
    local_mean = np.array([float(by_window[w]["average_precision"]) for w in windows])
    local_low = np.array([float(by_window[w]["average_precision_ci_low"]) for w in windows])
    local_high = np.array([float(by_window[w]["average_precision_ci_high"]) for w in windows])
    x = np.arange(1, len(windows) + 1)
    ax_a.errorbar(x, local_mean, yerr=np.vstack([local_mean - local_low, local_high - local_mean]),
                  color=LOCAL, lw=1.45, marker="o", ms=3.8, mec="white", mew=0.4,
                  capsize=2.0, elinewidth=0.8, label="Local windows", zorder=3)
    ax_a.errorbar([0], [global_mean], yerr=[[global_mean - global_low],
                  [global_high - global_mean]], fmt="D", ms=4.1, color=GLOBAL,
                  mfc="white", mec=GLOBAL, capsize=2.2, lw=0.9,
                  label="Whole protein", zorder=4)
    for xi, w in zip(x, windows):
        if float(by_window[w]["p_average_precision"]) < 0.05:
            ax_a.text(xi, local_mean[xi - 1] + 0.018, "*", color=LOCAL,
                      ha="center", va="bottom", fontsize=8, fontweight="bold")
    ax_a.set_xticks(np.arange(0, len(windows) + 1))
    ax_a.set_xticklabels(["Global", *[rf"$\pm${w}" for w in windows]], rotation=35,
                         ha="right", rotation_mode="anchor")
    ax_a.set_ylabel("Group-balanced average precision")
    ax_a.set_ylim(0.34, 0.80)
    ax_a.set_title("A   Average precision across windows", loc="left", fontweight="bold")
    ax_a.grid(axis="y", color="#E5E5E5", lw=0.45)
    ax_a.legend(frameon=False, loc="upper right", fontsize=5.8, handlelength=1.3)
    ax_a.spines[["top", "right"]].set_visible(False)

    # Compare each local window with the family-balanced global mean.
    global_recovery = np.array([
        np.mean(list(family_values(rows, "global", f"recovery_{k}").values()))
        for k in range(1, 6)
    ])
    delta = np.array([[float(by_window[w][f"recovery_{k}"]) - global_recovery[k - 1]
                       for k in range(1, 6)] for w in windows])
    cmap = LinearSegmentedColormap.from_list("gain", ["#B85C4A", "#FAFAFA", "#21618C"])
    lim = max(0.18, float(np.max(np.abs(delta))))
    image = ax_b.imshow(delta, cmap=cmap, vmin=-lim, vmax=lim,
                        interpolation="nearest", aspect="auto")
    ax_b.set_xticks(np.arange(5)); ax_b.set_xticklabels([1, 2, 3, 4, 5])
    ax_b.set_yticks(np.arange(len(windows)))
    ax_b.set_yticklabels([rf"$\pm${w}" for w in windows])
    ax_b.set_xlabel("Number of profile neighbors, k")
    ax_b.set_ylabel("Active-site window")
    ax_b.set_title("B   Recovery gain over global", loc="left", fontweight="bold")
    for i in range(delta.shape[0]):
        for j in range(delta.shape[1]):
            value = delta[i, j]
            ax_b.text(j, i, f"{value:+.2f}", ha="center", va="center", fontsize=4.9,
                      color="white" if abs(value) > lim * 0.58 else "#222222")
    cbar = fig.colorbar(image, ax=ax_b, orientation="horizontal", fraction=0.065,
                        pad=0.18, ticks=[-lim, 0, lim])
    cbar.set_label("Local − global recovery", labelpad=1)
    cbar.ax.tick_params(length=2, pad=1)

    sizes = collections.defaultdict(list)
    for row in size_rows:
        sizes[int(row["window"])].append(float(row["residue_count"]))
    medians = np.array([np.median(sizes[w]) for w in windows])
    low = np.array([np.quantile(sizes[w], 0.25) for w in windows])
    high = np.array([np.quantile(sizes[w], 0.75) for w in windows])
    x = np.arange(len(windows))
    ax_c.fill_between(x, low, high, color=LOCAL, alpha=0.15, lw=0)
    ax_c.plot(x, medians, color=LOCAL, marker="o", ms=3.8, lw=1.5,
              mec="white", mew=0.4)
    ax_c.set_xticks(x); ax_c.set_xticklabels([rf"$\pm${w}" for w in windows], rotation=35,
                                             ha="right", rotation_mode="anchor")
    ax_c.set_xlabel("Active-site window")
    ax_c.set_ylabel("Resolved residues")
    ax_c.set_title("C   Fragment extent", loc="left", fontweight="bold")
    ax_c.grid(axis="y", color="#E5E5E5", lw=0.45)
    ax_c.spines[["top", "right"]].set_visible(False)

    fig.subplots_adjust(left=0.09, right=0.985, top=0.88, bottom=0.27, wspace=0.48)
    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"omega_window_sensitivity.{suffix}", dpi=450, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

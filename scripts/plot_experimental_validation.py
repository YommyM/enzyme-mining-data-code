#!/usr/bin/env python3
"""Create the manuscript validation figure from candidate-level reaction records."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental_validation"
OUT = ROOT / "figures" / "experimental_validation"


def rows(name):
    with (DATA / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def ec_components(label):
    prefix = "2.6.1."
    values = label.removeprefix(prefix).split("/")
    return {prefix + value for value in values}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lal = rows("lal_products.csv")
    ta = rows("omega_ta_summary.csv")
    reactions = rows("omega_ta_reactions.csv")

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8,
        "axes.labelsize": 8, "axes.titlesize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig, axes = plt.subplots(
        1, 2, figsize=(7.1, 3.65), constrained_layout=True,
        gridspec_kw={"width_ratios": [0.9, 1.75]},
    )

    # Panel a: product breadth of the three experimentally chosen Lal candidates.
    ax = axes[0]
    colors = {"GHU": "#477A9F", "MBD": "#B85C4A", "MCQ": "#D89555"}
    ordered = sorted(lal, key=lambda row: (row["candidate"], float(row["response"])))
    labels = [f'{row["candidate"]}  {row["product"]}' for row in ordered]
    values = [float(row["response"]) / 1000 for row in ordered]
    y = np.arange(len(ordered))
    ax.barh(
        y, values, color=[colors[row["candidate"]] for row in ordered],
        edgecolor="black", linewidth=0.45,
    )
    ax.set_xscale("log")
    ax.set_xlim(9, 450)
    ax.set_xticks([10, 100, 400])
    ax.set_xticklabels(["10", "100", "400"])
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel(r"MS response ($\times 10^3$; arbitrary units)")
    ax.set_title("A   Lal dipeptide products", loc="left", fontweight="bold")

    # Panel b: local-structure EC hypotheses versus candidate-specific substrate tests.
    ax = axes[1]
    candidates = [row["candidate"] for row in ta]
    ec_order = [
        "2.6.1.13", "2.6.1.11", "2.6.1.22", "2.6.1.19",
        "2.6.1.18", "2.6.1.40", "2.6.1.44", "2.6.1.96",
    ]
    reaction_lookup = {}
    for row in reactions:
        ec = row["tested_ec_class"].replace("-associated", "")
        reaction_lookup[(row["candidate"], ec)] = row["result"]

    relation_color = {
        "concordant": "#6F8FA6",
        "discordant": "#D89555",
        "annotation_insufficient": "#A9A9A9",
    }
    relation_text = {"concordant": "C", "discordant": "D", "annotation_insufficient": "U"}
    detected = "#B85C4A"
    not_detected = "#FFFFFF"
    not_tested = "#E1E1E1"
    predicted_edge = "#315D7A"

    for yi, row in enumerate(ta):
        relation = row["uniprot_relation"]
        ax.add_patch(Rectangle((-1.5, yi - 0.5), 1, 1,
                               facecolor=relation_color[relation], edgecolor="white", linewidth=0.8))
        ax.text(-1.0, yi, relation_text[relation], ha="center", va="center",
                color="white" if relation != "annotation_insufficient" else "black",
                fontsize=6.5, fontweight="bold")
        predicted = ec_components(row["inferred_ec"])
        for xi, ec in enumerate(ec_order):
            result = reaction_lookup.get((row["candidate"], ec))
            face = detected if result == "detected" else not_detected if result == "not_detected" else not_tested
            is_predicted = ec in predicted
            ax.add_patch(Rectangle(
                (xi - 0.5, yi - 0.5), 1, 1, facecolor=face,
                edgecolor=predicted_edge if is_predicted else "white",
                linewidth=1.6 if is_predicted else 0.8,
            ))
            if result == "detected":
                ax.text(xi, yi, "+", ha="center", va="center", color="white",
                        fontsize=8, fontweight="bold")
            elif result == "not_detected":
                ax.text(xi, yi, "−", ha="center", va="center", color="black",
                        fontsize=8, fontweight="bold")

    inferred = [row["inferred_ec"].removeprefix("2.6.1.") for row in ta]
    ax.set_yticks(np.arange(len(candidates)))
    ax.set_yticklabels([f"{candidate}   {label}" for candidate, label in zip(candidates, inferred)])
    for label, candidate in zip(ax.get_yticklabels(), candidates):
        if candidate in {"TA7", "TA8", "TA9", "TA10"}:
            label.set_fontweight("bold")
            label.set_color("#8F493D")
    ax.invert_yaxis()
    ax.set_xticks([-1] + list(range(len(ec_order))))
    ax.set_xticklabels(
        ["UniProt\nrelation"] + [ec.removeprefix("2.6.1.") for ec in ec_order],
        rotation=0,
    )
    ax.xaxis.tick_bottom()
    ax.xaxis.set_label_position("bottom")
    ax.set_xlabel("Tested EC-associated substrate class")
    ax.set_xlim(-1.55, len(ec_order) - 0.45)
    ax.set_ylim(len(candidates) - 0.5, -0.5)
    ax.set_title("B   Nearest-local-neighbor hypotheses and reaction tests",
                 loc="left", fontweight="bold")
    ax.tick_params(axis="both", direction="in", top=True, right=True, width=0.8, length=3)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

    legend = [
        Patch(facecolor=detected, edgecolor="black", linewidth=0.4, label="Detected"),
        Patch(facecolor=not_detected, edgecolor="black", linewidth=0.4, label="Not detected"),
        Patch(facecolor=not_tested, edgecolor="white", label="Not tested"),
        Patch(facecolor=not_tested, edgecolor=predicted_edge, linewidth=1.6, label="Predicted EC component"),
        Patch(facecolor=relation_color["concordant"], label="C: concordant"),
        Patch(facecolor=relation_color["discordant"], label="D: discordant"),
        Patch(facecolor=relation_color["annotation_insufficient"], label="U: annotation insufficient"),
    ]
    ax.legend(handles=legend, frameon=False, fontsize=5.5, ncol=4,
              loc="upper center", bbox_to_anchor=(0.50, -0.16),
              columnspacing=0.75, handlelength=1.05, handletextpad=0.35)

    axes[0].tick_params(direction="in", top=True, right=True, width=0.8, length=3)
    for spine in axes[0].spines.values():
        spine.set_linewidth(0.8)
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"experimental_validation_summary.{suffix}", dpi=400, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

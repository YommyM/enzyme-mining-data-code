#!/usr/bin/env python3
"""Plot Lal NR retrieval, whole-protein structural similarity, and selected candidates."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "lals_reference"
OUT = ROOT / "figures" / "lal_uniprot_branch"


def matrix(name):
    rows = list(csv.reader((DATA / name).open()))
    return rows[0][1:], np.array([[float(value) for value in row[1:]] for row in rows[1:]])


def add_box(ax, xy, width, height, text, color, fontsize=6.2):
    ax.add_patch(plt.Rectangle(xy, width, height, facecolor=color,
                               edgecolor="#252525", linewidth=0.65))
    ax.text(xy[0] + width / 2, xy[1] + height / 2, text,
            ha="center", va="center", color="white", fontsize=fontsize,
            fontweight="bold", linespacing=1.18)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ids, global_tm = matrix("global_tm_similarity_nr_103.csv")
    metadata = {row["ID"]: row for row in csv.DictReader(
        (DATA / "nr_structure_set_103.csv").open())}
    if set(ids) != set(metadata):
        raise ValueError("Lal matrix and metadata contain different proteins")

    distance = np.clip(1 - global_tm, 0, 1)
    np.fill_diagonal(distance, 0)
    order = leaves_list(linkage(squareform(distance, checks=False), method="average"))
    roles = np.array([[0 if metadata[ids[i]]["role"] == "reference" else 1 for i in order]])
    selected_ligase_cluster = {"Q80WS1", "Q8IXN7", "Q9ULI2", "Q6PFX8",
                      "GHU07852", "MBD5089657", "MCQ2799571"}
    functional_group = np.array([[1 if ids[i] in selected_ligase_cluster else 0 for i in order]])

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8,
        "axes.labelsize": 8, "axes.titlesize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig = plt.figure(figsize=(7.1, 3.65))
    outer = fig.add_gridspec(1, 3, width_ratios=[0.88, 1.28, 0.93])

    # a, retrieval logic
    ax = fig.add_subplot(outer[0])
    add_box(ax, (0.12, 2.48), 0.76, 0.48, "29 characterized Lals\n8 sequence clusters", "#477A9F")
    add_box(ax, (0.18, 1.76), 0.64, 0.43, "9 profile HMMs", "#526A70")
    add_box(ax, (0.08, 0.91), 0.84, 0.57, "NR retrieval\n$E<10^{-5}$; score $>50$", "#6E858A")
    add_box(ax, (0.08, 0.02), 0.84, 0.60, "Deduplication and\nPfam-stratified sampling\n74 candidates", "#D89555", 5.9)
    for y0, y1 in [(2.47, 2.20), (1.75, 1.49), (0.90, 0.63)]:
        ax.annotate("", (0.50, y1), (0.50, y0),
                    arrowprops={"arrowstyle": "-|>", "lw": 0.8, "color": "#4A4A4A"})
    ax.set_xlim(0, 1); ax.set_ylim(-0.05, 3.13); ax.set_axis_off()
    ax.set_title("A   Sequence-guided retrieval", loc="left", fontweight="bold")

    # b, whole-protein structural similarity
    sub = outer[1].subgridspec(
        7, 1, height_ratios=[0.095, 0.035, 0.035, 1, 0.065, 0.07, 0.12], hspace=0.045
    )
    title_ax = fig.add_subplot(sub[0]); title_ax.set_axis_off()
    title_ax.set_title("B   Global structure ($n=103$)",
                       loc="left", fontweight="bold", fontsize=8, pad=5)
    strip_role = fig.add_subplot(sub[1])
    strip_role.imshow(roles, aspect="auto", cmap=ListedColormap(["#477A9F", "#D89555"]),
                      vmin=0, vmax=1)
    strip_role.set_axis_off()
    strip_group = fig.add_subplot(sub[2])
    strip_group.imshow(functional_group, aspect="auto",
                       cmap=ListedColormap(["#E8E8E8", "#B13F4A"]), vmin=0, vmax=1)
    strip_group.set_axis_off()
    ax_hm = fig.add_subplot(sub[3])
    image = ax_hm.imshow(global_tm[np.ix_(order, order)], vmin=0, vmax=1,
                         cmap="RdBu", interpolation="nearest", aspect="equal")
    ax_hm.set_xticks([]); ax_hm.set_yticks([])
    for spine in ax_hm.spines.values():
        spine.set_linewidth(0.7)
    colorbar_ax = fig.add_subplot(sub[5])
    colorbar = fig.colorbar(image, cax=colorbar_ax, orientation="horizontal", ticks=[0, 0.5, 1])
    colorbar.ax.set_title("Structural similarity", fontsize=6.2, pad=2)
    colorbar.ax.tick_params(direction="in", length=2, pad=1)
    legend_ax = fig.add_subplot(sub[6]); legend_ax.set_axis_off()
    legend_ax.legend(
        handles=[
            Patch(facecolor="#477A9F", label="Reference"),
            Patch(facecolor="#D89555", label="Candidate"),
            Patch(facecolor="#B13F4A", label="Selected cluster"),
            Patch(facecolor="#E8E8E8", edgecolor="#AAAAAA", label="Other"),
        ],
        loc="center", ncol=2, frameon=False, fontsize=5.7,
        handlelength=1.0, columnspacing=0.8, handletextpad=0.35,
    )

    # c, candidate-to-reference functional link used for experimental selection
    ax = fig.add_subplot(outer[2])
    add_box(ax, (0.03, 1.18), 0.39, 0.68, "Dipeptide-\nforming ligases\n(references)", "#477A9F", 5.1)
    candidate_boxes = [
        ("GHU07852\nbroad automatic\nannotation", 2.24),
        ("MBD5089657\nbroad automatic\nannotation", 1.20),
        ("MCQ2799571\nbroad automatic\nannotation", 0.16),
    ]
    for text_label, y in candidate_boxes:
        add_box(ax, (0.52, y), 0.46, 0.68, text_label, "#D89555", 5.1)
        ax.annotate("", (0.51, y + 0.34), (0.43, 1.52),
                    arrowprops={"arrowstyle": "-|>", "lw": 0.65, "color": "#6A6A6A",
                                "connectionstyle": "arc3,rad=0.0"})
    ax.text(0.50, -0.02, "Three candidates selected\nfor experimental testing",
            ha="center", va="top", fontsize=5.2, color="#333333")
    ax.set_xlim(0, 1); ax.set_ylim(-0.18, 3.13); ax.set_axis_off()
    ax.set_title("C   Candidate selection", loc="left", fontweight="bold")

    fig.subplots_adjust(left=0.025, right=0.985, bottom=0.08, top=0.91, wspace=0.22)
    fig.canvas.draw()
    heat_position = ax_hm.get_position()
    for strip in (strip_role, strip_group):
        strip_position = strip.get_position()
        strip.set_position([heat_position.x0, strip_position.y0, heat_position.width, strip_position.height])
    cb_position = colorbar_ax.get_position()
    colorbar_ax.set_position([heat_position.x0, 0.105, heat_position.width, 0.020])
    legend_ax.set_position([heat_position.x0, 0.010, heat_position.width, 0.055])
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"lal_uniprot_branch.{suffix}", dpi=400)
    plt.close(fig)


if __name__ == "__main__":
    main()

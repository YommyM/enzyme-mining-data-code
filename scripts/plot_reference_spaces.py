#!/usr/bin/env python3
"""Supporting Figure S1: sequence organization used to build profile HMMs."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import numpy as np
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "reference_spaces"
CMAP = LinearSegmentedColormap.from_list(
    "sequence_similarity", ["#FAFAFA", "#DCE8EF", "#9DBFD1", "#4E87A6", "#174A68"]
)
GROUP_COLORS = ["#4A4A4A", "#477A9F", "#B85C4A", "#D89555",
                "#6C8A63", "#8A6E9B", "#8B7355", "#6E8B8B"]


def read_matrix(path):
    rows = list(csv.reader(path.open()))
    return rows[0][1:], np.array([[float(value) for value in row[1:]] for row in rows[1:]])


def grouped_order(ids, groups, values):
    order = []
    for group in sorted(set(groups.values()), key=lambda x: int(x.replace("HMM", ""))):
        positions = [i for i, protein_id in enumerate(ids) if groups[protein_id] == group]
        if len(positions) > 2:
            distance = np.clip(1 - values[np.ix_(positions, positions)], 0, 1)
            np.fill_diagonal(distance, 0)
            within = leaves_list(linkage(squareform(distance, checks=False), method="average"))
            positions = [positions[i] for i in within]
        order.extend(positions)
    return order


def draw(fig, spec, values, ids, groups, title, panel, vmax_groups):
    order = grouped_order(ids, groups, values)
    ordered_ids = [ids[i] for i in order]
    grid = spec.subgridspec(2, 1, height_ratios=[0.055, 1], hspace=0.025)
    strip = fig.add_subplot(grid[0])
    group_numbers = np.array([[int(groups[x].replace("HMM", "")) - 1 for x in ordered_ids]])
    strip.imshow(group_numbers, aspect="auto", interpolation="nearest",
                 cmap=ListedColormap(GROUP_COLORS[:vmax_groups]), vmin=0,
                 vmax=vmax_groups - 1)
    strip.set_axis_off()
    strip.set_title(f"{panel}   {title}", loc="left", fontweight="bold", pad=5)

    ax = fig.add_subplot(grid[1])
    image = ax.imshow(values[np.ix_(order, order)], vmin=0, vmax=1, cmap=CMAP,
                      interpolation="nearest", aspect="equal", rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    ordered_groups = [groups[x] for x in ordered_ids]
    boundaries = [i - 0.5 for i in range(1, len(ordered_groups))
                  if ordered_groups[i] != ordered_groups[i - 1]]
    for boundary in boundaries:
        ax.axhline(boundary, color="white", lw=0.45)
        ax.axvline(boundary, color="white", lw=0.45)
    for spine in ax.spines.values():
        spine.set_linewidth(0.65)
    return image, ax, strip


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lals = ROOT / "data" / "lals_reference"
    omega = ROOT / "data" / "omega_reference"
    lal_ids, lal_seq = read_matrix(lals / "reference_sequence_identity_29.csv")
    lal_groups = {}
    build = ROOT / "data" / "lals_reference" / "profile_hmms"
    for number in range(1, 9):
        for line in (build / f"class{number}.fasta").read_text().splitlines():
            if line.startswith(">"):
                lal_groups[line[1:].split()[0]] = f"HMM{number}"

    metadata = list(csv.DictReader((omega / "reference_metadata_58.csv").open()))
    omega_groups = {row["ID"]: row["HMM_group"] for row in metadata}
    omega_ids, omega_seq = read_matrix(omega / "sequence_similarity_58.csv")

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7.0, "axes.titlesize": 8.2,
        "xtick.labelsize": 6.2, "ytick.labelsize": 6.2,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig = plt.figure(figsize=(7.1, 3.55))
    outer = fig.add_gridspec(1, 2, wspace=0.22)
    panels = [
        draw(fig, outer[0], lal_seq, lal_ids, lal_groups,
             "Lal sequence identity", "A", 8),
        draw(fig, outer[1], omega_seq, omega_ids, omega_groups,
             r"$\omega$-TA sequence similarity", "B", 6),
    ]
    fig.subplots_adjust(left=0.06, right=0.985, top=0.90, bottom=0.18)
    fig.canvas.draw()
    for _, heat_ax, strip_ax in panels:
        heat = heat_ax.get_position(); strip = strip_ax.get_position()
        strip_ax.set_position([heat.x0, heat.y1 + 0.008, heat.width, strip.height])
    left = panels[0][1].get_position().x0
    right = panels[1][1].get_position().x1
    bottom = min(panel[1].get_position().y0 for panel in panels)
    cax = fig.add_axes([left, bottom - 0.10, right - left, 0.018])
    cbar = fig.colorbar(panels[-1][0], cax=cax, orientation="horizontal",
                        ticks=[0, 0.25, 0.5, 0.75, 1])
    cbar.set_label("Sequence similarity", labelpad=1)
    cbar.ax.tick_params(direction="in", length=2, pad=1)
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"reference_spaces.{suffix}", dpi=450, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

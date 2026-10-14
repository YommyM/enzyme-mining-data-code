#!/usr/bin/env python3
"""Figure 3: global/local omega-TA matrices and functional benchmark."""

from pathlib import Path
import csv
import itertools
import re
import warnings

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, Normalize
from matplotlib.patches import Patch
import numpy as np
from scipy.cluster.hierarchy import ClusterWarning, dendrogram, leaves_list, linkage
from scipy.spatial.distance import cdist


ROOT = Path(__file__).resolve().parents[1]
ARTICLE = ROOT
DATA = ROOT / "data" / "omega_reference"
OUT = ROOT / "figures" / "omega_local_advantage"
HEATMAP = LinearSegmentedColormap.from_list(
    "structural_similarity", ["#FAFAFA", "#E6EDF2", "#A9C4D3", "#5E91AD", "#245A78", "#102F45"]
)
GLOBAL_COLOR = "#7A7A7A"
LOCAL_COLOR = "#0072B2"
SEED = 20260923


def original_matrix(kind, expected_ids):
    """Read the frozen manuscript matrix instead of thousands of raw outputs."""
    filename = {"global": "global_tm_original_workflow_58.csv",
                "local": "local_tm_original_workflow_window_3_58.csv"}[kind]
    rows = list(csv.reader((DATA / filename).open()))
    if rows[0][1:] != expected_ids or [r[0] for r in rows[1:]] != expected_ids:
        raise ValueError(f"Identifier order mismatch: {filename}")
    matrix = np.asarray([[float(v) for v in r[1:]] for r in rows[1:]])
    return matrix


def write_matrix(path, ids, matrix):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ID", *ids])
        for identifier, row in zip(ids, matrix):
            writer.writerow([identifier, *[f"{value:.5f}" for value in row]])


def hierarchy(matrix):
    # This is algebraically the same operation used in the original notebook:
    # linkage(1 - tm_df.values, method="average").
    with warnings.catch_warnings():
        # SciPy warns because this square observation matrix resembles a
        # precomputed distance matrix; treating rows as observations is the
        # intentional behavior of the original notebook.
        warnings.simplefilter("ignore", ClusterWarning)
        tree = linkage(np.clip(1.0 - matrix, 0.0, 1.0), method="average")
    return tree, leaves_list(tree)


def functional_memberships(curation_rows, level_rows):
    """Return evidence-based, potentially nonexclusive functional labels."""
    memberships = {}
    terminal = "terminal-amino/diamino--2-oxoglutarate aminotransferase"
    gaba = "GABA/5-aminovalerate--2-oxoglutarate aminotransferase"
    selective = "S-selective aromatic/aliphatic amine transaminase (Fold I)"
    broad = "broad-spectrum amine transaminase (Fold I)"
    for row in curation_rows:
        identifier = row["structure_id"].lower()
        labels = [level_rows[identifier]["substrate_class"]]
        # The S8EY38 parent is experimentally annotated for both ornithine and
        # GABA transamination.  All three deposited structures retain both.
        if "ornithine/GABA" in row["protein_function"]:
            labels.extend([terminal, gaba])
        # PDB 5VIU is assigned ArgD, while its mapped UniProt record also
        # carries the OAT-associated EC 2.6.1.13 annotation.
        if identifier == "5viu":
            labels.append(terminal)
        # S selectivity and broad amine scope overlap rather than form mutually
        # exclusive biochemical classes for these characterized enzymes.
        if row["stereochemical_scope"] == "S-selective" and (
            "broad" in row["detailed_biochemical_subgroup"]
            or "broad" in row["protein_function"]
        ):
            labels.extend([selective, broad])
        memberships[identifier] = tuple(dict.fromkeys(labels))
    return memberships


def draw_matrix(fig, spec, matrix, memberships, fold_labels, label_to_index,
                function_colors, title, panel):
    tree, order = hierarchy(matrix)
    grid = spec.subgridspec(
        3, 3, height_ratios=[0.105, 0.22, 1], width_ratios=[0.060, 0.030, 1],
        hspace=0.012, wspace=0.022,
    )
    title_ax = fig.add_subplot(grid[0, :])
    title_ax.set_axis_off()
    title_ax.set_title(f"{panel}   {title}", loc="left", fontweight="bold", fontsize=8.2, pad=4)

    dendro_ax = fig.add_subplot(grid[1, 2])
    dendrogram(tree, ax=dendro_ax, no_labels=True, color_threshold=0,
               above_threshold_color="#555555", link_color_func=lambda _: "#555555")
    dendro_ax.set_axis_off()

    ordered_function = []
    for i in order:
        labels = list(memberships[i])[:2]
        if len(labels) == 1:
            labels *= 2
        ordered_function.append([label_to_index[labels[0]], label_to_index[labels[1]]])
    function_ax = fig.add_subplot(grid[2, 0])
    function_ax.imshow(np.asarray(ordered_function), aspect="auto", interpolation="nearest",
                       cmap=ListedColormap(function_colors), vmin=0,
                       vmax=len(function_colors) - 1)
    function_ax.set_xticks([0.5])
    function_ax.set_xticklabels(["Substrate"], rotation=90, fontsize=5.0)
    function_ax.set_yticks([])

    ordered_fold = np.array([[fold_labels[i] - 1] for i in order])
    fold_ax = fig.add_subplot(grid[2, 1])
    fold_ax.imshow(ordered_fold, aspect="auto", interpolation="nearest",
                   cmap=ListedColormap(["#B8BDC5", "#343A40"]), vmin=0, vmax=1)
    fold_ax.set_xticks([0])
    fold_ax.set_xticklabels(["Fold"], rotation=90, fontsize=5.2)
    fold_ax.set_yticks([])
    for annotation_ax in (function_ax, fold_ax):
        annotation_ax.tick_params(length=0, pad=1, labeltop=True, labelbottom=False)
        for spine in annotation_ax.spines.values():
            spine.set_visible(False)

    heat_ax = fig.add_subplot(grid[2, 2])
    image = heat_ax.imshow(matrix[np.ix_(order, order)], norm=Normalize(vmin=0.2, vmax=1.0),
                           cmap=HEATMAP, interpolation="nearest", rasterized=True,
                           aspect="equal")
    heat_ax.set_xticks([])
    heat_ax.set_yticks([])
    for spine in heat_ax.spines.values():
        spine.set_linewidth(0.65)
    return image, heat_ax, (function_ax, fold_ax)


def profile_benchmark(matrix, ids, rows, memberships, primary):
    roles = {row["structure_id"].lower(): row["benchmark_role"] for row in rows}
    parents = {row["structure_id"].lower(): row["parent_entity"] for row in rows}
    candidates = [i for i, identifier in enumerate(ids)
                  if roles[identifier.lower()] in {"query_and_candidate", "candidate_only"}]
    queries = [i for i, identifier in enumerate(ids)
               if roles[identifier.lower()] == "query_and_candidate"]
    distance = cdist(1.0 - matrix, 1.0 - matrix, metric="euclidean")
    records = []
    for query in queries:
        qid = ids[query].lower()
        available = [j for j in candidates if j != query and
                     parents[ids[j].lower()] != parents[qid]]
        ranked = sorted(available, key=lambda j: (distance[query, j], ids[j].lower()))
        relevant = np.array([
            bool(set(memberships[qid]) & set(memberships[ids[j].lower()]))
            for j in ranked
        ], dtype=int)
        total_relevant = int(relevant.sum())
        record = {"structure_id": ids[query], "family": primary[qid]}
        for k in range(1, 6):
            denominator = min(k, total_relevant)
            record[f"recovery_{k}"] = float(relevant[:k].sum() / denominator) if denominator else np.nan
        precision = np.cumsum(relevant) / np.arange(1, len(relevant) + 1)
        record["average_precision"] = float((precision * relevant).sum() / total_relevant)
        records.append(record)
    return records


def family_summary(records, key):
    families = list(dict.fromkeys(row["family"] for row in records))
    return {
        family: float(np.mean([row[key] for row in records if row["family"] == family]))
        for family in families
    }


def bootstrap_ci(values, draws=20000):
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(SEED)
    means = np.mean(rng.choice(values, size=(draws, len(values)), replace=True), axis=1)
    return tuple(np.quantile(means, [0.025, 0.975]))


def exact_sign_flip(differences):
    differences = np.asarray(differences, dtype=float)
    observed = abs(float(np.mean(differences)))
    null = [abs(float(np.mean(differences * np.asarray(signs))))
            for signs in itertools.product((-1, 1), repeat=len(differences))]
    return float(np.mean(np.asarray(null) >= observed - 1e-12))


def draw_quantification(fig, spec_c, spec_d, global_records, local_records,
                        range_rows=None):
    families = list(family_summary(global_records, "average_precision"))
    ks = np.arange(1, 6)
    series = {}
    for name, records in (("Global", global_records), ("Local ±3", local_records)):
        means, lows, highs = [], [], []
        for k in ks:
            summary = family_summary(records, f"recovery_{k}")
            values = [summary[family] for family in families]
            low, high = bootstrap_ci(values)
            means.append(np.mean(values)); lows.append(low); highs.append(high)
        series[name] = (np.asarray(means), np.asarray(lows), np.asarray(highs))

    ax_c = fig.add_subplot(spec_c)
    for name, color, marker in (("Global", GLOBAL_COLOR, "o"), ("Local ±3", LOCAL_COLOR, "s")):
        mean, low, high = series[name]
        ax_c.plot(ks, mean, color=color, marker=marker, lw=1.5, ms=3.8, label=name, zorder=3)
        ax_c.fill_between(ks, low, high, color=color, alpha=0.13, lw=0)
    g3 = family_summary(global_records, "recovery_3")
    l3 = family_summary(local_records, "recovery_3")
    p3 = exact_sign_flip([l3[f] - g3[f] for f in families])
    ax_c.text(0.98, 0.04, f"paired-group sign-flip at k = 3\n$p$ = {p3:.3f}",
              transform=ax_c.transAxes, ha="right", va="bottom", fontsize=5.8)
    ax_c.set_title("C   Compatible substrate-group neighbors", loc="left", fontweight="bold", fontsize=8.2)
    ax_c.set_xlabel("Number of profile neighbors, k")
    ax_c.set_ylabel("Group-balanced recovery")
    ax_c.set_xticks(ks)
    ax_c.set_ylim(0.25, 0.84)
    ax_c.grid(axis="y", color="#D9D9D9", lw=0.45)
    ax_c.spines[["top", "right"]].set_visible(False)
    ax_c.legend(frameon=False, loc="upper right", fontsize=6.2, handlelength=1.6)

    ax_d = fig.add_subplot(spec_d)
    global_ap = family_summary(global_records, "average_precision")
    local_ap = family_summary(local_records, "average_precision")
    if range_rows is not None:
        windows = [3, 5, 7, 10, 15, 20]
        by_window = {int(row["representation"].replace("window_", "")): row
                     for row in range_rows}
        x = np.arange(1, len(windows) + 1)
        local_means = np.array([float(by_window[w]["average_precision"]) for w in windows])
        local_low = np.array([float(by_window[w]["average_precision_ci_low"]) for w in windows])
        local_high = np.array([float(by_window[w]["average_precision_ci_high"]) for w in windows])
        global_mean = float(np.mean(list(global_ap.values())))
        global_low, global_high = bootstrap_ci(list(global_ap.values()))

        local_error = np.vstack([local_means - local_low, local_high - local_means])
        ax_d.errorbar(x, local_means, yerr=local_error, color=LOCAL_COLOR, lw=1.45,
                      marker="o", ms=3.8, mec="white", mew=0.45, capsize=2.0,
                      elinewidth=0.8, label="Local windows", zorder=3)
        ax_d.errorbar([0], [global_mean],
                      yerr=[[global_mean - global_low], [global_high - global_mean]],
                      fmt="D", ms=4.2, color=GLOBAL_COLOR, mfc="white", mec=GLOBAL_COLOR,
                      capsize=2.4, lw=1.0, label="Whole protein", zorder=4)
        significant = [w for w in windows if float(by_window[w]["p_average_precision"]) < 0.05]
        for xi, w in zip(x, windows):
            if w in significant:
                ax_d.text(xi, local_means[xi - 1] + 0.018, "*", color=LOCAL_COLOR,
                          ha="center", va="bottom", fontsize=8.0, fontweight="bold",
                          clip_on=True)
        ax_d.set_title("D   Dependence on active-site window", loc="left",
                       fontweight="bold", fontsize=8.2)
        ax_d.set_ylabel("Family-balanced average precision")
        ax_d.set_xticks(np.arange(0, len(windows) + 1))
        ax_d.set_xticklabels(["Global", *[rf"$\pm${w}" for w in windows]], rotation=35,
                             ha="right", rotation_mode="anchor")
        ax_d.set_xlim(-0.45, len(windows) + 0.35)
        ax_d.set_ylim(0.34, 0.80)
        ax_d.grid(axis="y", color="#E4E4E4", lw=0.45)
        ax_d.spines[["top", "right"]].set_visible(False)
        ax_d.legend(frameon=False, loc="upper right", fontsize=5.8, handlelength=1.4)
        return p3, exact_sign_flip([local_ap[f] - global_ap[f] for f in families])

    short = {
        "N-acyl diamino-acid aminotransferase": "N-acyl diamino",
        "GABA/5-aminovalerate--2-oxoglutarate aminotransferase": "GABA/5AV",
        "GABA--pyruvate aminotransferase": "GABA–Pyr",
        "terminal-amino/diamino--2-oxoglutarate aminotransferase": "terminal/diamine",
        "beta-amino-acid aminotransferase": "β-amino acid",
        "noncarboxylate acidic omega-amino aminotransferase": "P/S ω-amino",
        "S-selective aromatic/aliphatic amine transaminase (Fold I)": "S-selective",
        "broad-spectrum amine transaminase (Fold I)": "broad amine",
    }
    order = sorted(families, key=lambda family: local_ap[family] - global_ap[family])
    y = np.arange(len(order))
    for position, family in enumerate(order):
        ax_d.plot([global_ap[family], local_ap[family]], [position, position],
                  color="#C2C2C2", lw=0.85, zorder=1)
    ax_d.scatter([global_ap[f] for f in order], y, color=GLOBAL_COLOR, s=17,
                 edgecolor="white", linewidth=0.35, label="Global", zorder=2)
    ax_d.scatter([local_ap[f] for f in order], y, color=LOCAL_COLOR, s=17, marker="s",
                 edgecolor="white", linewidth=0.35, label="Local ±3", zorder=3)
    p_ap = exact_sign_flip([local_ap[f] - global_ap[f] for f in families])
    ax_d.set_title("D   Ranking by substrate group", loc="left",
                   fontweight="bold", fontsize=8.2)
    ax_d.set_xlabel("Average precision")
    ax_d.set_yticks(y)
    ax_d.set_yticklabels([short.get(f, f) for f in order], fontsize=4.5)
    ax_d.tick_params(axis="y", pad=1)
    ax_d.set_xlim(0, 1.02)
    ax_d.grid(axis="x", color="#E2E2E2", lw=0.45)
    ax_d.spines[["top", "right"]].set_visible(False)
    legend = ax_d.legend(frameon=False, loc="upper left", fontsize=5.4,
                         handlelength=1.0, title=rf"paired $p$ = {p_ap:.3f}")
    legend.get_title().set_fontsize(5.2)
    return p3, p_ap


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader((DATA / "reference_function_reclassification.csv").open()))
    level_rows = {row["structure_id"].lower(): row
                  for row in csv.DictReader((DATA / "reference_function_levels.csv").open())}
    ids = [row["structure_id"] for row in rows]
    primary = {identifier.lower(): level_rows[identifier.lower()]["substrate_class"] for identifier in ids}
    memberships_by_id = functional_memberships(rows, level_rows)
    memberships = [memberships_by_id[identifier.lower()] for identifier in ids]

    global_tm = original_matrix("global", ids)
    local_tm = original_matrix("local", ids)
    write_matrix(OUT / "global_tm_original_workflow_58.csv", ids, global_tm)
    write_matrix(OUT / "local_tm_original_workflow_window_3_58.csv", ids, local_tm)

    benchmark_global = profile_benchmark(global_tm, ids, rows, memberships_by_id, primary)
    benchmark_local = profile_benchmark(local_tm, ids, rows, memberships_by_id, primary)
    benchmark_path = OUT / "original_workflow_multilabel_profile_benchmark.csv"
    with benchmark_path.open("w", newline="") as handle:
        fields = ["representation", "structure_id", "family", *[f"recovery_{k}" for k in range(1, 6)],
                  "average_precision"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for representation, records in (("global", benchmark_global), ("local_window_3", benchmark_local)):
            for record in records:
                writer.writerow({"representation": representation, **record})

    display_classes = list(dict.fromkeys(label for labels in memberships for label in labels))
    label_to_index = {label: index for index, label in enumerate(display_classes)}
    # Use one canonical left-to-right order for every multifunctional strip;
    # otherwise the same pair can be reversed when a different annotation was
    # encountered first for a given protein.
    memberships = [
        tuple(sorted(labels, key=lambda label: label_to_index[label]))
        for labels in memberships
    ]
    palette = [
        "#0072B2", "#D55E00", "#009E73", "#E69F00", "#6A3D9A", "#56B4E9",
        "#CC79A7", "#8C510A", "#1B9E77", "#E7298A", "#66A61E", "#A6761D",
        "#7570B3", "#E6AB02", "#666666", "#222222",
    ]
    function_colors = palette[:len(display_classes)]
    fold4 = {line[1:].split()[0].lower()
             for line in (DATA / "foldtype_4.fasta").read_text().splitlines()
             if line.startswith(">")}
    fold_labels = [2 if identifier.lower() in fold4 else 1 for identifier in ids]

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7.0, "axes.labelsize": 7.0,
        "axes.titlesize": 8.0, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig = plt.figure(figsize=(7.15, 5.15))
    outer = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.02], wspace=0.34)
    quantitative = outer[0, 2].subgridspec(2, 1, hspace=0.46)
    panels = [
        draw_matrix(fig, outer[0, 0], global_tm, memberships, fold_labels, label_to_index,
                    function_colors, "Global structure", "A"),
        draw_matrix(fig, outer[0, 1], local_tm, memberships, fold_labels, label_to_index,
                    function_colors, "Local active-site structure (±3)", "B"),
    ]
    p3, p_ap = draw_quantification(fig, quantitative[0], quantitative[1],
                                    benchmark_global, benchmark_local)

    legend_names = {
        "GABA/5-aminovalerate--2-oxoglutarate aminotransferase": "GABA/5-aminovalerate–2OG TA",
        "S-selective aromatic/aliphatic amine transaminase (Fold I)": "S-selective aromatic/aliphatic amine TA",
        "broad-spectrum amine transaminase (Fold I)": "broad-spectrum amine TA",
        "noncarboxylate acidic omega-amino aminotransferase": "phosphonate/sulfonate ω-amino TA",
        "terminal-amino/diamino--2-oxoglutarate aminotransferase": "terminal-amino/diamine–2OG TA",
    }
    handles = [Patch(facecolor=function_colors[i], edgecolor="none",
                     label=legend_names.get(label, label))
               for i, label in enumerate(display_classes)]
    handles += [Patch(facecolor="#B8BDC5", label="Fold I"),
                Patch(facecolor="#343A40", label="Fold IV")]
    fig.subplots_adjust(left=0.055, right=0.99, top=0.955, bottom=0.30)

    fig.canvas.draw()
    for _, heat_ax, annotation_axes in panels:
        heat_position = heat_ax.get_position()
        for annotation_ax in annotation_axes:
            strip_position = annotation_ax.get_position()
            annotation_ax.set_position([strip_position.x0, heat_position.y0,
                                        strip_position.width, heat_position.height])

    # Place the common structural-similarity scale and the function legend as
    # a compact block directly beneath the two matrix panels.
    left = panels[0][1].get_position().x0
    right = panels[1][1].get_position().x1
    bottom = min(panel[1].get_position().y0 for panel in panels)
    cbar_ax = fig.add_axes([left, bottom - 0.060, right - left, 0.011])
    cbar = fig.colorbar(panels[-1][0], cax=cbar_ax, orientation="horizontal",
                        ticks=[0.2, 0.4, 0.6, 0.8, 1.0])
    cbar.set_label("Structural similarity", labelpad=1)
    cbar.ax.tick_params(direction="in", length=2, pad=1)
    fig.legend(handles=handles, loc="upper center",
               bbox_to_anchor=((left + right) / 2, bottom - 0.096), ncol=4,
               frameon=False, fontsize=4.35, handlelength=0.78, columnspacing=0.62,
               handletextpad=0.28)

    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"omega_local_advantage.{suffix}", dpi=450, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure 3 written to {OUT}")
    print(f"k=3 exact family sign-flip p={p3:.6f}; AP p={p_ap:.6f}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Replot TA chromatograms from vector paths embedded in the source PPTX.

The EMF objects were converted to SVG without tracing a raster image.  Blue signal
paths are read directly, mapped to the original 0--65 min time axis, normalized
within each trace, and drawn with common axes and line widths.  Normalization makes
the result suitable for qualitative peak/retention-time evidence, not quantitation.
"""

from pathlib import Path
import re
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
SVG = ROOT / "data" / "experimental_validation" / "extracted_plots" / "ta" / "svg"
OUT = ROOT / "figures" / "experimental_validation"

BLACK = "#262626"
GREY = "#777777"
LIGHT = "#d9d9d9"
BLUE = "#587c95"
RED = "#a85f59"
ROW_SPACING = 1.30
TRACE_HEIGHT = 0.62

# Tick-label centers recovered from the vector export: 10 min = 131.57 pt and
# 20 min = 209.86 pt.  The relationship is shared by the supplied TA panels.
X_ZERO = 53.28
X_PER_MIN = 7.829

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7.5,
    "axes.titlesize": 7.5,
    "axes.labelsize": 7.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def _points(path_d):
    return [(float(x), float(y)) for x, y in
            re.findall(r"[ML]\s*([-0-9.]+)\s+([-0-9.]+)", path_d)]


def extract_traces(path):
    """Return bottom-to-top source traces as lists of vector polyline segments."""
    root = ET.parse(path).getroot()
    candidates = []
    for element in root.iter():
        style = element.attrib.get("style", "")
        if "stroke:rgb(0%,0%,100%)" not in style:
            continue
        pts = _points(element.attrib.get("d", ""))
        if len(pts) < 10:
            continue
        xs = [p[0] for p in pts]
        if max(xs) - min(xs) < 0.5:
            continue
        median_y = sorted(p[1] for p in pts)[len(pts) // 2]
        candidates.append((median_y, pts))

    # Instrument exports contain either four evenly spaced traces (whole-cell)
    # or three traces (purified enzyme). Large gaps robustly identify the rows.
    centers = []
    for median_y, _ in sorted(candidates):
        if not centers or median_y - centers[-1][-1] > 25:
            centers.append([median_y])
        else:
            centers[-1].append(median_y)
    row_centers = [sum(group) / len(group) for group in centers]
    traces = [[] for _ in row_centers]
    for median_y, pts in candidates:
        index = min(range(len(row_centers)), key=lambda i: abs(median_y - row_centers[i]))
        traces[index].append(pts)

    output = []
    for segments in traces:
        baseline = min(y for segment in segments for _, y in segment)
        amplitude = max(y - baseline for segment in segments for _, y in segment)
        amplitude = amplitude or 1.0
        converted = []
        for segment in segments:
            xy = [((x - X_ZERO) / X_PER_MIN, (y - baseline) / amplitude)
                  for x, y in segment]
            xy = [(x, y) for x, y in xy if 0 <= x <= 65]
            if len(xy) >= 2:
                converted.append(xy)
        output.append(converted)
    # SVG source coordinates run opposite to the visual top-to-bottom order.
    return list(reversed(output))


def extract_slide10_ta4_traces(path):
    """Recover the native-vector TA4/8a panel embedded directly in slide 10.

    Unlike the neighboring TA3 panel, this chromatogram is stored as PowerPoint
    drawing paths rather than as a separate EMF media object.  Coordinates below
    are calibrated from the product and donor standards in the rendered vector
    slide (L-glutamate, 34.3 min; donor 8a, 46.9 min).
    """
    root = ET.parse(path).getroot()
    row_baselines = [286.7, 340.9, 396.2, 451.5]
    raw = [[] for _ in row_baselines]
    matrix_pattern = re.compile(
        r"matrix\(([-0-9.]+),([-0-9.]+),([-0-9.]+),"
        r"([-0-9.]+),([-0-9.]+),([-0-9.]+)\)"
    )
    for element in root.iter():
        if "stroke:rgb(0%,0%,100%)" not in element.attrib.get("style", ""):
            continue
        pts = _points(element.attrib.get("d", ""))
        if len(pts) < 2:
            continue
        match = matrix_pattern.search(element.attrib.get("transform", ""))
        if match:
            a, b, c, d, tx, ty = map(float, match.groups())
            pts = [(a * x + c * y + tx, b * x + d * y + ty) for x, y in pts]
        pts = [(x, y) for x, y in pts if 494.8 <= x <= 895.2]
        if len(pts) < 2 or max(x for x, _ in pts) - min(x for x, _ in pts) < 0.2:
            continue
        # Peaks point upward in rendered-slide coordinates, so max(y) estimates
        # the local baseline for each segmented path.
        path_baseline = max(y for _, y in pts)
        index = min(range(4), key=lambda i: abs(path_baseline - row_baselines[i]))
        if abs(path_baseline - row_baselines[index]) > 2.2:
            continue
        if min(y for _, y in pts) < row_baselines[index] - 28:
            continue
        raw[index].append(pts)

    output = []
    x_zero = 495.0
    x_per_min = 6.14603
    for segments in raw:
        baseline = max(y for segment in segments for _, y in segment)
        amplitude = baseline - min(y for segment in segments for _, y in segment)
        amplitude = amplitude or 1.0
        converted = []
        for segment in segments:
            xy = [((x - x_zero) / x_per_min, (baseline - y) / amplitude)
                  for x, y in segment]
            xy = [(x, y) for x, y in xy if 0 <= x <= 65]
            if len(xy) >= 2:
                converted.append(xy)
        output.append(converted)
    return output


def trace_height(segments, retention_time, window=0.35):
    nearby = [y for segment in segments for x, y in segment
              if abs(x - retention_time) <= window]
    return max(nearby) if nearby else 0.0


def peak_label(ax, segments, row_offset, retention_time, text, *, color=BLACK,
               dx=0.0, text_height=0.82, align="center"):
    peak_y = row_offset + TRACE_HEIGHT * trace_height(segments, retention_time)
    label_y = row_offset + text_height
    ax.plot([retention_time], [peak_y], marker="o", ms=2.0, color=color,
            markeredgewidth=0, zorder=4)
    ax.annotate(text, xy=(retention_time, peak_y),
                xytext=(retention_time + dx, label_y),
                ha=align, va="bottom", fontsize=6.2, color=color,
                linespacing=0.92,
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                          edgecolor="none", alpha=0.94),
                arrowprops=dict(arrowstyle="-", color=color, lw=0.55,
                                shrinkA=1.5, shrinkB=1.5),
                zorder=5, clip_on=False)


def draw_trace_stack(ax, svg_path, title, result, reaction, labels, analytes,
                     trace_data=None):
    traces = trace_data if trace_data is not None else extract_traces(svg_path)
    if len(traces) != len(labels):
        raise ValueError(f"{svg_path.name}: extracted {len(traces)} traces, expected {len(labels)}")
    offsets = [ROW_SPACING * i for i in reversed(range(len(traces)))]
    for offset, segments, label in zip(offsets, traces, labels):
        ax.axhline(offset, color=LIGHT, lw=0.45, zorder=0)
        for segment in segments:
            xs, ys = zip(*segment)
            ax.plot(xs, [offset + TRACE_HEIGHT * y for y in ys], color=BLUE,
                    lw=0.75, solid_joinstyle="round", solid_capstyle="round")
        ax.text(-1.2, offset + 0.08, label, ha="right", va="center",
                fontsize=6.4, color=GREY, clip_on=False)

    # Identify molecules using retention times established by the donor/product
    # standards in the same source panel. Only target compounds are annotated.
    donor = analytes["donor"]
    product = analytes["product"]
    donor_row = analytes.get("donor_row", 2)
    product_row = analytes.get("product_row", 3)
    peak_label(ax, traces[donor_row], offsets[donor_row], donor[1], donor[0], color=BLUE,
               dx=1.2, align="left")
    peak_label(ax, traces[product_row], offsets[product_row], product[1], product[0], color=BLUE,
               dx=-1.2, align="right")
    if result == "detected":
        peak_label(ax, traces[0], offsets[0], product[1], product[0], color=RED,
                   dx=-1.1, text_height=0.92, align="right")

    result_color = RED if result == "detected" else GREY
    ax.set_title(f"{title}  |  {result}", loc="left", color=result_color,
                 fontweight="bold", pad=6)
    ax.text(0.0, -0.28, reaction, transform=ax.transAxes, ha="left", va="top",
            color=BLACK, fontsize=6.8)
    ax.set_xlim(0, 65)
    ax.set_ylim(-0.15, offsets[0] + 1.20)
    ax.set_xticks([0, 30, 60])
    ax.set_yticks([])
    ax.set_xlabel("Retention time (min)", labelpad=1)
    ax.tick_params(axis="x", direction="in", length=3, width=0.6, colors=BLACK)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.65)
    ax.spines["bottom"].set_color(BLACK)


def save(fig, stem):
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"{stem}.{suffix}", dpi=400, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)


def whole_cell():
    panels = [
        (12, "TA2  |  EC 2.6.1.11", "not detected",
         "N-acetylornithine + 2-oxoglutarate",
         {"donor": ("N-acetyl-\nornithine", 40.0), "product": ("L-glutamate", 34.3)}),
        (14, "TA10  |  EC 2.6.1.11", "detected",
         "N-acetylornithine + 2-oxoglutarate",
         {"donor": ("N-acetyl-\nornithine", 40.0), "product": ("L-glutamate", 34.3)}),
        (28, "TA5  |  EC 2.6.1.18", "not detected",
         "beta-alanine + pyruvate",
         {"donor": ("β-Alanine", 42.3), "product": ("L-alanine", 38.0)}),
        (32, "TA5  |  EC 2.6.1.44", "detected",
         "L-alanine + glyoxylate",
         {"donor": ("L-alanine", 37.9), "product": ("Glycine", 36.7)}),
    ]
    labels = ["Reaction", "Cell control", "Donor std.", "Product std."]
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 5.65))
    fig.subplots_adjust(left=0.15, right=0.985, top=0.95, bottom=0.16,
                        hspace=0.78, wspace=0.25)
    for i, (ax, (number, title, result, reaction, analytes)) in enumerate(zip(axes.flat, panels)):
        draw_trace_stack(ax, SVG / f"image{number}.svg", title, result, reaction,
                         labels, analytes)
        ax.text(-0.12, 1.06, chr(ord("a") + i), transform=ax.transAxes,
                fontweight="bold", fontsize=8.5, color=BLACK)
    fig.text(0.15, 0.025,
             "Signals were recovered from vector paths embedded in the supplied PPTX. Each trace is normalized independently;\n"
             "peak height is therefore not quantitative. Common retention-time axes and line widths are used throughout.",
             fontsize=6.6, color=GREY)
    save(fig, "si_ta_whole_cell_hplc_replotted")


def dual_function():
    """Show only candidate--reaction records backed by distinct source objects."""
    panels = [
        (17, "TA3  |  EC 2.6.1.19", "detected",
         "4-aminobutanoate + 2-oxoglutarate",
         {"donor": ("4-Aminobutanoate", 47.1),
          "product": ("L-glutamate", 34.3),
          "donor_row": 3, "product_row": 2}),
        (19, "TA4  |  EC 2.6.1.19", "detected",
         "4-aminobutanoate + 2-oxoglutarate",
         {"donor": ("4-Aminobutanoate", 47.1),
          "product": ("L-glutamate", 34.3),
          "donor_row": 3, "product_row": 2}),
        (22, "TA3  |  EC 2.6.1.22-related", "detected",
         "(S)-3-amino-2-methylpropanoate\n+ 2-oxoglutarate",
         {"donor": ("(S)-3-amino-2-\nmethylpropanoate", 46.9),
          "product": ("L-glutamate", 34.3),
          "donor_row": 3, "product_row": 2}),
        ("slide10_ta4", "TA4  |  EC 2.6.1.22-related", "detected",
         "(S)-3-amino-2-methylpropanoate\n+ 2-oxoglutarate",
         {"donor": ("(S)-3-amino-2-\nmethylpropanoate", 46.9),
          "product": ("L-glutamate", 34.3),
          "donor_row": 3, "product_row": 2}),
    ]
    labels = ["Reaction", "Cell control", "Product std.", "Donor std."]
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 6.10))
    fig.subplots_adjust(left=0.15, right=0.985, top=0.95, bottom=0.20,
                        hspace=0.88, wspace=0.26)
    for i, (ax, (source, title, result, reaction, analytes)) in enumerate(zip(axes.flat, panels)):
        if source == "slide10_ta4":
            path = SVG / "slide10_full.svg"
            traces = extract_slide10_ta4_traces(path)
        else:
            path = SVG / f"image{source}.svg"
            traces = None
        draw_trace_stack(ax, path, title, result, reaction, labels, analytes,
                         trace_data=traces)
        ax.text(-0.12, 1.06, chr(ord("a") + i), transform=ax.transAxes,
                fontweight="bold", fontsize=8.5, color=BLACK)
    fig.text(
        0.15, 0.025,
        "Distinct vector records are shown for each candidate--reaction pair. Each trace is normalized independently; peak height is qualitative.",
        fontsize=6.5, color=GREY,
    )
    save(fig, "si_ta_dual_function_hplc_replotted")


def purified():
    panels = [
        (37, "TA1  |  EC 2.6.1.13-related", "detected",
         "L-ornithine + 2-oxobutyrate", ["Reaction", "Product std.", "Donor std."],
         {"donor": ("L-ornithine", 42.9), "product": ("2-Aminobutyrate", 53.7),
          "donor_row": 2, "product_row": 1}),
        (41, "TA6  |  EC 2.6.1.96 activity", "detected",
         "4-aminobutanoate + pyruvate", ["Reaction", "Donor std.", "Product std."],
         {"donor": ("4-Aminobutanoate", 47.1), "product": ("L-alanine", 38.4),
          "donor_row": 1, "product_row": 2}),
    ]
    fig, axes = plt.subplots(2, 1, figsize=(7.25, 4.7))
    fig.subplots_adjust(left=0.15, right=0.985, top=0.93, bottom=0.17, hspace=0.84)
    for i, (ax, (number, title, result, reaction, labels, analytes)) in enumerate(zip(axes, panels)):
        draw_trace_stack(ax, SVG / f"image{number}.svg", title, result, reaction,
                         labels, analytes)
        ax.text(-0.12, 1.06, chr(ord("a") + i), transform=ax.transAxes,
                fontweight="bold", fontsize=8.5, color=BLACK)
    fig.text(0.15, 0.025,
             "Vector-derived, independently normalized qualitative traces. The supplied TA6--TA9 slide contains a TA6\n"
             "reaction chromatogram; candidate-level outcomes for TA6--TA9 are reported separately in Table S6.",
             fontsize=6.6, color=GREY)
    save(fig, "si_ta_purified_hplc_replotted")


if __name__ == "__main__":
    dual_function()
    whole_cell()
    purified()

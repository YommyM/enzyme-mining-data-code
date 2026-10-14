#!/usr/bin/env python3
"""Build publication-style SI figures from the experimental assets embedded in PPTX.

The script never uses rendered PowerPoint slides.  LC--HRMS/HPLC plots and gels are
the original embedded image objects extracted from the supplied decks.  It only
trims blank canvas and places the unchanged analytical panels in a consistent grid.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image, ImageChops


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental_validation" / "extracted_plots"
LAL = DATA / "lal" / "png"
TA = DATA / "ta"
OUT = ROOT / "figures" / "experimental_validation"

BLACK = "#242424"
GREY = "#6b6b6b"
LIGHT = "#e7e7e7"
BLUE = "#5f7f99"
RED = "#a65f5b"
ORANGE = "#b68752"

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


def trim_white(path, threshold=10, pad=8):
    """Remove only near-white outer canvas; do not alter the plotted data."""
    image = Image.open(path).convert("RGB")
    diff = ImageChops.difference(image, Image.new("RGB", image.size, "white")).convert("L")
    bbox = diff.point(lambda x: 0 if x < threshold else 255).getbbox()
    if bbox is None:
        return image
    left, top, right, bottom = bbox
    return image.crop((max(0, left - pad), max(0, top - pad),
                       min(image.width, right + pad), min(image.height, bottom + pad)))


def show_asset(ax, path, *, trim=True):
    image = trim_white(path) if trim else Image.open(path).convert("RGB")
    ax.imshow(image, interpolation="none")
    ax.set_axis_off()
    return ax


def panel_label(ax, label):
    ax.text(-0.02, 1.03, label, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8.5, fontweight="bold", color=BLACK, clip_on=False)


def save(fig, stem):
    OUT.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"{stem}.{suffix}", dpi=400, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)


def lal_analytical_evidence():
    # One representative product per candidate; the complete product list is tabulated.
    rows = [
        ("GHU", "Ser–Cys", "11,903", 31, 34),
        ("MBD", "Gln–Arg", "223,274", 35, 36),
        ("MCQ", "Gly–Ala", "350,524", 45, 46),
    ]
    fig = plt.figure(figsize=(7.25, 4.35))
    gs = fig.add_gridspec(len(rows) + 1, 3, width_ratios=[1.15, 3.0, 3.0],
                          height_ratios=[0.24] + [1] * len(rows), hspace=0.16, wspace=0.06)

    headers = ["Candidate / product", "Extracted-ion chromatogram", "HRMS spectrum"]
    for col, title in enumerate(headers):
        ax = fig.add_subplot(gs[0, col]); ax.set_axis_off()
        ax.text(0 if col == 0 else 0.5, 0.45, title,
                ha="left" if col == 0 else "center", va="center",
                fontweight="bold", color=BLACK)

    letters = "abcdefghij"
    group_colors = {"GHU": BLUE, "MBD": ORANGE, "MCQ": RED}
    for i, (enzyme, product, response, lc, ms) in enumerate(rows):
        label_ax = fig.add_subplot(gs[i + 1, 0]); label_ax.set_axis_off()
        color = group_colors[enzyme]
        label_ax.add_patch(Rectangle((0.00, 0.12), 0.045, 0.76, color=color,
                                     transform=label_ax.transAxes, clip_on=False))
        label_ax.text(0.09, 0.66, f"{letters[i]}  {enzyme}", fontweight="bold",
                      color=BLACK, transform=label_ax.transAxes, va="center")
        label_ax.text(0.09, 0.41, product, color=BLACK,
                      transform=label_ax.transAxes, va="center")
        label_ax.text(0.09, 0.18, f"response {response}", color=GREY,
                      fontsize=6.8, transform=label_ax.transAxes, va="center")
        for col, image_no in ((1, lc), (2, ms)):
            ax = fig.add_subplot(gs[i + 1, col])
            show_asset(ax, LAL / f"image{image_no}.png")
            for spine in ax.spines.values():
                spine.set_visible(True); spine.set_linewidth(0.45); spine.set_color(LIGHT)
    fig.text(0.04, 0.012,
             "One representative product is shown per candidate. Response = arbitrary instrument units, not concentration.",
             color=GREY, fontsize=6.8)
    save(fig, "si_lal_analytical_evidence")


def ta_expression():
    src = TA / "source"
    fig = plt.figure(figsize=(7.25, 7.3))
    gs = fig.add_gridspec(3, 4, height_ratios=[1.35, 1, 1.35], hspace=0.34, wspace=0.16)

    ax = fig.add_subplot(gs[0, :])
    image = Image.open(src / "image7.jpg").convert("RGB")
    # Remove photographed blank gel area while retaining every sample lane.
    image = image.crop((70, 350, 1215, 900))
    ax.imshow(image, interpolation="none"); ax.set_axis_off(); panel_label(ax, "a")
    ax.set_title("Initial expression screen (TA1--TA10; source gel)", loc="left", pad=3)

    panels = [
        ("b", "Lower-temperature expression", "image8.png"),
        ("c", "C41 host", "image9.png"),
        ("d", "BL21(DE3)pLysS host", "image10.png"),
    ]
    # Explicit spans keep the three troubleshooting gels at readable widths.
    for slc, item in zip((gs[1, 0:2], gs[1, 2], gs[1, 3]), panels):
        letter, title, name = item
        ax = fig.add_subplot(slc); show_asset(ax, src / name, trim=False)
        panel_label(ax, letter); ax.set_title(title, loc="left", pad=3)

    ax = fig.add_subplot(gs[2, 0:2])
    show_asset(ax, src / "image11.png", trim=False); panel_label(ax, "e")
    ax.set_title("Purified soluble candidates", loc="left", pad=3)

    ax = fig.add_subplot(gs[2, 2:]); ax.set_axis_off(); panel_label(ax, "f")
    ax.set_title("Assay routing", loc="left", pad=3)
    outcomes = [
        ("TA1, TA6--TA9", "soluble", BLUE),
        ("TA2--TA5, TA10", "predominantly insoluble", ORANGE),
        ("TA1, TA6--TA9", "purified-enzyme assays", RED),
        ("TA2--TA5, TA10", "whole-cell assays", GREY),
    ]
    for i, (candidates, outcome, color) in enumerate(outcomes):
        y = 0.82 - i * 0.21
        ax.add_patch(Rectangle((0.01, y - 0.055), 0.025, 0.11, color=color,
                               transform=ax.transAxes))
        ax.text(0.06, y + 0.025, candidates, transform=ax.transAxes,
                fontweight="bold", color=BLACK, va="center")
        ax.text(0.06, y - 0.045, outcome, transform=ax.transAxes,
                color=GREY, va="center")
    save(fig, "si_ta_expression")


def ta_whole_cell_hplc():
    png = TA / "png"
    panels = [
        ("TA2", "EC 2.6.1.11", "N-acetylornithine + 2-oxoglutarate", "not detected", 12),
        ("TA10", "EC 2.6.1.11", "N-acetylornithine + 2-oxoglutarate", "detected", 14),
        ("TA5", "EC 2.6.1.18", "beta-alanine + pyruvate", "not detected", 28),
        ("TA5", "EC 2.6.1.44", "L-alanine + glyoxylate", "detected", 32),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 5.25))
    fig.subplots_adjust(left=0.055, right=0.985, top=0.965, bottom=0.15,
                        hspace=0.58, wspace=0.14)
    for i, (ax, (cand, ec, reaction, result, image_no)) in enumerate(zip(axes.flat, panels)):
        show_asset(ax, png / f"image{image_no}.png")
        panel_label(ax, chr(ord("a") + i))
        outcome_color = RED if result == "detected" else GREY
        ax.set_title(f"{cand}  |  {ec}  |  {result}", loc="left", pad=3,
                     color=outcome_color, fontweight="bold")
        ax.text(0.0, -0.08, reaction, transform=ax.transAxes, ha="left", va="top",
                color=BLACK, fontsize=6.8)
    fig.text(0.055, 0.022,
             "Whole-cell conditions: 50 mg mL$^{-1}$ wet cells, 10 mM acceptor, 100 mM donor, 0.5 mM PLP,\n"
             "100 mM KPi pH 8.0, 30 $^\circ$C, 3 h. Traces are unchanged embedded exports; only blank margins were removed.",
             color=GREY, fontsize=6.6)
    save(fig, "audit_ta_whole_cell_source_panels")


def ta_purified_hplc():
    png = TA / "png"
    panels = [
        ("TA1", "EC 2.6.1.13-related", "L-ornithine + 2-oxobutyrate", "detected", 37,
         "2 mg mL$^{-1}$ enzyme; 10/100 mM acceptor/donor; 1 h"),
        ("TA6--TA9", "EC 2.6.1.96", "4-aminobutanoate + pyruvate", "all detected", 41,
         "3 mg mL$^{-1}$ enzyme; 2/10 mM acceptor/donor; overnight"),
    ]
    fig, axes = plt.subplots(2, 1, figsize=(7.25, 5.1))
    fig.subplots_adjust(left=0.055, right=0.985, top=0.94, bottom=0.12, hspace=0.44)
    for i, (ax, (cand, ec, reaction, result, image_no, condition)) in enumerate(zip(axes, panels)):
        show_asset(ax, png / f"image{image_no}.png")
        panel_label(ax, chr(ord("a") + i))
        ax.set_title(f"{cand}  |  {ec}  |  {result}", loc="left", pad=3,
                     color=RED, fontweight="bold")
        ax.text(0.0, -0.08, f"{reaction}   |   {condition}", transform=ax.transAxes,
                ha="left", va="top", color=BLACK, fontsize=6.8)
    fig.text(0.055, 0.025,
             "Both assays used 0.5 mM PLP in 100 mM KPi, pH 8.0, at 30 $^\circ$C. The source deck reports qualitative\n"
             "detection only; no conversion values or replicate statistics were supplied.",
             color=GREY, fontsize=6.6)
    save(fig, "audit_ta_purified_source_panels")


def main():
    lal_analytical_evidence()
    ta_expression()


if __name__ == "__main__":
    main()

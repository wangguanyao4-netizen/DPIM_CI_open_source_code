"""Reproduce the nonlinear core-tail coverage mechanism figure from formal data."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parent
RUN_ROOT = ROOT.parent / "data/publication"
OUTPUT_ROOT = None
M_VALUE = 50000
METHODS = ("Student t", "percentile bootstrap", "bootstrap-t")
R_VALUES = (16, 128)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def pchip_slopes(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Fritsch-Carlson slopes for shape-preserving cubic interpolation."""
    h = np.diff(x)
    delta = np.diff(y) / h
    slopes = np.zeros_like(y)

    for i in range(1, len(y) - 1):
        if delta[i - 1] * delta[i] <= 0:
            slopes[i] = 0.0
        else:
            w1 = 2.0 * h[i] + h[i - 1]
            w2 = h[i] + 2.0 * h[i - 1]
            slopes[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i])

    slopes[0] = ((2.0 * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (h[0] + h[1])
    if np.sign(slopes[0]) != np.sign(delta[0]):
        slopes[0] = 0.0
    elif np.sign(delta[0]) != np.sign(delta[1]) and abs(slopes[0]) > 3.0 * abs(delta[0]):
        slopes[0] = 3.0 * delta[0]

    slopes[-1] = ((2.0 * h[-1] + h[-2]) * delta[-1] - h[-1] * delta[-2]) / (h[-1] + h[-2])
    if np.sign(slopes[-1]) != np.sign(delta[-1]):
        slopes[-1] = 0.0
    elif np.sign(delta[-1]) != np.sign(delta[-2]) and abs(slopes[-1]) > 3.0 * abs(delta[-1]):
        slopes[-1] = 3.0 * delta[-1]
    return slopes


def pchip(x: np.ndarray, y: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    slopes = pchip_slopes(x, y)
    idx = np.searchsorted(x, x_new, side="right") - 1
    idx = np.clip(idx, 0, len(x) - 2)
    h = x[idx + 1] - x[idx]
    t = (x_new - x[idx]) / h
    h00 = 2.0 * t**3 - 3.0 * t**2 + 1.0
    h10 = t**3 - 2.0 * t**2 + t
    h01 = -2.0 * t**3 + 3.0 * t**2
    h11 = t**3 - t**2
    return h00 * y[idx] + h10 * h * slopes[idx] + h01 * y[idx + 1] + h11 * h * slopes[idx + 1]




def load_formal_data():
    raise RuntimeError("Use replay.py; a validated CSV adapter is required.")


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman"],
            "mathtext.fontset": "stix",
            "font.size": 9.0,
            "text.color": "black",
            "axes.labelcolor": "black",
            "axes.edgecolor": "black",
            "xtick.color": "black",
            "ytick.color": "black",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def region_spans(y: np.ndarray, density: np.ndarray) -> dict[str, float]:
    threshold = 0.05 * float(np.max(density))
    core = np.flatnonzero(density >= threshold)
    if core.size == 0 or core[0] == 0 or core[-1] == len(y) - 1:
        raise RuntimeError("The formal response grid does not bracket the 5% core-density boundary.")
    return {
        "threshold": threshold,
        "left_tail": float(y[core[0] - 1]),
        "left_core": float(y[core[0]]),
        "right_core": float(y[core[-1]]),
        "right_tail": float(y[core[-1] + 1]),
    }


def add_background(ax: plt.Axes, x_min: float, x_max: float, spans: dict[str, float]) -> None:
    tail_color = "#f7dddd"
    transition_color = "#fff0d7"
    core_color = "#e9f3e7"
    ax.axvspan(x_min, spans["left_tail"], color=tail_color, alpha=0.66, zorder=0)
    ax.axvspan(spans["left_tail"], spans["left_core"], color=transition_color, alpha=0.72, zorder=0)
    ax.axvspan(spans["left_core"], spans["right_core"], color=core_color, alpha=0.72, zorder=0)
    ax.axvspan(spans["right_core"], spans["right_tail"], color=transition_color, alpha=0.72, zorder=0)
    ax.axvspan(spans["right_tail"], x_max, color=tail_color, alpha=0.66, zorder=0)
    for boundary in (
        spans["left_tail"],
        spans["left_core"],
        spans["right_core"],
        spans["right_tail"],
    ):
        ax.axvline(boundary, color="#e58b31", linewidth=0.75, linestyle=(0, (4, 4)), zorder=1)


def arrow_axes(ax: plt.Axes, x_min: float, x_max: float, y_max: float) -> None:
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="both", which="both", length=0)
    ax.annotate(
        "",
        xy=(x_max, 0),
        xytext=(x_min, 0),
        arrowprops={"arrowstyle": "-|>", "color": "black", "lw": 0.85, "mutation_scale": 9},
        annotation_clip=False,
    )
    ax.annotate(
        "",
        xy=(x_min, y_max),
        xytext=(x_min, 0),
        arrowprops={"arrowstyle": "-|>", "color": "black", "lw": 0.85, "mutation_scale": 9},
        annotation_clip=False,
    )


def save_source_data(
    y: np.ndarray,
    density: np.ndarray,
    coverage_rows: list[dict[str, object]],
    spans: dict[str, float],
) -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    with (OUTPUT_ROOT / "nonlinear_reference_density.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["y", "reference_smoothed_density", "core_threshold", "region"])
        for yi, pi in zip(y, density, strict=True):
            if yi <= spans["left_tail"] or yi >= spans["right_tail"]:
                region = "tail"
            elif yi < spans["left_core"] or yi > spans["right_core"]:
                region = "transition"
            else:
                region = "core"
            writer.writerow([f"{yi:.15g}", f"{pi:.15g}", f"{spans['threshold']:.15g}", region])

    with (OUTPUT_ROOT / "nonlinear_formula_experiment_error.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["method", "R", "B", "M", "y", "coverage", "formula_prediction", "absolute_error"]
        )
        for row in coverage_rows:
            experiment = float(row["coverage"])
            formula = float(row["formula_prediction"])
            writer.writerow(
                [
                    row["method"],
                    row["R"],
                    row["B"],
                    row["M"],
                    row["y"],
                    f"{experiment:.15g}",
                    f"{formula:.15g}",
                    f"{abs(formula - experiment):.15g}",
                ]
            )


def plot() -> None:
    configure_style()
    y, density, coverage_rows, config = load_formal_data()
    spans = region_spans(y, density)
    save_source_data(y, density, coverage_rows, spans)

    x_range = float(y[-1] - y[0])
    x_min = float(y[0] - 0.045 * x_range)
    x_max = float(y[-1] + 0.045 * x_range)
    x_dense = np.linspace(y[0], y[-1], 900)
    density_dense = np.maximum(pchip(y, density, x_dense), 0.0)

    fig = plt.figure(figsize=(10.6, 6.25), facecolor="white")
    grid = fig.add_gridspec(3, 1, height_ratios=[2.55, 1.30, 2.30], hspace=0.04)
    ax_top = fig.add_subplot(grid[0])
    ax_note = fig.add_subplot(grid[1], sharex=ax_top)
    ax_bottom = fig.add_subplot(grid[2], sharex=ax_top)

    # Panel (a): nonlinear reference smoothed density and selected ordinates.
    top_max = 1.12 * float(np.max(density))
    add_background(ax_top, x_min, x_max, spans)
    ax_top.plot(x_dense, density_dense, color="black", linewidth=1.35, zorder=3)
    selected = np.rint(np.array([0, 3, 6, 9, 13, 18, 24]) * (len(y) - 1) / 24).astype(int)
    # Place labels on the open side of the local density slope so glyphs never cross the curve or stems.
    label_alignments = [(-4, 10), (4, -17), (4, -17), (4, 10), (4, 10), (4, 10), (-4, 10)]
    for number, (index, offset) in enumerate(zip(selected, label_alignments, strict=True), start=1):
        xi = float(y[index])
        yi = float(density[index])
        ax_top.vlines(xi, 0, yi, color="#333333", linewidth=0.65, zorder=2)
        ax_top.plot(xi, 0, "o", color="black", markersize=3.4, clip_on=False, zorder=4)
        ax_top.plot(xi, yi, "o", markerfacecolor="white", markeredgecolor="black", markersize=4.0, zorder=4)
        ax_top.annotate(
            rf"$\widehat{{p}}_h(y_{{{number}}})$",
            xy=(xi, yi),
            xytext=offset,
            textcoords="offset points",
            ha="left" if offset[0] >= 0 else "right",
            va="bottom" if offset[1] >= 0 else "top",
            fontsize=8.2,
            color="black",
        )
        ax_top.text(xi, -0.045 * top_max, rf"$y_{{{number}}}$", ha="center", va="top", fontsize=8.3)

    ax_top.text((x_min + spans["left_tail"]) / 2, 0.32 * top_max, "tail", ha="center", color="black")
    ax_top.text((spans["left_core"] + spans["right_core"]) / 2, 0.075 * top_max, "core", ha="center", color="black")
    ax_top.text((spans["right_tail"] + x_max) / 2, 0.32 * top_max, "tail", ha="center", color="black")
    ax_top.set_xlim(x_min, x_max)
    ax_top.set_ylim(-0.09 * top_max, top_max)
    ax_top.set_xticks([])
    ax_top.set_yticks([0])
    ax_top.set_yticklabels(["0"])
    arrow_axes(ax_top, x_min, x_max, top_max)
    ax_top.text(
        x_min - 0.055 * x_range,
        0.50 * top_max,
        r"$p_h(y)$",
        ha="center",
        va="center",
        rotation=90,
        rotation_mode="anchor",
        fontsize=12,
    )
    ax_top.text(x_max, -0.045 * top_max, r"$y$", ha="right", va="top", fontsize=10)
    ax_top.text(x_min - 0.055 * x_range, 1.01 * top_max, "(a)", fontsize=11, fontweight="bold")

    # Explanatory strip; all wording remains black while arrows encode the regimes.
    ax_note.set_ylim(0, 1)
    ax_note.axis("off")
    for boundary in (
        spans["left_tail"],
        spans["left_core"],
        spans["right_core"],
        spans["right_tail"],
    ):
        ax_note.axvline(boundary, color="#e58b31", linewidth=0.75, linestyle=(0, (4, 4)))
    left_mid = 0.5 * (x_min + spans["left_tail"])
    core_mid = 0.5 * (spans["left_core"] + spans["right_core"])
    right_mid = 0.5 * (spans["right_tail"] + x_max)
    tail_text = "$|\\gamma_n(y)|\\;\u2191$\n$|\\kappa_n(y)|\\;\u2191$\npilot moment uncertainty\nhigher order remainder"
    ax_note.text(left_mid, 0.61, tail_text, ha="center", va="center", linespacing=1.35, color="black")
    ax_note.text(right_mid, 0.61, tail_text, ha="center", va="center", linespacing=1.35, color="black")
    ax_note.text(core_mid, 0.68, "core region\nsecond-order prediction", ha="center", va="center", linespacing=1.35, color="black")
    ax_note.annotate("", xy=(spans["right_core"] - 0.03 * x_range, 0.20), xytext=(spans["left_core"] + 0.03 * x_range, 0.20), arrowprops={"arrowstyle": "<|-|>", "color": "#217a35", "lw": 0.9, "mutation_scale": 9})
    ax_note.annotate("", xy=(x_min + 0.03 * x_range, 0.08), xytext=(spans["left_tail"] - 0.02 * x_range, 0.08), arrowprops={"arrowstyle": "-|>", "color": "#c92323", "lw": 0.8, "linestyle": (0, (4, 4)), "mutation_scale": 8})
    ax_note.annotate("", xy=(x_max - 0.03 * x_range, 0.08), xytext=(spans["right_tail"] + 0.02 * x_range, 0.08), arrowprops={"arrowstyle": "-|>", "color": "#c92323", "lw": 0.8, "linestyle": (0, (4, 4)), "mutation_scale": 8})

    # Panel (b): actual formal M=5000 factorized-formula errors.
    add_background(ax_bottom, x_min, x_max, spans)
    colors = {"Student t": "#155aa6", "percentile bootstrap": "#cf2d2d", "bootstrap-t": "#15864b"}
    styles = {"Student t": "-", "percentile bootstrap": (0, (5, 3)), "bootstrap-t": "-."}
    widths = {R_VALUES[0]: 1.35, R_VALUES[1]: 0.70}
    alphas = {R_VALUES[0]: 1.0, R_VALUES[1]: 0.72}
    max_error = 0.0
    for method in METHODS:
        for r_value in R_VALUES:
            rows = [row for row in coverage_rows if row["method"] == method and int(row["R"]) == r_value]
            rows.sort(key=lambda row: float(row["y"]))
            x_values = np.array([float(row["y"]) for row in rows])
            errors = np.array([abs(float(row["formula_prediction"]) - float(row["coverage"])) for row in rows])
            max_error = max(max_error, float(np.max(errors)))
            curve = np.maximum(pchip(x_values, errors, x_dense), 0.0)
            ax_bottom.plot(
                x_dense,
                curve,
                color=colors[method],
                linestyle=styles[method],
                linewidth=widths[r_value],
                alpha=alphas[r_value],
                solid_capstyle="round",
                zorder=3 if r_value == R_VALUES[0] else 2,
            )

    bottom_max = max(0.03, 1.50 * max_error)
    ax_bottom.set_ylim(-0.07 * bottom_max, bottom_max)
    ax_bottom.set_xticks(y[selected])
    ax_bottom.set_xticklabels([rf"$y_{{{i}}}$" for i in range(1, 8)])
    tick_values = np.array([0.0, 0.01, 0.02, 0.03])
    tick_values = tick_values[tick_values <= bottom_max + 1e-12]
    ax_bottom.set_yticks(tick_values)
    ax_bottom.set_yticklabels(["0" if value == 0 else f"{value:.2f}" for value in tick_values])
    arrow_axes(ax_bottom, x_min, x_max, bottom_max)
    ax_bottom.text(
        -0.060,
        0.50,
        r"$|\widehat C_m^{\mathrm{TC}}(y)-\widehat C_m^{\mathrm{exp}}(y)|$",
        transform=ax_bottom.transAxes,
        ha="center",
        va="center",
        rotation=90,
        rotation_mode="anchor",
        fontsize=11.0,
        clip_on=False,
    )
    ax_bottom.text(x_max, -0.045 * bottom_max, r"$y$", ha="right", va="top", fontsize=10)
    ax_bottom.text(x_min - 0.055 * x_range, 1.01 * bottom_max, "(b)", fontsize=11, fontweight="bold")

    method_handles = [
        Line2D([0], [0], color=colors[method], linestyle=styles[method], linewidth=1.25, label={"Student t": r"$t$", "percentile bootstrap": r"$p$", "bootstrap-t": r"$bt$"}[method])
        for method in METHODS
    ]
    method_legend = ax_bottom.legend(
        handles=method_handles,
        loc="upper center",
        bbox_to_anchor=(0.41, 0.98),
        ncol=3,
        frameon=False,
        handlelength=3.0,
        columnspacing=1.8,
        fontsize=10.0,
    )
    ax_bottom.add_artist(method_legend)
    r_handles = [
        Line2D([0], [0], color="black", linewidth=widths[r_value], alpha=alphas[r_value], label=rf"$R={r_value}$")
        for r_value in R_VALUES
    ]
    ax_bottom.legend(
        handles=r_handles,
        loc="upper center",
        bbox_to_anchor=(0.66, 0.98),
        frameon=False,
        handlelength=3.0,
        fontsize=10.0,
    )

    # Shared x limits are useful for region alignment, but only panel (b) owns x tick labels.
    ax_top.tick_params(axis="x", labelbottom=False)
    ax_note.tick_params(axis="x", labelbottom=False)

    fig.subplots_adjust(left=0.105, right=0.985, top=0.975, bottom=0.085)
    stem = OUTPUT_ROOT / "nonlinear_real_data_core_tail_mechanism"
    pdf_path = stem.with_suffix(".pdf")
    fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.025)
    plt.close(fig)

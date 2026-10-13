"""Create Euler-beam and Kirchhoff-plate core-tail figures from formal data."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from plot_nonlinear_real_data_mechanism import (
    METHODS,
    add_background,
    arrow_axes,
    configure_style,
    pchip,
    read_csv,
    region_spans,
)


ROOT = Path(__file__).resolve().parent
RUN_ROOT = ROOT.parent / "data/publication"
M_VALUE = 30000


def load_model_data(model):
    raise RuntimeError("Use replay.py; a validated CSV adapter is required.")


MODEL_SETTINGS = {
    "linear": {
        "label": "Euler beam",
        "R": (32, 128),
        "B": 399,
        "output": None,
        "stem": "linear_real_data_core_tail_mechanism",
    },
    "plate": {
        "label": "Kirchhoff plate",
        "R": (16, 128),
        "B": 999,
        "output": None,
        "stem": "plate_real_data_core_tail_mechanism",
    },
}




def selected_locations(coverage_rows: list[dict[str, str]]) -> tuple[np.ndarray, np.ndarray]:
    y = np.array(sorted({float(row["y"]) for row in coverage_rows}))
    indices = np.unique(np.rint(np.linspace(0, len(y) - 1, 7)).astype(int))
    if len(indices) != 7:
        raise RuntimeError("Exactly seven display locations are required.")
    return y, indices


def label_offset(x: np.ndarray, density: np.ndarray, index: int) -> tuple[int, int]:
    if index == 0 or index == len(x) - 1:
        return (-4, 10)
    if density[index] < 0.20 * np.max(density):
        return (-4, 10)
    local_slope = density[index + 1] - density[index - 1]
    return (4, -17) if local_slope > 0 else (4, 10)


def write_source_data(
    output_root: Path,
    model: str,
    density_y: np.ndarray,
    density: np.ndarray,
    display_y: np.ndarray,
    spans: dict[str, float],
    coverage_rows: list[dict[str, str]],
) -> None:
    output_root.mkdir(parents=True, exist_ok=True)
    display_density = pchip(density_y, density, display_y)
    with (output_root / f"{model}_reference_density.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["y", "reference_smoothed_density", "core_threshold"])
        for yi, pi in zip(density_y, density, strict=True):
            writer.writerow([f"{yi:.15g}", f"{pi:.15g}", f"{spans['threshold']:.15g}"])

    with (output_root / f"{model}_display_locations.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["y", "interpolated_reference_density", "region"])
        for yi, pi in zip(display_y, display_density, strict=True):
            if yi <= spans["left_tail"] or yi >= spans["right_tail"]:
                region = "tail"
            elif yi < spans["left_core"] or yi > spans["right_core"]:
                region = "transition"
            else:
                region = "core"
            writer.writerow([f"{yi:.15g}", f"{pi:.15g}", region])

    with (output_root / f"{model}_formula_experiment_error.csv").open(
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


def normalize_svg(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8")


def plot_model(model: str) -> None:
    settings = MODEL_SETTINGS[model]
    density_y, density, coverage_rows, metadata = load_model_data(model)
    display_y, selected = selected_locations(coverage_rows)
    display_density = np.maximum(pchip(density_y, density, display_y), 0.0)
    spans = region_spans(display_y, display_density)
    write_source_data(
        settings["output"], model, density_y, density, display_y, spans, coverage_rows
    )

    x_range = float(display_y[-1] - display_y[0])
    x_min = float(display_y[0] - 0.045 * x_range)
    x_max = float(display_y[-1] + 0.045 * x_range)
    x_dense = np.linspace(display_y[0], display_y[-1], 900)
    density_dense = np.maximum(pchip(density_y, density, x_dense), 0.0)

    fig = plt.figure(figsize=(10.6, 6.25), facecolor="white")
    grid = fig.add_gridspec(3, 1, height_ratios=[2.55, 1.30, 2.30], hspace=0.04)
    ax_top = fig.add_subplot(grid[0])
    ax_note = fig.add_subplot(grid[1], sharex=ax_top)
    ax_bottom = fig.add_subplot(grid[2], sharex=ax_top)

    top_max = 1.12 * float(np.max(density_dense))
    add_background(ax_top, x_min, x_max, spans)
    ax_top.plot(x_dense, density_dense, color="black", linewidth=1.35, zorder=3)
    for number, index in enumerate(selected, start=1):
        xi = float(display_y[index])
        yi = float(display_density[index])
        offset = label_offset(display_y, display_density, int(index))
        ax_top.vlines(xi, 0, yi, color="#333333", linewidth=0.65, zorder=2)
        ax_top.plot(xi, 0, "o", color="black", markersize=3.4, clip_on=False, zorder=4)
        ax_top.plot(
            xi,
            yi,
            "o",
            markerfacecolor="white",
            markeredgecolor="black",
            markersize=4.0,
            zorder=4,
        )
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
        ax_top.text(
            xi,
            -0.045 * top_max,
            rf"$y_{{{number}}}$",
            ha="center",
            va="top",
            fontsize=8.3,
        )

    ax_top.text(
        (x_min + spans["left_tail"]) / 2,
        0.32 * top_max,
        "tail",
        ha="center",
        color="black",
    )
    ax_top.text(
        (spans["left_core"] + spans["right_core"]) / 2,
        0.075 * top_max,
        "core",
        ha="center",
        color="black",
    )
    ax_top.text(
        (spans["right_tail"] + x_max) / 2,
        0.32 * top_max,
        "tail",
        ha="center",
        color="black",
    )
    ax_top.set_xlim(x_min, x_max)
    ax_top.set_ylim(-0.09 * top_max, top_max)
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
    ax_top.text(
        x_min - 0.055 * x_range,
        1.01 * top_max,
        "(a)",
        fontsize=11,
        fontweight="bold",
    )

    ax_note.set_ylim(0, 1)
    ax_note.axis("off")
    for boundary in (
        spans["left_tail"],
        spans["left_core"],
        spans["right_core"],
        spans["right_tail"],
    ):
        # Keep the regime boundary visible without drawing it through annotation glyphs.
        ax_note.plot(
            [boundary, boundary],
            [0.00, 0.34],
            color="#e58b31",
            linewidth=0.75,
            linestyle=(0, (4, 4)),
        )
        ax_note.plot(
            [boundary, boundary],
            [0.88, 1.00],
            color="#e58b31",
            linewidth=0.75,
            linestyle=(0, (4, 4)),
        )
    left_mid = 0.5 * (x_min + spans["left_tail"])
    core_mid = 0.5 * (spans["left_core"] + spans["right_core"])
    right_mid = 0.5 * (spans["right_tail"] + x_max)
    tail_text = "$|\\gamma_n(y)|\\;\u2191$\n$|\\kappa_n(y)|\\;\u2191$\npilot moment uncertainty\nhigher order remainder"
    ax_note.text(
        left_mid,
        0.61,
        tail_text,
        ha="center",
        va="center",
        linespacing=1.35,
        color="black",
    )
    ax_note.text(
        right_mid,
        0.61,
        tail_text,
        ha="center",
        va="center",
        linespacing=1.35,
        color="black",
    )
    ax_note.text(
        core_mid,
        0.68,
        "core region\nsecond-order prediction",
        ha="center",
        va="center",
        linespacing=1.35,
        color="black",
    )
    ax_note.annotate(
        "",
        xy=(spans["right_core"] - 0.03 * x_range, 0.20),
        xytext=(spans["left_core"] + 0.03 * x_range, 0.20),
        arrowprops={
            "arrowstyle": "<|-|>",
            "color": "#217a35",
            "lw": 0.9,
            "mutation_scale": 9,
        },
    )
    ax_note.annotate(
        "",
        xy=(x_min + 0.03 * x_range, 0.08),
        xytext=(spans["left_tail"] - 0.02 * x_range, 0.08),
        arrowprops={
            "arrowstyle": "-|>",
            "color": "#c92323",
            "lw": 0.8,
            "linestyle": (0, (4, 4)),
            "mutation_scale": 8,
        },
    )
    ax_note.annotate(
        "",
        xy=(x_max - 0.03 * x_range, 0.08),
        xytext=(spans["right_tail"] + 0.02 * x_range, 0.08),
        arrowprops={
            "arrowstyle": "-|>",
            "color": "#c92323",
            "lw": 0.8,
            "linestyle": (0, (4, 4)),
            "mutation_scale": 8,
        },
    )

    add_background(ax_bottom, x_min, x_max, spans)
    colors = {
        "Student t": "#155aa6",
        "percentile bootstrap": "#cf2d2d",
        "bootstrap-t": "#15864b",
    }
    styles = {
        "Student t": "-",
        "percentile bootstrap": (0, (5, 3)),
        "bootstrap-t": "-.",
    }
    widths = {settings["R"][0]: 1.35, settings["R"][1]: 0.70}
    alphas = {settings["R"][0]: 1.0, settings["R"][1]: 0.72}
    max_error = 0.0
    for method in METHODS:
        for r_value in settings["R"]:
            rows = [
                row
                for row in coverage_rows
                if row["method"] == method and int(row["R"]) == r_value
            ]
            rows.sort(key=lambda row: float(row["y"]))
            x_values = np.array([float(row["y"]) for row in rows])
            errors = np.array(
                [
                    abs(float(row["formula_prediction"]) - float(row["coverage"]))
                    for row in rows
                ]
            )
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
                zorder=3 if r_value == settings["R"][0] else 2,
            )

    bottom_max = 1.50 * max_error
    ax_bottom.set_ylim(-0.07 * bottom_max, bottom_max)
    ax_bottom.set_xticks(display_y[selected])
    ax_bottom.set_xticklabels([rf"$y_{{{i}}}$" for i in range(1, 8)])
    ax_bottom.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))
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
    ax_bottom.text(
        x_min - 0.055 * x_range,
        1.01 * bottom_max,
        "(b)",
        fontsize=11,
        fontweight="bold",
    )

    method_handles = [
        Line2D(
            [0],
            [0],
            color=colors[method],
            linestyle=styles[method],
            linewidth=1.50,
            label={
                "Student t": r"$t$",
                "percentile bootstrap": r"$p$",
                "bootstrap-t": r"$bt$",
            }[method],
        )
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
        Line2D(
            [0],
            [0],
            color="black",
            linewidth=widths[r_value],
            alpha=alphas[r_value],
            label=rf"$R={r_value}$",
        )
        for r_value in settings["R"]
    ]
    ax_bottom.legend(
        handles=r_handles,
        loc="upper center",
        bbox_to_anchor=(0.66, 0.98),
        frameon=False,
        handlelength=3.0,
        fontsize=10.0,
    )

    ax_top.tick_params(axis="x", labelbottom=False)
    ax_note.tick_params(axis="x", labelbottom=False)
    fig.subplots_adjust(left=0.105, right=0.985, top=0.975, bottom=0.085)
    stem = settings["output"] / settings["stem"]
    pdf_path = stem.with_suffix(".pdf")
    fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.025)
    plt.close(fig)

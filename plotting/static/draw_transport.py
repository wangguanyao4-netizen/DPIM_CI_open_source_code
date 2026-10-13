"""Rebuild the probability-weighted DPIM schematic as vector artwork."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon


TEXT_COLOR = "#111111"
EDGE_COLOR = "#747474"
ORANGE = "#E99A00"
BLUE = "#1769E8"
MAGENTA = "#D61AA4"


def clip_half_plane(
    polygon: np.ndarray, normal: np.ndarray, offset: float, tolerance: float = 1e-11
) -> np.ndarray:
    """Clip a convex polygon to normal.dot(x) <= offset."""
    if len(polygon) == 0:
        return polygon

    output: list[np.ndarray] = []
    previous = polygon[-1]
    previous_value = float(np.dot(normal, previous) - offset)
    previous_inside = previous_value <= tolerance

    for current in polygon:
        current_value = float(np.dot(normal, current) - offset)
        current_inside = current_value <= tolerance
        if current_inside != previous_inside:
            direction = current - previous
            denominator = float(np.dot(normal, direction))
            if abs(denominator) > tolerance:
                fraction = (offset - float(np.dot(normal, previous))) / denominator
                output.append(previous + fraction * direction)
        if current_inside:
            output.append(current)
        previous = current
        previous_inside = current_inside

    return np.asarray(output, dtype=float)


def bounded_voronoi(points: np.ndarray, bounds: tuple[float, float, float, float]) -> list[np.ndarray]:
    """Construct Voronoi cells by clipping the plotting rectangle by all bisectors."""
    xmin, xmax, ymin, ymax = bounds
    rectangle = np.array(
        [[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]], dtype=float
    )
    cells: list[np.ndarray] = []
    squared_norms = np.einsum("ij,ij->i", points, points)
    for index, point in enumerate(points):
        cell = rectangle.copy()
        for other_index, other in enumerate(points):
            if other_index == index:
                continue
            normal = other - point
            offset = 0.5 * (squared_norms[other_index] - squared_norms[index])
            cell = clip_half_plane(cell, normal, offset)
            if len(cell) == 0:
                break
        cells.append(cell)
    return cells


def gaussian(x: np.ndarray, center: float, width: float, amplitude: float) -> np.ndarray:
    return amplitude * np.exp(-0.5 * ((x - center) / width) ** 2)


def draw_arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str,
    curvature: float,
) -> None:
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        connectionstyle=f"arc3,rad={curvature}",
        mutation_scale=13,
        linewidth=1.25,
        color=color,
        shrinkA=1,
        shrinkB=1,
        zorder=8,
    )
    ax.add_patch(arrow)


def draw_axis_arrow(
    ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=15,
            linewidth=1.25,
            color=TEXT_COLOR,
            shrinkA=0,
            shrinkB=0,
            zorder=10,
        )
    )


def build_figure(output_dir: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "STIXGeneral",
            "font.size": 13,
            "mathtext.fontset": "stix",
            "text.color": TEXT_COLOR,
            "axes.labelcolor": TEXT_COLOR,
            "axes.edgecolor": TEXT_COLOR,
            "xtick.color": TEXT_COLOR,
            "ytick.color": TEXT_COLOR,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    fig, ax = plt.subplots(figsize=(16.8, 9.4), constrained_layout=False)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.set_xlim(0, 16.8)
    ax.set_ylim(0, 9.4)
    ax.axis("off")

    # ------------------------------------------------------------------
    # Input probability space and bounded Voronoi cells.
    # ------------------------------------------------------------------
    bounds = (0.80, 5.85, 2.20, 7.20)
    rng = np.random.default_rng(20260726)
    columns, rows = 6, 6
    gx = np.linspace(1.18, 5.48, columns)
    gy = np.linspace(2.62, 6.78, rows)
    points = np.array([(x, y) for y in gy for x in gx], dtype=float)
    points += rng.normal(scale=[0.15, 0.14], size=points.shape)
    points[:, 0] = np.clip(points[:, 0], bounds[0] + 0.20, bounds[1] - 0.20)
    points[:, 1] = np.clip(points[:, 1], bounds[2] + 0.20, bounds[3] - 0.20)
    cells = bounded_voronoi(points, bounds)

    targets = np.array([[3.95, 5.62], [2.63, 4.04], [4.44, 3.48]])
    selected = [int(np.argmin(np.linalg.norm(points - target, axis=1))) for target in targets]
    highlight_colors = {selected[0]: "#F7D77E", selected[1]: "#B9D0F5", selected[2]: "#EAAAD5"}
    pastel = ["#F0F3E8", "#ECEAF2", "#F4EBE5", "#E9F0EF", "#F2F2F2"]

    for index, cell in enumerate(cells):
        if len(cell) < 3:
            continue
        face = highlight_colors.get(index, pastel[index % len(pastel)])
        alpha = 0.78 if index in highlight_colors else 0.68
        ax.add_patch(
            Polygon(
                cell,
                closed=True,
                facecolor=face,
                edgecolor=EDGE_COLOR,
                linewidth=0.62,
                alpha=alpha,
                zorder=1,
            )
        )

    auxiliary = rng.uniform(
        [bounds[0] + 0.08, bounds[2] + 0.08],
        [bounds[1] - 0.08, bounds[3] - 0.08],
        size=(56, 2),
    )
    ax.scatter(auxiliary[:, 0], auxiliary[:, 1], s=10, color="#D7191C", linewidths=0, zorder=4)
    ax.scatter(points[:, 0], points[:, 1], s=17, color=TEXT_COLOR, linewidths=0, zorder=5)

    draw_axis_arrow(ax, (0.65, 2.05), (6.20, 2.05))
    draw_axis_arrow(ax, (0.65, 2.05), (0.65, 7.48))
    ax.text(6.34, 1.92, r"$\Theta_1$", fontsize=19, ha="left", va="center", color=TEXT_COLOR)
    ax.text(0.65, 7.65, r"$\Theta_2$", fontsize=19, ha="center", va="bottom", color=TEXT_COLOR)

    label_offsets = [(0.13, 0.44), (-0.42, 0.34), (-0.12, 0.38)]
    cell_labels = []
    symbols = [("i", "P_i", r"\theta_i"), ("j", "P_j", r"\theta_j"), ("k", "P_k", r"\theta_k")]
    for selected_index, offset, (_, p_symbol, theta_symbol) in zip(selected, label_offsets, symbols):
        px, py = points[selected_index]
        p_label = ax.text(
            px + offset[0],
            py + offset[1] + 0.20,
            rf"${p_symbol}$",
            fontsize=18,
            ha="center",
            va="center",
            color=TEXT_COLOR,
            zorder=7,
        )
        theta_label = ax.text(
            px + offset[0],
            py + offset[1] - 0.18,
            rf"${theta_symbol}$",
            fontsize=18,
            ha="center",
            va="center",
            color=TEXT_COLOR,
            zorder=7,
        )
        cell_labels.append((selected_index, p_label, theta_label))

    # Test complete glyph bounding boxes, not just the text anchor, against each cell.
    from matplotlib.path import Path as MplPath
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for index, p_label, theta_label in cell_labels:
        polygon = cells[index]
        inside = MplPath(np.vstack([polygon, polygon[0]]))
        if index == selected[1]:
            # Keep the blue labels side by side below its representative point.
            for label, position in zip((p_label, theta_label), ((2.65, 4.23), (2.94, 4.23))):
                label.set_fontsize(16)
                label.set_position(position)
                box = label.get_window_extent(renderer).transformed(ax.transData.inverted())
                corners = np.array([[box.x0,box.y0],[box.x1,box.y0],
                                    [box.x1,box.y1],[box.x0,box.y1]])
                assert inside.contains_points(corners).all(), f"Blue label leaves its cell: {corners}"
            continue
        sizes = []
        for label in (p_label, theta_label):
            box = label.get_window_extent(renderer).transformed(ax.transData.inverted())
            sizes.append((box.width + 0.035, box.height + 0.035))
        candidates = []
        for x0 in np.linspace(polygon[:, 0].min(), polygon[:, 0].max(), 65):
            for y0, dx, dy in [(y, dx, dy) for y in np.linspace(polygon[:, 1].min(), polygon[:, 1].max(), 65)
                              for dx, dy in [(0.19, 0), (0, 0.16)]]:
                centers = [(x0 - dx, y0 + dy), (x0 + dx, y0 - dy)]
                valid = True
                clearance = []
                for (cx, cy), (w, h) in zip(centers, sizes):
                    corners = np.array([[cx-w/2,cy-h/2],[cx+w/2,cy-h/2],
                                        [cx+w/2,cy+h/2],[cx-w/2,cy+h/2]])
                    if not inside.contains_points(corners).all():
                        valid = False
                        break
                    px, py = points[index]
                    arrow_y = py + 0.20 * (cx - px)
                    if cx+w/2 >= px+0.10 and cy-h/2-0.035 <= arrow_y <= cy+h/2+0.035:
                        valid = False
                        break
                    distances = np.maximum(np.abs(np.vstack([points, auxiliary]) - [cx,cy]) - [w/2,h/2], 0)
                    clearance.append(np.linalg.norm(distances, axis=1).min())
                if valid:
                    candidates.append((min(clearance), x0, y0, dx, dy))
        assert candidates, f"No valid label placement for cell {index}"
        _, cx, cy, dx, dy = max(candidates)
        p_label.set_position((cx - dx, cy + dy))
        theta_label.set_position((cx + dx, cy - dy))

    ax.text(
        0.85,
        1.10,
        r"$P_i=\int_{\Omega_{\theta_i}}p_{\Theta}(\theta)\,\mathrm{d}\theta,"
        r"\qquad \widehat{P}_i=\frac{1}{N_{\mathrm{V}}}"
        r"\sum_{j=1}^{N_{\mathrm{V}}}\mathbf{1}\!\left\{\theta_j^{\mathrm{aux}}\in\Omega_{\theta_i}\right\}$",
        fontsize=17,
        ha="left",
        va="center",
        color=TEXT_COLOR,
    )

    # ------------------------------------------------------------------
    # Response mapping and kernel reconstruction.
    # ------------------------------------------------------------------
    ax.text(7.32, 6.88, r"$Y=g(\Theta)$", fontsize=19, ha="center", va="center", color=TEXT_COLOR)

    x_axis_left, x_axis_right, y_axis = 9.20, 16.20, 2.05
    draw_axis_arrow(ax, (x_axis_left - 0.28, y_axis), (x_axis_right, y_axis))
    draw_axis_arrow(ax, (x_axis_left, 1.78), (x_axis_left, 7.38))
    ax.text(x_axis_right - 0.02, 1.78, r"$y$", fontsize=19, ha="right", va="top", color=TEXT_COLOR)
    ax.text(9.20, 8.00, "pdf", fontsize=17, ha="center", va="bottom", color=TEXT_COLOR)
    ax.text(9.20, 7.60, r"$p_h(y)$", fontsize=18, ha="center", va="bottom", color=TEXT_COLOR)

    x = np.linspace(x_axis_left, x_axis_right - 0.15, 1000)
    centers = np.array([10.15, 12.35, 13.90])
    widths = np.array([0.36, 0.50, 0.50])
    amplitudes = np.array([3.55, 2.05, 1.15])
    curve_colors = [ORANGE, BLUE, MAGENTA]

    gray_centers = np.array([9.72, 10.48, 11.10, 11.62, 12.02, 12.62, 13.15, 13.70, 14.25, 14.80, 15.22])
    gray_amplitudes = np.array([0.55, 0.78, 0.72, 0.76, 0.62, 0.68, 0.60, 0.57, 0.52, 0.48, 0.42])
    for center, amplitude in zip(gray_centers, gray_amplitudes):
        ax.plot(x, y_axis + gaussian(x, center, 0.34, amplitude), color="#C8C8C8", linewidth=0.85, zorder=2)

    component_curves = []
    for center, width, amplitude, color in zip(centers, widths, amplitudes, curve_colors):
        component = y_axis + gaussian(x, center, width, amplitude)
        component_curves.append(component)
        mask = component > y_axis + 0.015
        ax.plot(x[mask], component[mask], color=color, linewidth=1.45, zorder=5)
        peak_y = y_axis + amplitude
        ax.plot([center, center], [y_axis - 0.08, peak_y], color=color, linewidth=0.75, linestyle=(0, (2, 2)), alpha=0.78, zorder=3)

    broad = (
        gaussian(x, 10.05, 0.50, 2.55)
        + gaussian(x, 11.18, 0.86, 1.32)
        + gaussian(x, 12.35, 0.54, 1.62)
        + gaussian(x, 13.55, 0.82, 1.05)
        + gaussian(x, 14.55, 0.90, 0.62)
    )
    density = y_axis + broad
    ax.plot(x, density, color=TEXT_COLOR, linewidth=1.75, zorder=6)

    for center, label in zip(centers, [r"$g(\theta_i)$", r"$g(\theta_j)$", r"$g(\theta_k)$"]):
        ax.plot([center, center], [y_axis - 0.10, y_axis + 0.02], color=TEXT_COLOR, linewidth=0.75, zorder=8)
        ax.text(center, 1.72, label, fontsize=17, ha="center", va="top", color=TEXT_COLOR)

    formula = r"$Q_n(y)=\sum_{i=1}^{n}P_iK_h\!\left\{y-g(\theta_i)\right\}$"
    ax.text(14.05, 7.15, formula, fontsize=18, ha="center", va="center", color=TEXT_COLOR)

    label_x = 14.28
    label_y = float(np.interp(label_x, x, density)) + 0.42
    ax.annotate(
        r"$p_h(y)$",
        xy=(14.28, float(np.interp(14.28, x, density))),
        xytext=(14.78, label_y),
        fontsize=18,
        ha="left",
        va="center",
        color=TEXT_COLOR,
        arrowprops={"arrowstyle": "-", "color": TEXT_COLOR, "linewidth": 1.0},
    )

    arrow_starts = [
        tuple(points[selected[0]] + np.array([0.10, 0.02])),
        tuple(points[selected[1]] + np.array([0.10, 0.00])),
        tuple(points[selected[2]] + np.array([0.10, 0.00])),
    ]
    # Each arrow terminates on its own colored kernel, not on a different component.
    arrow_ends = [(float(c - 0.35*w), float(y_axis + a*np.exp(-0.5*0.35**2)))
                  for c,w,a in zip(centers,widths,amplitudes)]
    curvatures = [-0.28, -0.12, 0.18]
    for start, end, color, curvature in zip(arrow_starts, arrow_ends, curve_colors, curvatures):
        draw_arrow(ax, start, end, color, curvature)

    output_dir.mkdir(parents=True, exist_ok=True)
    base = output_dir / "Fig01a_DPIM_probability_transport"
    metadata = {
        "Title": "Probability-weighted DPIM response transport schematic",
        "Author": "Rebuilt with Python and Matplotlib",
        "Subject": "Voronoi probability weights and kernel density reconstruction",
    }
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.05, metadata=metadata)
    plt.close(fig)

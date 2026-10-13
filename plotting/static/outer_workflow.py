#!/usr/bin/env python3
"""Draw the DPIM outer-sample workflow from a formal randomized curve pool."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
from reportlab.graphics import renderPDF
from reportlab.graphics.shapes import Circle, Drawing, Line, Path as RLPath, Polygon, Rect, String
from reportlab.lib.colors import Color, HexColor
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont


WIDTH = 1180
HEIGHT = 650

def configure_fonts(roman_font, italic_font):
    # Explicit user-installed files preserve fonts without bundling them.
    pdfmetrics.registerFont(TTFont("PaperRoman", str(roman_font)))
    pdfmetrics.registerFont(TTFont("PaperItalic", str(italic_font)))


def load_pool(path: Path) -> np.ndarray:
    pool = np.loadtxt(path, delimiter=",")
    if pool.ndim != 2:
        raise ValueError("The curve pool must be a two-dimensional matrix.")
    return pool


def smooth_curves(y: np.ndarray, curves: np.ndarray, points_per_interval: int = 14) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate actual ordinates without altering their amplitudes."""
    kernel = np.array([0.0, 1.0, 0.0])
    vertices = np.vstack([
        np.convolve(np.pad(curve, (1, 1), mode="edge"), kernel, mode="valid")
        for curve in curves
    ])
    slopes = np.gradient(vertices, y, axis=1)
    dense_parts = []
    curve_parts = []
    for interval in range(len(y) - 1):
        t = np.linspace(0.0, 1.0, points_per_interval, endpoint=False)
        h = y[interval + 1] - y[interval]
        h00 = 2.0 * t**3 - 3.0 * t**2 + 1.0
        h10 = t**3 - 2.0 * t**2 + t
        h01 = -2.0 * t**3 + 3.0 * t**2
        h11 = t**3 - t**2
        dense_parts.append(y[interval] + t * h)
        curve_parts.append(
            vertices[:, interval, None] * h00
            + h * slopes[:, interval, None] * h10
            + vertices[:, interval + 1, None] * h01
            + h * slopes[:, interval + 1, None] * h11
        )
    dense_y = np.concatenate([*dense_parts, y[-1:]])
    dense_curves = np.concatenate([*curve_parts, vertices[:, -1, None]], axis=1)
    return dense_y, np.maximum(dense_curves, 0.0)


def data_to_page(x: np.ndarray, y: np.ndarray, box: tuple[float, float, float, float],
                 xlim: tuple[float, float], ylim: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    left, bottom, right, top = box
    xp = left + (x - xlim[0]) / (xlim[1] - xlim[0]) * (right - left)
    yp = bottom + (y - ylim[0]) / (ylim[1] - ylim[0]) * (top - bottom)
    return xp, yp


def add_polyline(drawing: Drawing, x: np.ndarray, y: np.ndarray, color: Color,
                 width: float, opacity: float = 1.0, dash: list[float] | None = None) -> None:
    path = RLPath()
    path.moveTo(float(x[0]), float(y[0]))
    for xi, yi in zip(x[1:], y[1:]):
        path.lineTo(float(xi), float(yi))
    path.strokeColor = color
    path.strokeWidth = width
    path.fillColor = None
    path.strokeOpacity = opacity
    if dash:
        path.strokeDashArray = dash
    drawing.add(path)


def arrow(drawing: Drawing, x1: float, y1: float, x2: float, y2: float,
          color: Color = HexColor("#17191C"), width: float = 1.25) -> None:
    drawing.add(Line(x1, y1, x2, y2, strokeColor=color, strokeWidth=width))
    direction = np.array([x2 - x1, y2 - y1], dtype=float)
    direction /= np.linalg.norm(direction)
    normal = np.array([-direction[1], direction[0]])
    tip = np.array([x2, y2])
    base = tip - 9.0 * direction
    p1 = base + 4.0 * normal
    p2 = base - 4.0 * normal
    drawing.add(Polygon([tip[0], tip[1], p1[0], p1[1], p2[0], p2[1]],
                        fillColor=color, strokeColor=color))


def boxed_text(drawing: Drawing, x: float, y: float, value: str, font_name: str,
               font_size: float, color: Color, padding: float = 4.0) -> None:
    width = pdfmetrics.stringWidth(value, font_name, font_size)
    drawing.add(Rect(x - padding, y - padding * 0.7, width + 2 * padding,
                     font_size + 1.7 * padding, fillColor=HexColor("#FFFFFF"),
                     strokeColor=None, fillOpacity=0.94))
    drawing.add(String(x, y, value, fontName=font_name, fontSize=font_size, fillColor=color))


def write_latex_overlay(output_dir: Path, x_star: float, label_pages: list[float],
                        q_mean_page: float) -> Path:
    labels = [
        r"$Q_1(y^*)$", r"$Q_2(y^*)$", r"$Q_r(y^*)$",
        r"$Q_s(y^*)$", r"$Q_{R-1}(y^*)$", r"$Q_R(y^*)$",
    ]
    label_nodes = "\n".join(
        rf"\node[mathlabel] at (1032,{page - 5:.4f}) {{{label}}};"
        for page, label in zip(label_pages, labels)
    )
    tex = rf"""\documentclass[border=0pt]{{standalone}}
\usepackage{{fontspec}}
\setmainfont{{Latin Modern Roman}}
\usepackage{{amsmath,bm,graphicx,tikz}}
\begin{{document}}
\begin{{tikzpicture}}[x=1bp,y=1bp]
\tikzset{{
  textlabel/.style={{anchor=base west,text=black,font=\fontsize{{17bp}}{{20bp}}\selectfont\itshape}},
  mathlabel/.style={{anchor=base west,text=black,fill=white,fill opacity=.94,text opacity=1,
                    inner sep=2.5bp,font=\fontsize{{15bp}}{{18bp}}\selectfont}},
  clearlabel/.style={{anchor=base west,text=black,fill=white,fill opacity=.94,text opacity=1,
                     inner sep=3bp,font=\fontsize{{17bp}}{{20bp}}\selectfont\itshape}}
}}
\node[anchor=south west,inner sep=0] at (0,0)
  {{\includegraphics[width=1180bp,height=650bp]{{real_data_outer_sample_workflow_background.pdf}}}};
\node[textlabel] at (80,585) {{pdf}};
\node[textlabel] at (879,102) {{$y$}};
\node[anchor=base west,text=black,font=\fontsize{{16bp}}{{19bp}}\selectfont] at ({x_star - 27:.4f},95) {{$y=y^*$}};
\node[clearlabel] at (205,515) {{$P$ curves}};
\node[clearlabel] at (655,515) {{$R$ selected curves}};
\node[textlabel] at (940,585) {{outer sample at $y^*$}};
{label_nodes}
\node[mathlabel] at (906,{q_mean_page + 14:.4f}) {{$\overline{{Q}}_R^{{\,*}}$}};
\node[anchor=base west,text=black,font=\fontsize{{13bp}}{{16bp}}\selectfont\itshape] at (897,22) {{Student t}};
\node[anchor=base west,text=black,font=\fontsize{{13bp}}{{16bp}}\selectfont\itshape] at (976,22) {{Percentile}};
\node[anchor=base west,text=black,font=\fontsize{{13bp}}{{16bp}}\selectfont\itshape] at (1065,22) {{Bootstrap-t}};
\node[anchor=base west,text=black,font=\fontsize{{16bp}}{{20bp}}\selectfont] at (76,58)
  {{$Q_r(y)=\displaystyle\sum_{{i=1}}^{{n}}P_{{r,i}}K_h\!\left\{{y-g(\bm{{\theta}}_{{r,i}})\right\}}$}};
\end{{tikzpicture}}
\end{{document}}
"""
    tex_path = output_dir / "real_data_outer_sample_workflow_latex.tex"
    tex_path.write_text(tex, encoding="utf-8")
    return tex_path


def compile_latex_figure(output_dir: Path, tex_path: Path) -> tuple[Path, Path, Path]:
    xelatex = shutil.which("xelatex")
    dvisvgm = shutil.which("dvisvgm")
    if not xelatex:
        raise FileNotFoundError(f"XeLaTeX was not found at {xelatex}")
    build_dir = output_dir / "latex_build"
    build_dir.mkdir(exist_ok=True)
    command = [
        str(xelatex), "-interaction=nonstopmode", "-halt-on-error",
        f"-output-directory={build_dir}", tex_path.name,
    ]
    try:
        for _ in range(2):
            subprocess.run(command, cwd=output_dir, check=True, capture_output=True, text=True)
        compiled_pdf = build_dir / f"{tex_path.stem}.pdf"
        final_pdf = output_dir / "real_data_outer_sample_workflow.pdf"
        shutil.copy2(compiled_pdf, final_pdf)
        final_svg = output_dir / "real_data_outer_sample_workflow.svg"
        if dvisvgm:
            subprocess.run(
                [str(dvisvgm), "--pdf", "--exact", "--font-format=woff",
                 f"--output={final_svg.name}", final_pdf.name],
                cwd=output_dir, check=True, capture_output=True, text=True,
            )
        final_png = output_dir / "real_data_outer_sample_workflow.png"
        return final_pdf, final_svg, final_png
    finally:
        shutil.rmtree(build_dir, ignore_errors=True)


def draw_figure(y: np.ndarray, pool: np.ndarray, config: dict, r_count: int,
                y_star: float, output_dir: Path, pool_path: Path,
                config_path: Path, source_mat: Path | None) -> dict:
    if pool.shape[1] != len(y):
        raise ValueError(f"Pool has {pool.shape[1]} ordinates but config has {len(y)} y values.")
    if r_count > pool.shape[0]:
        raise ValueError("R cannot exceed the available curve count.")

    display_rank_amplitude = 0.0
    display_shape_retention = 1.0
    dense_y, dense_pool = smooth_curves(y, pool)
    r_curves = pool[:r_count]
    y_star_index = int(np.argmin(np.abs(y - y_star)))
    actual_y_star = float(y[y_star_index])
    dense_y_star_index = int(np.argmin(np.abs(dense_y - actual_y_star)))
    r_actual_dense = dense_pool[:r_count]
    actual_order = np.argsort(r_actual_dense[:, dense_y_star_index])
    rank_scores = np.empty(r_count)
    rank_scores[actual_order] = np.linspace(-1.0, 1.0, r_count)
    r_dense_mean = r_actual_dense.mean(axis=0, keepdims=True)
    r_dense = (
        r_dense_mean
        + display_shape_retention * (r_actual_dense - r_dense_mean)
        + display_rank_amplitude * rank_scores[:, None] * r_dense_mean
    )
    r_dense = np.maximum(r_dense, 0.0)
    q_values = r_dense[:, dense_y_star_index]

    order = np.argsort(q_values)
    selected_local = np.array([
        order[-1], order[-2], order[round(0.66 * (r_count - 1))],
        order[round(0.33 * (r_count - 1))], order[1], order[0],
    ])
    if len(selected_local) != 6:
        raise RuntimeError("Quantile selection did not produce six distinct curves.")
    selected_local = selected_local[np.argsort(q_values[selected_local])[::-1]]

    highlight_colors = [
        HexColor("#A51C45"), HexColor("#E58B00"), HexColor("#4C8C2B"),
        HexColor("#00A6B2"), HexColor("#1B56A8"), HexColor("#8D4BC2"),
    ]
    pale_blue = HexColor("#AFC4E3")
    axis_color = HexColor("#000000")
    mean_color = HexColor("#4B5563")

    drawing = Drawing(WIDTH, HEIGHT)
    drawing.add(Rect(0, 0, WIDTH, HEIGHT, fillColor=HexColor("#FFFFFF"), strokeColor=None))

    left_box = (68.0, 132.0, 885.0, 578.0)
    xlim = (float(y.min()), float(y.max()))
    ylim = (0.0, float(np.quantile(pool, 0.999) * 1.13))

    arrow(drawing, left_box[0], left_box[1], left_box[2] + 10, left_box[1], axis_color)
    arrow(drawing, left_box[0], left_box[1], left_box[0], left_box[3] + 12, axis_color)
    drawing.add(String(left_box[0] + 12, left_box[3] + 7, "density",
                       fontName="Times-Italic", fontSize=17, fillColor=axis_color))
    drawing.add(String(left_box[2] - 6, left_box[1] - 30, "y",
                       fontName="Times-Italic", fontSize=18, fillColor=axis_color))

    xp, _ = data_to_page(dense_y, np.zeros_like(dense_y), left_box, xlim, ylim)
    for curve in dense_pool:
        _, yp = data_to_page(dense_y, curve, left_box, xlim, ylim)
        add_polyline(drawing, xp, yp, pale_blue, 0.55, opacity=0.30)

    r_mean = r_dense.mean(axis=0)
    _, mean_yp = data_to_page(dense_y, r_mean, left_box, xlim, ylim)
    add_polyline(drawing, xp, mean_yp, mean_color, 1.15, opacity=0.90, dash=[4, 3])

    x_star_page, _ = data_to_page(np.array([actual_y_star]), np.array([0.0]), left_box, xlim, ylim)
    x_star_page = float(x_star_page[0])
    drawing.add(Line(x_star_page, left_box[1] - 8, x_star_page, left_box[3] - 5,
                     strokeColor=axis_color, strokeWidth=1.1, strokeDashArray=[6, 4]))
    drawing.add(String(x_star_page - 27, left_box[1] - 37, "y = y*",
                       fontName="Times-Italic", fontSize=16, fillColor=axis_color))

    for local_index, color in zip(selected_local, highlight_colors):
        _, yp = data_to_page(dense_y, r_dense[local_index], left_box, xlim, ylim)
        add_polyline(drawing, xp, yp, color, 1.25, opacity=0.98)

    boxed_text(drawing, 205, 515, "P curves",
               "PaperItalic", 17, axis_color)
    boxed_text(drawing, 655, 515, "R selected curves",
               "PaperItalic", 17, axis_color)

    right_x = 1010.0
    right_bottom = 120.0
    right_top = 575.0
    arrow(drawing, right_x, right_bottom, right_x, right_top, axis_color)
    drawing.add(String(940, right_top + 7, "outer sample at y*",
                       fontName="Times-Italic", fontSize=16, fillColor=axis_color))

    q_min = float(q_values.min())
    q_max = float(q_values.max())
    q_pad = max(0.06 * (q_max - q_min), 1e-4)
    q_low, q_high = q_min - q_pad, q_max + q_pad

    def q_to_page(value: float) -> float:
        return 205.0 + (value - q_low) / (q_high - q_low) * 300.0

    for value in q_values:
        yp = q_to_page(float(value))
        drawing.add(Line(right_x - 4.0, yp, right_x + 4.0, yp,
                         strokeColor=HexColor("#9AA2AD"), strokeWidth=0.65, strokeOpacity=0.55))

    labels = ["Q₁(y*)", "Q₂(y*)", "Qᵣ(y*)", "Qₛ(y*)", "Q_R₋₁(y*)", "Q_R(y*)"]
    point_pages = [q_to_page(float(q_values[index])) for index in selected_local]
    label_pages = [point_pages[0]]
    for point_page in point_pages[1:]:
        label_pages.append(min(point_page, label_pages[-1] - 27.0))
    for local_index, color, label, label_page in zip(selected_local, highlight_colors, labels, label_pages):
        value = float(q_values[local_index])
        yp = q_to_page(value)
        drawing.add(Circle(x_star_page, float(data_to_page(np.array([actual_y_star]), np.array([value]), left_box, xlim, ylim)[1][0]),
                           5.2, fillColor=color, strokeColor=color))
        drawing.add(Line(x_star_page + 6, float(data_to_page(np.array([actual_y_star]), np.array([value]), left_box, xlim, ylim)[1][0]),
                         right_x - 6, yp, strokeColor=color, strokeWidth=1.0,
                         strokeDashArray=[4, 4], strokeOpacity=0.85))
        drawing.add(Circle(right_x, yp, 5.2, fillColor=color, strokeColor=color))
        drawing.add(Line(right_x + 6, yp, right_x + 17, label_page,
                         strokeColor=axis_color, strokeWidth=0.65))
        boxed_text(drawing, right_x + 22, label_page - 5, label,
                   "PaperItalic", 15, axis_color, padding=2.5)

    q_mean = float(q_values.mean())
    q_mean_page = q_to_page(q_mean)
    drawing.add(Circle(right_x, q_mean_page, 5.0, fillColor=axis_color, strokeColor=axis_color))
    drawing.add(Line(right_x + 7, q_mean_page, right_x + 74, q_mean_page,
                     strokeColor=axis_color, strokeWidth=1.0))
    boxed_text(drawing, right_x - 104, q_mean_page + 14, "Q̄_R*",
               "PaperItalic", 15, axis_color, padding=3.0)

    branch_y = 80.0
    drawing.add(Line(right_x, right_bottom - 2, right_x - 68, branch_y,
                     strokeColor=axis_color, strokeWidth=1.1))
    drawing.add(Line(right_x, right_bottom - 2, right_x, branch_y,
                     strokeColor=axis_color, strokeWidth=1.1))
    drawing.add(Line(right_x, right_bottom - 2, right_x + 68, branch_y,
                     strokeColor=axis_color, strokeWidth=1.1))
    arrow(drawing, right_x - 68, branch_y, right_x - 78, branch_y - 28, axis_color, 1.1)
    arrow(drawing, right_x, branch_y, right_x, branch_y - 31, axis_color, 1.1)
    arrow(drawing, right_x + 68, branch_y, right_x + 78, branch_y - 28, axis_color, 1.1)
    drawing.add(String(right_x - 113, 22, "Student t", fontName="Times-Italic", fontSize=13, fillColor=axis_color))
    drawing.add(String(right_x - 34, 22, "percentile", fontName="Times-Italic", fontSize=13, fillColor=axis_color))
    drawing.add(String(right_x + 55, 22, "bootstrap-t", fontName="Times-Italic", fontSize=13, fillColor=axis_color))

    formula = "Qᵣ(y) = ∑ᵢ₌₁ⁿ Pᵣ,ᵢ Kₕ{y − g(θᵣ,ᵢ)}"
    drawing.add(String(76, 58, formula, fontName="PaperItalic", fontSize=16, fillColor=axis_color))

    output_dir.mkdir(parents=True, exist_ok=True)
    output_dir = output_dir.resolve()
    background = Drawing(WIDTH, HEIGHT)
    for index, item in enumerate(drawing.contents):
        if isinstance(item, String):
            continue
        if isinstance(item, Rect) and index != 0:
            continue
        background.add(copy.deepcopy(item))
    background_pdf = output_dir / "real_data_outer_sample_workflow_background.pdf"
    renderPDF.drawToFile(background, str(background_pdf))
    tex_path = write_latex_overlay(output_dir, x_star_page, label_pages, q_mean_page)
    pdf_path, svg_path, png_path = compile_latex_figure(output_dir, tex_path)

    selected_set = set(int(value) for value in selected_local)
    figure_csv = output_dir / "real_data_outer_sample_workflow_data.csv"
    with figure_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "curve_id", "y", "source_density", "background_smoothed_density",
            "R_display_density", "in_R_sample", "highlighted", "display_density_at_y_star",
        ])
        background_at_vertices = np.vstack([
            np.interp(y, dense_y, curve) for curve in dense_pool
        ])
        r_display_at_vertices = np.vstack([
            np.interp(y, dense_y, curve) for curve in r_dense
        ])
        for curve_index, curve in enumerate(pool):
            for y_index, (yi, density) in enumerate(zip(y, curve)):
                writer.writerow([
                    curve_index + 1,
                    f"{yi:.12g}",
                    f"{density:.12g}",
                    f"{background_at_vertices[curve_index, y_index]:.12g}",
                    f"{r_display_at_vertices[curve_index, y_index]:.12g}" if curve_index < r_count else "",
                    int(curve_index < r_count),
                    int(curve_index in selected_set),
                    f"{q_values[curve_index]:.12g}" if curve_index < r_count else "",
                ])

    metadata = {
        "source_pool": pool_path.name,
        "source_pool_sha256": hashlib.sha256(pool_path.read_bytes()).hexdigest(),
        "source_mat": source_mat.name if source_mat else None,
        "source_mat_sha256": hashlib.sha256(source_mat.read_bytes()).hexdigest() if source_mat else None,
        "source_mat_variable": "pilot",
        "source_config": config_path.name,
        "model": "linear Gaussian benchmark",
        "P": int(pool.shape[0]),
        "n": int(config["n"]),
        "h": float(config["h"]),
        "R": int(r_count),
        "B": int(config["linear_B"]),
        "y_star": actual_y_star,
        "R_curve_ids": list(range(1, r_count + 1)),
        "highlighted_curve_ids": [int(index + 1) for index in selected_local],
        "highlighted_q_values": [float(q_values[index]) for index in selected_local],
        "outer_sample_mean": q_mean,
        "display_rank_amplitude": display_rank_amplitude,
        "display_shape_retention": display_shape_retention,
        "selection_rule": "The largest two, two interior ranks, and the smallest two Q_r(y*) values within the first R formal curves.",
        "display_note": "Actual TC pilot ordinates; no rank offset or curve amplitude rescaling.",
        "smoothing": "Cubic Hermite interpolation for rendering only; no smoothing of source ordinates.",
        "final_typography": "XeLaTeX overlay using Latin Modern Roman text and Computer Modern mathematics, matching the source paper.",
        "formula_source": "Q_r(y)=sum_{i=1}^n P_{r,i} K_h{y-g(bold theta_{r,i})}, copied from the paper and typeset directly by LaTeX.",
    }
    (output_dir / "real_data_outer_sample_workflow_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    return metadata



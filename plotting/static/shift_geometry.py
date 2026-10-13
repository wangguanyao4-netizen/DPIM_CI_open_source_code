"""Vector illustration of first and last common random shifts on one response map."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
HELPER = ROOT / 'method_geometry.py'
spec = importlib.util.spec_from_file_location('vector_helpers', HELPER)
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)
np, plt = v.np, v.plt
from matplotlib.patches import FancyArrowPatch, Polygon
from matplotlib.path import Path as MPath
from matplotlib.text import Text
from mpl_toolkits.mplot3d import proj3d

COLORS = ['#285e87', '#a06034']
FILLS = ['#e4edf4', '#f2e8df']
SHIFTS = [np.array([.16, .08]), np.array([-.15, .18])]
LABELS = ['1', 'R']
MARKERS = ['o', '^']
SELECTED = 7
H = .17
ARROW_WIDTH = .78
ARROW_HEAD = 10.5
ARROW_ALPHA = .90


def response(theta):
    return (theta[..., 0] + .55 * np.sin(2 * np.pi * theta[..., 1])
            + .15 * theta[..., 1])


def periodic_path(p, shift):
    """Split the unwrapped line at every unit-square boundary crossing."""
    times = [0., 1.]
    for j in range(2):
        if shift[j] != 0:
            for boundary in [0., 1.]:
                t = (boundary - p[j]) / shift[j]
                if 1e-10 < t < 1 - 1e-10:
                    times.append(float(t))
    times = sorted(set(times))
    segments = []
    for a, b in zip(times[:-1], times[1:]):
        cell = np.floor(p + (.5 * (a + b)) * shift)
        segments.append((p + a * shift - cell, p + b * shift - cell))
    return segments


def sample_bezier_path(path, samples_per_segment=180):
    """Sample every cubic segment in a Matplotlib path."""
    vertices=np.asarray(path.vertices)
    if (len(vertices)-1) % 3 != 0:
        raise ValueError('Expected one or more cubic Bezier segments')
    sampled=[]
    for start in range(0, len(vertices)-1, 3):
        p0, p1, p2, p3=vertices[start:start+4]
        tt=np.linspace(.04, .96, samples_per_segment)[:, None]
        sampled.append((1-tt)**3*p0 + 3*(1-tt)**2*tt*p1
                       + 3*(1-tt)*tt**2*p2 + tt**3*p3)
    return np.vstack(sampled)


def build():
    fig = plt.figure(figsize=(12.2, 6.25))
    paper = fig.add_axes([0, 0, 1, 1])
    paper.set(xlim=(0, 1), ylim=(0, 1)); paper.axis('off')
    text = lambda x, y, s, size=12, **kw: v.text(paper, x, y, s, size=size, **kw)
    cells = v.periodic_cells()
    weights = np.array([v.area(c) for c in cells])
    points = [(v.POINTS + s) % 1 for s in SHIFTS]
    outputs = [response(u) for u in points]
    yy = np.linspace(-.8, 2, 1200)
    parts = [weights[:, None] * np.exp(-.5 * ((yy - y[:, None]) / H)**2)
             / (H * np.sqrt(2 * np.pi)) for y in outputs]
    ymax = max(p.sum(axis=0).max() for p in parts) * 1.17
    text(.040, .955, r'$(a)$', 12, ha='left')
    text(.282, .955, r'$(b)$', 12, ha='left')
    text(.755, .955, r'$(c)$', 12, ha='left')
    text(.143, .955, 'Random shifts', 12)
    text(.505, .955, r'Deterministic response surface, $Y=g(\boldsymbol{\theta})$', 12)
    text(.866, .955, 'DPIM density estimates', 12)
    rows = [.575, .120]
    input_axes=[]
    density_axes=[]
    shift_arrow_data=[]
    for k, (label, row) in enumerate(zip(LABELS, rows)):
        ax = fig.add_axes([.055, row, .190, .304])
        ax.set(xlim=(0, 1), ylim=(0, 1), aspect='equal',
               xticks=[0, 1], yticks=[0, 1])
        ax.tick_params(direction='out', length=2.5, width=.6, labelsize=8, pad=2)
        for spine in ax.spines.values():
            spine.set_linewidth(.7)
        ax.set_xlabel(r'$\theta_1$', fontsize=10, labelpad=-2)
        ax.set_ylabel(r'$\theta_2$', fontsize=10, labelpad=-5)
        shifted_areas=[]
        for i, cell in enumerate(cells):
            fragments=v.pieces(cell, SHIFTS[k])
            shifted_areas.append(sum(v.area(poly) for poly in fragments))
            for poly in fragments:
                ax.add_patch(Polygon(
                    poly, closed=True,
                    facecolor='#e7e7e7' if i == SELECTED else '#fafafa',
                    edgecolor='#858585', linewidth=.52, zorder=0))
        assert np.allclose(shifted_areas, weights, atol=1e-12, rtol=0)
        for i, original in enumerate(v.POINTS):
            segments=periodic_path(original, SHIFTS[k])
            for a, b in segments:
                ax.add_patch(FancyArrowPatch(
                    a, b, arrowstyle='-|>', color=COLORS[k],
                    linewidth=ARROW_WIDTH, mutation_scale=ARROW_HEAD,
                    shrinkA=4, shrinkB=6, alpha=ARROW_ALPHA, zorder=3))
                shift_arrow_data.append(
                    (ax, np.asarray(a), np.asarray(b), k, i))
            assert np.allclose(
                sum((b-a for a, b in segments)), SHIFTS[k])
        ax.scatter(v.POINTS[:, 0], v.POINTS[:, 1], s=18, facecolors='white',
                   edgecolors='#666666', linewidths=.7, zorder=4)
        ax.scatter(*v.POINTS[SELECTED], s=39, facecolors='white',
                   edgecolors='black', linewidths=.75, zorder=5)
        ax.scatter(points[k][:, 0], points[k][:, 1], s=21, c=COLORS[k],
                   marker=MARKERS[k], zorder=5)
        ax.scatter(*points[k][SELECTED], s=45, c=COLORS[k], marker=MARKERS[k],
                   edgecolors='black', linewidths=.7, zorder=6)
        v.text(ax, .02, 1.035, rf'$r={label}$', size=11,
               ha='left', va='bottom', clip_on=False)
        midpoint=(v.POINTS[SELECTED]+.5*SHIFTS[k]) % 1
        normal=np.array([-SHIFTS[k][1], SHIFTS[k][0]])
        normal=normal/np.linalg.norm(normal)
        delta_pos=midpoint + [1, -1][k]*[.125, .120][k]*normal
        v.text(ax, *delta_pos, rf'$\Delta_{{{label}}}$', size=11)
        den = fig.add_axes([.770, row+.033, .200, .225])
        den.set(xlim=(-.8, 2), ylim=(0, ymax), xticks=[], yticks=[])
        den.spines[['right', 'top']].set_visible(False)
        den.spines[['bottom', 'left']].set_linewidth(.7)
        den.fill_between(yy, 0, parts[k][SELECTED], color=FILLS[k], zorder=1)
        den.plot(yy, parts[k][SELECTED], color=COLORS[k], lw=.72,
                 linestyle=(0, (3, 2)), zorder=2)
        den.plot(yy, parts[k].sum(axis=0), color=COLORS[k], lw=1.35)
        den.scatter(outputs[k], np.zeros(len(weights)), c=COLORS[k], s=8,
                    marker=MARKERS[k], clip_on=False, zorder=3)
        den.scatter(outputs[k][SELECTED], 0, c=COLORS[k], s=40, marker=MARKERS[k],
                    edgecolors='black', linewidths=.6, clip_on=False, zorder=4)
        v.text(den, 1.91, .87*ymax, rf'$Q_{{{label}}}(y)$', size=12, ha='right')
        v.text(den, 1.98, -.11*ymax, '$y$', size=10, ha='right', va='top')
        input_axes.append(ax)
        density_axes.append(den)
    text(.150, .505, r'$\vdots$', 18)
    text(.870, .505, r'$\vdots$', 18)

    surf = fig.add_axes([.270, .105, .465, .775], projection='3d')
    surf.set_facecolor('white')
    a, b = np.meshgrid(np.linspace(0, 1, 35), np.linspace(0, 1, 35))
    z = response(np.stack([a, b], axis=-1))
    surf.plot_surface(a, b, z, rstride=2, cstride=2, color='#dfe1e2', alpha=.88,
                      edgecolor='#a5aaac', linewidth=.22, antialiased=True, shade=True)
    surf.contour(a, b, z, levels=7, zdir='z', offset=-.65,
                 colors='#b5b8ba', linewidths=.45)
    surf.set_proj_type('persp', focal_length=.82)
    surf.view_init(elev=29, azim=-132)
    surf.set(xlim=(0, 1), ylim=(0, 1), zlim=(-.65, 1.8),
             xticks=[0, .5, 1], yticks=[0, .5, 1], zticks=[-.5, .5, 1.5])
    surf.set_xticklabels(['0', '0.5', '1'])
    surf.set_yticklabels(['0', '0.5', '1'])
    surf.set_zticklabels(['-0.5', '0.5', '1.5'])
    surf.tick_params(labelsize=8, pad=0, width=.6)
    surf.set_xlabel(r'$\theta_1$', fontsize=11, labelpad=4)
    surf.set_ylabel(r'$\theta_2$', fontsize=11, labelpad=4)
    surf.set_zlabel('$Y$', fontsize=11, labelpad=5)
    surf.set_box_aspect((1.22, 1.0, .92)); surf.grid(True)
    for axis in [surf.xaxis, surf.yaxis, surf.zaxis]:
        axis.pane.fill=True
        axis.pane.set_facecolor((.975, .975, .975, 1))
        axis.pane.set_edgecolor('#777777')
        axis._axinfo['grid']['color'] = (.80, .80, .80, .65)
        axis._axinfo['grid']['linewidth'] = .35
    # Reproject the exact surface coordinates as a vector overlay so markers
    # remain legible without moving them away from Y = g(theta).
    fig.canvas.draw()
    overlay = fig.add_axes([0, 0, 1, 1], facecolor='none')
    overlay.set(xlim=(0, 1), ylim=(0, 1)); overlay.axis('off')
    to_fig = fig.transFigure.inverted()
    projected_sets=[]
    arrow_paths=[]
    arrow_targets=[]
    for k in range(2):
        projected=[]
        for p, y in zip(points[k], outputs[k]):
            sx, sy, _ = proj3d.proj_transform(p[0], p[1], y, surf.get_proj())
            projected.append(to_fig.transform(surf.transData.transform((sx, sy))))
        projected=np.asarray(projected)
        projected_sets.append(projected)
        start=to_fig.transform(input_axes[k].transData.transform(points[k][SELECTED]))
        middle=projected[SELECTED]
        end=to_fig.transform(density_axes[k].transData.transform((outputs[k][SELECTED], 0)))
        if k == 0:
            routes=[
                [start,
                 (start[0]+.060, start[1]+.120),
                 (middle[0]-.100, middle[1]+.100),
                 middle],
                [middle,
                 (.735, middle[1]),
                 (.735, end[1]-.055),
                 (.770, end[1]-.055),
                 (.820, end[1]-.055),
                 (end[0], end[1]-.055),
                 end],
            ]
        else:
            routes=[
                [start,
                 (start[0]+.030, start[1]+.160),
                 (middle[0]-.100, middle[1]+.100),
                 middle],
                [middle,
                 (.735, middle[1]),
                 (.735, end[1]-.055),
                 (.770, end[1]-.055),
                 (.820, end[1]-.055),
                 (end[0], end[1]-.055),
                 end],
            ]
        for vertices in routes:
            path=MPath(vertices, [MPath.MOVETO] + [MPath.CURVE4]*(len(vertices)-1))
            arrow_paths.append(path)
            arrow_targets.append((np.asarray(vertices[0]), np.asarray(vertices[-1])))
            overlay.add_patch(FancyArrowPatch(
                path=path, transform=overlay.transData, arrowstyle='-|>',
                mutation_scale=ARROW_HEAD, shrinkA=7, shrinkB=17,
                linewidth=ARROW_WIDTH, color=COLORS[k],
                alpha=ARROW_ALPHA, zorder=12,
                joinstyle='round', capstyle='round'))
        overlay.scatter(projected[:, 0], projected[:, 1], s=25,
                        c=COLORS[k], marker=MARKERS[k], edgecolors='black',
                        linewidths=.35, zorder=20)
        overlay.scatter(projected[SELECTED, 0], projected[SELECTED, 1], s=62,
                        c=COLORS[k], marker=MARKERS[k], edgecolors='black',
                        linewidths=.75, zorder=21)
    paper.scatter([.395], [.048], s=22, c=COLORS[0], marker=MARKERS[0])
    text(.408, .048, '$r=1$', 9.5, ha='left')
    paper.scatter([.475], [.048], s=24, c=COLORS[1], marker=MARKERS[1])
    text(.488, .048, '$r=R$', 9.5, ha='left')
    paper.plot([.645, .675], [.048, .048], color='black', lw=1.15)
    text(.683, .048, 'density', 9.2, ha='left')
    paper.fill_between([.780, .810], [.037, .037], [.059, .059], color='#e8ecef')
    paper.plot([.780, .810], [.048, .048], color='black', lw=.72,
               linestyle=(0, (3, 2)))
    text(.818, .048, 'selected kernel', 9.2, ha='left')

    # Reject any correspondence curve that crosses text or an unrelated marker.
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    arrow_text_overlaps=[]
    for arrow_index, path in enumerate(arrow_paths):
        samples=sample_bezier_path(path, samples_per_segment=240)
        for artist in fig.findobj(match=Text):
            label=artist.get_text().strip()
            if not label or not artist.get_visible():
                continue
            bbox=artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())
            pad_x, pad_y=.0025, .0035
            inside=((samples[:, 0] >= bbox.x0-pad_x)
                    & (samples[:, 0] <= bbox.x1+pad_x)
                    & (samples[:, 1] >= bbox.y0-pad_y)
                    & (samples[:, 1] <= bbox.y1+pad_y))
            if np.any(inside):
                hit=samples[np.flatnonzero(inside)[0]]
                arrow_text_overlaps.append(
                    dict(arrow=arrow_index, text=label,
                         text_bounds=[float(bbox.x0), float(bbox.y0),
                                      float(bbox.x1), float(bbox.y1)],
                         curve_point=[float(hit[0]), float(hit[1])]))
    assert not arrow_text_overlaps, arrow_text_overlaps

    shift_arrow_text_overlaps=[]
    shift_arrow_marker_overlaps=[]
    shift_samples=[]
    for arrow_index, (ax, a, b, k, point_index) in enumerate(shift_arrow_data):
        a_fig=to_fig.transform(ax.transData.transform(a))
        b_fig=to_fig.transform(ax.transData.transform(b))
        ts=np.linspace(.08, .92, 180)[:, None]
        samples=(1-ts)*a_fig + ts*b_fig
        shift_samples.append(samples)
        for artist in fig.findobj(match=Text):
            label=artist.get_text().strip()
            if not label or not artist.get_visible():
                continue
            bbox=artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())
            inside=((samples[:, 0] >= bbox.x0-.0025)
                    & (samples[:, 0] <= bbox.x1+.0025)
                    & (samples[:, 1] >= bbox.y0-.0035)
                    & (samples[:, 1] <= bbox.y1+.0035))
            if np.any(inside):
                shift_arrow_text_overlaps.append(
                    dict(arrow=arrow_index, text=label))
        local_markers=np.vstack([v.POINTS, points[k]])
        for marker_index, center in enumerate(local_markers):
            if (np.linalg.norm(center-a) < 1e-8
                    or np.linalg.norm(center-b) < 1e-8):
                continue
            center_fig=to_fig.transform(ax.transData.transform(center))
            if np.min(np.linalg.norm(samples-center_fig, axis=1)) < .0065:
                shift_arrow_marker_overlaps.append(
                    dict(arrow=arrow_index, marker=marker_index))
    assert not shift_arrow_text_overlaps, shift_arrow_text_overlaps
    assert not shift_arrow_marker_overlaps, shift_arrow_marker_overlaps

    visible_text=[]
    for artist in fig.findobj(match=Text):
        label=artist.get_text().strip()
        if not label or not artist.get_visible():
            continue
        bbox=artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())
        visible_text.append((label, bbox))
    text_text_overlaps=[]
    for left in range(len(visible_text)):
        label_a, box_a=visible_text[left]
        for right in range(left+1, len(visible_text)):
            label_b, box_b=visible_text[right]
            overlap_x=min(box_a.x1, box_b.x1)-max(box_a.x0, box_b.x0)
            overlap_y=min(box_a.y1, box_b.y1)-max(box_a.y0, box_b.y0)
            if overlap_x > .0015 and overlap_y > .0015:
                text_text_overlaps.append(dict(first=label_a, second=label_b))
    assert not text_text_overlaps, text_text_overlaps

    marker_centers=[]
    for k in range(2):
        for p in v.POINTS:
            marker_centers.append(to_fig.transform(input_axes[k].transData.transform(p)))
        for p in points[k]:
            marker_centers.append(to_fig.transform(input_axes[k].transData.transform(p)))
        marker_centers.extend(projected_sets[k])
        for y in outputs[k]:
            marker_centers.append(
                to_fig.transform(density_axes[k].transData.transform((y, 0))))
    marker_centers=np.asarray(marker_centers)
    arrow_marker_overlaps=[]
    for arrow_index, path in enumerate(arrow_paths):
        samples=sample_bezier_path(path, samples_per_segment=220)
        allowed=arrow_targets[arrow_index]
        for marker_index, center in enumerate(marker_centers):
            if any(np.linalg.norm(center-endpoint) < 1e-8 for endpoint in allowed):
                continue
            if np.min(np.linalg.norm(samples-center, axis=1)) < .0065:
                arrow_marker_overlaps.append(
                    dict(arrow=arrow_index, marker=marker_index))
    assert not arrow_marker_overlaps, arrow_marker_overlaps

    arrow_density_curve_overlaps=[]
    for k, arrow_index in enumerate([1, 3]):
        samples=sample_bezier_path(arrow_paths[arrow_index], samples_per_segment=260)
        density_curves=[]
        for values in [parts[k].sum(axis=0), parts[k][SELECTED]]:
            density_curves.append(np.asarray([
                to_fig.transform(density_axes[k].transData.transform((x, value)))
                for x, value in zip(yy, values)
            ]))
        for curve_index, curve in enumerate(density_curves):
            minimum=min(
                np.min(np.linalg.norm(curve-point, axis=1)) for point in samples)
            if minimum < .0035:
                arrow_density_curve_overlaps.append(
                    dict(arrow=arrow_index, density_curve=curve_index,
                         minimum_distance=float(minimum)))
    assert not arrow_density_curve_overlaps, arrow_density_curve_overlaps

    checks=[]
    for k in range(2):
        wide = np.linspace(-4, 5, 6000)
        integrals = np.trapezoid(weights[:, None] * np.exp(-.5*((wide-outputs[k][:, None])/H)**2)
                                / (H*np.sqrt(2*np.pi)), wide, axis=1)
        assert np.allclose(integrals, weights, atol=1e-12, rtol=0)
        assert abs(integrals.sum()-1) < 1e-12
        checks.append(dict(r=LABELS[k],shift=SHIFTS[k].tolist(),points=points[k].tolist(),
                           responses=outputs[k].tolist(),density_integral=float(integrals.sum())))
    return fig, dict(repetitions=checks, weights=weights.tolist(), h=H,
                     response='g(theta) = theta1 + 0.55 sin(2 pi theta2) + 0.15 theta2',
                     interpretation='Illustrative geometry, not numerical-example data',
                     common_shift_all_points=True, intermediate_repetitions='vertical ellipsis',
                     displayed_shift_arrows='all representative points',
                     shifted_periodic_cells_drawn=True,
                     shifted_cell_areas_preserved=True,
                     response_markers_projected_at_exact_surface_coordinates=True,
                     selected_point_arrows='theta_i^(r) to g(theta_i^(r)) to Y_i^(r)',
                     curved_arrowheads_visible=True,
                     arrow_text_overlaps=arrow_text_overlaps,
                     arrow_marker_overlaps=arrow_marker_overlaps,
                     arrow_density_curve_overlaps=arrow_density_curve_overlaps,
                     shift_arrow_text_overlaps=shift_arrow_text_overlaps,
                     shift_arrow_marker_overlaps=shift_arrow_marker_overlaps,
                     text_text_overlaps=text_text_overlaps,
                     shaded_component='selected probability-weighted kernel')





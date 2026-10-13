"""Maintainer-only extraction of approved, lightweight publication evidence.

Requires an explicit historical archive. Never runs simulations or changes it.
Ordinary users should run replay.py instead; it needs only packaged CSV files.
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1]
M30 = "coverage_M30000_manuscript_20260927/publication_data/"
N500 = "nonlinear_reference_500M_20260927/publication_data/"
REFRESH = "paper_refresh_30000_20260926/"
SELECTED = "manuscript_selected_nonlinear500M_20260927/publication_data/"
REPORT = "reports/manuscript_structure_tables_20261012/"
TC = "tc_recalculation_20260922/plotting/"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def extract_nodes(path, names):
    text = path.read_text(encoding="utf-8-sig")
    lines = text.splitlines(keepends=True)
    result = []
    found = set()
    for node in ast.parse(text).body:
        name = getattr(node, "name", None)
        if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
        if name in names:
            result.append("".join(lines[node.lineno - 1:node.end_lineno]))
            found.add(name)
    if found != set(names):
        raise ValueError(f"Missing source nodes: {set(names) - found}")
    return "\n\n".join(result) + "\n"


def prepare(archive, include_sources=False):
    records = []

    def inputs(paths, selection):
        return [{"archive_relative_path": p, "sha256": sha(archive / p),
                 "selection": selection} for p in paths]

    def record(dest, paths, selection):
        p = PACKAGE / dest
        records.append({"path": dest, "sha256": sha(p), "bytes": p.stat().st_size,
                        "inputs": inputs(paths, selection)})

    def copy(source, dest):
        target = PACKAGE / dest
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(archive / source, target)
        record(dest, [source], "byte-for-byte copy of approved source")

    def csv_out(name, rows, paths, selection):
        dest = "data/publication/" + name
        write_csv(PACKAGE / dest, rows)
        record(dest, paths, selection)

    primary_path = M30 + "nested_M_pointwise_coverage.csv"
    final_path = N500 + "nested_M_pointwise_coverage.csv"
    rows = [r for r in read_csv(archive / primary_path)
            if r["model"] in ("linear", "plate") and int(r["M"]) == 30000]
    nonlinear = read_csv(archive / final_path)
    assert len(rows) == 300 and len(nonlinear) == 150
    assert all(r["model"] == "nonlinear" and int(r["M"]) == 50000
               and int(r["P"]) == 50000 and int(r["reference_N"]) == 500000000
               for r in nonlinear)
    for r in nonlinear:
        r["dataset"] = "primary"
    csv_out("coverage.csv", rows + nonlinear, [primary_path, final_path],
            "linear/plate M30000 only; discard ALL old nonlinear rows; append nonlinear P50000 M50000 reference_N500000000")
    comparisons = M30 + "nested_M_all_comparisons.csv"
    controls = [r for r in read_csv(archive / comparisons)
                if r["dataset"] == "euler_pool_control" and int(r["M"]) == 30000]
    assert len(controls) == 150
    csv_out("euler_pool_control.csv", controls, [comparisons],
            "dataset=euler_pool_control, M=30000; separate from primary evidence")

    summary_paths = [M30 + "coverage_summary.csv", N500 + "coverage_summary.csv"]
    summary = [r for r in read_csv(archive / summary_paths[0])
               if r["model"] in ("linear", "plate") and int(r["M"]) == 30000]
    for r in read_csv(archive / summary_paths[1]):
        if int(r["reference_N"]) == 500000000:
            r.update(model="nonlinear", M="50000", P="50000", dataset="primary")
            r["mean_absolute_nominal_error"] = r.pop("mean_abs_nominal_error")
            summary.append(r)
    csv_out("coverage_summary.csv", summary, summary_paths,
            "linear/plate M30000 including labelled Euler control; nonlinear reference_N500000000 only")

    for model in ("linear", "plate"):
        copy(REFRESH + f"publication_data/reference_{model}/reference_density.csv",
             f"data/publication/reference_{model}.csv")
    for name in ("reference_density", "reference_coverage", "reference_n_density"):
        copy(N500 + name + ".csv", "data/publication/" + name.replace("reference", "reference_nonlinear", 1) + ".csv")
    copy(SELECTED + "n_uncentered_coverage_results.csv", "data/publication/internal_n.csv")
    copy(REFRESH + "results/plate_mc_rqmc_recommended_R_intervals.csv", "data/publication/intervals_plate.csv")
    # Final common CSV keeps the independently generated endpoints; only the
    # nonlinear reference was updated. Check the Euler block was not changed.
    old = read_csv(archive / (REFRESH + "results/mc_rqmc_recommended_R_intervals.csv"))
    final = read_csv(archive / (SELECTED + "mc_rqmc_recommended_R_intervals.csv"))
    csv_out("intervals_linear_nonlinear.csv",
            [r for r in old if r["model"] == "linear"] + [r for r in final if r["model"] == "nonlinear"],
            [REFRESH + "results/mc_rqmc_recommended_R_intervals.csv", SELECTED + "mc_rqmc_recommended_R_intervals.csv"],
            "Euler from original full-precision CSV; nonlinear from final500M CSV; retains independent endpoints, no re-estimation")

    rp = REFRESH + "results/recommended_R_before_intervals.csv"
    np = "nonlinear_P50000_M50000_20260927/publication_data/recommended_R.csv"
    fields = ["model", "scheme", "method", "n", "B", "R_recommended", "Lambda", "P"]
    recommendations = [{k: r[k] for k in fields} for r in read_csv(archive / rp)
                       if not (r["model"] == "nonlinear" and r["scheme"] == "rqmc")]
    for r in read_csv(archive / np):
        recommendations.append(dict(model="nonlinear", scheme="rqmc",
            method=r["method"].replace("t distribution", "Student t"), n="768", B="399",
            R_recommended=r["recommended_R"], Lambda=r["Lambda"], P=r["P"]))
    csv_out("recommended_R.csv", recommendations, [rp, np],
            "replace nonlinear RQMC rows with P50000 recommendation; MC remains P30000; NOT new independent intervals")
    for name in ("saved_trial_intervals", "table_interval_values"):
        copy(REPORT + "data/" + name + ".csv", "data/publication/" + name + ".csv")
    for name in ("summary", "summary_paper_signed_difference"):
        copy("DPIM_CI_paper_used_only_20260804/03_results/finite_B/" + name + ".csv",
             "data/publication/finite_B_" + name + ".csv")

    # A figure-only 256-by-201 slice, not a pool for coverage recounting.
    from scipy.io import loadmat
    import numpy as np
    pool = "tc_recalculation_20260922/results/paper_scale_original_grid_20260922/linear_rqmc_n768/pilot_complete.mat"
    q = loadmat(archive / pool)["pilot"][:256]
    assert q.shape == (256, 201)
    target = PACKAGE / "data/publication/outer_workflow_curves.csv"
    np.savetxt(target, q, delimiter=",", fmt="%.17g")
    record("data/publication/outer_workflow_curves.csv", [pool], "pilot[0:256, :] verbatim float64 ordinates; no random draws")

    expected = REPORT + "compile/Wang_DPIM_Unified_Framework_Euler_Corrected/DPIM_CI_full_integrated_figures/formal_current/"
    pdfs = sorted((archive / expected).glob("*.pdf"))
    assert len(pdfs) == 17
    for p in pdfs:
        copy(expected + p.name, "data/publication/expected/" + p.name)

    if include_sources:
        mappings = {
            "plotting/matlab/plot_nested_M_paper_style.m": "reference_expansion_20260925/plot_nested_M_paper_style.m",
            "plotting/matlab/plot_intervals_pdf_label_only.m": REPORT + "figure_labels/interval_density/plot_intervals_pdf_label_only.m",
            "plotting/plot_linear_plate_real_data_mechanisms.py": TC + "plot_linear_plate_real_data_mechanisms.py",
            "plotting/plot_nonlinear_real_data_mechanism.py": TC + "plot_nonlinear_real_data_mechanism.py",
            "plotting/static/draw_transport.py": "deliverables/workflow_figure_revision_20260909/draw_transport.py",
            "plotting/static/euler_model.tex": "DPIM_CI_paper_used_only_20260804/02_figure_sources/mechanical_models/linear_euler_bernoulli_beam_model.tex",
            "plotting/static/draw_plate.py": "DPIM_CI_paper_used_only_20260804/02_figure_sources/mechanical_models/draw_plate_sfem_model.py",
            "plotting/static/outer_workflow.py": TC + "plot_real_data_outer_sample_workflow.py",
            "plotting/static/outer_overlay.tex": "deliverables/euler_analytic_full_manuscript_20260923/review_figures/workflow/real_data_outer_sample_workflow_latex.tex",
            "plotting/static/method_geometry.py": "figures/dpim_method_drafts_20260928/draw_method_figures.py",
            "plotting/static/framework_geometry.py": "figures/dpim_method_drafts_v2_20260928/draw_method_figures.py",
            "plotting/static/shift_geometry.py": "figures/dpim_method_drafts_v6_20260928/draw_first_last_shift.py",
            "plotting/static/framework_layout.py": "reports/framework_figure_fix_20260928/original_layout.py",
            "plotting/static/shift_layout.py": "reports/mechanism_arrow_tips_20260928/workflow.py",
        }
        for dest, source in mappings.items():
            copy(source, dest)
        selections = [
            ("plotting/internal_n.py", TC + "update_n_coverage_revision_20260630.py",
             ["METHODS", "METHOD_LABEL", "METHOD_COLOR", "N_VALUES", "configure_plots", "plot_n_coverage"]),
            ("plotting/labels.py", "deliverables/euler_analytic_full_manuscript_20260923/review_figures/redraw_python_labels.py", ["save_labels"]),
            ("plotting/spacing.py", "coverage_M30000_manuscript_20260927/render_mechanisms.py", ["spaced_save", "fix_annotation_separators"]),
            ("plotting/current_labels.py", REPORT + "redraw_mechanism_labels.py", ["geometry", "change_labels"]),
        ]
        for dest, source, names in selections:
            (PACKAGE / dest).write_text(extract_nodes(archive / source, names), encoding="utf-8")
            record(dest, [source], "mechanical AST extraction: " + ", ".join(names))
    provenance = PACKAGE / "provenance/plot_sources.json"
    provenance.parent.mkdir(exist_ok=True)
    old_records = json.loads(provenance.read_text()).get("files", []) if provenance.exists() else []
    by_path = {r["path"]: r for r in old_records}
    by_path.update({r["path"]: r for r in records})
    provenance.write_text(json.dumps({
        "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
        "path_policy": "path is package-relative; archive_relative_path is historical provenance, NEVER a runtime dependency",
        "data_scope": "publication summaries, density ordinates, archived interval endpoints and a 256-curve illustration slice; NOT raw pools or M-trial archives",
        "generation_command": "python plotting/prepare_publication.py --archive-root <approved-local-archive>",
        "new_simulations": 0, "files": sorted(by_path.values(), key=lambda r: r["path"]),
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Staged {len(records)} items; no experiments executed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-root", required=True, type=Path)
    parser.add_argument("--bootstrap-sources", action="store_true",
                        help="One-time maintainer extraction; refuses to overwrite adapted sources")
    args = parser.parse_args()
    if args.bootstrap_sources and (PACKAGE / "plotting/internal_n.py").exists():
        parser.error("Source bootstrap already exists; do not overwrite adapted code")
    prepare(args.archive_root.resolve(), args.bootstrap_sources)

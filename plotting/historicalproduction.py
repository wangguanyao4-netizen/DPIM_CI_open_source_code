"""Explicit adapters for archived plotters. Importing this file runs nothing.

These commands have NOT been executed during public-release assembly.
No raw pools, experiments, reference estimation, or manuscript are invoked.
"""
from __future__ import annotations

import argparse
import importlib
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "data/publication"


def python_figures(data, output):
    import matplotlib
    matplotlib.use("Agg")
    import pandas as pd
    from matplotlib.figure import Figure
    import labels
    import current_labels
    import spacing
    import internal_n as ns
    import plot_linear_plate_real_data_mechanisms as lp
    import plot_nonlinear_real_data_mechanism as nl

    cov = pd.read_csv(data / "coverage.csv")
    native = Figure.savefig
    current = {}

    def final_save(fig, *args, **kwargs):
        if not getattr(fig, "_current_labels", False) and len(fig.axes) == 3:
            current_labels.change_labels(fig)
            fig._current_labels = True
        return native(fig, *args, **kwargs)

    labels.save_original = spacing.spaced_save(final_save)

    def dispatch(fig, *args, **kwargs):
        if current.get("model") == "nonlinear":
            kwargs["pad_inches"] = max(float(kwargs.get("pad_inches", .025)), .10)
        if getattr(fig, "_labels_applied", False):
            return labels.save_original(fig, *args, **kwargs)
        fig._labels_applied = True
        return labels.save_labels(fig, *args, **kwargs)

    Figure.savefig = dispatch
    try:
        lp.configure_style()
        lp.M_VALUE = 30000
        for model, stem in [("linear", "Fig03b_Euler_core_tail_mechanism"),
                            ("plate", "Fig08b_plate_core_tail_mechanism")]:
            current["model"] = model
            labels.analytic = model == "linear"
            ref = pd.read_csv(data / f"reference_{model}.csv")
            frame = cov[cov.model == model]
            assert len(frame) == 150 and (frame.M == 30000).all()
            lp.MODEL_SETTINGS[model].update(output=output / model, stem=stem)
            lp.load_model_data = lambda unused, r=ref, d=frame: (
                r.y.to_numpy(), r.reference_density.to_numpy(), d.to_dict("records"), {"P": 30000})
            lp.plot_model(model)
            shutil.move(output / model / (stem + ".pdf"), output / (stem + ".pdf"))
        current["model"] = "nonlinear"
        labels.analytic = False
        nl.configure_style()
        nl.M_VALUE = 50000
        nl.OUTPUT_ROOT = output / "nonlinear"
        ref = pd.read_csv(data / "reference_nonlinear_density.csv")
        frame = cov[cov.model == "nonlinear"]
        assert len(frame) == 150 and (frame.M == 50000).all()
        assert (frame.P == 50000).all() and (frame.reference_N == 500000000).all()
        nl.load_formal_data = lambda: (ref.y.to_numpy(), ref.reference_density.to_numpy(),
            frame.to_dict("records"), dict(nonlinear_B=399, h=.08, n=768,
                                          P=50000, N_V=50000, N_ref=500000000))
        nl.plot()
        shutil.move(nl.OUTPUT_ROOT / "nonlinear_real_data_core_tail_mechanism.pdf",
                    output / "Fig10b_nonlinear_core_tail_mechanism.pdf")
        current.clear()
        ns.configure_plots()
        ns.FIG_DIR = ns.OUT_DIR = output
        frame = pd.read_csv(data / "internal_n.csv")
        assert len(frame) == 1458 and (frame.M == 30000).all()
        frame["method"] = frame.method.replace({"Student t": "t distribution"})
        for model, stem in [("linear", "Fig04_Euler_internal_sample_size"),
                            ("nonlinear", "Fig11_nonlinear_internal_sample_size")]:
            ns.plot_n_coverage(frame, model, stem)
    finally:
        Figure.savefig = native
        labels.save_original = native


def static_figure(kind, data, output, args):
    sys.path.insert(0, str(ROOT / "static"))
    if kind == "transport":
        importlib.import_module("draw_transport").build_figure(output)
    elif kind == "plate":
        importlib.import_module("draw_plate").draw_pdf(output / "Fig06_Kirchhoff_plate_model.pdf")
    elif kind in ("framework", "shift"):
        importlib.import_module(kind + "_layout").draw(output)
    elif kind == "euler":
        exe = shutil.which("pdflatex")
        if not exe:
            raise RuntimeError("pdflatex with standalone/TikZ must be on PATH")
        subprocess.run([exe, "-no-shell-escape", "-interaction=nonstopmode", "-halt-on-error",
            "-jobname=Fig02_Euler_beam_model", str(ROOT / "static/euler_model.tex")],
            cwd=output, check=True)
    elif kind == "outer":
        import numpy as np
        import pandas as pd
        if not args.font_roman or not args.font_italic:
            raise ValueError("Supply the original Times Roman/Italic font files with --font-roman and --font-italic")
        module = importlib.import_module("outer_workflow")
        module.configure_fonts(args.font_roman, args.font_italic)
        path = data / "outer_workflow_curves.csv"
        y = pd.read_csv(data / "reference_linear.csv").y.to_numpy()
        config = dict(n=768, h=.08, linear_B=399)
        module.draw_figure(y, np.loadtxt(path, delimiter=","), config, 32, 0.,
                           output, path, data / "reference_linear.csv", None)
        # The accepted, later label overlay is supplied separately. The earlier
        # automatic overlay is NOT claimed to be the current approved figure.
        shutil.copyfile(ROOT / "static/outer_overlay.tex", output / "outer_overlay.tex")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kind", required=True,
        choices=["python", "tables", "transport", "plate", "euler", "shift", "framework", "outer"])
    parser.add_argument("--data-dir", type=Path, default=DATA)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--font-roman", type=Path)
    parser.add_argument("--font-italic", type=Path)
    args = parser.parse_args()
    data = args.data_dir.resolve()
    output = (args.output_dir or ROOT.parent / "outputs" /
              ("historicalproduction_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f"))).resolve()
    if output == data or output in data.parents or data in output.parents:
        parser.error("Output must be separate from the input data tree")
    output.mkdir(parents=True, exist_ok=False)
    if args.kind == "python":
        python_figures(data, output)
    elif args.kind == "tables":
        from tables import export_tables
        export_tables(data, output)
    else:
        static_figure(args.kind, data, output, args)


if __name__ == "__main__":
    main()

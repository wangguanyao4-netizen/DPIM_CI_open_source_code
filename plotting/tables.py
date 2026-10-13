"""Archived first-trial table formatter; no experiment or manuscript writes."""
import csv
from pathlib import Path


def export_tables(data_dir, OUT):
    rows = list(csv.DictReader((data_dir / 'table_interval_values.csv').open(newline='', encoding='utf-8-sig')))
    INDICES = [1, 5, 9, 13, 17, 21, 25]
    METHODS = ['Student t', 'percentile bootstrap', 'bootstrap t']
    selected = []
    for model, label in [('linear', 'Euler beam'), ('plate', 'Kirchhoff plate'), ('nonlinear', 'nonlinear response')]:
        data = [r for r in rows if r['model'] == model]
        assert len(data) == 42 and {int(r['trial']) for r in data} == {1}
        r_values = sorted({int(r['R']) for r in data})
        lookup = {(int(r['location_index']), int(r['R']), r['method']): r for r in data}
        target = r'\(p_h(y)\)' if model == 'linear' else r'\(\widehat p_h^{\mathrm{ref}}(y)\)'
        yhead = r'\(y\)' if model != 'plate' else r'\(y\;(10^{-4}\,\mathrm m)\)'
        lines = [r'\vspace{0.6em}', r'\begin{minipage}{\linewidth}',
                 rf'\captionof{{table}}{{Pointwise 95\% {label} density confidence intervals from the first trial}}',
                 rf'\label{{tab:{model}-R-interval-endpoints}}', r'\centering',
                 r'\begingroup\small', r'\setlength{\tabcolsep}{3pt}',
                 r'\renewcommand{\arraystretch}{1.03}', r'\begin{tabular}{cccccc}', r'\toprule',
                 yhead + ' & ' + target + r' & \(R\) & Student t & percentile bootstrap & bootstrap t \\',
                 r'\midrule']
        for ix in INDICES:
            for ir, R in enumerate(r_values):
                first = lookup[ix, R, METHODS[0]]
                scale_y = 1e4 if model == 'plate' else 1
                scale_p = 1e-4 if model == 'plate' else 1
                y = float(first['y']) * scale_y
                ref = float(first['reference_density']) * scale_p
                cells = [f'{y:.5f}' if ir == 0 else '', f'{ref:.5f}' if ir == 0 else '', str(R)]
                for method in METHODS:
                    r = lookup[ix, R, method]
                    lower, upper = float(r['lower']) * scale_p, float(r['upper']) * scale_p
                    assert lower <= upper
                    cells.append(rf'\([{lower:.5f},\,{upper:.5f}]\)')
                    selected.append(r)
                lines.append(' & '.join(cells) + r' \\')
            if ix != INDICES[-1]:
                lines.append(r'\addlinespace[2pt]')
        lines.extend([r'\bottomrule', r'\end{tabular}', r'\endgroup'])
        if model == 'plate':
            lines.extend([r'\par\smallskip', r'{\small Densities and interval endpoints are expressed in \(10^4\,\mathrm m^{-1}\).}'])
        lines.append(r'\end{minipage}')
        (OUT / f'{model}_interval_endpoints.tex').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    assert len(selected) == 126
    with (OUT / 'table_interval_values.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(selected)

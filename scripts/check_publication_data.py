"""Validate archived publication data; do not generate any new trial."""
from pathlib import Path
import csv
import math
import statistics

DATA = Path(__file__).resolve().parents[1] / 'data/publication'


def read(name):
    with (DATA/name).open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def close(a, b, tol=2e-12):
    assert math.isclose(float(a),float(b),rel_tol=tol,abs_tol=tol), (a,b)


def main():
    cov = read('coverage.csv')
    ctrl = read('euler_pool_control.csv')
    assert len(cov) == 450 and len(ctrl) == 150
    for model, n, B, M, rvalues in [('linear',768,399,30000,{32,128}), ('plate',200,999,30000,{16,128}), ('nonlinear',768,399,50000,{16,128})]:
        rows = [r for r in cov if r['model']==model]
        assert len(rows)==150 and len({r['y'] for r in rows})==25
        assert {int(r['R']) for r in rows} == rvalues
        assert all(int(r['n'])==n and int(r['B'])==B and int(r['M'])==M for r in rows)
        if model == 'linear':
            assert all(r['experiment_type']=='fresh_independent' for r in rows)
            for r in rows:
                y=float(r['y'])
                close(r['truth'],math.exp(-y*y/(2*(1+.08**2)))/math.sqrt(2*math.pi*(1+.08**2)))
        else:
            assert all(int(float(r['P'])) == M for r in rows)
        if model == 'nonlinear':
            assert all(int(r['reference_N'])==500000000 for r in rows)
    assert all(r['dataset']=='euler_pool_control' for r in ctrl)
    for r in cov+ctrl:
        close(r['coverage'],int(r['hit_count'])/int(r['M']))
        assert 0 <= float(r['cp95_lower']) <= float(r['coverage']) <= float(r['cp95_upper']) <= 1

    for summary in read('coverage_summary.csv'):
        data = ctrl if summary['dataset']=='euler_pool_control' else cov
        rows = [r for r in data if r['model']==summary['model'] and r['method']==summary['method'] and r['R']==summary['R']]
        assert len(rows)==25
        close(summary['mean_coverage'], statistics.mean(float(r['coverage']) for r in rows))
        close(summary['mean_formula'], statistics.mean(float(r['formula_prediction']) for r in rows))
        close(summary['rmse'], math.sqrt(statistics.mean((float(r['coverage'])-float(r['formula_prediction']))**2 for r in rows)))

    full = read('saved_trial_intervals.csv')
    selected = read('table_interval_values.csv')
    assert len(full)==450 and len(selected)==126
    assert all(r in full for r in selected)
    for r in full:
        assert int(r['trial'])==1 and float(r['lower'])<=float(r['upper'])
        match = [x for x in cov if x['model']==r['model'] and x['R']==r['R'] and math.isclose(float(x['y']),float(r['y']),abs_tol=1e-12)]
        assert len(match)==3
        close(r['reference_density'],match[0]['truth'])
    internal = read('internal_n.csv')
    assert len(internal)==1458 and {int(r['n']) for r in internal} == {192,384,768}
    assert all(int(r['M'])==30000 and int(r['R'])==128 for r in internal)
    rec = read('recommended_R.csv')
    assert len(rec)==18
    for name in ('intervals_linear_nonlinear.csv','intervals_plate.csv'):
        rows=read(name)
        for r in rows:
            recommendation = [x for x in rec if all(x[k]==r[k] for k in ('model','scheme','method','n','B'))]
            assert len(recommendation)==1
            assert int(r['R'])==int(recommendation[0]['R_recommended'])
            assert float(r['ci_lower']) <= float(r['ci_upper'])
        groups={(r['model'],r['scheme'],r['method']) for r in rows}
        assert all(sum((r['model'],r['scheme'],r['method'])==g for r in rows)==201 for g in groups)
    for summary in read('gf_summary.csv'):
        design=read(summary['design_path'])
        dims=[k for k in design[0] if k.startswith('u')]
        close(sum(float(r['probability_weight']) for r in design),1)
        discrepancies=[]
        for dim in dims:
            before=0.; maximum=0.
            for row in sorted(design,key=lambda r:float(r[dim])):
                after=before+float(row['probability_weight'])
                maximum=max(maximum,abs(before-float(row[dim])),abs(after-float(row[dim])))
                before=after
            discrepancies.append(maximum)
        close(summary['GF'],max(discrepancies))
    print('PASS: primary/control counts, analytical Euler target, coverage and RMSE arithmetic, first-trial references, recommended-R curves, internal-n budgets and three original GF values.')
    print('No new bootstrap trials, density estimates or finite element solutions were computed.')


if __name__ == '__main__':
    main()

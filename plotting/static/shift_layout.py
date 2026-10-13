"""Final historical arrow-tip correction, without paper workflow or QA imports."""
from pathlib import Path
import importlib.util
ROOT = Path(__file__).resolve().parent

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def draw(output_dir):
    d=load('mechanism',ROOT/'shift_geometry.py')
    fig,geometry=d.build()
    for t in fig.findobj(d.Text):
        if t.get_text()==r'Deterministic response surface, $Y=g(\boldsymbol{\theta})$':
            t.set_text(r'Deterministic response surface, $Y=g(T(u))$')
        elif t.get_text() in [r'$\theta_1$',r'$\theta_2$']:
            t.set_text(r'$u_1$' if t.get_text()==r'$\theta_1$' else r'$u_2$')
            t.set_rotation(0)
    surf=next(a for a in fig.axes if a.name=='3d')
    surf.xaxis.set_rotate_label(False);surf.yaxis.set_rotate_label(False)
    surf.xaxis.label.set_rotation(0);surf.yaxis.label.set_rotation(0)
    for ax in fig.axes:
        if ax.get_ylabel()==r'$u_2$' and ax is not surf:
            ax.yaxis.label.set_rotation(0);ax.yaxis.labelpad=4
    overlay=fig.axes[-1]
    arrows=[p for p in overlay.patches if isinstance(p,d.FancyArrowPatch)]
    audit=[]
    for index,old in enumerate(arrows):
        vertices=old._path_original.vertices.copy()
        a,b,c,e=vertices[-4:]
        # Subdivide the last cubic exactly; every retained point follows the
        # original curve. Only the final 10 pt approaching the marker is removed.
        tt=d.np.linspace(0,1,20001)[:,None]
        samples=(1-tt)**3*a+3*(1-tt)**2*tt*b+3*(1-tt)*tt**2*c+tt**3*e
        points_per_unit=fig.get_size_inches()*72
        distances=d.np.linalg.norm((samples-e)*points_per_unit,axis=1)
        k=d.np.flatnonzero(distances>=10)[-1]
        t=float(tt[k,0])
        ab=(1-t)*a+t*b;bc=(1-t)*b+t*c;ce=(1-t)*c+t*e
        abc=(1-t)*ab+t*bc;bce=(1-t)*bc+t*ce
        endpoint=(1-t)*abc+t*bce
        trimmed=d.np.vstack([vertices[:-3],ab,abc,endpoint])
        assert d.np.array_equal(trimmed[:-3],vertices[:-3])
        path=d.MPath(trimmed,[d.MPath.MOVETO]+[d.MPath.CURVE4]*(len(trimmed)-1))
        old.remove()
        overlay.add_patch(d.FancyArrowPatch(path=path,transform=overlay.transData,
            arrowstyle='-|>',mutation_scale=d.ARROW_HEAD,linewidth=d.ARROW_WIDTH,
            color=d.COLORS[index//2],alpha=d.ARROW_ALPHA,zorder=12,
            joinstyle='round',capstyle='round'))
        audit.append(dict(arrow=index,original_curve_preserved=True,
            target_clearance_pt=float(d.np.linalg.norm((endpoint-e)*points_per_unit))))
    fig.canvas.draw()
    d.plt.rcParams.update({'pgf.texsystem':'pdflatex','pgf.rcfonts':False,
        'pgf.preamble':d.plt.rcParams['text.latex.preamble']})
    target=Path(output_dir)/'DPIM_random_shift_response_density_mechanism.pdf'
    fig.savefig(target,backend='pgf');d.plt.close(fig)
    return target

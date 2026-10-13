"""Local spacing repair; all original module centers and routes are retained."""
from pathlib import Path
import importlib.util
import os

ROOT=Path(__file__).resolve().parent


def draw(output_dir):
    spec=importlib.util.spec_from_file_location('original',ROOT/'framework_geometry.py')
    d=importlib.util.module_from_spec(spec);spec.loader.exec_module(d)
    # The response is written without the extra intermediate symbol Y_i^(r).
    # The density expression follows eq:dpim-estimator verbatim.
    replacements={
        r'$Y_i^{(r)}=g(T(u_i^{(r)}))$':r'$g(T(u_i^{(r)}))$',
        r'$Q_r(y)=\sum_i P_i^0 K_h(y-Y_i^{(r)})$':r'$Q_r(y)=\sum_{i=1}^{n}P_i^0H_y(u_i^{(r)})$'}
    def module(ax,x,y,lines,color,w=2.63,h=.7,size=11.5):
        density=lines[0]=='Density calculation'
        h=1.10 if density else max(h,.98)
        patch=d.FancyBboxPatch((x-w/2,y-h/2),w,h,
            boxstyle='round,pad=0.01,rounding_size=0.06',
            facecolor=color,edgecolor='black',linewidth=1.05,zorder=3)
        ax.add_patch(patch)
        labels=[]
        for i,s in enumerate(lines):
            offset=(.31 if i==0 else -.12) if density else (.255 if i==0 else -.235)
            labels.append(d.txt(ax,x,y+offset,replacements.get(s,s),size,zorder=4))
        d.MODULES.append((patch,labels))
    d.module=module
    fig=d.architecture()
    from matplotlib.patches import FancyArrowPatch
    for arrow in fig.axes[0].patches:
        if not isinstance(arrow,FancyArrowPatch) or arrow._posA_posB is None:continue
        a,b=arrow._posA_posB
        if abs(a[0]-b[0])>1e-8 or b[1]<=a[1]:continue
        mid=(a[1]+b[1])/2
        aligned=[p for p,_ in d.MODULES if abs(p.get_x()+p.get_width()/2-a[0])<1e-8]
        below=[p for p in aligned if p.get_y()+p.get_height()/2<mid]
        above=[p for p in aligned if p.get_y()+p.get_height()/2>mid]
        y0=max(p.get_y()+p.get_height() for p in below)+.035 if below else a[1]
        y1=min(p.get_y() for p in above)-.035 if above else b[1]
        assert y1>y0
        arrow.set_positions((a[0],y0),(b[0],y1))
        arrow.shrinkA=0;arrow.shrinkB=0
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    for patch,labels in d.MODULES:
        a,b=[t.get_window_extent(renderer) for t in labels]
        assert a.y0-b.y1>3,(labels[0].get_text(),a.y0-b.y1)
    # Give the original repetition enclosures a small top margin.
    for p in fig.axes[0].patches:
        if isinstance(p,d.FancyBboxPatch) and abs(p.get_width()-3.4)<1e-6:
            p.set_height(4.2)
    d.plt.rcParams.update({'pgf.texsystem':'pdflatex','pgf.rcfonts':False,
        'pgf.preamble':d.plt.rcParams['text.latex.preamble']})
    target=Path(output_dir)/'DPIM_unified_statistical_framework.pdf'
    fig.savefig(target,backend='pgf',bbox_inches='tight',pad_inches=.045)
    d.plt.close(fig)
    return target

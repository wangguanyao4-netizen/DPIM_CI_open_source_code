"""Vector method diagrams; illustrative geometry, not manuscript experiment data."""
from __future__ import annotations

import hashlib
import json
from itertools import product
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle
from matplotlib.path import Path as MPath
import numpy as np

ROOT = Path(__file__).resolve().parent
BLACK = '#000000'
BLUE = '#3b6483'
RUST = '#a36b4b'
GREEN = '#65847a'
EDGE = '#82909a'
FILLS = ['#e3ebf0', '#eee2ce', '#dce9e2', '#f2f3f3', '#f2f3f3', '#f2f3f3', '#f2f3f3', '#f2f3f3', '#f2f3f3']
POINTS = np.array([[.77,.63],[.18,.26],[.47,.37],[.14,.76],[.43,.86],
                   [.82,.14],[.48,.08],[.08,.03],[.76,.93]])
SHIFTS = [np.array([.34,.0]), np.array([.13,.29])]
plt.rcParams.update({'text.usetex':True, 'font.family':'serif',
    'text.latex.preamble':r'\usepackage{lmodern}\usepackage{amsmath,amssymb}\usepackage{bm}',
    'font.size':12, 'text.color':BLACK, 'axes.labelcolor':BLACK,
    'xtick.color':BLACK, 'ytick.color':BLACK, 'axes.edgecolor':BLACK,
    'svg.fonttype':'path', 'pdf.compression':6, 'savefig.facecolor':'white'})


def clip(poly, normal, offset):
    if len(poly) == 0:
        return np.empty((0,2))
    output = []
    a = poly[-1]
    va = a@normal-offset
    for b in poly:
        vb = b@normal-offset
        if (va <= 1e-12) != (vb <= 1e-12):
            output.append(a+(b-a)*va/(va-vb))
        if vb <= 1e-12:
            output.append(b)
        a, va = b, vb
    return np.asarray(output).reshape((-1,2))


def unit_clip(poly):
    for n, b in [(np.array([-1.,0]),0), (np.array([1.,0]),1),
                 (np.array([0.,-1]),0), (np.array([0.,1]),1)]:
        poly = clip(poly,n,b)
    return poly


def area(poly):
    if len(poly) < 3:
        return 0.
    return abs(np.dot(poly[:,0],np.roll(poly[:,1],1))-np.dot(poly[:,1],np.roll(poly[:,0],1)))/2


def periodic_cells():
    cells = []
    for v in POINTS:
        poly = np.array([[-1.,-1],[2.,-1],[2.,2],[-1.,2]])
        for delta in product([-1,0,1],repeat=2):
            for other in POINTS+delta:
                if np.linalg.norm(other-v) < 1e-10:
                    continue
                poly = clip(poly, other-v, (other@other-v@v)/2)
        cells.append(poly)
    return cells


def pieces(poly, shift):
    out = []
    for delta in product([-2,-1,0,1],repeat=2):
        part = unit_clip(poly+shift+delta)
        if area(part) > 1e-10:
            out.append(part)
    return out


def text(ax,x,y,s,size=12,ha='center',va='center',**kw):
    return ax.text(x,y,s,fontsize=size,ha=ha,va=va,color=BLACK,**kw)


def arrow(ax,a,b,color=BLACK,lw=1.1,**kw):
    p=FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=11,
                      linewidth=lw,color=color,shrinkA=2,shrinkB=2,**kw)
    ax.add_patch(p)
    return p


def route(ax, vertices, color=BLACK, lw=1.1):
    codes=[MPath.MOVETO]+[MPath.LINETO]*(len(vertices)-1)
    ax.add_patch(FancyArrowPatch(path=MPath(vertices,codes),arrowstyle='-|>',
                               mutation_scale=11,linewidth=lw,color=color,fill=False))


def canvas(height):
    fig=plt.figure(figsize=(11.8,height))
    ax=fig.add_axes([0,0,1,1])
    ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    return fig,ax


def torus(fig,box,cells,shift,highlight_wrap=False):
    ax=fig.add_axes(box)
    ax.set(xlim=(-.03,1.03),ylim=(-.03,1.03),aspect='equal');ax.axis('off')
    for i,cell in enumerate(cells):
        for poly in pieces(cell,shift):
            ax.add_patch(Polygon(poly,closed=True,facecolor=FILLS[i],edgecolor=EDGE,lw=.65))
    ax.add_patch(Rectangle((0,0),1,1,fill=False,edgecolor=BLACK,lw=1.0))
    u=(POINTS+shift)%1
    ax.scatter(u[:,0],u[:,1],s=19,c=BLACK,zorder=4)
    for i in [0,1,2]:
        x,y=u[i]
        label=rf'$P_{i+1}^0$'
        text(ax,x+.034 if x<.9 else x-.04,y-.067 if y>.84 else y+.056,label,size=11,
             ha='left' if x<.9 else 'right',bbox=dict(facecolor=FILLS[i],edgecolor='none',pad=.5))
    text(ax,-.015,-.035,'0',size=10,ha='right',va='top')
    text(ax,1,-.035,'1',size=10,va='top')
    text(ax,-.02,1,'1',size=10,ha='right')
    if highlight_wrap:
        y=POINTS[0,1]
        ax.scatter([POINTS[0,0]],[y],s=27,facecolors='white',edgecolors=BLUE,lw=1,zorder=5)
        arrow(ax,(POINTS[0,0]+.025,y),(1.005,y),color=BLUE,lw=1.3)
        arrow(ax,(-.012,y),(u[0,0]-.014,y),color=BLUE,lw=1.3)
        text(ax,.86,y-.10,r'$\Delta_1$',size=11)
    return ax


def build_mechanism():
    cells=periodic_cells()
    weights=np.array([area(p) for p in cells])
    assert abs(weights.sum()-1)<1e-12
    deviations=[]
    for shift in [np.zeros(2)]+SHIFTS:
        areas=np.array([sum(area(p) for p in pieces(c,shift)) for c in cells])
        deviations.append(float(np.max(abs(areas-weights))))
        assert np.allclose(areas,weights,rtol=0,atol=1e-12)
        grid=(np.arange(83)+.413)/83
        sample=np.array(list(product(grid,grid)))
        cover=np.zeros(len(sample),dtype=int)
        delta=abs(sample[:,None,:]-(POINTS+shift)[None,:,:]%1)
        expected=np.argmin(np.sum(np.minimum(delta,1-delta)**2,axis=2),axis=1)
        for i,cell in enumerate(cells):
            for part in pieces(cell,shift):
                inside=MPath(part).contains_points(sample)
                cover+=inside
                assert np.all(expected[inside]==i)
        assert np.all(cover==1), 'Periodic partition has overlaps or holes.'
    fig,ax=canvas(6.8)
    xs=[.043,.376,.709]
    titles=[r'\textbf{(a)} Initial representative points',r'\textbf{(b)} First outer repetition',r'\textbf{(c)} Second outer repetition']
    for x,title in zip(xs,titles):
        text(ax,x,.935,title,size=13,ha='left')
    text(ax,.5,.873,r'One fixed initial set; independent shifts $\Delta_1$ and $\Delta_2$',size=12)
    # Equal physical axis lengths retain the unit-square geometry.
    boxes=[[x,.393,.247,.4285] for x in xs]
    torus(fig,boxes[0],cells,np.zeros(2))
    torus(fig,boxes[1],cells,SHIFTS[0],True)
    torus(fig,boxes[2],cells,SHIFTS[1])
    formulas=[r'$v_i,\ V_i^0,\ P_i^0$',r'$u_i^{(1)}=(v_i+\Delta_1)\bmod 1$',r'$u_i^{(2)}=(v_i+\Delta_2)\bmod 1$']
    for x,s in zip(xs,formulas):text(ax,x+.1235,.354,s,size=12)
    text(ax,.1665,.306,r'$P_i^0\geq0,\quad\sum_iP_i^0=1$',size=12)
    text(ax,.4995,.306,r'Periodic boundary crossing',size=11)
    text(ax,.8325,.306,r'Translated cells; unchanged probabilities',size=11)
    ax.plot([.04,.96],[.264,.264],color='#c5cbd0',lw=.7)
    text(ax,.043,.229,r'\textbf{(d)} From each shifted point set to a density estimate',ha='left',size=13)
    text(ax,.157,.149,r'$u_i^{(r)}\ \xrightarrow{\ T\ }\ \bm\theta_i^{(r)}\ \xrightarrow{\ g\ }\ Y_i^{(r)}$',size=13)
    text(ax,.157,.087,'Deterministic structural response',size=11)
    arrow(ax,(.285,.147),(.35,.147),color=EDGE)
    text(ax,.499,.146,r'$Q_r(y)=\displaystyle\sum_{i=1}^{n}P_i^0K_h\!\left(y-Y_i^{(r)}\right)$',size=13)
    text(ax,.499,.087,'Fixed weights; varying responses',size=11)
    arrow(ax,(.648,.147),(.708,.147),color=EDGE)
    ca=fig.add_axes([.732,.072,.217,.135])
    y=np.linspace(-.9,1.8,250)
    for j,(shift,color) in enumerate(zip(SHIFTS,[BLUE,RUST])):
        u=(POINTS+shift)%1
        response=u[:,0]+.55*np.sin(2*np.pi*u[:,1])+.15*u[:,1]
        h=.17
        q=np.exp(-.5*((y[:,None]-response)/h)**2)@weights/(np.sqrt(2*np.pi)*h)
        ca.plot(y,q,color=color,lw=1.35,label=rf'$Q_{j+1}(y)$')
    ca.spines[['top','right','left']].set_visible(False)
    ca.spines['bottom'].set_linewidth(.7)
    ca.set(xticks=[],yticks=[],xlim=(y[0],y[-1]),ylim=(0,None))
    ca.legend(frameon=False,fontsize=10,ncol=2,loc='upper center',bbox_to_anchor=(.5,1.36),handlelength=1.5,columnspacing=.9)
    ca.text(1.02,0,r'$y$',transform=ca.transAxes,fontsize=11)
    text(ax,.5,.031,r'$\mathbb{E}\!\left[Q_r(y)\mid\mathcal D_n\right]=p_h(y)$',size=12)
    return fig,dict(cell_volumes=weights.tolist(),maximum_volume_error=max(deviations),
                    shifts=[s.tolist() for s in SHIFTS],illustrative_model='g(u)=u1+0.55*sin(2*pi*u2)+0.15*u2',
                    illustrative_kernel_h=.17,initial_points=POINTS.tolist())


def box(ax,x,y,w,h,title,lines,accent=BLUE):
    ax.add_patch(Rectangle((x,y),w,h,facecolor='white',edgecolor='#9ba6ae',lw=.8))
    ax.plot([x,x+w],[y+h,y+h],color=accent,lw=2)
    text(ax,x+w/2,y+h-.035,title,size=13)
    step=.030
    first=y+h-.077
    for i,line in enumerate(lines):text(ax,x+w/2,first-i*step,line,size=12)


def build_framework():
    fig,ax=canvas(6.65)
    text(ax,.042,.958,r'\textbf{Density calculation and statistical accuracy in one framework}',ha='left',size=15)
    box(ax,.042,.743,.269,.142,'Deterministic structural model',
        [r'$\bm\Theta=T(u),\quad Y=g(\bm\Theta)$','Prescribed model and smoothing $h$'],RUST)
    box(ax,.365,.743,.269,.142,'Fixed representative points',
        [r'$\mathcal D_n=\{v_i,V_i^0,P_i^0\}_{i=1}^{n}$','Probabilities assigned once'],GREEN)
    box(ax,.688,.743,.269,.142,'Independent random shifts',
        [r'$u_i^{(r)}=(v_i+\Delta_r)\bmod1$','Same weights in every repetition'],BLUE)
    arrow(ax,(.637,.813),(.685,.813),color=EDGE)
    route(ax,[(.1765,.74),(.1765,.688),(.351,.688)],RUST)
    route(ax,[(.8225,.74),(.8225,.688),(.648,.688)],BLUE)
    text(ax,.5,.687,r'$Q_r(y)=\displaystyle\sum_{i=1}^{n}P_i^0K_h\!\left\{y-g(T(u_i^{(r)}))\right\}$',size=14)
    text(ax,.5,.626,'Repeated DPIM density estimates',size=12)
    route(ax,[(.5,.608),(.5,.575),(.1765,.575),(.1765,.553)],BLUE)
    box(ax,.042,.397,.269,.153,'Preliminary outer sample',
        [r'$Q_1^{\mathrm{pilot}},\ldots,Q_{L_{\mathrm{pilot}}}^{\mathrm{pilot}}$',r'Independent preliminary shifts'],BLUE)
    box(ax,.365,.397,.269,.153,'Moments of the density estimator',
        ['Variance, skewness, excess kurtosis',r'$\widehat v_2(y),\quad\widehat\gamma_n(y),\quad\widehat\kappa_n(y)$'],BLUE)
    box(ax,.688,.397,.269,.153,'Unified coverage approximation',
        [r'Student t, percentile, bootstrap t',r'$C_m^{(2)}(R,B,n;y)$',r'Edgeworth expansion; finite $B$'],BLUE)
    arrow(ax,(.313,.474),(.362,.474),BLUE)
    arrow(ax,(.636,.474),(.685,.474),BLUE)
    text(ax,.5,.35,'Deterministic responses enter the coverage error through these moments',size=12)
    arrow(ax,(.8225,.394),(.8225,.271),BLUE)
    box(ax,.688,.104,.269,.164,'Recommended outer sample size',
        [r'Prescribed $\varepsilon_{\mathrm c}$ and response region',r'$R_{\min,m}(n,B;\mathcal Y_c)$','Second order asymptotic criterion'],BLUE)
    box(ax,.365,.104,.269,.164,'New independent outer sample',
        [r'$Q_1(y),\ldots,Q_R(y)$','At the recommended $R$','Independent of the preliminary sample'],GREEN)
    box(ax,.042,.104,.269,.164,'Density and pointwise intervals',
        [r'$\bar Q_R(y),\quad I_t(y),\ I_p(y),\ I_{bt}(y)$',r'Coverage target: $p_h(y)$'],RUST)
    arrow(ax,(.685,.184),(.637,.184),BLUE)
    arrow(ax,(.362,.184),(.314,.184),GREEN)
    text(ax,.5,.043,r'Prescribed $n$, $h$, $B$ and confidence level $1-\alpha$; fixed model throughout',size=12)
    return fig







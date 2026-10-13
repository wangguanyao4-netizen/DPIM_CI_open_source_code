"""Minimal periodic geometry and a Transformer-style DPIM architecture diagram."""
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
v.ROOT = ROOT
np, plt = v.np, v.plt
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle
from matplotlib.path import Path as MPath


def txt(ax, x, y, s, size=12, **kw):
    return v.text(ax, x, y, s, size=size, **kw)


def geometry():
    fig = plt.figure(figsize=(11.4, 7.0))
    outer = fig.add_axes([0, 0, 1, 1])
    outer.set(xlim=(0, 1), ylim=(0, 1)); outer.axis('off')
    cells = v.periodic_cells()
    weights = np.array([v.area(c) for c in cells])
    shifts = [np.zeros(2), np.array([.64, .15]), np.array([.14, .72])]
    geometry_errors = []
    grid = (np.arange(83) + .413) / 83
    sample = np.array([(x, y) for x in grid for y in grid])
    xs = [.045, .375, .705]
    # Both realizations start from the same initial design, not from one another.
    for x in xs[1:]:
        v.route(outer,[(xs[0]+.125,.907),(xs[0]+.125,.946),
                      (x+.125,.946),(x+.125,.907)],lw=.95)
    txt(outer,.50,.978,r'$\Delta_1$',13)
    txt(outer,.83,.978,r'$\Delta_2$',13)
    txt(outer,.667,.978,r'$\Delta_1\perp\Delta_2$',11)
    responses = []
    for j, delta in enumerate(shifts):
        ax = fig.add_axes([xs[j], .49, .25, .40715])
        ax.set(xlim=(-.015, 1.015), ylim=(-.015, 1.015), aspect='equal')
        ax.axis('off')
        counts = np.zeros(len(sample), dtype=int)
        calculated = []
        for i, cell in enumerate(cells):
            fragments = v.pieces(cell, delta)
            calculated.append(sum(v.area(p) for p in fragments))
            for p in fragments:
                counts += MPath(p).contains_points(sample)
                ax.add_patch(Polygon(p, closed=True, facecolor='#dce9ef' if i == 2 else 'white',
                                     edgecolor='#7b7b7b', linewidth=.7))
        assert np.all(counts == 1)
        err = float(np.max(abs(weights - calculated)))
        geometry_errors.append(err)
        assert err < 1e-12 and abs(sum(calculated)-1) < 1e-12
        ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, edgecolor='black', linewidth=1.05))
        points = (v.POINTS + delta) % 1
        ax.scatter(points[:, 0], points[:, 1], s=15, color='black', zorder=4)
        x, y = points[2]
        txt(ax, x+.04, y+.075, r'$v_i$' if j == 0 else rf'$u_i^{{({j})}}$', 12, ha='left')
        if j == 0:
            txt(ax,x+.055,y-.075,r'$V_i^0$',12)
        elif j == 1:
            txt(ax,.165,.422,r'$P_i^0$',11)
        else:
            txt(ax,.53,.19,r'$P_i^0$',11)
        if j:
            start = v.POINTS[2]
            axis = int(np.argmin((1-start)/delta))
            crossing = start + delta*((1-start[axis])/delta[axis])
            wrapped = crossing.copy(); wrapped[axis] = 0
            ax.scatter(*start,s=19,facecolors='white',edgecolors='#3b6483',linewidths=.8,zorder=5)
            v.arrow(ax,start,crossing,color='#3b6483',lw=.9,linestyle='--',zorder=2)
            v.arrow(ax,wrapped,points[2],color='#3b6483',lw=.9,linestyle='--',zorder=2)
        txt(ax, -.01, -.035, '$0$', 10, ha='right', va='top')
        txt(ax, 1, -.035, '$1$', 10, va='top')
        txt(ax, -.025, 1, '$1$', 10, ha='right')
        txt(outer,xs[j]-.027,.916,rf'\textbf{{({chr(97+j)})}}',13)
        if j:
            responses.append(points[:,0]+.55*np.sin(2*np.pi*points[:,1])+.15*points[:,1])
    txt(outer,.17,.422,r'$P_i^0=|V_i^0|$',13)
    for j in [1,2]:
        txt(outer,xs[j]+.125,.432,rf'$u_i^{{({j})}}=(v_i+\Delta_{j})\bmod1$',12)
        v.arrow(outer,(xs[j]+.125,.401),(xs[j]+.125,.337),lw=1)
        txt(outer,xs[j]+.16,.37,r'$g\circ T,\ K_h$',11,ha='left')
    w_ax = fig.add_axes([.065,.145,.21,.19])
    colors = ['#b5b5b5']*len(weights); colors[2] = '#709caf'
    w_ax.bar(np.arange(1,10),weights,color=colors,width=.5,edgecolor='black',linewidth=.5)
    w_ax.set(xlim=(.5,9.5),ylim=(0,max(weights)*1.16),xticks=[1,3,6,9],yticks=[])
    w_ax.set_xlabel(r'$i$',fontsize=12,labelpad=2)
    w_ax.set_ylabel(r'$P_i^0$',fontsize=12,rotation=0,labelpad=12)
    w_ax.spines[['right','top']].set_visible(False)
    w_ax.tick_params(axis='x',labelsize=10,length=2)
    txt(outer,.17,.368,r'$\sum_i P_i^0=1$',12)
    density_checks = []
    for j, response in enumerate(responses,1):
        d_ax = fig.add_axes([xs[j]+.01,.145,.23,.19])
        h = .17
        yy = np.linspace(-.8,2.0,700)
        components = weights[:,None]*np.exp(-.5*((yy-response[:,None])/h)**2)/(h*np.sqrt(2*np.pi))
        total = components.sum(axis=0)
        for component in components:
            d_ax.plot(yy,component,color='#b8b8b8',lw=.6)
        d_ax.fill_between(yy,0,components[2],color='#dce9ef')
        d_ax.plot(yy,components[2],color='#3b6483',lw=.8)
        d_ax.plot(yy,total,color='black',lw=1.2)
        d_ax.set(xlim=(-.8,2.0),ylim=(0,total.max()*1.23),xticks=[],yticks=[])
        d_ax.spines[['right','top']].set_visible(False)
        d_ax.scatter([response[2]],[0],s=15,color='black',clip_on=False,zorder=5)
        txt(d_ax,response[2],-.08*total.max(),rf'$Y_i^{{({j})}}$',11,va='top')
        txt(d_ax,1.99,-.08*total.max(),'$y$',12,ha='right',va='top')
        txt(d_ax,1.97,total.max()*1.12,rf'$Q_{j}(y)$',13,ha='right')
        wide = np.linspace(-4,5,5000)
        values = weights[:,None]*np.exp(-.5*((wide-response[:,None])/h)**2)/(h*np.sqrt(2*np.pi))
        integral = float(np.trapezoid(values.sum(axis=0),wide))
        assert abs(integral-1)<1e-12
        density_checks.append(integral)
    txt(outer,.5,.045,r'$Q_r(y)=\sum_{i=1}^{n}P_i^0K_h(y-Y_i^{(r)}),\qquad Y_i^{(r)}=g(T(u_i^{(r)}))$',14)
    return fig, {'shifts': [s.tolist() for s in shifts], 'maximum_volume_error': max(geometry_errors),
                 'illustrative_density_integrals':density_checks,
                 'illustrative_response':'g(T(u)) = u1 + 0.55*sin(2*pi*u2) + 0.15*u2; h = 0.17',
                 'grid_points_per_panel': len(sample), 'grid_coverage_exactly_once': True}


PINK = '#f8dddd'
BLUE = '#cde8ef'
YELLOW = '#f2efbd'
PEACH = '#f8e2c1'
GREEN = '#d7ead5'
LILAC = '#deddea'


MODULES = []


def module(ax, x, y, lines, color, w=2.63, h=.7, size=11.5):
    patch = FancyBboxPatch((x-w/2, y-h/2), w, h,
                 boxstyle='round,pad=0.01,rounding_size=0.06',
                 facecolor=color, edgecolor='black', linewidth=1.05, zorder=3)
    ax.add_patch(patch)
    labels = []
    for i, s in enumerate(lines):
        labels.append(txt(ax, x, y+(len(lines)-1-2*i)*.145, s, size, zorder=4))
    MODULES.append((patch, labels))


def up(ax, x, y0, y1):
    v.arrow(ax, (x,y0), (x,y1), lw=1.05, zorder=2)


def architecture():
    fig = plt.figure(figsize=(7.8, 8.2))
    ax = fig.add_axes([.025, .02, .95, .96])
    ax.set(xlim=(0,10), ylim=(0,11.6)); ax.axis('off')
    left, right = 2.6, 7.45
    for x, repeat in [(left,r'$L_{\mathrm{pilot}}\times$'), (right,r'$R\times$')]:
        ax.add_patch(FancyBboxPatch((x-1.7,2.15),3.4,4.0,
                     boxstyle='round,pad=0.015,rounding_size=0.17',
                     facecolor='#f3f3f3',edgecolor='black',linewidth=1.1,zorder=0))
        # The enclosure denotes independent repetitions, not sequential layers.
        txt(ax, x-1.85 if x == left else x+1.85,4.15,repeat,11,
            ha='right' if x == left else 'left')
        txt(ax,x,.33,r'$v_i,\ P_i^0,\quad i=1,\ldots,n$',12)
        module(ax,x,1.13,['Representative points','Fixed probabilities'],PINK,h=.76)
        up(ax,x,.55,.73)
        up(ax,x,1.53,2.54)
        module(ax,x,2.96,['Random shift',r'$u_i^{(r)}=(v_i+\Delta_r)\bmod1$'],PEACH,h=.82,size=12)
        up(ax,x,3.38,3.99)
        module(ax,x,4.41,['Deterministic response',r'$Y_i^{(r)}=g(T(u_i^{(r)}))$'],BLUE,h=.82,size=12)
        up(ax,x,4.83,5.24)
        module(ax,x,5.64,['Density calculation',r'$Q_r(y)=\sum_i P_i^0 K_h(y-Y_i^{(r)})$'],BLUE,h=.78,size=10.8)
        # A genuine bypass supplies unchanged weights to the kernel sum.
        v.route(ax,[(x+1.33,1.13),(x+1.51,1.13),(x+1.51,5.64),(x+1.33,5.64)],lw=.95)
        txt(ax,x+1.73,1.62,r'$P_i^0$',10,ha='left')
    up(ax,left,6.06,6.69)
    module(ax,left,7.06,['Moment estimates',r'$\widehat v_2,\ \widehat\gamma_n,\ \widehat\kappa_n$'],GREEN,h=.72)
    up(ax,left,7.43,7.98)
    module(ax,left,8.36,['Coverage approximation',r'$C_m^{(2)}(R,B,n;y)$'],PEACH,h=.72)
    up(ax,left,8.73,9.29)
    module(ax,left,9.66,['Outer sample size selection',r'$R_{\min,m}$'],YELLOW,h=.72,size=11.5)
    txt(ax,left,10.65,'Preliminary calculation',13)
    # R sets the right-hand repetition count; it is not density data.
    v.route(ax,[(left+1.33,9.66),(5.15,9.66),(5.15,4.15),(right-1.73,4.15)],lw=1.05)
    txt(ax,5.36,4.48,r'$R$',12)
    up(ax,right,6.06,6.69)
    module(ax,right,7.06,['Independent outer sample',r'$Q_1(y),\ldots,Q_R(y)$'],GREEN,h=.72)
    up(ax,right,7.43,7.98)
    module(ax,right,8.36,['Mean and standard error',r'$\bar Q_R(y),\quad S_{Q,R}(y)/\sqrt{R}$'],LILAC,h=.72,size=12)
    up(ax,right,8.73,9.29)
    module(ax,right,9.66,['Pointwise intervals',r'$I_t(y),\ I_p(y),\ I_{bt}(y)$'],GREEN,h=.72)
    # Bootstrap quantiles require the outer sample itself, not just mean/SE.
    v.route(ax,[(right+1.33,7.06),(right+1.63,7.06),(right+1.63,9.66),(right+1.33,9.66)],lw=.95)
    up(ax,right,10.04,10.52)
    txt(ax,right,10.83,r'Inference for $p_h(y)$',13)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for patch, labels in MODULES:
        bounds = patch.get_window_extent(renderer)
        for label in labels:
            b = label.get_window_extent(renderer)
            assert b.x0 > bounds.x0+3 and b.x1 < bounds.x1-3, label.get_text()
            assert b.y0 > bounds.y0+2 and b.y1 < bounds.y1-2, label.get_text()
    return fig





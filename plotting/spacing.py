def spaced_save(native):
    def save(fig,*args,**kwargs):
        # The original label adapter is invoked once per figure, including
        # the first export; later formats use the same resolved geometry.
        if not getattr(fig,'_spacing_fixed',False) and len(fig.axes)==3:
            fix_annotation_separators(fig)
            fig._spacing_fixed=True
        return native(fig,*args,**kwargs)
    return save


def fix_annotation_separators(fig):
    import numpy as np
    from matplotlib.lines import Line2D
    ax=fig.axes[1];fig.canvas.draw();renderer=fig.canvas.get_renderer()
    textboxes=[t.get_window_extent(renderer).padded(4) for t in ax.texts if t.get_text()]
    lines=[l for l in ax.lines if str(l.get_color()).lower()=='#e58b31']
    xs=sorted(set(float(l.get_xdata()[0]) for l in lines))
    assert len(xs)==4
    for line in lines:line.remove()
    transform=ax.get_xaxis_transform();inverse=ax.transAxes.inverted()
    for x in xs:
        screen_x=transform.transform((x,0))[0];segments=[(0.,1.)]
        for box in textboxes:
            if box.x0<=screen_x<=box.x1:
                lo=inverse.transform((screen_x,box.y0))[1]
                hi=inverse.transform((screen_x,box.y1))[1]
                next_segments=[]
                for a,b in segments:
                    if hi<=a or lo>=b:next_segments.append((a,b))
                    else:
                        if a<lo:next_segments.append((a,min(b,lo)))
                        if hi<b:next_segments.append((max(a,hi),b))
                segments=next_segments
        for lo,hi in segments:
            if hi-lo<.025:continue
            line=Line2D([x,x],[lo,hi],transform=transform,color='#e58b31',linewidth=.75,linestyle=(0,(4,4)))
            ax.add_line(line)
            p=transform.transform([[x,lo],[x,hi]])
            for box in textboxes:
                assert not (box.x0<screen_x<box.x1 and max(p[:,1].min(),box.y0)<min(p[:,1].max(),box.y1)-1e-6)
    fig._separator_audit=dict(boundaries=len(xs),textboxes=len(textboxes),overlaps=0)


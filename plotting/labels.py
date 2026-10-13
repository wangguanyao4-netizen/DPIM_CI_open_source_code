"""Accepted terminology adapter; no filesystem or manuscript side effects."""
from matplotlib.figure import Figure
from matplotlib.text import Text

save_original = Figure.savefig
analytic = False


def save_labels(fig, *args, **kwargs):
    for obj in fig.findobj(match=Text):
        value = obj.get_text().replace('bootstrap-t', 'bootstrap t')
        value = value.replace('second-order prediction', 'second order prediction')
        value = value.replace('pilot moment uncertainty', 'moment estimation error')
        value = value.replace(r'\widehat C_m^{\mathrm{TC}}', r'C_m^{\mathrm{formula}}')
        value = value.replace(r'\;\uparrow', '').replace('\\;\u2191','')
        if analytic:
            value = value.replace(r'\widehat{p}_h', 'p_h')
        obj.set_text(value)
        if value == 'core' and obj.axes is not None:
            lo, hi = obj.axes.get_xlim()
            x, y = obj.get_position()
            obj.set_position((x + .035*(hi-lo), y))
    return save_original(fig, *args, **kwargs)


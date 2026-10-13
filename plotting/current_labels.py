"""Current full method names and pdf label, with geometry assertions."""
import numpy as np
from matplotlib.legend import Legend


def geometry(fig):
    arrays = []
    for ax in fig.axes:
        arrays.extend([np.asarray(ax.get_xlim()), np.asarray(ax.get_ylim())])
        for line in ax.lines:
            arrays.extend([np.asarray(line.get_xdata()), np.asarray(line.get_ydata())])
        for collection in ax.collections:
            arrays.append(np.asarray(collection.get_offsets()))
            arrays.extend(path.vertices.copy() for path in collection.get_paths())
    return [array.copy() for array in arrays]


def change_labels(fig):
    before = geometry(fig)
    top, _, bottom = fig.axes
    targets = [t for t in top.texts if t.get_rotation() == 90]
    assert len(targets) == 1
    targets[0].set_text('pdf')
    legends = [obj for obj in bottom.get_children() if isinstance(obj, Legend)]
    methods = next(leg for leg in legends if [t.get_text() for t in leg.get_texts()] == ['$t$', '$p$', '$bt$'])
    for label, value in zip(methods.get_texts(), ['Student t', 'percentile bootstrap', 'bootstrap t']):
        label.set_text(value)
    methods.set_bbox_to_anchor((.43, .98))
    rlegend = next(leg for leg in legends if leg is not methods)
    for anchor in [.77, .87, .93, .72]:
        rlegend.set_bbox_to_anchor((anchor, .98))
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        a, b = [leg.get_window_extent(renderer) for leg in [methods, rlegend]]
        boundaries = [line.get_transform().transform(line.get_xydata())[0, 0]
                      for line in bottom.lines if str(line.get_color()).lower() == '#e58b31']
        if not a.overlaps(b) and not any(b.x0 <= x <= b.x1 for x in boundaries):
            break
    else:
        raise AssertionError('No clear position for the R legend')
    assert not a.overlaps(b), 'Method and R legends overlap'
    # The source curves use the lower two thirds; legends retain the existing top band.
    for leg, box in [(methods, a), (rlegend, b)]:
        for line in bottom.lines:
            if len(line.get_xdata()) < 20:
                continue
            points = line.get_transform().transform(line.get_xydata())
            assert not np.any((points[:, 0] >= box.x0) & (points[:, 0] <= box.x1)
                              & (points[:, 1] >= box.y0) & (points[:, 1] <= box.y1)), 'Legend covers data'
    after = geometry(fig)
    assert len(before) == len(after)
    assert all(np.array_equal(a, b, equal_nan=True) for a, b in zip(before, after))
    return dict(numerical_arrays_unchanged=len(before), legend_overlap=False, legend_data_overlap=False)


'''Helper for plotting reblocking plots.'''

# copyright: (c) 2015 James Spencer
# license: modified BSD license; see LICENSE for further details.

import matplotlib.pyplot as plt
import mpl_toolkits
import mpl_toolkits.axisartist as mpl_aa
import matplotlib.tight_layout as mpl_tl

def plot_reblocking(block_info, plotfile=None, plotshow=True):
    '''Plot the reblocking data.

Parameters
----------
block_info : :class:`pandas.DataFrame`
    Reblocking data (i.e. the first item of the tuple returned by ``reblock``).
plotfile : string
    If not null, save the plot to the given filename.  If '-', then show the
    plot interactively.  See also ``plotshow``.
plotshow : bool
    If ``plotfile`` is not given or is '-', then show the plot interactively.

Returns
-------
fig : :class:`matplotlib.figure.Figure`
    plot of the reblocking data.
'''

    # See http://matplotlib.org/examples/axes_grid/demo_parasite_axes2.html.

    # Create host axes.  Must plot to here first and then clone...
    fig = plt.figure()
    host = mpl_toolkits.axes_grid1.host_subplot(111, axes_class=mpl_aa.Axes)

    offset = -90 # distance between y axes.
    axes = []
    for (i, col) in enumerate(block_info.columns.get_level_values(0).unique()):
        if i == 0:
            ax = host
        else:
            ax = host.twinx()
            # Create a new y axis a little to the left of the current figure.
            new_fixed_axis = ax.get_grid_helper().new_fixed_axis
            ax.axis["left"] = new_fixed_axis(loc='left', axes=ax,
                                             offset=(i*offset,0))
            # Hide right-hand side ticks as the multiple y axes just interfere
            # with each other.
            ax.axis["right"].toggle(all=False)

        block = block_info.index.values
        std_err = block_info.ix[:,(col, 'standard error')].values
        if 'standard error error' in block_info[col]:
            std_err_err = block_info.ix[:,(col, 'standard error error')].values
        else:
            std_err_err = 0*std_err
        line = ax.errorbar(block, std_err, std_err_err, marker='o', label=col)

        # There should only be (at most) one non-null value for optimal block.
        opt = block_info.xs((col, 'optimal block'), axis=1, copy=False)
        opt = opt[opt != ''].index.values
        if opt:
            opt = opt[0]
            ax.annotate('', (block[opt], std_err[opt]-std_err_err[opt]),
                         xytext=(0, -30), textcoords='offset points',
                         arrowprops=dict(
                             arrowstyle="->", color=line[0].get_color()
                       ),)

        ax.set_ylabel('%s standard error' % col, labelpad=0)
        ax.set_xlim((-0.1,len(block)-0.9))
        axes.append(ax)

    plt.xlabel('Reblock iteration')
    plt.legend(loc=2)

    plt.tight_layout()

    if plotfile == '-' or (not plotfile and plotshow):
        plt.show()
    elif plotfile:
        plt.savefig(plotfile)

    return fig

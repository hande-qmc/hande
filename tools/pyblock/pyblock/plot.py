'''Helper for plotting reblocking plots.'''

# copyright: (c) 2015 James Spencer
# license: modified BSD license; see LICENSE for further details.

import matplotlib.pyplot as plt

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

    fig = plt.figure()
    data_sets = block_info.columns.get_level_values(0).unique()

    for (i, col) in enumerate(data_sets):

        ax = fig.add_subplot(len(data_sets), 1, i+1)

        # There should only be (at most) one non-null value for optimal block.
        opt = block_info[block_info[(col,'optimal block')] != ''].index.values
        if opt:
            opt = opt[0]

        std_err = block_info[(col, 'standard error')]
        if 'standard error error' in block_info[col]:
            std_err_err = block_info[(col, 'standard error error')]
        else:
            std_err_err = 0*std_err
        line = ax.errorbar(block_info.index, std_err, std_err_err, marker='o',
                           label=col)

        if opt:
            ax.annotate('', (opt, std_err[opt]-std_err_err[opt]),
                         xytext=(0, -20), textcoords='offset points',
                         arrowprops=dict(
                             arrowstyle="->", #color=line[0].get_color()
                             linewidth=1.2*line[0].get_linewidth(),
                       ),)

        ax.legend(loc=2)
        ax.set_ylabel('standard error')
        ax.set_xlabel('Reblock iteration')

    size = fig.get_size_inches()
    fig.set_size_inches(size[0], size[1]*len(data_sets))
    fig.tight_layout()

    if plotfile == '-' or (not plotfile and plotshow):
        plt.show()
    elif plotfile:
        fig.savefig(plotfile)

    return fig

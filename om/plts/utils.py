"""Utilities for plotting."""

from os.path import join as pjoin

###################################################################################################
###################################################################################################

def set_axis_spines(ax, lw=None):
    """Helper function for setting axis spines"""

    # Set the top and right side frame & ticks off
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Set linewidth of remaining spines
    if lw:
        ax.spines['left'].set_linewidth(lw)
        ax.spines['bottom'].set_linewidth(lw)


def save_figure(save_out, save_name, fig_info):
    """Helper function for saving out plots."""

    if save_out:
        plt.savefig(save_name, format=f_info.format, bbox_inches=f_info.bbox, dpi=f_info.dpi)
        save_name = pjoin(fig_info.save_path, save_name + fig_info.format)

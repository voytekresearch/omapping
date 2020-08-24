"""Fig info defines a class to store figure settings for plots created in omegamappin."""

###################################################################################################
###################################################################################################

class FigInfo(object):
    """Object to hold settings to save figures.

    Attributes
    ----------
    t_fs : int
        Font size for figure title.
    sp_fs : int
        Font size for figure sub-plots.
    ax_fs : int
        Font size for axis labels.
    ti_fs : int
        Font size for ticks.
    ax_lw : float
        Line width.
    add_title : boolean
        Whether to add titles.
    title : str
        Title to add to plots.
    vis_opac : float
        Opacity to use for figures.
    save_path : str
        Path to save out figures to.
    format : {'svg', 'pdf'}
        Format to save out figure as.
    bbox : {'tight'}
        Setting for boundary box.
    dpi : int
        DPI to save out the figure with.
    """

    def __init__(self):

        # Default Settings - font sizes
        self.t_fs = 36         # Title font size
        self.sp_fs = 28        # Subplot title font size
        self.ax_fs = 24        # Axis font size
        self.ti_fs = 16        # Ticks font size

        # Default Settings - other settings
        self.ax_lw = 3

        # Default Settings - what to add to plot
        self.add_title = False

        # Plot Information
        self.title = 'Group'
        self.vis_opac = 0.005

        # Save Information
        self.save_path = ("/Users/tom/Documents/Research/"
                          "1-Projects/MEGmapping/4-Figures/MegAnalysis/")
        self.format = 'pdf'
        self.bbox = 'tight'
        self.dpi = 600

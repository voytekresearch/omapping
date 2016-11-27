

##
##
##

class FigInfo(object):
    """Object to hold settings to save figures.

    Attributes
    ----------
    t_fs : int
        Font size for figure title.
    sp_fs : int
        xx
    ax_fs : int
        xx
    ti_fs : int
        xx
    ax_lw : float
        Line width.
    add_title : boolean
        Whether to add titles
    title : str
        xx
    vis_opac : float
        xx
    save_path : str
        xx
    format : {'.svg', '.pdf'}
        Format to save out figure as.
    bbox : {'tight'}
        Setting for ...
    dpi : int
        DPI to save out the figure with.
    """

    def __init__(self):

        # Default Settings - font sizes
        self.t_fs = 22         # Title font size
        self.sp_fs = 20        # Subplot title font size
        self.ax_fs = 20        # Axis font size
        self.ti_fs = 14        # Ticks font size

        # Default Settings - other settings
        self.ax_lw = 2.5

        # Default Settings - what to add to plot
        self.add_title = False

        # Plot Information
        self.title = 'Group'
        self.vis_opac = 0.005

        # Save Information
        self.save_path = ("/Users/thomasdonoghue/Documents/Research/"
                          "1-Projects/OMEGA/4-Figures/MegAnalysis/")
        self.format = 'svg'
        self.bbox = 'tight'
        self.dpi = 150
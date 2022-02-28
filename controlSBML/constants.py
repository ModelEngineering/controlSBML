"""Constants for Project."""

from docstring_expander.kwarg import Kwarg
import matplotlib.pyplot
import os


# Legend specification
class LegendSpec():

    def __init__(self, names, crd=(1.15, 1), loc="upper right"):
        """

        Parameters
        ----------
        names: list-srs - names of the legends
        crd: (float, float) - coordinate of the legend
        loc: str - position of the legend
        """
        self.names = list(names)
        self.crd = crd
        self.loc = loc


################ DIRECTORIES #################
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
for _ in range(1):
  PROJECT_DIR = os.path.dirname(PROJECT_DIR)
CODE_DIR = os.path.join(PROJECT_DIR, "controlSBML")
TEST_DIR = os.path.join(PROJECT_DIR, "tests")
DATA_DIR = os.path.join(PROJECT_DIR, "data")
BIOMODELS_ZIP_FILENAME = "biomodels.zip"

# Constants
TIME = "time"
START_TIME = 0  # Default start time
END_TIME = 5  # Default endtime
POINTS_PER_TIME = 10

# Keyword options
O_AX = "ax"
O_END_TIME = "end_time"
O_FIGURE = "figure"
O_FIGSIZE = "figsize"
O_IS_PLOT = "is_plot"
O_LEGEND_CRD = "legend_crd"  # Legend coordinate
O_LEGEND_SPEC = "legend_spec"
O_POINTS_PER_TIME = "points_per_time"
O_START_TIME = "start_time"
O_SUPTITLE = "suptitle"
O_TITLE = "title"
O_XLABEL = "xlabel"
O_XLIM = "xlim"
O_XTICKLABELS = "xticklabels"
O_YLABEL = "ylabel"
O_YLIM = "ylim"
O_YTICKLABELS = "yticklabels"

# Default values of options
SIM_DCT = dict(
      start_time=START_TIME,  # Start time of the simulation
      end_time=END_TIME,      # End time of the simulation
      points_per_time=POINTS_PER_TIME,    # Number of points in the simulation
      )
# Options for a single plot
PLOT_DCT = dict(
      legend_spec=None,     # LegendSpec
      legend_crd=None,      # coordinate for legend
      ylim=None,           # maximum and minimum value of y
      xlim=None,
      xlabel="",           
      ylabel="",           
      title="",             # plot title
      ax=None,              # axis to plot
      xticklabels=None,
      yticklabels=None,
      )
# Options for the full figure
FIG_DCT = dict(
      is_plot=True,         # Is a figure generated
      figure=None,          # Figure option
      figsize=(10, 10),     # Size of the figure
      suptitle="",          # Title for the figure
      )

DEFAULT_DCTS = [PLOT_DCT, FIG_DCT, SIM_DCT]


# Must maintain this in correspondence with SIM_DCT, PLOT_DCT, FIG_DCT
KWARGS = [
    #SIMULATION OPTIONS
    Kwarg(O_END_TIME, default=10, dtype=float, doc="end time of simulation"),
    Kwarg(O_POINTS_PER_TIME, default=10, dtype=float,
          doc="number of simulation points per time period"),
    Kwarg(O_START_TIME, default=0, dtype=float, doc="when simulation begins"),
    #PLOT OPTIONS
    Kwarg(O_AX, default=None, dtype=matplotlib.pyplot.axes, doc="Plotting axis"),
    Kwarg(O_LEGEND_SPEC, default=None, dtype=LegendSpec, doc="Position of the legend"),
    Kwarg(O_LEGEND_CRD, default=None, dtype=tuple, doc="Coordinate position of the legend"),
    Kwarg(O_TITLE, default="", dtype=str, doc="Plot title"),
    Kwarg(O_XLABEL, default="", dtype=str, doc="x-axis label"),
    Kwarg(O_XLIM, default=None, dtype=(float, float), doc="Lower and upper values of x axis"),
    Kwarg(O_XTICKLABELS, default=None, dtype=list, doc="x-axis tic marks"),
    Kwarg(O_YLIM, default=None, dtype=(float, float), doc="Lower and upper values of y axis"),
    Kwarg(O_YLABEL, default="", dtype=str, doc="y-axis label"),
    Kwarg(O_YTICKLABELS, default=None, dtype=list, doc="y-axis tic marks"),
    # FIGURE OPTIONS
    Kwarg(O_FIGURE, default=None, dtype=matplotlib, doc="Figure option"),
    Kwarg(O_FIGSIZE, default=None, dtype=(float, float), doc="widith, height"),
    Kwarg(O_IS_PLOT, default=True, dtype=bool, doc="Do the plot"),
    Kwarg(O_SUPTITLE, default="", dtype=str, doc="Figure title"),
    ]
SIM_KWARGS = list(SIM_DCT.keys())
PLOT_KWARGS = list(PLOT_DCT.keys())
FIG_KWARGS = list(FIG_DCT.keys())
ALL_KWARGS = []
for kwargs in [SIM_KWARGS, PLOT_KWARGS, FIG_KWARGS]:
    ALL_KWARGS.extend(kwargs)

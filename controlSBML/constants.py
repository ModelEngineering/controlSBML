"""Constants for Project."""

from docstring_expander.kwarg import Kwarg
import matplotlib.pyplot
import os


# Legend specification
class LegendSpec():

    def __init__(self, names, crd=(1.4, 1), loc="upper right"):
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
O_FIGSIZE = "figsize"
O_IS_PLOT = "is_plot"
O_LEGEND_SPEC = "legend_spec"
O_POINTS_PER_TIME = "points_per_time"
O_START_TIME = "start_time"
O_SUPTITLE = "suptitle"
O_TITLE = "title"
O_XLABEL = "xlabel"
O_XTICKLABELS = "xticklabels"
O_YTICKLABELS = "yticklabels"
O_YLABEL = "ylabel"
O_YLIM = "ylim"

# Default values of options
SIM_DCT = dict(
      start_time=START_TIME,  # Start time of the simulation
      end_time=END_TIME,      # End time of the simulation
      points_per_time=POINTS_PER_TIME,    # Number of points in the simulation
      )
# Options for a single plot
PLOT_DCT = dict(
      ylim=None,           # maximum and minimum value of y
      xlim=None,
      xlabel="",           
      ylabel="",           
      title="",             # plot title
      legend_spec=None,     # LegendSpec
      ax=None,              # axis to plot
      xticklabels=None,
      yticklabels=None,
      )
# Options for the full figure
FIG_DCT = dict(
      is_plot=True,         # Is a figure generated
      figsize=(10, 10),     # Size of the figure
      suptitle="",          # Title for the figure
      )

DEFAULT_DCTS = [PLOT_DCT, FIG_DCT, SIM_DCT]


# Must maintain this in correspondence with SIM_DCT, PLOT_DCT, FIG_DCT
KWARGS = [
    #SIMULATION OPTIONS
    Kwarg("end_time", default=10, dtype=float, doc="end time of simulation"),
    Kwarg("points_per_time", default=10, dtype=float, doc="number of simulation points per time period"),
    Kwarg("start_time", default=0, dtype=float, doc="when simulation begins"),
    #PLOT OPTIONS
    Kwarg("ax", default=None, dtype=matplotlib.pyplot.axes, doc="Plotting axis"),
    Kwarg("legend_spec", default=None, dtype=LegendSpec, doc="position of the legend"),
    Kwarg("title", default="", dtype=str, doc="Plot title"),
    Kwarg("xlabel", default="", dtype=str, doc="x-axis label"),
    Kwarg("xlim", default=None, dtype=(float, float), doc="lower and upper values of x axis"),
    Kwarg("xticklabels", default=None, dtype=list, doc="x-axis tic marks"),
    Kwarg("ylim", default=0, dtype=(float, float), doc="lower and upper values of y axis"),
    Kwarg("ylabel", default="", dtype=str, doc="y-axis label"),
    Kwarg("yticklabels", default=None, dtype=list, doc="y-axis tic marks"),
    # FIGURE OPTIONS
    Kwarg("figsize", default=(10, 10), dtype=(float, float), doc="widith, height"),
    Kwarg("is_plot", default=True, dtype=bool, doc="do the plot"),
    Kwarg("suptitle", default=0, dtype=str, doc="figure title"),
    ]
SIM_KWARGS = list(SIM_DCT.keys())
PLOT_KWARGS = list(PLOT_DCT.keys())
FIG_KWARGS = list(FIG_DCT.keys())
ALL_KWARGS = []
for kwargs in [SIM_KWARGS, PLOT_KWARGS, FIG_KWARGS]:
    ALL_KWARGS.extend(kwargs)

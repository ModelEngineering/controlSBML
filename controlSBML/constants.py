"""Constants for Project."""

from docstring_expander.kwarg import Kwarg
import matplotlib.pyplot
import pandas as pd
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
A_DF = None  # Use default value of A matrix
B_DF = None  # Use default value of B matrix
C_DF = None  # Use default value of C matrix
END_TIME = 5  # Default endtime
EVENT = "event"
INPUT = "input"
OUTPUT = "output"
OUT_STATE = "out_state"
PARAMS = "params"
POINTS_PER_TIME = 10
START_TIME = 0  # Default start time
STATE = "state"
STEP_VAL = 1  # Multiplier used for simulation input
TIME = "time"

# Keyword options
O_A_DF = "A_df"
O_AX = "ax"
O_B_DF = "B_df"
O_C_DF = "C_df"
O_END_TIME = "end_time"
O_FIGURE = "figure"
O_FIGSIZE = "figsize"
O_STEP_VAL = "step_val"
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
      A_df=A_DF,
      B_df=B_DF,
      C_df=C_DF,
      step_val=STEP_VAL,
      start_time=START_TIME,
      end_time=END_TIME,
      points_per_time=POINTS_PER_TIME,
      )
# Options for a single plot
PLOT_DCT = dict(
      legend_spec=None,
      legend_crd=None,
      ylim=None,
      xlim=None,
      xlabel="",           
      ylabel="",           
      title="",
      ax=None,
      xticklabels=None,
      yticklabels=None,
      )
# Options for the full figure
FIG_DCT = dict(
      is_plot=True,
      figure=None,
      figsize=(10, 10),
      suptitle="",
      )

DEFAULT_DCTS = [PLOT_DCT, FIG_DCT, SIM_DCT]


# Must maintain this in correspondence with SIM_DCT, PLOT_DCT, FIG_DCT
KWARGS = [
    #SIMULATION OPTIONS
    Kwarg(O_A_DF, default=10, dtype=pd.DataFrame, doc="Linear system A matrix"),
    Kwarg(O_AX, default=None, dtype=matplotlib.pyplot.axes, doc="Plotting axis"),
    Kwarg(O_B_DF, default=10, dtype=pd.DataFrame, doc="Linear system B matrix"),
    Kwarg(O_C_DF, default=10, dtype=pd.DataFrame, doc="Linear system C matrix"),
    Kwarg(O_END_TIME, default=10, dtype=float, doc="end time of simulation"),
    Kwarg(O_POINTS_PER_TIME, default=10, dtype=float,
          doc="number of simulation points per time period"),
    Kwarg(O_START_TIME, default=0, dtype=float, doc="when simulation begins"),
    #PLOT OPTIONS
    Kwarg(O_LEGEND_SPEC, default=None, dtype=LegendSpec,
          doc="Position of the legend"),
    Kwarg(O_LEGEND_CRD, default=None, dtype=tuple,
          doc="Coordinate position of the legend"),
    Kwarg(O_STEP_VAL, default=10, dtype=float, doc="value of step input"),
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

# TimeSeries
MS_IN_SEC = 1000
SEC_IN_MS = 1.0/MS_IN_SEC
TIMESERIES_INDEX_NAME = "miliseconds"

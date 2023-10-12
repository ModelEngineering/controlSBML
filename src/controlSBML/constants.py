"""Constants for Project."""
import collections
from docstring_expander.kwarg import Kwarg
import matplotlib.pyplot
import numpy as np
import pandas as pd
import os


DEFAULT_NUM_STEP = 5

# Classes
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

FitterResult = collections.namedtuple("FitterResult",
            "transfer_function staircase_arr staircase_name parameters rms_residuals stderr nfev redchi time_series input_name output_name")


################ DIRECTORIES #################
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
for _ in range(2):
  PROJECT_DIR = os.path.dirname(PROJECT_DIR)
CODE_DIR = os.path.join(PROJECT_DIR, "controlSBML")
PLOT_DIR = PROJECT_DIR
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
IN = "in"
IS_PLOT = False
OUT = "out"
OUTPUT = "output"
OUT_STATE = "out_state"
PARAMS = "params"
POINTS_PER_TIME = 10
START_TIME = 0  # Default start time
STATE = "state"
STEP_VAL = 1  # Multiplier used for simulation input
TIME = "time"

# Keyword options
O_AX = "ax"
O_AX2 = "ax2"
O_END_TIME = "end_time"
O_FIGURE = "figure"
O_FIGSIZE = "figsize"
O_STEP_VAL = "step_val"
O_IS_PLOT = "is_plot"
O_LEGEND_CRD = "legend_crd"  # Legend coordinate
O_LEGEND_SPEC = "legend_spec"
O_POINTS_PER_TIME = "points_per_time"
O_PREDICTED = "predicted"
O_STAIRCASE = "staircase"
O_START_TIME = "start_time"
O_SUPTITLE = "suptitle"
O_TITLE = "title"
O_WRITEFIG = "writefig"  # Write the figure
O_XLABEL = "xlabel"
O_XLIM = "xlim"
O_XTICKLABELS = "xticklabels"
O_YLABEL = "ylabel"
O_YLIM = "ylim"
O_YTICKLABELS = "yticklabels"

# Default values of options
SIM_DCT = dict(
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
      ax2=None,
      xticklabels=None,
      yticklabels=None,
      )
# Options for the full figure
FIG_DCT = dict(
      is_plot=True,
      figure=None,
      figsize=(10, 10),
      suptitle="",
      writefig=False,
      )

DEFAULT_DCTS = [PLOT_DCT, FIG_DCT, SIM_DCT]


# Must maintain this in correspondence with SIM_DCT, PLOT_DCT, FIG_DCT
KWARGS = [
    #SIMULATION OPTIONS
    Kwarg(O_AX, default=None, dtype=matplotlib.pyplot.axes, doc="Plotting axis"),
    Kwarg(O_AX2, default=None, dtype=matplotlib.pyplot.axes, doc="Plotting 2nd axis"),
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
    Kwarg(O_IS_PLOT, default=True, dtype=bool, doc="Display the plot"),
    Kwarg(O_WRITEFIG, default=False, dtype=bool, doc="Write the plot"),
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
TIMES = np.linspace(0, 5, 50)

# URLs
WOLF_URL = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000206.2?filename=BIOMD0000000206_url.xml"
METFORMIN_URL = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000001039.5?filename=Zake2021_Metformin%2BMice%2BIV.xml"
MTOR_URL = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000823.2?filename=Varusai2018.xml"

# Transfer function building
DEFAULT_NUM_NUMERATOR = 2
DEFAULT_NUM_DENOMINATOR = 3
DEFAULT_DEVIATION = 0.5

# Plot line colors
PREDICTED_COLOR = "blue"
SIMULATED_COLOR = "green"
INPUT_COLOR = "red"

COMMENT_STR = "//"
DEFAULT_NUM_STEP = 5
DEFAULT_INITIAL_VALUE = 0
DEFAULT_FINAL_VALUE = 10
DEFAULT_POINT_PER_STEP = 10

# Types
TYPE_PARAMETER = "parameter"
TYPE_BOUNDARY_SPECIES = "boundary_species"
TYPE_FLOATING_SPECIES = "floating_species"
TYPE_REACTION = "reaction"
TYPE_ASSIGNMENT = "assignment"   # assignment rule
ANTIMONY_TYPES = [TYPE_PARAMETER, TYPE_BOUNDARY_SPECIES, TYPE_FLOATING_SPECIES, TYPE_REACTION, TYPE_ASSIGNMENT]
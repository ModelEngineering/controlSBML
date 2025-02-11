"""Constants for Project."""
import collections
from docstring_expander.kwarg import Kwarg # type: ignore
import matplotlib.pyplot
import numpy as np
import pandas as pd # type: ignore
import os


##########################
# Classes
##########################
class NoiseSpec(object):
    # Specifications for a noise model
    def __init__(self, sine_amp=0, sine_freq=0, random_mag=0, random_std=0, offset=None, slope=0):
        """_summary_

        Args:
            sine_amp (int, optional): Amplitude of the sine wave. Defaults to 0.
            sine_freq (int, optional): Frequency of the sine wave. Defaults to 0.
            random_mag (int, optional): Magnitude scalinging of the log-normal random variable. Defaults to 0.
            random_std (int, optional): Standard deviation of the normal for the log-normal random variable added to the function. Defaults to 0.
            offset (int, optional): Constant offset. Defaults to  sine_amp.
            slope (int, optional): Slope w.r.t. time. Defaults to 0.
        """
        self.sine_amp = sine_amp
        self.sine_freq = sine_freq
        self.random_mag = random_mag
        self.random_std = random_std
        if offset is None:
            if (np.isclose(sine_amp, 0) and np.isclose(slope, 0)):
                offset = sine_amp
            else:
                offset = 0
        self.offset = offset
        self.slope = slope

    def __str__(self):
        result = f"NoiseSpec(sine_amp={self.sine_amp}, sine_freq={self.sine_freq}, random_mag={self.random_mag}"
        result += f" random_std={self.random_std}, dc_gain={self.offset}, slope={self.slope})"
        return result
    
    def __eq__(self, other):
        return self.__dict__ == other.__dict__

class DisturbanceSpec(NoiseSpec):
    # Specifications for a disturbance model
    pass


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

class FitterResult(object):
    ATTRS = ["transfer_function", "staircase_arr", "staircase_name", "parameters", "rms_residuals",
             "stderr", "nfev", "redchi", "time_series", "input_name", "output_name", "antimony_builder"]

    def __init__(self, transfer_function=None, staircase_arr=None, staircase_name=None, parameters=None, rms_residuals=None,
                 stderr=None, nfev=None, redchi=None, time_series=None, input_name=None, output_name=None, antimony_builder=None):
        self.transfer_function = transfer_function
        self.staircase_arr = staircase_arr
        self.staircase_name = staircase_name
        self.parameters = parameters
        self.rms_residuals = rms_residuals
        self.stderr = stderr
        self.nfev = nfev
        self.redchi = redchi
        self.time_series = time_series
        self.input_name = input_name
        self.output_name = output_name
        self.antimony_builder = antimony_builder

    def copy(self):
        kwargs = {a: getattr(self, a) for a in self.ATTRS}
        return FitterResult(**kwargs)
    
    def equals(self, other):
        for attr in self.ATTRS:
            if getattr(self, attr) != getattr(other, attr):
                return False
        return True


################ DIRECTORIES #################
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
for _ in range(2):
  PROJECT_DIR = os.path.dirname(PROJECT_DIR)
CODE_DIR = os.path.join(PROJECT_DIR, "controlSBML")
PLOT_DIR = PROJECT_DIR
TEST_DIR = os.path.join(PROJECT_DIR, "tests")
MODEL_DIR = os.path.join(PROJECT_DIR, "models")
NOTEBOOK_DIR = os.path.join(PROJECT_DIR, "notebooks")
DATA_DIR = os.path.join(PROJECT_DIR, "data")
BIOMODELS_ZIP_FILENAME = "biomodels.zip"
MTOR_PATH = os.path.join(MODEL_DIR, "Varusai2018.xml")
WOLF_PATH = os.path.join(MODEL_DIR, "BIOMD0000000206.xml")

# Constants
END_TIME = 5.0  # Default endtime
EVENT = "event"
FITTER_METHOD = "fitter_method"
INPUT = "input"
IN = "in"
IS_PLOT = "is_plot"
MAX = "max"
MIN = "min"
MSE = "mse"
OUT = "out"
OUTPUT = "output"
OUT_STATE = "out_state"
PARAMS = "params"
POINTS_PER_TIME = 10
REASON = "reason"
SCORE = "score"
START_TIME = 0.0  # Default start time
STATE = "state"
STEP_VAL = 1  # Multiplier used for simulation input
TIME = "time"

# Keyword options
O_AX = "ax"
O_AX2 = "ax2"
O_DISTURBANCE_SPEC = "disturbance_spec"
O_END_TIME = "end_time"
O_FIGURE = "figure"
O_FIGSIZE = "figsize"
O_FINAL_VALUE = "final_value"
O_FITTER_METHOD = "fitter_method"
O_KP_SPEC = "kP_spec"
O_KI_SPEC = "kI_spec"
O_KD_SPEC = "kD_spec"
O_KF_SPEC = "kF_spec"
O_INITIAL_VALUE = "initial_value"
O_INPUT_NAME = "input_name"
O_INPUT_NAMES = "input_names"
O_IS_GREEDY = "is_greedy"
O_IS_FIXED_INPUT_SPECIES = "is_fixed_input_species"
O_IS_STEADY_STATE = "is_steady_state"
O_IS_PLOT = "is_plot"
O_LEGEND = "legend"
O_LEGEND_CRD = "legend_crd"  # Legend coordinate
O_LEGEND_SPEC = "legend_spec"
O_MARKERS = "markers"
O_MARKERSIZE = "markersize"
O_NOISE_SPEC = "noise_spec"
O_NUM_POINT = "num_point"
O_NUM_PROCESS = "num_process"
O_NUM_RESTART = "num_restart"
O_NUM_STEP = "num_step"
O_OUTPUT_NAME = "output_name"
O_OUTPUT_NAMES = "output_names"
O_POINTS_PER_TIME = "points_per_time"
O_PREDICTED = "predicted"
O_SELECTIONS = "selections"   # Selections for the plot
O_SETPOINT = "setpoint"
O_SIGN = "sign"
O_STAIRCASE = "staircase"
O_START_TIME = "start_time"
O_STEP_VAL = "step_val"
O_SUPTITLE = "suptitle"
O_TIMES = "times"
O_TITLE = "title"
O_WRITEFIG = "writefig"  # Write the figure
O_XLABEL = "xlabel"
O_XLABEL_ANGLE = "xlabel_angle"
O_XLIM = "xlim"
O_XTICKLABELS = "xticklabels"
O_YLABEL = "ylabel"
O_YLIM = "ylim"
O_YTICKLABELS = "yticklabels"

# Control parameters
CP_KP = "kP"
CP_KI = "kI"
CP_KF = "kF"
CP_KD = "kD"
CONTROL_PARAMETERS = [CP_KP, CP_KI, CP_KD, CP_KF]
KP_SPEC = "kP_spec"
KI_SPEC = "kI_spec"
KD_SPEC = "kD_spec"
KF_SPEC = "kF_spec"
CONTROL_PARAMETER_SPECS = [KP_SPEC, KI_SPEC, KD_SPEC, KF_SPEC]

# URLs
WOLF_URL = "https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL3352181362/2/BIOMD0000000206_url.xml"
METFORMIN_URL = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000001039.5?filename=Zake2021_Metformin%2BMice%2BIV.xml"
MTOR_URL = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000823.2?filename=Varusai2018.xml"

# Transfer function building
DEFAULT_NUM_ZERO = 1
DEFAULT_NUM_POLE = 3
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

# TimeSeries
MS_IN_SEC = 1000.0
SEC_IN_MS = 1.0/MS_IN_SEC
TIMESERIES_INDEX_NAME = "miliseconds"
TIMES = np.linspace(DEFAULT_INITIAL_VALUE, DEFAULT_FINAL_VALUE, DEFAULT_POINT_PER_STEP*DEFAULT_FINAL_VALUE)

# Types
TYPE_PARAMETER = "parameter"
TYPE_BOUNDARY_SPECIES = "boundary_species"
TYPE_FLOATING_SPECIES = "floating_species"
TYPE_REACTION = "reaction"
TYPE_ASSIGNMENT = "assignment"   # assignment rule
VALID_INPUT_TYPES = [TYPE_BOUNDARY_SPECIES, TYPE_FLOATING_SPECIES, TYPE_PARAMETER, TYPE_ASSIGNMENT]
VALID_OUTPUT_TYPES = [TYPE_BOUNDARY_SPECIES, TYPE_FLOATING_SPECIES, TYPE_PARAMETER, TYPE_ASSIGNMENT, TYPE_REACTION]
ANTIMONY_TYPES = list(set(VALID_INPUT_TYPES + VALID_OUTPUT_TYPES))

# Staircase
DEFAULT_NUM_STEP = 5
DEFAULT_INITIAL_VALUE = 0
DEFAULT_FINAL_VALUE = 10
DEFAULT_NUM_POINT = 100

# Fitter methods
FITTER_METHOD_GPZ = "gpz"
FITTER_METHOD_POLY = "poly"
DEFAULT_FITTER_METHOD = FITTER_METHOD_POLY


# Default values of options
SIM_DCT = dict(
      step_val=STEP_VAL,
      start_time=START_TIME,
      end_time=END_TIME,
      points_per_time=POINTS_PER_TIME,
      times=TIMES,
      initial_value=DEFAULT_INITIAL_VALUE,
      final_value=DEFAULT_FINAL_VALUE,
      )
# Options for a single plot
PLOT_DCT:dict = dict(
      legend=None,
      legend_spec=None,
      legend_crd=None,
      markers=None,
      markersize=8,
      ylim=None,
      xlim=None,
      xlabel="",
      ylabel="",
      title="",
      ax=None,
      ax2=None,
      xticklabels=None,
      yticklabels=None,
      xlabel_angle=0,
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
    Kwarg(O_MARKERS, default=None, dtype=list, doc="list of markers for plots"),
    Kwarg(O_MARKERSIZE, default=8, dtype=float, doc="size of markers"),
    Kwarg(O_POINTS_PER_TIME, default=10, dtype=float,
          doc="number of simulation points per time period"),
    Kwarg(O_INITIAL_VALUE, default=DEFAULT_INITIAL_VALUE, dtype=float, doc="initial value of step input"),
    Kwarg(O_FINAL_VALUE, default=DEFAULT_FINAL_VALUE, dtype=float, doc="final value of step input"),
    Kwarg(O_TIMES, default=np.linspace(0, 10, 100), dtype=np.array, doc="times of simulations"),
    Kwarg(O_START_TIME, default=0, dtype=float, doc="when simulation begins"),
    #PLOT OPTIONS
    Kwarg(O_LEGEND, default=None, dtype=list,
          doc="List of names for plots"),
    Kwarg(O_LEGEND_SPEC, default=None, dtype=LegendSpec,
          doc="Position of the legend"),
    Kwarg(O_LEGEND_CRD, default=None, dtype=tuple,
          doc="Coordinate position of the legend"),
    Kwarg(O_STEP_VAL, default=10, dtype=float, doc="value of step input"),
    Kwarg(O_TITLE, default="", dtype=str, doc="Plot title"),
    Kwarg(O_XLABEL, default="", dtype=str, doc="x-axis label"),
    Kwarg(O_XLABEL_ANGLE, default=0, dtype=int, doc="x-axis label angle"),
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

# Design results
DESIGN_RESULT_SUCCESS = "Design successful."
DESIGN_RESULT_CANNOT_SIMULATE = "No design. Cannot simulate the closed loop system." 
DESIGN_RESULT_OUTPUT_TOO_SMALL = "No design. Output is too small."
DESIGN_RESULT_OUTPUT_TOO_LARGE = "No design. Output is too large."
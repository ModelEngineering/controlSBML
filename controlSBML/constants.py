"""Constants for Project."""
from controlSBML.options import Options

import collections
import os


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
SIM_OPTS = Options(dict(
      start_time=START_TIME,  # Start time of the simulation
      end_time=END_TIME,      # End time of the simulation
      points_per_time=POINTS_PER_TIME,    # Number of points in the simulation
      ))
# Options for a single plot
PLOT_OPTS = Options(dict(
      ylim=None,           # maximum and minimum value of y
      xlabel="",           
      ylabel="",           
      title="",             # plot title
      legend_spec=None,     # LegendSpec
      ax=None,              # axis to plot
      xticklabels=None,
      yticklabels=None,
      ))
# Options for the full figure
FIG_OPTS = Options(dict(
      is_plot=True,         # Is a figure generated
      figsize=(10, 10),     # Size of the figure
      suptitle="",          # Title for the figure
      ))

OPTS_LST = [PLOT_OPTS, FIG_OPTS, SIM_OPTS]

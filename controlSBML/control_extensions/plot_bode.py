"""Constructs bode plots for MIMO models."""

import controlSBML.constants as cn
from controlSBML.option_management.option_manager import OptionManager

import control
from docstring_expander.expander import Expander
import matplotlib.pyplot as plt
import numpy as np

NUM_FREQ = 40


@Expander(cn.KWARGS, cn.ALL_KWARGS)
def plotBode(sys, names=None, is_magnitude=True, is_phase=True, is_plot=True, **kwargs):
    """
    Constructs bode plots for a MIMO system.

    Parameters
    ----------
    sys: control.system
    names: list-str
        Names of the states
    is_magnitude: bool
        Do magnitude plots
    is_phase: bool
        Do phase plots
    is_plot: bool
        Display plots
    #@expand
    """
    # Calculate magnitudes and phases for al inputs and outputs
    freqs = np.array([np.pi*n/NUM_FREQ for n in range(1, NUM_FREQ+1)])
    all_mag, all_phase, angle = sys.frequency_response(freqs)
    # Construct the matrix of plots
    num_state = np.shape(sys.A)[0]
    if names is None:
        names = [str(n) for n in range(num_state)]
    mgr = OptionManager(kwargs)
    figure, axes = plt.subplots(num_state, num_state,
         figsize=mgr.fig_opts[cn.O_FIGSIZE])
    for idx1 in range(num_state):
        for idx2 in range(num_state):
            # Create the plot data
            arr = np.array([10e-5 if np.isclose(v, 0) else v
                  for v in all_mag[idx1, idx2]])
            db_arr = 20*np.log10(arr)
            db_arr = np.array([v if v > -10 else -80 for v in db_arr])
            log_angle = np.log(np.array(angle))
            # Construct the plot
            new_mgr = mgr.copy()
            ax = axes[idx1, idx2]
            new_mgr.plot_opts[cn.O_AX] = ax
            title ="%s->%s" % (names[idx1], names[idx2])
            new_mgr.plot_opts.set(cn.O_TITLE, default=title)
            new_mgr.setYlim(db_arr, is_override=False)
            ax.plot(log_angle, db_arr)
            ax.plot([log_angle[0], log_angle[-1]], [0, 0],
                   linestyle="dashed", color="grey")
            if idx1 < num_state - 1:
                new_mgr.plot_opts.set(cn.O_XTICKLABELS, override=[])
            else:
                new_mgr.plot_opts.set(cn.O_XLABEL, default="log rads")
            if idx2 == 0:
                new_mgr.plot_opts.set(cn.O_YLABEL, default="db")
            new_mgr.doPlotOpts()
    mgr.doFigOpts()
        

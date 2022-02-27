"""Creates a table of transfer functions for a StateSpace (MIMO) model."""

import controlSBML.constants as cn
from controlSBML.option_management.option_manager import OptionManager

import control
from docstring_expander.expander import Expander
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NUM_FREQ = 100  # Number of frequences in radians/sec


class StateSpaceTF(object):

    def __init__(self, mimo_sys, input_names=None, output_names=None):
        """
        Parameters
        ----------
        mimo_sys: control.StateSpace
            control.StateSpace
        input_names: list-str
            names of the inputs
        output_names: list-str
            names of the outputs
        """
        self.dataframe = self.ss2tf(mimo_sys, input_names=input_names,
              output_names=output_names)

    def __str__(self):
        stgs = ["(input, output)\n\n"]
        indents = "            "
        for inp in self.dataframe.index:
            for out in self.dataframe.columns:
                pfx = "(%s, %s):  " % (inp, out)
                stg = str(self.dataframe.loc[inp, out])[:-1]  # Exclude nl
                stg = stg.replace("\n", "", 1)
                stg = stg.replace("\n", "\n" + indents)
                stgs.append(pfx + stg + "\n")
        return ("\n").join(stgs)

    @staticmethod
    def getSystemShape(sys):
        """
        Provides the number of states, number of inputs, and number of outputs.

        Parameters
        ----------
        sys: control.StateSpace
        
        Returns
        -------
        int: num states
        int: num inputs
        int: num outputs
        """
        return sys.nstates, sys.ninputs, sys.noutputs

    @classmethod
    def ss2tf(cls, mimo_sys, input_names=None, output_names=None):
        """
        Creates a dataframe of transform functions for each combination
        of inputs and outputs.

        Parameters
        ----------
        mimo_sys: control.StateSpace

        Returns
        -------
        DataFrame
            index: input
            column: output
            value: control.TransferFunction
        """
        num_state, num_input, num_output = cls.getSystemShape(mimo_sys)
        A_mat = mimo_sys.A
        if input_names is None:
            input_names = [str(n) for n in range(1, num_input+1)]
        else:
            input_names = list(input_names)
        if output_names is None:
            output_names = [str(n) for n in range(1, num_output+1)]
        else:
            output_names = list(output_names)
        # Construct matrices for a 1-input, 1-output state space model
        B_base_mat = np.reshape(np.repeat(0, num_state), (num_state, 1))
        C_base_mat = np.reshape(np.repeat(0, num_state), (1, num_state))
        D_mat = np.repeat(0, 1)
        # Construct the dataframe entries
        dct = {n: [] for n in output_names}
        for inp_idx in range(num_input):  # Index for output state
            for out_idx in range(num_output):  # Index for input sate
                # Construct the SISO system
                B_mat = np.array(B_base_mat)
                B_mat[inp_idx, 0] = 1
                C_mat = np.array(C_base_mat)
                C_mat[0, out_idx] = 1
                new_sys = control.StateSpace(A_mat, B_mat, C_mat, D_mat)
                siso_tf = control.ss2tf(new_sys)
                dct[output_names[out_idx]].append(siso_tf)
        #
        df = pd.DataFrame(dct, index=input_names)
        df.index.name = "Inputs"
        df.columns.name = "Outputs"
        return df


    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotBode(self, is_magnitude=True, is_phase=True, is_plot=True, **kwargs):
        """
        Constructs bode plots for a MIMO system. This is done by constructing n*n
        SISO systems where there n states.
    
        Parameters
        ----------
        is_magnitude: bool
            Do magnitude plots
        is_phase: bool
            Do phase plots
        is_plot: bool
            Display plots
        #@expand
        """
        # Calculate magnitudes and phases for al inputs and outputs
        freqs = np.array([n/NUM_FREQ for n in range(1, NUM_FREQ+1)])
        mgr = OptionManager(kwargs)
        figure, axes = plt.subplots(num_state, num_state,
             figsize=mgr.fig_opts[cn.O_FIGSIZE])
        for inp_name in self.dataframe.index:
            for out_name in self.dataframe.columns:
                siso_tf = self.dataframe.loc[inp_name, out_name]
                # Create the plot data
                mag_arr, phase_arr, angle_arr = control.bode(siso_tf, freqs)
                db_mag_arr = 20*np.log10(mag_arr)
                log_angle_arr = np.log10(angle_arr)
                # Construct the plot
                new_mgr = mgr.copy()
                ax = axes[inp_idx, out_idx]
                new_mgr.plot_opts[cn.O_AX] = ax
                title ="%s->%s" % (names[inp_idx], names[out_idx])
                new_mgr.plot_opts.set(cn.O_TITLE, default=title)
                new_mgr.setYlim(db_arr, is_override=False)
                ax.plot(log_angle_arr, db_mag_arr)
                ax.plot([log_angle_arr[0], log_angle_arr[-1]], [0, 0],
                       linestyle="dashed", color="grey")
                if inp_idx < num_state - 1:
                    new_mgr.plot_opts.set(cn.O_XTICKLABELS, override=[])
                else:
                    new_mgr.plot_opts.set(cn.O_XLABEL, default="log rads")
                if out_idx == 0:
                    new_mgr.plot_opts.set(cn.O_YLABEL, default="db")
                new_mgr.doPlotOpts()
        mgr.doFigOpts()
        

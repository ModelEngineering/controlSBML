"""Creates a table of transfer functions for a StateSpace (MIMO) model."""

import controlSBML.constants as cn
import controlSBML as ctl
from controlSBML.option_management.option_manager import OptionManager

import control
from docstring_expander.expander import Expander
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NUM_FREQ = 100  # Number of frequences in radians/sec
FREQ_RNG = [1e-2, 1e2]
LOW = 0
HIGH = 1


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
        self.input_names = list(self.dataframe.columns)
        self.output_names = list(self.dataframe.index)
        self.num_state, self.num_input, self.num_output = self.getSystemShape(
              mimo_sys)

    def __str__(self):
        stgs = ["(input, output)\n\n"]
        indents = "            "
        for inp in self.dataframe.columns:
            for out in self.dataframe.index:
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
            column: input
            index: output
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
        for out_idx, output_name in enumerate(output_names):
            for inp_idx, input_name in enumerate(input_names):
                # Construct the SISO system
                input_name = input_names[inp_idx]
                B_mat = np.array(B_base_mat)
                B_mat[inp_idx, 0] = 1
                C_mat = np.array(C_base_mat)
                C_mat[0, out_idx] = 1
                new_sys = control.StateSpace(A_mat, B_mat, C_mat, D_mat)
                siso_tf = control.ss2tf(new_sys)
                dct[output_name].append(siso_tf)
        #
        df = pd.DataFrame(dct).transpose()
        df.index = output_names
        df.columns = input_names
        df.index.name = "Outputs"
        df.columns.name = "Inputs"
        return df

    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotBode(self, is_magnitude=True, is_phase=True, **kwargs):
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
        freq_arr = np.array(range(NUM_FREQ))
        delta = (FREQ_RNG[HIGH] - FREQ_RNG[LOW])/NUM_FREQ
        freq_arr = (freq_arr + FREQ_RNG[LOW])*delta
        mgr = OptionManager(kwargs)
        legends = []
        for out_idx, out_name in enumerate(self.output_names):
            for inp_idx, inp_name in enumerate(self.input_names):
                siso_tf = self.dataframe.loc[out_name, inp_name]
                # Create the plot data
                _ = control.bode(siso_tf, freq_arr)
                # Construct the plot
                legend =" %s->%s" % (inp_name, out_name)
                legends.append(legend)
        mgr.plot_opts.set(cn.O_LEGEND_SPEC, default=ctl.LegendSpec(legends,
              crd=mgr.plot_opts[cn.O_LEGEND_CRD]))
        mgr.doPlotOpts()
        mgr.doFigOpts()

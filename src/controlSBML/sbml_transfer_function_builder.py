"""
Builds transfer functions for an SBML model from all inputs to all outputs.
"""

import controlSBML as ctl
import controlSBML.constants as cn
from controlSBML import util
import controlSBML.siso_transfer_function_builder as tfb
import controlSBML.simulate_system as ss
import controlSBML.timeseries as ts
from controlSBML.option_management.option_manager import OptionManager
from controlSBML.option_management.options import Options

from docstring_expander.expander import Expander
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd


class SBMLTransferFunctionBuilder(object):
    
    def __init__(self, ctlsb):
        """
        Construction of an array of control.TransferFunction

        Parameters
        ----------
        ctlsb: ControlSBML
        """
        self.ctlsb = ctlsb
        self.sys = ctlsb.makeNonlinearIOSystem("SBMLTransferFunctionBuilder")

    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotStaircaseResponse(self, staircase_spec=cn.StaircaseSpec(),
           ax2=None, is_steady_state=True, **options):
        """
        Plots the Nonlinear simulation response to a monotonic sequence of step inputs.

        Parameters
        ----------
        staircase_spec: StaircaseSpec
        ax2: Matplotlib.Axes (second y axis)
        is_steady_state: bool (initialize to steady state values)
        #@expand

        Returns
        -------
        dict
            key: str, str (input name, output name)
            value: PlotResult

        """
        mgr = OptionManager(options)
        nrows = self.sys.num_input
        ncols = self.sys.num_output
        fig, axes = plt.subplots(nrows, ncols, figsize=mgr.fig_opts[cn.O_FIGSIZE])
        fig = plt.figure(constrained_layout=True)
        grid_spec = gridspec.GridSpec(ncols=ncols, nrows=nrows, figure=fig)
        mgr.setFigure(fig)
        irow = 0
        icol = 0
        result_dct = {(i, o): [] for i, o in zip(self.ctlsb.input_names, self.sys.ctlsb.output_names)}
        for output_name in self.sys.output_names:
            for input_name in self.sys.input_names:
                if input_name == output_name:
                    continue
                sys = self.sys.getSubsystem(self.sys.name, [input_name], [output_name])
                irow = self.sys.input_names.index(input_name)
                icol = self.sys.output_names.index(output_name)
                ax = fig.add_subplot(grid_spec[irow, icol])
                siso_tfb = tfb.SISOTransferFunctionBuilder(sys)
                mgr.plot_opts[cn.O_AX] = ax
                mgr.plot_opts[cn.O_IS_PLOT] = False
                mgr.plot_opts[cn.O_TITLE] = "%s->%s" % (input_name, output_name)
                plot_result = siso_tfb.plotStaircaseResponse(staircase_spec=staircase_spec, option_mgr=mgr,
                        is_steady_state=is_steady_state)
                result_dct[(input_name, output_name)] = plot_result
                if icol < ncols - 1:
                    plot_result.ax2.set_ylabel("")
                    plot_result.ax2.set_yticklabels([])
                if irow < nrows - 1:
                    plot_result.ax.set_xlabel("")
                    plot_result.ax2.set_xticklabels([])
                if icol > 0:
                    plot_result.ax.set_yticklabels([])
                icol += 1
                if icol >= ncols:
                    icol = 0
                    irow += 1
                mgr.doPlotOpts()
                plot_result.ax.set_title("")
                plot_result.ax.set_title(mgr.plot_opts[cn.O_TITLE], y=0.2, pad=-24, fontsize=10)
        mgr.doFigOpts()
        return result_dct

    @Expander(cn.KWARGS, cn.ALL_KWARGS)    
    def fitTransferFunction(self, num_numerator, num_denominator, is_tf_only=True, **kwargs):
        """
        Constructs transfer functions for the NonlinearIOSystem.

        Parameters
        ----------
        num_numerator: int (number of numerator terms)
        num_denominator: int (number of denominator terms)
        is_tf_only: bool (Values are control.TransferFunction)
        kwargs: dict
            num_step: int (number of steps in staircase)
            initial_value: float (value for first step)
            final_value: float (value for final step)
        #@expand

        Returns
        -------
        DataFrame: 
            column names: str (output)
            index: str (input)
            values: tfb.FitterResult or control.TransferFunction
        """
        result_dct = {n: [] for n in self.sys.output_names}
        for output_name in self.sys.output_names:
            for input_name in self.sys.input_names:
                if input_name == output_name:
                    continue
                sys = self.sys.getSubsystem(self.sys.name, [input_name], [output_name])
                siso_tfb = tfb.SISOTransferFunctionBuilder(sys)
                value = siso_tfb.fitTransferFunction(num_numerator, num_denominator, **kwargs)
                if is_tf_only:
                    value = value.transfer_function
                result_dct[output_name].append(value)
        # Construct the output
        df = pd.DataFrame(result_dct)
        df.columns.name = "Outputs"
        df.index = list(self.sys.input_names)
        df.index.name = "Inputs"
        return df
    
    @classmethod
    def makeTransferFunctionBuilder(cls, *pargs, **kwargs):
        """
        Constructs transfer functions for SBML systems.

        Parameters
        ----------
        Same as for constructing ControlSBML
        """
        ctlsb = ctl.ControlSBML(*pargs, **kwargs)
        return cls(ctlsb)
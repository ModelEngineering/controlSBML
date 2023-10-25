"""
Builds transfer functions for an SBML model from all inputs to all outputs.

1. Staircase specification can be relative to steadystate and/or specific to an input.
"""

import controlSBML.constants as cn
import controlSBML.siso_transfer_function_builder as tfb
from controlSBML.staircase import Staircase
from controlSBML.option_management.option_manager import OptionManager
from controlSBML import msgs

from docstring_expander.expander import Expander
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd


class MIMOTransferFunctionBuilder(object):
    
    def __init__(self, ctlsb, **kwargs):
        """
        Construction of an array of control.TransferFunction

        Parameters
        ----------
        ctlsb: ControlSBML
        kwargs: dict (options for ctl.NonlinearIOSystem)
        """
        self.ctlsb = ctlsb
        self.sys = ctlsb.makeNonlinearIOSystem("SBMLTransferFunctionBuilder", **kwargs)

    @classmethod 
    def _makeResultDF(cls, result_dct, input_names):
        df = pd.DataFrame(result_dct)
        df.columns.name = "Outputs"
        df.index = list(input_names)
        df.index.name = "Inputs"
        return df
    
    def _makeStaircaseDct(self, staircase):
        if isinstance(staircase, dict):
            names = set(staircase.keys())
            diff = names.difference(self.sys.input_names)
            if len(diff) > 0:
                msg = f"Input names {diff} not in know system inputs: {self.sys.input_names}"
                msgs.error(msg)
            staircase_dct = staircase
        else:
            staircase_dct = {n: staircase for n in self.sys.input_names}
        return staircase_dct

    @Expander(cn.KWARGS, cn.SIM_KWARGS)
    def makeStaircaseResponse(self, staircase=Staircase(), is_steady_state=True, **kwargs):
        """
        Plots the Nonlinear simulation response to a monotonic sequence of step inputs.

        Parameters
        ----------
        staircase: Staircase or dict (key: str (input name), value: Staircase)
        is_steady_state: bool (initialize to steady state values)
        #@expand

        Returns
        -------
        Dataframe
            column names: str (output)
            index: str (input)
            values: Timeseries

        """
        mgr = OptionManager(kwargs)
        result_dct = {n: [] for n in self.sys.output_names}
        staircase_dct = self._makeStaircaseDct(staircase)
        for output_name in self.sys.output_names:
            for input_name in self.sys.input_names:
                # Adjust the staircase as required
                sys = self.sys.getSubsystem(self.sys.name, [input_name], [output_name])
                builder = tfb.SISOTransferFunctionBuilder(sys)
                response_df, _ = builder.makeStaircaseResponse(staircase=staircase_dct[input_name],
                                                            is_steady_state=is_steady_state,
                                                             **kwargs)
                result_dct[output_name].append(response_df)
        result_df = self._makeResultDF(result_dct, self.sys.input_names)
        return result_df
    
    @classmethod
    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def _plotMIMO(cls, dataframe, plotFunc, **options):
        """
        Plots the Nonlinear simulation response to a monotonic sequence of step inputs.

        Parameters
        ----------
        dataframe: Dataframe (rows are inputs; columns are outputs)
        plotFunc: function with the signature
            ARGUMENTS
                object in dataframe
                OptionManger (or None)
                Plot keyword arguments
            RETURNS
                PlotResult
        #@expand

        Returns
        -------
        Dataframe
            column names: str (output)
            index: str (input)
            values: PlotResult
        """
        mgr = OptionManager(options)
        input_names = list(dataframe.index)
        output_names = list(dataframe.columns)
        nrow = len(input_names)
        ncol = len(output_names)
        fig = plt.figure(constrained_layout=True)
        grid_spec = gridspec.GridSpec(ncols=ncol, nrows=nrow, figure=fig)
        mgr.setFigure(fig)
        irow = 0
        icol = 0
        result_dct = {n: [] for n in output_names} 
        for input_name in input_names:
            done_first_plot = False
            for output_name in output_names:
                entry = dataframe.loc[input_name, output_name]
                if entry is None:
                    plot_result = None
                else:
                    irow = input_names.index(input_name)
                    icol = output_names.index(output_name)
                    ax = fig.add_subplot(grid_spec[irow, icol])
                    ax2 = ax.twinx()
                    # Plot options
                    mgr.plot_opts[cn.O_AX] = ax
                    mgr.plot_opts[cn.O_AX2] = ax2
                    mgr.plot_opts[cn.O_XLABEL] = "time"
                    plot_result = plotFunc(entry, mgr=mgr)
                    if plot_result is not None:
                        if done_first_plot:
                            ax.set_ylabel("")
                            ax2.set_ylabel("")
                        else:
                            ax2.set_ylabel(input_name, color=cn.INPUT_COLOR)
                            done_first_plot = True
                        mgr.doPlotOpts()
                        if irow < nrow:
                            ax.set_xlabel("")
                            ax2.set_xlabel("")
                        else:
                            ax.set_xlabel("time")
                        plot_result.ax.legend([output_name])
                        ax.get_legend().set_visible(False)
                result_dct[output_name].append(plot_result)
        mgr.doFigOpts()
        result_df = cls._makeResultDF(result_dct, input_names)
        return result_df
    
    @classmethod
    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotStaircaseResponse(cls, response_df, **options):
        """
        Plots the Nonlinear simulation response to a monotonic sequence of step inputs.

        Parameters
        ----------
        response_df: Dataframe (see makeStaircaseResponse)
        #@expand

        Returns
        -------
        Dataframe
            column names: str (output)
            index: str (input)
            values: PlotResult
        """
        return cls._plotMIMO(response_df, tfb.SISOTransferFunctionBuilder.plotStaircaseResponse, **options)
    
    @Expander(cn.KWARGS, cn.SIM_KWARGS)
    def fitTransferFunction(self, num_numerator=cn.DEFAULT_NUM_NUMERATOR, 
                            num_denominator=cn.DEFAULT_NUM_DENOMINATOR, staircase=Staircase(), **sim_kwargs):
        """
        Constructs transfer functions for the NonlinearIOSystem.

        Parameters
        ----------
        num_numerator: int (number of numerator terms)
        num_denominator: int (number of denominator terms)
        staircase: Staircase or dict (key: input name, value: Staircase)
        #@expand

        Returns
        -------
        DataFrame: 
            column names: str (output)
            index: str (input)
            values: cn.FitterResult
        """
        fitter_dct = {n: [] for n in self.sys.output_names}
        staircase_dct = self._makeStaircaseDct(staircase)
        for output_name in self.sys.output_names:
            for input_name in self.sys.input_names:
                sys = self.sys.getSubsystem(self.sys.name, [input_name], [output_name])
                siso_tfb = tfb.SISOTransferFunctionBuilder(sys)
                # FIXME: Always getting the same transfer function
                fitter_result = siso_tfb.fitTransferFunction(num_numerator, num_denominator,
                                                             staircase=staircase_dct[input_name], **sim_kwargs)
                fitter_dct[output_name].append(fitter_result)
        # Construct the output
        fitter_df = self._makeResultDF(fitter_dct, self.sys.input_names)
        return fitter_df
    
    @classmethod
    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotFitTransferFunction(cls, fitter_result_df, **kwargs):
        """
        Plots the results of fitting a transfer function.

        Parameters
        ----------
        fitter_df: Dataframe (see fitTransferFunction)
        #@expand

        Returns
        -------
        Dataframe
            column names: str (output)
            index: str (input)
            values: PlotResult
        """
        return cls._plotMIMO(fitter_result_df, tfb.SISOTransferFunctionBuilder.plotFitterResult, **kwargs)
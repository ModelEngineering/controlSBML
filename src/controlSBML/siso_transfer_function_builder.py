"""
Builds a transfer function for a SISO System with options for the fitting method.

    plotStaircaseResponse: plots response to a staircase input to the transfer function

    TO DO
    1. Tests for fitting
"""

import controlSBML.constants as cn
import controlSBML as ctl
from controlSBML import msgs
from controlSBML import util
from controlSBML.option_management.option_manager import OptionManager
from controlSBML.staircase import Staircase
from controlSBML.poly_fitter import PolyFitter
from controlSBML.gpz_fitter import GPZFitter

from docstring_expander.expander import Expander # type: ignore
import numpy as np


MIN_ELAPSED_TIME = 1e-2
STAIRCASE = "staircase"
MIN_PARAMETER_VALUE = -1e6
MAX_PARAMETER_VALUE = 1e6


class SISOTransferFunctionBuilder(object):

    def __init__(self, sbml_system, input_name=None, output_name=None, fitter_method=cn.FITTER_METHOD_POLY):
        """
        Parameters
        ----------
        sys: SBMLSystem
        input_name: str
        output_name: str
        fitter_method: str (method for fitting the transfer function, either 'poly' or 'gpz')
        """
        #
        self.sbml_system = sbml_system
        self.input_name = input_name
        self.output_name = output_name
        if self.input_name is None:
            self.input_name = sbml_system.input_names[0]
        if self.output_name is None:
            self.output_name = sbml_system.output_names[0]
        self.fitter_method = fitter_method

    def copy(self):
        return SISOTransferFunctionBuilder(self.sbml_system, input_name=self.input_name, output_name=self.output_name,
                                           fitter_method=self.fitter_method)

    @Expander(cn.KWARGS, cn.SIM_KWARGS)
    def makeStaircaseResponse(self, staircase=Staircase(), mgr=None, times=None,
           is_steady_state=True, **kwargs):
        """
        Constructs the staircase response of the NonlinearIOSystem.
        Assumes a single output. If there is more than one,
        only the first is plotted. The operating point of the system is steady state.

        Parameters
        ----------
        staircase: Staircase (num_point will be adjusted as per options)
        is_steady_state: bool (initialize to steady state values)
        mgr: OptionManager
        times: np.array (times at which to simulate)
        #@expand

        Returns
        -------
        ctl.Timeseries
            index: time (ms)
            columns: <output_name>, staircase
        AntimonyBuilder
        """
        # Handle the options. If an option manager is specified, then caller handles figure generation.
        if mgr is None:
            mgr = OptionManager(kwargs)
        # Construct the time
        if times is None:
            start_time = mgr.options.get(cn.O_START_TIME)
            end_time = mgr.options.get(cn.O_END_TIME)
            points_per_time = staircase.num_point/(end_time - start_time)
            times = util.makeSimulationTimes(start_time=start_time, end_time=end_time, points_per_time=points_per_time)
        else:
            points_per_time = len(times)/(times[-1] - times[0])
            mgr.options[cn.START_TIME] = times[0]
            mgr.options[cn.END_TIME] = times[-1]
        mgr.options[cn.O_POINTS_PER_TIME] = points_per_time
        # Do the simulations
        result_ts, antimony_builder = self.sbml_system.simulateStaircase(self.input_name, self.output_name, times=times,
                                                initial_value=staircase.initial_value, num_step=staircase.num_step,
                                                final_value=staircase.final_value, is_steady_state=is_steady_state,
                                                inplace=False)
        staircase.setNumPoint(len(result_ts))
        staircase_name = "%s_%s" % (self.input_name, STAIRCASE)
        staircase_arr= staircase.staircase_arr
        result_ts[staircase_name] = staircase_arr
        del result_ts[self.input_name]
        return result_ts, antimony_builder

    @staticmethod 
    def setYAxColor(ax, position, color):
        # Set the colors of the labels, axes, and spines
        ax.tick_params(axis='y', labelcolor=color)
        ax.spines[position].set_color(color)
        ax.yaxis.label.set_color(color)
    
    @classmethod
    @Expander(cn.KWARGS, cn.PLOT_KWARGS)
    def plotStaircaseResponse(cls, response_ts, mgr=None, **kwargs):
        """
        Plots the response to a monotonic sequence of step inputs. Assumes a
        single input. Assumes a single output. If there is more than one,
        only the first is plotted. The operating point of the system is steady state.
        If option_mgr is not None, then the caller handles figure generation.

        Parameters
        ----------
        response_ts: Timeseries (Staircase response)
            columns: <output_name>, staircase_name
        ax2: Matplotlib.Axes (second y axis)
        is_steady_state: bool (initialize to steady state values)
        #@expand

        Returns
        -------
        util.PlotResult
        """
        # Handle the options. If an option manager is specified, then caller handles figure generation.
        if mgr is None:
            mgr = OptionManager(kwargs)
            is_fig = True
        else:
            is_fig = False
        #
        new_response_ts, staircase_name, output_name = cls._extractStaircaseResponseInformation(response_ts)
        staircase_ts = response_ts[staircase_name]
        response_ts = new_response_ts
        if mgr.plot_opts.asDict().get(cn.IS_PLOT, True):
            # Do the plots
            plot_result = util.plotOneTS(response_ts, mgr=mgr, colors=[cn.SIMULATED_COLOR])
            ax = plot_result.ax
            mgr.plot_opts.set(cn.O_AX, ax)
            mgr.plot_opts.set(cn.O_YLABEL, output_name)
            if mgr.plot_opts[cn.O_AX2] is None:
                ax2 = ax.twinx()
                mgr.plot_opts[cn.O_AX2] = ax2
            else:
                ax2 = mgr.plot_opts[cn.O_AX2]
            # Plot the staircase
            times = np.array(response_ts.index)/cn.MS_IN_SEC
            ax2.plot(times, staircase_ts, color=cn.INPUT_COLOR,
                linestyle="--")
            cls.setYAxColor(ax, "left", cn.SIMULATED_COLOR)
            cls.setYAxColor(ax2, "right", cn.INPUT_COLOR)
            ax2.set_ylabel(staircase_name)
            mgr.doPlotOpts()
            ax.legend([])
            if is_fig:
                mgr.doFigOpts()
        else:
            ax = None
            ax2 = None
        #
        return util.PlotResult(time_series=response_ts, ax=ax, ax2=ax2)
    
    @classmethod
    def _extractStaircaseResponseInformation(cls, timeseries):
        """Extracts informatin present in the staircase response.

        Args:
            timeseries (ctl.Timeseries):
                columns: output_name, staircase_name
                index: time (ms)

        Returns:
            ctl.Timeseries: timeseries without the staircase
            str: name of the staircase column
            str: name of output column
        """
        staircase_column, other_columns = cls._getStaircaseColumnName(timeseries)
        new_timeseries = timeseries[other_columns]
        return new_timeseries, staircase_column, other_columns[0]
    
    @staticmethod
    def _getStaircaseColumnName(timeseries):
        all_columns = set(timeseries.columns)
        other_columns = [c for c in timeseries.columns if STAIRCASE not in c]
        staircase_column = list(all_columns.difference(other_columns))[0]
        return staircase_column, other_columns

    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotTransferFunctionFit(self, num_zero=cn.DEFAULT_NUM_ZERO,
                            num_pole=cn.DEFAULT_NUM_POLE, staircase=Staircase(), 
                            fit_start_time=None, fit_end_time=None, **kwargs):
        """
        Constructs a transfer function for the System. This is done by first estimating the gain to the
        staircase input and then estimating the transients, zeroes and poles.

        Parameters
        ----------
        num_zero: int (number of zeros)
        num_pole: int (number of poles)
        staircase: Staircase
        fit_start_time: float (time at which fitting starts)
        fit_end_time: float (time at which fitting ends)
        kwargs: dict (options as described below)
        #@expand

        Returns
        -------
        FitterResult
            transfer_function: control.TransferFunction
            time_series: ("predicted" <name> <name>_staircase)
            parameters: lmfit.Parameters
            minimizer_result: lmfit.MinimizerResult
            stderr: dict (key: term, value: stderr)
            nfev: number of function evaluations
            redchi: float (reduced ChiSq)
            parameters: lmfit.Parameters
            antimony_builder: AntimonyBuilder
        """
        def adjustArray(arr):
            # Normalizes by mean and shapes
            new_arr = arr - np.mean(arr)
            new_arr = np.reshape(new_arr, (len(new_arr),))
            return new_arr
        #
        fitter_method = kwargs.get(cn.FITTER_METHOD, self.fitter_method)
        new_kwargs = dict(kwargs)
        if cn.FITTER_METHOD in new_kwargs:
            del new_kwargs[cn.FITTER_METHOD]
        mgr = OptionManager(new_kwargs)
        # Get the observational data
        new_staircase = staircase.copy()
        data_ts, antimony_builder = self.makeStaircaseResponse(staircase=new_staircase, **mgr.sim_opts)
        ms_times = util.cleanTimes(data_ts.index)
        _, input_name, _ = self._extractStaircaseResponseInformation(data_ts)
        data_ts.index = ms_times
        #  Construct the fitting data
        start_idx = 0
        end_idx = len(ms_times)
        if fit_start_time is not None:
            start_idx = np.sum(ms_times <= cn.MS_IN_SEC*fit_start_time)
        if fit_end_time is not None:
            end_idx = np.sum(ms_times <= cn.MS_IN_SEC*fit_end_time)
        if np.isclose(start_idx, end_idx):
            msgs.error("fit_start_time >= fit_end_time")
        sel_ms_times = ms_times[start_idx:end_idx]
        new_data_ts = data_ts.loc[sel_ms_times]
        # Do the fit
        if fitter_method == cn.FITTER_METHOD_POLY:
            fitter = PolyFitter(new_data_ts, input_name=input_name, output_name=self.output_name,
                                num_zero=num_zero, num_pole=num_pole)
        elif fitter_method == cn.FITTER_METHOD_GPZ:
            fitter = GPZFitter(new_data_ts, input_name=input_name, output_name=self.output_name,
                                num_zero=num_zero, num_pole=num_pole)
        else:
            raise ValueError("Unknown method: %s" % fitter_method)
        fitter.fit()
        #
        _, y_arr = fitter.simulateTransferFunction(fitter.transfer_function)
        df = new_data_ts.copy()
        df[cn.O_PREDICTED] = y_arr
        output_ts = ctl.Timeseries(df)
        residuals = output_ts[self.output_name] - output_ts[cn.O_PREDICTED].values
        rms_residuals = np.sqrt(np.mean(residuals**2))
        fitter_result = cn.FitterResult(
              input_name=self.input_name,
              output_name=self.output_name,
              transfer_function=fitter.transfer_function,
              time_series=output_ts,
              staircase_name=input_name,
              antimony_builder=antimony_builder,
              rms_residuals=rms_residuals)
        is_plot = kwargs.get(cn.IS_PLOT, True)
        if is_plot:
            new_kwargs = dict(mgr.plot_opts.asDict())
            new_kwargs.update(mgr.fig_opts.asDict())
            fitter.plot(**new_kwargs) 
        return fitter_result
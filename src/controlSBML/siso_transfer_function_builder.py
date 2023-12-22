"""
Builds a transfer function for a SISO NonlinearIOSystem

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

import control
from docstring_expander.expander import Expander
import lmfit
import numpy as np


MIN_ELAPSED_TIME = 1e-2
STAIRCASE = "staircase"
MIN_PARAMETER_VALUE = -10
MAX_PARAMETER_VALUE = 10
INITIAL_PARAMETER_VALUE = 0.1
NUMERATOR_PREFIX = "n"
DENOMINATOR_PREFIX = "d"
MAX_ABS_RESIDUAL = 10

global _mse_history


################## FUNCTIONS ####################
def _makeParameters(num_numerator, num_denominator, min_value=MIN_PARAMETER_VALUE, max_value=MAX_PARAMETER_VALUE,
                    initial_value=INITIAL_PARAMETER_VALUE):
    """
    Makes the parameters used to construct a transfer function.

    Parameters
    ----------
    num_numerator: int (number of numerator terms)
    num_denominator: int (number of denominator terms)
    min_value: float
    max_value: float
    initial_value: float (if None, choose random value)

    Returns
    -------
    lmfit.Parameter
    """
    pfit = lmfit.Parameters()
    for idx in range(num_numerator):
        if initial_value is None:
            initial_value = np.random.uniform(min_value, max_value)
        pfit.add(name='%s%d' % (NUMERATOR_PREFIX, idx),
              min=min_value,
              value=initial_value,
              max=max_value)
    for idx in range(num_denominator):
        if initial_value is None:
            initial_value = np.random.uniform(min_value, max_value)
        pfit.add(name='%s%d' % (DENOMINATOR_PREFIX, idx),
              min=min_value,
              value=initial_value,
              max=max_value)
    return pfit

def _makeTransferFunction(parameters):
    """
    Constructs a transfer function from a dictionary representation.

    Parameters
    ----------
    parameters: lmfit.Parameter
        parameters.valuesdict(): dict
            key=n<int>: numerator coefficient for int-th element
            key=d<int>: denominator coefficient for int-th element

    Returns
    -------
    control.TransferFunction
    """
    tf_dct = parameters.valuesdict()
    def makeVec(letter):
        """
        Creates a vector for the keys are a single letter followed
        by an index.

        Parameters
        ----------
        letter: char

        Returns
        -------
        np.array
        """
        dct = {int(k[1:]): v for k, v in tf_dct.items() if k[0] == letter}
        size = max([i for i in dct.keys()]) + 1
        arr = list(np.repeat(0, size))
        for idx, value in dct.items():
            arr[idx] = value
        return np.array(arr)
    #
    num_arr = makeVec(NUMERATOR_PREFIX)
    den_arr = makeVec(DENOMINATOR_PREFIX)
    tf = control.TransferFunction(num_arr, den_arr)
    new_tf = util.simplifyTransferFunction(tf)
    if len(new_tf.num[0][0]) > len(new_tf.den[0][0]):
        # Avoid improper transfer function
        new_tf = tf
    return new_tf

def _calculateTransferFunctionResiduals(parameters, data_in, data_out):
    """
    Computes the residuals for a transfer function.

    Parameters
    ----------
    parameters: lmfit.Parameters
        n<int>: numerator coefficient for int-th element
        d<int>: denominator coefficient
    data_in: (list-float, list-float) (times, input signal)
    data_out: array-float (calibration data)
    Returns
    -------
    float
    """
    def isDone(last_history:int=10, min_history:int=50, threshold:float=1e-3):
        """
        Checks if further progress is unlikely.

        Returns
        -------
        bool
            True if MSEs are not changing significantly
        """
        if len(_mse_history) < min_history:
            return False
        history_arr = np.array(_mse_history)
        median_mse = np.median(history_arr)
        credible_history_arr = history_arr[history_arr < median_mse]
        last_arr = credible_history_arr[-last_history:]
        metric = np.std(last_arr)/median_mse
        return metric < threshold
    #
    global _mse_history
    #
    times, inputs = data_in
    tf = _makeTransferFunction(parameters)
    times, y_arr = control.forced_response(tf, T=times, U=inputs)
    residuals = data_out - y_arr
    is_bads = [np.isnan(v) or np.isinf(v) or (v is None) for v in residuals]
    is_bad = any(is_bads)
    if is_bad:
        residuals = np.repeat(1e10, len(times))
    mse = np.sum(residuals**2)/len(residuals)
    _mse_history.append(mse)
    if isDone():
        # Force the minimizer to stop
        residuals = np.repeat(0, len(residuals))
    return residuals


################## CLASSES ####################
class SISOTransferFunctionBuilder(object):

    def __init__(self, sbml_system, input_name=None, output_name=None):
        """
        Parameters
        ----------
        sys: SBMLSystem
        input_name: str
        output_name: str
        """
        #
        self.sbml_system = sbml_system
        self.input_name = input_name
        self.output_name = output_name
        if self.input_name is None:
            self.input_name = sbml_system.input_names[0]
        if self.output_name is None:
            self.output_name = sbml_system.output_names[0]
        #

    def copy(self):
        return SISOTransferFunctionBuilder(self.sbml_system, input_name=self.input_name, output_name=self.output_name)

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

    @Expander(cn.KWARGS, cn.SIM_KWARGS)
    def fitTransferFunction(self, num_numerator=cn.DEFAULT_NUM_NUMERATOR,
                            num_denominator=cn.DEFAULT_NUM_DENOMINATOR, staircase=Staircase(), 
                            fit_start_time=None, fit_end_time=None, **kwargs):
        """
        Constructs a transfer function for the NonlinearIOSystem.

        Parameters
        ----------
        num_numerator: int (number of numerator terms)
        num_denominator: int (number of denominator terms)
        staircase: Staircase
        fit_start_time: float (time at which fitting starts)
        fit_end_time: float (time at which fitting ends)
        kwargs: dict (simulation options as described below)
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
        global _mse_history
        #
        # Get the calibration data
        new_staircase = staircase.copy()
        data_ts, antimony_builder = self.makeStaircaseResponse(staircase=new_staircase, **kwargs)
        ms_times = util.cleanTimes(data_ts.index)
        output_ts, staircase_column_name, _ = self._extractStaircaseResponseInformation(data_ts)
        data_ts.index = ms_times
        output_ts.index = ms_times
        new_staircase.name = staircase_column_name
        #  Construct the fitting data
        start_idx = 0
        end_idx = len(ms_times)
        if fit_start_time is not None:
            start_idx = np.sum(ms_times <= cn.MS_IN_SEC*fit_start_time)
        if fit_end_time is not None:
            end_idx = np.sum(ms_times <= cn.MS_IN_SEC*fit_end_time)
        if np.isclose(start_idx, end_idx):
            msgs.error("Start time is greater than end time")
        sel_ms_times = ms_times[start_idx:end_idx]
        sel_ts = data_ts.loc[sel_ms_times]
        staircase_arr = sel_ts[staircase_column_name].values
        sel_sec_times = sel_ms_times/cn.MS_IN_SEC
        data_in = (sel_sec_times, staircase_arr)
        data_out = sel_ts[self.output_name]
        # Do the fit
        parameters = _makeParameters(num_numerator, num_denominator)
        _mse_history = []   # History of mean squared errors
        mini = lmfit.Minimizer(_calculateTransferFunctionResiduals,
                               parameters, fcn_args=(data_in, data_out))
        minimizer_result = mini.leastsq()
        residuals = _calculateTransferFunctionResiduals(minimizer_result.params, data_in, data_out)
        max_abs_residual = np.max(np.abs(residuals))
        if max_abs_residual > MAX_ABS_RESIDUAL:
            msgs.warn("Possible numerical instability: max abs residual is %f" % max_abs_residual)
        rms_residuals = np.sqrt(np.mean((residuals**2)))
        stderr_dct = {k: v.stderr for k,v in minimizer_result.params.items()}
        transfer_function = _makeTransferFunction(minimizer_result.params)
        #
        _, y_arr = control.forced_response(transfer_function, T=sel_sec_times, U=staircase_arr, X0=0)
        output_ts = output_ts.loc[sel_ms_times]
        output_ts[cn.O_PREDICTED] = y_arr
        output_ts[staircase_column_name] = staircase_arr
        output_ts = ctl.Timeseries(output_ts)
        output_ts.index = sel_ms_times
        fitter_result = cn.FitterResult(
              input_name=self.input_name,
              output_name=self.output_name,
              transfer_function=transfer_function,
              stderr=stderr_dct,
              nfev=minimizer_result.nfev,
              rms_residuals = rms_residuals,
              redchi=minimizer_result.redchi,
              time_series=output_ts,
              staircase_arr=staircase_arr,
              staircase_name=staircase_column_name,
              parameters=minimizer_result.params,
              antimony_builder=antimony_builder,
              )
        return fitter_result
    
    @classmethod
    @Expander(cn.KWARGS, cn.PLOT_KWARGS)
    def plotFitterResult(cls, fitter_result, mgr=None, **kwargs):
        """
        Plots the results of fitting a transfer function for the NonlinearIOSystem.
        If an option manager is specified, then the caller handles figure generation.

        Parameters
        ----------
        fitter_result: FitterResult
        mgr: OptionManager
        kwargs: dict
        #@expand
        """
        # Initializations
        if mgr is None:
            is_fig = True
            mgr = OptionManager(kwargs)
        else:
            is_fig = False
        output_name = fitter_result.output_name
        staircase_arr = fitter_result.staircase_arr
        staircase_name = fitter_result.staircase_name
        transfer_function = fitter_result.transfer_function
        times = fitter_result.time_series.index
        #
        ts = fitter_result.time_series.drop(staircase_name, axis=1)
        util.plotOneTS(ts, mgr=mgr,
                       colors=[cn.SIMULATED_COLOR, cn.PREDICTED_COLOR],
                       markers=["o", ""])
        ax = mgr.plot_opts[cn.O_AX]
        mgr.plot_opts.set(cn.O_YLABEL, output_name)
        if mgr.plot_opts[cn.O_AX2] is None:
            ax2 = ax.twinx()
            mgr.plot_opts[cn.O_AX2] = ax2
        else:
            ax2 = mgr.plot_opts[cn.O_AX2]
        ax2.set_ylabel(staircase_name, color=cn.INPUT_COLOR)
        ax2.plot(times/cn.MS_IN_SEC, staircase_arr, color=cn.INPUT_COLOR, linestyle="--")
        latex = util.latexifyTransferFunction(transfer_function)
        if len(mgr.plot_opts[cn.O_TITLE]) == 0:
            #title = "%s->%s;  %s   " % (input_name, output_name, latex)
            title = latex
        else:
            title = mgr.plot_opts[cn.O_TITLE]
        cls.setYAxColor(ax, "left", cn.SIMULATED_COLOR)
        cls.setYAxColor(ax2, "right", cn.INPUT_COLOR)
        ax.set_title(title, y=0.2, pad=-14, fontsize=14, loc="right")
        mgr.doPlotOpts()
        ax.legend([output_name, cn.O_PREDICTED], loc="upper left")
        if is_fig:
            mgr.doFigOpts()
        return util.PlotResult(time_series=fitter_result.time_series, ax=ax, ax2=ax2)
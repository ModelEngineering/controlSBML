"""
Builds a transfer function for a SISO System

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

import control # type: ignore
from docstring_expander.expander import Expander # type: ignore
import lmfit # type: ignore
import numpy as np
from typing import Optional


MIN_ELAPSED_TIME = 1e-2
STAIRCASE = "staircase"
MIN_PARAMETER_VALUE = -1e6
MAX_PARAMETER_VALUE = 1e6
MIN_INITIAL_S_VALUE = -1
MAX_INITIAL_S_VALUE = 0
INITIAL_PARAMETER_VALUE = 0.1
NUMERATOR_PREFIX = "n"
DENOMINATOR_PREFIX = "d"
ZERO_PREFIX = "z"
POLE_PREFIX = "p"
MAX_ABS_RESIDUAL = 10

global _mse_history


################## FUNCTIONS ####################
def makeParameters(num_zero:int, num_pole:int, gain:float, min_value:float=MIN_PARAMETER_VALUE,
                    max_value:float=MAX_PARAMETER_VALUE,
                    initial_value:Optional[float]=None):
    """
    Makes the parameters used to use lmfit to search for a best fitting transfer function.

    Parameters
    ----------
    num_zero: int (number of zeroes in the transfer function)
    num_denominator: int (number of poles in the transfer function)
    gain: float (gain of the transfer function)
    min_value: float
    max_value: float
    initial_value: float (if None, choose random value)

    Returns
    -------
    lmfit.Parameter for zeros (begins with 'z') and poles (begins with 'p')
    """
    def getValue():
        if initial_value is None:
             value = np.random.uniform(MIN_INITIAL_S_VALUE, MAX_INITIAL_S_VALUE)
             return value
        else:
            return initial_value
    #   
    pfit = lmfit.Parameters()
    # Gain
    pfit.add(name="gain", min=min_value, value=gain, max=max_value)
    # Zeros
    for idx in range(num_zero):
        value = getValue()
        pfit.add(name='%s%d' % (ZERO_PREFIX, idx+1),
              min=min_value,
              value=value,
              max=max_value)
    # Poles
    for idx in range(num_pole):
        value = getValue()
        pfit.add(name='%s%d' % (POLE_PREFIX, idx+1),
              min=min_value,
              value=value,
              max=max_value)
    return pfit

def makeTransferFunction(parameters:lmfit.Parameter):
    """
    Constructs a transfer function from a dictionary representation.

    Parameters
    ----------
    parameters: lmfit.Parameter
        parameters.valuesdict(): dict
            key=n<int>: numerator coefficient for int-th element
            key=d<int>: denominator coefficient for int-th element
            "gain": float (gain of the transfer function)

    Returns
    -------
    control.TransferFunction
    """
    s = control.TransferFunction.s
    tf = control.TransferFunction([1], [1])
    for key, value in parameters.valuesdict().items():
        if key[0] == ZERO_PREFIX:
            tf *= (s - value)
        elif key[0] == POLE_PREFIX:
            tf *= 1/(s - value)
        elif key == "gain":
            continue
        else:
            import pdb; pdb.set_trace()
            raise ValueError("Unknown key in transfer function parameters: %s" % key)
    cur_gain = tf.dcgain()
    tf = parameters.valuesdict()["gain"]*tf/cur_gain
    return tf


def simulateTransferFunction(transfer_function:control.TransferFunction, 
                                times:np.ndarray, initial_value:float, inputs:np.ndarray)->np.ndarray:
    """
    Simulates the transfer function.

    Parameters
    ----------
    transfer_function: control.TransferFunction
    times: np.ndarray
    initial_value: float
    inputs: np.ndarray

    Returns
    -------
    np.ndarray
    """
    result_initial = control.initial_response(transfer_function, T=times, X0=initial_value)
    result_input = control.forced_response(transfer_function, T=times, U=inputs)
    y_arr_input = np.reshape(result_input.y, (len(times),)) 
    y_arr_initial = np.reshape(result_initial.y, (len(times),)) 
    y_arr = y_arr_input + y_arr_initial
    return y_arr

def _calculateTransferFunctionResiduals(parameters, data_in, data_out):
    """
    Computes the residuals for a transfer function.

    Parameters
    ----------
    parameters: lmfit.Parameters
        z<int>: zero
        p<int>: pole
        "gain" float (gain of the transfer function)
    data_in: (list-float, list-float) (times, input signal)
    data_out: array-float (calibration data)
    Returns
    -------
    float
    """
    global _mse_history
    _mse_history = []
    def isDone(last_history:int=10, min_history:int=50, threshold:float=1e-3):
        """
        Checks if further progress is unlikely.

        Returns
        -------
        bool
            True if MSEs are not changing significantly
        """
        global _mse_history
        if len(_mse_history) < min_history:  # type: ignore
            return False
        history_arr = np.array(_mse_history)  # type: ignore
        median_mse = np.median(history_arr)
        credible_history_arr = history_arr[history_arr < median_mse]
        last_arr = credible_history_arr[-last_history:]
        metric = np.std(last_arr)/median_mse
        return metric < threshold
    #
    times, inputs = data_in
    tf = makeTransferFunction(parameters)
    y_arr = simulateTransferFunction(tf, times, data_out[0], inputs)
    residuals = data_out - y_arr
    is_bads = [np.isnan(v) or np.isinf(v) or (v is None) for v in residuals]
    if False:
        is_neg_bads = [(v < 0) and b for b, v in zip(is_bads, residuals)]
        is_pos_bads = [(v > 0) and b for b, v in zip(is_bads, residuals)]
        residuals[is_pos_bads] = 1e6
        residuals[is_neg_bads] = -1e6
        other_bads = [b1 and (not b2) and (not b3) for b1, b2, b3 in zip(is_bads, is_neg_bads, is_pos_bads)]
        # For others, randomized negative and positive residuals
        randoms = np.random.randint(0, 2, np.sum(other_bads))
        randoms = (randoms - 1) + randoms
        residuals[other_bads] = randoms*1e6
    if any(is_bads):
        # randomized negative and positive residuals
        randoms = np.random.randint(0, 2, np.sum(is_bads))
        randoms = (randoms - 1) + randoms
        residuals[is_bads] = randoms*1e6
    mse = np.sum(residuals**2)/len(residuals)
    print(tf.poles(), tf.zeros(), tf.dcgain(), mse)
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

    @Expander(cn.KWARGS, cn.SIM_KWARGS)
    def fitTransferFunction(self, num_zero=cn.DEFAULT_NUM_ZERO,
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
        def adjustArray(arr):
            # Normalizes by mean and shapes
            new_arr = arr - np.mean(arr)
            new_arr = np.reshape(new_arr, (len(new_arr),))
            return new_arr
        #
        global _mse_history
        #
        # Get the observational data
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
        # Estimate the gain
        value_arr, idx_arr = new_staircase.makeEndStepInfo(start_idx=start_idx, end_idx=end_idx)
        full_data_arr = data_ts.values[:, 0]
        sel_idx_arr = np.array([n for n in idx_arr if (start_idx <= n) and (n <= end_idx)])
        output_arr = full_data_arr[sel_idx_arr]  # Get the ends of the steps
        adj_value_arr = adjustArray(value_arr)
        adj_output_arr = adjustArray(output_arr)
        gain = adj_output_arr.dot(adj_value_arr)/adj_value_arr.dot(adj_value_arr)
        # Do the fit
        parameters = makeParameters(num_zero, num_pole, gain)
        _mse_history = []   # History of mean squared errors
        out_arr = data_out.values
        minimizer_result = lmfit.minimize(_calculateTransferFunctionResiduals, parameters, args=(data_in, out_arr),
                                          method="differential_evolution")
        residuals = _calculateTransferFunctionResiduals(minimizer_result.params, data_in, out_arr)
        max_abs_residual = np.max(np.abs(residuals))
        if max_abs_residual > MAX_ABS_RESIDUAL:
            msgs.warn("Possible numerical instability: max abs residual is %f" % max_abs_residual)
        rms_residuals = np.sqrt(np.mean((residuals**2)))
        stderr_dct = {k: v.stderr for k,v in minimizer_result.params.items()}
        transfer_function = makeTransferFunction(minimizer_result.params)
        #
        y_arr = simulateTransferFunction(transfer_function, sel_sec_times, out_arr[0], staircase_arr)
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
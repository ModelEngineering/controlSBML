"""
Builds a transfer function for a SISO NonlinearIOSystem

    plotStaircaseResponse: plots response to a staircase input to the transfer function

    TO DO
    1. Tests for fitting
"""

import controlSBML.constants as cn
from controlSBML import msgs
from controlSBML import util
import controlSBML.simulate_system as ss
import controlSBML.timeseries as ts
from controlSBML.option_management.option_manager import OptionManager
from controlSBML.option_management.options import Options

import collections
import control
from docstring_expander.expander import Expander
import lmfit
import numpy as np
import pandas as pd


MIN_ELAPSED_TIME = 1e-2
STAIRCASE = "staircase"
MIN_PARAMETER_VALUE = -10
MAX_PARAMETER_VALUE = 10
INITIAL_PARAMETER_VALUE = 0.1
NUMERATOR_PREFIX = "n"
DENOMINATOR_PREFIX = "d"


################## FUNCTIONS ####################
def makeParameters(num_numerator, num_denominator):
    """
    Makes the parameters used to construct a transfer function.

    Parameters
    ----------
    num_numerator: int (number of numerator terms)
    num_denominator: int (number of denominator terms)

    Returns
    -------
    lmfit.Parameter
    """
    pfit = lmfit.Parameters()
    for idx in range(num_numerator):
        pfit.add(name='%s%d' % (NUMERATOR_PREFIX, idx),
              min=MIN_PARAMETER_VALUE,
              value=INITIAL_PARAMETER_VALUE,
              max=MAX_PARAMETER_VALUE)
    for idx in range(num_denominator):
        pfit.add(name='%s%d' % (DENOMINATOR_PREFIX, idx),
              min=MIN_PARAMETER_VALUE,
              value=INITIAL_PARAMETER_VALUE,
              max=MAX_PARAMETER_VALUE)
    return pfit

def makeTransferFunction(parameters):
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
    return control.TransferFunction(num_arr, den_arr)

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
    times, inputs = data_in
    tf = makeTransferFunction(parameters)
    times, y_arr = control.forced_response(tf, T=times, U=inputs)
    residuals = data_out - y_arr
    is_bads = [np.isnan(v) or np.isinf(v) or (v is None) for v in residuals]
    is_bad = any(is_bads)
    if is_bad:
        residuals = np.repeat(1e10, len(times))
    return residuals


################## CLASSES ####################
class SISOTransferFunctionBuilder(object):

    def __init__(self, sys, input_name=None, output_name=None):
        """
        Parameters
        ----------
        sys: NonlinearIOSystem with single input and single output
        input_name: str
        output_name: str
        """
        if len(sys.input_labels) == 0:
            raise ValueError("System must have an input.")
        if len(sys.output_labels) == 0:
            raise ValueError("System must have an output.")
        #
        self.sys = sys
        self.input_name = input_name
        self.output_name = output_name
        if self.input_name is None:
            self.input_name = sys.input_labels[0]
        if self.output_name is None:
            self.output_name = sys.output_labels[0]

    @staticmethod
    def _makeStaircase(num_point, num_step, initial_value, final_value):
        """
        A staircase is a sequence of steps of the same magnitude and duration.

        Parameters
        ----------
        num_point: number of points in the staircase
        num_step: int (number of steps in the stair response.
        start_value: float (initial values of the inputs)
        final_value: float (ending values of the inputs)

        Returns
        -------
        np.ndarray
        """
        steps = []  # Steps in the staircase
        num_point_in_step = int(num_point/num_step)
        for num in range(num_step):
            steps.extend(list(np.repeat(num, num_point_in_step)))
        num_added_point = num_point - len(steps)
        if num_added_point < 0:
            raise RuntimeError("Negative residual count")
        elif num_added_point > 0:
            steps.extend(list(np.repeat(num_step - 1, num_added_point)))
        staircase_arr = np.array(steps)
        # Rescale
        staircase_arr = staircase_arr*(final_value - initial_value)/(num_step - 1)
        staircase_arr += initial_value
        #
        return staircase_arr

    @staticmethod
    def getStaircaseArr(time_series):
        """
        Extracts the staircase array.

        Parameters
        ----------
        time_series: Timeseries

        Returns
        -------
        np.array
        """
        columns = list(time_series.columns)
        sel_columns = [c for c in columns if STAIRCASE in c]
        if len(sel_columns) != 1:
            raise ValueError("Invalid staircase timeseries")
        column = sel_columns[0]
        return time_series[column].values


    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotStaircaseResponse(self, final_value=0, num_step=5, initial_value=0,
           ax2=None, **kwargs):
        """
        Plots the response to a monotonic sequence of step inputs. Assumes a
        single input. Assumes a single output. If there is more than one,
        only the first is plotted.

        Parameters
        ----------
        num_step: int (number of steps in staircase)
        initial_value: float (value for first step)
        final_value: float (value for final step)
        ax2: Matplotlib.Axes (second y axis)
        #@expand

        Returns
        -------
        util.PlotResult
        """
        # Handle the options
        mgr = OptionManager(kwargs)
        start_time = mgr.options.get(cn.O_START_TIME)
        end_time = mgr.options.get(cn.O_END_TIME)
        points_per_time = mgr.options.get(cn.O_POINTS_PER_TIME)
        is_plot = mgr.options.get(cn.O_IS_PLOT)
        # Construct the staircase inputs
        num_point = points_per_time*(end_time - start_time) + 1
        staircase_arr = self._makeStaircase(num_point, num_step, initial_value,
               final_value)
        # Restructure
        result_ts = ss.simulateSystem(self.sys, u_vec=staircase_arr,
               start_time=start_time, output_names=[self.output_name],
               end_time=end_time, points_per_time=points_per_time)
        staircase_name = "%s_%s" % (self.input_name, STAIRCASE)
        result_ts[staircase_name] = staircase_arr
        # Do the plots
        ax = None
        ax2 = None
        if is_plot:
            plot_opts = Options(mgr.plot_opts, cn.DEFAULT_DCTS)
            # Plot the output
            output_ts = result_ts.copy()
            del output_ts[staircase_name]
            column_names = list(result_ts)
            revised_opts = Options(plot_opts, cn.DEFAULT_DCTS)
            revised_opts.set(cn.O_WRITEFIG, False)
            revised_opts.set(cn.O_IS_PLOT,  False)
            plot_result = util.plotOneTS(output_ts, **revised_opts)
            ax = plot_result.ax
            if ax2 is None:
                ax2 = ax.twinx()
            # Plot the staircase
            times = np.array(result_ts.index)/cn.MS_IN_SEC
            ax2.plot(times, result_ts[staircase_name], color="red",
                  linestyle="--")
            ax2.set_ylabel(staircase_name, color="red")
            ax2.legend([])
            mgr.doFigOpts()
        #
        return util.PlotResult(time_series=result_ts, ax=ax, ax2=ax2)

    @Expander(cn.KWARGS, cn.SIM_KWARGS)
    def fitTransferFunction(self, num_numerator, num_denominator, **kwargs):
        """
        Constructs a transfer function for the NonlinearIOSystem.

        Parameters
        ----------
        num_numerator: int (number of numerator terms)
        num_denominator: int (number of denominator terms)
        kwargs: dict
            num_step: int (number of steps in staircase)
            initial_value: float (value for first step)
            final_value: float (value for final step)
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
        """
        FitterResult = collections.namedtuple("FitterResult",
            "transfer_function parameters minimizer_result stderr nfev redchi time_series")
        # Initialize staircase arguments
        final_value = kwargs.get("final_value", 0)
        initial_value = kwargs.get("initial_value", 0)
        num_step = kwargs.get("num_step", 5)
        new_kwargs = dict(kwargs)
        new_kwargs["final_value"] = final_value
        new_kwargs["initial_value"] = initial_value
        new_kwargs["num_step"] = num_step
        # Get the calibration data
        plot_result = self.plotStaircaseResponse(is_plot=False,
              **new_kwargs)
        data_ts = plot_result.time_series
        times = data_ts.times
        times_diff = np.diff(times)
        if not np.allclose(times_diff, times_diff[0]):
            import pdb; pdb.set_trace()
            pass
        staircase_arr = self.getStaircaseArr(data_ts)
        data_in = (times, staircase_arr)
        data_out = data_ts[self.output_name]
        # Do the fit
        parameters = makeParameters(num_numerator, num_denominator)
        mini = lmfit.Minimizer(_calculateTransferFunctionResiduals,
                               parameters, fcn_args=(data_in, data_out))
        minimizer_result = mini.leastsq()
        stderr_dct = {k: v.stderr for k,v in minimizer_result.params.items()}
        transfer_function = makeTransferFunction(minimizer_result.params)
        #
        time_series = data_ts.copy()
        _, y_arr = control.forced_response(transfer_function,
              T=times, U=staircase_arr)
        time_series["predicted"] = y_arr
        fitter_result = FitterResult(transfer_function=transfer_function,
              minimizer_result=minimizer_result,
              stderr=stderr_dct,
              nfev=minimizer_result.nfev,
              redchi=minimizer_result.redchi,
              time_series=time_series,
              parameters=minimizer_result.params)
        return fitter_result


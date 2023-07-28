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
from controlSBML.option_management.option_manager import OptionManager
from controlSBML.option_management.options import Options
from controlSBML.staircase import Staircase

import collections
import control
from docstring_expander.expander import Expander
import lmfit
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sympy


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
    tf = control.TransferFunction(num_arr, den_arr)
    return util.simplifyTransferFunction(tf)

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
        sys: NonlinearIOSystem
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
    def plotStaircaseResponse(self, staircase=Staircase(),
           ax2=None, is_steady_state=True, option_mgr=None, **kwargs):
        """
        Plots the response to a monotonic sequence of step inputs. Assumes a
        single input. Assumes a single output. If there is more than one,
        only the first is plotted. The operating point of the system is steady state.

        Parameters
        ----------
        staircase: Staircase (num_point will be adjusted as per options)
        ax2: Matplotlib.Axes (second y axis)
        is_steady_state: bool (initialize to steady state values)
        #@expand

        Returns
        -------
        util.PlotResult
        """
        # Handle the options
        if option_mgr is None:
            option_mgr = OptionManager(kwargs)
            is_fig = True
        else:
            is_fig = False
        #
        start_time = option_mgr.options.get(cn.O_START_TIME)
        end_time = option_mgr.options.get(cn.O_END_TIME)
        points_per_time = option_mgr.options.get(cn.O_POINTS_PER_TIME)
        is_plot = option_mgr.options.get(cn.O_IS_PLOT)
        # Construct the staircase inputs
        staircase.num_point = (end_time-start_time)*points_per_time + 1
        staircase_arr = staircase.staircase_arr
        # Do the simulations
        if is_steady_state:
            success = self.sys.setSteadyState()
            if not success:
                msgs.warn("Could not set steady state.")
        result_ts = ss.simulateSystem(self.sys, u_vec=staircase_arr,
               start_time=start_time, output_names=[self.output_name],
               is_steady_state=True,
               end_time=end_time, points_per_time=points_per_time)
        staircase_name = "%s_%s" % (self.input_name, STAIRCASE)
        result_ts[staircase_name] = staircase_arr
        ax = None
        ax2 = None
        if is_plot:
            # Do the plots
            ax2 = None
            plot_opts = Options(option_mgr.plot_opts, cn.DEFAULT_DCTS)
            # Plot the output
            output_ts = result_ts.copy()
            del output_ts[staircase_name]
            revised_opts = Options(plot_opts, cn.DEFAULT_DCTS)
            revised_opts.set(cn.O_IS_PLOT,  False)
            revised_opts.set(cn.O_FIGURE, option_mgr.fig_opts[cn.O_FIGURE])
            plot_result = util.plotOneTS(output_ts, **revised_opts)
            ax = plot_result.ax
            ax.legend([self.input_name], loc="lower left")
            if ax2 is None:
                ax2 = ax.twinx()
            # Plot the staircase
            times = np.array(result_ts.index)/cn.MS_IN_SEC
            ax2.plot(times, result_ts[staircase_name], color="red",
                linestyle="--")
            ax2.set_ylabel(staircase_name, color="red")
            ax2.legend([staircase_name], loc="lower right")
            option_mgr.doPlotOpts()
            if is_fig:
                option_mgr.doFigOpts()
        #
        return util.PlotResult(time_series=result_ts, ax=ax, ax2=ax2, staircase=staircase)
    
    @staticmethod
    def _extractFromTimeseries(timeseries):
        new_timeseries = timeseries.copy()
        all_columns = set(new_timeseries.columns)
        other_columns = [c for c in new_timeseries.columns if STAIRCASE not in c]
        staircase_column = list(all_columns.difference(other_columns))[0]
        new_timeseries = new_timeseries[other_columns]
        return new_timeseries, staircase_column

    @Expander(cn.KWARGS, cn.SIM_KWARGS)
    def fitTransferFunction(self, num_numerator, num_denominator, staircase=Staircase(), 
                            **kwargs):
        """
        Constructs a transfer function for the NonlinearIOSystem.

        Parameters
        ----------
        num_numerator: int (number of numerator terms)
        num_denominator: int (number of denominator terms)
        staircase: Staircase
        kwargs: dict
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
        # Initialize staircase arguments
        new_kwargs = dict(kwargs)
        # Get the calibration data
        new_staircase = staircase.copy()
        plot_result = self.plotStaircaseResponse(is_plot=False, staircase=new_staircase, **kwargs)
        data_ts = plot_result.time_series
        time_series, staircase_column_name = self._extractFromTimeseries(data_ts)
        new_staircase.name = staircase_column_name
        times = data_ts.times
        times_diff = np.diff(times)
        if not np.allclose(times_diff, times_diff[0]):
            import pdb; pdb.set_trace()
            pass
        staircase_arr = plot_result.staircase.staircase_arr
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
        _, y_arr = control.forced_response(transfer_function,
              T=times, U=staircase_arr)
        time_series[cn.O_PREDICTED] = y_arr
        fitter_result = cn.FitterResult(
              transfer_function=transfer_function,
              stderr=stderr_dct,
              nfev=minimizer_result.nfev,
              redchi=minimizer_result.redchi,
              time_series=time_series,
              staircase=new_staircase,
              parameters=minimizer_result.params)
        return fitter_result
    

    @Expander(cn.KWARGS, cn.PLOT_KWARGS)
    def plotFit(self, fitter_result, is_plot=True, **kwargs):
        """
        Plots the results of fitting a transfer function for the NonlinearIOSystem.

        Parameters
        ----------
        fitter_result: FitterResult
        kwargs: dict
        #@expand
        """
        # Initializations
        mgr = OptionManager(kwargs)
        new_kwargs = util.setNoPlot(kwargs)
        staircase = fitter_result.staircase
        staircase_arr = staircase.staircase_arr
        transfer_function = fitter_result.transfer_function
        times = fitter_result.time_series.times
        #
        util.plotOneTS(fitter_result.time_series, mgr=mgr, **new_kwargs)
        ax = plt.gca()
        ax.legend([self.input_name, cn.O_PREDICTED], loc="lower right")
        ax2 = ax.twinx()
        ax2.set_ylabel(staircase.name, color="red")
        ax2.plot(times, staircase_arr, color="red", linestyle="--")
        numr = transfer_function.num[0][0]
        denr = transfer_function.den[0][0]
        latex = r'$\frac{%s}{%ss + %s}$' % (numr[0], denr[0], denr[1])
        latex = util.latexifyTransferFunction(transfer_function)
        if len(mgr.plot_opts[cn.O_TITLE]) == 0:
            title = "%s->%s;  %s   " % (self.input_name, self.output_name, latex)
        else:
            title = mgr.plot_opts[cn.O_TITLE]
        ax.set_title(title, y=0.2, pad=-24, fontsize=14, loc="right")
        ax.legend([self.output_name, cn.O_PREDICTED], loc="lower right")
        mgr.doPlotOpts()
        if is_plot:
            mgr.doFigOpts()
        return fitter_result


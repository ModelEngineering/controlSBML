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
    
    @Expander(cn.KWARGS, cn.SIM_KWARGS)
    def makeStaircaseResponse(self, staircase=Staircase(), mgr=None,
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
        #@expand

        Returns
        -------
        ctl.Timeseries
            index: time (ms)
            columns: <output_name>, staircase
        """
        # Handle the options. If an option manager is specified, then caller handles figure generation.
        if mgr is None:
            mgr = OptionManager(kwargs)
        #
        start_time = mgr.options.get(cn.O_START_TIME)
        end_time = mgr.options.get(cn.O_END_TIME)
        points_per_time = mgr.options.get(cn.O_POINTS_PER_TIME)
        # Construct the staircase inputs
        staircase.num_point = (end_time-start_time)*points_per_time + 1
        staircase_arr = staircase.staircase_arr
        # Do the simulations
        if is_steady_state:
            success = self.sys.setSteadyState()
            if not success:
                msgs.warn("Could not find a steady state. Using current state.")
        result_ts = ss.simulateSystem(self.sys, u_vec=staircase_arr,
               start_time=start_time, output_names=[self.output_name],
               is_steady_state=True,
               end_time=end_time, points_per_time=points_per_time)
        staircase_name = "%s_%s" % (self.input_name, STAIRCASE)
        result_ts[staircase_name] = staircase_arr
        return result_ts
    
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
        # Set the colors of the labels, axes, and spines
        def setYAxColor(ax, position, color):
            ax.tick_params(axis='y', labelcolor=color)
            ax.spines[position].set_color(color)
            ax.yaxis.label.set_color(color)
        #
        AX_COLOR = "green"
        AX2_COLOR = "red"
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
        plot_result = util.plotOneTS(response_ts, mgr=mgr, color=AX_COLOR)
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
        ax2.plot(times, staircase_ts, color=AX2_COLOR,
            linestyle="--")
        setYAxColor(ax, "left", AX_COLOR)
        setYAxColor(ax2, "right", AX2_COLOR)
        ax2.set_label(staircase_name)
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
    def fitTransferFunction(self, num_numerator, num_denominator, staircase=Staircase(), 
                            **kwargs):
        """
        Constructs a transfer function for the NonlinearIOSystem.

        Parameters
        ----------
        num_numerator: int (number of numerator terms)
        num_denominator: int (number of denominator terms)
        staircase: Staircase
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
        """
        # Get the calibration data
        new_staircase = staircase.copy()
        data_ts = self.makeStaircaseResponse(staircase=new_staircase, **kwargs)
        output_ts, staircase_column_name, _ = self._extractStaircaseResponseInformation(data_ts)
        new_staircase.name = staircase_column_name
        times = data_ts.times
        times_diff = np.diff(times)
        if not np.allclose(times_diff, times_diff[0]):
            import pdb; pdb.set_trace()
            pass
        staircase_arr = new_staircase.staircase_arr
        data_in = (times, staircase_arr)
        data_out = data_ts[self.output_name]
        # Do the fit
        parameters = makeParameters(num_numerator, num_denominator)
        mini = lmfit.Minimizer(_calculateTransferFunctionResiduals,
                               parameters, fcn_args=(data_in, data_out))
        minimizer_result = mini.leastsq()
        residuals = _calculateTransferFunctionResiduals(minimizer_result.params, data_in, data_out)
        rms_residuals = np.sqrt(np.mean(residuals**2))
        stderr_dct = {k: v.stderr for k,v in minimizer_result.params.items()}
        transfer_function = makeTransferFunction(minimizer_result.params)
        #
        _, y_arr = control.forced_response(transfer_function, T=times, U=staircase_arr, X0=0)
        output_ts[cn.O_PREDICTED] = y_arr
        fitter_result = cn.FitterResult(
              input_name=self.input_name,
              output_name=self.output_name,
              transfer_function=transfer_function,
              stderr=stderr_dct,
              nfev=minimizer_result.nfev,
              rms_residuals = rms_residuals,
              redchi=minimizer_result.redchi,
              time_series=output_ts,
              staircase=new_staircase,
              parameters=minimizer_result.params)
        return fitter_result
    
    @classmethod
    @Expander(cn.KWARGS, cn.PLOT_KWARGS)
    def plotFitTransferFunction(cls, fitter_result, mgr=None, **kwargs):
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
        staircase = fitter_result.staircase
        staircase_arr = staircase.staircase_arr
        transfer_function = fitter_result.transfer_function
        times = fitter_result.time_series.times
        #
        util.plotOneTS(fitter_result.time_series, mgr=mgr)
        ax = mgr.plot_opts[cn.O_AX]
        if mgr.plot_opts[cn.O_AX2] is None:
            ax2 = ax.twinx()
            mgr.plot_opts[cn.O_AX2] = ax2
        else:
            ax2 = mgr.plot_opts[cn.O_AX2]
        ax2.set_ylabel(staircase.name, color="red")
        ax2.plot(times, staircase_arr, color="red", linestyle="--")
        latex = util.latexifyTransferFunction(transfer_function)
        if len(mgr.plot_opts[cn.O_TITLE]) == 0:
            #title = "%s->%s;  %s   " % (input_name, output_name, latex)
            title = latex
        else:
            title = mgr.plot_opts[cn.O_TITLE]
        ax.set_title(title, y=0.2, pad=-14, fontsize=14, loc="right")
        mgr.doPlotOpts()
        ax.legend([output_name, cn.O_PREDICTED], loc="upper left")
        if is_fig:
            mgr.doFigOpts()
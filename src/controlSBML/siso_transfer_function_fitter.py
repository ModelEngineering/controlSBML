"""
Efficient fitting of transfer functions to data. The fit proceeds in the following steps:
1. Fit the gain of the transfer function.
2. Fit the poles of the transfer function.
3. Fit the zeros of the transfer function.

Todo:
1. Selective specify minimimum and max values for poles, zeros
"""

import controlSBML.constants as cn
from controlSBML import util
from controlSBML.option_management.option_manager import OptionManager
from controlSBML.timeseries import Timeseries

import control # type: ignore
from docstring_expander.expander import Expander # type: ignore
import itertools
import lmfit # type: ignore
import numpy as np
import scipy  # type: ignore
from typing import Optional, Tuple, Union


MIN_PARAMETER_VALUE = -1e6
MAX_PARAMETER_VALUE = 1e6
MIN_INITIAL_S_VALUE = -1
MAX_INITIAL_S_VALUE = 0
DEFAULT_MIN_POLE_VALUE = -5
DEFAULT_MAX_POLE_VALUE = -0.01
DEFAULT_MIN_ZERO_VALUE = -100
DEFAULT_MAX_ZERO_VALUE = DEFAULT_MAX_POLE_VALUE
ZERO_PREFIX = "z"
POLE_PREFIX = "p"


class SISOTransferFunctionFitter(object):

    def __init__(self, timeseries:Timeseries, num_zero:Optional[int]=0, num_pole:Optional[int]=2,
                 input_name:Optional[str]=None, output_name:Optional[str]=None,
                 min_pole_value:float=DEFAULT_MIN_POLE_VALUE,
                 max_pole_value:float=DEFAULT_MAX_POLE_VALUE,
                 min_zero_value:float=DEFAULT_MIN_ZERO_VALUE,
                 max_zero_value:float=DEFAULT_MAX_ZERO_VALUE,
                 ):
        """
        Parameters
        ----------
        timeseries: Timeseries (input has "staircase" in name; time is index; output is a column)
        num_zero: int (number of zeroes in the transfer function)
        num_denominator: int (number of poles in the transfer function)
        """
        self.input_name, self.output_name = self._extractNames(timeseries, input_name, output_name)
        self.in_arr = timeseries[self.input_name].values
        self.out_arr = timeseries[self.output_name].values
        dt = np.round(scipy.stats.mode(np.diff(timeseries.index))[0])/cn.MS_IN_SEC
        self.times = np.array([dt*n for n in range(len(timeseries.index))])  # Must be evenly spaced
        self.length = len(self.times)
        self.num_zero:int = num_zero  # type: ignore
        self.num_pole:int = num_pole # type: ignore
        self.min_zero_value = self.setValue(min_zero_value, DEFAULT_MIN_ZERO_VALUE)
        self.max_zero_value = self.setValue(max_zero_value, DEFAULT_MAX_ZERO_VALUE)
        self.min_pole_value = self.setValue(min_pole_value, DEFAULT_MIN_POLE_VALUE)
        self.max_pole_value = self.setValue(max_pole_value, DEFAULT_MAX_POLE_VALUE)
        # Outputs of calculations
        self.dcgain:Union[None, float] = None
        self.poles:Union[None, np.ndarray[float]] = None
        self.zeros:Union[None, np.ndarray[float]] = None
        self.transfer_function:Union[None, control.TransferFunction] = None

    @staticmethod
    def setValue(value, default):
        if value is None:
            return default
        else:
            return value

    @staticmethod
    def _extractNames(timeseries:Timeseries, input_name:Union[str, None], output_name:Union[str, None])->Tuple[str, str]:
        """
        Extracts the names of the input and output from the timeseries.

        Parameters
        ----------
        timeseries: Timeseries

        Returns
        -------
        Tuple[str, str]
        """
        columns = timeseries.columns
        if len(columns) != 2:
            raise ValueError("Timeseries must have two columns")
        # Find the input name
        if input_name is None:
            input_names = [n for n in columns if "staircase" in n]
            if len(input_names) != 1:
                raise ValueError("Timeseries must have one column with 'staircase' in name")
            input_name = input_names[0]
        # Find the output name
        if output_name is None:
            output_names = [n for n in columns if n != input_name]
            if len(output_names) != 1:
                raise ValueError("Timeseries must have one column that is not the input")
            output_name = output_names[0]
        else:
            if output_name not in columns:
                raise ValueError("Output name not in Timeseries")
        return input_name, output_name

    def makeParameters(self, zeros:Optional[np.ndarray[complex]]=None,
                       poles:Optional[np.ndarray[complex]]=None)->lmfit.Parameters:
        """
        Makes the parameters used to use lmfit to search for a best fitting transfer function.

        Parameters
        ----------
        zeros: np.ndarray[float] (None - use default or create random)
        poles: np.ndarray[float]

        Returns
        -------
        lmfit.Parameter for zeros (begins with 'z') and poles (begins with 'p')
        """
        def getValues(values, default_values, min_values, max_values):
            if values is not None:
                # Specified a list of values
                results = values
            else:
                # No values specified. Use existing if present
                if default_values is not None:
                    results = default_values
                else:
                    results = []
                    for min_value, max_value in zip(min_values, max_values):
                        results.append(self.uniformFromLogspace(min_value, max_value, 1))
            return results
        #
        zeros = getValues(zeros, self.zeros, self.min_zero_value, self.max_zero_value)
        poles = getValues(poles, self.poles, self.min_pole_value, self.max_pole_value)
        pfit = lmfit.Parameters()
        # Zeros
        for idx, zero in enumerate(zeros):  # type: ignore
            pfit.add(name='%s%d' % (ZERO_PREFIX, idx+1),
                min=self.min_zero_value,
                value=zero,
                max=self.max_zero_value)
        # Poles
        for idx, pole in enumerate(poles): # type: ignore
            pfit.add(name='%s%d' % (POLE_PREFIX, idx+1),
                min=self.min_pole_value,
                value=pole,
                max=self.max_pole_value)
        return pfit

    def makeTransferFunction(self, parameters:lmfit.Parameter):
        """
        Constructs a transfer function from a parameter representation. Gain must have been calculated previously.

        Parameters
        ----------
        parameters: lmfit.Parameter
            name=n<int>: numerator coefficient for int-th element
            name=d<int>: denominator coefficient for int-th element

        Returns
        -------
        control.TransferFunction
        """
        if self.dcgain is None:
            raise ValueError("dcgain must be calculated before calling makeTransferFunction")
        #
        s = control.TransferFunction.s
        tf = control.TransferFunction([1], [1])
        for key, value in parameters.valuesdict().items():
            if key[0] == ZERO_PREFIX:
                tf *= (s - value)
            elif key[0] == POLE_PREFIX:
                tf *= 1/(s - value)
            else:
                raise ValueError("Unknown key in transfer function parameters: %s" % key)
        cur_gain = tf.dcgain()
        tf = self.dcgain*tf/cur_gain
        return tf

    def simulateTransferFunction(self, transfer_function:control.TransferFunction)->Tuple[np.ndarray, np.ndarray]:
        """
        Simulates the transfer function.

        Parameters
        ----------
        transfer_function: control.TransferFunction

        Returns
        -------
        np.ndarray - times
        np.ndarray - output
        """
        result_input = control.forced_response(transfer_function, T=self.times, U=self.in_arr)
        y_arr_input = np.reshape(result_input.y, (self.length,)) 
        return self.times, y_arr_input
    
    def fit(self):
        """
        Estimates the gain of the transfer function.

        Returns
        -------
        float
        """
        self.fitDCGain()
        self.fitPoles()
        self.fitZeros()
        parameters = self.makeParameters(self.poles, self.zeros)
        self.transfer_function = self.makeTransferFunction(parameters)
    
    def fitDCGain(self):
        """
        Estimates the gain of the transfer function.

        Returns
        -------
        float
        """
        adj_out_arr = self.out_arr - np.mean(self.out_arr)
        adj_in_arr = self.in_arr - np.mean(self.in_arr)
        self.dcgain = adj_out_arr.dot(adj_in_arr)/adj_in_arr.dot(adj_in_arr)
    
    def fitPoles(self, num_itr:int=10)->None:
        """
        Estimates the poles of the transfer function.

        Parameters
        ----------
        num_itr: int (number of iterations in the grid search)

        Returns
        -------
        np.ndarray
        """
        self.poles = None  # type: ignore
        best_mse: float = 1e6
        best_poles = np.repeat(None, self.num_pole)  # type: ignore
        for poles in self.getNext(self.min_pole_value, self.max_pole_value, num_itr, self.num_pole):
            parameters = self.makeParameters([], poles)  # type: ignore
            tf = self.makeTransferFunction(parameters)
            mse = self._calculateResiduals(tf)
            if mse < best_mse:
                best_mse = mse
                best_poles = poles
        self.poles = np.array(best_poles)

    # Fix: Handle a mix of positive and negative zeros
    def fitZeros(self, num_itr:int=10)->None:
        """
        Estimates the zeros of the transfer function.

        Returns
        -------
        np.ndarray
        """
        if self.poles is None:
            raise ValueError("Poles must be calculated before calling fitZeros")
        def search(min_value, max_value):
            self.zeros = None  # type: ignore
            best_mse = 1e6
            best_zeros = np.repeat(None, self.num_zero)  # type: ignore
            for zeros in self.getNext(min_value, max_value, num_itr, self.num_zero):
                parameters = self.makeParameters(zeros, self.poles)  # type: ignore
                tf = self.makeTransferFunction(parameters)
                mse = self._calculateResiduals(tf)
                if mse < best_mse:
                    best_mse = mse
                    best_zeros = zeros
            return best_zeros, best_mse
        #
        if (self.min_zero_value<0 and self.max_zero_value>0):
            max_zero1 = -0.01
            best_zeros1, best_mse1 = search(self.min_zero_value, max_zero1)
            min_zero2 = 0.01
            best_zeros2, best_mse2 = search(min_zero2, self.max_zero_value)
            if best_mse1 < best_mse2:
                self.zeros = best_zeros1
            else:
                self.zeros = best_zeros2
        else:
            self.zeros, _ = search(self.min_zero_value, self.max_zero_value)
    
    def plot(self, **kwargs):
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
            mgr = OptionManager(kwargs)
        #
        ts = Timeseries({cn.TIME: self.times, self.output_name: self.out_arr, self.input_name: self.in_arr})
        util.plotOneTS(ts, mgr=mgr,
                       colors=[cn.SIMULATED_COLOR, cn.PREDICTED_COLOR],
                       markers=["o", ""])
        ax = mgr.plot_opts[cn.O_AX]
        mgr.plot_opts.set(cn.O_YLABEL, self.output_name)
        if mgr.plot_opts[cn.O_AX2] is None:
            ax2 = ax.twinx()
            mgr.plot_opts[cn.O_AX2] = ax2
        else:
            ax2 = mgr.plot_opts[cn.O_AX2]
        ax2.set_ylabel(self.input_name, color=cn.INPUT_COLOR)
        ax2.plot(self.times/cn.MS_IN_SEC, self.in_arr, color=cn.INPUT_COLOR, linestyle="--")
        latex = util.latexifyTransferFunction(self.transfer_function)
        if len(mgr.plot_opts[cn.O_TITLE]) == 0:
            #title = "%s->%s;  %s   " % (input_name, output_name, latex)
            title = latex
        else:
            title = mgr.plot_opts[cn.O_TITLE]
        self.setYAxColor(ax, "left", cn.SIMULATED_COLOR)
        self.setYAxColor(ax2, "right", cn.INPUT_COLOR)
        ax.set_title(title, y=0.2, pad=-14, fontsize=14, loc="right")
        ax.legend([self.output_name, cn.O_PREDICTED], loc="upper left")
        mgr.doPlotOpts()
        mgr.doFigOpts()

    @staticmethod 
    def setYAxColor(ax, position, color):
        # Set the colors of the labels, axes, and spines
        ax.tick_params(axis='y', labelcolor=color)
        ax.spines[position].set_color(color)
        ax.yaxis.label.set_color(color)
    
    def _calculateResiduals(self, transfer_function:control.TransferFunction)->float:
        """
        Computes the residuals for a transfer function specified by parameters.

        Parameters
        ----------
        transfer_function: control.TransferFunction

        Returns
        -------
        float (MSE)
        """
        _, y_arr = self.simulateTransferFunction(transfer_function)
        residuals = self.out_arr - y_arr
        is_bads = [np.isnan(v) or np.isinf(v) or (v is None) for v in residuals]
        if any(is_bads):
            mse = 1e6
        else:
            mse = np.sum(residuals**2)/len(residuals)
        return mse
    
    @staticmethod
    def getNext(initial_value:float, final_value:float, num_increment:int,
                num_entry:int)->itertools.combinations_with_replacement:
        """
        Returns an multiple value array, incrementing the next value in even increments in log space.
        initial_value and final_value must be the same sign.

        Parameters
        ----------
        initial_value: float
        final_value: float
        num_increment: int (number of increments of the values)
        num_entry: int (number of entries in the array)

        Returns
        -------
        itertools.combinations_with_replacement
        """
        if (initial_value*final_value <= 0):
            raise ValueError("initial_value and final_value must be the same sign and non-zero.")
        sign = np.sign(initial_value)
        log_values = np.linspace(np.log10(sign*initial_value), np.log10(sign*final_value), num_increment)
        possible_values = [sign*10**v for v in log_values]
        # Generate evenly spaced values in log space
        return itertools.combinations_with_replacement(possible_values, num_entry)

    @staticmethod
    def uniformFromLogspace(min_value:float, max_value:float, num:int)->np.ndarray:
        """
        Returns an array of values evenly spaced in log space.

        Parameters
        ----------
        min_value: float
        max_value: float
        num: int

        Returns
        -------
        np.ndarray
        """
        if (min_value*max_value <= 0):
            raise ValueError("min_value and max_value must be the same sign and non-zero.")
        sign = np.sign(min_value)
        log_values = np.random.uniform(np.log10(sign*min_value), np.log10(sign*max_value), num)
        values = np.array([sign*10**v for v in log_values])
        return values
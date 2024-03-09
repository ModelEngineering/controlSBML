"""
Abstract class for fitting a SISO transfer function.
"""

import controlSBML.constants as cn
from controlSBML import util
from controlSBML.option_management.option_manager import OptionManager
from controlSBML.timeseries import Timeseries

import control # type: ignore
import numpy as np
import pandas as pd  # type: ignore
import scipy  # type: ignore
from typing import Optional, Tuple, Union, Iterator


MIN_PARAMETER_VALUE = -1e6
MAX_PARAMETER_VALUE = 1e6
MIN_INITIAL_S_VALUE = -1
MAX_INITIAL_S_VALUE = 0
DEFAULT_MIN_POLE_VALUE = -5
DEFAULT_MAX_POLE_VALUE = -0.01
DEFAULT_MIN_ZERO_VALUE = -100
DEFAULT_MAX_ZERO_VALUE = DEFAULT_MAX_POLE_VALUE
DEFAULT_NUM_ITR = 10


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
        num_zero: int (number of zeros)
        num_pole: int (number of poles)
        input_name: str (name of input column)
        output_name: str (name of output column)
        min_pole_value: float
        max_pole_value: float
        min_zero_value: float
        max_zero_value: float
        """
        self.timeseries = timeseries
        self.input_name, self.output_name = self._extractNames(timeseries, input_name, output_name)
        self.in_arr = timeseries[self.input_name].values
        self.out_arr = timeseries[self.output_name].values
        self.max_out_value = np.max(self.out_arr)
        self.dt = np.round(scipy.stats.mode(np.diff(timeseries.index))[0])/cn.MS_IN_SEC
        self.times = np.array([self.dt*n for n in range(len(timeseries.index))])  # Must be evenly spaced
        self.times = np.reshape(self.times, (len(self.times),))
        self.length = len(self.times)
        self.num_zero:int = num_zero  # type: ignore
        self.num_pole:int = num_pole # type: ignore
        self.min_zero_value = self._setValue(min_zero_value, DEFAULT_MIN_ZERO_VALUE)
        self.max_zero_value = self._setValue(max_zero_value, DEFAULT_MAX_ZERO_VALUE)
        self.min_pole_value = self._setValue(min_pole_value, DEFAULT_MIN_POLE_VALUE)
        self.max_pole_value = self._setValue(max_pole_value, DEFAULT_MAX_POLE_VALUE)
        # Outputs of calculations
        self.transfer_function:Union[None, control.TransferFunction] = None

    def copy(self):
        fitter = SISOTransferFunctionFitter(self.timeseries, num_zero=self.num_zero, num_pole=self.num_pole,
                    input_name=self.input_name, output_name=self.output_name,
                    min_pole_value=self.min_pole_value, max_pole_value=self.max_pole_value,
                    min_zero_value=self.min_zero_value, max_zero_value=self.max_zero_value)
        fitter.transfer_function = self.transfer_function
        return fitter

    @staticmethod
    def _setValue(value, default):
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

    def simulateTransferFunction(self,
                transfer_function:Optional[control.TransferFunction]=None)->Tuple[np.ndarray, np.ndarray]:
        """
        Simulates the transfer function. Special handling of negative DCGain to create inverse relationships.

        Parameters
        ----------
        transfer_function: control.TransferFunction

        Returns
        -------
        np.ndarray - times
        np.ndarray - output
        """
        if transfer_function is None:
            transfer_function = self.transfer_function
        sign = np.sign(transfer_function.dcgain())  # type: ignore
        result_input = control.forced_response(sign*transfer_function, T=self.times, U=self.in_arr)
        y_arr_output = np.reshape(result_input.y, (self.length,)) 
        if sign < 0:
            y_arr_output = self.max_out_value - y_arr_output
        return self.times, y_arr_output

    def _calculateTransferFunctionResiduals(self, transfer_function:control.TransferFunction)->float:
        """
        Computes the residuals for a transfer function specified by parameters.

        Parameters
        ----------
        transfer_function: control.TransferFunction

        Returns
        -------
        float (root of the mean square of residuals)
        """
        _, y_arr = self.simulateTransferFunction(transfer_function)
        residuals = self.out_arr - y_arr
        is_bads = [np.isnan(v) or np.isinf(v) or (v is None) for v in residuals]
        if any(is_bads):
            mse = 1e6
        else:
            rmse = np.sqrt(np.mean(residuals**2))
        return rmse
    
    def fit(self, **kwargs)->None:
        """
        Estimates the gain of the transfer function.

        Parameters
        ----------
        kwargs: options for method

        Returns
        -------
        float
        """
        raise NotImplementedError("fit must be implemented in a subclass")
    
    def plot(self, **kwargs):
        """
        Plots the results of fitting a transfer function.
        If an option manager is specified, then the caller handles figure generation.

        Parameters
        ----------
        fitter_result: FitterResult
        mgr: OptionManager
        kwargs: dict
        #@expand
        """
        # Initializations
        mgr = OptionManager(kwargs)
        #
        ts = Timeseries(pd.DataFrame({cn.TIME: self.times, self.output_name: self.out_arr,
                                      cn.O_PREDICTED: self.simulateTransferFunction()[1]}))
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
        ax2.plot(self.times, self.in_arr, color=cn.INPUT_COLOR, linestyle="--")
        latex = util.latexifyTransferFunction(self.transfer_function)
        title = mgr.plot_opts[cn.O_TITLE]
        self._setYAxColor(ax, "left", cn.SIMULATED_COLOR)
        self._setYAxColor(ax2, "right", cn.INPUT_COLOR)
        ax.set_title(title)
        ax.set_title(latex, y=0.2, x=0.5, pad=-14, fontsize=14, loc="right")
        ax.legend([self.output_name, cn.O_PREDICTED], loc="upper left")
        mgr.doPlotOpts()
        mgr.doFigOpts()

    @staticmethod 
    def _setYAxColor(ax, position, color):
        # Set the colors of the labels, axes, and spines
        ax.tick_params(axis='y', labelcolor=color)
        ax.spines[position].set_color(color)
        ax.yaxis.label.set_color(color)
    
    @staticmethod
    def _uniformFromLogspace(min_value:float, max_value:float, num:int)->np.ndarray:
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
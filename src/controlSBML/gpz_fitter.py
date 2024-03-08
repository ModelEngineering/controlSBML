"""
Efficient fitting of transfer functions to data. The fit proceeds in the following steps:
1. Fit the gain of the transfer function.
2. Fit the poles of the transfer function.
3. Fit the zeros of the transfer function.

Public methods:
    copy: Make a copy of the fitter
    fit: Fit the transfer function to the data

Todo:
1. Selective specify minimimum and max values for poles, zeros
2. Allow for a mix of positive and negative poles and zeros
3. Handle poles and zeros in the RHS of the s-plane
"""

import controlSBML.constants as cn
from controlSBML.siso_transfer_function_fitter import SISOTransferFunctionFitter
from controlSBML.option_management.option_manager import OptionManager
from controlSBML.timeseries import Timeseries

import control # type: ignore
from docstring_expander.expander import Expander # type: ignore
import itertools
import lmfit # type: ignore
import numpy as np
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
ZERO_PREFIX = "z"
POLE_PREFIX = "p"
KWARGS = ["num_zero", "num_pole", "input_name", "output_name", "min_pole_value", "max_pole_value",
          "min_zero_value", "max_zero_value", "num_itr"]


class GPZFitter(SISOTransferFunctionFitter):

    def __init__(self, *pargs, num_itr:int=DEFAULT_NUM_ITR, **kwargs):
        """
        Parameters
        ----------
        See SISOTransferFunctionFitter
        """
        diff = set(kwargs.keys()).difference(KWARGS)
        if len(diff) > 0:
            raise ValueError("Unknown keyword arguments: %s" % diff)
        self.num_itr = num_itr
        super().__init__(*pargs, **kwargs)
        # Outputs of calculations
        self.dcgain:Union[None, float] = None
        self.poles:Union[None, np.ndarray[float]] = None
        self.zeros:Union[None, np.ndarray[float]] = None

    def copy(self):
        fitter = self.copy()
        # Outputs of calculations
        fitter.dcgain = self.dcgain
        fitter.poles = self.poles
        fitter.zeros = self.zeros
        fitter.transfer_function = self.transfer_function
        return fitter

    def _makeParameters(self, zeros:Optional[np.ndarray[complex]]=None,
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
        def getValues(values, default_values, min_value, max_value, num_value):
            if values is not None:
                # Specified a list of values
                results = values
            else:
                # No values specified. Use existing if present
                if default_values is not None:
                    results = default_values
                else:
                    results = self._uniformFromLogspace(min_value, max_value, num_value)
            return results
        #
        zeros = getValues(zeros, self.zeros, self.min_zero_value, self.max_zero_value, self.num_zero)
        poles = getValues(poles, self.poles, self.min_pole_value, self.max_pole_value, self.num_pole)
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

    def _makeTransferFunction(self, parameters:lmfit.Parameter):
        """
        Constructs a transfer function from a parameter representation. Gain must have been calculated previously.

        Parameters
        ----------
        parameters: lmfit.Parameter
            name=z<int>: value of a zero
            name=p<int>: value of a pole

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
    
    def fit(self, **kwargs)->None:
        """
        Estimates the gain of the transfer function.

        Parameters
        ----------
        num_itr: int (number of iterations)

        Returns
        -------
        float
        """
        if "num_itr" in kwargs:
            num_itr = kwargs["num_itr"]
        else:
            num_itr = self.num_itr
        self._fitDCGain()
        self._fitPoles(num_itr)
        self._fitZeros(num_itr)
        parameters = self._makeParameters(self.zeros, self.poles)
        self.transfer_function = self._makeTransferFunction(parameters)
    
    def _fitDCGain(self):
        """
        Estimates the gain of the transfer function.

        Returns
        -------
        float
        """
        adj_out_arr = self.out_arr - np.mean(self.out_arr)
        adj_in_arr = self.in_arr - np.mean(self.in_arr)
        self.dcgain = adj_out_arr.dot(adj_in_arr)/adj_in_arr.dot(adj_in_arr)
    
    def _fitPoles(self, num_itr:Optional[int]=None)->None:
        """
        Estimates the poles of the transfer function.

        Returns
        -------
        np.ndarray
        """
        if num_itr is None:
            num_itr = self.num_itr
        self.poles = None  # type: ignore
        best_mse: float = 1e6
        best_poles = np.repeat(None, self.num_pole)  # type: ignore
        for poles in self._getNext(self.min_pole_value, self.max_pole_value, num_itr, self.num_pole):
            parameters = self._makeParameters([], poles)  # type: ignore
            tf = self._makeTransferFunction(parameters)
            mse = self._calculateTransferFunctionResiduals(tf)
            if mse < best_mse:
                best_mse = mse
                best_poles = poles
        self.poles = np.array(best_poles)

    # Fix: Handle a mix of positive and negative zeros
    def _fitZeros(self, num_itr:Optional[int]=None)->None:
        """
        Estimates the zeros of the transfer function.

        Parameters
        ----------
        num_itr: int (number of iterations)

        Returns
        -------
        np.ndarray
        """
        if self.poles is None:
            raise ValueError("Poles must be calculated before calling fitZeros")
        if num_itr is None:
            num_itr = self.num_itr
        # Establish the baseline for comparison
        base_tf = self._makeTransferFunction(self._makeParameters([], self.poles))   # type: ignore
        base_mse = self._calculateTransferFunctionResiduals(base_tf)
        self.zeros = None  # type: ignore
        #
        def search(min_value, max_value):
            best_mse = 1e6
            for zeros in self._getNext(min_value, max_value, num_itr, self.num_zero):
                adj_zeros = [v for v in zeros if v is not None]
                if len(adj_zeros)> len(self.poles):
                    continue
                parameters = self._makeParameters(adj_zeros, self.poles)  # type: ignore
                tf = self._makeTransferFunction(parameters)
                mse = self._calculateTransferFunctionResiduals(tf)
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
                best_zeros = best_zeros1
                best_mse = best_mse1
            else:
                best_zeros = best_zeros2
                best_mse = best_mse2
        else:
            best_zeros, best_mse = search(self.min_zero_value, self.max_zero_value)
        if best_mse < base_mse:
            self.zeros = np.array(best_zeros)
        else:
            self.zeros = np.array([])   # type: ignore
    
    @staticmethod
    def _getNext(initial_value:float, final_value:float, num_increment:int,
                num_entry:int)->Iterator[np.ndarray]:
        """
        Returns an multiple value array, incrementing the next value in even increments in log space.
        initial_value and final_value must be the same sign. The size of the array may be less than num_entry.

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
        possible_values.append(None)
        # Generate evenly spaced values in log space
        for combination in itertools.combinations_with_replacement(possible_values, num_entry):
            refined_combination = np.array([v for v in combination if v is not None])
            yield refined_combination
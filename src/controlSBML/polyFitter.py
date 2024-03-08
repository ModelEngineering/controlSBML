"""
Fits a transfer function by fitting the numerator and denominator polynomials.
"""

import controlSBML.constants as cn
from controlSBML.siso_transfer_function_fitter import SISOTransferFunctionFitter

import control # type: ignore
import lmfit # type: ignore
import numpy as np
from typing import Optional


MIN_PARAMETER_VALUE = -1e6
MAX_PARAMETER_VALUE = 1e6
INITIAL_PARAMETER_VALUE = 0.1
NUMERATOR_PREFIX = "n"
DENOMINATOR_PREFIX = "d"


class polyFitter(SISOTransferFunctionFitter):

    def __init__(self, *args, **kwargs):
        """
        """
        super().__init__(*args, **kwargs)
        self.num_numerator = self.num_zero + 1
        self.num_denominator = self.num_pole + 1
    
    def copy(self):
        fitter = self.copy()
        # Outputs of calculations
        fitter.num_numerator = self.num_numerator
        fitter.num_denominator = self.num_denominator
        return fitter
    
    def makeParameters(self, num_numerator:Optional[int]=None, num_denominator:Optional[int]=None,
                       initial_value:Optional[float]=None):
        """
        Makes the parameters used to use lmfit to search for a best fitting transfer function.

        Parameters
        ----------
        num_numerator_term: int (number of terms in the zero polynomial)
        num_denominator: int (number of poles in the transfer function)
        initial_value: float (initial value for the parameters)

        Returns
        -------
        lmfit.Parameter for zeros (begins with 'z') and poles (begins with 'p')
        """
        def getValue():
            if initial_value is None:
                value = np.random.uniform(MIN_PARAMETER_VALUE, MAX_PARAMETER_VALUE)
                return value
            else:
                return initial_value
        #   
        pfit = lmfit.Parameters()
        if num_numerator is None:
            num_numerator = self.num_numerator
        if num_denominator is None:
            num_denominator = self.num_denominator
        #
        # Numerator
        for idx in range(num_numerator):
            value = getValue()
            pfit.add(name='%s%d' % (NUMERATOR_PREFIX, idx+1),
                min=MIN_PARAMETER_VALUE,
                value=value,
                max=MAX_PARAMETER_VALUE)
        # Poles
        for idx in range(num_numerator):
            value = getValue()
            pfit.add(name='%s%d' % (DENOMINATOR_PREFIX, idx+1),
                min=MIN_PARAMETER_VALUE,
                value=value,
                max=MAX_PARAMETER_VALUE)
        return pfit

    def fit(self, **kwargs):
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
        num_numerator = kwargs.get("num_numerator", self.num_denominator)
        num_denominator = kwargs.get("num_denominator", self.num_denominator)
        initial_value = kwargs.get("initial_value", self.initial_value)
        parameters = self.makeParameters(num_numerator, num_denominator, initial_value)
        minimizer_result = lmfit.minimize(self._calculateResiduals, parameters,
                                          method="differential_evolution")
        self.transfer_function = self.makeTransferFunction(minimizer_result.params)
    
    def _calculateResiduals(self, parameters):
        """
        Computes the residuals for a transfer function specified by parameters.

        Parameters
        ----------
        parameters: lmfit.Parameters

        Returns
        -------
        float (MSE)
        """
        transfer_function = self._makeTransferFunction(parameters)
        return super()._calculateTransferFunctionResiduals(transfer_function)

    def _makeTransferFunction(self, parameters:lmfit.Parameter):
        """
        Constructs a transfer function from a parameter representation. Gain must have been calculated previously.

        Parameters
        ----------
        parameters: lmfit.Parameter
            name=n<int>: numerator coefficient for the s**n
            name=d<int>: denominator coefficient for s**n

        Returns
        -------
        control.TransferFunction
        """
        s = control.TransferFunction.s
        numerator = 0
        denominator = 0
        for key, value in parameters.valuesdict().items():
            spower = int(key[1:])
            if key[0] == NUMERATOR_PREFIX:
                numerator += value*s**spower
            elif key[0] == DENOMINATOR_PREFIX:
                numerator += value*s**spower
            else:
                raise ValueError("Unknown key in transfer function parameters: %s" % key)
        tf = numerator/denominator
        return tf
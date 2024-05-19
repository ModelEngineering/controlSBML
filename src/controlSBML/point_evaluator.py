"""Evaluates a single design point."""

import controlSBML.constants as cn
from controlSBML.parallel_search import Evaluator
from controlSBML.sbml_system import SBMLSystem

import numpy as np
from typing import Tuple, Union


##################################################################
class PointEvaluator(Evaluator):
    # Evaluates a point in the design space
    def __init__(self, sbml_system:SBMLSystem, input_name:str, output_name:str, setpoint:float, 
                 sign:float, times:np.array, is_greedy:bool=False):
        """
        Args:
            sbml_system (SBMLSystem): Should be a copy of the original system so it is light weight
            input_name (str): system input
            output_name (str): system output
            setpoint (float):
            times (np.array):
            is_greedy (bool): if True, then a greedy search is done to find a feasible system
        """
        self.sbml_system = sbml_system
        if sbml_system.isInitialized():
            raise ValueError("sbml_system should be a copy of the original system so it is light weight")
        self.input_name = input_name
        self.output_name = output_name
        self.setpoint = setpoint
        self.sign = sign
        self.times = times
        self.siso_evaluator = None
        self.is_greedy = is_greedy

    def initialize(self)->None:
        pass

    def evaluate(self, **candidate_dct:dict)->dict:
        """
        Calculates MSE for the specified points.

        Args:
            candidate_dct: dict (key: name of parameter, value: value of parameter)
        Returns:
            dict
                score: float (MSE)
                cadidate_dct keys, values
        """
        # FIXME
#        if self.is_greedy:
#            new_dct = self._searchForFeasibleClosedLoopSystem(**candidate_dct)
#            if new_dct is None:
#                result_dct[cn.SCORE] = None
#                return result_dct
#        else:
#            new_dct = dict(candidate_dct)
        new_dct = dict(candidate_dct)
        result_dct = dict(new_dct)
        reason, residual_mse = self._calculateMse(**new_dct)   # type: ignore
        result_dct[cn.SCORE] = residual_mse # type: ignore
        result_dct[cn.REASON] = reason  # type: ignore
        return result_dct

    def _calculateMse(self, **parameter_dct:dict)->Tuple[str, object]:
        """
        Attempts to calculate the mean squared error of the closed loop system. Reports if system is unstable.

        Args:
            parameter_dct: dict: {name: value for each of kP, kI, kF}

        Returns:
            str (description of design result)
            float (mean squared error)
        """
        SIZE_MARGIN = 1e6
        _ = self.sbml_system.roadrunner  # Initialize roadrunner
        max_output = SIZE_MARGIN*self.sbml_system._max_value_dct[self.output_name]
        min_output = self.sbml_system._min_value_dct[self.output_name]/SIZE_MARGIN
        try:
            response_ts, builder = self.sbml_system.simulateSISOClosedLoop(setpoint=self.setpoint,
                        input_name=self.input_name, output_name=self.output_name,
                        times=self.times, sign=self.sign,
                        is_steady_state=self.sbml_system.is_steady_state, inplace=False,
                        **parameter_dct)  # type: ignore
        except Exception as error:
            if "CVODE" in str(error):
                return cn.DESIGN_RESULT_CANNOT_SIMULATE, None
            else:
                raise ValueError(str(error))
        # Check for large outputs
        if response_ts is None:
            return cn.DESIGN_RESULT_CANNOT_SIMULATE, None
        outputs = response_ts[self.output_name].values
        max_value = np.max([np.max(outputs), np.abs(np.min(outputs))])
        if False:
            # Disable these checks
            if max_value > max_output:
                return cn.DESIGN_RESULT_OUTPUT_TOO_LARGE, None
            min_value = np.min([np.max(outputs), np.abs(np.min(outputs))])
            if min_value < min_output:
                return cn.DESIGN_RESULT_OUTPUT_TOO_SMALL, None
        #
        residuals = self.setpoint - response_ts[self.output_name].values
        mse = np.mean(residuals**2)
        return cn.DESIGN_RESULT_SUCCESS, mse
    
    def _searchForFeasibleClosedLoopSystem(self, max_iteration:int=10, **parameter_dct)->Union[dict, None]:
        """
        Does a greedy search to find values of the parameters that are minimally stable.

        Args:
            max_iteration: int (maximum number of iterations)
            parameter_dct: dict (name: value)
                key: name of parameter
                value: value of parameter or None
        Returns:
            dict: {name: value} or None (no stable result found)
                key: name of parameter
                value: value of parameter or None
        """
        MINIMAL_FACTOR = 0.01
        factor = 0.5
        def mult(dct, factor):
            new_dct = {}
            for name, value in dct.items():
                if value is None:
                    new_dct[name] = None
                else:
                    new_dct[name] = factor*value
            return new_dct
        # Iterate to find values
        dct = dict(parameter_dct)
        last_stable_dct = None
        for _ in range(max_iteration):
            is_feasible, _ = self._calculateMse(**dct)
            if is_feasible:
                if factor < MINIMAL_FACTOR:
                    break
                else:
                    # Try a bigger factor
                    factor = 1 + factor
                    last_stable_dct = dict(dct)
            else:
                if last_stable_dct is not None:
                    dct = dict(last_stable_dct)
                    break
                else:
                    factor = factor/2
            dct = mult(dct, factor)
        # Return the result
        if last_stable_dct is None:
            return None
        return dct
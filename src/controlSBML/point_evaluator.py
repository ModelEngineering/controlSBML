"""Evaluates a single design point."""

import controlSBML.constants as cn
from controlSBML.parallel_search import Evaluator
from controlSBML.sbml_system import SBMLSystem

import numpy as np


##################################################################
class PointEvaluator(Evaluator):
    # Evaluates a point in the design space
    def __init__(self, sbml_system:SBMLSystem, input_name:str, output_name:str, setpoint:float, times:np.array,
                 is_greedy:bool=False):
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
        #
        if self.is_greedy:
            new_dct = self._searchForFeasibleClosedLoopSystem(**candidate_dct)
            if new_dct is None:
                result_dct[cn.SCORE] = None
                return result_dct
        else:
            new_dct = dict(candidate_dct)
        result_dct = dict(new_dct)
        _, residual_mse = self._calculateMse(**new_dct)
        result_dct[cn.SCORE] = residual_mse
        return result_dct

    def _calculateMse(self, max_output:float=1e6, min_output:float=0, **parameter_dct:dict)->(bool, float):
        """
        Attempts to calculate the mean squared error of the closed loop system. Reports if system is unstable.

        Args:
            max_output: float (maximum output)
            min_output: float (minimum output)
            parameter_dct: dict: {name: value for eack of kp, ki, kf}

        Returns:
            bool (successful simulation)
            float (mean squared error)
        """
        try:
            response_ts, _ = self.sbml_system.simulateSISOClosedLoop(setpoint=self.setpoint,
                        input_name=self.input_name, output_name=self.output_name,
                        times=self.times,
                        is_steady_state=self.sbml_system.is_steady_state, inplace=False,
                        **parameter_dct)
        except Exception as error:
            if "CVODE" in str(error):
                return False, None
            else:
                raise ValueError(str(error))
        # Check for large outputs
        outputs = response_ts[self.output_name].values
        max_value = np.max([np.max(outputs), np.abs(np.min(outputs))])
        min_value = np.min([np.max(outputs), np.abs(np.min(outputs))])
        if (max_value > max_output) or (min_value < min_output):
            return False, None
        # Check for negative values
        for column in response_ts.columns:
            if np.any(response_ts[column].values < 0):
                return False, None
        #
        residuals = self.setpoint - response_ts[self.output_name].values
        mse = np.mean(residuals**2)
        return True, mse
    
    def _searchForFeasibleClosedLoopSystem(self, max_iteration:int=10, **parameter_dct)->dict:
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
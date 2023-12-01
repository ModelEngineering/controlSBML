"""Evaluates a closed loop design."""

from controlSBML import util
from controlSBML.sbml_system import SBMLSystem
import controlSBML.constants as cn

import numpy as np
from typing import List


class SISODesignEvaluator:
    # Evaluates designs and remembers the best one

    def __init__(self, system:SBMLSystem, input_name, output_name, setpoint:float=1, times:List[float]=cn.TIMES):
        self.system = system
        self.setpoint = setpoint
        self.times = times
        self.input_name = input_name
        self.output_name = output_name
        #
        self.kp = None
        self.ki = None
        self.kf = None
        self.residual_mse = None

    def evaluate(self, kp:float=None, ki:float=None, kf:float=None):
        """
        Evaluates the closed loop system. Updates the control parameters if the design is better.

        Args:
            kp: float (proportional gain)
            ki: float (integral gain)
            kf: float (feedforward gain)
        Returns:
            bool (successful simulation)
        """
        def update(mse, kp=kp, ki=ki, kf=kf):
            self.kp = kp
            self.ki = ki
            self.kf = kf
            self.residual_mse = mse
        #
        value_dct = {cn.CP_KP: kp, cn.CP_KI: ki, cn.CP_KF: kf}
        is_feasible, mse = self._calculateMse(**value_dct)
        if not is_feasible:
            return False
        if self.residual_mse is None:
            update(mse, **value_dct)
        elif mse < self.residual_mse:
            update(mse, **value_dct)
        return True
    
    def _calculateMse(self, max_output:float=1e6, min_output:float=0, **parameter_dct:dict):
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
            response_ts, _ = self.system.simulateSISOClosedLoop(setpoint=self.setpoint,
                        input_name=self.input_name, output_name=self.output_name,
                        times=self.times,
                        is_steady_state=self.system.is_steady_state, inplace=False,
                        **parameter_dct)
        except Exception:
            return False, None
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
"""Evaluates a closed loop design."""

from controlSBML import util
from controlSBML.sbml_system import SBMLSystem
import controlSBML.constants as cn

import numpy as np
import pandas as pd
import os
from typing import List

# Columns: kp, ki, kf, mse


class SISODesignEvaluator:
    # Evaluates designs and remembers the best one

    def __init__(self, system:SBMLSystem, input_name, output_name, setpoint:float=1, times:List[float]=cn.TIMES,
                 save_path=None):
        self.system = system
        self.setpoint = setpoint
        self.times = times
        self.input_name = input_name
        self.output_name = output_name
        self.save_path = save_path
        #
        self.kp = None
        self.ki = None
        self.kf = None
        self.residual_mse = None

    def getEvaluatorResults(self)->pd.DataFrame:
        """
        Returns the design results if they have been saved.

        Returns:
            pd.DataFrame
                columns: kp, ki, kf, mse
        """
        if (self.save_path is not None) and (os.path.isfile(self.save_path)):
            df = pd.read_csv(self.save_path)
            return df
        else:
            return None

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
        if (self.save_path is not None) and (os.path.isfile(self.save_path)):
            df = pd.read_csv(self.save_path)
            result_dct = df.to_dict(orient="list")
        else:
            result_dct = {cn.CP_KP: [], cn.CP_KI: [], cn.CP_KF: [], cn.MSE: []}
        def update(mse, kp=kp, ki=ki, kf=kf):
            self.kp = kp
            self.ki = ki
            self.kf = kf
            self.residual_mse = mse
        #
        # Evaluate the design
        value_dct = {cn.CP_KP: kp, cn.CP_KI: ki, cn.CP_KF: kf}
        is_feasible, mse = self._calculateMse(**value_dct)
        # Save the results
        if self.save_path is not None:
            result_dct[cn.CP_KP].append(kp)
            result_dct[cn.CP_KI].append(ki)
            result_dct[cn.CP_KF].append(kf)
            result_dct[cn.MSE].append(mse)
            df = pd.DataFrame(result_dct)
            df.to_csv(self.save_path, index=False)
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
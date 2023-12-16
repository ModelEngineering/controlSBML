"""Evaluates a closed loop design. Saves the best design using MSE criteria. Preserves all designs."""

from controlSBML import util
from controlSBML.sbml_system import SBMLSystem
import controlSBML.constants as cn

import numpy as np
import pandas as pd
import os
from typing import List


class EvaluatorResult(dict):
    # Container of evaluation results

    def __init__(self, kp_spec=False, ki_spec=False, kf_spec=False):
        self[cn.MSE] = []
        self.attrs = [cn.MSE]
        if kp_spec:
            self[cn.CP_KP] = []
            self.attrs.append(cn.CP_KP)
        if ki_spec:
            self[cn.CP_KI] = []
            self.attrs.append(cn.CP_KI)
        if kf_spec:
            self[cn.CP_KF] = []
            self.attrs.append(cn.CP_KF)

    def __len__(self):
        keys = list(self.keys())
        if len(keys) == 0:
            return 0
        return len(self[keys[0]])

    def add(self, **kwargs):
        for key, value in kwargs.items():
            self[key].append(value)

    @classmethod
    def makeFromDataframe(self, dataframe):
        """
        Creates a result object from a dataframe

        Args:
            dataframe (dataframe): columns: kp, ki, kf, mse

        Returns:
            EvaluatorResult
        """
        kp_spec = cn.CP_KP in dataframe.columns
        ki_spec = cn.CP_KI in dataframe.columns
        kf_spec = cn.CP_KF in dataframe.columns
        evaluator_result = EvaluatorResult(kp_spec=kp_spec, ki_spec=ki_spec, kf_spec=kf_spec)
        for _, row in dataframe.iterrows():
            kwargs = {cn.MSE: row[cn.MSE]}
            if kp_spec:
                kwargs[cn.CP_KP] = row[cn.CP_KP]
            if ki_spec:
                kwargs[cn.CP_KI] = row[cn.CP_KI]
            if kf_spec:
                kwargs[cn.CP_KF] = row[cn.CP_KF]
            evaluator_result.add(**kwargs)
        return evaluator_result
    
    def writeCsv(self, file_path):
        df = self.getDataframe()
        df.to_csv(file_path, index=False)

    def getDataframe(self):
        return pd.DataFrame(self)

    def copy(self):
        evaluator_result = EvaluatorResult()
        evaluator_result[cn.MSE] = list(self[cn.MSE])
        if cn.CP_KP in self:
            evaluator_result[cn.CP_KP] = list(self[cn.CP_KP])
        if cn.CP_KI in self:
            evaluator_result[cn.CP_KI] = list(self[cn.CP_KI])
        if cn.CP_KF in self:
            evaluator_result[cn.CP_KF] = list(self[cn.CP_KF])
        return evaluator_result

    @classmethod
    def merge(cls, evaluator_results):
        """
        Merges multiple evaluator results.

        Args:
            evaluator_results (list): list of EvaluatorResult objects

        Returns:
            EvaluatorResult
        """
        template_result = evaluator_results[0]
        merged_result = {k: [] for k in template_result.keys()}
        for key in merged_result.keys():
            for results in evaluator_results:
                merged_result[key].extend(results[key])
        return merged_result


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
        # Best results
        self.kp = None
        self.ki = None
        self.kf = None
        self.residual_mse = None
        # All results
        self.evaluator_result = None

    @classmethod
    def merge(cls, evaluators):
        """
        Merges multiple evaluators.

        Args:
            evaluators (_type_): _description_
        """
        def update(evaluator, merged_evaluator):
            merged_evaluator.kp = evaluator.kp
            merged_evaluator.ki = evaluator.ki
            merged_evaluator.kf = evaluator.kf
            merged_evaluator.residual_mse = evaluator.residual_mse
        #
        merged_evaluator = evaluators[0]
        #merged_evaluator = merged_evaluator.copy(is_set_outputs=True)
        merged_evaluator = merged_evaluator.copy()
        results = [e.evaluator_result for e in evaluators]
        merged_evaluator.evaluator_result = EvaluatorResult.merge(results)
        df = pd.DataFrame(merged_evaluator.evaluator_result)
        ser = df[cn.MSE].dropna()
        min_mse = np.min(ser)
        idx = ser[ser == min_mse].index[0]
        if cn.CP_KP in df.columns:
            merged_evaluator.kp = df[cn.CP_KP].values[idx]
        if cn.CP_KI in df.columns:
            merged_evaluator.ki = df[cn.CP_KI].values[idx]
        if cn.CP_KF in df.columns:
            merged_evaluator.kf = df[cn.CP_KF].values[idx]
        if False:
            for evaluator in evaluators:
                if evaluator.residual_mse is None:
                    continue
                if merged_evaluator.residual_mse is None:
                    update(evaluator, merged_evaluator)
                    continue
                if evaluator.residual_mse < merged_evaluator.residual_mse:
                    update(evaluator, merged_evaluator)
                    continue
        return merged_evaluator
        
    def copy(self, is_set_outputs:bool=True):
        """
        Args:
            is_set_outputs (bool, optional): copy the outputs. Defaults to True.
        """
        evaluator = SISODesignEvaluator(self.system.copy(), self.input_name, self.output_name, self.setpoint, self.times,
                                   self.save_path)
        if evaluator.evaluator_result is not None:
            evaluator.evaluator_result = self.evaluator_result.copy()
        else:
            evaluator.evaluator_result = None
        if is_set_outputs:
            evaluator.kp = self.kp
            evaluator.ki = self.ki
            evaluator.kf = self.kf
            evaluator.residual_mse = self.residual_mse
        return evaluator

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
        def update(**kwargs):
            for key, value in kwargs.items():
                setattr(self, key, value)
        #
        # Initialization
        if self.evaluator_result is None:
            if (self.save_path is not None) and (os.path.isfile(self.save_path)):
                df = pd.read_csv(self.save_path)
                self.evaluator_result = EvaluatorResult.makeFromDataframe(df)
            else:
                self.evaluator_result = EvaluatorResult(kp_spec=kp, ki_spec=ki, kf_spec=kf)
        # Evaluate the design
        value_dct = {cn.CP_KP: kp, cn.CP_KI: ki, cn.CP_KF: kf}
        original_value_dct = dict(value_dct)
        for key, value in original_value_dct.items():
            if value is None:
                del value_dct[key]
        is_feasible, residual_mse = self.calculateMse(**value_dct)
        # Save the results
        self.evaluator_result.add(**value_dct, mse=residual_mse)
        self.evaluator_result.writeCsv(self.save_path)
        if not is_feasible:
            return False
        if self.residual_mse is None:
            update(residual_mse=residual_mse, **value_dct)
        elif residual_mse < self.residual_mse:
            update(residual_mse=residual_mse, **value_dct)
        return True
    
    def calculateMse(self, max_output:float=1e6, min_output:float=0, **parameter_dct:dict):
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
    
    @classmethod
    def makeFromDataframe(cls, system:SBMLSystem, input_name:str, output_name:str, dataframe:pd.DataFrame,
                setpoint:float=1, times:List[float]=cn.TIMES,
                save_path=None):
        """
        Creates a SISODesignEvaluator object from a dataframe

        Args:
            dataframe (dataframe): columns: kp, ki, kf, mse

        Returns:
            EvaluatorResult
        """
        evaluator = cls(system, input_name, output_name,
                setpoint=setpoint, times=times, save_path=save_path)
        evaluator_result = EvaluatorResult.makeFromDataframe(dataframe)
        df = dataframe.sort_values(cn.MSE)
        if cn.CP_KP in df.columns:
            evaluator.kp = df[cn.CP_KP].values[0]
        if cn.CP_KI in df.columns:
            evaluator.ki = df[cn.CP_KI].values[0]
        if cn.CP_KF in df.columns:
            evaluator.kf = df[cn.CP_KF].values[0]
        evaluator.residual_mse = df[cn.MSE].values[0]
        return evaluator
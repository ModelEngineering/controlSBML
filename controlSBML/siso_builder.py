"""Builds and Evaluates SISO Systems"""

"""
The SISO system has:
  * A ControlSBML SISO System
  * PID Controller
  * Filter (with dcgain = 1)

The design parameters of the system are:
  * PID parameters: kp, ki, kd
  * Filter parameter: constant

Evaluations provided are:
  * Find dcgain of ControlSBML system with many inputs and one output
  * DC gain of the closed loop system
  * Evaluate the closed loop system for disturbances and noise
"""

import controlSBML as ctl
from controlSBML import msgs

import control
import numpy as np
import pandas as pd


class SISOBuilder(object):

    def __init__(self, ctlsb):
        self.ctlsb = ctlsb

    def evaluateControllability(self, times, input_names=None, output_names=None):
        """
        Evaluates the controllability of the inputs on the outputs.
        If no input (output) is specified, then all inputs (outputs)
        are considered for self.ctlsb.

        Parameters
        ----------
        times: list-float
             times at which dc gain is evaluated
        input_names: list-str
        output_names: list-str
        
        Returns
        -------
        dict
            key: time
            value: pd.DataFrame
                 index: input_names
                 columns: output_names
                 value: dcgain of output for input
        """
        if input_names is None:
            input_names = self.ctlsb.input_names
        if output_names is None:
            output_names = self.ctlsb.output_names
        dct = {t: {} for t in times}
        for output_name in output_names:
            for input_name in input_names:
                ctlsb = ctl.ControlSBML(self.ctlsb.model_reference,
                      input_names=[input_name],
                      output_names=[output_names],
                      is_reduced=self.ctlsb.is_reduced)
                 dct[time][output_name] = []
                for time in times:
                     tf = ctlsb.makeTransferFunction(time)
                     dct[time][output_name].append(tf.dcgain())
        # Construct the DataFrames
        result_dct = {}
        for time in dct.keys():
            result_dct[time] = pd.DataFrame(dct[time], index=input_names)
        #
        return result_dct

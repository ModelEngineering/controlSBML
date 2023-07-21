"""
Builds transfer functions for an SBML model from all inputs to all outputs.
"""

import controlSBML as ctl
import controlSBML.constants as cn
from controlSBML import util
import controlSBML.siso_transfer_function_builder as tfb
import controlSBML.simulate_system as ss
import controlSBML.timeseries as ts
from controlSBML.option_management.option_manager import OptionManager
from controlSBML.option_management.options import Options

import collections
import control
from docstring_expander.expander import Expander
import numpy as np
import pandas as pd


class SBMLTransferFunctionBuilder(object):
    
    def __init__(self, ctlsb):
        """
        Construction of an array of control.TransferFunction

        Parameters
        ----------
        ctlsb: ControlSBML
        """
        self.ctlsb = ctlsb
        self.sys = ctlsb.makeNonlinearIOSystem("SBMLTransferFunctionBuilder")
        
    def fitTransferFunction(self, num_numerator, num_denominator, is_tf_only=True, **kwargs):
        """
        Constructs transfer functions for the NonlinearIOSystem.

        Parameters
        ----------
        num_numerator: int (number of numerator terms)
        num_denominator: int (number of denominator terms)
        is_tf_only: bool (Values are control.TransferFunction)
        kwargs: dict
            num_step: int (number of steps in staircase)
            initial_value: float (value for first step)
            final_value: float (value for final step)
        #@expand

        Returns
        -------
        DataFrame: 
            column names: str (output)
            index: str (input)
            values: tfb.FitterResult or control.TransferFunction
        """
        result_dct = {n: [] for n in self.sys.output_names}
        for output_name in self.sys.output_names:
            for input_name in self.sys.input_names:
                sys = self.sys.getSubsystem(self.sys.name, [input_name], [output_name])
                siso_tfb = tfb.SISOTransferFunctionBuilder(sys)
                value = siso_tfb.fitTransferFunction(num_numerator, num_denominator, **kwargs)
                if is_tf_only:
                    value = value.transfer_function
                result_dct[output_name].append(value)
        # Construct the output
        df = pd.DataFrame(result_dct)
        df.columns.name = "Outputs"
        df.index = list(self.sys.input_names)
        df.index.name = "Inputs"
        return df
    
    @classmethod
    def makeBuilder(cls, *pargs, **kwargs):
        """
        Constructs transfer functions for SBML systems.

        Parameters
        ----------
        Same as for constructing ControlSBML
        """
        ctlsb = ctl.ControlSBML(*pargs, **kwargs)
        return cls(ctlsb)
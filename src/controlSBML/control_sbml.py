"""
Top level API

Usage:

    # Examine an SBML model
    ctlsb = ControlSBML(model)
    ctlsb.setTimes(np.linspace(0, 100, 1000))
    ctlsb.plotModel()
    # Define the system
    ctlsb.getPossibleInputs()
    ctlsb.getPossibleOutputs()
    ctlsb.makeSystem(input_names=["S1"], output_names=["S3"])
    # Evaluate controllability --do the inputs affect the outputs?
    ctlsb.setStaircase(input_value=0, final_value=10, num_step=5)
    ctlsb.plotStaircaseResponse() # Displays the outputs associated with the inputs
    # Identify the system by fitting to the input staircase
    fitter_result = ctlsb.fitTransferFunction()
    ctlsb.plotFitterResult(fitter_result)
    # Construct a SISO closed loop system and do a trail and error design
    ctlsb.makeSISOClosedLoopSystem(kp=1, ki=0.1, kf=0.5, setpoint=3)
    _ = ctlsb.simulateClosedLoopSystem()  # Plots the results by default
    ctlsb.makeSISOClosedLoopSystem(kp=0.8, ki=0.2, kf=0.1, setpoint=5)
    _ = ctlsb.simulateClosedLoopSystem()
    # Automatically design
    ctlsb.design(kp=True, ki=True, min_value=0, max_value=10, num_iteration=20, setpoint=4)
    _ = ctlsb.simulateClosedLoopSystem()
    # Explore deviations on the design
    param_dct = ctlsb.get()  # Returns a dictionary of parameter values
    ctlsb.makeSISOClosedLoopSystem(setpoint=4, **param_dct)  # How well does it do with a different setpoint
    new_param_dct = {n: v*1.1 for n, v param_dct.items()}
    ctlsb.makeSISOClosedLoopSystem(setpoint=5, **new_param_dct)  # Do further trail and error


"""

from controlSBML.sbml_system import SBMLSystem
from controlSBML.siso_transfer_function_builder import SISOTransferFunctionBuilder
from controlSBML.siso_closed_loop_designer import SISOClosedLoopDesigner
from controlSBML.timeseries import Timeseries
from controlSBML.staircase import Staircase
from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML import util
import controlSBML.constants as cn

import numpy as np


class ControlSBML(object):

    def __init__(self, model_reference):
        """
        model_reference: str
            string, SBML file or Roadrunner object
        is_fixed_input_species: bool (input species are fixed)
        is_steady_state: bool (start the simulation at steady state)
        """
        # First initializations
        self.model_reference = model_reference
        self.roadrunner = makeRoadrunner(self.model_reference)
        self.setTimes()   # self.times
        self.sbml_system = None
        # Staircase descriptions
        self.staircase_initial_value = cn.DEFAULT_INITIAL_VALUE
        self.staircase_num_step = cn.DEFAULT_NUM_STEP
        self.staircase_final_value = cn.DEFAULT_FINAL_VALUE

    def makeSystem(self, input_names=None, output_names=None, is_fixed_input_species=True, is_steady_state=False):
        """
        Creates an SBML system for the specified inputs and outputs.

        Args:
            input_names: list-str
            output_names: list-str
            is_fixed_input_species: bool (if True, a species used as input is treated as fixed)
            is_fixed_steady_state: bool (if True, simulations are started in steady state)
        """
        self.sbml_system = SBMLSystem(self.roadrunner, input_names=input_names, output_names=output_names,
                                      is_fixed_input_species=is_fixed_input_species, is_steady_state=is_steady_state)

    def plotModel(self, **kwargs):
        """
        Plots the original model.

        Args:
            kwargs: plot options
        
        Returns:
            Timeseries
        """
        data = self.roadrunner.simulate(self.times[0], self.times[-1], len(self.times))
        ts = Timeseries(data)
        util.plotOneTS(ts, **kwargs)
        return ts
        
    def get(self):
        """
        Provides the parameters of the closed loop system
        """
        if self.sbml_system is None:
            raise ValueError("Must use 'makeSystem' before 'get'.")
        return self.sbml_system.get()
    
    def setTimes(self, times=None):
        if times is None:
            times = np.linspace(cn.START_TIME, cn.END_TIME, cn.POINTS_PER_TIME*(cn.END_TIME-cn.START_TIME))
        self.times = times

    def copy(self):
        ctlsb = ControlSBML(self.model_reference)
        ctlsb.times = self.times
        if self.sbml_system is not None:
            ctlsb.makeSystem(input_names=self.sbml_system.input_names, output_names=self.sbml_system.output_names,
                             is_fixed_input_species=self.sbml_system.is_fixed_input_species,
                             is_steady_state=self.sbml_system.is_steady_state)
        ctlsb.staircase_initial_value = self.staircase_initial_value
        ctlsb.staircase_num_step = self.staircase_num_step
        ctlsb.staircase_final_value = self.staircase_final_value

    @property
    def input_staircase(self):
        """
        Creates a staircase to use for system evaluations.
        """
        return Staircase(initial_value=self.staircase_initial_value, final_value=self.staircase_final_value,
                         num_step=self.staircase_num_step, num_point=len(self.times))
    
    def setStaircase(self, initial_value=cn.DEFAULT_INITIAL_VALUE, final_value=cn.DEFAULT_FINAL_VALUE, num_step=cn.DEFAULT_NUM_STEP):
        self.staircase_initial_value = initial_value
        self.staircase_final_value = final_value
        self.staircase_num_step = num_step
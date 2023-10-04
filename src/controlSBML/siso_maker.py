"""
Makes a SISO model of a system with plots and defaults for staircase and closed loop designs.

Usage:
    from controlSBML.siso_maker import SISOMaker
    maker = SISOMaker(model_id, model)
    maker.makeStaircase()
    maker.makeClosedLoop()

    To run on all of BioModels:
        SISOMaker.runBiomodels(start=1)
"""

from controlSBML import constants as cn
from controlSBML.sbml_system import SBMLSystem
from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML.iterate_biomodels import iterateBiomodels
from controlSBML.make_roadrunner import makeRoadrunner

import matplotlib.pyplot as plt
import numpy as np
import os
import tellurium as te

INPUT_COLOR = "red"
OUTPUT_COLOR = "blue"
STAIRCASE_SFX = "_staircase"
CLOSED_LOOP_SFX = "_closed_loop"
PNG_EXT = ".png"


class SISOMaker(object):

    def __init__(self, model_id, model, input_name=None, output_name=None, end_time=20):
        """
        Args:
            model_id: str (identifier of model)
            model: model reference
            input_name: str
            output_name: str
            end_time: float
        """
        self.model_id = model_id
        self.model = model
        self.end_time = end_time
        #
        try:
            self.roadrunner = makeRoadrunner(model)
        except Exception as error:
            print("Error for %s: %s" % (self.filename, error))
            return
        self.input_name, self.output_name = self._makeInputOutput(input_name, output_name)
        self.system, self.times = self._makeSystem()
        self.times = None

    def _makeInputOutput(self, input_name, output_name):
        def getValidNames(names):
            return [n for n in names if not n in ["at", "in"]]
        #
        if (input_name is None) or (output_name is None):
            names = self.roadrunner.getFloatingSpeciesIds()
            names.extend(self.roadrunner.getGlobalParameterIds())
            names.extend(self.roadrunner.getAssignmentRuleIds())
            valid_names = getValidNames(names)
        if input_name is None:
            if len(valid_names) < 1:
                print("No valid_names for %s" % self.filename)
                return
            new_input_name = valid_names[0]
        if output_name is None:
            if len(valid_names) < 2:
                print("Only one floating species for %s" % filename)
                return
            new_output_name = valid_names[1]
        return new_input_name, new_output_name

    def _makeSystem(self):
        """
        Creates a SISO system for the model.

        Returns:
            SBMLSystem
            list-float (times)
        """
        system = SBMLSystem(self.roadrunner, [self.input_name], [self.output_name], is_fixed_input_species=False)
        times = np.linspace(0, self.end_time, 10*self.end_time)
        return system, times
    
    def makeStaircase(self):
        pass

    def makeClosedLoop(self):
        pass

    @staticmethod
    def _checkFile(filename, _, suffix=STAIRCASE_SFX):
        """Checks if the file has already been processed.
        Args:
            filename: str (filename with extension)
            extension: str (extension to check for)
        Returns:
            str
        """
        base_filename = iterateBiomodels._makeBasefilename(filename)
        plot_name = base_filename + suffix + PNG_EXT
        if os.path.isfile(plot_name):
            return ("Already processed")
        else:
            return ""
        
    def _checkStaircase(self):
        return self._checkFile(suffix=STAIRCASE_SFX)
        
    def _checkClosedloop(self):
        return self._checkFile(suffix=CLOSED_LOOP_SFX)

    def makeStaircase(self, initial_value=0, final_value=10):
        """
        Creates a staircase response and a plot. Creates the file
        <filename>_staircase.png.

        Args:
            initial_value (int, optional): _description_. Defaults to 0.
            final_value (int, optional): _description_. Defaults to 10.
        """
        try:
            ts = self.system.simulateStaircase(self.input_name, self.output_name, times=self.times,
                                               final_value=final_value, num_step=5, is_steady_state=False)
        except Exception as error:
            print("Error in staricase for %s: %s" % (self.filename, error))
            return
        # Make the plot
        plot_filename = self.model_id + STAIRCASE_SFX + PNG_EXT
        inputs = ts[self.input_name].values
        outputs = ts[self.output_name].values
        times = ts.index/1000
        _, ax = plt.subplots(1)
        ax.scatter(times, outputs, color=OUTPUT_COLOR)
        ax2 = ax.twinx()
        ax2.plot(times, inputs, color=INPUT_COLOR)
        ax2.set_ylabel(self.input_name, color=INPUT_COLOR)
        ax.set_ylabel(self.output_name, color=OUTPUT_COLOR)
        ax.set_xlabel("time")
        ax.set_title(self.model_id)
        plt.savefig(plot_filename)

    def makeClosedLoop(self, setpoint=1):
        """
        Creates a close loop system and a companion step response plot.
        <filename>_closedloop.png.
        """
        try:
            ts = self.system.simulateSISOClosedLoop(input_name=self.input_name, output_name=self.output_name,
                                                    setpoint=setpoint, times=self.times,
                                           kp=1, ki=None, kf=None, reference=1, is_steady_state=False)
        except Exception as error:
            print("Error in closed loop for %s: %s" % (self.filename, error))
            return
        # Make the plot
        plot_filename = self.model_id + CLOSED_LOOP_SFX + PNG_EXT
        outputs = ts[self.output_name].values
        inputs = np.repeat(reference, len(outputs))
        times = ts.index/1000
        _, ax = plt.subplots(1)
        ax.scatter(times, outputs, color=OUTPUT_COLOR)
        ax.plot(times, inputs, color=INPUT_COLOR)
        ax.set_ylabel(self.input_name, color=INPUT_COLOR)
        ax.set_ylabel(self.output_name, color=OUTPUT_COLOR)
        ax.set_xlabel("time")
        ax.set_title(self.model_id)
        plt.savefig(plot_filename)

    @classmethod
    def runBiomodels(cls, start=0, end=1e5, is_report=True, end_time=20):
        """
        Verifies files in BioModels
        """
        verifier = cls(filename, model, is_report=is_report, end_time=end_time)
        checkFunctions = [verifier._checkStaircase, verifier._checkClosedloop]
        iterator = iterateBiomodels(start=start, end=end, is_report=is_report, checkerFunctions=checkFunctions)
        for filename, model in iterator:
            verifier = cls(filename, model)
            verifier._makeSystem()
            verifier.makeStaircase()
            verifier.makeClosedLoop()


if __name__ == '__main__':
    SISOMaker.runBiomodels(start=1)
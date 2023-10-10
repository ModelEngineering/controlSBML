"""
SISOMaker creates a default SBMLSystem for an SBML model to verify that the SBMLSystem works for the SBML.
This consists of: (a) create an SBMLSystem; (b) create a staircase response,
(c) and construct a closed loop system and simulate it.

Usage:
    from controlSBML.siso_maker import SISOMaker
    SBMLMaker.run(model_reference)

    To run on all of BioModels:
        SISOMaker.runBiomodels(start=1)
"""

from controlSBML import constants as cn
from controlSBML.sbml_system import SBMLSystem
from controlSBML.staircase import Staircase
from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML.iterate_biomodels import iterateBiomodels
from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML.siso_transfer_function_builder import SISOTransferFunctionBuilder
from controlSBML import util

import matplotlib.pyplot as plt
import numpy as np
import os

INPUT_COLOR = "red"
OUTPUT_COLOR = "blue"
STAIRCASE_SFX = "_staircase"
CLOSED_LOOP_SFX = "_closed_loop"
PNG_EXT = ".png"


class SISOMaker(object):

    def __init__(self, model, model_id="model", input_name=None, output_name=None, end_time=20):
        """
        Args:
            model: model reference (URL, SBML, Antimony, XML)
            model_id: str (identifier of model)
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
        self.input_name = SBMLSystem.makeInputName(input_name, self.roadrunner)
        self.output_name = SBMLSystem.makeOutputName(output_name, self.roadrunner, input_names=[self.input_name])
        self.system, self.times = self._makeSystem()
        self.builder = SISOTransferFunctionBuilder(self.system)

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

    @staticmethod
    def _checkFile(filename, _, suffix=STAIRCASE_SFX):
        """
        Checks if a file exists for the analysis done.

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

    def makeStaircase(self, initial_value=0, final_value=10, is_show=False, **kwargs):
        """
        Creates a staircase response and a plot. Creates the file
        <filename>_staircase.png.

        Args:
            initial_value (int, optional): _description_. Defaults to 0.
            final_value (int, optional): _description_. Defaults to 10.
        """
        staircase = Staircase(initial_value=initial_value, final_value=final_value)
        try:
            ts = self.builder.makeStaircaseResponse(staircase=staircase)
        except Exception as error:
            print("Error in staricase for %s: %s" % (self.model_id, error))
            return
        # Make the plot
        is_plot, new_kwargs = util.setNoPlot(kwargs)
        self.builder.plotStaircaseResponse(ts, title=self.model_id, **new_kwargs)
        plot_filename = self.model_id + STAIRCASE_SFX + PNG_EXT
        plt.savefig(plot_filename)
        if is_plot:
            plt.show()

    def makeClosedLoop(self, setpoint=1, **kwargs):
        """
        Creates a close loop system and a companion step response plot.
        <filename>_closedloop.png.
        """
        try:
            ts = self.system.simulateSISOClosedLoop(input_name=self.input_name, output_name=self.output_name,
                                                    setpoint=setpoint,
                                           kp=1, ki=None, kf=None, is_steady_state=False)
        except Exception as error:
            print("Error in closed loop for %s: %s" % (self.model_id, error))
            return
        # Make the plot
        is_plot, new_kwargs = util.setNoPlot(kwargs)
        self.system.plotSISOClosedLoop(ts, setpoint=setpoint, figsize=(5,5), title=self.model_id, **new_kwargs)
        plot_filename = self.model_id + CLOSED_LOOP_SFX + PNG_EXT
        plt.savefig(plot_filename)
        if is_plot:
            plt.show()

    @classmethod
    def runModel(cls, model, **kwargs):
        """
        Constructs staircase and closed loop plots for a model.
            model__staircase.png, model_closedloop.png
        Args:
            model: model reference (URL, SBML, Antimony, XML)
            kwargs: keyword arguments for SISOMaker
        """
        maker = cls(model, **kwargs)
        maker._makeSystem()
        maker.makeStaircase()
        maker.makeClosedLoop()
    
    @classmethod
    def runBiomodels(cls, start=0, end=1e5, is_report=True, end_time=20):
        """
        Makes models for all models in BioModels
        """
        maker = cls(filename, model, is_report=is_report, end_time=end_time)
        checkFunctions = [maker._checkStaircase, maker._checkClosedloop]
        iterator = iterateBiomodels(start=start, end=end, is_report=is_report, checkerFunctions=checkFunctions)
        for filename, model in iterator:
            cls.runModel(model, model_id=filename)
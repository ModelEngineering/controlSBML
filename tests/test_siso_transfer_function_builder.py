import controlSBML.constants as cn
import controlSBML as ctl
import controlSBML.siso_transfer_function_builder as stb
from controlSBML.sbml_system import SBMLSystem
from controlSBML.staircase import Staircase
import controlSBML.util as util
import helpers

import copy
import matplotlib.pyplot as plt
import numpy as np
import os
import unittest
import shutil
import tellurium as te
import tempfile


IGNORE_TEST = True
IS_PLOT = True
PLOT_PATH = helpers.setupPlotting(__file__)
END_TIME = 5
DT = 0.01
POINTS_PER_TIME = int(1.0 / DT)
NUM_TIME = int(POINTS_PER_TIME*END_TIME) + 1
TIMES = [n*DT for n in range(0, NUM_TIME)]
LINEAR_MDL = """
// Illustrate Antimony File
model *linear()
J0:  -> S1; k1
J1: S1 -> S2; S1
J2: S2 -> ; S2

k1 = 5
S1 = 10
S2 = 0
end
"""
INPUT_NAME = "S1"
OUTPUT_NAME = "S2"
INITIAL_VALUE = 2
FINAL_VALUE = 15
STAIRCASE= Staircase(initial_value=INITIAL_VALUE, final_value=FINAL_VALUE, num_step=5)
SYSTEM = SBMLSystem(LINEAR_MDL, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME], is_fixed_input_species=True)
BUILDER = stb.SISOTransferFunctionBuilder(SYSTEM)
RESPONSE_TS, _ = BUILDER.makeStaircaseResponse(staircase=STAIRCASE, times=np.linspace(0, END_TIME, NUM_TIME))


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

    def testMakeParameters(self):
        if IGNORE_TEST:
            return
        def test(parameters, prefix, expected_count):
            names = [n for n in parameters.valuesdict().keys()
                  if n[0] == prefix]
            self.assertEqual(len(names), expected_count)
        #
        parameters = stb._makeParameters(3, 4)
        test(parameters, stb.NUMERATOR_PREFIX, 3)
        test(parameters, stb.DENOMINATOR_PREFIX, 4)
        self.assertTrue(len(parameters.valuesdict()) == 7)

    def testMakeTransferFunction(self):
        if IGNORE_TEST:
            return
        parameters = stb._makeParameters(2, 2)
        tf = stb._makeTransferFunction(parameters)
        self.assertTrue(tf.poles()[0] == -1)
        self.assertTrue(tf.dcgain() == 1)

    def testCalculateTransferFunctionResiduals(self):
        if IGNORE_TEST:
            return
        times = list(RESPONSE_TS.index)
        data_in = (times, RESPONSE_TS["S1_staircase"].values)
        data_out = RESPONSE_TS["S2"].values
        parameters = stb._makeParameters(3, 3)
        residuals = stb._calculateTransferFunctionResiduals(parameters, data_in,
              data_out)
        self.assertEqual(len(residuals), len(times))
        self.assertTrue("float" in str(residuals.dtype))


#############################
# Tests
#############################
class TestSBMLSystem(unittest.TestCase):

    def setUp(self):
        if IGNORE_TEST:
            return
        self.init()
        self.removeFiles()

    def tearDown(self):
        plt.close()
        self.removeFiles()

    def init(self, do_simulate_on_update=True):
        if IS_PLOT:
            cn.PLOT_DIR = cn.TEST_DIR
        else:
            cn.PLOT_DIR= tempfile.mkdtemp()
        self.system = copy.deepcopy(SYSTEM)
        self.initial_value = INITIAL_VALUE
        self.final_value = FINAL_VALUE
        self.builder = BUILDER
        self.response_ts = RESPONSE_TS

    def removeFiles(self):
        if IS_PLOT:
            return
        for ffile in os.listdir(cn.PLOT_DIR):
            if ("figure_" in ffile) and (".pdf") in ffile:
                path = os.path.join(cn.PLOT_DIR, ffile)
                if os.path.isfile(path) and IGNORE_TEST:
                    os.remove(path)
        # FIXME: This won't work in windows
        if IS_PLOT and ("var" in cn.PLOT_DIR):
            shutil.rmtree(cn.PLOT_DIR)
        #
        if os.path.isfile(PLOT_PATH):
            os.remove(PLOT_PATH)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.init()
        self.assertTrue(isinstance(self.builder, stb.SISOTransferFunctionBuilder ))

    def testMakeStaircaseResponse(self):
        if IGNORE_TEST:
            return
        self.init()
        diff = set(self.response_ts.columns) - set(["S2", "S1_staircase"])
        self.assertEqual(len(diff), 0)
        self.assertGreater(len(self.response_ts), 0)
        self.assertEqual(self.response_ts["S1_staircase"].values[0], self.initial_value)
        self.assertEqual(self.response_ts["S1_staircase"].values[-1], self.final_value)

    def testPlotStaircaseResponse(self):
        if IGNORE_TEST:
            return
        self.init()
        plot_result = self.builder.plotStaircaseResponse(self.response_ts, is_plot=IS_PLOT)
        self.assertTrue(isinstance(plot_result, util.PlotResult))

    def testFitTransferFunction(self):
        if IGNORE_TEST:
            return
        self.init()
        system = SBMLSystem(LINEAR_MDL, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME], is_fixed_input_species=True)
        builder = stb.SISOTransferFunctionBuilder(system)
        fitter_result = builder.fitTransferFunction(num_numerator=4, num_denominator=4,
              end_time=100)
        if IS_PLOT:
            builder.plotFitterResult(fitter_result, is_plot=IS_PLOT)
        self.assertTrue(isinstance(fitter_result.time_series, ctl.Timeseries))
        self.assertLess(fitter_result.rms_residuals, 0.2)

    def testFitTransferFunctionTimes(self):
        if IGNORE_TEST:
            return
        self.init()
        times = np.linspace(0, 100, 1000)
        system = SBMLSystem(LINEAR_MDL, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME], is_fixed_input_species=True)
        builder = stb.SISOTransferFunctionBuilder(system)
        fitter_result = builder.fitTransferFunction(num_numerator=4, num_denominator=4,
              end_time=100, fit_start_time=0, fit_end_time=40, times=times, staircase=STAIRCASE)
        if IS_PLOT:
            builder.plotFitterResult(fitter_result, is_plot=IS_PLOT)
        self.assertTrue(isinstance(fitter_result.time_series, ctl.Timeseries))
        self.assertLess(fitter_result.rms_residuals, 0.1)

    def testPlotFitTransferFunction(self):
        if IGNORE_TEST:
            return
        system = SBMLSystem(LINEAR_MDL, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME], is_fixed_input_species=True)
        builder = stb.SISOTransferFunctionBuilder(system)
        fitter_result = builder.fitTransferFunction(num_numerator=4, num_denominator=4,
              end_time=50)
        builder.plotFitterResult(fitter_result, is_plot=IS_PLOT, figsize=(5,5))

    def testFitTransferFunction2(self):
        if IGNORE_TEST:
            return
        system = SBMLSystem(cn.WOLF_URL,
              input_names=["na"], output_names=["s6"])
        builder = stb.SISOTransferFunctionBuilder(system)
        staircase = Staircase(initial_value=50, final_value=100)
        fitter_result = builder.fitTransferFunction(1, 2, staircase=staircase,
              end_time=5)
        builder.plotFitterResult(fitter_result, is_plot=IS_PLOT)
        self.assertTrue(isinstance(fitter_result.time_series, ctl.Timeseries))


if __name__ == '__main__':
  unittest.main()
import controlSBML as ctl
import controlSBML.constants as cn
import controlSBML.siso_transfer_function_builder as stb
from controlSBML.staircase import Staircase
import helpers

import control
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import unittest
import shutil
import tellurium as te
import tempfile


IGNORE_TEST = False
IS_PLOT = True
PLOT_PATH = helpers.setupPlotting(__file__)
END_TIME = 5
DT = 0.01
POINTS_PER_TIME = int(1.0 / DT)
NUM_TIME = int(POINTS_PER_TIME*END_TIME) + 1
TIMES = [n*DT for n in range(0, NUM_TIME)]

LINEAR_MDL = """
J0: $S0 -> S1; k1*S0
J1: S1 -> S2; S1
J2: S2 -> S3; S2

k1 = 10
S0 = 1
S1 = 10
S2 = 0
S3 = 0
"""
LINEAR_MDL = """
J0:  -> S1; k1
J1: S1 -> S2; S1
J2: S2 -> ; S2

k1 = 5
S1 = 10
S2 = 0
"""
rr = te.loada(LINEAR_MDL)
INPUT_NAME = "S1"
OUTPUT_NAME = "S2"
ctlsb = ctl.ControlSBML(LINEAR_MDL, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME])
LINEAR_SYS = ctlsb.makeNonlinearIOSystem("LINEAR_SYS")
builder = stb.SISOTransferFunctionBuilder(LINEAR_SYS)
plot_response = builder.plotStaircaseResponse(staircase=Staircase(final_value=10), is_plot=False)
LINEAR_TS = plot_response.time_series


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
        parameters = stb.makeParameters(3, 4)
        test(parameters, stb.NUMERATOR_PREFIX, 3)
        test(parameters, stb.DENOMINATOR_PREFIX, 4)
        self.assertTrue(len(parameters.valuesdict()) == 7)

    def testMakeTransferFunction(self):
        if IGNORE_TEST:
            return
        parameters = stb.makeParameters(2, 2)
        tf = stb.makeTransferFunction(parameters)
        self.assertTrue(tf.poles()[0] == -1)
        self.assertTrue(tf.dcgain() == 1)

    def testCalculateTransferFunctionResiduals(self):
        if IGNORE_TEST:
            return
        times = LINEAR_TS.times
        data_in = (times, LINEAR_TS["S1_staircase"].values)
        data_out = LINEAR_TS[OUTPUT_NAME].values
        parameters = stb.makeParameters(3, 3)
        residuals = stb._calculateTransferFunctionResiduals(parameters, data_in,
              data_out)
        self.assertEqual(len(residuals), len(times))
        self.assertTrue("float" in str(residuals.dtype))


#############################
# Tests
#############################
class TestNonlinearIOSystem(unittest.TestCase):

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
        self.ctlsb = ctl.ControlSBML(LINEAR_MDL,
              input_names=[INPUT_NAME], output_names=[OUTPUT_NAME])
        self.sys = ctl.NonlinearIOSystem("test_sys", self.ctlsb,
               do_simulate_on_update=do_simulate_on_update)
        self.builder = stb.SISOTransferFunctionBuilder(self.sys)

    def removeFiles(self):
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

    def testFitTransferFunction(self):
        if IGNORE_TEST:
            return
        self.init()
        fitter_result = self.builder.fitTransferFunction(4, 4,
              end_time=100)
        self.assertTrue(isinstance(fitter_result.time_series, ctl.Timeseries))
        if IS_PLOT:
            fitter_result.time_series["staircase"] = fitter_result.staircase.staircase_arr
            ctl.plotOneTS(fitter_result.time_series, writefig=PLOT_PATH)
    
    def testPlotFit(self):
        if IGNORE_TEST:
            return
        self.init()
        fitter_result = self.builder.fitTransferFunction(4, 4,
              end_time=50)
        self.builder.plotFit(fitter_result, is_plot=IS_PLOT)

    def testFitTransferFunction2(self):
        if IGNORE_TEST:
            return
        ctlsb = ctl.ControlSBML("https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000206.2?filename=BIOMD0000000206_url.xml",
              input_names=["s5"], output_names=["s6"])
        builder = ctlsb.makeSISOTransferFunctionBuilder()
        staircase = Staircase(initial_value=0, final_value=10)
        fitter_result = builder.fitTransferFunction(3, 3, staircase=staircase,
              end_time=5)
        builder.plotFit(fitter_result, is_plot=IS_PLOT)
        self.assertTrue(isinstance(fitter_result.time_series, ctl.Timeseries))


if __name__ == '__main__':
  unittest.main()

import controlSBML.constants as cn  # type: ignore
import controlSBML as ctl
import controlSBML.siso_transfer_function_builder as stb  # type: ignore
from controlSBML.sbml_system import SBMLSystem # type: ignore
from controlSBML.staircase import Staircase # type: ignore
import controlSBML.util as util # type: ignore
import helpers

import copy
import matplotlib.pyplot as plt
import numpy as np
import os
import unittest
import shutil
import tellurium as te # type: ignore
import tempfile


IGNORE_TEST = False
IS_PLOT = False
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
NEG_MODEL = """
// Negative relationship between S0 and S3
model *neg()
$S1 -> S2; k1*S1
S2 -> S3; k2*S2
S3 ->; k3*S3
S3 -> ; k4*S0
S1 = 10
$S0 = 10
S3 = 10
k1 = 1
k2 = 2
k3 = 2
k4 = 2
end
"""
INPUT_NAME = "S1"
OUTPUT_NAME = "S2"
INITIAL_VALUE = 2
FINAL_VALUE = 15
STAIRCASE= Staircase(initial_value=INITIAL_VALUE, final_value=FINAL_VALUE, num_step=5)
SYSTEM = SBMLSystem(LINEAR_MDL, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME], is_fixed_input_species=True)
if not IGNORE_TEST:
    BUILDER = stb.SISOTransferFunctionBuilder(SYSTEM)
    RESPONSE_TS, _ = BUILDER.makeStaircaseResponse(staircase=STAIRCASE, times=np.linspace(0, END_TIME, NUM_TIME))


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
        fitter_result = builder.plotTransferFunctionFit(is_plot=IS_PLOT)
        self.assertTrue(isinstance(fitter_result.time_series, ctl.Timeseries))
        self.assertLess(fitter_result.rms_residuals, 1)

    def testPlotFitTransferFunction(self):
        if IGNORE_TEST:
            return
        system = SBMLSystem(LINEAR_MDL, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME], is_fixed_input_species=True)
        builder = stb.SISOTransferFunctionBuilder(system)
        _ = builder.plotTransferFunctionFit(num_zero=3, num_pole=3,
              end_time=500, is_plot=IS_PLOT, figsize=(5,5))

    def testFitTransferFunction2(self):
        if IGNORE_TEST:
            return
        system = SBMLSystem(cn.WOLF_PATH,
              input_names=["na"], output_names=["s6"])
        builder = stb.SISOTransferFunctionBuilder(system)
        staircase = Staircase(initial_value=50, final_value=100)
        fitter_result = builder.plotTransferFunctionFit(1, 2, staircase=staircase,
              end_time=5, is_plot=IS_PLOT)
        self.assertTrue(isinstance(fitter_result.time_series, ctl.Timeseries))

    def testFitTransferFunctionDecrease(self):
        if IGNORE_TEST:
            return
        system = SBMLSystem(NEG_MODEL,
              input_names=["S0"], output_names=["S3"])
        builder = stb.SISOTransferFunctionBuilder(system)
        staircase = Staircase()
        fitter_result = builder.plotTransferFunctionFit(num_zero=1, num_pole=2, staircase=staircase,
                                                    fit_start_time=2, fit_end_time=10, is_plot=IS_PLOT)

    def testFitTransferFunctionBug(self):
        if IGNORE_TEST:
            return
        url = URL = "https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1304160000/2/BIOMD0000000449_url.xml"
        system = SBMLSystem(url,
              input_names=["IR"], output_names=["GLUT4"], is_fixed_input_species=True)
        builder = stb.SISOTransferFunctionBuilder(system, fitter_method=cn.FITTER_METHOD_POLY)
        staircase = Staircase(initial_value=10, final_value=25)
        _ = builder.plotTransferFunctionFit(num_zero=1, num_pole=3, staircase=staircase,
                                            times=np.linspace(0, 1000, 10000),
              fit_start_time=100, fit_end_time=1000, is_plot=IS_PLOT)


if __name__ == '__main__':
  unittest.main()
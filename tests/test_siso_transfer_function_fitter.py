import controlSBML.constants as cn  # type: ignore
from controlSBML.siso_transfer_function_fitter import SISOTransferFunctionFitter  # type: ignore
from controlSBML.siso_transfer_function_builder import SISOTransferFunctionBuilder  # type: ignore
from controlSBML.sbml_system import SBMLSystem # type: ignore
from controlSBML.staircase import Staircase # type: ignore
import controlSBML.util as util # type: ignore
import helpers

import control  # type: ignore
import matplotlib.pyplot as plt
import numpy as np
import unittest
import tellurium as te # type: ignore


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
INPUT_STAIRCASE = "S1_staircase"
OUTPUT_NAME = "S2"
INITIAL_VALUE = 2
FINAL_VALUE = 15
STAIRCASE= Staircase(initial_value=INITIAL_VALUE, final_value=FINAL_VALUE, num_step=5)
SYSTEM = SBMLSystem(LINEAR_MDL, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME], is_fixed_input_species=True)
BUILDER = SISOTransferFunctionBuilder(SYSTEM)
RESPONSE_TS, _ = BUILDER.makeStaircaseResponse(staircase=STAIRCASE, times=np.linspace(0, END_TIME, NUM_TIME))


#############################
# Tests
#############################
class TestSISOTransferFunctionFitter(unittest.TestCase):

    def setUp(self):
        self.fitter = SISOTransferFunctionFitter(RESPONSE_TS)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.fitter, SISOTransferFunctionFitter))
        self.assertEqual(self.fitter.input_name, INPUT_STAIRCASE)
        self.assertEqual(self.fitter.output_name, OUTPUT_NAME)

    def testPlot(self):
        #if IGNORE_TEST:
        #    return
        self.fitter.transfer_function = control.TransferFunction([1], [1, 1])
        self.fitter.plot(is_plot=IS_PLOT)

    def testUniformFromLogspace(self):
        if IGNORE_TEST:
            return
        def test(min_val, max_val, num):
            values = self.fitter._uniformFromLogspace(min_val, max_val, num)
            self.assertEqual(len(values), num)
            self.assertTrue(all([min_val <= v <= max_val for v in values]))
        #
        test(1, 10, 50)
        test(1, 10, 10)
        test(1, 10, 100)


if __name__ == '__main__':
  unittest.main()
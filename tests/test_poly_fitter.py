import controlSBML.constants as cn  # type: ignore
from controlSBML.poly_fitter import PolyFitter  # type: ignore
from controlSBML.siso_transfer_function_builder import SISOTransferFunctionBuilder  # type: ignore
from controlSBML.sbml_system import SBMLSystem # type: ignore
from controlSBML.staircase import Staircase # type: ignore
from controlSBML.timeseries import Timeseries # type: ignore
from controlSBML.control_sbml import ControlSBML # type: ignore

import matplotlib.pyplot as plt  # type: ignore
import control  # type: ignore
import numpy as np
import unittest


IGNORE_TEST = False
IS_PLOT = False
END_TIME = 10
TIMES = [0.01*n for n in range(100*END_TIME)]
NUM_TIME = len(TIMES)
LINEAR_MDL = """
// Simple Antimony File
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
NUM_ZERO = 3
NUM_POLE = 4
STAIRCASE= Staircase(initial_value=INITIAL_VALUE, final_value=FINAL_VALUE, num_step=5)
SYSTEM = SBMLSystem(LINEAR_MDL, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME], is_fixed_input_species=True)
BUILDER = SISOTransferFunctionBuilder(SYSTEM)
RESPONSE_TS, _ = BUILDER.makeStaircaseResponse(staircase=STAIRCASE, times=TIMES)


#############################
# Tests
#############################
class TestPolyFitter(unittest.TestCase):

    def setUp(self):
        self.num_zero = NUM_ZERO
        self.num_pole = NUM_POLE
        self.ts = Timeseries(RESPONSE_TS.copy())
        self.fitter = PolyFitter(self.ts, num_zero=self.num_zero, num_pole=self.num_pole)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.fitter, PolyFitter))

    def testMakeParameters(self):
        if IGNORE_TEST:
            return
        def test(num_zero, num_pole, is_zero_None=False, is_pole_None=False):
            def addOne(num):
                if num is None:
                    return None
                else:
                    return num + 1
            #
            fitter = PolyFitter(self.ts, num_zero=num_zero, num_pole=num_pole)
            if is_zero_None:
                num_zero = None
            if is_pole_None:
                num_pole = None
            parameters = fitter._makeParameters(addOne(num_zero), addOne(num_pole))
            num_n = len([n for n in parameters.valuesdict().keys() if n[0] == 'n'])
            num_d = len([n for n in parameters.valuesdict().keys() if n[0] == 'd'])
            if (num_zero is not None) and (num_pole is not None):
                self.assertEqual(len(parameters.valuesdict()), num_n + num_d)
                self.assertEqual(num_n, addOne(num_zero))
                self.assertEqual(num_d, addOne(num_pole))
        #
        test(2, 3, is_zero_None=False)
        test(2, 3, is_zero_None=True)
        test(2, 3, is_pole_None=True)
        test(2, 3, is_zero_None=True, is_pole_None=True)
        test(3, 4)
        test(0, 1)

    def testZeroGreaterThanPole(self):
        if IGNORE_TEST:
            return
        with self.assertRaises(ValueError):
            PolyFitter(self.ts, num_zero=4, num_pole=3)

    def testMakeTransferFunction(self):
        if IGNORE_TEST:
            return
        def test(num_zero, num_pole):
            fitter = PolyFitter(self.ts, num_zero=num_zero, num_pole=num_pole)
            parameters = fitter._makeParameters()
            tf = fitter._makeTransferFunction(parameters)
            self.assertTrue(isinstance(tf, control.TransferFunction))
        #
        test(3, 4)
        test(0, 1)
        test(1, 1)

    def testCalculateResiduals(self):
        if IGNORE_TEST:
            return
        parameters = self.fitter._makeParameters()
        rmse = np.sqrt(np.sum(self.fitter._calculateResiduals(parameters)**2))
        self.assertTrue(isinstance(rmse, float))

    def testFit(self):
        if IGNORE_TEST:
            return
        for _ in range(5):
            self.fitter.fit()
            self.fitter.plot(is_plot=IS_PLOT)
            self.assertTrue(isinstance(self.fitter.transfer_function, control.TransferFunction))
            mse = self.fitter._calculateNormalizedMSE(self.fitter.transfer_function)
            if mse < 2:
                self.assertLess(mse, 2)
                break
        self.assertLess(mse, 2)

    def testFitTransferFunctionBug(self):
        if IGNORE_TEST:
            return
        URL = "https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1304160000/2/BIOMD0000000449_url.xml"
        ctlsb = ControlSBML(URL,
              input_name="IR", output_name="GLUT4", is_fixed_input_species=True, is_plot=False)
        result = ctlsb.plotStaircaseResponse(initial_value=10, final_value=25, times=np.linspace(0, 1.5, 150),
                                            is_plot=False)
        plt.close()
        fitter = PolyFitter(result.timeseries, num_pole=2, num_zero=0)
        fitter.fit()
        fitter.plot(title="FitTransferFunctionBug", is_plot=IS_PLOT)


if __name__ == '__main__':
  unittest.main()
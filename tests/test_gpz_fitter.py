import controlSBML.constants as cn  # type: ignore
from controlSBML.gpz_fitter import GPZFitter  # type: ignore
from controlSBML.siso_transfer_function_builder import SISOTransferFunctionBuilder  # type: ignore
from controlSBML.sbml_system import SBMLSystem # type: ignore
from controlSBML.staircase import Staircase # type: ignore
from controlSBML.timeseries import Timeseries # type: ignore
from controlSBML.control_sbml import ControlSBML # type: ignore

import matplotlib.pyplot as plt
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
STAIRCASE= Staircase(initial_value=INITIAL_VALUE, final_value=FINAL_VALUE, num_step=5)
SYSTEM = SBMLSystem(LINEAR_MDL, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME], is_fixed_input_species=True)
BUILDER = SISOTransferFunctionBuilder(SYSTEM)
RESPONSE_TS, _ = BUILDER.makeStaircaseResponse(staircase=STAIRCASE, times=TIMES)


#############################
# Tests
#############################
class TestGPZFitter(unittest.TestCase):

    def setUp(self):
        self.num_zero = 3
        self.num_pole = 4
        self.ts = Timeseries(RESPONSE_TS.copy())
        self.fitter = GPZFitter(self.ts, num_zero=self.num_zero, num_pole=self.num_pole)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.fitter, GPZFitter))

    def testMakeParameters(self):
        if IGNORE_TEST:
            return
        def test(num_zero, num_pole, is_zero_None=False, is_pole_None=False):
            fitter = GPZFitter(self.ts, num_zero=num_zero, num_pole=num_pole)
            if is_zero_None:
                zeros = None
            else:
                zeros = range(num_zero)
            if is_pole_None:
                poles = None
            else:
                poles = range(num_pole)
            parameters = fitter._makeParameters(zeros, poles)
            num_z = len([n for n in parameters.valuesdict().keys() if n[0] == 'z'])
            num_p = len([n for n in parameters.valuesdict().keys() if n[0] == 'p'])
            self.assertEqual(len(parameters.valuesdict()), num_z + num_p)
            self.assertEqual(num_z, num_zero)
            self.assertEqual(num_p, num_pole)
        #
        test(2, 3, is_zero_None=True)
        test(2, 3, is_pole_None=True)
        test(2, 3, is_zero_None=True, is_pole_None=True)
        test(3, 4)
        test(0, 1)

    def testMakeTransferFunction(self):
        if IGNORE_TEST:
            return
        def test(num_zero, num_pole, dcgain):
            fitter = GPZFitter(self.ts, num_zero=num_zero, num_pole=num_pole)
            zeros = -np.array(range(num_zero))
            poles = -1*np.array(range(num_pole))
            parameters = fitter._makeParameters(zeros, poles)
            fitter.dcgain = dcgain
            tf = fitter._makeTransferFunction(parameters)
            self.assertTrue(isinstance(tf, control.TransferFunction))
            tf = fitter._makeTransferFunction(parameters)
            for pole in tf.poles():
                self.assertTrue(any([np.isclose(pole, v)]) for n, v in parameters.valuesdict().values() if n[0] == 'p')
            for zero in tf.zeros():
                self.assertTrue(any([np.isclose(zero, v)]) for n, v in parameters.valuesdict().values() if n[0] == 'z')
            self.assertTrue(np.isclose(tf.dcgain(), dcgain))
        #
        test(3, 4, -15)
        test(0, 1, 10)

    def testCalculateTransferFunctionResiduals(self):
        if IGNORE_TEST:
            return
        dcgain = 1
        self.fitter.dcgain = dcgain
        parameters = self.fitter._makeParameters()
        tf = self.fitter._makeTransferFunction(parameters)
        rmse = self.fitter._calculateNormalizedMSE(tf)
        self.assertTrue(isinstance(rmse, float))

    def testFitDCGain(self):
        if IGNORE_TEST:
            return
        self.fitter._fitDCGain()
        self.plot()
        self.assertTrue(isinstance(self.fitter.dcgain, float))

    def testFitPoles(self):
        if IGNORE_TEST:
            return
        ts, _ = BUILDER.makeStaircaseResponse(staircase=STAIRCASE, times=np.linspace(0, 5, 50))
        fitter = GPZFitter(ts, num_zero=0, num_pole=2)
        fitter._fitDCGain()
        fitter._fitPoles()
        self.plot(fitter)
        self.assertTrue(isinstance(fitter.poles, np.ndarray))

    def testFitZeros(self):
        if IGNORE_TEST:
            return
        ts, _ = BUILDER.makeStaircaseResponse(staircase=STAIRCASE, times=np.linspace(0, 5, 50))
        fitter = GPZFitter(ts, num_zero=2, num_pole=2)
        fitter._fitDCGain()
        fitter._fitPoles()
        fitter._fitZeros()
        self.plot(fitter)
        self.assertTrue(isinstance(fitter.zeros, np.ndarray))

    def plot(self, in_fitter=None, title=""):
        if not IS_PLOT:
            return
        plt.close()
        if in_fitter is None:
            fitter = self.fitter.copy()
        else:
            fitter = in_fitter.copy()
        if fitter.dcgain is None:
            raise ValueError("dcgain is None")
        if fitter.poles is None:
            fitter.pole = []
            fitter.zeros = []
        if fitter.zeros is None:
            fitter.zeros = []
        #
        if fitter.transfer_function is None:
            fitter.fit()
        _, predictions = fitter._simulateTransferFunction()
        plt.plot(fitter.times, fitter.in_arr, linestyle='--')
        plt.plot (fitter.times, fitter.out_arr)
        plt.plot (fitter.times, predictions)
        plt.legend(["input", "output", "prediction"])
        plt.title(title)
        plt.title(str(fitter.transfer_function), x=0.7, y=0.1)
        plt.show()

    def testFitTransferFunctionBug(self):
        #if IGNORE_TEST:
        #    return
        URL = "https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1304160000/2/BIOMD0000000449_url.xml"
        ctlsb = ControlSBML(URL,
              input_name="IR", output_name="GLUT4", is_fixed_input_species=True, is_plot=False)
        ts, _ = ctlsb.plotStaircaseResponse(initial_value=10, final_value=25, times=np.linspace(0, 1.5, 150),
                                            is_plot=False)
        plt.close()
        fitter = GPZFitter(ts, num_itr=5)
        fitter.fit()
        fitter.plot(title="FitTransferFunctionBug", is_plot=IS_PLOT)
        self.assertLess(len(fitter.poles), 3)
        self.assertEqual(len(fitter.zeros), 0)

    def testGetNext(self):
        if IGNORE_TEST:
            return
        def count(num_increment, num_entry):
            results = list(self.fitter._getNext(1, 10, num_increment, num_entry))
            self.assertTrue(any([len(r) == 0 for r in results]))
            trues = [len(r) <= num_entry for r in results]
            self.assertTrue(all(trues))
            return len(results)
        #
        self.assertEqual(count(3, 2), 10)
        self.assertEqual(count(3, 1), 4)
        self.assertEqual(count(2, 1), 3)


if __name__ == '__main__':
  unittest.main()
import controlSBML.constants as cn  # type: ignore
from controlSBML.siso_transfer_function_fitter import SISOTransferFunctionFitter  # type: ignore
from controlSBML.siso_transfer_function_builder import SISOTransferFunctionBuilder  # type: ignore
from controlSBML.sbml_system import SBMLSystem # type: ignore
from controlSBML.staircase import Staircase # type: ignore
from controlSBML.timeseries import Timeseries # type: ignore

import matplotlib.pyplot as plt
import control
import numpy as np
import unittest


IGNORE_TEST = True
IS_PLOT = True
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
class TestSISOTransferFunctionFitter(unittest.TestCase):

    def setUp(self):
        self.num_zero = 3
        self.num_pole = 4
        self.ts = Timeseries(RESPONSE_TS.copy())
        self.fitter = SISOTransferFunctionFitter(self.ts, num_zero=self.num_zero, num_pole=self.num_pole)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.fitter, SISOTransferFunctionFitter))
        self.assertEqual(len(self.fitter.min_pole_values), self.num_pole)

    def testUniformFromLogspace(self):
        if IGNORE_TEST:
            return
        def test(min_val, max_val, num):
            values = self.fitter.uniformFromLogspace(min_val, max_val, num)
            self.assertEqual(len(values), num)
            self.assertTrue(all([min_val <= v <= max_val for v in values]))
        #
        test(1, 10, 50)
        test(1, 10, 10)
        test(1, 10, 100)

    def testMakeParameters(self):
        if IGNORE_TEST:
            return
        def test(num_zero, num_pole, is_zero_None=False, is_pole_None=False):
            fitter = SISOTransferFunctionFitter(self.ts, num_zero=num_zero, num_pole=num_pole)
            if is_zero_None:
                zeros = None
            else:
                zeros = range(num_zero)
            if is_pole_None:
                poles = None
            else:
                poles = range(num_pole)
            parameters = fitter.makeParameters(zeros, poles)
            num_z = len([n for n in parameters.valuesdict().keys() if n[0] == 'z'])
            num_p = len([n for n in parameters.valuesdict().keys() if n[0] == 'p'])
            self.assertEqual(len(parameters.valuesdict()), num_zero + num_pole)
            self.assertEqual(num_z, num_zero)
            self.assertEqual(num_p, num_pole)
        #
        test(2, 1, is_pole_None=True)
        test(2, 1, is_zero_None=True)
        test(2, 1, is_zero_None=True, is_pole_None=True)
        test(3, 4)
        test(0, 1)
        test(1, 0)

    def testMakeTransferFunction(self):
        if IGNORE_TEST:
            return
        def test(num_zero, num_pole, dcgain):
            fitter = SISOTransferFunctionFitter(self.ts, num_zero=num_zero, num_pole=num_pole)
            zeros = -np.array(range(num_zero))
            poles = -1*np.array(range(num_pole))
            parameters = fitter.makeParameters(zeros, poles)
            fitter.dcgain = dcgain
            tf = fitter.makeTransferFunction(parameters)
            self.assertTrue(isinstance(tf, control.TransferFunction))
            tf = fitter.makeTransferFunction(parameters)
            for pole in tf.poles():
                self.assertTrue(any([np.isclose(pole, v)]) for n, v in parameters.valuesdict().values() if n[0] == 'p')
            for zero in tf.zeros():
                self.assertTrue(any([np.isclose(zero, v)]) for n, v in parameters.valuesdict().values() if n[0] == 'z')
            self.assertTrue(np.isclose(tf.dcgain(), dcgain))
        #
        test(3, 4, -15)
        test(0, 1, 10)
        test(1, 0, 10)

    def testCalculateTransferFunctionResiduals(self):
        if IGNORE_TEST:
            return
        dcgain = 1
        self.fitter.dcgain = dcgain
        parameters = self.fitter.makeParameters()
        tf = self.fitter.makeTransferFunction(parameters)
        rmse = self.fitter._calculateResiduals(tf)
        self.assertTrue(isinstance(rmse, float))

    def testFitDCGain(self):
        if IGNORE_TEST:
            return
        self.fitter.fitDCGain()
        if IS_PLOT:
            ts, _ = BUILDER.makeStaircaseResponse(staircase=STAIRCASE, times=np.linspace(0, 100, 1000))
            fitter = SISOTransferFunctionFitter(ts, num_zero=0, num_pole=0)
            fitter.fitDCGain()
            parameters = fitter.makeParameters([], [])
            tf = fitter.makeTransferFunction(parameters)
            prediction_times, predictions = fitter.simulateTransferFunction(tf)
            plt.plot(fitter.times, fitter.in_arr, linestyle='--')
            plt.plot(fitter.times, fitter.out_arr)
            plt.plot(prediction_times, predictions)
            plt.title("gain: {dcgain}".format(dcgain=self.fitter.dcgain))
            plt.legend(["input", "output", "prediction"])
            plt.show()
        self.assertTrue(isinstance(self.fitter.dcgain, float))

    def testFitPoles(self):
        if IGNORE_TEST:
            return
        ts, _ = BUILDER.makeStaircaseResponse(staircase=STAIRCASE, times=np.linspace(0, 5, 50))
        fitter = SISOTransferFunctionFitter(ts, num_zero=0, num_pole=2)
        fitter.fitDCGain()
        fitter.fitPoles()
        if IS_PLOT:
            parameters = fitter.makeParameters([], fitter.poles)
            tf = fitter.makeTransferFunction(parameters)
            times, predictions = fitter.simulateTransferFunction(tf)
            import pdb; pdb.set_trace()
            plt.plot(fitter.times, fitter.in_arr, linestyle='--')
            plt.plot (fitter.times, fitter.out_arr)
            plt.plot (times, predictions)
            plt.legend(["input", "output", "prediction"])
            plt.title("Fit dcgain, poles")
            plt.show()
        self.assertTrue(isinstance(self.fitter.poles, np.ndarray))

    def testFitPoles(self):
        #if IGNORE_TEST:
        #    return
        ts, _ = BUILDER.makeStaircaseResponse(staircase=STAIRCASE, times=np.linspace(0, 5, 50))
        fitter = SISOTransferFunctionFitter(ts, num_zero=2, num_pole=2, min_zero_value=-100)
        fitter.fitDCGain()
        fitter.fitPoles()
        fitter.fitZeros()
        if IS_PLOT:
            parameters = fitter.makeParameters(fitter.zeros, fitter.poles)
            tf = fitter.makeTransferFunction(parameters)
            _, predictions = fitter.simulateTransferFunction(tf)
            import pdb; pdb.set_trace()
            plt.plot(fitter.times, fitter.in_arr, linestyle='--')
            plt.plot (fitter.times, fitter.out_arr)
            plt.plot (fitter.times, predictions)
            plt.legend(["input", "output", "prediction"])
            plt.title("Fit dcgain, poles and zeros")
            plt.show()
        self.assertTrue(isinstance(self.fitter.poles, np.ndarray))


if __name__ == '__main__':
  unittest.main()
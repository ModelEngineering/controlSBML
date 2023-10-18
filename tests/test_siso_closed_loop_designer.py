import controlSBML.siso_closed_loop_designer as cld
from controlSBML.siso_transfer_function_builder import SISOTransferFunctionBuilder
import controlSBML as ctl
import controlSBML.util as util
import helpers

import copy
import control
import matplotlib.pyplot as plt
import numpy as np
import sympy
import unittest

IGNORE_TEST = False
IS_PLOT = False
helpers.setupPlotting(__file__)
MODEL = """
model *main_model()
S0 -> S1; k0*S0
S1 -> S2; k1*S1

k0 = 1
k1 = 1
S0 = 0
S1 = 0
S2 = 0
end
"""
MODEL2 = """
model *main2_model()
S0 -> S1; k0*S0
S1 -> S2; k1*S1
S2 -> ; k2*S2

k0 = 1
k1 = 1
k2 = 2
S0 = 10
S1 = 0
S2 = 0
end
"""
LINEAR_MDL = """
model *linear_model()
species S3

 -> S1; k0
S1 -> S2; k1*S1
S2 -> S3; k2*S2

S1 = 0
S2 = 0
S3 = 0
k0 = 0
k1 = 1
k2 = 2
end
"""
# Construct a transfer function for the model. This is a linear model, and so it should be accurate.
INPUT_NAME = "S0"
OUTPUT_NAME = "S2"
SYSTEM = ctl.SBMLSystem(MODEL, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME])
builder = SISOTransferFunctionBuilder(SYSTEM, input_name=INPUT_NAME, output_name=OUTPUT_NAME)
staircase = ctl.Staircase(final_value=15, num_step=5)
fitter_result = builder.fitTransferFunction(num_numerator=2, num_denominator=3, staircase=staircase)
if False:
    builder.plotFitTransferFunction(fitter_result, figsize=(5,5))
TRANSFER_FUNCTION = fitter_result.transfer_function
PARAMETER_DCT = {p: n+1 for n, p in enumerate(cld.PARAM_NAMES)}
TIMES = np.linspace(0, 20, 200)


#############################
# Tests
#############################
class TestSISOClosedLoopDesigner(unittest.TestCase):

    def setUp(self):
        self.sys_tf = copy.deepcopy(TRANSFER_FUNCTION)
        self.system = copy.deepcopy(SYSTEM)
        self.designer = cld.SISOClosedLoopDesigner(self.system, self.sys_tf, times=TIMES)

    def testGetSet(self):
        if IGNORE_TEST:
            return
        self.designer.set(**PARAMETER_DCT)
        for name in cld.PARAM_NAMES:
            self.assertEqual(getattr(self.designer, name), PARAMETER_DCT[name])
        #
        dct = self.designer.get()
        for name in cld.PARAM_NAMES:
            self.assertEqual(dct[name], PARAMETER_DCT[name])

    def testCalculateClosedLoopTf(self):
        if IGNORE_TEST:
            return
        sys_tf = control.tf([1], [1, 2])
        closed_loop_tf_kp = cld._calculateClosedLoopTf(sys_tf=sys_tf, kp=3)
        closed_loop_tf_ki = cld._calculateClosedLoopTf(sys_tf=sys_tf, ki=3)
        _, ys_kp = control.step_response(closed_loop_tf_kp, TIMES)
        _, ys_ki = control.step_response(closed_loop_tf_ki, TIMES)
        self.assertTrue(ys_kp[-1] < ys_ki[-1])
        self.assertTrue(np.isclose(ys_ki[-1], 1, atol=0.01))

    def test_closed_loop_tf(self):
        if IGNORE_TEST:
            return
        sys_tf = control.tf([1], [1, 1])
        designer = cld.SISOClosedLoopDesigner(sys_tf)
        with self.assertRaises(ValueError):
            _ = designer.closed_loop_tf()
        #
        designer.kf = 4
        closed_loop_tf = designer.closed_loop_tf
        numr = np.array(closed_loop_tf.num[0][0])
        self.assertTrue(np.allclose(numr, [20000, 80000, 0]))
        denr = np.array(closed_loop_tf.den[0][0])
        self.assertTrue(np.allclose(denr, [1.00000e+00, 1.00050e+04, 1.30004e+05, 4.00000e+04]))

    def testDesign1(self):
        if IGNORE_TEST:
            return
        def checkParams(names):
            for name in names:
                self.assertIsNotNone(getattr(designer, name))
            other_names = set(cld.PARAM_NAMES) - set(names)
            for name in other_names:
                self.assertIsNone(getattr(designer, name))
        designer = cld.SISOClosedLoopDesigner(SYSTEM, self.sys_tf)
        designer.design(kp=True)
        checkParams(["kp"])

    def testDesign2(self):
        if IGNORE_TEST:
            return
        numr = np.array([0.38613466, 0.16315592])
        denr = np.array([2.04313198, 1.77120743, 0.49351805])
        sys_tf = control.tf(numr, denr)
        designer1 = cld.SISOClosedLoopDesigner(SYSTEM, sys_tf)

    def testSimulate(self):
        if IGNORE_TEST:
            return
        def calcDiff(arr):
            return np.abs(arr[-1] - 1)
        #
        sys_tf = control.tf([1], [1, 1])
        designer = cld.SISOClosedLoopDesigner(SYSTEM, sys_tf)
        designer.set(kp=20)
        _, prediction1s = designer.simulate()
        designer.set(kp=20, ki=50)
        _, prediction2s = designer.simulate()
        self.assertLess(calcDiff(prediction2s), calcDiff(prediction1s))

    def testPlot(self):
        if IGNORE_TEST:
            return
        self.designer.set(**PARAMETER_DCT)
        self.designer.plot(is_plot=IS_PLOT, markers=["", ""])
        self.designer.set(kp=10, ki=5)
        self.designer.plot(is_plot=IS_PLOT)

    def makeDesigner(self, end_time=200):
        system = ctl.SBMLSystem(MODEL2, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME], is_fixed_input_species=True)
        times = np.linspace(0, end_time, 10*end_time)
        builder = SISOTransferFunctionBuilder(system, input_name=INPUT_NAME, output_name=OUTPUT_NAME)
        staircase = ctl.Staircase(final_value=15, num_step=5)
        fitter_result = builder.fitTransferFunction(num_numerator=2, num_denominator=3, staircase=staircase)
        if False:
            builder.plotFitTransferFunction(fitter_result, figsize=(5,5))
        designer = cld.SISOClosedLoopDesigner(system, fitter_result.transfer_function, times=times)
        return designer

    def testPlot2(self):
        if IGNORE_TEST:
            return
        designer = self.makeDesigner()
        designer.design(kp=True, ki=True)
        designer.plot(is_plot=IS_PLOT, markers=["", ""])
        # Check for stability
        self.assertTrue(util.isTransferFunctionStable(designer.closed_loop_tf))
    
    def testEvaluate2(self):
        if IGNORE_TEST:
            return
        def test(kp, ki):
            designer = self.makeDesigner(end_time=10)
            designer.design(kp=kp, ki=ki)
            designer.evaluate(is_plot=IS_PLOT)
            return designer.residual_rmse
        #
        rmse1 = test(kp=True, ki=False)
        rmse2 = test(kp=True, ki=True)
        self.assertLess(rmse2, rmse1)

    def test_closed_loop_tf(self):
        # Checks that the closed loop transfer function is calculated correctly
        if IGNORE_TEST:
            return
        # Setup
        s, kp, ki = sympy.symbols("s kp ki")
        # System transfer function
        sys_tf = control.tf([1], [1, 1])
        systf = 1/(s + 1)
        # Symbolic calculation of transfer function
        ctltf = kp + ki/s
        looptf = sympy.simplify(systf*ctltf)
        cltf = sympy.simplify(looptf/(1 + looptf))
        #
        designer = cld.SISOClosedLoopDesigner(self.system, sys_tf, step_size=5)
        designer.set(kp=2, ki=3)
        closed_loop_tf = designer.closed_loop_tf
        func1 = lambda x: float(closed_loop_tf(x))
        cltf_nums = cltf.subs({kp: 2, ki: 3})
        func2 = lambda x: sympy.N(cltf_nums.subs({s: x}))
        result = ctl.util.compareSingleArgumentFunctions(func1, func2, 0, 100)
        self.assertTrue(result)

    # FIXME: Bad fit
    def testBug1(self):
        if IGNORE_TEST:
            return
        # Setup
        system = ctl.SBMLSystem(LINEAR_MDL, input_names=["k0"], output_names=["S3"], is_fixed_input_species=True,
                                is_steady_state=False)
        linear_bldr = SISOTransferFunctionBuilder(system)
        linear_staircase = ctl.Staircase(initial_value=0, final_value=10, num_step=5)
        fitter_result = linear_bldr.fitTransferFunction(num_numerator=2, num_denominator=3, 
                                                    staircase=linear_staircase, fit_start_time=20,
                                                start_time=0, end_time=200)
        linear_bldr.plotFitTransferFunction(fitter_result, figsize=(5,5), is_plot=IS_PLOT)
        linear_tf = fitter_result.transfer_function
        numr = linear_tf.num[0][0]
        denr = linear_tf.den[0][0]
        #
        designer = cld.SISOClosedLoopDesigner(system, linear_tf, step_size=5)
        designer.set(kp=2, ki=3)
        if IS_PLOT:
            designer.plot(figsize=(5,5))  # Have options for a period



#############################
class TestHistory(unittest.TestCase):

    def setUp(self):
        self.sys_tf = TRANSFER_FUNCTION
        self.designer = cld.SISOClosedLoopDesigner(SYSTEM, self.sys_tf)
        self.history = cld._History(self.designer)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertEqual(len(self.history), 0)
        diff = set(cld.PARAM_NAMES) - set(self.history._dct.keys())
        self.assertEqual(len(diff), 0)

    def addItems(self, count):
        for _ in range(count):
            self.history.add()
    
    def testAdd(self):
        if IGNORE_TEST:
            return
        self.addItems(1)
        self.assertEqual(len(self.history), 1)
        self.addItems(1)
        self.assertEqual(len(self.history), 2)

    def testReport(self):
        if IGNORE_TEST:
            return
        COUNT = 3
        self.addItems(COUNT)
        df = self.history.report()
        self.assertEqual(len(df), COUNT)
    
    def testGet(self):
        if IGNORE_TEST:
            return
        COUNT = 1
        self.addItems(COUNT)
        designer = self.history.get(0)
        self.assertTrue(isinstance(designer, cld.SISOClosedLoopDesigner))

    def testClear(self):
        if IGNORE_TEST:
            return
        self.addItems(3)
        self.history.clear()
        self.assertEqual(len(self.history), 0)


if __name__ == '__main__':
    unittest.main()

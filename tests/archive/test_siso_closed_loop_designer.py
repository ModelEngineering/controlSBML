import controlSBML.siso_closed_loop_designer as scld
import controlSBML as ctl
import helpers

import control
import numpy as np
import sympy
import unittest

IGNORE_TEST = False
IS_PLOT = False
helpers.setupPlotting(__file__)
MODEL = """
S0 -> S1; k0*S0
S1 -> S2; k1*S1

k0 = 1
k1 = 1
S0 = 0
S1 = 0
S2 = 0
"""
LINEAR_MDL = """
species S3

 -> S1; k0
S1 -> S2; k1*S1
S2 -> S3; k2*S2
S3 -> ; k3*S3

S1 = 0
S2 = 0
S3 = 0
k0 = 0
k1 = 1
k2 = 2
k3 = 3
"""
# Construct a transfer function for the model. This is a linear model, and so it should be accurate.
CTLSB = ctl.ControlSBML(MODEL, input_names=["S0"], output_names=["S2"])
builder = CTLSB.makeSISOTransferFunctionBuilder(is_fixed_input_species=True)
staircase = ctl.Staircase(final_value=15, num_step=5)
fitter_result = builder.fitTransferFunction(num_numerator=2, num_denominator=3, staircase=staircase)
TRANSFER_FUNCTION = fitter_result.transfer_function
PARAMETER_DCT = {p: n+1 for n, p in enumerate(scld.PARAM_NAMES)}
TIMES = np.linspace(0, 20, 200)


#############################
# Tests
#############################
class TestSISOClosedLoopDesigner(unittest.TestCase):

    def setUp(self):
        self.sys_tf = TRANSFER_FUNCTION
        self.designer = scld.SISOClosedLoopDesigner(self.sys_tf)
        self.ctlsb = CTLSB

    def testGetSet(self):
        if IGNORE_TEST:
            return
        self.designer.set(**PARAMETER_DCT)
        for name in ["kp", "ki", "kd", "kf"]:
            self.assertEqual(getattr(self.designer, name), PARAMETER_DCT[name])
        #
        dct = self.designer.get()
        for name in ["kp", "ki", "kd", "kf"]:
            self.assertEqual(dct[name], PARAMETER_DCT[name])

    def testCalculateClosedLoopTf(self):
        if IGNORE_TEST:
            return
        sys_tf = control.tf([1], [1, 2])
        closed_loop_tf_kp = scld._calculateClosedLoopTf(open_loop_transfer_function=sys_tf, kp=3)
        closed_loop_tf_ki = scld._calculateClosedLoopTf(open_loop_transfer_function=sys_tf, ki=3)
        _, ys_kp = control.step_response(closed_loop_tf_kp, TIMES)
        _, ys_ki = control.step_response(closed_loop_tf_ki, TIMES)
        self.assertTrue(ys_kp[-1] < ys_ki[-1])
        self.assertTrue(np.isclose(ys_ki[-1], 1, atol=0.01))
        #
        # Check degree of the closed loop function
        closed_loop_tf_kd = scld._calculateClosedLoopTf(open_loop_transfer_function=sys_tf, kd=3)
        self.assertGreater(len(closed_loop_tf_kd.den[0][0]), len(closed_loop_tf_kp.den[0][0]))

    def test_closed_loop_tf(self):
        if IGNORE_TEST:
            return
        sys_tf = control.tf([1], [1, 1])
        designer = scld.SISOClosedLoopDesigner(sys_tf)
        with self.assertRaises(ValueError):
            _ = designer.closed_loop_tf()
        #
        designer.kd = 2
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
            other_names = set(scld.PARAM_NAMES) - set(names)
            for name in other_names:
                self.assertIsNone(getattr(designer, name))
        sys_tf = control.tf([1], [1, 1])
        designer = scld.SISOClosedLoopDesigner(sys_tf)
        designer.design(kp=1.0)
        checkParams(["kp"])
        #
        designer.design(kd=True, kf=True)
        checkParams(["kd", "kf"])

    def testDesign2(self):
        if IGNORE_TEST:
            return
        numr = np.array([0.38613466, 0.16315592])
        denr = np.array([2.04313198, 1.77120743, 0.49351805])
        sys_tf = control.tf(numr, denr)
        designer1 = scld.SISOClosedLoopDesigner(sys_tf)
        designer1.design(kp_spec=True, ki_spec=True, residual_precision=2)
        designer2 = scld.SISOClosedLoopDesigner(sys_tf)
        designer2.design(kp=True, ki=True, residual_precision=5)
        self.assertLess(designer1.kp, designer2.kp)

    def testSimulate(self):
        if IGNORE_TEST:
            return
        def calcDiff(arr):
            return np.abs(arr[-1] - 1)
        #
        sys_tf = control.tf([1], [1, 1])
        designer = scld.SISOClosedLoopDesigner(sys_tf)
        designer.set(kp=20)
        _, prediction1s = designer.simulateTransferFunction()
        designer.set(kp=20, ki=50)
        _, prediction2s = designer.simulateTransferFunction()
        self.assertLess(calcDiff(prediction2s), calcDiff(prediction1s))

    def testPlot(self):
        if IGNORE_TEST:
            return
        self.designer.set(**PARAMETER_DCT)
        times = np.linspace(0, 50, 500)
        self.designer.plot(times=times, is_plot=IS_PLOT, markers=["", ""])
        self.designer.set(kp=10, ki=5)
        self.designer.plot(times=times, is_plot=IS_PLOT)
        designer = self.designer.history.get(1)
        designer.set(kp=1)
        designer.plot(times=times, is_plot=IS_PLOT, markers=["", ""])
    
    def testEvaluatePlotResult(self):
        if IGNORE_TEST:
            return
        self.designer.set(kp=1)
        self.designer.evaluateNonlinearIOSystemClosedLoop(self.ctlsb, is_plot=IS_PLOT)

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
        designer = ctl.SISOClosedLoopDesigner(sys_tf, times=TIMES, step_size=5)
        designer.set(kp=2, ki=3)
        closed_loop_tf = designer.closed_loop_tf
        func1 = lambda x: float(closed_loop_tf(x))
        cltf_nums = cltf.subs({kp: 2, ki: 3})
        func2 = lambda x: sympy.N(cltf_nums.subs({s: x}))
        self.assertTrue(ctl.util.compareSingleArgumentFunctions(func1, func2, 0, 100))

    # FIXME: Finish tests
    def testBug1(self):
        if IGNORE_TEST:
            return
        # Setup
        ctlsb = ctl.ControlSBML(LINEAR_MDL, input_names=["k0"], output_names=["S3"])
        linear_bldr = ctlsb.makeSISOTransferFunctionBuilder(is_fixed_input_species=True)
        linear_staircase = ctl.Staircase(initial_value=0, final_value=10, num_step=5)
        fitter_result = linear_bldr.fitTransferFunction(num_numerator=3, num_denominator=4, 
                                                    staircase=linear_staircase, fit_start_time=20,
                                                start_time=0, end_time=100)
        linear_tf = fitter_result.transfer_function
        numr = linear_tf.num[0][0]
        denr = linear_tf.den[0][0]
        #
        designer = ctl.SISOClosedLoopDesigner(linear_tf, times=TIMES, step_size=5)
        designer.set(kp=2, ki=3)
        if IS_PLOT:
            designer.plot(figsize=(5,5), times=TIMES)  # Have options for a period



#############################
class TestHistory(unittest.TestCase):

    def setUp(self):
        self.sys_tf = TRANSFER_FUNCTION
        self.designer = scld.SISOClosedLoopDesigner(self.sys_tf)
        self.history = scld._History(self.designer)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertEqual(len(self.history), 0)
        diff = set(scld.PARAM_NAMES) - set(self.history._dct.keys())
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
        self.assertTrue(isinstance(designer, scld.SISOClosedLoopDesigner))

    def testClear(self):
        if IGNORE_TEST:
            return
        self.addItems(3)
        self.history.clear()
        self.assertEqual(len(self.history), 0)


if __name__ == '__main__':
    unittest.main()

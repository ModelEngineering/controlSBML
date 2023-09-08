from controlSBML import util
from controlSBML.timeseries import Timeseries
import controlSBML.siso_closed_loop_designer as scld
import controlSBML.constants as cn
import controlSBML as ctl
import helpers

import control
import pandas as pd
import numpy as np
import tellurium as te
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
# Construct a transfer function for the model. This is a linear model, and so it should be accurate.
ctlsb = ctl.ControlSBML(MODEL, input_names=["S0"], output_names=["S2"])
builder = ctlsb.makeSISOTransferFunctionBuilder(is_fixed_input_species=True)
staircase = ctl.Staircase(final_value=15, num_step=5)
fitter_result = builder.fitTransferFunction(num_numerator=2, num_denominator=3, staircase=staircase)
TRANSFER_FUNCTION = fitter_result.transfer_function
if False:
    builder.plotFitTransferFunction(fitter_result)
PARAMETER_DCT = {p: n+1 for n, p in enumerate(scld.PARAM_NAMES)}


#############################
# Tests
#############################
class TestSISOClosedLoopDesigner(unittest.TestCase):

    def setUp(self):
        self.sys_tf = TRANSFER_FUNCTION
        self.designer = scld.SISOClosedLoopDesigner(self.sys_tf)

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
        sys_tf = control.tf([1], [1, 1])
        designer = scld.SISOClosedLoopDesigner(sys_tf)
        with self.assertRaises(ValueError):
            _ = designer.closed_loop_tf()
        #
        designer.kd = 2
        designer.kf = 4
        closed_loop_tf = designer.closed_loop_tf
        numr = np.array(closed_loop_tf.num[0][0])
        self.assertTrue(np.allclose(numr, [2, 8, 0]))
        denr = np.array(closed_loop_tf.den[0][0])
        self.assertTrue(np.allclose(denr, [1, 13, 4]))

    def testDesign1(self):
        if IGNORE_TEST:
            return
        sys_tf = control.tf([1], [1, 1])
        designer = scld.SISOClosedLoopDesigner(sys_tf)
        designer.design(kd=True, kf=True)
        for name in ["kd", "kf"]:
            self.assertIsNotNone(getattr(designer, name))
        for name in ["ki", "kp"]:
            self.assertIsNone(getattr(designer, name))

    def testSimulate(self):
        if IGNORE_TEST:
            return
        def calcDiff(arr):
            return np.abs(arr[-1] - 1)
        #
        sys_tf = control.tf([1], [1, 1])
        designer = scld.SISOClosedLoopDesigner(sys_tf)
        designer.set(kp=20)
        prediction1s = designer.simulate()
        designer.set(kp=20, ki=50)
        prediction2s = designer.simulate()
        self.assertLess(calcDiff(prediction2s), calcDiff(prediction1s))

    def testPlot(self):
        if IGNORE_TEST:
            return
        self.designer.set(**PARAMETER_DCT)
        self.designer.plot(is_plot=IS_PLOT)


if __name__ == '__main__':
    unittest.main()
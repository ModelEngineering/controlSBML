from controlSBML import util
from controlSBML.timeseries import Timeseries
from controlSBML.siso_closed_loop_designer import SISOClosedLoopDesigner
import controlSBML.constants as cn
import controlSBML as ctl
import helpers

import control
import pandas as pd
import numpy as np
import tellurium as te
import unittest

IGNORE_TEST = True
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
if IS_PLOT:
    builder.plotFitTransferFunction(fitter_result)


#############################
# Tests
#############################
class TestSISOClosedLoopDesigner(unittest.TestCase):

    def setUp(self):
        self.sys_tf = TRANSFER_FUNCTION
        self.designer = SISOClosedLoopDesigner(self.sys_tf)

    def testDesignSimple(self):
        #if IGNORE_TEST:
        #    return
        tf = control.tf([1], [1, 1])
        designer = SISOClosedLoopDesigner(tf)
        designer.design(kp=1)
        import pdb; pdb.set_trace()

    def testDesign(self):
        if IGNORE_TEST:
            return
        self.designer.design(ki=1)
        import pdb; pdb.set_trace()

if __name__ == '__main__':
    unittest.main()
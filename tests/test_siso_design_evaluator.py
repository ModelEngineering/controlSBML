from controlSBML.grid import Grid, Axis
from controlSBML.siso_design_evaluator import SISODesignEvaluator

from controlSBML import util
from controlSBML.sbml_system import SBMLSystem
import controlSBML.constants as cn

import pandas as pd
import numpy as np
import os
import tellurium as te
import unittest


IGNORE_TEST = False
IS_PLOT = False
TIMES = np.linspace(0, 100, 1000)
MODEL_UNSTABLE = """
model *simple()
S1 -> S2; k1*S1

S1 = 10
S2 = 0
k1 = 0.1
end
"""
MODEL_STABLE = """
model *simple()
S1 -> S2; k1*S1
S2 -> ; k2*S2

S1 = 10
S2 = 0
k1 = 0.1
k2 = 1
end
"""
INPUT_NAME = "S1"
OUTPUT_NAME = "S2"
SYSTEM_STABLE = SBMLSystem(MODEL_STABLE, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME])
SYSTEM_UNSTABLE = SBMLSystem(MODEL_UNSTABLE, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME])


#############################
# Tests
#############################
class TestSISOClosedLoopDesigner(unittest.TestCase):

    def setUp(self):
        self.evaluator_stable = SISODesignEvaluator(SYSTEM_STABLE, INPUT_NAME, OUTPUT_NAME, times=TIMES)
        self.evaluator_unstable = SISODesignEvaluator(SYSTEM_UNSTABLE, INPUT_NAME, OUTPUT_NAME, times=TIMES)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.evaluator_stable, SISODesignEvaluator))

    def testCalculateMse(self):
        if IGNORE_TEST:
            return
        parameter_dct = {cn.CP_KP: 1, cn.CP_KI: 1, cn.CP_KF: 1}
        is_feasible, mse = self.evaluator_unstable._calculateMse(**parameter_dct)
        self.assertFalse(is_feasible)
        #
        is_feasible, mse = self.evaluator_stable._calculateMse(**parameter_dct)
        self.assertTrue(is_feasible)
        self.assertTrue(mse > 0)

    def testEvaluate(self):
        if IGNORE_TEST:
            return
        parameter_dct = {cn.CP_KP: 1, cn.CP_KI: 1, cn.CP_KF: 1}
        self.evaluator_stable.evaluate(**parameter_dct)
        self.assertTrue(self.evaluator_stable.residual_mse > 0)
        #
        parameter_dct = {cn.CP_KP: 1, cn.CP_KI: 1, cn.CP_KF: 1}
        self.evaluator_unstable.evaluate(**parameter_dct)
        self.assertTrue(self.evaluator_unstable.residual_mse is None)


if __name__ == '__main__':
    unittest.main()
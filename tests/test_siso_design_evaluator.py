from controlSBML.grid import Grid, Axis
from controlSBML.siso_design_evaluator import SISODesignEvaluator, EvaluatorResult

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
FILE_SAVE = os.path.join(cn.TEST_DIR, "siso_design_evaluator.csv")
DF = pd.DataFrame({cn.CP_KP: [1, 2], cn.CP_KI: [1, 2], cn.MSE: [1, 2]})


#############################
def remove():
    if os.path.isfile(FILE_SAVE):
        os.remove(FILE_SAVE)


#############################
# Tests
#############################
class TestSISOEvaluatorResult(unittest.TestCase):

    def setUp(self):
        self.evaluator_result = EvaluatorResult(kp_spec=True, ki_spec=True, kf_spec=False)
        remove()

    def tearDown(self):
        remove()

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.evaluator_result, EvaluatorResult))
        self.assertTrue(cn.MSE in self.evaluator_result)
        self.assertTrue(cn.CP_KP in self.evaluator_result)
        self.assertTrue(cn.CP_KI in self.evaluator_result)
        self.assertFalse(cn.CP_KF in self.evaluator_result)

    def testAdd(self):
        if IGNORE_TEST:
            return
        self.evaluator_result.add(kp=1, ki=1, mse=1)
        self.evaluator_result.add(kp=1, ki=1, mse=1)
        self.assertTrue(len(self.evaluator_result[cn.MSE]) == 2)

    def testMakeFromDataframe(self):
        if IGNORE_TEST:
            return
        evaluator_result = EvaluatorResult.makeFromDataframe(DF)
        self.assertTrue(len(evaluator_result[cn.MSE]) == 2)

    def testWriteCsv(self):
        if IGNORE_TEST:
            return
        evaluator_result = EvaluatorResult.makeFromDataframe(DF)
        evaluator_result.writeCsv(FILE_SAVE)
        self.assertTrue(os.path.isfile(FILE_SAVE))

    def testCopy(self):
        if IGNORE_TEST:
            return
        evaluator_result = EvaluatorResult.makeFromDataframe(DF)
        copy_result = evaluator_result.copy()
        for key, value in evaluator_result.items():
            if value is not None:
                self.assertTrue(np.allclose(value, copy_result[key]))
            else:
                self.assertIsNone(copy_result[key])

    def testMerge(self):
        if IGNORE_TEST:
            return
        evaluator_result = EvaluatorResult.makeFromDataframe(DF)
        merged_result = self.evaluator_result.merge([evaluator_result, evaluator_result])
        self.assertEqual(len(merged_result[cn.MSE]), 2*len(DF))


class TestSISODesignEvaluator(unittest.TestCase):

    def setUp(self):
        self.evaluator_stable = SISODesignEvaluator(SYSTEM_STABLE, INPUT_NAME, OUTPUT_NAME, times=TIMES, save_path=FILE_SAVE)
        self.evaluator_unstable = SISODesignEvaluator(SYSTEM_UNSTABLE, INPUT_NAME, OUTPUT_NAME, times=TIMES, save_path=FILE_SAVE)
        remove()

    def tearDown(self):
        remove()

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

    def testGetDesignResults(self):
        if IGNORE_TEST:
            return
        parameter_dct = {cn.CP_KP: 1, cn.CP_KI: 1, cn.CP_KF: 1}
        self.evaluator_stable.evaluate(**parameter_dct)
        df = self.evaluator_stable.evaluator_result.getDataframe()
        self.assertTrue(df is not None)
        self.assertTrue(len(df) == 1)

    def testEvaluate(self):
        if IGNORE_TEST:
            return
        self.assertFalse(os.path.isfile(FILE_SAVE))
        parameter_dct = {cn.CP_KP: 1, cn.CP_KI: 1, cn.CP_KF: 1}
        self.evaluator_stable.evaluate(**parameter_dct)
        self.assertTrue(os.path.isfile(FILE_SAVE))
        self.assertTrue(self.evaluator_stable.residual_mse > 0)
        #
        parameter_dct = {cn.CP_KP: 1.1, cn.CP_KI: 1, cn.CP_KF: 1}
        self.evaluator_stable.evaluate(**parameter_dct)
        self.assertIsNotNone(self.evaluator_stable.evaluator_result)
        self.assertEqual(len(self.evaluator_stable.evaluator_result[cn.CP_KP]), 2)
        #
        parameter_dct = {cn.CP_KP: 1, cn.CP_KI: 1, cn.CP_KF: 1}
        self.evaluator_unstable.evaluate(**parameter_dct)
        self.assertTrue(self.evaluator_unstable.residual_mse is None)


if __name__ == '__main__':
    unittest.main()
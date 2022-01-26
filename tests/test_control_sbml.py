from controlSBML.control_sbml import ControlSBML
import helpers

import numpy as np
import pandas as pd
import os
import unittest
import tellurium as te


IGNORE_TEST = False
IS_PLOT = False
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
ANTIMONY_FILE = os.path.join(TEST_DIR, "Model_antimony.ant")
LINEAR_MODEL = """
S0 -> 2 S0; S0
S0 -> S1; S0
S1 -> S2; S1
S2 -> S3; S2

S0 = 1
S1 = 10
S2 = 0
S3 = 0
"""


#############################
# Tests
#############################
class TestControlSBML(unittest.TestCase):

    def setUp(self):
      # Cannot modify self.control
      self.ctlsb = ControlSBML(ANTIMONY_FILE)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue("RoadRunner" in str(type(self.ctlsb.roadrunner)))

    def testConstructWithRoadrunner(self):
        if IGNORE_TEST:
            return
        model = te.loada(helpers.TEST_PATH_1)
        ctlsb1 = ControlSBML(model)
        ctlsb2 = ControlSBML(helpers.TEST_PATH_1)
        diff = set(ctlsb1.get().values()).symmetric_difference(
              ctlsb2.get().values())
        self.assertEqual(len(diff), 0)

    def testGet(self):
        if IGNORE_TEST:
            return
        S0 = "S0"
        dct = self.ctlsb.get(S0)
        self.assertEqual(dct, self.ctlsb.roadrunner[S0])
        #
        dct = self.ctlsb.get([S0])
        self.assertEqual(dct[S0], self.ctlsb.roadrunner[S0])
        #
        dct = self.ctlsb.get()
        self.assertGreater(len(dct), 0)

    def testSet(self):
        if IGNORE_TEST:
            return
        VALUE = 2
        self.ctlsb.set({"S0": VALUE})
        value = self.ctlsb.get("S0")
        self.assertTrue(np.isclose(VALUE, value))

    def testCopyEquals(self):
        if IGNORE_TEST:
            return
        ctlsb = self.ctlsb.copy()
        self.assertTrue(ctlsb.equals(self.ctlsb))
        ctlsb.antimony = ""
        self.assertFalse(ctlsb.equals(self.ctlsb))

    def testMkStateSpace(self):
        if IGNORE_TEST:
            return
        sys = self.ctlsb.mkStateSpace()
        self.assertEqual(sys.nstates, 3)

    def testMkInitialState(self):
        if IGNORE_TEST:
            return
        X0 = self.ctlsb.mkInitialState()
        self.assertEqual(np.shape(X0), (3,))
        dX0 = np.matmul(self.ctlsb.jacobian, X0)
        self.assertEqual(len(X0), len(dX0))

    def testSetTime(self):
        if IGNORE_TEST:
            return
        def isEqual(jac1, jac2, val_bl):
            mat = jac1 - jac2
            diff = sum((mat.flatten())**2)
            self.assertEqual(np.isclose(diff, 0), val_bl)
        #
        self.ctlsb.setTime(0)
        jac_0 = self.ctlsb.jacobian
        self.ctlsb.setTime(5)
        jac_5 = self.ctlsb.jacobian
        isEqual(jac_0, jac_5, False)
        #
        self.ctlsb.setTime(0)
        jac_00 = self.ctlsb.jacobian
        isEqual(jac_00, jac_0, True)

    def testSimulateLinearSystem(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(LINEAR_MODEL)
        linear_df = ctlsb.simulateLinearSystem()
        rr = te.loada(LINEAR_MODEL)
        data = rr.simulate()
        columns = [c[1:-1] for c in data.colnames]
        rr_df = pd.DataFrame(data, columns=columns)
        for column in rr_df.columns[1:]:
            linear_arr = linear_df[column].values
            rr_arr = rr_df[column].values
            self.assertLess(np.abs(linear_arr[-1] - rr_arr[-1]), 0.1)
        
        


if __name__ == '__main__':
  unittest.main()

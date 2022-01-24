from controlSBML import constants as cn
from controlSBML import control_sbml
from controlSBML.control_sbml import ControlSBML
from tests import helpers

import copy
import numpy as np
import os
import unittest
import tellurium as te


IGNORE_TEST = True
IS_PLOT = True
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
ANTIMONY_FILE = os.path.join(TEST_DIR, "Model_antimony.ant")


#############################
# Tests
#############################
class TestControlSBML(unittest.TestCase):

    def setUp(self):
      # Cannot modify self.control
      self.ctlsb = ControlSBML(ANTIMONY_FILE)
  
    def testConstructor(self):
        # TESTING
        self.assertTrue("RoadRunner" in type(self.ctlsb.roadrunner))
  
    def testConstructWithRoadrunner(self):
        if IGNORE_TEST:
            return
        model = te.loada(helpers.TEST_PATH_1)
        ctlsb1 = SymmathSBML(model)
        ctlsb2 = SymmathSBML(helpers.TEST_PATH_1)
        diff = set(ctlsb1.get().values()).symmetric_difference(ctlsb2.get().values())
        self.assertGreater(len(diff), 0)

    def testGet(self):
        if IGNORE_TEST:
            return
        dct = self.ctlsb.get("S0")
        import pdb; pdb.set_trace()

    def testSet(self):
        if IGNORE_TEST:
            return
        VALUE = 2
        self.ctlsb.set("S0", VALUE)
        value = self.ctlsb_get("S0")
        self.assertTrue(np.isclose(VALUE, value)
        import pdb; pdb.set_trace()

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
        x0 = self.ctlsb.mkInitialState()
        self.assertEqual(np.shape(x0), (3, 1))
    

if __name__ == '__main__':
  unittest.main()

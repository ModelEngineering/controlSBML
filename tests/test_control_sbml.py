from controlSBML.control_sbml import ControlSBML
from controlSBML import control_sbml
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
MODEL_FILE = os.path.join(TEST_DIR, "BIOMD0000000823.xml")
LINEAR_MDL = """
J0: $S0 -> S1; k0*$S0
J1: S1 -> S2; k1*S1
J2: S2 -> S3; k2*S2

S0 = 1
S1 = 10
S2 = 0
S3 = 0
k0 = 1
k1 = 1
k2 = 1
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

    def testMakeSISOTransferFunctionBuilder(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(MODEL_FILE,
            input_names=["IR"], output_names=["mTORC1_DEPTOR"])
        builder = ctlsb.makeSISOTransferFunctionBuilder()
        self.assertTrue("Builder" in str(type(builder)))
    
    def testMakeNonlinearIOSystem(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(LINEAR_MDL, input_names=["S1"])
        non_sys = ctlsb.makeNonlinearIOSystem("tst")
        self.assertTrue("NonlinearIOSystem" in str(type(non_sys)))


if __name__ == '__main__':
  unittest.main()

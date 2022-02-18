from controlSBML.control_analysis import ControlAnalysis
from controlSBML import control_analysis
import helpers

import numpy as np
import pandas as pd
import os
import unittest
import tellurium as te


IGNORE_TEST = False
IS_PLOT = False
if IS_PLOT:
    import matplotlib
    matplotlib.use('TkAgg')

HTTP_FILE = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000206.2?filename=BIOMD0000000206_url.xml"
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
ANTIMONY_FILE = os.path.join(TEST_DIR, "Model_antimony.ant")
LINEAR_MDL = """
$S0 -> S1; $S0
S1 -> S2; S1
S2 -> S3; S2

S0 = 1
S1 = 10
S2 = 0
S3 = 0
"""
NONLINEAR_MDL = """
S0 -> 2 S0; S0
S0 -> S1; S0
S1 -> S2; S1*S1
S2 -> S3; S2*S1

S0 = 1
S1 = 10
S2 = 0
S3 = 0
"""
NONLINEAR1_MDL = """
//S0 -> 2 S0; S0
$S0 -> S1; $S0
S1 -> S2; k2*S1
S2 -> ; k3*S2*S1

k1 = 1;
k2 = 1
k3 = 1
k4 = 1
S0 = 1
k0 = 1
S1 = 10
S2 = 0
S3 = 0
"""


#############################
# Tests
#############################
class TestControlAnalysis(unittest.TestCase):

    def setUp(self):
      # Cannot modify self.control
      self.ctlsb = ControlAnalysis(ANTIMONY_FILE)

    def testSimulateLinearSystemRoadrunner(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlAnalysis(LINEAR_MDL)
        linear_df = ctlsb.simulateLinearSystem()
        rr = te.loada(LINEAR_MDL)
        rr_df = ctlsb.simulateRoadrunner()
        for column in rr_df.columns[1:]:
            linear_arr = linear_df[column].values
            rr_arr = rr_df[column].values
            self.assertLess(np.abs(linear_arr[-1] - rr_arr[-1]), 0.1)


if __name__ == '__main__':
  unittest.main()

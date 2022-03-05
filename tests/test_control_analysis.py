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
J0: $S0 -> S1; $S0
J1: S1 -> S2; S1
J2: S2 -> S3; S2

$S0 = 0
S1 = 10
S2 = 0
S3 = 0
"""
LONG_LINEAR_MDL = """
J0: S0 -> S1; S0
J1: S1 -> S2; S1
J2: S2 -> S3; S2
J3: S2 -> S3; S3

S0 = 10
S1 = 0
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
        if IGNORE_TEST:
          return
        self.ctlsb = ControlAnalysis(ANTIMONY_FILE)

    def testSimulateLinearSystemRoadrunner(self):
        if IGNORE_TEST:
          return
        ctlsb = ControlAnalysis(LINEAR_MDL)
        linear_df = ctlsb.simulateLinearSystem(step_val=0)
        rr_df = ctlsb.simulateRoadrunner()
        for column in linear_df.columns[1:]:
            linear_arr = linear_df[column].values
            rr_arr = rr_df[column].values
            self.assertLess(np.abs(linear_arr[-1] - rr_arr[-1]), 0.1)

    def testReducedStoich(self):
        if IGNORE_TEST:
          return
        ctlsb = ControlAnalysis(LONG_LINEAR_MDL)
        #import pdb; pdb.set_trace()

    def testinearApproximationNonzeroInput(self):
        if IGNORE_TEST:
          return
        step_val = 2
        ctlsb = ControlAnalysis(LINEAR_MDL, input_names=["J0"])
        ctlsb.setTime(2)
        ctlsb.set({"S0": step_val})
        rr_df = ctlsb.simulateRoadrunner()
        linear_df = ctlsb.simulateLinearSystem(step_val=2)
        squared_difference = ((rr_df-linear_df)**2).sum().sum()
        self.assertTrue(np.isclose(squared_difference, 0))

    def testSimulateLinearSystem2(self):
        if IGNORE_TEST:
          return
        ctlsb = ControlAnalysis(
              "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000823.2?filename=Varusai2018.xml")
        timepoint = 1
        linear_df = ctlsb.simulateLinearSystem(end_time=10, timepoint=timepoint)
        self.assertTrue(isinstance(linear_df, pd.DataFrame))


if __name__ == '__main__':
  unittest.main()

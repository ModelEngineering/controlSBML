from controlSBML.control_analysis import ControlAnalysis
from controlSBML import control_analysis
import helpers

import matplotlib.pyplot as plt
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

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
ANTIMONY_FILE = os.path.join(TEST_DIR, "Model_antimony.ant")
MODEL_FILE = os.path.join(TEST_DIR, "BIOMD0000000823.xml")
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

    def plot(self, df1, df2):
        for name in df1.columns:
            randoms = np.random.rand(len(df1[name]))
            plt.plot(df1.index, df1[name]+0.1*randoms)
            plt.plot(df2.index, df2[name])
        plt.show()

    def testSimulateLinearSystemRoadrunner(self):
        if IGNORE_TEST:
          return
        ctlsb = ControlAnalysis(LINEAR_MDL)
        linear_df = ctlsb.simulateLinearSystem(step_val=0)
        rr_df = ctlsb.simulateRoadrunner()
        ams = np.sqrt(np.sum((rr_df.to_numpy() - linear_df.to_numpy())**2))
        ams = ams/len(rr_df)
        self.assertLess(np.abs(ams), 1e-4)
        if IS_PLOT:
            self.plot(linear_df, rr_df)

    def testReducedStoich(self):
        if IGNORE_TEST:
          return
        ctlsb = ControlAnalysis(LONG_LINEAR_MDL)
        #import pdb; pdb.set_trace()

    def testinearApproximationNonzeroInput(self):
        if IGNORE_TEST:
          return
        step_val = 2
        ctlsb = ControlAnalysis(LINEAR_MDL, input_names=["S1"])
        ctlsb.setTime(2)
        ctlsb.set({"S0": step_val})
        rr_df = ctlsb.simulateRoadrunner()
        linear_df = ctlsb.simulateLinearSystem(step_val=2)
        squared_difference = (
              (rr_df.df - linear_df.df)**2).sum().sum()
        self.assertTrue(np.isclose(squared_difference, 0))

    def testSimulateLinearSystem2(self):
        if IGNORE_TEST:
          return
        ctlsb = ControlAnalysis(MODEL_FILE)
        time = 1
        linear_df = ctlsb.simulateLinearSystem(end_time=10, time=time)
        self.assertTrue(isinstance(linear_df, pd.DataFrame))


if __name__ == '__main__':
  unittest.main()

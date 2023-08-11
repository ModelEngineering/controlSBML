from controlSBML.control_analysis import ControlAnalysis
import controlSBML.constants as cn
import helpers

import control
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import unittest
import tellurium as te


IGNORE_TEST = False
IS_PLOT = False
helpers.setupPlotting(__file__)
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
LINEAR2_MDL = """
J0: $S0 -> S1; $S0
J1: S1 -> S2; S1
J2: S2 -> S3; S2
J3: S3 -> S4; S3

S0 = 10
S1 = 0
S2 = 0
S3 = 0
S4 = 0
"""
#
LINEAR3_MDL = """
species $S0, S1, S2, S3

  -> S1; S0
J1a: S2 -> S1; k2*S2
J1b: S1 -> S2; S1
J2: S2 -> S3; S2

k = 10
k2 = 2
S1 = 5
S2 = 0
S3 = 0
S0 = k
"""
LINEAR3_MDL_S0 = 10
LINEAR3_MDL_k2 = 2
LINEAR3_MDL_S1 = 5
#
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

    def testSimulateRoadrunner(self):
        # Compares simulation result with expectation for a known transfer function when there is an initial condition.
        if IGNORE_TEST:
          return
        ctlsb = ControlAnalysis(LINEAR3_MDL, input_names=["S0"], output_names=["S1", "S2"])
        ctlsb.set({"S0": LINEAR3_MDL_S0, "k2": LINEAR3_MDL_k2, "S1": LINEAR3_MDL_S1})
        rr_ts = ctlsb.simulateRoadrunner()
        diff = set(["S1", "S2", "S3"]).symmetric_difference(rr_ts.columns)
        self.assertEqual(len(diff), 0)

    def testSimulateLinearSystemIdealized(self):
        # Compares simulation result with expectation for a known transfer function when there is an initial condition.
        if IGNORE_TEST:
          return
        ctlsb = ControlAnalysis(LINEAR3_MDL, input_names=["S0"], output_names=["S1", "S2"])
        ctlsb.set({"S0": LINEAR3_MDL_S0, "k2": LINEAR3_MDL_k2, "S1": LINEAR3_MDL_S1})
        rr_ts = ctlsb.simulateRoadrunner()
        rr_ts = rr_ts[["S1"]]
        # Transfer function
        tf = control.TransferFunction([1, 1+ LINEAR3_MDL_k2], [1, 2 + LINEAR3_MDL_k2, 1])
        times = rr_ts.index/cn.MS_IN_SEC
        linear_ts = ctlsb.simulateLinearSystem(tf, times=times,
                                                step_val=LINEAR3_MDL_S0, input_name="S0", output_name="S1")["S0"]
        ams = np.sqrt(np.sum((rr_ts.to_numpy() - linear_ts.to_numpy())**2))
        ams = ams/len(rr_ts)
        self.assertLess(np.abs(ams), 1e-4)

    # FIXME: (1) Not getting correct answer for fixed species
    #       (2) Not correctly handling fixed rates with initial conditions.
    def testSimulateLinearSystemFixedSpecies(self):
        if IGNORE_TEST:
          return
        return
        ctlsb = ControlSBML(LINEAR2_MDL, input_names=["S0"], output_names=["S1", "S2", "S3"])
        rr_ts = ctlsb.simulateRoadrunner()
        # No additional external input
        transfer_function_df = ctlsb.makeMIMOTransferFunctionDF(staircase=Staircase(final_value=0),
                                                                is_steady_state=False, is_fixed_input_species=True)
        linear_ts = ctlsb.simulateLinearSystem(transfer_function_df, step_val=10)
        import pdb; pdb.set_trace()
        linear_ts.columns = rr_ts.columns # Ensure have the same columns
        ams = np.sqrt(np.sum((rr_ts.to_numpy() - linear_ts.to_numpy())**2))
        ams = ams/len(rr_ts)
        self.assertLess(np.abs(ams), 1e-4)


if __name__ == '__main__':
  unittest.main()

from controlSBML.control_plot import ControlPlot  # type: ignore
from controlSBML.control_sbml import ControlSBML # type: ignore

import control  # type: ignore
import pandas as pd  # type: ignore
import os
import unittest


IGNORE_TEST = False
IS_PLOT = False
#helpers.setupPlotting(__file__)

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
ANTIMONY_FILE = os.path.join(TEST_DIR, "Model_antimony.ant")
MODEL_FILE = os.path.join(TEST_DIR, "tests/BIOMD0000000206.xml")
LINEAR_MDL = """
J0: $S0 -> S1; $S0
J1: S1 -> S2; S1
J2: S2 -> S3; S2

$S0 = 0
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
S2 = 5
S3 = 5
S0 = k
"""
LINEAR3_MDL_S0 = 10
LINEAR3_MDL_k2 = 2
LINEAR3_MDL_S1 = 5


#############################
# Tests
#############################
class TestControlPlot(unittest.TestCase):

    def setUp(self):
      # Cannot modify self.control
      self.ctlsb = ControlPlot(ANTIMONY_FILE)

    def testPlotTrueModel(self):
        if IGNORE_TEST:
          return
        self.ctlsb.plotTrueModel(is_plot=IS_PLOT, ylabel="values",
              end_time=10, title="title", figsize=(5, 10))

    def testPlotAccuracy(self):
        #if IGNORE_TEST:
        #  return
        ctlsb = ControlPlot(LINEAR3_MDL, input_names=["S0", "S2"], output_names=["S1", "S3"])
        tf = control.TransferFunction([1, 1+ LINEAR3_MDL_k2], [1, 2 + LINEAR3_MDL_k2, 1])
        transfer_function_df = pd.DataFrame({"S1": [tf, tf], "S3": [tf, tf]})
        transfer_function_df.index = ["S0", "S2"]
        ctlsb.plotAccuracy(transfer_function_df, suptitle="Test", is_plot=IS_PLOT, figsize=(5,5))

    def testPlotBode(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(NONLINEAR_MDL, input_names=["S1", "S2"],
              output_names=["S0", "S2"])
        ctlsb.plotBode(is_plot=IS_PLOT)


if __name__ == '__main__':
  unittest.main()

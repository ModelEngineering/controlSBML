from controlSBML.control_plot import ControlPlot
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
class TestControlPlot(unittest.TestCase):

    def setUp(self):
      # Cannot modify self.control
      self.ctlsb = ControlPlot(ANTIMONY_FILE)

    def testPlotLinearApproximation(self):
        if IGNORE_TEST:
          return
        ctlsb = ControlPlot(NONLINEAR_MDL)
        ctlsb.setTime(2)
        A_df = ctlsb.jacobian
        ctlsb.plotLinearApproximation(A_mat=A_df, suptitle="Test",
              is_plot=IS_PLOT, figsize=(5,5))

    def testPlotTrueModel(self):
        if IGNORE_TEST:
          return
        self.ctlsb.plotTrueModel(is_plot=IS_PLOT, ylabel="values",
              end_time=10, title="title")

    def testPlotAccuracy(self):
        if IGNORE_TEST:
          return
        ctlsb = ControlPlot(NONLINEAR_MDL)
        ctlsb.plotAccuracy(figsize=(5, 5), is_plot=IS_PLOT)
        self.ctlsb.plotAccuracy(NONLINEAR_MDL,
              [0, 1, 2, 3], suptitle="Test", is_plot=IS_PLOT)


if __name__ == '__main__':
  unittest.main()

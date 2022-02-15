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

    def testMakeStateSpace(self):
        if IGNORE_TEST:
            return
        sys = self.ctlsb.makeStateSpace()
        self.assertEqual(sys.nstates, 3)

    def testMkInitialState(self):
        if IGNORE_TEST:
            return
        def test(mdl):
            ctlsb = ControlSBML(mdl)
            ctlsb = ControlSBML(ctlsb.roadrunner)
            num_species = ctlsb.roadrunner.model.getNumFloatingSpecies()  \
                  + ctlsb.roadrunner.model.getNumBoundarySpecies()
            X0 = ctlsb.current_state
            self.assertEqual(len(X0), num_species)
            jacobian = ctlsb.jacobian
            self.assertTrue(isinstance(jacobian, pd.DataFrame))
            self.assertEqual(len(jacobian), len(jacobian.columns))
            dX0 = np.matmul(ctlsb.jacobian.values, X0.values)
            self.assertEqual(len(X0), len(dX0))
        #
        test(NONLINEAR1_MDL)
        test(LINEAR_MDL)

    def testSetTime(self):
        if IGNORE_TEST:
            return
        def isEqual(jac1, jac2, val_bl):
            mat = jac1 - jac2
            diff = sum((mat.values.flatten())**2)
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

    def testSimulateLinearSystemRoadrunner(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(LINEAR_MDL)
        linear_df = ctlsb.simulateLinearSystem()
        rr = te.loada(LINEAR_MDL)
        rr_df = ctlsb.simulateRoadrunner()
        for column in rr_df.columns[1:]:
            linear_arr = linear_df[column].values
            rr_arr = rr_df[column].values
            self.assertLess(np.abs(linear_arr[-1] - rr_arr[-1]), 0.1)

    def testEvaluateAccuracy(self):
        if IGNORE_TEST:
          return
        self.ctlsb.evaluateAccuracy(NONLINEAR_MDL,
              [0, 1, 2, 3], suptitle="Test", is_plot=IS_PLOT)

    def testPlotLinearApproximation(self):
        if IGNORE_TEST:
          return
        ctlsb = ControlSBML(NONLINEAR_MDL)
        ctlsb.setTime(2)
        A_mat = ctlsb.jacobian
        ctlsb.plotLinearApproximation(A_mat, suptitle="Test", is_plot=IS_PLOT,
              figsize=(5,5))

    def testPlotTrueModel(self):
        if IGNORE_TEST:
          return
        self.ctlsb.plotTrueModel(is_plot=IS_PLOT, ylabel="values")

    def testParseOptions(self):
        if IGNORE_TEST:
          return
        def isSameDct(dct1, dct2):
            diff = set(dct1.keys()).symmetric_difference(dct2.keys())
            self.assertEqual(len(diff), 0)
            diff = set(dct1.values()).symmetric_difference(dct2.values())
            self.assertEqual(len(diff), 0)
        #
        kwargs = dict(control_sbml.PLOT_OPTS)
        kwargs.update(control_sbml.FIG_OPTS)
        kwargs.update(control_sbml.SIM_OPTS)
        plot_opts, fig_opts, sim_opts = self.ctlsb._parseOpts(**kwargs)
        isSameDct(fig_opts, control_sbml.FIG_OPTS)
        isSameDct(plot_opts, control_sbml.PLOT_OPTS)
        isSameDct(sim_opts, control_sbml.SIM_OPTS)


if __name__ == '__main__':
  unittest.main()

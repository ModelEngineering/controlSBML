from controlSBML.control_base import ControlBase
import helpers

import control
import numpy as np
import pandas as pd
import os
import unittest
import tellurium as te


IGNORE_TEST = False
IS_PLOT = False

HTTP_FILE = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000206.2?filename=BIOMD0000000206_url.xml"
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
ANTIMONY_FILE = os.path.join(TEST_DIR, "Model_antimony.ant")
INPUT_NAMES = ["J0"]
OUTPUT_NAMES = ["S3", "S2"]
END_TIME = 20
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
SIMPLE_MDL = """
S1 -> S2; 1
S1 = 10; S2 = 0;
"""
LINEAR_MDL_INPUT_NAMES = ["J0", "J1", "J2"]
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
REDUCABLE_MDL = """
J0: S0 -> S1; S0
J1: S1 -> S2; S1
J2: S2 -> S3; S2
J3: S2 -> S3; S3
J0a: S0a -> S1a; S0a
J1a: S1a -> S2a; S1a
J2a: S2a -> S3a; S2a
J3a: S2a -> S3a; S3a
S0a = 10
S1a = 0
S2a = 0
S3a = 0
S0 = 10
S1 = 0
S2 = 0
S3 = 0
"""

OUTPUT_NAMES_REDUCABLE = [
  "S0a", "S1a", "S2a", "S3a", "S0", "S1", "S2", "S3"]

SPECIES_NAMES = ["S0", "S1", "S2"]

def setList(lst, default):
    if lst is None:
        return default
    return lst


#############################
# Tests
#############################
class TestControlBase(unittest.TestCase):

    def setUp(self):
        # Cannot modify self.control
        if IGNORE_TEST:
          return
        self.init()

    def init(self):
        # Cannot modify self.control
        self.ctlsb = ControlBase(ANTIMONY_FILE)

    def testConstructor(self):
        if IGNORE_TEST:
          return
        self.assertTrue("RoadRunner" in str(type(self.ctlsb.roadrunner)))
        for lst in [self.ctlsb.state_names, self.ctlsb.output_names]:
            diff = set(SPECIES_NAMES).symmetric_difference(lst)
            self.assertEqual(len(diff), 0)

    def testMakeCDF(self):
        if IGNORE_TEST:
          return
        C_df = self.ctlsb._makeCDF()
        num_row = len(self.ctlsb.output_names)
        num_col = self.ctlsb.num_state
        self.assertEqual(np.shape(C_df.values), (num_row, num_col))

    def testConstructWithRoadrunner(self):
        if IGNORE_TEST:
          return
        model = te.loada(helpers.TEST_PATH_1)
        ctlsb1 = ControlBase(model)
        ctlsb2 = ControlBase(helpers.TEST_PATH_1)
        diff = set(ctlsb1.get().values()).symmetric_difference(
              ctlsb2.get().values())
        self.assertEqual(len(diff), 0)

    def test_roadrunner_namespace(self):
        if IGNORE_TEST:
          return
        dct = self.ctlsb.roadrunner_namespace
        not_contained = set(["S1", "S2", "S0"]).difference(dct.keys())
        self.assertEqual(len(not_contained), 0)

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

    def testJacobian(self):
        if IGNORE_TEST:
          return
        self.assertTrue(isinstance(self.ctlsb.jacobian_df, pd.DataFrame))

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

    def testMakeStateSpace1(self):
        if IGNORE_TEST:
          return
        self.init()
        sys = self.ctlsb.makeStateSpace(A_mat=self.ctlsb.jacobian_df.values)
        self.assertEqual(sys.nstates, 3)

    def testMakeStateSpace2(self):
        if IGNORE_TEST:
          return
        ctlsb = ControlBase(SIMPLE_MDL)
        sys = ctlsb.makeStateSpace()
        self.assertEqual(sys.nstates, 1)
        times = [0.1*v for v in range(50)]
        X0 = ctlsb.state_ser.values
        _, y_vals = control.forced_response(sys, T=times, X0=X0)
        self.assertEqual(len(y_vals), 2)  # S1, S2 are outputs
        self.assertEqual(len(y_vals[0]), len(times))

    def testMakeStateSpaceReducable(self):
        if IGNORE_TEST:
          return
        ctlsb = ControlBase(REDUCABLE_MDL, output_names=OUTPUT_NAMES_REDUCABLE)
        sys = ctlsb.makeStateSpace()
        self.assertGreater(np.shape(sys.C)[0], np.shape(sys.A)[0])

    def _simulate(self, u_val, mdl=LINEAR_MDL, input_names=None,
          output_names=None, end_time=END_TIME):
        """
        Simulates the linear model with an input.

        Parameters
        ----------
        u_val: float
            level of step applied to inputs.
        mdl: model_reference

        Returns
        -------
        array-float
        """
        input_names = setList(input_names, INPUT_NAMES)
        output_names = setList(output_names, OUTPUT_NAMES)
        ctlsb = ControlBase(mdl,
              input_names=input_names, output_names=output_names)
        sys = ctlsb.makeStateSpace()
        times = [0.1*v for v in range(end_time)]
        X0 = ctlsb.state_ser.values
        U = np.repeat(u_val, len(times))
        _, y_vals = control.forced_response(sys, T=times, X0=X0, U=U)
        return y_vals

    def testMakeStateSpaceLeveledInputs(self):
        if IGNORE_TEST:
          return
        y_dct = {u: self._simulate(u) for u in [0, 1, 2]}
        keys = y_dct.keys()
        for key in range(len(keys)-1):
            y1_vals = y_dct[key]
            y2_vals = y_dct[key+1]
            for idx in range(2):
                # Compare current and enxt
                self.assertGreater(y2_vals[idx][-1], y1_vals[idx][-1])

    def testMakeStateSpaceZeroInput(self):
        if IGNORE_TEST:
          return
        input_names = ["J0"]
        ctlsb = ControlBase(LINEAR_MDL,
              input_names=input_names, output_names=["S3", "S2"])
        sys = ctlsb.makeStateSpace(A_mat=ctlsb.jacobian_df.values)
        num_state, _ = np.shape(ctlsb.jacobian_df.values)
        self.assertEqual(np.shape(sys.B), (num_state, 1))
        self.assertEqual(np.shape(sys.C), (2, num_state))
        # Simulate the system
        END = 30
        times = [0*v for v in range(END)]
        X0 = ctlsb.state_ser.values
        U = np.repeat(2, END)
        times, y_val1s = control.forced_response(sys, T=times, X0=X0, U=U)
        times, y_val2s = control.forced_response(sys, T=times, X0=X0)
        for idx in range(2):
            self.assertEqual(np.sum(y_val1s[idx, :] - y_val2s[idx, :])**2, 0)

    def test_state_ser(self):
        if IGNORE_TEST:
          return
        def test(mdl):
            ctlsb = ControlBase(mdl)
            num_species = len(ctlsb.state_names)
            X0 = ctlsb.state_ser
            self.assertEqual(len(X0), num_species)
            jacobian_df = ctlsb.jacobian_df
            self.assertTrue(isinstance(jacobian_df, pd.DataFrame))
            self.assertEqual(len(jacobian_df), len(jacobian_df.columns))
            dX0 = np.matmul(ctlsb.jacobian_df.values, X0.values)
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
        jac_0 = self.ctlsb.jacobian_df
        self.ctlsb.setTime(5)
        jac_5 = self.ctlsb.jacobian_df
        isEqual(jac_0, jac_5, False)
        #
        self.ctlsb.setTime(0)
        jac_00 = self.ctlsb.jacobian_df
        isEqual(jac_00, jac_0, True)

    def testMakeBMatrix(self):
        if IGNORE_TEST:
          return
        ctlsb = ControlBase(LINEAR_MDL, input_names=LINEAR_MDL_INPUT_NAMES)
        B_df = ctlsb._makeBDF()
        self.assertEqual(np.shape(B_df.values), (3, 3))
        self.assertTrue(all(B_df.values.flatten()
              == ctlsb.full_stoichiometry_df.values.flatten()))
        #
        ctlsb = ControlBase(LINEAR_MDL, input_names=["J0", "J2"])
        B_df = ctlsb._makeBDF()
        self.assertEqual(np.shape(B_df.values), (3, 2))

    def testMakeUserError(self):
        if IGNORE_TEST:
          return
        with self.assertRaises(ValueError):
            _ = ControlBase(LINEAR_MDL, input_names=["J0", "K2"])
        with self.assertRaises(ValueError):
            _ = ControlBase(LINEAR_MDL, output_names=["S1", "SS2"])

    def testMakeNonlinearIOSystem(self):
        if IGNORE_TEST:
          return
        ctlsb = ControlBase(LINEAR_MDL, input_names=["J0"])
        effector_dct = {"J0": "S0"}
        non_sys = ctlsb.makeNonlinearIOSystem("tst", effector_dct=effector_dct)
        self.assertTrue("NonlinearIOSystem" in str(type(non_sys)))


if __name__ == '__main__':
  unittest.main()

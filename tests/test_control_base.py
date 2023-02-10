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
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
BIOMD1015 = os.path.join(TEST_DIR, "Jarrah2014.xml")
MODEL_FILE = os.path.join(TEST_DIR, "BIOMD0000000823.xml")
ANTIMONY_FILE = os.path.join(TEST_DIR, "Model_antimony.ant")
REACTION_NAMES = ["J1"]
OUTPUT_NAMES = ["S3", "S2"]
END_TIME = 20
LINEAR2_MDL = """
$S1 -> S2; k1*$S1
J1: S2 -> S3; k2*S2
J2: S3 -> S2; k3*S3
J3: S2 -> ; k4*S2

k1 = 1
k2 = 2
k3 = 3
k4 = 4
$S1 = 10
S2 = 0
S3 = 0
S4 = 0
"""
LINEAR3_MDL = """
S1 -> S2; k1*S1
J1: S1 -> S3; k2*S1
J2: S2 -> S3; k3*S2
J3: S3 -> ; k3*S3

k1 = 1
k2 = 2
k3 = 3
k4 = 4
S1 = 10
S2 = 0
S3 = 0
S4 = 0
"""
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
LINEAR_MDL_REACTION_NAMES = ["J0", "J1", "J2"]
LINEAR_MDL_SPECIES_NAMES = ["S1", "S2"]
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
        self.ctlsb = ControlBase(ANTIMONY_FILE, is_reduced=True)

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
        self.init()
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

    def testAdd(self):
        if IGNORE_TEST:
          return
        self.init()
        dct = {"S1": 1, "S2": 2}
        names = dct.keys()
        cur_dct = self.ctlsb.get(names)
        self.ctlsb.add(dct)
        new_dct = self.ctlsb.get(names)
        trues = [new_dct[n] == cur_dct[n] + dct[n] for n in names]
        self.assertTrue(all(trues))

    def testCopyEquals(self):
        if IGNORE_TEST:
            return
        self.init()
        ctlsb = self.ctlsb.copy()
        # FIXME: takes a long time
        self.assertTrue(ctlsb.equals(self.ctlsb, is_quick_check=True))
        ctlsb.antimony = ""
        self.assertFalse(ctlsb.equals(self.ctlsb, is_quick_check=True))

    def testMakeStateSpace1(self):
        if IGNORE_TEST:
          return
        self.init()
        sys = self.ctlsb.makeStateSpace(A_mat=self.ctlsb.jacobian_df.values)
        self.assertEqual(sys.nstates, 3)
        #
        sys1 = self.ctlsb.makeStateSpace(time=4)
        self.assertFalse(sys.A[0][0] == sys1.A[0][0])
        sys2 = self.ctlsb.makeStateSpace(time=0)
        self.assertTrue(np.isclose(sys.A[0][0], sys2.A[0][0]))

    def testMakeStateSpace2(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlBase(LINEAR_MDL, input_names=["S1"],
              output_names=["S2"])
        sys = ctlsb.makeStateSpace()
        self.assertEqual(sys.nstates, 3)
        times = [0.1*v for v in range(50)]
        X0 = ctlsb.state_ser.values
        _, y_vals = control.forced_response(sys, T=times, X0=X0)
        self.assertGreater(y_vals[-1], y_vals[0])

    def testMakeStateSpace3(self):
        if IGNORE_TEST:
          return
        ctlsb = ControlBase(BIOMD1015)
        sys = ctlsb.makeStateSpace()
        self.assertTrue("StateSpace" in str(type(sys)))

    def testMakeStateSpaceSpeciesInput(self):
        if IGNORE_TEST:
          return
        ctlsb = ControlBase(LINEAR2_MDL, input_names=["S2"],
              output_names=["S3"])
        sys = ctlsb.makeStateSpace()
        tf = control.ss2tf(sys)
        num = tf.num[0][0]
        self.assertTrue(np.isclose(num[0], 2))
        den = tf.den[0][0]
        self.assertEqual(den[0], 1)

    def testMakeStateSpaceReducable(self):
        if IGNORE_TEST:
          return
        ctlsb = ControlBase(REDUCABLE_MDL, output_names=OUTPUT_NAMES_REDUCABLE,
              is_reduced=True)
        sys = ctlsb.makeStateSpace()
        self.assertGreater(np.shape(sys.C)[0], np.shape(sys.A)[0])

    def testMakeStateSpaceReducableWithOutputFlux(self):
        if IGNORE_TEST:
            return
        def test(time, is_equal):
            ctlsb = ControlBase(REDUCABLE_MDL, output_names=["S1", "S2"],
                  is_reduced=True)
            ctlsb.setTime(time)
            sys = ctlsb.makeStateSpace()
            is_equal = all([x == y for x, y in zip(sys.C[0, :], sys.C[1,:])])
            self.assertEqual(is_equal, is_equal)
        #
        test(0, False)
        test(1, True)

    def _simulate(self, u_val, mdl=LINEAR_MDL, input_names=["S1"],
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
        self.init()
        input_names = ["S1"]
        ctlsb = ControlBase(LINEAR_MDL,
              input_names=input_names, output_names=["S1", "S2", "S3"])
        sys = ctlsb.makeStateSpace(A_mat=ctlsb.jacobian_df.values)
        num_state, _ = np.shape(ctlsb.jacobian_df.values)
        self.assertEqual(np.shape(sys.B), (num_state, 1))
        self.assertEqual(np.shape(sys.C), (3, num_state))
        # Simulate the system
        END = 50
        times = [0.1*v for v in range(END)]
        X0 = ctlsb.state_ser.values
        U = np.repeat(2, END)
        times, y_val1s = control.forced_response(sys, T=times, X0=X0, U=U)
        times, y_val2s = control.forced_response(sys, T=times, X0=X0)
        for time in times:
            total = sum([y_val2s[i, int(time)] for i in range(3)])
            self.assertTrue(np.isclose(total, X0[0]))
        #
        time_start = int(times[0])
        time_end = int(times[-1])
        total_end = sum([y_val1s[i, time_end] for i in range(3)])
        total_start = sum([y_val1s[i, time_start] for i in range(3)])
        self.assertGreater(total_end - total_start, 0.2)

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

    def testMakeBMatrixReactionInputs(self):
        if IGNORE_TEST:
            return
        with self.assertRaises(ValueError):
            # Reactions cannot be an input
            ctlsb = ControlBase(LINEAR_MDL, 
                 input_names=LINEAR_MDL_REACTION_NAMES)

    def testMakeBMatrixSpeciesInputs(self):
        if IGNORE_TEST:
            return
        for input_names in [ None, ["S2"], ["S1", "S2"]]:
            ctlsb = ControlBase(LINEAR_MDL, input_names=input_names)
            B_df = ctlsb._makeBDF()
            trues = [B_df[c].sum() == 1 for c in B_df.columns]
            self.assertTrue(all(trues))
            self.assertEqual(len(B_df), len(ctlsb.state_names))
        #
        ctlsb = ControlBase(LINEAR_MDL, input_names=LINEAR_MDL_SPECIES_NAMES)
        B_df = ctlsb._makeBDF()
        self.assertEqual(np.shape(B_df.values),
              (3, len(LINEAR_MDL_SPECIES_NAMES)))

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
        ctlsb = ControlBase(LINEAR_MDL, input_names=["S1"])
        non_sys = ctlsb.makeNonlinearIOSystem("tst")
        self.assertTrue("NonlinearIOSystem" in str(type(non_sys)))

    def testMakeTransferFunction(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlBase(LINEAR_MDL, input_names=["S3"], output_names=["S1"])
        tf = ctlsb.makeTransferFunction()
        self.assertTrue("TransferFunction" in str(type(tf)))
        with(self.assertRaises(ValueError)):
            ctlsb = ControlBase(LINEAR_MDL, output_names=["S1"])
            tf = ctlsb.makeTransferFunction()

    def testReduceTransferFunction(self):
        if IGNORE_TEST:
          return
        self.init()
        def test(num_numerator_zero, num_denominator_zero):
            numerator_ply = list(np.repeat(0, num_numerator_zero))
            numerator_ply.insert(0, 1)
            denominator_ply = list(np.repeat(0, num_denominator_zero))
            denominator_ply.insert(0, 1)
            denominator_ply.insert(0, 1)
            tf = control.TransferFunction(numerator_ply, denominator_ply)
            new_tf = self.ctlsb.reduceTransferFunction(tf)
            dcgain = new_tf.dcgain()
            if num_numerator_zero == num_denominator_zero:
                self.assertTrue(np.isclose(new_tf.dcgain(), 1.0))
            elif num_numerator_zero > num_denominator_zero:
                self.assertTrue(np.isclose(new_tf.dcgain(), 0))
            elif num_numerator_zero < num_denominator_zero:
                self.assertEqual(new_tf.dcgain(), np.inf)
        #
        test(0, 0)
        test(1, 1)
        test(0, 1)
        test(1, 2)
        test(1, 0)
        test(2, 1)

    def testMakeTransferFunction2(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlBase(LINEAR3_MDL,
            input_names=["S1"], output_names=["S3"])
        tf = ctlsb.makeTransferFunction(time=1, atol=1e-3)
        dcgain = tf.dcgain()
        self.assertLess(np.abs(dcgain -  0.33), 0.01)

    def testMakeFluxJacobian(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlBase(MODEL_FILE,
            input_names=["IR"], output_names=["mTORC1_DEPTOR"])
        df_0 = ctlsb.makeFluxJacobian(0)
        df_2 = ctlsb.makeFluxJacobian(2)
        df = df_0 - df_2
        df = df.applymap(lambda v: np.abs(v))
        max_val = df.max().max()
        self.assertGreater(max_val, 0.5)

    def testMakeSISOTransferFunctionBuilder(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlBase(MODEL_FILE,
            input_names=["IR"], output_names=["mTORC1_DEPTOR"])
        builder = ctlsb.makeSISOTransferFunctionBuilder()
        self.assertTrue("Builder" in str(type(builder)))
        #
        self.init()
        with self.assertRaises(ValueError):
            builder = self.ctlsb.makeSISOTransferFunctionBuilder()


if __name__ == '__main__':
  unittest.main()

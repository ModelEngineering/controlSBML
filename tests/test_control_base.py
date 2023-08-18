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
MODEL_FILE = os.path.join(TEST_DIR, "BIOMD0000000823.xml")
BIOMD1015 = os.path.join(TEST_DIR, "Jarrah2014.xml")
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

SPECIES_NAMES = ["S0", "S1", "S2", "SS0"]

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
        diff = set(SPECIES_NAMES).symmetric_difference(self.ctlsb.input_names)
        self.assertEqual(len(diff), 0)
        #
        diff = set(SPECIES_NAMES).difference(self.ctlsb.output_names)
        self.assertGreater(len(diff), 0)

    # FIXME: Delete test? 
    def deprecatedTestConstructor(self):
        if IGNORE_TEST:
          return
        self.assertTrue("RoadRunner" in str(type(self.ctlsb.roadrunner)))
        for lst in [self.ctlsb.species_names, self.ctlsb.output_names]:
            diff = set(SPECIES_NAMES).symmetric_difference(lst)
            self.assertEqual(len(diff), 0)

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

    def testSetSteadyState(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlBase(NONLINEAR1_MDL)
        ser1 = ctlsb.species_ser
        ctlsb.setSteadyState()
        ser2 = ctlsb.species_ser
        ctlsb.setSteadyState()
        ser3 = ctlsb.species_ser
        self.assertFalse(ser1.equals(ser2))
        self.assertFalse(ser2.equals(ser3))

    def test_species_ser(self):
        if IGNORE_TEST:
          return
        def test(mdl):
            ctlsb = ControlBase(mdl)
            num_species = len(ctlsb.species_names)
            X0 = ctlsb.species_ser
            self.assertEqual(len(X0), num_species)
            jacobian_df = ctlsb.jacobian_df
            self.assertTrue(isinstance(jacobian_df, pd.DataFrame))
            self.assertEqual(len(jacobian_df), len(jacobian_df.columns))
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

    def testMakeUserError(self):
        if IGNORE_TEST:
          return
        with self.assertRaises(ValueError):
            _ = ControlBase(LINEAR_MDL, input_names=["J0", "K2"])
        with self.assertRaises(ValueError):
            _ = ControlBase(LINEAR_MDL, output_names=["S1", "SS2"])


if __name__ == '__main__':
  unittest.main()

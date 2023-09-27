import controlSBML.constants as cn
from controlSBML.sbml_system import SBMLSystem
from controlSBML import util

import matplotlib.pyplot as plt
import numpy as np
import unittest


IGNORE_TEST = True
IS_PLOT = True
IMPROPER_LINEAR_MDL = """
// Illustrate Antimony File

aa := S1
bb := S4

S1 -> S2; k1*S1
J1: S2 -> S3; k2*S2
J2: S3 -> S2; k3*S3
J3: S2 -> S4; k4*S2

k1 = 1
k2 = 1
k3 = 1
k4 = 1
S1 = 10
S2 = 0
S3 = 0
"""
LINEAR_MDL = "model *main_model()\n" + IMPROPER_LINEAR_MDL + "\nend"
SPECIES_NAMES = ["S1", "S2", "S3"]


#############################
# Tests
#############################
class TestSBMLSystem(unittest.TestCase):

    def setUp(self):
        self.system = SBMLSystem(LINEAR_MDL, ["S1"], ["S3"], is_fixed_input_species=True)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.system.antimony, str))
        self.assertGreater(len(self.system.reaction_names), 0)
        self.assertGreater(len(self.system.floating_species_names), 0)

    def testImproperModel(self):
        if IGNORE_TEST:
            return
        with self.assertRaises(SystemExit):
            SBMLSystem(IMPROPER_LINEAR_MDL, ["S1"], ["S3"], is_fixed_input_species=True)
                          
    def testConstructor2(self):
        if IGNORE_TEST:
                return
        with self.assertRaises(ValueError):
            _ = SBMLSystem(LINEAR_MDL, ["S9"], ["S3"])
        with self.assertRaises(ValueError):
            _ = SBMLSystem(LINEAR_MDL, ["S1"], ["S9"])

    def testConstructor3(self):
        if IGNORE_TEST:
                return
        self.system = SBMLSystem(LINEAR_MDL, ["aa"], ["bb"], is_fixed_input_species=True)

    def testGet(self):
        if IGNORE_TEST:
            return
        for name in SPECIES_NAMES:
            self.assertEqual(self.system.get(name), self.system.roadrunner[name])

    def testSet(self):
        if IGNORE_TEST:
            return
        for name in SPECIES_NAMES:
            self.system.set({name: 100})
            self.assertEqual(self.system.get(name), 100)

    def testSimulate1(self):
        if IGNORE_TEST:
            return
        initial_s1 = self.system.get("S1")
        ts = self.system.simulate(0, 100)
        self.assertTrue(np.isclose(ts["S4"].values[-1], 10))
       
    def testSimulateSISOClosedLoop(self):
        if IGNORE_TEST:
            return
        def test(is_fixed_input_species):
            system = SBMLSystem(LINEAR_MDL, ["S1"], ["S3"], is_fixed_input_species=is_fixed_input_species)
            reference = 5
            ts = system.simulateSISOClosedLoop("S1", "S3", kp=2, ki=0.8, kf=0.5, reference=reference, end_time=20, num_point=1000)
            self.assertGreater(len(ts), 0)
            if is_fixed_input_species:
                tolerance = 0.3
            else:
                tolerance = 0.8
            self.assertLess(np.abs(ts["S3"].values[-1] - reference), tolerance)
            if IS_PLOT:
                for column in ts.columns:
                    if column not in ["time", "S1", "S3"]:
                        del ts[column]
                if is_fixed_input_species:
                    title = "Boundary species"
                else:
                    title = "Boundary reaction"
                util.plotOneTS(ts, ax2=0, figsize=(8,8), title=title)
                plt.show()
        #
        test(False)
        test(True)

    def testSimulateStaircase(self):
        if IGNORE_TEST:
            return
        times = np.linspace(0, 50, 500)
        def test(is_fixed_input_species):
            system = SBMLSystem(LINEAR_MDL, ["S1"], ["S3"], is_fixed_input_species=True)
            ts = system.simulateStaircase("S1", "S3", times=times, final_value=10, num_step=5, is_steady_state=False)
            self.assertGreater(len(ts), 0)
            variance = np.var(ts["S3"])
            self.assertFalse(np.isclose(variance, 0))
            if IS_PLOT:
                util.plotOneTS(ts, ax2=0, figsize=(8,8), ylim=[0, 10])
                plt.show()
        #
        test(True)
        test(False)

    # Test staircase on a closedloop system


if __name__ == '__main__':
  unittest.main()
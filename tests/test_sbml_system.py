import controlSBML.constants as cn
from controlSBML.sbml_system import SBMLSystem
from controlSBML import util

import matplotlib.pyplot as plt
import numpy as np
import unittest


IGNORE_TEST = False
IS_PLOT = False
LINEAR_MDL = """
// Illustrate Antimony File

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

    def testConstructor2(self):
       if IGNORE_TEST:
            return
       with self.assertRaises(ValueError):
           _ = SBMLSystem(LINEAR_MDL, ["S9"], ["S3"])
       with self.assertRaises(ValueError):
           _ = SBMLSystem(LINEAR_MDL, ["S1"], ["S9"])

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
        ts = self.system.simulateSISOClosedLoop("S1", "S3", kp=60, ki=1, kf=None, reference=5, end_time=20, num_point=1000)
        self.assertGreater(len(ts), 0)
        if IS_PLOT:
            del ts["integral_control_error_S1_S3"]
            if "filter_S1_S3" in ts.columns:
                del ts["filter_S1_S3"]
            del ts["S4"]
            util.plotOneTS(ts, ax2=0, figsize=(8,8))
            plt.show()

    def testSimulateStaircase(self):
        if IGNORE_TEST:
            return
        times = np.linspace(0, 50, 500)
        system = SBMLSystem(LINEAR_MDL, ["S1"], ["S3"], is_fixed_input_species=True)
        ts = system.simulateStaircase("S1", "S3", times=times, final_value=10, num_step=5, is_steady_state=False)
        self.assertGreater(len(ts), 0)
        if IS_PLOT:
            util.plotOneTS(ts, ax2=0, figsize=(8,8), ylim=[0, 10])
            plt.show()


if __name__ == '__main__':
  unittest.main()
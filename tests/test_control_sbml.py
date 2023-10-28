from controlSBML.control_sbml import ControlSBML

import numpy as np
import unittest


IGNORE_TEST = True
IS_PLOT = False
LINEAR_MDL = """
model *main_model()
// Illustrate Antimony File
species S1, S2, S3, S4

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
end
"""
SPECIES_NAMES = ["S1", "S2", "S3", "S4"]


#############################
# Tests
#############################
class TestControlSBML(unittest.TestCase):

    def setUp(self):
        self.ctlsb = ControlSBML(LINEAR_MDL)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue("RoadRunner" in str(type(self.ctlsb.roadrunner)))

    def testSetTimes(self):
        if IGNORE_TEST:
            return
        max_value = 10
        self.ctlsb.setTimes(np.linspace(0, max_value, 100))
        times = self.ctlsb.getTimes()
        self.assertEqual(times[-1], max_value)

    def testSetSBMLSystem(self):
        #if IGNORE_TEST:
        #    return
        self.ctlsb.setSystem(input_names=["S1"], output_names=["S2"], is_fixed_input_species=False, is_steady_state=False)
        system = self.ctlsb.getSBMLSystem()
        self.assertTrue("S1" in system.input_names)


if __name__ == '__main__':
  unittest.main()
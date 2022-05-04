import controlSBML.constants as cn
import controlSBML as ctl
from controlSBML.siso_closed_loop_system import SISOClosedLoopSystem

import control
import pandas as pd
import numpy as np
import tellurium as te
import unittest


IGNORE_TEST = False
IS_PLOT = False
SIZE = 10
if IS_PLOT:
    import matplotlib
    matplotlib.use('TkAgg')
TIMES = ctl.makeSimulationTimes()
MODEL = """
S0 -> 2 S0; S0
S0 -> S1; S0
S1 -> S2; S1*S1
S2 -> S3; S2*S1

S0 = 1
S1 = 10
S2 = 0
S3 = 0
"""
        

#############################
# Tests
#############################
class TestSISOClosedLoopSystem(unittest.TestCase):

    def setUp(self):
        self.ctlsb = ctl.ControlSBML(MODEL, input_names=["S0"],
              output_names=["S3"])
        self.siso = SISOClosedLoopSystem(self.ctlsb)

    def testConstructor(self):
        self.assertTrue("controller" in self.siso.component_names)


if __name__ == '__main__':
  unittest.main()

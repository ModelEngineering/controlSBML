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
HTTP_FILE = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000206.2?filename=BIOMD0000000206_url.xml"
HTTP_FILE = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000823.2?filename=Varusai2018.xml"
        
if IS_PLOT:
    import matplotlib
    matplotlib.use('TkAgg')
TIMES = ctl.makeSimulationTimes()
MODEL1 = """
S0 -> S1; S0
S0a -> S1; S0a
S1 -> S2; S1*S1
S2 -> S3; S2*S1
S3 ->; S3

S0a = 1
S0 = 1
S1 = 10
S2 = 0
S3 = 0
"""
MODEL2 = """
S0 -> S1; S0
S2 -> S3; S2

S0 = 10
S1 = 0
S2 = 10
S3 = 0
"""
        

#############################
# Tests
#############################
class TestSISOClosedLoopSystem(unittest.TestCase):

    def setUp(self):
        self.ctlsb = ctl.ControlSBML(MODEL1, input_names=["S0", "S2"],
              output_names=["S1", "S3"])
        self.siso = SISOClosedLoopSystem(self.ctlsb)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue("controller" in self.siso.component_names)

    def testEvaluateControllability1(self):
        if IGNORE_TEST:
          return
        ctlsb = ctl.ControlSBML(MODEL2, input_names=["S0", "S2"],
              output_names=["S1", "S3"])
        siso = SISOClosedLoopSystem(ctlsb)
        times = [0, 1]
        dct = siso.evaluateControllability(times)
        self.assertEqual(len(times), len(dct))
        df = dct[times[0]]
        self.assertTrue(df.loc["S0", :].equals(df.loc["S2", :]))

    def testEvaluateControllability2(self):
        if IGNORE_TEST:
          return
        times = [0, 1, 2, 3]
        dct = self.siso.evaluateControllability(times)
        self.assertEqual(len(dct), len(times))

    def testEvaluateControllability3(self):
        if IGNORE_TEST:
          return
        ctlsb = ctl.ControlSBML(HTTP_FILE)
        state_names = ctlsb.state_names
        ctlsb = ctl.ControlSBML(HTTP_FILE, input_names=state_names[0:3],
              output_names=state_names[-3:])
        siso = SISOClosedLoopSystem(ctlsb)
        times = [0, 1, 2, 3]
        dct = siso.evaluateControllability(times)
        self.assertEqual(len(times), len(dct))


if __name__ == '__main__':
  unittest.main()

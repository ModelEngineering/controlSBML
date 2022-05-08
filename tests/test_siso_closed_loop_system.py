import controlSBML.constants as cn
import controlSBML as ctl
from controlSBML.siso_closed_loop_system import SISOClosedLoopSystem

import control
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import tellurium as te
import unittest


IGNORE_TEST = False
IS_PLOT = False
SIZE = 10
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
BIOMD823 = os.path.join(TEST_DIR, "BIOMD0000000823.xml")
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
MODEL3 = """
S0 -> S1; S0
S1 -> S2; S1

S0 = 10
S1 = 0
S2 = 0
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
        ctlsb = ctl.ControlSBML(BIOMD823)
        state_names = ctlsb.state_names
        ctlsb = ctl.ControlSBML(BIOMD823, input_names=state_names[0:3],
              output_names=state_names[-3:])
        siso = SISOClosedLoopSystem(ctlsb)
        times = [0, 1, 2, 3]
        dct = siso.evaluateControllability(times)
        self.assertEqual(len(times), len(dct))

    def testMakeClosedLoopSystem1(self):
        if IGNORE_TEST:
          return
        ctlsb = ctl.ControlSBML(MODEL2, input_names=["S2"],
              output_names=["S3"])
        siso = SISOClosedLoopSystem(ctlsb)
        siso.makeClosedLoopSystem("cl_sys", kp=20)
        times = ctl.makeSimulationTimes(end_time=100)
        result = control.input_output_response(siso.closed_loop_system,
              times, U=1)
        if IS_PLOT:
            plt.plot(result.t.flatten(), result.y.flatten())
            plt.show()
        outputs = result.y.flatten()
        self.assertGreater(outputs[-1], 0.7)

    def testMakeClosedLoopSystem2(self):
        if IGNORE_TEST:
          return
        ctlsb = ctl.ControlSBML(BIOMD823)
        state_names = ctlsb.state_names
        ctlsb = ctl.ControlSBML(BIOMD823, input_names=["pAkt"],
              output_names=["pDEPTOR"])
        siso = SISOClosedLoopSystem(ctlsb)
        closed_loop_outputs=["sum_R_F.in2", "sum_R_F.out", 
              "sum_U_D.out", "sum_Y_N.out"]
        siso.makeClosedLoopSystem("cl_sys", kp=1000, ki=2, kf=None,
            noise_amp=0, noise_frq=20,
            closed_loop_outputs=closed_loop_outputs)
        ts = siso.makeStepResponse(time=1, step_size=1, end_time=300,
              points_per_time=2)
        ts.columns = ["input", "e(t)", "system.in", "output"]
        if IS_PLOT:
            ctl.plotOneTS(ts, xlabel="time", ylim=[-5, 5])
        self.assertGreater(len(ts), 0)
        df = siso.factory.report()
        if IS_PLOT:
            plt.plot(df.index, df["sum_R_F.out"])
            plt.xlabel("time")
            plt.ylabel("e(t)")
            plt.show()
        # Do a second simulation to check the log
        ts = siso.makeStepResponse(time=1, step_size=1, end_time=300,
              points_per_time=2)
        ts.columns = ["input", "e(t)", "system.in", "output"]
        if IS_PLOT:
            ctl.plotOneTS(ts, xlabel="time", ylim=[-5, 5])
        df2 = siso.report()
        self.assertTrue(df.equals(df2))


if __name__ == '__main__':
  unittest.main()

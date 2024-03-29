import controlSBML.constants as cn
import controlSBML as ctl
from controlSBML import siso_closed_loop_system
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
MODEL3a = """

-> S0; k0
S0 -> S1; k1*S0
S1 -> S2; k2*S1
S2 -> S1; k3*S2
S2 -> ; k4*S2

S0 = 0
S1 = 0
S2 = 0
S3 = 0
k0 = 0.5
k1 = 1
k2 = 2
k3 = 3
k4 = 4
"""
BIOMD1015 = os.path.join(TEST_DIR, "Jarrah2014.xml")
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
S1 -> S3; S1

S0 = 0
S1 = 0
S2 = 0
S3 = 0
"""
MODEL2A = """
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
        self.assertTrue("Factory" in str(type(self.siso.factory)))

    def testMakeClosedLoopSystem1(self):
        if IGNORE_TEST:
            return
        ctlsb = ctl.ControlSBML(MODEL3, input_names=["S0"],
              output_names=["S2"])
        siso = SISOClosedLoopSystem(ctlsb, is_fixed_input_species=True)
        siso.makePIDClosedLoopSystem(kp=0.2, closed_loop_outputs=["S1", "S2"])
        times = ctl.makeSimulationTimes(end_time=100)
        inputs = np.repeat(0, len(times))
        inputs[0] = 1
        result = control.input_output_response(siso.closed_loop_system,
              times, U=inputs)
        if IS_PLOT:
            plt.plot(result.t.flatten(), result.y[0].flatten())
            plt.plot(result.t.flatten(), result.y[1].flatten())
            plt.show()
        self.assertEqual(np.shape(result.y)[0], 2)  # Should have 2 outputs
        max_val = np.max(result.y[1,:])
        max_ratio = result.y[1,-1]/max_val
        self.assertLess(max_ratio, 0.01)  # S2 should go to 0

    def testMakePIDClosedLoopSystem2(self):
        # FIXME: very inefficient test
        if True:
          return
        if IGNORE_TEST:
          return
        ctlsb = ctl.ControlSBML(BIOMD823, input_names=["pAkt"],
              output_names=["pDEPTOR"])
        siso = SISOClosedLoopSystem(ctlsb)
        closed_loop_outputs=["cl_input.out", "sum_F_R.out", 
              "sum_D_U.out", "cl_output.out"]
        siso.makePIDClosedLoopSystem(kp=10, ki=0, kf=None,
            noise_amp=0, noise_frq=20,
            closed_loop_outputs=closed_loop_outputs)
        ts = siso.makeStepResponse(time=1, step_size=1, end_time=300,
              points_per_time=2)
        if IS_PLOT:
            ctl.plotOneTS(ts, xlabel="time", ylim=[-5, 5])
        self.assertGreater(len(ts), 0)
        df = siso.factory.report()
        self.assertGreater(len(df), 0)
        if IS_PLOT:
            plt.plot(df.index, df["sum_F_R.out"])
            plt.xlabel("time")
            plt.ylabel("e(t)")
            plt.show()

    def testMakeFullStatelosedLoopSystem2(self):
        if IGNORE_TEST:
            return
        ctlsb = ctl.ControlSBML(BIOMD823, input_names=["pAkt"],
              output_names=["pDEPTOR"])
        siso = SISOClosedLoopSystem(ctlsb)
        closed_loop_outputs=["sum_F_R.out", 
              "sum_D_U.out", "cl_output.out"]
        siso.makePIDClosedLoopSystem(kp=1, ki=0, kf=None,
            noise_amp=0, noise_frq=20,
            closed_loop_outputs=closed_loop_outputs)
        ts = siso.makeResponse(time=0, step_size=1, end_time=5, period=2,
              points_per_time=10)
        ts.columns = ["input", "e(t)", "system.in", "output"]
        self.assertGreater(len(ts), 0)
        self.assertEqual(sum(sum(np.isnan(ts.to_numpy()))), 0)

    # FIXME: makeFullStateClosedLoopSystem is broken
    def testMakeFullState(self):
        if IGNORE_TEST:
          return
        return
        ctlsb = ctl.ControlSBML(MODEL2, input_names=["S0"],
              output_names=["S3"])
        siso = SISOClosedLoopSystem(ctlsb)
        siso.makeFullStateClosedLoopSystem(poles=-2)
        times = ctl.makeSimulationTimes(end_time=10)
        ts = siso.makeStepResponse(step_size=2, end_time=10)
        if IS_PLOT:
            ctl.plotOneTS(ts, figsize=(5,5))
            plt.show()
        self.assertGreater(np.abs(ts["S3"].values[-1]), 0)

    # TODO: More tests of fullstate filters
    def testMakeFullStateControllerWithFilters(self):
        # TODO: make test more efficient
        if True:
            return
        if IGNORE_TEST:
          return
        ctlsb = ctl.ControlSBML(MODEL2, input_names=["S0"],
              output_names=["S3"])
        siso = SISOClosedLoopSystem(ctlsb)
        siso.makeFullStateClosedLoopSystem(poles=-2, 
              kf=-1, noise_amp=0.1, noise_frq=20)
        ts = siso.makeStepResponse(step_size=2, end_time=200)
        siso1 = SISOClosedLoopSystem(ctlsb)
        siso1.makeFullStateClosedLoopSystem(poles=-2, 
              kf=0, noise_amp=0.1, noise_frq=20)
        ts1 = siso1.makeStepResponse(step_size=2, end_time=200)
        diff = np.sum((ts1.values.flatten() - ts.values.flatten())**2)
        self.assertGreater(diff, 0)
        if IS_PLOT:
            ctl.plotOneTS(ts, figsize=(5,5), ylim=[0, 10])
            plt.show()
            ctl.plotOneTS(ts1, figsize=(5,5), ylim=[0, 10])
            plt.show()
        self.assertGreater(np.abs(ts["S3"].values[-1]), 0)
        df = siso.factory.report()
        trues = ["fltr.%s_in" % n in df.columns for n in ["S1", "S3"]]
        self.assertTrue(all(trues))

    def testMakePerturbation(self):
        if IGNORE_TEST:
          return
        PERTURB = "perturb"
        SUM = "sum"
        SUM_INPUT_NAMES = ["left", "right"]
        AMP = 1
        FRQ = 5
        dct = self.siso._makePerturbation(PERTURB, SUM, SUM_INPUT_NAMES,
              amp=AMP, frq=FRQ)
        trues = [k in [PERTURB, SUM, siso_closed_loop_system.CONNECTIONS]
              for k in dct.keys()]
        self.assertTrue(all(trues))

    def testMakeDisturbanceNoiseCLinputoutput(self):
        if IGNORE_TEST:
          return
        assembly = self.siso._makeDisturbanceNoiseCLinputoutput()
        def test(attribute, data_type, sub_data_type=None):
            value = assembly.__getattribute__(attribute)
            self.assertTrue(isinstance(value, data_type))
            if sub_data_type is not None:
                sub_value = value[0]
                self.assertTrue(isinstance(sub_value, sub_data_type))
        #
        test("con", list, sub_data_type=list)
        test("sys", list, sub_data_type=control.iosys.NonlinearIOSystem)
        test("inp", list, sub_data_type=str)
        test("out", list, sub_data_type=str)


if __name__ == '__main__':
  unittest.main()

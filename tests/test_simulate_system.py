import controlSBML as ctl
import controlSBML.constants as cn
from controlSBML import util

import control
import control
import os
import pandas as pd
import numpy as np
import unittest


IGNORE_TEST = False
IS_PLOT = False
if IS_PLOT:
    import matplotlib
    matplotlib.use('TkAgg')
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_FILE = os.path.join(TEST_DIR, "BIOMD0000000823.xml")
NONLINEAR_MDL = """
species E1;
J0:  -> S1; E1
J1: S1 -> S2; S1*S1
J2: S2 -> S3; S2*S1

E1 = 10
S1 = 1
S2 = 2
S3 = 3
"""
INITIAL_VALUES = [10, 1, 2, 3]
INPUT_NAMES = ["E1"]
OUTPUT_NAMES = ["S1", "S2"]

##### HELPERS ###
# Defined only in terms of control system
xeq = [4]
def outfcn(t, x, u, _):
    # State is accumulated error
    return xeq - u
CONTROLLER = control.NonlinearIOSystem(
  None,
  outfcn,
  inputs=['in'],
  outputs=['out'], name='controller')
# Control sbml
ctlsb = ctl.ControlSBML(NONLINEAR_MDL,
      input_names=INPUT_NAMES, output_names=OUTPUT_NAMES)
CTL_SYS = ctl.NonlinearIOSystem("CTL_SYS", ctlsb)
# Interconnect
CLOSED_OUTPUTS = list(CTL_SYS.outlist)
CLOSED_OUTPUTS.append("controller.out")
INTERCONNECT = control.interconnect(
  [CTL_SYS, CONTROLLER],       # systems
  connections=[
    ['CTL_SYS.E1', 'controller.out'],
    ['controller.in',  'CTL_SYS.S1'],
  ],
  inplist=["controller.in"],
  outlist=CLOSED_OUTPUTS,
)
MTOR_OUTPUT_NAMES = ["mTORC1_DEPTOR", "mTORC2_DEPTOR"]


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

    def setUp(self):
        self.times = util.makeSimulationTimes()
        self.ctl_sys = CTL_SYS
        self.sys = CONTROLLER
        self.interconnect = INTERCONNECT

    def testMakeStateVectorCtl(self):
        if IGNORE_TEST:
            return
        x_vec = ctl.simulate_system.makeStateVector(self.ctl_sys)
        trues = [np.isclose(x, y) for x, y in zip(x_vec, INITIAL_VALUES)]
        self.assertTrue(all(trues))

    def testMakeStateVectorControlNonlinearSystem(self):
        if IGNORE_TEST:
          return
        x_vec = ctl.simulate_system.makeStateVector(self.sys)
        ssq = np.sum(x_vec)**2
        self.assertTrue(np.isclose(ssq, 0))

    def testMakeStateVectorInterconnectedSystem(self):
        if IGNORE_TEST:
          return
        x_vec = ctl.simulate_system.makeStateVector(self.interconnect)
        self.assertEqual(len(x_vec), 4)

    def testSimulateSystemCtlNonlinearIOSystem(self):
        if IGNORE_TEST:
          return
        ts = ctl.simulateSystem(self.ctl_sys)
        self.assertTrue("Timeseries" in str(type(ts)))

    def testSimulateSystemControlNonlinearIOSystem(self):
        if IGNORE_TEST:
          return
        ts = ctl.simulateSystem(self.sys, output_names=["out"])
        self.assertTrue("Timeseries" in str(type(ts)))
        #
        ts = ctl.simulateSystem(self.sys)
        self.assertTrue("Timeseries" in str(type(ts)))

    def testSimulateSystemControlInterconnectedSystem(self):
        if IGNORE_TEST:
          return
        ts = ctl.simulateSystem(self.interconnect)
        self.assertTrue(all([isinstance(c, int) for c in ts.columns]))
        self.assertTrue("Timeseries" in str(type(ts)))
        #
        ts = ctl.simulateSystem(self.interconnect, output_names=CLOSED_OUTPUTS)
        self.assertTrue("Timeseries" in str(type(ts)))
        util.plotOneTS(ts, figsize=(5,5), title="InterconnectedSystem",
              is_plot=IS_PLOT)


if __name__ == '__main__':
  unittest.main()

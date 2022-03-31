import controlSBML as ctl 
import controlSBML.constants as cn
from controlSBML import util

import control
import pandas as pd
import numpy as np
import unittest


IGNORE_TEST = False
IS_PLOT = False
if IS_PLOT:
    import matplotlib
    matplotlib.use('TkAgg')
NONLINEAR_MDL = """
J0:  -> S1; $S0
J1: S1 -> S2; S1*S1
J2: S2 -> S3; S2*S1

$S0 = 10
S1 = 1
S2 = 2
S3 = 3
"""
INITIAL_VALUES = [1, 2, 3]
INPUT_NAMES = ["J0"]
OUTPUT_NAMES = ["S1", "S2"]

##### HELPERS ###
# Defined only in terms of control system
xeq = [4]
def updfcn(t, x, u, _):
    # x[0] - current error
    # x[1] - integral of error
    err, acc_err = x
    new_err = xeq - u
    acc_err += new_err
    return new_err, acc_err
def outfcn(t, x, u, _):
    # State is accumulated error
    control_error = x[0] + x[1] 
    control_error = min(control_error, 25)
    return control_error
CONTROLLER = control.NonlinearIOSystem(
  updfcn,
  outfcn,
  states=["error", "accumulated_error"],
  inputs=['in'],
  outputs=['out'], name='controller')
# Control sbml
ctlsb = ctl.ControlSBML(NONLINEAR_MDL,
      input_names=INPUT_NAMES, output_names=OUTPUT_NAMES)
CTL_SYS = ctl.NonlinearIOSystem("CTL_SYS", ctlsb, effector_dct={"J0": "S0"})
# Interconnect
CLOSED_OUTPUTS = list(CTL_SYS.outlist)
CLOSED_OUTPUTS.append("controller.out")
INTERCONNECT = control.interconnect(
  [CTL_SYS, CONTROLLER],       # systems
  connections=[
    ['CTL_SYS.J0', 'controller.out'],
    ['controller.in',  'CTL_SYS.S1'],
  ],
  inplist=["controller.in"],
  outlist=CLOSED_OUTPUTS,
)


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
        self.assertEqual(len(x_vec), 5)

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
        ts = ctl.simulateSystem(self.interconnect, output_names=CLOSED_OUTPUTS)
        self.assertTrue("Timeseries" in str(type(ts)))
        util.plotOneTS(ts, figsize=(5,5), title="InterconnectedSystem",
              is_plot=IS_PLOT)
        #
        ts = ctl.simulateSystem(self.interconnect)
        self.assertTrue(all([isinstance(c, int) for c in ts.columns]))
        self.assertTrue("Timeseries" in str(type(ts)))


if __name__ == '__main__':
  unittest.main()

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

    def makeMtor(self, timepoint=0, input_names=None, output_names=None):
        # Creates a NonlinearIOSystem named "mtro"
        if output_names is None:
            output_names = MTOR_OUTPUT_NAMES
        if input_names is None:
            input_names=["v1", "v11"]
        ctlsb = ctl.ControlSBML(
            "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000823.2?filename=Varusai2018.xml",
            input_names=input_names,
            output_names=output_names,
            is_reduced=True)
        ctlsb.setTime(timepoint)
        return ctlsb.makeNonlinearIOSystem("mtor")

    def testMtorBug1(self):
        if IGNORE_TEST:
          return
        non_linear_mtor = self.makeMtor()
        ts = ctl.simulateSystem(non_linear_mtor)
        self.assertEqual(len(ts.columns), len(MTOR_OUTPUT_NAMES))

    def testMtorBug2(self):
        if IGNORE_TEST:
          return
        mtor = self.makeMtor(timepoint=0)
        input_names = ["v6"]
        output_names = ["mTORC1_DEPTOR", "pAkt"]
        mtor = self.makeMtor(input_names=input_names, output_names=output_names)
        xeq = [100]  # Desired concentration for mTORC1_DEPTOR
        def outfcn(t, x, u, _):
            # State is accumulated error
            new_err = xeq[0] - u[0]
            return -30*new_err
        controller = control.NonlinearIOSystem(
          None,
          outfcn,
          inputs=['in'],
          outputs=['out'], name='controller')
        # Create the closed loop system
        closed_outputs = list(mtor.outlist)
        closed_outputs.append('controller.out')
        #closed_outputs.append("controller.out")  # Make this visible as well
        mtor_closed = control.interconnect(
          [mtor, controller],       # systems
          connections=[
            ['mtor.v6', 'controller.out'],
            ['controller.in',  'mtor.mTORC1_DEPTOR'],
          ],
          inplist=["controller.in"],
          outlist=closed_outputs,
        )
        initial_x_vec = ctl.makeStateVector(mtor_closed)
        ts = ctl.simulateSystem(mtor_closed,
              output_names=closed_outputs,
              initial_x_vec=initial_x_vec, end_time=200)
        self.assertGreater(len(ts), 0)


if __name__ == '__main__':
  unittest.main()

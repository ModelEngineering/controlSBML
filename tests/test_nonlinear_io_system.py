import controlSBML as ctl

import control
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import unittest
import tellurium as te


IGNORE_TEST = False
IS_PLOT = False
if IS_PLOT:
    import matplotlib
    matplotlib.use('TkAgg')
EFFECTOR_DCT = {"J0": "E_J0"}
END_TIME = 5
DT = 0.01
POINTS_PER_TIME = int(1.0 / DT)
NUM_TIME = int(POINTS_PER_TIME*END_TIME) + 1
TIMES = [n*DT for n in range(0, NUM_TIME)]

LINEAR_MDL = """
J0: $S0 -> S1; $S0
J1: S1 -> S2; S1
J2: S2 -> S3; S2

$S0 = 0
S1 = 10
S2 = 0
S3 = 0
"""
NONLINEAR_MDL = """
#J0: ->  S0; 1 + E_J0
J0: S0 -> S1; S0 + E_J0
J2: S1 -> S2; S1*S1
J3: S2 -> $S3; S2*S1

S0 = 10
S1 = 0
S2 = 0
E_J0 = 0
$S3 = 3
"""


#############################
# Tests
#############################
class TestNonlinearIOSystem(unittest.TestCase):

    def setUp(self):
        if IGNORE_TEST:
            return
        self.init()

    def init(self):
        self.ctlsb = ctl.ControlSBML(NONLINEAR_MDL,
              input_names=["J0"], output_names=["S1", "S2"])
        self.sys = ctl.NonlinearIOSystem("test_sys", self.ctlsb,
               effector_dct=EFFECTOR_DCT)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.ctlsb, ctl.ControlSBML))

    def runSimulation(self, states=None):
        u_vec = np.repeat(0.0, len(self.ctlsb.input_names))
        self.ctlsb.setTime(TIMES[0])
        x_vec = self.ctlsb.state_ser.values
        dct = {n: [] for n in self.ctlsb.state_names}
        for time in TIMES:
            dx_vec = self.sys.updfcn(time, x_vec, u_vec, {})
            new_x_vec = x_vec + DT*dx_vec
            [dct[n].append(v) for n, v in
                  zip(self.ctlsb.state_names, new_x_vec)]
            x_vec = new_x_vec
        if states is None:
           newDct = dct
        else:
           newDct = {k: v for k, v in dct.items() if k in states}
        df = pd.DataFrame(dct)
        return df

    def testCtlsbUpdfcn(self):
        if IGNORE_TEST:
          return
        self.init()
        df = self.runSimulation()
        self._checkWithSimulation(df)

    def testOutfcn(self):
        if IGNORE_TEST:
          return
        self.init()
        time = 3
        self.ctlsb.setTime(time)
        x_vec = self.ctlsb.state_ser.values
        u_vec = np.repeat(0, len(self.ctlsb.input_names))
        out_vec = self.sys.outfcn(time, x_vec, u_vec, None)
        self.assertEqual(len(out_vec), len(self.ctlsb.output_names))

    def testOutlist(self):
        if IGNORE_TEST:
          return
        self.init()
        names = self.sys.outlist
        self.assertTrue(all([self.sys.name in n for n in names]))
        self.assertEqual(len(names), self.ctlsb.num_output)

    def testGetStateSer(self):
        if IGNORE_TEST:
          return
        self.init()
        ser_0 = self.sys.getStateSer()
        ser_1 = self.sys.getStateSer(time=1)
        self.assertGreater(ser_0.loc["S0"], ser_1.loc["S1"])

    def testInplist(self):
        if IGNORE_TEST:
          return
        self.init()
        names = self.sys.inplist
        self.assertTrue(all([self.sys.name in n for n in names]))
        self.assertEqual(len(names), self.ctlsb.num_input)

    def _checkWithSimulation(self, df):
        df_rr = self.ctlsb.simulateRoadrunner(start_time=0, end_time=END_TIME,
               points_per_time=POINTS_PER_TIME)     
        df_rr = df_rr[df.columns]
        ssq = 0.0
        for column in df.columns:
            ssq += np.sum((df[column].values - df_rr[column].values)**2)
        # Calculate error per value
        rms_per_value = np.sqrt(ssq/(len(df.columns)*len(df)))
        self.assertLess(rms_per_value, 1e-1)

    def runInputOutputResponse(self, u_val, end_time=END_TIME):
        u_val = np.array(u_val)
        #
        x_vec = self.ctlsb.state_ser.values
        u_vec = np.repeat(u_val, len(self.ctlsb.input_names))
        u_vecs = np.array([np.array(u_vec) for _ in range(NUM_TIME)])
        u_vecs = np.reshape(u_vecs, (1, NUM_TIME))
        t, y = control.input_output_response(self.sys, TIMES, u_vecs, x_vec)
        df = pd.DataFrame(y, self.ctlsb.output_names)
        return df.transpose()

    def testWithInputOutputResponse(self):
        if IGNORE_TEST:
          return
        self.init()
        df = self.runInputOutputResponse(0)
        self._checkWithSimulation(df)
        # Output increases with inpu
        df1 = self.runInputOutputResponse(1)
        df2 = self.runInputOutputResponse(10)
        diff = (df2 - df1).sum().sum()
        self.assertGreater(diff, 10)
        


if __name__ == '__main__':
  unittest.main()

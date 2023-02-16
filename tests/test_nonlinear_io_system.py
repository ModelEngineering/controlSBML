import controlSBML as ctl
import controlSBML.constants as cn

import control
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import unittest
import shutil
import tellurium as te
import tempfile


IGNORE_TEST = False
IS_PLOT = False
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
species E_J0;
#J0: ->  S0; 1 + E_J0
J0: S0 -> S1; k1*S0*E_J0
J2: S1 -> S2; k2*S1*S1
J3: S2 -> $S3; k3*S2*S1

S0 = 10
S1 = 0
S2 = 0
E_J0 = 1
$S3 = 3
k1 =1; k2=2; k3=3
"""

# Temporarily change the plot path
if IS_PLOT:
    cn.PLOT_DIR = cn.TEST_DIR
else:
    cn.PLOT_DIR= tempfile.mkdtemp()


#############################
# Tests
#############################
class TestNonlinearIOSystem(unittest.TestCase):

    def setUp(self):
        if IGNORE_TEST:
            return
        self.init()
        self.removeFiles()

    def tearDown(self):
        plt.close()
        self.removeFiles()

    def init(self, do_simulate_on_update=True):
        if IS_PLOT:
            cn.PLOT_DIR = cn.TEST_DIR
        else:
            cn.PLOT_DIR= tempfile.mkdtemp()
        self.ctlsb = ctl.ControlSBML(NONLINEAR_MDL,
              input_names=["E_J0"], output_names=["S1", "S2"])
        self.sys = ctl.NonlinearIOSystem("test_sys", self.ctlsb,
               do_simulate_on_update=do_simulate_on_update)

    def removeFiles(self):
        for ffile in os.listdir(cn.PLOT_DIR):
            if ("figure_" in ffile) and (".pdf") in ffile:
                path = os.path.join(cn.PLOT_DIR, ffile)
                if os.path.isfile(path) and IGNORE_TEST:
                    os.remove(path)
        # FIXME: this won't work in Windows
        if IS_PLOT and ("var" in cn.PLOT_DIR):
            shutil.rmtree(cn.PLOT_DIR)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.init()
        self.assertTrue(isinstance(self.ctlsb, ctl.ControlSBML))

    def runSimulation(self, states=None):
        u_vec = np.repeat(0.0, len(self.ctlsb.input_names))
        self.sys.setTime(TIMES[0])
        x_vec = self.sys.makeStateSer().values
        dct = {n: [] for n in self.sys.state_names}
        for idx, time in enumerate(TIMES):
            dx_vec = self.sys.updfcn(time, x_vec, u_vec, {})
            new_x_vec = x_vec + DT*dx_vec
            [dct[n].append(v) for n, v in
                  zip(self.sys.state_names, new_x_vec)]
            x_vec = new_x_vec
        if states is None:
           newDct = dct
        else:
           newDct = {k: v for k, v in dct.items() if k in states}
        df = pd.DataFrame(dct)
        return df

    def testUpdfcn(self):
        if IGNORE_TEST:
          return
        self.init()
        df = self.runSimulation()
        self._checkWithSimulation(df)
        #
        self.init(do_simulate_on_update=True)
        df = self.runSimulation()
        self._checkWithSimulation(df)

    def testOutfcn(self):
        if IGNORE_TEST:
          return
        self.init()
        time = 3
        self.sys.setTime(time)
        x_vec = self.sys.makeStateSer().values
        u_vec = np.repeat(0, len(self.sys.input_names))
        out_vec = self.sys.outfcn(time, x_vec, u_vec, None)
        self.assertEqual(len(out_vec), self.sys.num_output)

    def testOutlist(self):
        if IGNORE_TEST:
          return
        self.init()
        names = self.sys.outlist
        self.assertTrue(all([self.sys.name in n for n in names]))
        self.assertEqual(len(names), self.sys.num_output)

    def testGetStateSer(self):
        if IGNORE_TEST:
          return
        self.init()
        ser_0 = self.sys.makeStateSer()
        ser_1 = self.sys.makeStateSer(time=1)
        self.assertGreater(ser_0.loc["S0"], ser_1.loc["S1"])

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
        x_vec = self.sys.makeStateSer().values
        u_vec = np.repeat(u_val, self.sys.num_input)
        u_vecs = np.array([np.array(u_vec) for _ in range(NUM_TIME)])
        u_vecs = np.reshape(u_vecs, (1, NUM_TIME))
        t, y = control.input_output_response(self.sys, TIMES, u_vecs, x_vec)
        df = pd.DataFrame(y, self.sys.output_names)
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
        self.assertGreater(df2["S1"].values[1], df1["S1"].values[1])

    def testWithInputOutputResponseWithoutEffector(self):
        if IGNORE_TEST:
            return
        self.init()
        self.sys = ctl.NonlinearIOSystem("test_sys", self.ctlsb)
        df = self.runInputOutputResponse(0)


if __name__ == '__main__':
  unittest.main()

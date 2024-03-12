import controlSBML as ctl
import controlSBML.constants as cn
import controlSBML.util as util
import helpers
from controlSBML.option_management.option_manager import OptionManager

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
#helpers.setupPlotting(__file__)

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
J0: ->  S0; 1 + E_J0
J1: S0 -> S1; k1*S0
J2: S1 -> S2; k2*S1*S1
J3: S2 -> $S3; k3*S2*S1

S0 = 10
S1 = 0
S2 = 0
E_J0 = 0
$S3 = 3
k1 =1; k2=2; k3=3
"""
LINEAR_MDL2 = """
J0:  -> S1; k0
J0a: -> S1; S0
J0b: S0 ->; S0
J1: S1 -> S2; k1*S1
J2: S2 -> S3; k2*S2
J3: S3 -> S4; k3*S3
J4: S4 -> ; k4*S4

S0 = 0
k0 = 50
k1 = 1
k2 = 2
k3 = 3
k4 = 4
S1 = 1
S2 = 1
S3 = 1
S4 = 1
"""
INPUT_NAME2 = "S2"
OUTPUT_NAME2 = "S4"

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
        for time in TIMES:
            dx_vec = self.sys.updfcn(time, x_vec, u_vec, {})
            new_x_vec = x_vec + DT*dx_vec
            [dct[n].append(v) for n, v in
                  zip(self.sys.state_names, new_x_vec)]
            x_vec = new_x_vec
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
        common_columns = list(df.columns.intersection(df_rr.columns))
        df_rr = df_rr[common_columns]
        ssq = 0.0
        for column in common_columns:
            ssq += np.sum((df[column].values - df_rr[column].values)**2)
        # Calculate error per value
        rms_per_value = np.sqrt(ssq/(len(df.columns)*len(df)))
        self.assertLess(rms_per_value, 1e-1)
        if IS_PLOT:
            mgr = OptionManager({})
            num_col = len(common_columns)
            _, axes = plt.subplots(num_col, 1)
            for i, column in enumerate(common_columns):
                axes[i].scatter(df_rr[column], df[column])
                axes[i].set_xlabel("simulated")
                axes[i].set_ylabel("predicted")
                axes[i].set_title(column)
                max_val = max(df_rr[column].max(), df[column].max())
                axes[i].plot([0, max_val], [0, max_val], color="red", linestyle="--")
            mgr.doFigOpts()

    def runInputOutputResponse(self, u_val, end_time=END_TIME, sys=None):
        if sys is None:
            sys = self.sys
        u_val = np.array(u_val)
        #
        x_vec = sys.makeStateSer().values
        u_vec = np.repeat(u_val, sys.num_input)
        u_vecs = np.array([np.array(u_vec) for _ in range(NUM_TIME)])
        u_vecs = np.reshape(u_vecs, (1, NUM_TIME))
        t, y = control.input_output_response(sys, TIMES, U=u_vecs, X0=x_vec)
        df = pd.DataFrame(y, sys.output_names)
        df = df.transpose()
        df.index = t
        return df

    def testWithInputOutputResponse(self):
        if IGNORE_TEST:
            return
        self.init()
        df = self.runInputOutputResponse(0)
        J1_flux_0 = self.ctlsb.roadrunner["J1"]
        self._checkWithSimulation(df)
        # Output increases with inpu
        _ = self.runInputOutputResponse(10)
        J1_flux_10 = self.ctlsb.roadrunner["J1"]
        self.assertGreater(J1_flux_10, J1_flux_0)

    def testLogger(self):
        if IGNORE_TEST:
            return
        def test(df):
            self.assertTrue(isinstance(df, pd.DataFrame))
            self.assertGreater(len(df), 0)
            true_dct = {c: c in df.columns for c in sys.state_names}
            self.assertTrue(all(true_dct.values()))
        #
        if IS_PLOT:
            cn.PLOT_DIR = cn.TEST_DIR
        else:
            cn.PLOT_DIR= tempfile.mkdtemp()
        ctlsb = ctl.ControlSBML(NONLINEAR_MDL,
              input_names=["E_J0"], output_names=["S1", "S2"])
        sys = ctl.NonlinearIOSystem("test_sys", ctlsb, is_log=True)
        _ = self.runInputOutputResponse(0, sys=sys)
        df1 = sys.updfcn_logger.report()
        test(df1)
        df2 = sys.outfcn_logger.report()
        test(df2)
       
    def testWithInputOutputResponseWithoutEffector(self):
        if IGNORE_TEST:
            return
        self.init()
        self.sys = ctl.NonlinearIOSystem("test_sys", self.ctlsb)
        df = self.runInputOutputResponse(0)

    def testAccuracy(self):
        if IGNORE_TEST:
            return
        ctlsb = ctl.ControlSBML(LINEAR_MDL2,
              input_names=[INPUT_NAME2], output_names=[OUTPUT_NAME2])
        sys = ctl.NonlinearIOSystem("linear_mdl2", ctlsb)
        x_vec = sys.makeStateSer().values
        times = np.linspace(0, 1000, 10000)
        times, y_vals = util.control.input_output_response(sys, T=times, X0=x_vec, U=50)
        self.assertLess(np.abs(y_vals[-1] - 25), 0.1)   # S4 = k2/k4*U at steady state

    def testFixedRateInput(self):
        # Simulating with roadrunner should produce the same results as simulating with controlSBML
        if IGNORE_TEST:
            return
        mgr = OptionManager({})
        times = np.linspace(0, 5, 50)
        # Simulate in road runner
        rr = te.loada(LINEAR_MDL2)
        data = rr.simulate(times[0], times[-1], len(times))
        rr_vals = np.array(data["[S3]"])
        # Simulate in controlSBML
        ctlsb = ctl.ControlSBML(LINEAR_MDL2,
              input_names=["S1"], output_names=["S3"])
        U = ctlsb.get("k0")
        ctlsb.set({"k0": 0})
        dct = {"S1": False}
        sys = ctl.NonlinearIOSystem("linear_mdl2", ctlsb,
                                    is_fixed_input_species=dct)
        x_vec = sys.makeStateSer().values
        times, y_vals = util.control.input_output_response(sys, T=times, X0=x_vec, U=U)
        rms = (np.sum((rr_vals - y_vals)**2)/len(rr_vals))**0.5
        self.assertLess(rms, 1e-2)
        if IS_PLOT:
            plt.scatter(rr_vals, y_vals)
            plt.xlabel("Roadrunner")
            plt.ylabel("controlSBML")
            plt.plot([0, max(y_vals)], [0, max(y_vals)], color="red", linestyle="--")
            mgr.doPlotOpts()
            mgr.doFigOpts()

    def testFixedConcentrationInput(self):
        # Simulating with roadrunner should produce the same results as simulating with controlSBML
        if IGNORE_TEST:
            return
        mgr = OptionManager({})
        times = np.linspace(0, 5, 50)
        # Simulate in road runner
        rr = te.loada(LINEAR_MDL2)
        data = rr.simulate(times[0], times[-1], len(times))
        rr_vals = np.array(data["[S3]"])
        # Simulate in controlSBML
        ctlsb = ctl.ControlSBML(LINEAR_MDL2,
              input_names=["S0"], output_names=["S3"])
        U = ctlsb.get("k0")
        ctlsb.set({"k0": 0})
        dct = {"S0": True}
        sys = ctl.NonlinearIOSystem("linear_mdl2", ctlsb,
                                    is_fixed_input_species=dct)
        x_vec = sys.makeStateSer().values
        times, y_vals = util.control.input_output_response(sys, T=times, X0=x_vec, U=U)
        if IS_PLOT:
            plt.scatter(rr_vals, y_vals)
            plt.xlabel("Roadrunner")
            plt.ylabel("controlSBML")
            plt.plot([0, max(y_vals)], [0, max(y_vals)], color="red", linestyle="--")
            mgr.doPlotOpts()
            mgr.doFigOpts()
        rms = (np.sum((rr_vals - y_vals)**2)/len(rr_vals))**0.5
        self.assertLess(rms, 1e-2)

    def testSetSteadyState(self):
        if IGNORE_TEST:
            return
        ctlsb = ctl.ControlSBML(LINEAR_MDL2,
              input_names=["S1"], output_names=["S3"])
        sys = ctl.NonlinearIOSystem("test_sys", ctlsb)
        initial_value = ctlsb.get("S3")
        sys.setSteadyState()
        final_value = ctlsb.get("S3")
        self.assertGreater(final_value, initial_value)


if __name__ == '__main__':
  unittest.main()
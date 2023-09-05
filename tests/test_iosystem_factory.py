from controlSBML.iosystem_factory import IOSystemFactory
from controlSBML import *
import controlSBML as ctl
import controlSBML.constants as cn
import controlSBML.util as util

import control
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import tellurium as te
import unittest


IGNORE_TEST = True
IS_PLOT = False
SIZE = 20
if IS_PLOT:
    import matplotlib
    matplotlib.use('TkAgg')
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PLOT_PATH = os.path.join(TEST_DIR, "plot.pdf")
TIMES = ctl.makeSimulationTimes(end_time=50)
MODEL_BUG = """
S1 -> S2; k1*S1
J1: S2 -> S3; k2*S2
J2: S3 -> S2; k3*S3
J3: S2 -> ; k4*S2

k1 = 1
k2 = 2
k3 = 3
k4 = 4
S1 = 10
S2 = 0
S3 = 0
S4 = 0
"""
MODEL = """
J1: -> S0; 5
S0 -> S1; k0*S0
S1 -> S2; k1*S1
S2 -> S1; k2*S2
S2 -> ; k3*S2

k0 = 0.5
k1 = 1
k2 = 2
k3 = 3
S0 = 5
S1 = 0
S2 = 0
"""

def makeSystemForModel():
    A = [ 
    [-0.5, 0, 0],
    [0.5, -1, 2],
    [0, 1, -5],
        ]
    B = [[5], [0], [0]]
    C = np.eye(3)
    model_sys = control.ss(A, B, C, 0)
    X0=[5, 0, 0]
    if IS_PLOT:
        times, yvals = control.step_response(model_sys, np.linspace(0, 10, 100), X0=X0)
        legends = []
        _, axes = plt.subplots(2)
        for idx in range(yvals.shape[0]):
            ys = np.reshape(yvals[idx], (100,1))
            legends.append("S%d" % idx)
            axes[0].plot(times, ys)
        axes[0].legend(legends)
        rr = te.loada(MODEL)
        arr = rr.simulate(0, 10, 100)
        df = pd.DataFrame(arr, columns=["time", "S0", "S1", "S2"])
        df.plot(x="time", y=["S0", "S1", "S2"], ax=axes[1])
        plt.savefig("plot.pdf")
    return model_sys, X0



#############################
# Tests
#############################
class TestIOSystemFactory(unittest.TestCase):

    def setUp(self):
        self.factory = IOSystemFactory()

    def runController(self, name="controller", is_log=False, U=None, **kwargs):
        times = list(range(SIZE))
        if U is None:
           U = times
        factory = IOSystemFactory(is_log=is_log)
        controller = factory.makePIDController(name, **kwargs)
        return factory, control.input_output_response(controller, T=times, U=U)

    # TODO: More tests for integral and differential control
    def testMakePIDController(self):
        #if IGNORE_TEST:
        #  return
        kp = 2
        factory, result = self.runController(kp=kp)
        trues = [r == kp*( t) for t, r in zip(result.t, result.outputs)]
        self.assertTrue(all(trues))
        #
        _, result_ki = self.runController(ki=kp)
        _, result_kd = self.runController(kd=kp)
        self.assertGreater(result_ki.y[0][-1], result_kd.y[0][-1])
        U = np.array(range(SIZE//2))
        U = np.concatenate([U, -U])
        _, result_kd1 = self.runController(kd=kp, is_nonnegative_output=True, U=U)
        tot = np.sum(result_kd1.y[SIZE//2:])
        self.assertEqual(tot, 0)

    def testMakeSinusoid(self):
        if IGNORE_TEST:
          return
        sys = self.factory.makeSinusoid("sine", amplitude=10, frequency=20)
        result = control.input_output_response(sys, T=TIMES)
        self.assertTrue(len(result.y) > 0)
        self.assertTrue(np.var(result.y) > 0)

    def testMakeAdder1(self):
        if IGNORE_TEST:
          return
        adder = self.factory.makeAdder("adder", num_input=3)
        times = [0, 1, 2, 3, 4]
        u_arr = np.array([[t, t, t] for t in times])
        u_arr = u_arr.transpose()
        result = control.input_output_response(adder, T=times, U=u_arr)
        trues = [r == 3*t for t, r in zip(times, result.outputs.flatten())]
        self.assertTrue(all(trues))

    def testMakeAdder2(self):
        # Check on use of minus sign
        if IGNORE_TEST:
          return
        input_names = ["a", "b", "c", "-d"]
        num_input = len(input_names)
        adder = self.factory.makeAdder("adder", input_names=input_names)
        self.assertEqual(adder.ninputs, num_input)
        input_values = [0, 1, 2, 3, 4]
        u_arr = np.array([np.repeat(t, len(input_names)) for t in input_values])
        u_arr = u_arr.transpose()
        result = control.input_output_response(adder, T=input_values, U=u_arr)
        trues = [r == (num_input-2)*t for t, r in zip(input_values, result.outputs.flatten())]
        self.assertTrue(all(trues))

    def testMakeFilter(self):
        if IGNORE_TEST:
          return
        sys = self.factory.makeFilter("filter", 1)
        length = len(TIMES)
        mean = 5
        U = np.random.normal(mean, 1, length)
        result = control.input_output_response(sys, T=TIMES, U=U)
        y_values = result.y.flatten()
        times = result.t.flatten()
        self.assertTrue(len(result.y) > 0)
        lin_sys = sys.linearize(x0=0, u0=0)
        self.assertTrue(lin_sys.dcgain() == 1)
        if IS_PLOT:
             self.makePlot(sys, times=TIMES, U=U)

    def makePlot(self, sys, times=TIMES, U=None):
        if U is None:
            result = control.input_output_response(sys, T=times)
        else:
            result = control.input_output_response(sys, T=times, U=U)
        plt.plot(result.t, result.outputs.flatten())
        plt.xlabel("time")
        plt.savefig(PLOT_PATH)

    def testMakeFilterZeroInput(self):
        if IGNORE_TEST:
          return
        sys = self.factory.makeFilter("filter", 0)
        length = len(TIMES)
        results = control.input_output_response(sys, T=TIMES, U=1)
        outputs = results.outputs.flatten()
        self.assertEqual(np.mean(outputs), 1)
        self.assertEqual(np.var(outputs), 0)

    def testMakeArbitrarySignal(self):
        if IGNORE_TEST:
          return
        sys = self.factory.makeArbitrarySignal("arbitrary", signal_function=lambda t: t)
        result = control.input_output_response(sys, T=TIMES)
        self.assertTrue(len(result.y) > 0)
        self.assertTrue(np.var(result.y) > 0)
        if IS_PLOT:
            self.makePlot(sys)

    def testMakeStep(self):
        if IGNORE_TEST:
          return
        step_size = 3
        sys = self.factory.makeStep("step", step_size=step_size)
        result = control.input_output_response(sys, T=TIMES)
        self.assertTrue(np.var(result.y) == 0)
        self.assertTrue(result.y.flatten()[0] == step_size)

    def testMakeMultiplier(self):
        if IGNORE_TEST:
          return
        factor = 3
        sys = self.factory.makeMultiplier("multiplier", factor)
        U = np.repeat(1, len(TIMES)).flatten()
        result = control.input_output_response(sys, T=TIMES.flatten(), U=U)
        self.assertTrue(np.var(result.y) == 0)
        self.assertTrue(result.y.flatten()[0] == factor)

    def testPassthru(self):
        if IGNORE_TEST:
          return
        sys = self.factory.makePassthru("passthru")
        result = control.input_output_response(sys, T=TIMES.flatten(), U=TIMES)
        trues = [x == y for x, y in zip(TIMES, result.outputs.flatten())]
        self.assertTrue(all(trues))

    def testMakeLogReport(self):
        if IGNORE_TEST:
          return
        def test(is_log=False):
            factory, _ = self.runController(is_log=is_log)
            df = factory.report()
            if not is_log:
               self.assertIsNone(df)
            else:
                self.assertGreater(len(df), 0)
                self.assertEqual(len(df.columns), 4)
            return df
        #
        _ = test(is_log=False)
        df = test(is_log=True)

    # TODO: Needs more tests
    def testMakeFullStateController(self):
        if IGNORE_TEST:
            return
        model_sys, X0 = makeSystemForModel()
        state_names = ["S0", "S1", "S2"]
        controller = self.factory.makeFullStateController("controller",
              model_sys, state_names, dcgain=1.0, poles=-10)
        times = [0.1*n for n in range(50)]
        result = control.input_output_response(controller, T=times, U=1)
        outputs = result.outputs[0]
        self.assertEqual(len(times), len(outputs))

    def testBug1(self):
        if IGNORE_TEST:
          return
        # Elements of the system
        kp = 1
        ki = 1
        kd = 0
        dt = 1/cn.POINTS_PER_TIME
        factory = ctl.IOSystemFactory(dt=dt)
        # Create the elements of the feedback loop
        noise = factory.makeSinusoid("noise", amplitude=0, frequency=20)
        disturbance = factory.makeSinusoid("disturbance", 0, 2)
        
        ctlsb = ctl.ControlSBML(MODEL_BUG, input_names=["S1"], output_names=["S3"])
        system = ctlsb.makeNonlinearIOSystem("system")
        controller = factory.makePIDController("controller", kp=kp, ki=ki, kd=kd)
        fltr = factory.makePassthru("fltr")
        
        sum_N_Y = factory.makeAdder("sum_N_Y")
        sum_D_U = factory.makeAdder("sum_D_U")
        sum_F_R = factory.makeAdder("sum_F_R")
        
        # Create the closed loop system
        closed_loop = control.interconnect(
          [noise, disturbance, sum_N_Y, sum_F_R, sum_D_U, system, fltr, controller ], 
          connections=[
            ['controller.in', 'sum_F_R.out'],    # e(t)
            ['sum_D_U.in1', 'controller.out'],   # u(t)
            ['sum_D_U.in2', 'disturbance.out'],  # d(t)
            ['system.S1',   'sum_D_U.out'],
            ['sum_N_Y.in1', 'system.S3'],        # y(t)
            ['sum_N_Y.in2', 'noise.out'],        # n(t)
            ['fltr.in',     'sum_N_Y.out'],
            ['sum_F_R.in1', '-fltr.out'],
          ],
          inplist=["sum_F_R.in2"],
          #outlist=["sum_F_R.in2", "sum_N_Y.out", 'system.S2', 'system.S3'],
          outlist=['system.S3'],
        )
        
        times = ctl.makeSimulationTimes(0, 50, 1/dt)
        X0 = ctl.makeStateVector(closed_loop)
        U = 1
        result = control.input_output_response(closed_loop, T=times, U=U, X0=X0)
        plt.plot([result.t[0], result.t[-1]], [U, U])
        plt.plot(result.t, result.outputs.flatten())
        #plt.plot(result.t, result.outputs[2].flatten())
        #plt.plot(result.t, result.outputs[3].flatten())
        plt.ylim([0, 2])
        legends = ["reference", "output"]
        plt.legend(legends)
        if IS_PLOT:
            plt.savefig(PLOT_PATH)

    def testMakeStateFilter(self):
        if IGNORE_TEST:
          return
        filter_name = "state_filter"
        state_names = ["a", "b", "c"]
        KF = -1
        sys, connections = self.factory.makeStateFilter(
              filter_name, KF, "in_sys", "ot_sys", state_names)
        results = control.input_output_response(sys, T=TIMES, U=1)
        self.assertEqual(np.shape(results.outputs)[0], len(state_names))
        #
        kfs = np.repeat(KF, len(state_names))
        sys, connections = self.factory.makeStateFilter(
              filter_name+"1", kfs, "in_sys", "ot_sys", state_names)
        result2s = control.input_output_response(sys, T=TIMES, U=1)
        total = np.sum((results.outputs - result2s.outputs)**2)
        self.assertTrue(np.isclose(total, 0))
        #
        kfs = [-1 - n for n in range(len(state_names))]
        sys, connections = self.factory.makeStateFilter(
              filter_name+"2", kfs, "in_sys", "ot_sys", state_names)
        result3s = control.input_output_response(sys, T=TIMES, U=1)
        self.assertEqual(np.shape(result3s.outputs)[0], len(state_names))
        if IS_PLOT:
            ts = ctl.timeresponse2Timeseries(result3s)
            ctl.plotOneTS(ts)
            plt.savefig(PLOT_PATH)

    def testMakeStateFilterZeroFilter(self):
        if IGNORE_TEST:
          return
        filter_name = "state_filter"
        state_names = ["a", "b", "c"]
        kf = 0
        sys, connections = self.factory.makeStateFilter(
              filter_name, kf, "in_sys", "ot_sys", state_names)
        results = control.input_output_response(sys, T=TIMES, U=1)
        self.assertEqual(np.shape(results.outputs)[0], len(state_names))
        outputs = results.outputs.flatten()
        self.assertEqual(np.mean(outputs), 1)
        self.assertEqual(np.var(outputs), 0)
    
    def testMakeArbitrarySignalTimeshift(self):
        if IGNORE_TEST:
            return
        start_time = 1
        end_time = 5
        def shiftedRamp(time):
           if time > start_time:
              return time - start_time
           else:
              return 0
        #
        factory = IOSystemFactory()
        sys = factory.makeArbitrarySignal("shifted_ramp", shiftedRamp, start_time=start_time, end_time=end_time)
        result = control.input_output_response(sys, T=TIMES)
        y_values = result.y.flatten()
        plt.plot(result.t, result.outputs.flatten())
        max_value = end_time - start_time
        self.assertEqual(max(y_values), max_value)
        self.assertEqual(min(y_values), 0)
        if IS_PLOT:
            plt.savefig(PLOT_PATH)


if __name__ == '__main__':
  unittest.main()

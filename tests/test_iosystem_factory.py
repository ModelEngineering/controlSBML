from controlSBML.iosystem_factory import IOSystemFactory
import controlSBML as ctl
import controlSBML.constants as cn

import control
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False
if IS_PLOT:
    import matplotlib
    matplotlib.use('TkAgg')
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
-> S0; 5
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


#############################
# Tests
#############################
class TestIOSystemFactory(unittest.TestCase):

    def setUp(self):
        self.factory = IOSystemFactory()

    def runController(self, name="controller", **kwargs):
        factory = IOSystemFactory()
        controller = factory.makePIDController(name, **kwargs)
        times = list(range(20))
        return factory, control.input_output_response(controller, T=times, U=times)

    # TODO: More tests for integral and differential control
    def testMakePIDController(self):
        if IGNORE_TEST:
          return
        kp = 2
        factory, result = self.runController(kp=kp)
        trues = [r == kp*( t) for t, r in zip(result.t, result.outputs)]
        self.assertTrue(all(trues))
        #
        _, result_ki = self.runController(ki=kp)
        _, result_kd = self.runController(kd=kp)
        self.assertGreater(result_ki.y[0][-1], result_kd.y[0][-1])

    def testMakeSinusoid(self):
        if IGNORE_TEST:
          return
        sys = self.factory.makeSinusoid("sine", 10, 20)
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
        if IGNORE_TEST:
          return
        input_names = ["a", "b", "c", "d"]
        num_input = len(input_names)
        adder = self.factory.makeAdder("adder", input_names=input_names)
        self.assertEqual(adder.ninputs, num_input)
        times = [0, 1, 2, 3, 4]
        u_arr = np.array([np.repeat(t, len(input_names)) for t in times])
        u_arr = u_arr.transpose()
        result = control.input_output_response(adder, T=times, U=u_arr)
        trues = [r == num_input*t for t, r in zip(times, result.outputs.flatten())]
        self.assertTrue(all(trues))

    def testMakeFilter(self):
        if IGNORE_TEST:
          return
        sys = self.factory.makeFilter("filter", -1)
        length = len(TIMES)
        mean = 5
        U = np.random.normal(mean, 1, length)
        result = control.input_output_response(sys, T=TIMES, U=U)
        y_values = result.y.flatten()
        times = result.t.flatten()
        self.assertTrue(len(result.y) > 0)
        self.assertLess(np.abs(np.mean(y_values[-10:]) - mean), mean/10)
        lin_sys = sys.linearize(x0=0, u0=0)
        self.assertTrue(lin_sys.dcgain() == 1)
        if IS_PLOT:
             plt.scatter(result.t.flatten(), result.y.flatten())
             plt.xlabel("time")
             plt.ylabel("filter output")
             plt.show()

    def testMakeConstant(self):
        if IGNORE_TEST:
          return
        constant = 3
        sys = self.factory.makeConstant("constant", constant)
        result = control.input_output_response(sys, T=TIMES)
        self.assertTrue(np.var(result.y) == 0)
        self.assertTrue(result.y.flatten()[0] == constant)

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
        def test():
            factory, _ = self.runController()
            df = factory.report()
            self.assertGreater(len(df), 0)
            self.assertEqual(len(df.columns), 4)
        #
        test()

    def testMakeFullStateController(self):
        if IGNORE_TEST:
          return
        ctlsb = ctl.ControlSBML(MODEL, input_names=["S0"], output_names=["S2"])
        controller = self.factory.makeFullStateController("controller",
              ctlsb, factor=1.0, poles=-10, time=1)
        times = [0, 1, 2, 3, 4]
        U = np.array([(1, 1, 1,) for _ in range(len(times))])
        U = U.transpose()
        result = control.input_output_response(controller, T=times, U=U)
        outputs = result.outputs[0]
        self.assertEqual(len(times), len(outputs))
        self.assertEqual(len(set(outputs)), 1)

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
        noise = factory.makeSinusoid("noise", 0, 20)
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
            plt.show()


if __name__ == '__main__':
  unittest.main()

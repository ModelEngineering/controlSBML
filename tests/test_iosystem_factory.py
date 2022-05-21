from controlSBML.iosystem_factory import IOSystemFactory
import controlSBML as ctl

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
TIMES = ctl.makeSimulationTimes(0, 5, 50)
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
        self.factory = IOSystemFactory(
              callback_fcn=ctl.IOSystemFactory_CALLBACK_REPORT)

    def runController(self, name="controller",
          callback_fcn=ctl.IOSystemFactory_CALLBACK_REPORT, **kwargs):
        factory = IOSystemFactory(callback_fcn=callback_fcn)
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
        self.assertGreater(len(self.factory.callback_log), 0)

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
        self.assertGreater(len(self.factory.callback_log), 0)

    def testMakeAdder2(self):
        if IGNORE_TEST:
          return
        input_names = ["a", "b", "c", "d"]
        num_input = len(input_names)
        adder = self.factory.makeAdder("adder", input_names=input_names)
        self.assertEquals(adder.ninputs, num_input)
        times = [0, 1, 2, 3, 4]
        u_arr = np.array([np.repeat(t, len(input_names)) for t in times])
        u_arr = u_arr.transpose()
        result = control.input_output_response(adder, T=times, U=u_arr)
        trues = [r == num_input*t for t, r in zip(times, result.outputs.flatten())]
        self.assertTrue(all(trues))
        self.assertGreater(len(self.factory.callback_log), 0)

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
        self.assertGreater(len(self.factory.callback_log), 0)
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
        self.assertGreater(len(self.factory.callback_log), 0)

    def testMakeMultiplier(self):
        if IGNORE_TEST:
          return
        factor = 3
        sys = self.factory.makeMultiplier("multiplier", factor)
        U = np.repeat(1, len(TIMES)).flatten()
        result = control.input_output_response(sys, T=TIMES.flatten(), U=U)
        self.assertTrue(np.var(result.y) == 0)
        self.assertTrue(result.y.flatten()[0] == factor)
        self.assertGreater(len(self.factory.callback_log), 0)

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
        self.assertEquals(len(times), len(outputs))
        self.assertEquals(len(set(outputs)), 1)


if __name__ == '__main__':
  unittest.main()

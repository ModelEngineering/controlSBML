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
        _, result = self.runController(kp=kp)
        trues = [r == kp*( t) for t, r in zip(result.t, result.outputs)]
        self.assertTrue(all(trues))
        #
        _, result_ki = self.runController(ki=kp)
        _, result_kd = self.runController(kd=kp)
        trues = [v1 <= v2 for v1, v2 in 
              zip(result_kd.y.flatten(), result_ki.y.flatten())]
        self.assertTrue(all(trues))

    def testMakeSinusoid(self):
        if IGNORE_TEST:
          return
        sys = self.factory.makeSinusoid("sine", 10, 20)
        result = control.input_output_response(sys, T=TIMES)
        self.assertTrue(len(result.y) > 0)
        self.assertTrue(np.var(result.y) > 0)

    def testMakeAdder(self):
        if IGNORE_TEST:
          return
        adder = self.factory.makeAdder("adder", num_input=3)
        times = [0, 1, 2, 3, 4]
        u_arr = np.array([[t, t, t] for t in times])
        u_arr = u_arr.transpose()
        result = control.input_output_response(adder, T=times, U=u_arr)
        trues = [r == 3*t for t, r in zip(times, result.outputs.flatten())]
        self.assertTrue(all(trues))

    def testMakeFilter(self):
        if IGNORE_TEST:
          return
        sys = self.factory.makeFilter("filter", -1)
        U = range(len(TIMES)) + 0.1*np.random.rand(len(TIMES))
        result = control.input_output_response(sys, T=TIMES, U=U)
        self.assertTrue(len(result.y) > 0)
        lin_sys = sys.linearize(x0=0, u0=0)
        self.assertTrue(lin_sys.dcgain() == 1)
        df = self.factory.report()
        if IS_PLOT:
             plt.scatter(df["in"], df["out"])
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


if __name__ == '__main__':
  unittest.main()

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
        self.factory = IOSystemFactory(
              callback_fcn=ctl.IOSystemFactor_CALLBACK_REPORT)

    def runController(self, name="controller",
          callback_fcn=ctl.IOSystemFactor_CALLBACK_REPORT, **kwargs):
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
        trues = [v1 <= v2 for v1, v2 in 
              zip(result_kd.y.flatten(), result_ki.y.flatten())]
        self.assertTrue(all(trues))
        self.assertGreater(len(factory.callback_log), 0)

    def testMakeSinusoid(self):
        if IGNORE_TEST:
          return
        sys = self.factory.makeSinusoid("sine", 10, 20)
        result = control.input_output_response(sys, T=TIMES)
        self.assertTrue(len(result.y) > 0)
        self.assertTrue(np.var(result.y) > 0)
        self.assertGreater(len(self.factory.callback_log), 0)

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

    def testCallout(self):
        if IGNORE_TEST:
          return
        global updfcn_cnt, outfcn_cnt
        updfcn_cnt = 0
        outfcn_cnt = 0
        def callback_fcn(name, time, x_vec, u_vec, dct, out_vec):
            global updfcn_cnt, outfcn_cnt
            if "upd" in name:
                updfcn_cnt += 1
            if "out" in name:
                outfcn_cnt += 1
            return "controller"
        #
        factory, _ = self.runController(callback_fcn=callback_fcn)
        self.assertGreater(updfcn_cnt, outfcn_cnt)


if __name__ == '__main__':
  unittest.main()

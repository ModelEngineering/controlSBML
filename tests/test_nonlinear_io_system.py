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

############### TEST HELPERS #################
def predprey_rhs(t, x, u, params):
    # Parameter setup
    a = params.get('a', 3.2)
    b = params.get('b', 0.6)
    c = params.get('c', 50.)
    d = params.get('d', 0.56)
    k = params.get('k', 125)
    r = params.get('r', 1.6)

    # Map the states into local variable names
    H = x[0]
    L = x[1]

    # Compute the control action (only allow addition of food)
    u_0 = u if u > 0 else 0

    # Compute the discrete updates
    dH = (r + u_0) * H * (1 - H/k) - (a * H * L)/(c + H)
    dL = b * (a * H *  L)/(c + H) - d * L

    return [dH, dL]


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
              input_names=["J0"], output_names=["S2"])
        self.sys = ctl.NonlinearIOSystem(None, None, ctlsb=self.ctlsb,
              inputs=self.ctlsb.input_names, outputs=self.ctlsb.output_names,
              states=self.ctlsb.state_names, name="test_sys",
              effector_dct=EFFECTOR_DCT)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.sys.state_call_dct, dict))

    def testWrapNonTellurium(self):
        if IGNORE_TEST:
          return
        io_predprey = ctl.NonlinearIOSystem(
            predprey_rhs, None, inputs=('u'), outputs=('H', 'L'),
            states=('H', 'L'), name='predprey')
        # Simulate the system
        X0 = [25, 20]                 # Initial H, L
        T = np.linspace(0, 70, 500)   # Simulation 70 years of time
        
        # Simulate the system
        t, y = control.input_output_response(io_predprey, T, 0, X0)
        self.assertGreater(len(t), 0)
        self.assertGreater(len(y), 0)
        self.assertTrue(isinstance(io_predprey.state_call_df, pd.DataFrame))
        
        # Plot the response
        if IS_PLOT:
            plt.figure(1)
            plt.plot(t, y[0])
            plt.plot(t, y[1])
            plt.legend(['Hare', 'Lynx'])
            plt.show()

    def testCtlsbUpdfcn(self):
        if IGNORE_TEST:
          return
        self.init()
        num_time = 51
        x_vec = self.ctlsb.state_ser.values
        u_vec = np.repeat(1, len(self.ctlsb.input_names))
        times = [n*0.1 for n in range(0, num_time)]
        dct = {n: [] for n in self.ctlsb.state_names}
        for time in times:
            _ = self.sys.ctlsbUpdfcn(time, x_vec, u_vec, {})
            self.sys._is_first_state_call = False
            states = [self.ctlsb.get(n) for n in self.ctlsb.state_names]
            [dct[self.ctlsb.state_names[i]].append(states[i])
                  for i in range(self.ctlsb.num_state)]
            x_vec = self.ctlsb.state_ser.values
        df = pd.DataFrame(dct)
        last_value = df.values[-1,:]
        self.assertGreater(last_value[-1], last_value[0])


if __name__ == '__main__':
  unittest.main()

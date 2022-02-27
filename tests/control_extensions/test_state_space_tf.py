from controlSBML.control_extensions.state_space_tf import StateSpaceTF

import control
import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = True
IS_PLOT = True
if IS_PLOT:
    import matplotlib
    matplotlib.use('TkAgg')
# x''' + 3*x'' + 2x' + x = 0
SIZE = 2
A1_MAT = np.array([
    [0, 1, 0],
    [0, 0, 1],
    [-3, -2, -1]
    ])
A_MAT = [[1., -2], [3, -4]]
B_MAT = np.identity(SIZE)
C_MAT = np.identity(SIZE)
D_MAT = 0*np.identity(SIZE)


#############################
# Tests
#############################
class TestStateSpaceTF(unittest.TestCase):

    def setUp(self):
        self.sys = control.StateSpace(A_MAT, B_MAT, C_MAT, D_MAT)
        self.ss_tf = StateSpaceTF(self.sys)

    def testStr(self):
        if IGNORE_TEST:
          return
        if IGNORE_TEST:
            print(self.ss_tf)
        self.assertTrue("s" in str(self.ss_tf))
        self.assertTrue("-----" in str(self.ss_tf))

    def testConstructor(self):
        if IGNORE_TEST:
          return
        self.assertTrue(isinstance(self.ss_tf.dataframe, pd.DataFrame))

    def testGetSystemShape(self):
        if IGNORE_TEST:
          return
        num_state, num_input, num_output = self.ss_tf.getSystemShape(self.sys)
        self.assertEqual(num_state, self.sys.nstates)
        self.assertEqual(num_input, len(self.ss_tf.dataframe))
        self.assertEqual(num_output, len(self.ss_tf.dataframe.columns))

    def testPlotBode(self):
        # TESTING
        self.ss_tf.plotBode(is_plot=IS_PLOT, legend_crd=(0.9, 1))


if __name__ == '__main__':
  unittest.main()

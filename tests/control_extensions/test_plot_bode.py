from controlSBML.control_extensions.plot_bode import plotBode

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
A_MAT = np.array([
    [0, 1, 0],
    [0, 0, 1],
    [-3, -2, -1]
    ])


#############################
# Tests
#############################
class TestFunction(unittest.TestCase):

    def setUp(self):
        self.sys = control.StateSpace(A_MAT, 0*A_MAT, 0*A_MAT, 0*A_MAT)

    def testFunction(self):
        # TESTING
        plotBode(self.sys)


if __name__ == '__main__':
  unittest.main()

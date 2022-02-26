from controlSBML.option_management.option_manager import OptionManager
import controlSBML.constants as cn

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import unittest


IGNORE_TEST = False
IS_PLOT = False
DCT = dict(cn.PLOT_DCT)
DCT.update(cn.SIM_DCT)
if IS_PLOT:
    import matplotlib
    matplotlib.use('TkAgg')



#############################
# Tests
#############################
class TestOptionManager(unittest.TestCase):

    def setUp(self):
        # Cannot modify self.control
        dct = dict(DCT)
        _, self.ax = plt.subplots(1)
        dct[cn.O_AX] = self.ax
        dct[cn.O_XLABEL] = "xlabel"
        dct[cn.O_IS_PLOT] = IS_PLOT
        self.option_mgr = OptionManager(dct, cn.DEFAULT_DCTS)

    def tearDown(self):
        plt.close()

    def testConstructor(self):
        if IGNORE_TEST:
          return
        diff = set(cn.PLOT_DCT.keys()).symmetric_difference(
              self.option_mgr.plot_opts)
        self.assertEqual(len(diff), 0)

    def testDoPlotOptsFigOpts(self):
        if IGNORE_TEST:
          return
        self.ax.plot(range(5), range(5))
        self.option_mgr.doPlotOpts()
        self.option_mgr.doFigOpts()



if __name__ == '__main__':
  unittest.main()

from controlSBML.option_management.option_manager import OptionManager
import controlSBML.constants as cn

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import unittest


IGNORE_TEST = True
IS_PLOT = True
DCT = dict(cn.PLOT_DCT)
DCT.update(cn.SIM_DCT)
if IS_PLOT:
    import matplotlib
    matplotlib.use('TkAgg')

# Temporarily change the plot path
cn.PLOT_DIR = os.path.join(cn.PLOT_DIR, "tests")



#############################
# Tests
#############################
class TestOptionManager(unittest.TestCase):

    def setUp(self):
        # Cannot modify self.control
        self.removeFiles()
        dct = dict(DCT)
        _, self.ax = plt.subplots(1)
        dct[cn.O_AX] = self.ax
        dct[cn.O_XLABEL] = "xlabel"
        dct[cn.O_IS_PLOT] = IS_PLOT
        self.option_mgr = OptionManager(dct)

    def tearDown(self):
        plt.close()
        self.removeFiles()

    def removeFiles(self):
        for ffile in os.listdir(cn.PLOT_DIR):
            if ("figure_" in ffile) and (".pdf") in ffile:
                path = os.path.join(cn.PLOT_DIR, ffile)
                if os.path.isfile(path):
                    os.remove(path)

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

    def testWritefig(self):
        # TESTING
        def test(writefig, idx=0):
            dct = dict(DCT)
            _, ax = plt.subplots(1)
            dct[cn.O_AX] = ax
            dct[cn.O_XLABEL] = "xlabel"
            dct[cn.O_WRITEFIG] = writefig
            dct[cn.O_IS_PLOT] = False
            option_mgr = OptionManager(dct)
            ax.plot(range(5), range(5))
            option_mgr.doPlotOpts()
            option_mgr.doFigOpts()
            #
            if isinstance(writefig, bool):
                filename = "figure_%d.pdf" % idx
                path = os.path.join(cn.PLOT_DIR, filename)
                is_file = writefig
            else:
                path = writefig
                is_file = True
            if is_file:
                self.assertTrue(os.path.isfile(path))
            else:
                self.assertFalse(os.path.isfile(path))

        test(False)
        test(True, idx=0)
        test(True, idx=1)
        test(False, idx=2)
        path = os.path.join(cn.PLOT_DIR, "figure_10.pdf")
        test(path)
        test(False, idx=2)


if __name__ == '__main__':
  unittest.main()

from controlSBML.control_extensions import time_plots

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False
if IS_PLOT:
    import matplotlib
    matplotlib.use('TkAgg')


#############################
# Tests
#############################
class TestTimePlots(unittest.TestCase):

    def setUp(self):
        pass

    def testPlotDF(self):
        if IGNORE_TEST:
          return
        df = pd.DataFrame({"a": range(10)})
        df["b"] = 10*df["a"]
        time_plots.plot1DF(df, ylabel="values", xlabel="sec", is_plot=IS_PLOT)



if __name__ == '__main__':
  unittest.main()

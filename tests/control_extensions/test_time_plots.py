from controlSBML.control_extensions import time_plots

import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = True
IS_PLOT = True
if IS_PLOT:
    import matplotlib
    matplotlib.use('TkAgg')
DF = pd.DataFrame({"a": range(10)})
DF["b"] = 10*DF["a"]


#############################
# Tests
#############################
class TestTimePlots(unittest.TestCase):

    def setUp(self):
        pass

    def testPlotOneDF(self):
        if IGNORE_TEST:
          return
        time_plots.plotOneDF(DF, ylabel="values", xlabel="sec",
              is_plot=IS_PLOT)

    def testPlotOneDF(self):
        # TESTING
        df = DF.applymap(lambda v: 100*v)
        time_plots.plotManyDF(DF, df, ylabel="values", xlabel="sec",
              is_plot=IS_PLOT, names=["first", "second"])



if __name__ == '__main__':
  unittest.main()

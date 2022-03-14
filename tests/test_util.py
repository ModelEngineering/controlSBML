from controlSBML import util
from controlSBML.timeseries import Timeseries
import controlSBML.constants as cn

import pandas as pd
import numpy as np
import tellurium as te
import unittest


IGNORE_TEST = False
IS_PLOT = False
SIZE = 10
if IS_PLOT:
    import matplotlib
    matplotlib.use('TkAgg')
times = [1.0*n for n in range(SIZE)]
TS = Timeseries(pd.DataFrame({"a": range(SIZE)}), times=times)
TS["b"] = 10*TS["a"]


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

    def setUp(self):
        pass

    def testCacluateMatrixDistance(self):
        if IGNORE_TEST:
            return
        SHAPE = (4, 4)
        SIZE = SHAPE[0]*SHAPE[1]
        mat1 = np.array(range(SIZE))
        mat1 = np.reshape(mat1, SHAPE)
        dist = util.calculateMatrixDistance(mat1, mat1)
        self.assertTrue(np.isclose(dist, 0))
        #
        mat2 = 4*mat1
        dist1 = util.calculateMatrixDistance(mat1, mat2)
        self.assertGreater(dist1, 0)
        dist2 = util.calculateMatrixDistance(mat2, mat1)
        self.assertTrue(np.isclose(dist1, dist2))

    def testGetModel(self):
        if IGNORE_TEST:
          return
        model_str = util.getModel()
        rr = te.loada(model_str)
        self.assertTrue("roadrunner" in str(type(rr)))
 
    def testPlotOneTS(self):
        if IGNORE_TEST:
          return
        util.plotOneTS(TS, ylabel="values", xlabel="sec",
              is_plot=IS_PLOT)

    def testPlotManyTS(self):
        if IGNORE_TEST:
          return
        df = TS.applymap(lambda v: 100*v)
        ts = Timeseries(df, times=df.index)
        util.plotManyTS(TS, ts, ylabel="values", xlabel="sec",
              is_plot=IS_PLOT, names=["first", "second"])
        util.plotManyTS(TS, ts, ylabel="values", xlabel="sec",
              is_plot=IS_PLOT, names=["first", "second"], ncol=2)


if __name__ == '__main__':
  unittest.main()

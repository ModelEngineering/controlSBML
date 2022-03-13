from controlSBML import util

import pandas as pd
import numpy as np
import tellurium as te
import unittest


IGNORE_TEST = False
IS_PLOT = False
if IS_PLOT:
    import matplotlib
    matplotlib.use('TkAgg')
DF = pd.DataFrame({"a": range(10)})
DF["b"] = 10*DF["a"]


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
        util.plotOneTS(DF, ylabel="values", xlabel="sec",
              is_plot=IS_PLOT)

    def testPlotOneTS(self):
        if IGNORE_TEST:
          return
        df = DF.applymap(lambda v: 100*v)
        util.plotManyTS(DF, df, ylabel="values", xlabel="sec",
              is_plot=IS_PLOT, names=["first", "second"])


if __name__ == '__main__':
  unittest.main()

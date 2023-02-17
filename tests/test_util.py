from controlSBML import util
from controlSBML.timeseries import Timeseries
import controlSBML.constants as cn
import controlSBML as ctl

import control
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
MDL = "A->B; 1; A=0; B=0; TOT:=A+B"
MDL_RR = te.loada(MDL)
NAMED_ARRAY = MDL_RR.simulate()
MAT = np.array(range(10))
MAT = np.reshape(MAT, (2, 5))
DF = pd.DataFrame(NAMED_ARRAY, columns=NAMED_ARRAY.colnames)


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

    def testMakeSimulationTimes(self):
        if IGNORE_TEST:
          return
        times = util.makeSimulationTimes()
        self.assertTrue(isinstance(times, np.ndarray))
        #
        times = util.makeSimulationTimes(start_time=1, end_time=4)
        self.assertTrue(times[0] == 1)
        self.assertTrue(np.isclose(times[-1], 4))
        #
        time1s = util.makeSimulationTimes(start_time=1, end_time=4,
            points_per_time=100)
        self.assertGreater(len(time1s), len(times))

    def testPpMat(self):
        if IGNORE_TEST:
          return
        mat = np.array(range(10))
        result1 = util.ppMat(mat, column_names=["a"], is_print=IS_PLOT)
        mat = np.reshape(mat, (5,2))
        result2 = util.ppMat(mat, column_names=["a", "b"], is_print=IS_PLOT)

    def testMat2DF(self):
        if IGNORE_TEST:
          return
        for mat in [NAMED_ARRAY, DF]:
            df = util.mat2DF(mat)
            self.assertTrue(isinstance(df, pd.DataFrame))
            self.assertTrue(any(["A" in c for c in df.columns]))
        df = util.mat2DF(MAT)
        self.assertTrue(isinstance(df, pd.DataFrame))

    def testPlotMat(self):
        if IGNORE_TEST:
          return
        for mat in [MAT, NAMED_ARRAY, DF]:
            util.plotMat(mat, title="test", figsize=(5,5), is_plot=IS_PLOT)

    def testTimeresponse2Timeseries(self):
        if IGNORE_TEST:
          return
        SIZE = 10
        factory = ctl.IOSystemFactory()
        sys = factory.makePassthru("sys")
        times = list(range(SIZE))
        timeresponse = control.input_output_response(sys, times, U=times)
        ts = util.timeresponse2Timeseries(timeresponse)
        self.assertEqual(len(ts), SIZE)

    def testSetRoadrunner(self):
        if IGNORE_TEST:
            return
        rr = te.loada(MDL)
        value = 100
        util.setRoadrunnerValue(rr, {"A": value})
        self.assertTrue(np.isclose(value, rr["A"]))
        #
        with self.assertWarns(Warning):
            util.setRoadrunnerValue(rr, {"TOT": value})




if __name__ == '__main__':
  unittest.main()

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
        return
        # FIXE: Feature is deprecated.
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

    def testSimplifyTransferFunction(self):
        if IGNORE_TEST:
            return
        def test(tf, expected_num, expected_den):
            new_tf = util.simplifyTransferFunction(tf)
            den = new_tf.den[0][0]
            num = new_tf.num[0][0]
            self.assertTrue(all([v == e for v, e in zip(num, expected_num)]))
            self.assertTrue(all([v == e for v, e in zip(den, expected_den)]))
        #
        tf = control.TransferFunction([1], [1, 0.0002, 1])
        test(tf, [1], [1, 0, 1])
        #
        tf = control.TransferFunction([1], [1, 1])
        test(tf, [1], [1, 1])
        #
        tf = control.TransferFunction([1], [1, 2, 0.0002])
        test(tf, [1], [1, 2])
        #
        tf = control.TransferFunction([1, 0.00001], [1, 2, 0.0002])
        test(tf, [1], [1, 2])
        #
        tf = control.TransferFunction([1, 0.001, 0.00001], [1, 2, 0.0002])
        test(tf, [1], [1, 2])
        #
        tf = control.TransferFunction([2, 1, 0.001, 0.00001], [4, 1, 0.00002, 0.0002])
        test(tf, [2, 1], [4, 1])
        #
        tf = control.TransferFunction([2, 0.001, 1, 0.00001], [4, 1, 0.00002, 3, 0.0002])
        test(tf, [2, 0, 1], [4, 1, 0, 3])

    def testLatexifyTransferFunction(self):
        if IGNORE_TEST:
           return
        tf = control.TransferFunction([1], [1, 0.0002, 1])
        latex = util.latexifyTransferFunction(tf)
        self.assertTrue(latex.count("$"))
        self.assertTrue(latex.count("frac"))
        #
        tf = control.TransferFunction([0], [1, 0.0002, 1])
        latex = util.latexifyTransferFunction(tf)
        self.assertFalse("frac" in latex)
        #
        tf = control.TransferFunction([1, 0, 0], [1, 0.0002, 1])
        latex = util.latexifyTransferFunction(tf)

    def testCleanTimes(self):
        if IGNORE_TEST:
           return
        expected_diff = 1000
        TIMES = expected_diff*np.array([1, 2, 3, 4, 5]).astype(int)
        def test(times):
            result = util.cleanTimes(times)
            mean_diff = np.mean(np.diff(result))
            self.assertTrue(np.isclose(mean_diff, expected_diff))
            self.assertTrue(all([r == t for r, t in zip(result, TIMES)]))
        #
        test(TIMES)
        new_times = np.array(TIMES)
        new_times[2] += 2
        test(new_times)

    def testCalculateInitialValue(self):
        if IGNORE_TEST:
           return
        times = np.linspace(0, 10, 100)
        def test(tf):
            value = util.calculateInitialValue(tf)
            _, ys = control.step_response(tf, T=times)
            self.assertTrue(np.isclose(value, ys[0]))
        #
        tf = control.TransferFunction([1], [1, 0.0002, 1])
        test(tf)
        tf = control.TransferFunction([2, 2, 3], [1, 0.0002, 1])
        test(tf)

    def testCompareSingleArgumentFunctions(self):
        if IGNORE_TEST:
           return
        func1 = lambda x: x**2
        func2 = lambda x: x
        result = util.compareSingleArgumentFunctions(func1, func2, 0, 100)
        self.assertFalse(result)
        result = util.compareSingleArgumentFunctions(func1, func1, 0, 100)
        self.assertTrue(result)
        result = util.compareSingleArgumentFunctions(func2, func2, 0, 100)
        self.assertTrue(result)

    def testIsPositiveRealPart(self):
        if IGNORE_TEST:
           return
        def test(arr, is_stable):
            result = util.isPositiveRealPart(arr)
            self.assertEqual(result, is_stable)
        #
        test(np.array([-1+0j, 1+2j]), True)
        test(np.array([-1+0j, -1+2j]), False)
    
    def testIsStablePoles(self):
        if IGNORE_TEST:
           return
        def test(tf, is_stable):
            result = util.isStablePoles(tf)
            self.assertEqual(result, is_stable)
        #
        tf = control.TransferFunction([1, 1], [1, -1])
        test(tf, False)
        tf = control.TransferFunction([1, -1], [1, 1])
        test(tf, True)

    def testIsStableZeros(self):
        if IGNORE_TEST:
           return
        def test(tf, is_stable):
            result = util.isStableZeros(tf)
            self.assertEqual(result, is_stable)
        #
        tf = control.TransferFunction([1, 1], [1, -1])
        test(tf, True)
        tf = control.TransferFunction([1, -1], [1, 1])
        test(tf, False)

    def testRoundToDigits(self):
        if IGNORE_TEST:
           return
        def test(value, digits, expected):
            result = util.roundToDigits(value, digits)
            self.assertTrue(np.isclose(result, expected))
        #
        test(0.000223, 2, 0.00022)
        test(0.22, 1, 0.2)
        test(22, 1, 22)
        

if __name__ == '__main__':
    unittest.main()

from controlSBML.timeseries import TimeSeries
import controlSBML.constants as cn

import numpy as np
import pandas as pd
import unittest
import tellurium as te


IGNORE_TEST = False
IS_PLOT = False
if IS_PLOT:
    import matplotlib
    matplotlib.use('TkAgg')
EFFECTOR_DCT = {"J0": "E_J0"}
END_TIME = 5
COLUMNS = ["[a]", "b"]
SIZE = 10
DF = pd.DataFrame([range(SIZE), range(SIZE)]).transpose()
DF.columns = COLUMNS
TIMES = 10.0*np.array(range(SIZE))
DF.index = TIMES
MAT = DF.values
LINEAR_MDL = """
J0: $S0 -> S1; $S0
J1: S1 -> S2; S1
J2: S2 -> S3; S2

$S0 = 0
S1 = 10
S2 = 0
S3 = 0
"""
rr = te.loada(LINEAR_MDL)
NAMED_ARRAY = rr.simulate()


#############################
# Tests
#############################
class TestTimeseries(unittest.TestCase):

    def setUp(self):
        pass

    def _validate(self, ts):
        self.assertGreater(len(ts), 0)
        falses = [("[" in c) or ("]" in c) for c in ts.columns]
        self.assertFalse(any(falses))
        trues = [isinstance(v, int) for v in ts.index]
        self.assertTrue("a" in ts.columns)
        self.assertTrue(all(trues))
        self.assertFalse(np.isnan(np.sum(ts.values.flatten())))

    def testConstructorTS(self):
        if IGNORE_TEST:
          return
        ts = TimeSeries(MAT, times=TIMES, columns=COLUMNS)
        new_ts = TimeSeries(ts)
        self._validate(new_ts)

    def testConstructorMat(self):
        if IGNORE_TEST:
            return
        ts = TimeSeries(MAT, times=TIMES, columns=COLUMNS)
        self._validate(ts)

    def testConstructorNamedArray(self):
        if IGNORE_TEST:
          return
        ts = TimeSeries(NAMED_ARRAY)
        columns = list(ts.columns)
        columns[0] = "a"
        ts.columns = columns
        self._validate(ts)

    def testConstructorDFTimeColumn(self):
        if IGNORE_TEST:
          return
        df = DF.copy()
        df[cn.TIME] = df.index
        df.index = range(SIZE)
        ts = TimeSeries(df)
        self._validate(ts)

    def testConstructorDFNoTimeColumn(self):
        if IGNORE_TEST:
          return
        df = DF.copy()
        df.index = range(SIZE)
        times = list(DF.index)
        ts = TimeSeries(df, times=times)
        self._validate(ts)

    def testGetItem(self):
        if IGNORE_TEST:
          return
        ts = TimeSeries(MAT, times=TIMES, columns=COLUMNS)
        new_ts = ts["a"]
        self._validate(ts)

    def testGetItem1(self):
        if IGNORE_TEST:
          return
        ts = TimeSeries(MAT, times=TIMES, columns=COLUMNS)
        ts["c"] = 10*ts["a"]
        ts["d"] = range(SIZE)
        self._validate(ts)

    def testGetItem2(self):
        if IGNORE_TEST:
          return
        ts = TimeSeries(MAT, times=TIMES, columns=COLUMNS)
        new_ts = ts[["a", "b"]]
        self._validate(new_ts)


if __name__ == '__main__':
  unittest.main()

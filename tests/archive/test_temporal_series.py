from controlSBML.temporal_series import TemporalSeries

import numpy as np
import unittest


IGNORE_TEST = False
IS_PLOT = False
TIMES = [0, 1, 1.1, 2, 2.03, 2.5, 3]
VALUES = list(TIMES)


#############################
# Tests
#############################
class TestTemporalSeries(unittest.TestCase):

    def setUp(self):
        # Cannot modify self.control
        self.series = TemporalSeries(times=TIMES, values=VALUES)

    def testConstructor(self):
        if IGNORE_TEST:
          return
        self.assertTrue(np.allclose(self.series._times, TIMES))
        self.assertTrue(np.allclose(self.series._values, VALUES))

    def testLen(self):
        if IGNORE_TEST:
          return
        self.assertEqual(len(self.series), len(TIMES))
        self.series.add(10, 10)
        self.assertEqual(len(self.series), len(TIMES)+1)

    def testGetValue(self):
        if IGNORE_TEST:
          return
        # Out of range
        with self.assertRaises(ValueError):
            self.series.getValue(7.1)
        # Exact match
        self.assertEqual(self.series.getValue(1), 1)
        self.assertEqual(self.series.getValue(1.1), 1.1)
        self.assertEqual(self.series.getValue(2), 2)
        self.assertEqual(self.series.getValue(2.5), 2.5)
        # Interpolate
        self.assertEqual(self.series.getValue(1.05), 1.05)

    def testGetTimeValues(self):
        if IGNORE_TEST:
          return
        count = 5
        times = np.linspace(0, 5, 50)
        series = TemporalSeries(times=times, values=times)
        times, values = series.getTimesValues(count)
        self.assertEqual(len(values), count)
        self.assertEqual(len(times), count)
        self.assertTrue(np.allclose(times, values))
        variance = np.var(np.diff(times))
        self.assertTrue(np.isclose(variance, 0))



if __name__ == '__main__':
  unittest.main()
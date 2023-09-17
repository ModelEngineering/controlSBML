from controlSBML.temporal_series import TemporalSeries

import numpy as np
import unittest


IGNORE_TEST = False
IS_PLOT = False
TIMES = [1, 1.1, 2, 2.03, 2.5, 3]
VALUES = list(TIMES)


#############################
# Tests
#############################
class TestTemporalSeries(unittest.TestCase):

    def setUp(self):
        # Cannot modify self.control
        if IGNORE_TEST:
          return
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



if __name__ == '__main__':
  unittest.main()
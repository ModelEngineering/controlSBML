from controlSBML import logger as lg

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import unittest


IGNORE_TEST = False
IS_PLOT = False
SIZE = 10
LOGGER_NAME = "logger"
LOGGER2_NAME = "logger2"
ITEM_NAMES = ["in", "out"]
        

#############################
# Tests
#############################
class TestLogger(unittest.TestCase):

    def setUp(self):
        self.logger = lg.Logger(LOGGER_NAME, ITEM_NAMES)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertEqual(len(self.logger.dct[lg.TIME]), 0)

    def testAdd(self):
        if IGNORE_TEST:
          return
        item_values = np.repeat(1, len(ITEM_NAMES))
        self.logger.add(1, item_values)
        trues = [v[0] == 1 for v in self.logger.dct.values()]
        self.assertTrue(all(trues))

    def testReport(self):
        if IGNORE_TEST:
          return
        size = 10
        [self.logger.add(n, np.repeat(n, len(ITEM_NAMES))) for n in range(size)]
        df = self.logger.report()
        self.assertEqual(len(df.columns), len(ITEM_NAMES))
        self.assertEqual(len(df), size)
        self.assertEqual(df.index.name, lg.TIME)

    def testMerge(self):
        if IGNORE_TEST:
          return
        logger = lg.Logger(LOGGER2_NAME, ITEM_NAMES)
        size = 10
        [self.logger.add(n, np.repeat(n, len(ITEM_NAMES))) for n in range(size)]
        [logger.add(n, np.repeat(n, len(ITEM_NAMES))) for n in range(size)]
        merged_logger = self.logger.merge("new_logger", [logger])
        df = merged_logger.report()
        self.assertEqual(len(df.columns), 2*len(ITEM_NAMES))
        self.assertEqual(len(df), size)
        self.assertEqual(df.index.name, lg.TIME)

   
if __name__ == '__main__':
  unittest.main()

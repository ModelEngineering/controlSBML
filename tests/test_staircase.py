from controlSBML.staircase import Staircase
import helpers

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False
INITIAL_VALUE = 5
FINAL_VALUE = 10
NUM_STEP = 5
NUM_POINT = 100
PLOT_PATH = helpers.setupPlotting(__file__)
    

class TestStaircase(unittest.TestCase):

    def setUp(self):
         self.remove()
         if False:
             self.staircase = Staircase.makeRelativeStaircase(
                 center=10, fractional_deviation=0.1, num_step=5, num_point=100)
         
    def tearDown(self):
        self.remove()
         
    def remove(self):
        if os.path.isfile(PLOT_PATH):
            os.remove(PLOT_PATH)

    def testConstructor(self):
        if IGNORE_TEST:
            return

        def test(num_point, num_step, initial_value=0, final_value=5):
            staircase = Staircase(initial_value=initial_value,
                                    final_value=final_value,
                                    num_step=num_step,
                                    num_point=num_point)
            self.assertTrue(len(staircase.staircase_arr), num_point)
            num_distinct = len(set(staircase.staircase_arr))
            self.assertEqual(num_distinct, num_step + 1)
            self.assertEqual(staircase.staircase_arr[0], initial_value)
            self.assertEqual(staircase.staircase_arr[-1], final_value)
            _ = staircase.plot()
            plt.savefig(PLOT_PATH)
        #
        test(20, 4)
        test(19, 4)
        test(191, 17)
        test(91, 15)

    def testMakeRealtive(self):
        #if IGNORE_TEST:
        #    return
        staircase = Staircase.makeRelativeStaircase(5, 1.0)
        self.assertTrue(isinstance(staircase, Staircase))

    def testMakeRealtive(self):
        if IGNORE_TEST:
            return
        pass


if __name__ == '__main__':
    unittest.main()
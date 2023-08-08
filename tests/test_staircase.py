from controlSBML.staircase import Staircase
import helpers

import matplotlib.pyplot as plt
import numpy as np
import os
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
            if not IS_PLOT: 
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
        if IGNORE_TEST:
            return
        center = 5
        fractional_deviation = 2.0
        num_point = 100
        deviation = center*fractional_deviation
        staircase = Staircase.makeRelativeStaircase(center=center,
                                                    fractional_deviation=fractional_deviation,
                                                    num_step=10, num_point=num_point)
        self.assertEqual(len(staircase.staircase_arr), num_point)
        self.assertTrue(np.isclose(staircase.staircase_arr[-1], center + deviation))
        self.assertTrue(np.isclose(staircase.staircase_arr[0], center - deviation))
        self.assertTrue(isinstance(staircase, Staircase))

    def testPlot(self):
        if IGNORE_TEST:
            return
        staircase = Staircase.makeRelativeStaircase(5, 2)
        staircase.plot()
        plt.savefig(PLOT_PATH)


if __name__ == '__main__':
    unittest.main()
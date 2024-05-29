from controlSBML.parallel_coordinates import ParallelCoordinates  # type: ignore

import pandas as pd # type: ignore
import matplotlib.pyplot as plt # type: ignore
import numpy as np
import unittest


IGNORE_TEST = False
IS_PLOT = False
# Create data
columns = ['mpg', 'displacement', 'cylinders', 'horsepower', 'weight', 'acceleration']
dct = {}
factor = 1
# Make it easy to distinguish the columns
for column in columns:
    factor *= 10
    dct[column] = factor*np.random.rand(10)
DF = pd.DataFrame(dct)
NAN_IDXS = [2, 9]
for idx in NAN_IDXS:
    DF.loc[idx, 'mpg'] = np.nan
DF1 = pd.read_csv("tests/plot_parallel_coordinates.csv")
del DF1['Unnamed: 0']
NUM_CATEGORY = 5


#############################
# Tests
#############################
class TestParallelCoordinates(unittest.TestCase):

    def setUp(self):
        self.df = DF.copy()
        self.parallel = ParallelCoordinates(self.df, value_column='mpg', num_category=NUM_CATEGORY,
                                   is_plot=IS_PLOT, round_digit=4)
        
    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue("ParallelCoordinates" in str(type(self.parallel)))
        for column in self.parallel.columns:
            if column == "mpg":
                continue
            ser = self.parallel.scaled_df[column]
            self.assertTrue(column in self.parallel.columns)
            self.assertLessEqual(ser.min(), ser.max())

    def testSetAxisTicks(self):
        if IGNORE_TEST:
            return
        _, ax = plt.subplots(1)
        self.parallel._setAxisTicks("displacement", ax)

    def testMakeCategoriesAndLabels(self):
        if IGNORE_TEST:
            return
        categories, untrimmed_labels, trimmed_labels = self.parallel._makeCategoriesAndLabels()
        for idx in NAN_IDXS:
            self.assertEqual(categories[idx], NUM_CATEGORY-1)
        self.assertEqual(untrimmed_labels[-1], "nan")
        self.assertEqual(trimmed_labels[-1], "nan")

    def testPlot(self):
        if IGNORE_TEST:
            return
        self.parallel.plot()

    def testPlotParallelCoordinates(self):
        if IGNORE_TEST:
            return
        ParallelCoordinates.plotParallelCoordinates(self.df, value_column='mpg', num_category=NUM_CATEGORY,
                                   is_plot=IS_PLOT, round_digit=4)

    def testPlotParallelCoordinates1(self):
        #if IGNORE_TEST:
        #    return
        ParallelCoordinates.plotParallelCoordinates(DF1, value_column='score', num_category=NUM_CATEGORY,
                                   round_digit=4, is_plot=IS_PLOT)

        

if __name__ == '__main__':
    unittest.main()
from controlSBML.plot_parallel_coordinates import plotParallelCoordinates  # type: ignore

import pandas as pd # type: ignore
import numpy as np
import unittest


IGNORE_TEST = False
IS_PLOT = False
# Create data
columns = ['mpg', 'displacement', 'cylinders', 'horsepower', 'weight', 'acceleration']
dct = {}
for column in columns:
    dct[column] = np.random.rand(10)
DF = 10*pd.DataFrame(dct)
DF.loc[9, 'mpg'] = np.nan
DF.loc[2, 'mpg'] = np.nan
DF1 = pd.read_csv("tests/plot_parallel_coordinates.csv")
del DF1['Unnamed: 0']


#############################
# Tests
#############################
class TestFunction(unittest.TestCase):

    def setUp(self):
        self.df = DF.copy()

    def testPlotParallelCoordinates(self):
        if IGNORE_TEST:
            return
        plotParallelCoordinates(self.df, value_column='mpg', num_category=3,
                                   is_plot=IS_PLOT)
        
    def testPlotParallelCoordinates1(self):
        if IGNORE_TEST:
            return
        columns = ["kD", "kP", "kI", "kF"]
        plotParallelCoordinates(DF1, value_column='score', num_category=10, columns=columns,
                                   is_plot=IS_PLOT, round_digit=4)

        

if __name__ == '__main__':
    unittest.main()

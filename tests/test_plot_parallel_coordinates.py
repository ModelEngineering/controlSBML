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

        

if __name__ == '__main__':
    unittest.main()

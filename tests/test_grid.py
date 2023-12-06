from controlSBML.grid import Grid, Axis
from controlSBML.timeseries import Timeseries

import numpy as np
import unittest


IGNORE_TEST = False
IS_PLOT = False
NUM_COORDINATE = 3


#############################
# Tests
#############################
class TestAxis(unittest.TestCase):

    def setUp(self):
        self.axis = Axis("test", min_value=0, max_value=10, num_coordinate=NUM_COORDINATE, is_random=False)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.axis, Axis))

    def testGetGridCoordinates(self):
        if IGNORE_TEST:
            return
        expected = [0, 5, 10]
        actual = self.axis.coordinates
        self.assertEqual(expected, actual)

    def testGetPointCoordinate(self):
        if IGNORE_TEST:
            return
        min_value = 0
        max_value = 10
        axis = Axis("test", min_value=min_value, max_value=max_value, num_coordinate=2, is_random=False)
        actual = axis.getPointCoordinate(0)
        self.assertTrue(np.isclose(5, actual))
        #
        axis = Axis("test", min_value=min_value, max_value=max_value, num_coordinate=2, is_random=True)
        actual = axis.getPointCoordinate(0)
        self.assertGreater(actual, min_value)
        self.assertLess(actual, max_value)


class TestGrid(unittest.TestCase):

    def setUp(self):
        self.grid = Grid()

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.grid, Grid))

    def testAddAxis(self):
        if IGNORE_TEST:
            return
        parameter_name = "test"
        self.grid.addAxis(parameter_name, min_value=0, max_value=10, num_coordinate=NUM_COORDINATE)
        self.assertTrue(parameter_name in self.grid.axis_dct.keys())
        self.assertTrue(isinstance(self.grid.axis_dct[parameter_name], Axis))
        #
        with self.assertRaises(ValueError):
            self.grid.addAxis(parameter_name, min_value=1000, max_value=10, num_coordinate=NUM_COORDINATE)

    def test_num_point(self):
        if IGNORE_TEST:
            return
        self.grid.addAxis("test1", min_value=0, max_value=10, num_coordinate=NUM_COORDINATE)
        self.grid.addAxis("test2", min_value=0, max_value=10, num_coordinate=NUM_COORDINATE)
        self.assertEqual(self.grid.num_point, (NUM_COORDINATE-1)**2)

    def testRepr(self):
        if IGNORE_TEST:
            return
        self.grid.addAxis("test1", min_value=0, max_value=10, num_coordinate=NUM_COORDINATE)
        self.grid.addAxis("test2", min_value=0, max_value=10, num_coordinate=NUM_COORDINATE)
        self.assertEqual(str(self.grid).count("test"), 2)

    def testIteratePoint(self):
        if IGNORE_TEST:
            return
        self.grid.addAxis("test1", min_value=0, max_value=10, num_coordinate=5)
        self.grid.addAxis("test2", min_value=100, max_value=200, num_coordinate=NUM_COORDINATE)
        points = list(self.grid._iteratePoints())
        self.assertEqual(len(points), self.grid.num_point)

    def testPlot(self):
        if IGNORE_TEST:
            return
        for is_random in [True, False]:
            grid = Grid(is_random=is_random)
            grid.addAxis("parm1", min_value=0, max_value=10, num_coordinate=10)
            grid.addAxis("parm2", min_value=100, max_value=200, num_coordinate=5)
            grid.plot(is_plot=IS_PLOT, figsize=(5,5))
            #
            axis = grid.getAxis("parm1")
            axis.setMinValue(150)
            axis.setMaxValue(160)
            axis.setNumCoordinate(10)
            grid.plot(is_plot=IS_PLOT, figsize=(5,5))

if __name__ == '__main__':
    unittest.main()
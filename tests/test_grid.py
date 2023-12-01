from controlSBML.grid import Grid, Axis

import unittest

from controlSBML import util
from controlSBML.timeseries import Timeseries
import controlSBML.constants as cn
import controlSBML as ctl

import control
import pandas as pd
import numpy as np
import tellurium as te
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

    def testGetPointCoordinates(self):
        if IGNORE_TEST:
            return
        expected = [2.5, 7.5]
        axis = Axis("test", min_value=0, max_value=10, num_coordinate=NUM_COORDINATE, is_random=False)
        actual = axis.getPointCoordinates()
        self.assertEqual(expected, actual)
        #
        axis = Axis("test", min_value=0, max_value=10, num_coordinate=NUM_COORDINATE, is_random=True)
        coordinates = axis.coordinates
        actual = axis.getPointCoordinates()
        for idx in range(len(actual)):
            self.assertTrue(actual[idx] >= coordinates[idx])
            self.assertTrue(actual[idx] <= coordinates[idx + 1])


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
        points = list(self.grid.iteratePoints())
        self.assertEqual(len(points), self.grid.num_point)

    def testPlot(self):
        if IGNORE_TEST:
            return
        self.grid.addAxis("parm1", min_value=0, max_value=10, num_coordinate=10)
        self.grid.addAxis("parm2", min_value=100, max_value=200, num_coordinate=5)
        self.grid.plot(is_plot=IS_PLOT, figsize=(5,5))
        #
        axis = self.grid.getAxis("parm1")
        axis.min_value = 150
        axis.max_value = 160
        axis.num_coordinate = 10
        self.grid.generatePoints()
        self.grid.plot(is_plot=IS_PLOT, figsize=(5,5))

if __name__ == '__main__':
    unittest.main()
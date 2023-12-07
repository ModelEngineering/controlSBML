"""
Creates an n-dimensional Grid for exploring the values of n parameters.

A grid is a collection of n-dimensional cubes called grid elements.


     10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        x                x                 x
        x     G1         x       G2        x
        x            p   x                 x
        x                x                 x
P1   5  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        x                x                 x
        x     G3         x       G4        x
        x                x                 x
        x                x                 x
     0  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        0                10               20
                       P2

Above is a grid with the grid elements G1-G4.
P1 and P2 are parameters. The range of values for P1 is [0, 10]; this is referred to as the P1 axis. The P2 axis is [0, 20].
Coordinates are values on an axis. G1 is specified by the coordinaes: P1: [5, 10]; P2: [0, 10].
A point in the grid space is specified by a coordinate value for each parameter. The point "p" in G1 is P1: 7.5; P2: 5.

Usage:
    grid = Grid()
    grid.addAxis("P1", min_value=0, max_value=10, num_coordinate=3)
    grid.addAxis("P2", min_value=0, max_value=20, num_coordinate=3)
    print(grid_num_grid)  # 4
    for point in grid.pointIterator():
        print(point)   # Calculate using the point values

    # Change the coordinates for an axis
    axis = grid.getAxis("P1")
    axis.setMinValue(2)
    axis.setMaxValue(8)
    axis.setNumCoordinate(25)
"""
from controlSBML.option_management.option_manager import OptionManager
import controlSBML.constants as cn

import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typing

DEFAULT_MIN = 0
DEFAULT_MAX = 10
DEFAULT_NUM_COORDINATE = 3
NAME = "name"

class Point(dict):
    def __init__(self, **kwargs):
        """
        Specifies a point in a grid.

        Args:
            kwargs: dict (key: parameter name; value: float)
        """
        self.update(kwargs)

    def __repr__(self)->str:
        return "Point({})".format(", ".join(["{}: {}".format(k, v) for k, v in self.items()]))


class Axis:
    def __init__(self, parameter_name:str, min_value:float=None, max_value:float=None, num_coordinate:int=None,
                 is_random:bool=False, notifier=None):
        """
        Specifies the axis for a parameter.

        Args:
            parameter_parameter_name: str (parameter_name of parameter)
            min_value: float (minimum value for parameter)
            max_value: float (maximum value for parameter)
            num_coordinate: int (includes the first and last coordinate and so this is one more than the number of grid elements)
            is_random: bool (True: random points; False: use midpoint of grid element)
            notifier: function (called when an axis value is changed: no arguments, no return value)
        """
        self.parameter_name = parameter_name
        self.min_value = min_value
        self.max_value = max_value
        self.num_coordinate = num_coordinate
        self.is_random = is_random
        if notifier is None:
            notifier = lambda: None
        self.notifier = notifier

    def __repr__(self)->str:
        return "Axis({}, min_value={}, max_value={}, num_coordinate={})".format(self.parameter_name,
                                                                                self.min_value, self.max_value, self.num_coordinate)
    ######## Setters ########
    def setMinValue(self, min_value:float):
        """
        Sets the minimum value for the axis.

        Args:
            min_value: float (minimum value for parameter)
        """
        self.min_value = min_value
        self.notifier()

    def setMaxValue(self, max_value:float):
        """
        Sets the maximum value for the axis.

        Args:
            max_value: float (maximum value for parameter)
        """
        self.max_value = max_value
        self.notifier()

    def setNumCoordinate(self, num_coordinate:int):
        """
        Sets the number of coordinates for the axis.

        Args:
            num_coordinate: int (number of point coordinates for parameter)
        """
        self.num_coordinate = num_coordinate
        self.notifier() 

    @property
    def coordinates(self):
        """
        Finds the boundary coordinates for grids for this Axis

        Returns:
            _type_: _description_
        """
        return np.linspace(self.min_value, self.max_value, self.num_coordinate).tolist()
    
    def getPointCoordinate(self, idx)->float:
        """
        Finds the coordinate for a point on this axis.
        Args:
            idx: int (index of coordinate)

        Returns:
            float
        """
        if idx < 0 or idx >= self.num_coordinate - 1:
            raise ValueError("Invalid index: {}".format(idx))
        coordinates = self.coordinates
        if self.is_random:
            return np.random.uniform(coordinates[idx], coordinates[idx + 1])
        else:
            return (coordinates[idx] + coordinates[idx + 1])/2
   

class Grid(object):
    def __init__(self, axis_dct=None, min_value=DEFAULT_MIN, max_value=DEFAULT_MAX,
                 num_coordinate=DEFAULT_NUM_COORDINATE, is_random=True):
        """
        Creates a dictionary of values to use for iterations for a grid search of a parameter space.

        Args:
            axis_dict: dict (species the axis for each parameter)
                key: name of parameter
                value: AxisSpecifier
            default_min: float/dict (parameter name: value)
            default_max: float/dict (parameter name: value)
            default_num_coordinate: int (number of coordinates for a parameter; minimum is 2)
            is_random: bool (True: random points; False: use midpoint of grid element)
        """
        if axis_dct is None:
            axis_dct = {}
        self.axis_dct = axis_dct
        self.default_min = min_value
        self.default_max = max_value
        if num_coordinate < 2:
            raise ValueError("default_num_coordinate must be at least 2.")
        self.default_num_coordinate = num_coordinate
        self.is_random = is_random
        # Updated when an axis is added
        #   key: parameter name
        #   value: list of coordinates corresponding to points for the parameter
        self._points = None

    def notifyAxisChange(self):
        """
        Notifies the grid that an axis has changed.
        """
        self.recalculatePoints()

    def recalculatePoints(self):
        """
        Recalculates the points.
        """
        self._points = None

    def __repr__(self)->str:
        dct = {NAME: [], cn.MIN: [], cn.MAX: [], "num_coordinate": []}
        for parameter_name, axis in self.axis_dct.items():
            dct[NAME].append(parameter_name)
            dct[cn.MIN].append(axis.min_value)
            dct[cn.MAX].append(axis.max_value)
            dct["num_coordinate"].append(axis.num_coordinate)
        return str(pd.DataFrame(dct))

    def addAxis(self, parameter_name:str, min_value:float=None,
                          max_value:float=None, num_coordinate:int=None):
        """
        Creates an axis for a parameter.
        Args:
            parameter_name: str (name of parameter)
            min_value: float (minimum value for parameter)
            max_value: float (maximum value for parameter)
            num_coordinate: int (number of point coordinates for parameter)
        """
        if min_value is None:
            min_value = self.default_min
        if max_value is None:
            max_value = self.default_max
        if num_coordinate is None:
            num_coordinate = self.default_num_coordinate
        if parameter_name in self.axis_dct:
            raise ValueError("Parameter name already exists: {}".format(parameter_name))
        if min_value >= max_value:
            raise ValueError("min_value must be less than max_value: {} >= {}".format(min_value, max_value))
        self.axis_dct[parameter_name] = Axis(parameter_name, min_value=min_value, max_value=max_value,
                                             num_coordinate=num_coordinate, is_random=self.is_random,
                                             notifier=self.notifyAxisChange)
        self.recalculatePoints()
        
    @property
    def num_point(self)->int:
        """
        Finds the number of points in the grid.

        Returns:
            int: number of points
        """
        return len(self.points)
    
    @property
    def points(self)->typing.List[Point]:
        """
        Gets the list of points.

        Returns:
            list-dict: list of points
        """
        if self._points is None:
            self._points = list(self._iteratePoints())
        return self._points

    def _iteratePoints(self):
        """
        Creates an iterator that returns a point.

        Returns:
            Point
        """
        index_lists = [list(range(self.axis_dct[p].num_coordinate - 1)) for p in self.axis_dct.keys()]
        for indices in itertools.product(*index_lists):
            yield Point(**{a.parameter_name: a.getPointCoordinate(i) for a, i in zip(self.axis_dct.values(), indices)})

    def getAxis(self, parameter_name:str)->Axis:
        """
        Gets the Axis for a parameter.

        Args:
            parameter_name: str (name of parameter)

        Returns:
            Axis
        """
        if not parameter_name in self.axis_dct:
            raise ValueError("Parameter name not found: {}".format(parameter_name)) 
        return self.axis_dct[parameter_name]
    
    def plot(self, **kwargs):
        """
        Plots a two dimensional grid.

        Args:
            kwargs: dict: keyword plot arguments
        """
        XAXIS = 0
        YAXIS = 1
        if len(self.axis_dct) != 2:
            raise ValueError("Grid must have two dimensions.")
        mgr = OptionManager(kwargs)
        ax = mgr.getAx()
        ##_, ax = plt.subplots(1)
        # Plot the grid elements
        axes = list(self.axis_dct.values())
        # Do the x-axis
        x_coords = axes[XAXIS].coordinates
        yv = [axes[YAXIS].min_value, axes[YAXIS].max_value]
        for val in x_coords:
            xv = [val, val]
            ax.plot(xv, yv, linestyle="--", color="blue")
        # Do the y-axis
        y_coords = axes[YAXIS].coordinates
        xv = [axes[XAXIS].min_value, axes[XAXIS].max_value]
        for val in y_coords:
            yv = [val, val]
            ax.plot(xv, yv, linestyle="--", color="blue")
        # Do the points
        xv = []
        yv = []
        for point in self.points:
            xv.append(point[axes[XAXIS].parameter_name])
            yv.append(point[axes[YAXIS].parameter_name])
        ax.scatter(xv, yv, marker="o", color="red")
        ax.set_xlabel(axes[XAXIS].parameter_name)
        ax.set_ylabel(axes[YAXIS].parameter_name)
        kwargs.setdefault("is_plot", True)
        mgr.doPlotOpts()
        mgr.doFigOpts()
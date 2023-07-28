"""Create a staircase of values. Build multiple staircases."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DEFAULT_NUM_STEP = 5
DEFAULT_INITIAL_VALUE = 0
DEFAULT_FINAL_VALUE = 10
DEFAULT_POINT_PER_STEP = 10
DEFAULT_NUM_POINT = 100
C_NUM_POINT = "num_point"
C_NUM_STEP = "num_step"


class Staircase(object):

    def __init__(self, initial_value=DEFAULT_INITIAL_VALUE,
                 num_step=DEFAULT_NUM_STEP, final_value=DEFAULT_FINAL_VALUE,
                 num_point=DEFAULT_NUM_POINT, name=None):
        """
        Creates a staircase of values.

        Parameters
        ----------
        initial_value: float (value for first step)
        final_value: float (value for final step)
        num_step: int (number of steps in staircase)
        name: label for staircase
        """
        self.initial_value = initial_value
        self.num_step = num_step
        self.final_value = final_value
        self.num_point = num_point
        self.name = name
        self._updateState()
        #self.staircase_arr is dynamically

    def copy(self):
        return Staircase(initial_value=self.initial_value,
                         num_step=self.num_step,
                         final_value=self.final_value,
                         num_point=self.num_point,
                         name=self.name)

    def _updateState(self):
        self.num_level = self.num_step + 1
        self.point_per_level = int(self.num_point/self.num_level)

    @property
    def staircase_arr(self):
        self._updateState()
        #
        steps = []  # Steps in the staircase
        #
        for num in range(self.num_level):
            steps.extend(list(np.repeat(num, self.point_per_level)))
        staircase_arr = np.array(steps)
        # Rescale
        staircase_arr = staircase_arr*(self.final_value - self.initial_value)/(self.num_level - 1)
        staircase_arr += self.initial_value
        # Add more if needed
        num_add_more = self.num_point - len(staircase_arr)
        if num_add_more < 0:
            import pdb; pdb.set_trace()
            raise RuntimeError("Negative residual count")
        elif num_add_more > 0:
            final_values = np.repeat(self.final_value, num_add_more)
            staircase_arr = np.append(staircase_arr, final_values)
        #
        return staircase_arr
        

    @classmethod
    def makeRelativeStaircase(cls, center=5, fractional_deviation=0.1,
                              num_step=DEFAULT_NUM_STEP, num_point=DEFAULT_NUM_POINT):
        """
        Specifies the staircase relative to a center point.

        Parameters
        ----------
        center: float (center point)
        fractional_deviation: float (fractional deviation from center)
        num_step: int (number of steps in staircase)
        num_point: int (number of points in staircase)

        Returns
        -------
        Staircase
        """
        max_deviation = center*fractional_deviation
        cls(initial_value=center - max_deviation,
            final_value=center + max_deviation,
            num_step=num_step, num_point=num_point)
        
    def plot(self):
        """
        Plots the staircase.
        """
        staircase_ser = pd.Series(self.staircase_arr)
        fig, ax = plt.subplots()
        staircase_ser.plot(ax=ax)
        return fig


class MultiStaircase(object):
    """
    Builds and provides access to multiple staircases with the same number of steps and points.

    Usage:
        multi_starcases = MultiStaircase(num_point=20, num_step=3)
        multi_starcases.buildAbsolute(["S1", "S2"], initial_value=0, final_value=10)
        multi_starcases.buildRelative(["S4", "S5"], 5, 0.1)
        S1_staircase = multi_starcases.get("S1")
    """

    def __init__(self, **kwargs):
        """
        Args:
            kwargs: dict (default arguments for Staircase)
        """
        self.kwargs = kwargs
        #
        self.staircase_dct = {}  # key is name; value is Staircase

    def _makeKwargs(self, kwargs):
        new_kwargs = dict(kwargs)
        for key, value in self.kwargs.items():
            if not key in new_kwargs:
                new_kwargs[key] = value
        return new_kwargs

    def buildAbsolute(self, names, **kwargs):
        """
        Creates a staircase for all of the names.
        Args:
            names: list-str
            kwargs: keyword arguments for Staircase constructor
        """
        new_kwargs = self._makeKwargs(kwargs)
        for name in names:
            self.staircase_dct[name] = Staircase(**new_kwargs)

    def buildRelative(self, names, **kwargs):
        """
        Args:
            names: list-str
            kwargs: keyword arguments for makeRelativeStaircase 
        """
        new_kwargs = self._makeKwargs(kwargs)
        for name in names:
            self.staircase_dct[name] = Staircase.makeRelativeStaircase(**new_kwargs)

    def get(self, name):
        """
        Gets the staircase for the name.

        Args:
            name: str
        Returns:
            Staircase
        """
        return self.staircase_dct[name] 
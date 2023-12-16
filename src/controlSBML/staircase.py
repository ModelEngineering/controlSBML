"""Create a staircase of values. Build multiple staircases."""
import controlSBML.constants as cn

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


C_NUM_POINT = "num_point"
C_NUM_STEP = "num_step"


class Staircase(object):

    def __init__(self, initial_value=cn.DEFAULT_INITIAL_VALUE,
                 num_step=cn.DEFAULT_NUM_STEP, final_value=cn.DEFAULT_FINAL_VALUE,
                 num_point=cn.DEFAULT_NUM_POINT, name=None):
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
        self.setNumPoint(num_point)  # self.num_point
        self.name = name
        self._updateState()
        #self.staircase_arr is dynamically

    def __repr__(self):
        dct = {"initial_value": self.initial_value, "num_step": self.num_step,
               "final_value": self.final_value, "num_point": self.num_point}
        return str(dct)

    def setNumPoint(self, num_point):
        self.num_point = num_point

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
            raise RuntimeError("Negative residual count")
        elif num_add_more > 0:
            final_values = np.repeat(self.final_value, num_add_more)
            staircase_arr = np.append(staircase_arr, final_values)
        #
        return staircase_arr
        

    @classmethod
    def makeRelativeStaircase(cls, center=5, fractional_deviation=0.1,
                              num_step=cn.DEFAULT_NUM_STEP, num_point=cn.DEFAULT_NUM_POINT):
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
        return cls(initial_value=center - max_deviation,
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
"""Sequence of values with associated times, not necessarily evenly spaced."""

from controlSBML import msgs

import numpy as np


class TemporalSeries(object):

    def __init__(self, times=None, values=None):
        """
        Args:
            times: list-float (positive)
            values: list-float
        """
        if times is None:
            times = []
        if values is None:
            values = []
        if len(times) != len(values):
            raise ValueError("len(times) != len(values)")
        self._times = list(times)
        self._values = list(values)

    def __len__(self):
        return len(self._times)

    @property 
    def num_nonzero_times(self):
        return len([t for t in self._times if not np.isclose(t, 0)])
    
    def add(self, time, value):
        """
        Adds a temporal value.

        Args:
            time: float
            value: float
        """
        self._times.append(time)
        self._values.append(value)

    def getValue(self, time, times=None, values=None):
        """
        Gets the value at the specified time. May require interpolation.

        Args:
            time: float
        Returns:
            float
        """
        if times is None:
            times = self._times
        if values is None:
            values = self._values
        #
        min_time = np.min(times)
        max_time = np.max(times)
        if (time < min_time) or (time > max_time):
            raise ValueError("time must be in [%f, %f]" % (min_time, max_time))
        lower_time = np.max([t for t in times if t <= time])
        upper_time = np.min([t for t in times if t >= time])
        if np.isclose(lower_time, upper_time):
            # Exact match
            chosen_time = lower_time
            idx = times.index(chosen_time)
            value = values[idx]
        else:
            # Interpolate
            lower_idx = times.index(lower_time)
            lower_val = values[lower_idx]
            upper_idx = times.index(upper_time)
            upper_val = values[upper_idx]
            time_span = upper_time - lower_time
            span_position = (time - lower_time) / time_span
            value_span = upper_val - lower_val
            value = lower_val + span_position*value_span
        return value

    def getTimesValues(self, count, precision=8):
        """
        Gets evenly spaced times and their associated values that span the range of times accumulated so far.

        Args:
            count: int (positive)
            precision: int (Number of decimal places)

        Returns:
            array-float (times)
            array-float (values)
        """
        min_time = np.min(self._times)
        max_time = np.max(self._times)
        range_time = max_time - min_time
        step_time = range_time / count
        # Do steps as integers to avoid floating point errors
        mult = 10**precision
        step = int(mult*np.round(step_time, precision))
        # Average values for repeated times
        new_times = []
        new_values = []
        for time in set(self._times):
            new_times.append(time)
            idxs = [i for i, t in enumerate(self._times) if np.isclose(t, time)]
            new_value = np.mean([self._values[i] for i in idxs])
            new_values.append(new_value)
        # Calculate values
        if not np.isclose(step, 0):
            times = np.linspace(min_time, max_time, count)
            values = [self.getValue(t, times=new_times, values=new_values) for t in times]
            # upper = min(int(mult*max_time+count*step), mult*max_time)
            # if upper > 0:
            #     times = np.arange(int(mult*min_time), upper, step)/mult
            #     values = [self.getValue(t) for t in times]
            # else:
            #     is_calculate = False
        else:
            msgs.warn("Not enough data to calculate evenly spaced times.")
            times = np.repeat(self._times[-1], len(self))
            values = np.repeat(self._values[-1], len(self))
        return times, values
"""Time Series Object"""

"""
A Time Series (TS) object is a DataFrame structured as:
  - Index is time in integer ms
  - Columns are variable names
  - name attribute is "Timeseries"

Note that if complex artithematic is failing using Timeseries, then
use the ".df" and ".ser" properties to convert back to a pandas objects if needed.
Alternatively, if a pandas object is returned from an operation,
then use Timeseries or TimeseriesSer to reconstruct the object.
"""

import controlSBML as ctl
from controlSBML import util
import controlSBML.constants as cn

import numpy as np
import pandas as pd # type: ignore

############# FUNCTIONS ###############
def findCommonIndices(index1, index2):
    """
    Finds the indices common to both.

    Parameters
    ----------
    index1: list/index
    index2: list/index

    Returns
    -------
    sorted list
    """
    indices = list(set(index1).intersection(index2))
    indices.sort()
    return(indices)

def align(ts1, ts2):
    """
    Returns objects with the same indices.

    Parameters
    ----------
    ts1: Timeseries/TimeseriesSer
    ts2: Timeseries/TimeseriesSer

    Returns
    -------
    Timeseries/TimeseriesSer, Timeseries/Timeseries/Ser
    """
    common_indices = findCommonIndices(ts1.index, ts2.index)
    new_ts1 = ts1.loc[common_indices, :]
    new_ts2 = ts2.loc[common_indices, :]
    return new_ts1, new_ts2

def convertTime(times):
    """
    Converts float seconds to int ms.

    Parameters
    ----------
    times: list-float

    Returns
    -------
    list-int
    """
    # Check if this is an index in the correct units
    if "pandas.core.indexes" in str(type(times)):
        if times.name == cn.TIMESERIES_INDEX_NAME:
            return list(times)
    # Must convert
    arr = np.array(times).astype(float)
    new_arr = arr*cn.MS_IN_SEC
    new_arr = new_arr.astype(int)
    return new_arr

def _calcTimes(ts):
    # Calculates evenly spaced times from a Timeseries index
    times_int = list(ts.index)
    dt = float((times_int[1] - times_int[0])*cn.SEC_IN_MS)
    start = times_int[0]*cn.SEC_IN_MS
    end = times_int[-1]*cn.SEC_IN_MS + dt
    times = np.arange(start, end, dt)
    times_diff = np.diff(times)
    # Ensure that times are evenly spaced
    if not np.allclose(times_diff, times_diff[0]):
        raise RuntimeError("Times are not evenly spaced.")
    return times


############# CLASSES ###############
class TimeseriesSer(pd.Series):

    def __init__(self, ser, times=None):
        if times is None:
            times = ser.index
        new_times = convertTime(times)
        super().__init__(ser, index=new_times)
        self.index.name = cn.TIMESERIES_INDEX_NAME

    @property
    def ser(self):
        return pd.Series(self)

    @property
    def times(self):
        return _calcTimes(self)

    def align(self, other):
        """
        Returns objects with the same indices.

        Parameters
        ----------
        other: Timeseries/TimeseriesSer

        Returns
        -------
        Timeseries/TimeseriesSer, Timeseries/Timeseries/Ser
        """
        common_indices = findCommonIndices(self.index, other.index)
        new_ts1 = self.loc[common_indices]
        if isinstance(other, Timeseries):
            new_ts2 = other.loc[common_indices, :]
        else:
            new_ts2 = other.loc[common_indices]
        return new_ts1, new_ts2


class Timeseries(pd.DataFrame):

    def __init__(self, mat, times=None, columns=None):

        """
        Parameters
        ----------
        mat: DataFrame, numpy.darray, NamedArray
        times: list-float (time in seconds)
        columns: list-str

        Notes:
            1. Assigning a value to times overrides an existing time index.
            2. A column labelled "time" overrides an index.
            3. For DataFrames, if index.name is "milliseconds", then times
               are not converted.
        """
        # The following blocks create df and times
        if isinstance(mat, Timeseries):
            df = mat
            times = df.index
        elif isinstance(mat, pd.Series):
            raise ValueError("Use TimeseriesSer, not TimeSeries.")
        elif isinstance(mat, pd.DataFrame):
            if columns is None:
                mat_columns = list(mat.columns)
            else:
                mat_columns = columns
            df = pd.DataFrame(mat, index=mat.index, columns=mat_columns)
            if times is None:
                if cn.TIME in mat.columns:
                    times = df[cn.TIME]
                    del df[cn.TIME]
                else:
                    times = df.index
        #
        elif "NamedArray" in str(type(mat)):
            if columns is None:
                mat_columns = list(mat.colnames)
            else:
                mat_columns = columns
            df = pd.DataFrame(mat, columns=mat_columns)
            if times is None:
                if cn.TIME in mat.colnames:
                    times = df[cn.TIME]
                    del df[cn.TIME]
                else:
                    raise ValueError("No time information found.")
        #
        elif isinstance(mat, np.ndarray):
            if columns is None:
                raise ValueError("No column information")
            else:
                mat_columns = columns
            df = pd.DataFrame(mat)
            num_column = len(df.columns)
            try:
                new_columns = list(df.columns)
                for pos, column in enumerate(mat_columns):
                    new_columns[pos] = column
                df.columns = new_columns
            except:
                raise RuntimeError("Unable to set column names.")
            if times is None:
                raise ValueError("No time information found.")
        else:
            raise ValueError("Unsupported data container.")
        #
        df.index = convertTime(times)
        # Fix the columns if needed
        new_columns = [str(c)[1:-1]
              if str(c)[0] == "[" else c for c in df.columns]
        df.columns = new_columns
        super().__init__(df)
        self.index.name = cn.TIMESERIES_INDEX_NAME

    @property
    def times(self):
        return _calcTimes(self)

    def __getitem__(self, key):
        """
        Return a Timeseries object.

        Parameters
        ----------
        key: column of Timeseries

        Returns
        -------
        Timeseries
        """
        item = super().__getitem__(key)
        if isinstance(item, pd.Series):
            ts = TimeseriesSer(item)
            ts.columns = [key]
            return ts
        else:
            return self.__class__(item, times=item.index)

    @property
    def df(self):
        return pd.DataFrame(self)

    def align(self, other):
        """
        Returns objects with the same indices.

        Parameters
        ----------
        other: Timeseries/TimeseriesSer

        Returns
        -------
        Timeseries/TimeseriesSer, Timeseries/Timeseries/Ser
        """
        common_indices = findCommonIndices(self.index, other.index)
        new_self = self.loc[common_indices, :]
        if isinstance(other, pd.DataFrame):
            new_other = other.loc[common_indices, :]
        else:  # Series
            new_other = other.loc[common_indices]
        return new_self, new_other

    @staticmethod
    def mat2TS(mat, column_names=None, row_names=None):
        """
        Converts a numpy ndarray or array-like to a Timeseries.

        Parameters
        ----------
        mat: np.Array, NamedArray, DataFrame
        column_names: list-str
        row_names: list-str
        """
        df = util.mat2DF(mat, column_names=column_names, row_names=row_names)
        return Timeseries(df)

    def copy(self, **kwargs):
        return Timeseries(self)

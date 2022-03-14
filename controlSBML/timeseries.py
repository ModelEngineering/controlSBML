"""Time Series Object"""

"""
A Time Series (TS) object is a DataFrame structured as:
  - Index is time in integer ms
  - Columns are variable names
  - name attribute is "Timeseries"

Note that if complex artithematic is failing using Timeseries, then 
use the to_pandas method and convert back to a Timeseries if needed.
"""

import controlSBML.constants as cn

import numpy as np
import pandas as pd


class TimeseriesSer(pd.Series):

    def to_pandas(self):
        return pd.Series(self)

    @property
    def times(self):
        times = np.array(self.index)
        return cn.SEC_IN_MS*times


class Timeseries(pd.DataFrame):

    def __init__(self, mat, times=None, columns=None, use_index=False):

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
        is_milliseconds = False
        if isinstance(mat, Timeseries):
            df = mat
            times = list(df.index)
            is_milliseconds = True
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
                    is_milliseconds = mat.index.name == cn.TIMESERIES_INDEX_NAME
                    times = list(df.index)
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
            df = pd.DataFrame(mat, columns=mat_columns)
            if times is None:
                raise ValueError("No time information found.")
        else:
            raise ValueError("Unsupported data container.")
        #
        if not is_milliseconds:
            df.index = self._convertTime(times)
        # Fix the columns if needed
        new_columns = [c[1:-1] if c[0] == "[" else c for c in df.columns]
        df.columns = new_columns
        super().__init__(df)
        self.index.name = cn.TIMESERIES_INDEX_NAME

    @property
    def times(self):
        times = np.array(self.index)
        return cn.SEC_IN_MS*times

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
               
    @staticmethod
    def _convertTime(times):
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

    def to_pandas(self):
        return pd.DataFrame(self)

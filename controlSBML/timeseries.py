"""Time Series Object"""

"""
A Time Series (TS) object is a DataFrame structured as:
  - Index is time in integer ms
  - Columns are variable names
  - name attribute is "TimeSeries"
"""

import controlSBML.constants as cn

import numpy as np
import pandas as pd


MS_IN_SEC =1000
INDEX_NAME = "miliseconds"


class TimeSeriesSer(pd.Series):
    pass


class TimeSeries(pd.DataFrame):

    def __init__(self, mat, times=None, columns=None, use_index=False):

        """
        Parameters
        ----------
        mat: DataFrame, numpy.darray, NamedArray
        times: list-float/list-int
            float: seconds
           int: ms
        columns: list-str
        use_index: bool
           Use the dataframe index as the milliseconds time
        """
        # The following blocks create df and times
        if isinstance(mat, TimeSeries):
            df = mat
            times = list(df.index)
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
                elif use_index:
                    times = list(df.index)
                else:
                    raise ValueError("No time information found.")
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
        df.index = self._convertTime(times)
        # Fix the columns if needed
        new_columns = [c[1:-1] if c[0] == "[" else c for c in df.columns]
        df.columns = new_columns
        super().__init__(df)
        self.index.name = INDEX_NAME

    def __getitem__(self, key):
        """
        Return a TimeSeries object.

        Parameters
        ----------
        key: column of TimeSeries
        
        Returns
        -------
        TimeSeries
        """
        item = super().__getitem__(key)
        if isinstance(item, pd.Series):
            ts = TimeSeriesSer(item)
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
        # Check if already correct units
        if isinstance(times[0], int):
            return times
        # Convert units
        arr = np.array(times).astype(float)
        new_arr = arr*MS_IN_SEC
        new_arr = new_arr.astype(int)
        return new_arr

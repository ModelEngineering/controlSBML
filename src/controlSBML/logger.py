"""Logs time series events"""

import numpy as np
import pandas as pd

TIME = "time"
TIMEUS = "timeus"  # Time in microseconds
SEC_TO_US = 1e6


class Logger(object):

    def __init__(self, logger_name, item_names):
        """
        Parameters
        ----------
        logger_name: str
            name of logger
        item_names: list-str
        """
        self.logger_name = logger_name
        self.item_names = item_names
        self.dct = None
        self.dct = {n: [] for n in self.item_names}
        self.dct[TIME] = []
        false = any(["fltr.fltr" in n for n in item_names])
        if false:
            import pdb; pdb.set_trace()
            pass

    def initialize(self):
        """
        Initializes the log.
        """
        for key in self.dct.keys():
            self.dct[key] = []

    def __repr__(self):
        return self.logger_name + " logger"

    def add(self, time, item_values):
        """
        Add an entry in the log.

        Parameters
        ----------
        time: float
        item_values: list-obj
        """
        self.dct[TIME].append(time)
        [self.dct[n].append(v) for n, v in zip(self.item_names, item_values)]

    def report(self, is_filter=True):
        """
        Generate a report of the log.

        Parameters
        ----------
        is_filter: bool
            Filter times to eliminate repetitions

        Returns
        -------
        pd.DataFrame
            index: time
            columns: item_names
        """
        # Check the logging lencths
        len_dct = {n: len(self.dct[n]) for n in self.dct.keys()}
        lst = list(len_dct.values())
        var = np.var(lst)
        if var > 0:
            stgs = ["%s=%d" % (k, v) for k, v in len_dct.items()]
            stg = "\n  ".join(stgs)
            text = "Unequal logging lencths. Details follow:\n  %s" % stg
            raise ValueError(text)
        df = pd.DataFrame(self.dct)
        # Eliminate duplicate times by using millisecond granularity
        if is_filter:
            df[TIMEUS] = df[TIME].apply(lambda v: int(SEC_TO_US*v))
            dfg = df.groupby([TIMEUS])
            sers = []
            for idx, indices in dfg.groups.items():
                new_indices = list(indices)
                index = new_indices[-1]  # pick last one
                ser = df.loc[index, :]
                if isinstance(ser, pd.DataFrame):
                    ser = ser.reset_index()
                    ser = ser.loc[0, :]
                sers.append(ser)
            new_df = pd.DataFrame(sers)
            new_df[TIME] = new_df[TIMEUS].apply(lambda v: v/SEC_TO_US)
            del new_df[TIMEUS]
        else:
            new_df = df
        # Format the dataframe
        new_df = new_df.set_index(TIME)
        return new_df

    def merge(self, new_name, other_logs):
        """
        Combines multiple logs.

        Parameters
        ----------
        new_name: str
        other_logs: list-Logger
        
        Returns
        -------
        Logger
        """
        def combineName(logger_name, item_name):
            return "%s.%s" % (logger_name, item_name)
        def mergeNames(loggers):
            names = []
            for logger in loggers:
                names.extend([combineName(logger.logger_name, n)
                      for n in logger.item_names])
            return names
        #
        loggers = list(other_logs)
        loggers.append(self)
        new_item_names = mergeNames(loggers)
        new_logger = Logger(new_name, new_item_names)
        new_logger.dct[TIME] = self.dct[TIME]
        for logger in loggers:
            for item_name in logger.item_names:
                new_name = combineName(logger.logger_name, item_name)
                new_logger.dct[new_name] = logger.dct[item_name]
        return new_logger
            
        
        

   

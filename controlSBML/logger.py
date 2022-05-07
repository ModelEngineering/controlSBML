"""Logs time series events"""

import numpy as np
import pandas as pd

TIME = "time"


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
        self.dct = {n: [] for n in item_names}
        self.dct[TIME] = []

    def __repr__(self):
        return self.logger_name

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

    def report(self):
        """
        Generate a report of the log.

        Returns
        -------
        pd.DataFrame
            index: time
            columns: item_names
        """
        df = pd.DataFrame(self.dct)
        df = df.set_index(TIME)
        return df

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
            
        
        

   

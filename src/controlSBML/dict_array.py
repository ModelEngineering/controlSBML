"""Dictionary whose values are lists of the same length."""

"""
 Usage:
    dl = DictList(key1, key2, ...)   # Insert None to maintain length if a value is missing
    dl.append(key1=value1, key2=value2, ...)
    dl.addkey(key)  # Adds the key, inserting leading None values to maintain length
    len(dl)        # Returns the length of the lists
"""

import pandas as pd
import numpy as np
from typing import List

class DictArray(dict):

    def __init__(self, *args:List[str]):
        """
        Args:
            args (List[str]): List of keys
        """
        super().__init__()
        for key in self.keys():
            self[key] = []

    def __len__(self)->int:
        keys = list(self.keys())
        if len(keys) == 0:
            return 0
        first_key = keys[0]
        return len(self[first_key])
    
    def __eq__(self, other)->bool:
        IS_DEBUG = False
        if not isinstance(other, DictArray):
            if IS_DEBUG:
                print("1")
            return False
        if len(self) != len(other):
            if IS_DEBUG:
                print("2")
            return False
        for key in self.keys():
            if not key in other.keys():
                if IS_DEBUG:
                    print("3")
                return False
            for key in self.keys():
                values1 = self[key]
                values2 = other[key]
                for value1, value2 in zip(values1, values2):
                    if (value1 is None) and (value2 is None):
                        continue
                    if (isinstance(value1, float) and isinstance(value2, float)):
                        if np.isclose(value1, value2):
                            continue
                        else:
                            if IS_DEBUG:
                                print("4")
                            return False
                    if self._isNa(value1) and self._isNa(value2):
                        continue
                    if value1 != value2:
                        if IS_DEBUG:
                            print("5")
                        return False
        return True

    @staticmethod 
    def _isNa(value)->bool:
        if isinstance(value, float):
            if np.isnan(value):
                return True
        return value is None
    
    def copy(self)->'DictArray':
        dl = DictArray()
        for key, value in self.items():
            dl[key] = list(value)
        return dl
        
    
    def append(self, **kwargs)->None:
        """
        Appends values to list. If a new key is encountered include it (is_addkey=True)

        Args:
            kwargs: dict of key, value pairs
        """
        fill_list = list(np.repeat(None, len(self)))
        for key, value in kwargs.items():
            if not key in self.keys():
                self[key] = list(fill_list)
            self[key].append(value)
        for key in self.keys():
            if not key in kwargs.keys():
                self[key].append(None)

    @classmethod
    def makeFromDicts(cls, dcts:List[dict])->'DictArray':
        """
        Makes a DictList from a list of dicts.

        Args:
            dcts (List[dict]): List of dicts

        Returns:
            DictList: DictList
        """
        dl = DictArray()
        for dct in dcts:
            dl.append(**dct)
        return dl
    
    def makeDicts(self)->List[dict]:
        """
        Creates a list of dictionaries

        Returns:
            List[dict]:
        """
        dcts = []
        for i in range(len(self)):
            dct = {}
            for key in self.keys():
                dct[key] = self[key][i]
            dcts.append(dct)
        return dcts

    def getDataframe(self)->pd.DataFrame:
        return pd.DataFrame(self)

    @classmethod
    def merge(cls, dict_lists)->'DictArray':
        """
        Merges DictLists.

        Args:
            dict_lists (List[DictList]): List of DictLists

        Returns:
            DictList: merged DictList
        """
        merged_dl = DictArray()
        # Get the list of keys
        keys = []
        for dl in dict_lists:
            keys.extend(dl.keys())
        keys = list(set(keys))
        # Process the lists
        for dl in dict_lists:
            merged_fills = list(np.repeat(None, len(merged_dl)))
            dl_fills = list(np.repeat(None, len(dl)))
            for key in keys:
                if not key in merged_dl.keys():
                    merged_dl[key] = list(merged_fills)
                if key in dl.keys():
                    merged_dl[key].extend(dl[key])
                else:
                    merged_dl[key].extend(dl_fills)
        return merged_dl
    
    @classmethod
    def makeFromDataframe(cls, df:pd.DataFrame)->'DictArray':
        dct = df.to_dict(orient='list')
        return cls.merge([dct])
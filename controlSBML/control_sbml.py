"""Main class for SBMLControl"""

import control as c
import numpy as np
import pandas as pd
import tellurium as te

class SBMLControl(object):

    def __init__(file_path=None, url=None, model_str=None, is_sbml=True):
        """
        Parameters
        ----------
        file_path: str
        url: str
        model_str: str
        is_sbml: bool
            True: SBML
            False: Antimony
        """
        self.road_runner = None
        if is_sbml:
            pass
        else:
            if sring is not None:
                self.road_runner = te.loada(model_str)
        if self.road_runner is None:
            raise ValueError("Invalid model specification")

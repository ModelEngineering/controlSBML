"""Evaluator used with parallel search."""

from controlSBML.parallel_search import Evaluator
import controlSBML.constants as cn

import numpy as np
import tellurium as te


class TelluriumEvaluator(Evaluator):

    def __init__(self, model_url:str):
        self.model_url = model_url

    def initialize(self):
        self.rr = te.loadSBMLModel(self.model_url)

    def evaluate(self, **kwargs)->dict:
        self.rr.reset()
        for key, value in kwargs.items():
            #self.rr[key] = value
            pass
        data = self.rr.simulate(0, 10, 100)
        score = np.mean(data["[s6]"])
        return_dct = dict(kwargs)
        return_dct[cn.SCORE] = score
        return return_dct


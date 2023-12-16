from controlSBML.parallel_search import ParallelSearch, Evaluator
from controlSBML.dict_array import DictArray
from parallel_search_evaluator import TelluriumEvaluator
import controlSBML.constants as cn

import numpy as np
import time
import unittest


IGNORE_TEST = False
IS_PLOT = False
CANDIDATES = [{"kp": 1}, {"ki": 2, "kp": 3}, {"kf": 4, "kp": 5}]


class MyEvaluator(Evaluator):

    def __init__(self, wait_time:float=0):
        self.wait_time = wait_time

    def initialize(self):
        self.wait_time = self.wait_time

    def evaluate(self, **kwargs)->dict:
        if self.wait_time > 0:
            time.sleep(self.wait_time)
        score = np.prod([v for v in kwargs.values()])
        return_dct = dict(kwargs)
        return_dct["score"] = score
        [np.exp(0.001*i) for i in range(1000)]    # Use CPU
        return return_dct


#############################
# Tests
#############################
class TestParallelSearch(unittest.TestCase):

    def setUp(self):
        self.evaluator = MyEvaluator()
        self.search = ParallelSearch(self.evaluator, CANDIDATES, num_process=2)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(self.search is not None)

    def testEvaluateCandidates(self):
        if IGNORE_TEST:
            return
        return_dct = {}
        self.search._evaluateCandidates(1, self.evaluator, CANDIDATES, 1, return_dct)
        search_result = return_dct[1]
        self.assertTrue(isinstance(search_result, DictArray))
        self.assertEqual(len(CANDIDATES), len(search_result))
        score = np.prod([v for v in CANDIDATES[0].values() if v is not None])
        self.assertEqual(score, search_result[cn.SCORE][0])

    def doSearch(self, num_process, evaluator=None, candidates=None, expected_score=1):
        if candidates is None:
            candidates = CANDIDATES
        if evaluator is None:
            evaluator = MyEvaluator()
        search = ParallelSearch(evaluator, candidates, num_process=num_process)
        search.search()
        best_candidate = search.getBestCandidate()
        self.assertEqual(len(search._search_results_dict_array), len(candidates))
        self.assertTrue(best_candidate is not None)
        if expected_score is not None:
            self.assertEqual(best_candidate[cn.SCORE], expected_score)
        trues = [c in best_candidate.keys() for c in ["kp", "ki", "kf", "score"]]
        self.assertTrue(all(trues))

    def testSearchOneProcess(self):
        if IGNORE_TEST:
            return
        self.doSearch(1)

    def testSearchTwoProcess(self):
        if IGNORE_TEST:
            return
        self.doSearch(2)

    def testSearchManyProcess(self):
        if IGNORE_TEST:
            return
        candidates = []
        [candidates.extend(CANDIDATES) for _ in range(83)]
        evaluator = MyEvaluator(wait_time=0.1)
        self.doSearch(5, candidates=candidates, evaluator=evaluator)
    
    def testSearchTellurium(self):
        if IGNORE_TEST:
            return
        candidates = []
        [candidates.extend(CANDIDATES) for _ in range(100)]
        evaluator = TelluriumEvaluator(cn.WOLF_URL)
        self.doSearch(5, candidates=candidates, evaluator=evaluator, expected_score=None)

if __name__ == '__main__':
    unittest.main()
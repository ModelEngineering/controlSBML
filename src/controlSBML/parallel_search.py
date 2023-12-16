"""Finds the Candidate with the lowest score."""

"""
Parallel search implements a parallel search of a collection of candidate.
A candidate is a dict whose keys are parameter names and whose values are a value assigned to the parameter.
An Evaluator calculates a score for a candidate and returns a dictionary with the score and other information about the
candidate that produced the score. The evaluator should be constructed to be simple and light weight
(e.g., simple types such as str, float, array). It has an initialize method that creates heavy weight objects
(e.g., a Roadrunner instance).
A ParalellSearch is a singleton that is constructed with an Evaluator and a list of candidates. The search method
calculates the best search result.

Usage:
    from parallel_search import ParallelSearch
    evaluator = MyEvaluator(arguments)   # Construct the evaluator
    search = ParallelSearch(evaluator, candidates)
    search.search()   # Runs in the evaluators in parallel parallel
    best_result = search.getBestCandidate()  # Returns the candidate with the lowest score
    df = search.getSearchResuts()  # Returns a dataframe with the scores and other results from the evaluator
"""
import controlSBML.constants as cn
from controlSBML.dict_array import DictArray

import multiprocessing as mp
import numpy as np
import pandas as pd
import random
from typing import List
from tqdm import tqdm


class Evaluator(object):
    # The constructed object should be light weight and pickleable.
    
    def initialize(self)->None:
        """
        Initializes the evaluator. Called before evaluating candidates.
        """
        raise NotImplementedError("Subclass must implement")

    def evaluate(self, **kwargs)->dict:
        """
        Evaluates a candidate.

        Returns:
            dict:
                keys: parameters, "score", other
                values: float, str, int, bool
        """
        raise NotImplementedError("Subclass must implement")


class ParallelSearch(object):
    # Finds the candidate with the lowest score provided by an evaluator. Operates in parallel. Reports progress.
    # Finds the best candidate.

    def __init__(self, evaluator:Evaluator, candidates:List[dict], num_process:int=-1, is_report:bool=False):
        """
        Args:
            evaluator (Evaluator): evaluates a candidate
            candidates (list-dict): each dict is a candidate
            num_process (int): (default=-1, use all available processors)
            is_report (bool): True to print reports
        """
        self.evaluator = evaluator
        self.candidates = candidates
        if num_process == -1:
            num_process = mp.cpu_count()
        self.num_process = num_process
        self.is_report = is_report
        # Outputs
        self._search_results_dict_array = None
        self._best_candidate = None

    def getBestCandidate(self):
        if self._best_candidate is None:
            dict_list = self._search_results_dict_array.makeDicts()
            sorted_candidates = sorted(dict_list, key=lambda dct: dct[cn.SCORE])
            self._best_candidate = sorted_candidates[0]
        return self._best_candidate
    
    def getSearchResults(self):
        if self._search_results_dict_array is None:
            raise ValueError("Must run search() first")
        return self._search_results_dict_array.getDataframe()

    def search(self):
        num_process = min(len(self.candidates), self.num_process)
        if num_process == 1:
            # Run in a single process
            procnum = 0
            return_dct = {}
            self._evaluateCandidates(procnum, self.evaluator, self.candidates, num_process, return_dct)
            self._search_results_dict_array = return_dct[procnum]
        else:
            # Initialize for parallel calculations
            jobs = []
            num_process = min(num_process, len(self.candidates))
            random.shuffle(self.candidates)   # Process in random order
            num_candidate_per_process = int(len(self.candidates)//num_process)
            extra_candidates = len(self.candidates) % num_process
            manager = mp.Manager()
            return_dct = manager.dict()
            candidates = list(self.candidates)
            # Start the processes
            for procnum in range(num_process):
                if self.is_report:
                    print("**Starting process %d" % procnum)
                pos = min(num_candidate_per_process, len(candidates))
                if procnum < extra_candidates:
                    # Add one extra candidate to the first few processes
                    pos += 1
                these_candidates = candidates[:pos]
                candidates = candidates[pos:]
                if procnum == num_process - 1:
                    # Add the remaining candidates to the last process
                    these_candidates.extend(candidates)
                p = mp.Process(target=self._evaluateCandidates, args=(procnum, self.evaluator,
                      these_candidates, num_process, return_dct))
                jobs.append(p)
                p.start()
            # Wait for the processes to finish
            for proc in jobs:
                proc.join()
                proc.terminate()
            if self.is_report:
                print("**All processes finished")
            # Updates the candidates
            self._search_results_dict_array = DictArray.merge(return_dct.values())

    @classmethod
    def _evaluateCandidates(cls, procnum:int, evaluator:Evaluator,
                            candidates:List[dict], num_process:int, return_dct:dict)->None:
        """
        Evaluates a candidate in another process. Updates the score for each candidate. Runs in a separate process.

        Args:
            procnum (int): Process number
            num_process (int):  number of processes spawned
            candidates (List[dict])
        """
        evaluator.initialize()
        results = DictArray()
        if procnum == 0:
            # Show progress via process 0
            num_iteration = len(candidates)*num_process # Count for all proceses
            pos = 0  # Position in candidates
            for idx in tqdm(range(num_iteration)):
                residual_cnt = idx % num_process
                if residual_cnt == 0:
                    candidate_dct = candidates[pos]
                    result_dct = evaluator.evaluate(**candidate_dct)
                    results.append(**result_dct)
                    pos += 1
        else:
            for candidate_dct in candidates:
                result_dct = evaluator.evaluate(**candidate_dct)
                results.append(**result_dct)
        return_dct[procnum] = results
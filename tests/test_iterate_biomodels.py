from controlSBML.iterate_biomodels import iterateBiomodels

import os
import unittest
import tellurium as te


IGNORE_TEST = False
IS_PLOT = False


#############################
# Tests
#############################
class TestIterateBiomodels(unittest.TestCase):

    def setUp(self):
        pass

    def testIterateBiomodels(self):
        if IGNORE_TEST:
            return
        iterator = iterateBiomodels(is_report=True)
        ffiles = [f for f, _ in iterator]
        self.assertGreater(len(ffiles), 1000)

    def testIterateBiomodels2(self):
        if IGNORE_TEST:
            return
        iterator = iterateBiomodels(start=1, end=3, is_report=IS_PLOT)
        for ffile, content in iterator:
            rr = te.loads(content)
            self.assertTrue("RoadRunner" in str(type(rr)))

    def testCheckerFunction(self):
        if IGNORE_TEST:
            return
        def checker(filename, contents):
            if "BIOMD0000000001.xml" in filename:
                return "Model #1 is skipped"
            else:
                return ""
        iterator = iterateBiomodels(start=1, end=2, is_report=IS_PLOT, checkerFunctions=[checker])
        ffiles = [f for f, _ in iterator]
        self.assertEqual(len(ffiles), 1)


if __name__ == '__main__':
  unittest.main()

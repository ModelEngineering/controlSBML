from controlSBML.iterate_biomodels import iterateBiomodels

import os
import unittest
import tellurium as te


IGNORE_TEST = True
IS_PLOT = True


#############################
# Tests
#############################
class TestIterateBiomodels(unittest.TestCase):

    def setUp(self):
        pass

    def testIterateBiomodels(self):
        iterator = iterateBiomodels(is_report=True)
        ffiles = [f for f, _ in iterator]
        self.assertGreater(len(ffiles), 1000)

    def testIterateBiomodels2(self):
        iterator = iterateBiomodels(start=1, end=3, is_report=IS_PLOT)
        for ffile, content in iterator:
            rr = te.loads(content)
            self.assertTrue("RoadRunner" in str(type(rr)))



if __name__ == '__main__':
  unittest.main()

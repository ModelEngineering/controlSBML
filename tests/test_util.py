from controlSBML import util

import tellurium as te
import unittest
import numpy as np


IGNORE_TEST = False
IS_PLOT = False


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

    def setUp(self):
        pass

    def testCacluateMatrixDistance(self):
        if IGNORE_TEST:
            return
        SHAPE = (4, 4)
        SIZE = SHAPE[0]*SHAPE[1]
        mat1 = np.array(range(SIZE))
        mat1 = np.reshape(mat1, SHAPE)
        dist = util.calculateMatrixDistance(mat1, mat1)
        self.assertTrue(np.isclose(dist, 0))
        #
        mat2 = 4*mat1
        dist1 = util.calculateMatrixDistance(mat1, mat2)
        self.assertGreater(dist1, 0)
        dist2 = util.calculateMatrixDistance(mat2, mat1)
        self.assertTrue(np.isclose(dist1, dist2))

    def testGetModel(self):
        if IGNORE_TEST:
          return
        model_str = util.getModel()
        rr = te.loada(model_str)
        self.assertTrue("roadrunner" in str(type(rr)))
 


if __name__ == '__main__':
  unittest.main()

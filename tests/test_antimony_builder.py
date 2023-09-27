import controlSBML.constants as cn
from controlSBML import antimony_builder as ab

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import unittest
import tellurium as te


IGNORE_TEST = False
IS_PLOT = False
LINEAR_MDL = """
// Illustrate Antimony File
model *main_model()
S1 -> S2; k1*$S1
J1: S2 -> S3; k2*S2
J2: S3 -> S2; k3*S3
J3: S2 -> ; k4*S2

k1 = 1
k2 = 2
k3 = 3
k4 = 4
S1 = 10
S2 = 0
S3 = 0
end
"""
SYMBOL_DCT = {"S1": cn.TYPE_FLOATING_SPECIES, "S2": cn.TYPE_FLOATING_SPECIES, "S3": cn.TYPE_FLOATING_SPECIES}


#############################
# Tests
#############################
class TestAntimonyBuilder(unittest.TestCase):

    def setUp(self):
       self.builder = ab.AntimonyBuilder(LINEAR_MDL, SYMBOL_DCT)

    def check(self):
       rr = te.loada(str(self.builder))
       data = rr.simulate()
       self.assertTrue(len(data) > 0)
       if IS_PLOT:
          rr.plot()

    def testConstructor(self):
       if IGNORE_TEST:
           return
       self.assertTrue(isinstance(self.builder.antimony, str))

    def testMakeBoundarySpecies(self):
       if IGNORE_TEST:
           return
       self.builder.makeBoundarySpecies("S1")
       self.assertTrue("const" in self.builder.antimony_strs[-1])
       self.check()

    def testMakeBoundaryReaction(self):
       if IGNORE_TEST:
           return
       self.builder.makeBoundaryReaction("S1")
       self.assertTrue("->" in self.builder.antimony_strs[-2])
       self.check()

    def testMakeComment(self):
       if IGNORE_TEST:
           return
       self.builder.makeComment("comment")
       self.assertTrue("comment" in self.builder.antimony_strs[-1])

    def testMakeStaircase(self):
       if IGNORE_TEST:
           return
       self.builder.makeBoundarySpecies("S1")
       value_arr = self.builder.makeStaircase("S1", initial_value=2)
       self.assertTrue("at " in self.builder.antimony_strs[-1])
       self.assertEqual(len(value_arr), len(cn.TIMES))
       self.check()
    
    def testMakeSISClosedLoop(self):
       if IGNORE_TEST:
           return
       self.builder.makeBoundarySpecies("S1")
       self.builder.makeSISOClosedLoop("S1", "S3", kp=1, reference=4)
       self.check()
       

if __name__ == '__main__':
  unittest.main()
from controlSBML import constants as cn
from controlSBML import control_sbml
from controlSBML.control_sbml import ControlSBML
from tests import helpers

import copy
import numpy as np
import os
import unittest
import tellurium as te


IGNORE_TEST = False
IS_PLOT = False
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
ANTIMONY_FILE = os.path.join(TEST_DIR, "Model_antimony.ant")


#############################
# Tests
#############################
class TestControlSBML(unittest.TestCase):

    def setUp(self):
      # Cannot modify self.control
      self.ctl = ControlSBML(ANTIMONY_FILE)
  
    def testConstructor(self):
        if IGNORE_TEST:
            return
    

if __name__ == '__main__':
  unittest.main()

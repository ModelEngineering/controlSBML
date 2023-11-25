from controlSBML.option_set import OptionSet

import os
import unittest
import tellurium as te


IGNORE_TEST = False
IS_PLOT = False
OPTION_DCT = {
    "opt1": 1,
    "opt2": 2,
}


#############################
# Tests
#############################
class TestOptionSet(unittest.TestCase):

    def setUp(self):
        self.option_set = OptionSet(**OPTION_DCT)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertEqual(self.option_set.opt1, 1)
        self.assertEqual(self.option_set.opt2, 2)
        #
        diff = set(self.option_set.keys()).difference(set(OPTION_DCT.keys()))
        self.assertEqual(len(diff), 0)

    def testSetOptionSetDefaults(self):
        if IGNORE_TEST:
            return
        self.option_set.setOptionSet(opt1=3)
        self.assertEqual(self.option_set.opt1, 3)

    def testCopyOptionSet(self):
        if IGNORE_TEST:
            return
        new_option_set = self.option_set._copyOptionSet()
        self.assertEqual(new_option_set, self.option_set)

    def testUpdateOptionSet(self):
        if IGNORE_TEST:
            return
        new_option_set = self.option_set.getOptionSet(opt1=3)
        self.assertEqual(new_option_set.opt1, 3)
        self.assertEqual(self.option_set.opt1, 1)



if __name__ == '__main__':
  unittest.main()

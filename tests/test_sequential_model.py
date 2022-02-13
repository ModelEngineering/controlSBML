from controlSBML.sequential_model import SequentialModel

import tellurium as te
import unittest


IGNORE_TEST = True
IS_PLOT = True
NUM_REACTION = 4


#############################
# Tests
#############################
class TestSequentialModel(unittest.TestCase):

    def setUp(self):
        self.smodel = SequentialModel(NUM_REACTION)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertEqual(len(self.smodel.species_values), NUM_REACTION + 1)

    def testMakeKineticsFormula(self):
        if IGNORE_TEST:
          return
        num = 20
        formula_str = self.smodel._makeKineticsFormula(num)
        self.assertEqual(formula_str.count(str(num)), 2)

    def testGenerateDefault(self):
        if IGNORE_TEST:
          return
        model_str = self.smodel.generate()
        rr = te.loada(model_str)
        self.assertTrue("roadrunner" in str(type(rr)))

    def testGenerateDifferentKineticsAndInitialValues(self):
        # TESTING
        kinetics_values = [.1, .2, .3]
        species_values = [0, 1, 2, 3]
        smodel = SequentialModel(len(kinetics_values),
               kinetics_values=kinetics_values,
               species_values=species_values)
        model_str = smodel.generate()
        self.assertTrue("k0 = 0.1" in model_str)
        rr = te.loada(model_str)
        self.assertTrue("roadrunner" in str(type(rr)))


if __name__ == '__main__':
  unittest.main()

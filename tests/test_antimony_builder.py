import controlSBML.constants as cn
from controlSBML import antimony_builder as ab

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import unittest
import tellurium as te


IGNORE_TEST = False
IS_PLOT = False
MODEL_NAME = "main_model"
INCOMPLETE_LINEAR_MDL = """
S1 -> S2; k1*S1
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
"""
LINEAR_MDL = "model *%s()\n" % MODEL_NAME + INCOMPLETE_LINEAR_MDL + "end"
SYMBOL_DCT = {"S1": cn.TYPE_FLOATING_SPECIES, "S2": cn.TYPE_FLOATING_SPECIES, "S3": cn.TYPE_FLOATING_SPECIES}
rr = te.loadSBMLModel(cn.MTOR_URL)
MTOR_MDL = rr.getAntimony()
MTOR_NAME = "Varusai2018___Dynamic_modelling_of_the_mTOR_signalling_network_reveals_complex_emergent_behaviours_conferred_by_DEPTOR"


#############################
# Tests
#############################
class TestAntimonyBuilder(unittest.TestCase):

    def setUp(self):
        if IGNORE_TEST:
            return
        self.init()

    def init(self):
        if "builder" in dir(self):
            return
        self.builder = ab.AntimonyBuilder(LINEAR_MDL, symbol_dct=SYMBOL_DCT)

    def check(self, builder=None):
        if builder is None:
            builder = self.builder
        rr = te.loada(str(builder))
        data = rr.simulate(0,20, 2000, selections=["time", "S1", "S2", "S3"])
        self.assertTrue(len(data) > 0)
        if IS_PLOT:
            rr.plot()
        return data
    
    def testProperties(self):
        if IGNORE_TEST:
            return
        builder = ab.AntimonyBuilder(MTOR_MDL)
        self.assertGreater(len(builder.floating_species_names), 0)
        self.assertEqual(len(builder.boundary_species_names), 0)
        self.assertGreater(len(builder.reaction_names), 0)
        self.assertGreater(len(builder.parameter_names), 0)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.init()
        self.assertTrue(isinstance(self.builder.antimony, str))

    def testMakeBoundarySpecies(self):
        if IGNORE_TEST:
            return
        self.init()
        self.builder.makeBoundarySpecies("S1")
        self.assertTrue("const" in self.getStatement())
        self.check()

    def testMakeBoundaryReaction(self):
        if IGNORE_TEST:
            return
        self.init()
        self.builder.makeBoundaryReaction("S1")
        self.assertTrue("->" in self.getStatement(pos=2))
        self.check()

    def getStatement(self, pos=1, builder=None):
        if builder is None:
            builder = self.builder
        return builder.antimony_strs[builder.insert_pos-pos]

    def testMakeComment(self):
        if IGNORE_TEST:
            return
        self.init()
        self.builder.makeComment("comment")
        self.assertTrue("comment" in self.getStatement())

    def testMakeAdditionStatement(self):
        if IGNORE_TEST:
            return
        self.init()
        self.builder.makeAdditionStatement("S1", "S2", "S3")
        result = re.search("S1.*:=.*S2.*\+.*S3", self.getStatement())
        self.assertTrue(result)
        self.builder.makeAdditionStatement("S2", "S3", is_assignment=False)
        result = re.search("S2.* =.*S3", self.getStatement())
        self.assertTrue(result)

    def testMakeSinusoidSignal(self):
        if IGNORE_TEST:
            return
        ot_name = self.builder.makeSinusoidSignal(1, 2, suffix="_S1_S2")
        result = re.search("%s.*=.*1.*sin.*2" % ot_name, self.getStatement())
        self.assertTrue(result)
        self.builder.makeBoundarySpecies("S1")
        self.builder.makeAdditionStatement("S1", ot_name)
        self.check()

    def testMakeFilterElement(self):
        if IGNORE_TEST:
            return
        self.init()
        filter_in, filter_ot = self.builder.makeFilterElement(1.0, suffix="_S1_S3")
        sin_ot = self.builder.makeSinusoidSignal(1, 2, suffix="_S1_S2")
        self.builder.makeAdditionStatement(filter_in, sin_ot)
        self.check()

    def testMakeControlErrorSignal(self):
        if IGNORE_TEST:
            return
        self.init()
        def test(sign):
            signal_ot = self.builder.makeControlErrorSignal(-7, "S3", sign, suffix="_S1_S3")
            if sign == -1:
                result = re.search("%s.*:=.*7.*-.*S3" % signal_ot, self.getStatement())
            else:
                result = re.search("%s.*:=.*7.*\+.*S3" % signal_ot, self.getStatement())
            self.assertTrue(result)
            self.check()
        #
        test(-1)
        test(1)
    
    def testMakePIController(self):
        if IGNORE_TEST:
            return
        self.init()
        name_in, name_ot = self.builder.makePIControllerElement(kp=7, suffix="_S1_S3")
        self.builder.makeBoundarySpecies("S1")
        self.builder.makeAdditionStatement("S1", name_ot)
        self.builder.makeAdditionStatement(name_in, 3, "-"+"S3")
        self.check()

    def testMakePIControllerInputIsAParameter(self):
        if IGNORE_TEST:
            return
        builder = ab.AntimonyBuilder(LINEAR_MDL, symbol_dct=SYMBOL_DCT)
        name_in, name_ot = builder.makePIControllerElement(kp=7, suffix="_S1_S3")
        builder.makeAdditionStatement("k0", name_ot)
        builder.makeAdditionStatement(name_in, 3, "-"+"S3")
        self.check(builder=builder)

    def testMakeSISOClosedLoopSystem(self):
        if IGNORE_TEST:
            return
        self.init()
        self.builder.makeBoundarySpecies("S1")
        self.builder.makeSISOClosedLoopSystem("S1", "S3", kp=10, setpoint=5,
                           noise_amplitude=1, noise_frequency=20, disturbance_ampliude=2, disturbance_frequency=20)
        self.check()
        self.builder.makeSISOClosedLoopSystem("S1", "S3", kp=10, setpoint=5,
                           noise_amplitude=1, noise_frequency=20, disturbance_ampliude=2, disturbance_frequency=20,
                           initial_output_value=33)
        self.check()

    def testMakeStaircase(self):
        if IGNORE_TEST:
            return
        def test(is_fixed_input_species):
            builder = ab.AntimonyBuilder(LINEAR_MDL, symbol_dct=SYMBOL_DCT)
            if is_fixed_input_species:
                builder.makeBoundarySpecies("S1")
            else:
                builder.makeBoundaryReaction("S1")
            value_arr = builder.makeStaircase("S1", initial_value=2)
            self.assertTrue("at " in self.getStatement(builder=builder))
            self.assertEqual(len(value_arr), len(cn.TIMES))
            self.check(builder=builder)
            return builder
        #
        _ = test(True)
        _ = test(False)
    
    def testMakeSISClosedLoop(self):
        if IGNORE_TEST:
            return
        self.init()
        self.builder.makeBoundarySpecies("S1")
        self.builder.makeSISOClosedLoopSystem("S1", "S3", kp=10000, ki=1, setpoint=4)
        data = self.check()
        self.assertTrue(np.isclose(data["S3"][-1], 4, atol=0.01))

    def testCopyAndEqual(self):
        if IGNORE_TEST:
            return
        self.init()
        builder = self.builder.copy()
        self.assertTrue(builder == self.builder)
        #
        builder.makeBoundarySpecies("S1")
        self.assertFalse(builder == self.builder)
       

if __name__ == '__main__':
  unittest.main()
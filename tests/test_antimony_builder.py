import controlSBML.constants as cn # type: ignore
from controlSBML import antimony_builder as ab

import numpy as np
import re # type: ignore
import unittest
import tellurium as te # type: ignore


IGNORE_TEST = True
IS_PLOT = True
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
        self.builder = ab.AntimonyBuilder(LINEAR_MDL, symbol_dct=SYMBOL_DCT)

    def check(self, builder=None, times=None):
        if builder is None:
            builder = self.builder
        if times is None:
            times = np.linspace(0, 20, 2000)
        rr = te.loada(str(builder))
        selections = ["time", "S1", "S2", "S3"]
        if "setpoint_S1_S3" in rr.keys():
            selections.append("setpoint_S1_S3")
        if "noise_S1_S3_ot" in rr.keys():
            selections.append("noise_S1_S3_ot")
        data = rr.simulate(times[0], times[-1], len(times), selections=selections)
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
        #
        self.builder.makeAdditionStatement("S2", "S3", is_assignment=False)
        result = re.search("S2.* =.*S3", self.getStatement())
        self.assertTrue(result)
        #
        self.builder.makeAdditionStatement("S2", "S3", is_assignment=False, comment="comment")
        result = re.search("S2.* =.*S3.*# comment", self.getStatement())
        self.assertTrue(result)

    def testMakeNoiseElement(self):
        #if IGNORE_TEST:
        #    return
        self.init()
        noise_spec = cn.NoiseSpec(sine_amp=0.01, sine_freq=0.2, random_mag=0.3,
                                    random_std=0.004, offset=0.5, slope=0.0006)
        ot_name = self.builder.makeNoiseElement(noise_spec, suffix="_S1_S2")
        result = re.search("%s.*=.*1.*sin.*2.*3*lognormal.*4.*5.*6.*time" % ot_name, self.getStatement())
        self.assertTrue(result)
        self.builder.makeBoundarySpecies("S1")
        self.builder.makeAdditionStatement("S1", ot_name)
        import pdb; pdb.set_trace()
        self.check()

    def testMakeFilterElement(self):
        if IGNORE_TEST:
            return
        self.init()
        filter_in, filter_ot, calculation = self.builder.makeFilterElement(1.0, suffix="_S1_S3")
        noise_spec = cn.NoiseSpec(sine_amp=1, sine_freq=2)
        sin_ot = self.builder.makeNoiseElement(noise_spec, suffix="_S1_S2")
        self.builder.makeAdditionStatement(filter_in, sin_ot)
        self.check()

    def testMakeFilterElementNofilter(self):
        if IGNORE_TEST:
            return
        self.init()
        filter_in, filter_ot, calculation = self.builder.makeFilterElement(0, suffix="_S1_S3")
        noise_spec = cn.NoiseSpec(sine_amp=1, sine_freq=2)
        sin_ot = self.builder.makeNoiseElement(noise_spec, suffix="_S1_S2")
        self.builder.makeAdditionStatement(filter_in, sin_ot)
        self.check()

    def testMakeControlErrorSignal(self):
        if IGNORE_TEST:
            return
        self.init()
        def test(sign):
            signal_ot = self.builder.makeControlErrorSignal(-7, "S3", sign, suffix="_S1_S3")
            if sign == -1:
                result = re.search("-7.*S3", self.getStatement())
            else:
                result = re.search("S3.*-7", self.getStatement())
            self.assertTrue(result)
            self.check()
        #
        test(-1)
        test(1)

    def testMakePIDController(self):
        if IGNORE_TEST:
            return
        self.init()
        name_in, name_ot = self.builder.makePIDControllerElement("S3", kP=7, kD=5, suffix="_S1_S3")
        self.builder.makeBoundarySpecies("S1")
        self.builder.makeAdditionStatement("S1", name_ot)
        self.builder.makeAdditionStatement(name_in, 3, "-"+"S3")
        self.check()
    
    def testMakePIDController4(self):
        if IGNORE_TEST:
            return
        self.init()
        name_in, name_ot = self.builder.makePIDControllerElement("S3", kP=7, suffix="_S1_S3")
        self.builder.makeBoundarySpecies("S1")
        self.builder.makeAdditionStatement("S1", name_ot)
        self.builder.makeAdditionStatement(name_in, 3, "-"+"S3")
        self.check()

    def testMakePIDController3(self):
        # Filter without differential control
        if IGNORE_TEST:
            return
        self.init()
        name_in, name_ot = self.builder.makePIDControllerElement("S3", kP=7, suffix="_S1_S3")
        self.builder.makeBoundarySpecies("S1")
        self.builder.makeAdditionStatement("S1", name_ot)
        self.builder.makeAdditionStatement(name_in, 3, "-"+"S3")
        self.check()

    def testMakePIDController2(self):
        # Test interaction between filter and differential control
        if IGNORE_TEST:
            return
        self.init()
        name_in, name_ot = self.builder.makePIDControllerElement("S3", kP=7, kD=5, suffix="_S1_S3",
                    filter_derivative_calculation="-3*4")
        self.builder.makeBoundarySpecies("S1")
        self.builder.makeAdditionStatement("S1", name_ot)
        self.builder.makeAdditionStatement(name_in, 3, "-"+"S3")
        self.check()

    def testMakePIDControllerInputIsAParameter(self):
        if IGNORE_TEST:
            return
        builder = ab.AntimonyBuilder(LINEAR_MDL, symbol_dct=SYMBOL_DCT)
        name_in, name_ot = builder.makePIDControllerElement("S1", kP=7, suffix="_S1_S3")
        builder.makeAdditionStatement("k0", name_ot)
        builder.makeAdditionStatement(name_in, 3, "-"+"S3")
        self.check(builder=builder)

    def testClosedLoopSymbols(self):
        if IGNORE_TEST:
            return
        def test(name, suffix=None):
            if suffix is None:
                symbols = [s for s in self.builder.closed_loop_symbols if name in s]
            else:
                symbols = [s for s in self.builder.closed_loop_symbols if (name in s) and ("ot" in s)]
            self.assertTrue(len(symbols) > 0)
        #
        self.init()
        self.builder.makeBoundarySpecies("S1")
        noise_spec = cn.NoiseSpec(sine_amp=1, sine_freq=2)
        self.builder.makeSISOClosedLoopSystem("S1", "S3", kP=10, kI=1, kD=2, kF=10e5, setpoint=5, noise_spec=noise_spec,
                                              disturbance_spec=cn.NoiseSpec(sine_amp=2, sine_freq=3))
        self.assertGreater(len(self.builder.closed_loop_symbols), 0)
        for prefix in ["noise", "disturbance"]:
            test(prefix, suffix="out")
        for prefix in ["filter", "control_error", "controller"]:
            test(prefix, suffix="in")
            test(prefix, suffix="out")
        for prefix in ["kP", "kI", "kD", "setpoint"]:
            test(prefix)

    def testMakeSISOClosedLoopSystem(self):
        if IGNORE_TEST:
            return
        self.init()
        self.builder.makeBoundarySpecies("S1")
        noise_spec = cn.NoiseSpec(sine_amp=1, sine_freq=2)
        self.builder.makeSISOClosedLoopSystem("S1", "S3", kP=10, kI=1, setpoint=5, noise_spec=noise_spec,
                                              disturbance_spec=cn.NoiseSpec(sine_amp=2, sine_freq=3))
        self.check()
        #
        self.builder.initializeOutput()
        self.builder.makeBoundarySpecies("S1")
        noise_spec = cn.NoiseSpec(sine_amp=1, sine_freq=2)
        self.builder.makeSISOClosedLoopSystem("S1", "S3", kP=1, kI=0.1, kF=0.1, kD=2, setpoint=5,
                                              noise_spec=noise_spec)
        self.check()
        #
        self.builder.initializeOutput()
        self.builder.makeBoundarySpecies("S1")
        self.builder.makeSISOClosedLoopSystem("S1", "S3", kP=10, setpoint=5, kI=5, kD=2,
                                              noise_spec=noise_spec)
        self.check()
        #
        self.builder.initializeOutput()
        self.builder.makeBoundarySpecies("S1")
        self.builder.makeSISOClosedLoopSystem("S1", "S3", kP=10, setpoint=5, noise_spec=noise_spec,
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
        self.builder.makeSISOClosedLoopSystem("S1", "S3", kP=10, kI=1, setpoint=4)
        data = self.check(times=np.linspace(0, 100, 1000))
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

    def testGetClosedLoopSymbols(self):
        if IGNORE_TEST:
            return
        input_name = "S1"
        output_name = "S3"
        builder = ab.AntimonyBuilder(LINEAR_MDL, symbol_dct=SYMBOL_DCT)
        builder.makeBoundarySpecies("S1")
        noise_spec = cn.NoiseSpec(sine_amp=1, sine_freq=2, random_mag=3, random_std=4, offset=5, slope=0.0006)
        builder.makeSISOClosedLoopSystem(input_name, output_name, kP=10, kI=1, kF=10e5, setpoint=5, noise_spec=noise_spec,
                                              disturbance_spec=cn.NoiseSpec(sine_amp=2, sine_freq=3))
        def test(name):
            symbols = builder._getClosedLoopSymbols(input_name, output_name)
            self.assertGreaterEqual(len([s for s in symbols if name in s]), 1)
        #
        for name in ["kP", "kI", "kF", "noise", "disturbance", "setpoint"]:
            test(name)
       

if __name__ == '__main__':
  unittest.main()
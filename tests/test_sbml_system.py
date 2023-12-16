from controlSBML.sbml_system import SBMLSystem
from controlSBML import util
import controlSBML.constants as cn
from controlSBML.timeseries import Timeseries

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import unittest


IGNORE_TEST = False
IS_PLOT = False
IMPROPER_LINEAR_MDL = """
// Illustrate Antimony File
species S1, S2, S3, S4

aa := S1
bb := S4

S1 -> S2; k1*S1
J1: S2 -> S3; k2*S2
J2: S3 -> S2; k3*S3
J3: S2 -> S4; k4*S2

k1 = 1
k2 = 1
k3 = 1
k4 = 1
S1 = 10
S2 = 0
S3 = 0
"""
LINEAR_MDL = "model *main_model()\n" + IMPROPER_LINEAR_MDL + "\nend"
SPECIES_NAMES = ["S1", "S2", "S3", "S4"]


#############################
# Tests
#############################
class TestSBMLSystem(unittest.TestCase):

    def setUp(self):
        self.system = SBMLSystem(LINEAR_MDL, ["S1"], ["S3"], is_fixed_input_species=True)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.system.original_antimony, str))

    def testImproperModel(self):
        if IGNORE_TEST:
            return
        sbml_system = SBMLSystem(IMPROPER_LINEAR_MDL, ["S1"], ["S3"], is_fixed_input_species=True)
        with self.assertRaises(ValueError):
            _ = sbml_system.antimony_builder
                          
    def notestConstructor2(self):
        if IGNORE_TEST:
                return
        with self.assertRaises(ValueError):
            _ = SBMLSystem(LINEAR_MDL, ["S9"], ["S3"])
        with self.assertRaises(ValueError):
            _ = SBMLSystem(LINEAR_MDL, ["S1"], ["S9"])

    def testConstructor3(self):
        if IGNORE_TEST:
                return
        self.system = SBMLSystem(LINEAR_MDL, ["aa"], ["bb"], is_fixed_input_species=True)

    def testGet(self):
        if IGNORE_TEST:
            return
        for name in SPECIES_NAMES:
            self.assertEqual(self.system.get(name), self.system._roadrunner[name])

    def testSet(self):
        if IGNORE_TEST:
            return
        for name in SPECIES_NAMES:
            self.system.set({name: 100})
            self.assertEqual(self.system.get(name), 100)

    def testSimulate1(self):
        if IGNORE_TEST:
            return
        system = SBMLSystem(LINEAR_MDL, ["S1"], ["S4"], is_fixed_input_species=True)
        ts = system.simulate(start_time=0, end_time=100)
        self.assertTrue(np.isclose(ts["S4"].values[-1], 10))
       
    def testSimulateSISOClosedLoop(self):
        if IGNORE_TEST:
            return
        def test(is_fixed_input_species):
            system = SBMLSystem(LINEAR_MDL, ["S1"], ["S3"], is_fixed_input_species=is_fixed_input_species)
            setpoint = 5
            ts, _ = system.simulateSISOClosedLoop("S1", "S3", kp=2, ki=0.8, kf=0.5, setpoint=setpoint, end_time=200, num_point=1000)
            self.assertGreater(len(ts), 0)
            if is_fixed_input_species:
                tolerance = 0.1
            else:
                tolerance = 0.1
            self.assertLess(np.abs(ts["S3"].values[-1] - setpoint), tolerance)
            if IS_PLOT:
                for column in ts.columns:
                    if column not in ["time", "S1", "S3"]:
                        del ts[column]
                if is_fixed_input_species:
                    title = "Boundary species"
                else:
                    title = "Boundary reaction"
                util.plotOneTS(ts, ax2=0, figsize=(8,8), title=title)
                plt.show()
        #
        test(True)
        test(False)

    def testSimulateStaircase(self):
        if IGNORE_TEST:
            return
        times = np.linspace(0, 50, 500)
        def test(is_fixed_input_species):
            system = SBMLSystem(LINEAR_MDL, ["S1"], ["S3"], is_fixed_input_species=True)
            ts, _ = system.simulateStaircase("S1", "S3", times=times, final_value=10, num_step=5, is_steady_state=False)
            self.assertGreater(len(ts), 0)
            variance = np.var(ts["S3"])
            self.assertFalse(np.isclose(variance, 0))
            if IS_PLOT:
                util.plotOneTS(ts, ax2=0, figsize=(8,8), ylim=[0, 10])
                plt.show()
        #
        test(True)
        test(False)

    def testCopyAndEqual(self):
        if IGNORE_TEST:
            return
        system = self.system.copy()
        self.assertTrue(system == self.system)
        #
        times = np.linspace(0, 50, 500)
        _ = system.simulateStaircase("S1", "S3", times=times, final_value=10, num_step=5, is_steady_state=False)
        self.assertFalse(system == self.system)

    def testPlotSISOClosedLoop(self):
        #if IGNORE_TEST:
        #    return
        system = SBMLSystem(LINEAR_MDL, ["S1"], ["S3"], is_fixed_input_species=False)
        setpoint = 5
        ts, _ = system.simulateSISOClosedLoop(input_name="S1", output_name="S3", kp=2, ki=0.8, kf=0.5,
                                           setpoint=setpoint, end_time=100, num_point=1000)
        self.system.plotSISOClosedLoop(ts, setpoint, figsize=(5,5), title="Closed Loop", is_plot=IS_PLOT)

    def testPlotModel(self):
        if IGNORE_TEST:
            return
        system = SBMLSystem(LINEAR_MDL)
        ts = system.plotModel(is_plot=IS_PLOT)
        self.assertTrue(isinstance(ts, Timeseries))

    def testGetValidSymbolsInput(self):
        if IGNORE_TEST:
            return
        input_ser = self.system.getValidSymbols(is_input=True, is_str=False)
        self.assertTrue(isinstance(input_ser, pd.Series))
        self.assertTrue("S1" in input_ser.loc[cn.TYPE_FLOATING_SPECIES])
        self.assertTrue("k1" in input_ser.loc[cn.TYPE_PARAMETER])
        self.assertFalse("J1" in input_ser.index)
        #
        input_str = self.system.getValidInputs()
        self.assertTrue(isinstance(input_str, str))

    def testGetValidSymbolsOutput(self):
        if IGNORE_TEST:
            return
        output_ser = self.system.getValidSymbols(is_input=False, is_str=False)
        self.assertTrue(isinstance(output_ser, pd.Series))
        self.assertTrue("S1" in output_ser.loc[cn.TYPE_FLOATING_SPECIES])
        self.assertTrue("k1" in output_ser.loc[cn.TYPE_PARAMETER])
        self.assertFalse("J1" in output_ser.index)
        #
        out_str = self.system.getValidOutputs()
        self.assertTrue(isinstance(out_str, str))


if __name__ == '__main__':
  unittest.main()
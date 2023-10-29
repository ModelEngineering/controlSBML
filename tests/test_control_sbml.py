from controlSBML.control_sbml import ControlSBML
from controlSBML.siso_transfer_function_builder import SISOTransferFunctionBuilder
from controlSBML.timeseries import Timeseries
from controlSBML.antimony_builder import AntimonyBuilder

import numpy as np
import unittest


IGNORE_TEST = False
IS_PLOT = False
FIGSIZE = (5, 5)
LINEAR_MDL = """
model *main_model()
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
end
"""
SPECIES_NAMES = ["S1", "S2", "S3", "S4"]


#############################
# Tests
#############################
class TestControlSBML(unittest.TestCase):

    def setUp(self):
        self.ctlsb = ControlSBML(LINEAR_MDL)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue("RoadRunner" in str(type(self.ctlsb._roadrunner)))

    def testSetTimes(self):
        if IGNORE_TEST:
            return
        max_value = 10
        self.ctlsb.setTimes(np.linspace(0, max_value, 100))
        times = self.ctlsb.getTimes()
        self.assertEqual(times[-1], max_value)

    def testGetSetSystem(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSISOSystemSpecification(input_name="S1", output_name="S2",
                                              is_fixed_input_species=False, is_steady_state=False)
        system, builder = self.ctlsb.getSystem()
        self.assertTrue("S1" in system.input_names)
        self.assertTrue(isinstance(builder, SISOTransferFunctionBuilder))

    def testPlotModel(self):
        if IGNORE_TEST:
            return
        ts = self.ctlsb.plotModel(is_plot=IS_PLOT, figsize=FIGSIZE)
        self.assertTrue(isinstance(ts, Timeseries))

    def testPlotStaircaseResponse(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSISOSystemSpecification(input_name="S1", output_name="S3")
        ts, builder = self.ctlsb.plotStaircaseResponse(is_plot=IS_PLOT, figsize=FIGSIZE)
        self.assertTrue(isinstance(ts, Timeseries))
        self.assertTrue(isinstance(builder, AntimonyBuilder))

    def testPlotTransferFunctionFit(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSISOSystemSpecification(input_name="S1", output_name="S3")
        ts, builder = self.ctlsb.plotTransferFunctionFit(num_numerator=2, num_denominator=2, is_plot=IS_PLOT, figsize=FIGSIZE)
        self.assertTrue(isinstance(ts, Timeseries))
        self.assertTrue(isinstance(builder, AntimonyBuilder))

    def testPlotSISOClosedLoopSystem(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSISOSystemSpecification(input_name="S1", output_name="S3")
        ts, builder = self.ctlsb.plotSISOClosedLoopSystem(setpoint=3, is_plot=IS_PLOT, kp=1, figsize=FIGSIZE)
        self.assertTrue(isinstance(ts, Timeseries))
        self.assertTrue(isinstance(builder, AntimonyBuilder))
    
    def testPlotSISOClosedLoopDesign(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSISOSystemSpecification(input_name="S1", output_name="S3")
        _ = self.ctlsb.plotTransferFunctionFit(num_numerator=2, num_denominator=2, is_plot=False)
        ts, builder = self.ctlsb.plotSISOClosedLoopDesign(setpoint=5, sign=-1, kp=True, ki=True, figsize=FIGSIZE)
        self.assertTrue(isinstance(ts, Timeseries))
        self.assertTrue(isinstance(builder, AntimonyBuilder))
    
    def testCopy(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSISOSystemSpecification(input_name="S1", output_name="S3")
        _ = self.ctlsb.plotTransferFunctionFit(num_numerator=2, num_denominator=2, is_plot=False)
        ctlsb = self.ctlsb.copy()
        ts, builder = ctlsb.plotSISOClosedLoopDesign(setpoint=5, sign=-1, kp=True, ki=True, figsize=FIGSIZE, is_plot=IS_PLOT)
        self.assertTrue(isinstance(ts, Timeseries))
        self.assertTrue(isinstance(builder, AntimonyBuilder))


if __name__ == '__main__':
  unittest.main()
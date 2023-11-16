from controlSBML.control_sbml import ControlSBML
from controlSBML import constants as cn
from controlSBML.sbml_system import SBMLSystem
from controlSBML.siso_transfer_function_builder import SISOTransferFunctionBuilder
from controlSBML.timeseries import Timeseries
from controlSBML.antimony_builder import AntimonyBuilder

import control
import matplotlib.pyplot as plt
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
CTLSB = ControlSBML(LINEAR_MDL)


#############################
# Tests
#############################
class TestControlSBML(unittest.TestCase):

    def setUp(self):
        self.ctlsb = CTLSB.copy()

    def testConstructor(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(LINEAR_MDL)
        self.assertTrue("RoadRunner" in str(type(ctlsb._roadrunner)))

    def testSetOptions(self):
        if IGNORE_TEST:
            return
        max_value = 10
        self.ctlsb.setOptions(times=np.linspace(0, max_value, 100))
        times = self.ctlsb.getTimes()
        self.assertEqual(times[-1], max_value)

    def testGetOptions(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S2",
                                              is_fixed_input_species=False, is_steady_state=False)
        system, builder = self.ctlsb.getSystem()
        self.assertTrue("S1" in system.input_names)
        self.assertTrue(isinstance(builder, SISOTransferFunctionBuilder))

    def testSetSystem(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S2",
                                              is_fixed_input_species=False, is_steady_state=False)
        system, builder = self.ctlsb.getSystem()
        self.assertTrue("S1" in system.input_names)
        self.assertTrue(isinstance(builder, SISOTransferFunctionBuilder))

    def testPlotModel(self):
        if IGNORE_TEST:
            return
        ts = self.ctlsb.plotModel(is_plot=IS_PLOT, figsize=FIGSIZE)
        if IS_PLOT:
            plt.savefig("test_plot_model.png")
        self.assertTrue(isinstance(ts, Timeseries))

    def testPlotStaircaseResponse(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S3")
        ts, builder = self.ctlsb.plotStaircaseResponse(is_plot=IS_PLOT, figsize=FIGSIZE,
                                                       times=np.linspace(0, 100, 1000))
        self.assertTrue(isinstance(ts, Timeseries))
        self.assertTrue(isinstance(builder, AntimonyBuilder))
        if IS_PLOT:
            plt.savefig("test_plot_staircase_response.png")

    def testPlotTransferFunctionFit(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S3")
        ts, builder = self.ctlsb.plotTransferFunctionFit(num_numerator=1, num_denominator=3, is_plot=IS_PLOT, 
                                                         figsize=FIGSIZE, times=np.linspace(0, 100, 1000))
        self.assertTrue(isinstance(ts, Timeseries))
        self.assertTrue(isinstance(builder, AntimonyBuilder))

    def testPlotSISOClosedLoopSystem(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S3")
        ts, builder = self.ctlsb.plotClosedLoop(setpoint=3, is_plot=IS_PLOT, kp=1, figsize=FIGSIZE,
                                                          times=np.linspace(0, 100, 1000))
        self.assertTrue(isinstance(ts, Timeseries))
        self.assertTrue(isinstance(builder, AntimonyBuilder))
    
    def testPlotDesign(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S3")
        _ = self.ctlsb.plotTransferFunctionFit(num_numerator=2, num_denominator=2, is_plot=False)
        ctlsb = self.ctlsb.copy()
        ts, builder = ctlsb.plotDesign(setpoint=5, sign=-1, kp_spec=True, ki_spec=True, figsize=FIGSIZE, is_plot=IS_PLOT)
        self.assertTrue(isinstance(ts, Timeseries))
        self.assertTrue(isinstance(builder, AntimonyBuilder))

    def testEqualsCopy(self):
        if IGNORE_TEST:
            return
        ctlsb = self.ctlsb.copy()
        self.assertTrue(self.ctlsb.equals(ctlsb))
        ctlsb.setSystem(input_name="S1", output_name="S3")
        self.assertFalse(self.ctlsb.equals(ctlsb))

    def testGetters(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S3")
        self.ctlsb.plotTransferFunctionFit(is_plot=False)
        self.assertTrue(isinstance(self.ctlsb.getSystem()[0], SBMLSystem))
        self.assertTrue(isinstance(self.ctlsb.getSystem()[1], SISOTransferFunctionBuilder))
        self.assertTrue(isinstance(self.ctlsb.getAntimony(), str))
        self.assertTrue(isinstance(self.ctlsb.getInputName(), str))
        self.assertTrue(isinstance(self.ctlsb.getOutputName(), str))
        self.assertTrue(isinstance(self.ctlsb.getClosedLoopTransferFunction(), control.TransferFunction))
        self.assertTrue(isinstance(self.ctlsb.getOptions(), dict))
        self.assertTrue(isinstance(self.ctlsb.getFitterResult(), cn.FitterResult))
        self.assertTrue(isinstance(self.ctlsb.getTimes(), np.ndarray))


if __name__ == '__main__':
  unittest.main()
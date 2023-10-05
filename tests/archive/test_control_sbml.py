from controlSBML.control_sbml import ControlSBML
from controlSBML.staircase import Staircase
from controlSBML import util
import controlSBML.constants as cn
import helpers

import control
import os
import pandas as pd
import tellurium as te
import unittest


IGNORE_TEST = False
IS_PLOT = False

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
ANTIMONY_FILE = os.path.join(TEST_DIR, "Model_antimony.ant")
MODEL_FILE = os.path.join(TEST_DIR, "BIOMD0000000823.xml")
LINEAR_MDL = """
J0: $S0 -> S1; k0*$S0
J1: S1 -> S2; k1*S1
J2: S2 -> S3; k2*S2

S0 = 1
S1 = 10
S2 = 0
S3 = 0
k0 = 1
k1 = 1
k2 = 1
"""
helpers.setupPlotting(__file__)


#############################
# Tests
#############################
class TestControlSBML(unittest.TestCase):

    def setUp(self):
      # Cannot modify self.control
      self.ctlsb = ControlSBML(ANTIMONY_FILE)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue("RoadRunner" in str(type(self.ctlsb.roadrunner)))

    def testConstructWithRoadrunner(self):
        if IGNORE_TEST:
            return
        model = te.loada(helpers.TEST_PATH_1)
        ctlsb1 = ControlSBML(model)
        ctlsb2 = ControlSBML(helpers.TEST_PATH_1)
        diff = set(ctlsb1.get().values()).symmetric_difference(
              ctlsb2.get().values())
        self.assertEqual(len(diff), 0)

    def testMakeSISOTransferFunctionBuilder(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(MODEL_FILE,
            input_names=["IR"], output_names=["mTORC1_DEPTOR"])
        builder = ctlsb.makeSISOTransferFunctionBuilder()
        self.assertTrue("Builder" in str(type(builder)))
    
    def testMakeMIMOTransferFunctionBuilder(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(MODEL_FILE,
            input_names=["IR", "mTORC2_DEPTOR"], output_names=["mTORC1_DEPTOR"])
        builder = ctlsb.makeMIMOTransferFunctionBuilder()
        self.assertTrue("Builder" in str(type(builder)))

    def testMakeMIMOTransferFunctionDF(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(MODEL_FILE,
            input_names=["mTORC2_DEPTOR", "IR"], output_names=["mTORC1_DEPTOR"])
        transfer_function_df = ctlsb.makeMIMOTransferFunctionDF()
        self.assertTrue("DataFrame" in str(type(transfer_function_df)))
        trues = transfer_function_df.applymap(lambda x: isinstance(x, control.TransferFunction))
        self.assertTrue(trues.all().all())
    
    def testMakeNonlinearIOSystem(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(LINEAR_MDL, input_names=["S1"])
        non_sys = ctlsb.makeNonlinearIOSystem("tst")
        self.assertTrue("NonlinearIOSystem" in str(type(non_sys)))
    
    def testPlotMIMOStaircaseResponse(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(LINEAR_MDL)
        plot_result_df = ctlsb.plotMIMOStaircaseResponse(
            is_plot=IS_PLOT, figsize=(8,8),
            input_names=["S1"], output_names=["S2", "S3"],
            staircase=Staircase(initial_value=1, final_value=10, num_step=4),
            end_time=100)
        self.assertTrue(isinstance(plot_result_df, pd.DataFrame))
        self.assertTrue(all([isinstance(plot_result_df.loc[i, o], util.PlotResult) 
                             for i in plot_result_df.index for o in plot_result_df.columns]))
        
    def testPlotMIMOStaircaseResponse2(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(LINEAR_MDL)
        plot_result_df = ctlsb.plotMIMOStaircaseResponse(
            is_plot=IS_PLOT, figsize=(8,8),
            input_names=["k0"], output_names=["S2", "J1"],
            staircase=Staircase(initial_value=1, final_value=10, num_step=4),
            end_time=100)
        self.assertTrue(isinstance(plot_result_df, pd.DataFrame))
        self.assertTrue(all([isinstance(plot_result_df.loc[i, o], util.PlotResult) 
                             for i in plot_result_df.index for o in plot_result_df.columns]))

    def testFitMIMOTransferFunction(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(LINEAR_MDL)
        fitter_result_df = ctlsb.fitMIMOTransferFunction(
            num_numerator=1, num_denominator=2,
            is_fixed_input_species=True,
            input_names=["S1"], output_names=["S2", "S3"],
            staircase=Staircase(initial_value=1, final_value=10, num_step=4),
            end_time=100)
        self.assertTrue(isinstance(fitter_result_df, pd.DataFrame))
        self.assertTrue(all([isinstance(fitter_result_df.loc[i, o], cn.FitterResult) 
                             for i in fitter_result_df.index for o in fitter_result_df.columns]))
    
    def testPlotMIMOTransferFunction(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(LINEAR_MDL)
        plot_result_df = ctlsb.plotFitMIMOTransferFunction(
            num_numerator=1, num_denominator=2, is_plot=IS_PLOT, figsize=(8,8),
            is_fixed_input_species=True,
            input_names=["S1"], output_names=["S2", "S3"],
            staircase=Staircase(initial_value=1, final_value=10, num_step=4),
            end_time=100)
        self.assertTrue(isinstance(plot_result_df, pd.DataFrame))
        self.assertTrue(all([isinstance(plot_result_df.loc[i, o], util.PlotResult) 
                             for i in plot_result_df.index for o in plot_result_df.columns]))
        


if __name__ == '__main__':
  unittest.main()
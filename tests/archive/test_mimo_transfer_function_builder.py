import controlSBML as ctl
from controlSBML.control_sbml import ControlSBML
import controlSBML.constants as cn
import controlSBML.mimo_transfer_function_builder as tfb
import controlSBML.util as util
from controlSBML.staircase import Staircase
import helpers

import control
import numpy as np
import os
import pandas as pd
import unittest
import tellurium as te


IGNORE_TEST = False
IS_PLOT = False
END_TIME = 5
DT = 0.01
POINTS_PER_TIME = int(1.0 / DT)
NUM_TIME = int(POINTS_PER_TIME*END_TIME) + 1
TIMES = [n*DT for n in range(0, NUM_TIME)]

LINEAR_MDL = """
J0:  -> S1; k0
J1: S1 -> S2; k1*S1
J2: S2 -> S3; k2*S2
J3: S3 -> S4; k3*S3
J4: S4 -> ; k4*S4

k0 = 50
k1 = 1
k2 = 2
k3 = 3
k4 = 4
S1 = 1
S2 = 1
S3 = 1
S4 = 1
"""
INPUT_NAME = "S1"
OUTPUT_NAME = "S3"
INPUT_NAMES = ["S1", "S2"]
OUTPUT_NAMES = ["S3", "S4"]
#
#PLOT_PATH = helpers.setupPlotting(__file__)

#############################
# Tests
#############################
class TestSBMLTransferFunctionBuilder(unittest.TestCase):

    def setUp(self):
        self.remove()
        self.ctlsb = ControlSBML(LINEAR_MDL, input_names=INPUT_NAMES, output_names=OUTPUT_NAMES)
        self.builder = self.ctlsb.makeMIMOTransferFunctionBuilder(
            input_names=INPUT_NAMES, output_names=OUTPUT_NAMES)
        
    def tearDown(self):
        self.remove()
        
    def remove(self):
        return

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.builder, tfb.MIMOTransferFunctionBuilder))

    def checkTransferFunctions(self, result_df):
        self.assertTrue(isinstance(result_df, pd.DataFrame))
        result_arr = result_df.values
        result_arr = np.array([r.transfer_function for r in result_arr.flatten()])
        self.assertTrue(all([isinstance(tf, control.TransferFunction) for tf in result_arr]))
        return result_arr
    
    def checkTimeseries(self, result_df):
        self.assertTrue(isinstance(result_df, pd.DataFrame))
        for column in result_df.columns:
            for idx in result_df.index:
                ts = result_df.loc[idx, column]
                self.assertTrue(isinstance(ts, ctl.Timeseries))
    
    def testMakeStaircaseResponse(self):
        if IGNORE_TEST:
            return
        result_df = self.builder.makeStaircaseResponse(staircase=Staircase(final_value=10),
              end_time=100)
        self.checkTimeseries(result_df)
    
    def testFitTransferFunction(self):
        if IGNORE_TEST:
            return
        result_df = self.builder.fitTransferFunction(1, 2,
              staircase=Staircase(final_value=10),
              end_time=100)
        result_arr = self.checkTransferFunctions(result_df)
        self.assertGreater(result_arr[0].dcgain(), result_arr[1].dcgain())

    def testFitTransferFunction2(self):
        if IGNORE_TEST:
            return
        ctlsb = ctl.ControlSBML(cn.WOLF_URL, input_names=["at", "s1"], output_names=["s5", "s6"])
        builder = tfb.MIMOTransferFunctionBuilder(ctlsb, is_fixed_input_species=False)
        result_df = builder.fitTransferFunction(3, 3,
                staircase=Staircase(final_value=5),
                end_time=10)
        _ = self.checkTransferFunctions(result_df)
      
    def testPlotStaircaseResponse(self):
        if IGNORE_TEST:
           return
        response_df = self.builder.makeStaircaseResponse(staircase=Staircase(final_value=5), end_time=10)
        result_df = self.builder.plotStaircaseResponse(response_df, is_plot=IS_PLOT)
        self.checkPlotResultDF(result_df)

    def plotLinearTransferFunction(self, is_fixed_input_species):
        builder = self.ctlsb.makeMIMOTransferFunctionBuilder(
            is_fixed_input_species=is_fixed_input_species,
            input_names=INPUT_NAMES, output_names=OUTPUT_NAMES)
        response_df = builder.fitTransferFunction(4, 4,
                                                       staircase=Staircase(initial_value=1, final_value=6, num_step=10),
                                                       end_time=100)
        result_df = builder.plotFitTransferFunction(response_df, is_plot=IS_PLOT, figsize=(8,8))
        self.checkPlotResultDF(result_df)

    def testPlotFitTransferFunctionFixed(self):
        if IGNORE_TEST:
           return
        self.plotLinearTransferFunction(True)

    def testPlotFitTransferFunctionNotFixed(self):
        if IGNORE_TEST:
           return
        self.plotLinearTransferFunction(False)

    def checkPlotResultDF(self, result_df):
        self.assertTrue(isinstance(result_df, pd.DataFrame))
        for column in result_df.columns:
            for idx in result_df.index:
                plot_result = result_df.loc[idx, column]
                if plot_result is None:
                    continue
                self.assertTrue(isinstance(plot_result, util.PlotResult))

    def testPlotStaircaseResponse2(self):
        if IGNORE_TEST:
           return
        staircase_dct = {"at": Staircase(initial_value=-5, final_value=5),
                         "na": Staircase(initial_value=-10, final_value=10)}
        ctlsb = ControlSBML(cn.WOLF_URL, input_names=["at", "na"], output_names=["s6", "s5"])
        builder = ctlsb.makeMIMOTransferFunctionBuilder()
        response_df = builder.makeStaircaseResponse(staircase=staircase_dct)
        result_df = builder.plotStaircaseResponse(response_df, is_plot=IS_PLOT)
        self.checkPlotResultDF(result_df)

if __name__ == '__main__':
  unittest.main()

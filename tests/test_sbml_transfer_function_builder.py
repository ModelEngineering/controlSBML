import controlSBML as ctl
import controlSBML.constants as cn
import controlSBML.sbml_transfer_function_builder as tfb
import controlSBML.util as util
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
PLOT_PATH = helpers.setupPlotting(__file__)

#############################
# Tests
#############################
class TestSBMLTransferFunctionBuilder(unittest.TestCase):

    def setUp(self):
        self.remove()
        self.builder = tfb.SBMLTransferFunctionBuilder.makeTransferFunctionBuilder(
            LINEAR_MDL, input_names=INPUT_NAMES, output_names=OUTPUT_NAMES)
        
    def tearDown(self):
        self.remove()
        
    def remove(self):
        if os.path.isfile(PLOT_PATH):
            os.remove(PLOT_PATH)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.builder, tfb.SBMLTransferFunctionBuilder))
    
    def testFitTransferFunction(self):
        if IGNORE_TEST:
            return
        transfer_function_df = self.builder.fitTransferFunction(1, 2,
                                                                staircase_spec=cn.StaircaseSpec(final_value=10),
              end_time=100)
        self.assertTrue(isinstance(transfer_function_df, pd.DataFrame))
        tfs = transfer_function_df.values.flatten()
        self.assertTrue(all([isinstance(tf, control.TransferFunction) for tf in tfs]))

    def testFitTransferFunction2(self):
        if IGNORE_TEST:
            return
        ctlsb = ctl.ControlSBML("https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000206.2?filename=BIOMD0000000206_url.xml",
              input_names=["at"], output_names=["s5"])
        builder = ctlsb.makeSISOTransferFunctionBuilder()
        fitter_result = builder.fitTransferFunction(1, 2, 
                                                    staircase_spec=cn.StaircaseSpec(final_value=10),
              end_time=100)
        self.assertTrue(isinstance(fitter_result.time_series, ctl.Timeseries))
        if IS_PLOT:
            ctl.plotOneTS(fitter_result.time_series, writefig=True)
      
    def testPlotStaircaseResponse(self):
        if IGNORE_TEST:
           return
        builder = tfb.SBMLTransferFunctionBuilder.makeTransferFunctionBuilder(LINEAR_MDL)
        plot_result_dct = builder.plotStaircaseResponse(is_plot=IS_PLOT, end_time=100,
                                                         staircase_spec=cn.StaircaseSpec(final_value=10))
        self.checkPlotResultDct(plot_result_dct)

    def checkPlotResultDct(self, plot_result_dct):
        self.assertTrue(isinstance(plot_result_dct, dict))
        for plot_result in plot_result_dct.values():
            if not isinstance(plot_result, util.PlotResult):
                import pdb; pdb.set_trace()
            self.assertTrue(isinstance(plot_result, util.PlotResult))

    def testPlotStaircaseResponse2(self):
        if IGNORE_TEST:
           return
        url = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000206.2?filename=BIOMD0000000206_url.xml"
        builder = tfb.SBMLTransferFunctionBuilder.makeTransferFunctionBuilder(url, input_names=["at", "s5"], output_names=["s6"])
        plot_result_dct = builder.plotStaircaseResponse(is_plot=IS_PLOT, end_time=5,
                                                         staircase_spec=cn.StaircaseSpec(final_value=3))
        self.checkPlotResultDct(plot_result_dct)

if __name__ == '__main__':
  unittest.main()

import controlSBML as ctl
import controlSBML.constants as cn
import controlSBML.sbml_transfer_function_builder as tfb

import control
import numpy as np
import os
import pandas as pd
import unittest
import tellurium as te


IGNORE_TEST = True
IS_PLOT = True
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
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PLOT_PATH = os.path.join(TEST_DIR, "test_sbml_transfer_function_builder.pdf")
cn.DEFAULT_DCTS[1].update({cn.O_WRITEFIG: PLOT_PATH}) 
cn.FIG_DCT[cn.O_WRITEFIG] = PLOT_PATH

#############################
# Tests
#############################
class TestSBMLTransferFunctionBuilder(unittest.TestCase):

    def setUp(self):
        self.builder = tfb.SBMLTransferFunctionBuilder.makeTransferFunctionBuilder(
            LINEAR_MDL, input_names=INPUT_NAMES, output_names=OUTPUT_NAMES)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.builder, tfb.SBMLTransferFunctionBuilder))
    
    def testFitTransferFunction(self):
        if IGNORE_TEST:
            return
        transfer_function_df = self.builder.fitTransferFunction(1, 2, final_value=10,
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
        fitter_result = builder.fitTransferFunction(1, 2, final_value=10,
              end_time=100)
        self.assertTrue(isinstance(fitter_result.time_series, ctl.Timeseries))
        if IS_PLOT:
            ctl.plotOneTS(fitter_result.time_series, writefig=True)
      
    def testPlotStaircaseResponse(self):
        #if IGNORE_TEST:
        #   return
        #plot_result = self.builder.plotStaircaseResponse(is_plot=IS_PLOT, end_time=100,
        #                                                 staircase_spec=cn.StaircaseSpec(final_value=10))
        builder = tfb.SBMLTransferFunctionBuilder.makeTransferFunctionBuilder(LINEAR_MDL)
        plot_result = builder.plotStaircaseResponse(is_plot=IS_PLOT, end_time=10000,
                                                         staircase_spec=cn.StaircaseSpec(final_value=10))

if __name__ == '__main__':
  unittest.main()

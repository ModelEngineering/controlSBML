import controlSBML as ctl
import controlSBML.constants as cn
import controlSBML.sbml_transfer_function_builder as tfb

import control
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import unittest
import shutil
import tellurium as te
import tempfile


IGNORE_TEST = False
IS_PLOT = True
END_TIME = 5
DT = 0.01
POINTS_PER_TIME = int(1.0 / DT)
NUM_TIME = int(POINTS_PER_TIME*END_TIME) + 1
TIMES = [n*DT for n in range(0, NUM_TIME)]

LINEAR_MDL = """
J0: $S0 -> S1; k1*S0
J1: S1 -> S2; S1
J2: S2 -> S3; S2
J3: S3 -> S4; S3
J4: S4 -> S5; S4

k1 = 10
S0 = 1
S1 = 10
S2 = 0
S3 = 0
"""
INPUT_NAME = "S1"
OUTPUT_NAME = "S3"
INPUT_NAMES = ["S1", "S2"]
OUTPUT_NAMES = ["S3", "S4"]

#############################
# Tests
#############################
class TestMIMOTransferFunctionBuilder(unittest.TestCase):

    def setUp(self):
        self.builder = tfb.SBMLTransferFunctionBuilder.makeBuilder(
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


if __name__ == '__main__':
  unittest.main()

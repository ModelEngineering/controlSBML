import controlSBML as ctl
import controlSBML.constants as cn
import controlSBML.siso_transfer_function_builder as stb

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
IS_PLOT = False
END_TIME = 5
DT = 0.01
POINTS_PER_TIME = int(1.0 / DT)
NUM_TIME = int(POINTS_PER_TIME*END_TIME) + 1
TIMES = [n*DT for n in range(0, NUM_TIME)]

LINEAR_MDL = """
J0: $S0 -> S1; k1*S0
J1: S1 -> S2; S1
J2: S2 -> S3; S2

k1 = 10
S0 = 1
S1 = 10
S2 = 0
S3 = 0
"""
INPUT_NAME = "S1"
OUTPUT_NAME = "S3"
ctlsb = ctl.ControlSBML(LINEAR_MDL, input_names=["S1"], output_names=["S3"])
LINEAR_SYS = ctlsb.makeNonlinearIOSystem("LINEAR_SYS")
builder = stb.SISOTransferFunctionBuilder(LINEAR_SYS)
plot_response = builder.plotStaircaseResponse(final_value=10, is_plot=False)
LINEAR_TS = plot_response.time_series


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

    def testMakeParameters(self):
        if IGNORE_TEST:
            return
        def test(parameters, prefix, expected_count):
            names = [n for n in parameters.valuesdict().keys()
                  if n[0] == prefix]
            self.assertEqual(len(names), expected_count)
        #
        parameters = stb.makeParameters(3, 4)
        test(parameters, stb.NUMERATOR_PREFIX, 3)
        test(parameters, stb.DENOMINATOR_PREFIX, 4)
        self.assertTrue(len(parameters.valuesdict()) == 7)

    def testMakeTransferFunction(self):
        if IGNORE_TEST:
            return
        parameters = stb.makeParameters(2, 2)
        tf = stb.makeTransferFunction(parameters)
        self.assertTrue(tf.poles()[0] == -1)
        self.assertTrue(tf.dcgain() == 1)

    def testCalculateTransferFunctionResiduals(self):
        if IGNORE_TEST:
            return
        times = LINEAR_TS.times
        data_in = (times, LINEAR_TS["S1_staircase"].values)
        data_out = LINEAR_TS["S3"].values
        parameters = stb.makeParameters(3, 3)
        residuals = stb._calculateTransferFunctionResiduals(parameters, data_in,
              data_out)
        self.assertEqual(len(residuals), len(times))
        self.assertTrue("float" in str(residuals.dtype))


#############################
# Tests
#############################
class TestNonlinearIOSystem(unittest.TestCase):

    def setUp(self):
        if IGNORE_TEST:
            return
        self.init()
        self.removeFiles()

    def tearDown(self):
        plt.close()
        self.removeFiles()

    def init(self, do_simulate_on_update=True):
        if IS_PLOT:
            cn.PLOT_DIR = cn.TEST_DIR
        else:
            cn.PLOT_DIR= tempfile.mkdtemp()
        self.ctlsb = ctl.ControlSBML(LINEAR_MDL,
              input_names=[INPUT_NAME], output_names=[OUTPUT_NAME])
        self.sys = ctl.NonlinearIOSystem("test_sys", self.ctlsb,
               do_simulate_on_update=do_simulate_on_update)
        self.builder = stb.SISOTransferFunctionBuilder(self.sys)

    def removeFiles(self):
        for ffile in os.listdir(cn.PLOT_DIR):
            if ("figure_" in ffile) and (".pdf") in ffile:
                path = os.path.join(cn.PLOT_DIR, ffile)
                if os.path.isfile(path) and IGNORE_TEST:
                    os.remove(path)
        # FIXME: This won't work in windows
        if IS_PLOT and ("var" in cn.PLOT_DIR):
            shutil.rmtree(cn.PLOT_DIR)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.init()
        self.assertTrue(isinstance(self.builder, stb.SISOTransferFunctionBuilder ))

    def testMakeStaircase(self):
        if IGNORE_TEST:
            return
        self.init()
        def test(num_point, num_step, initial_value=0, final_value=5):
            result = self.builder._makeStaircase(num_point,
                  num_step, initial_value, final_value)
            self.assertTrue(len(result), num_point)
            num_distinct = len(set(result))
            self.assertEqual(num_distinct, num_step)
            self.assertEqual(result[0], initial_value)
            self.assertEqual(result[-1], final_value)
        #
        test, (20, 4)
        test(19, 4)
        test(191, 17)
        test(91, 15)

    def testPlotStaircaseResponse(self):
        if IGNORE_TEST:
            return
        self.init()
        staircase_name = "%s_staircase" % self.builder.input_name
        legend_crd = (.5, 1)
        def test(num_step, initial_value=0, final_value=20, start_time=0):
            plot_result = self.builder.plotStaircaseResponse(
                  final_value, num_step=num_step,
                  initial_value=initial_value, end_time=500, start_time=start_time,
                  writefig=True, figsize=(5,5), is_plot=IS_PLOT,
                  legend_crd=legend_crd)
            arr = plot_result.time_series[staircase_name].values
            num_distinct = len(set(arr))
            self.assertEqual(num_distinct, num_step)
            self.assertEqual(arr[0], initial_value)
            self.assertEqual(arr[-1], final_value)
            return plot_result
        #
        result = test(4, start_time=20)
        result = test(4)
        self.assertEqual(str(result), "")
        result = test(17)
        result = test(15)
        test(4, )

    def testPlotStaircaseResponseWithoutPlot(self):
        if IGNORE_TEST:
            return
        self.init()
        name = "S1"
        staircase_name = "%s_staircase" % name
        output_names = ["S2"]
        ctlsb = ctl.ControlSBML(LINEAR_MDL,
              input_names=["S1"], output_names=output_names)
        legend_spec = cn.LegendSpec(output_names, crd=(.5, 1))
        builder = ctlsb.makeSISOTransferFunctionBuilder(input_name=name,
              output_name=["S2"])
        def test(num_step, initial_value=0, final_value=11):
            plot_result = builder.plotStaircaseResponse(num_step=num_step,
                  initial_value=initial_value,
                  final_value=final_value, end_time=200,
                  is_plot=False,
                  legend_spec=legend_spec)
            arr = plot_result.time_series[staircase_name].values
            num_distinct = len(set(arr))
            self.assertEqual(num_distinct, num_step)
            self.assertEqual(arr[0], initial_value)
            self.assertEqual(arr[-1], final_value)
            return plot_result
        #
        result = test(4)
        self.assertEqual(str(result), "")
        result = test(17)
        result = test(15)
        test(4, )

    def testFitTransferFunction(self):
        if IGNORE_TEST:
            return
        self.init()
        fitter_result = self.builder.fitTransferFunction(4, 4, final_value=10,
              end_time=100)
        self.assertTrue(isinstance(fitter_result.time_series, ctl.Timeseries))
        if IS_PLOT:
            ctl.plotOneTS(fitter_result.time_series, writefig=True)

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

import controlSBML.siso_closed_loop_designer as cld
from controlSBML.control_sbml import ControlSBML
from controlSBML.siso_transfer_function_builder import SISOTransferFunctionBuilder
import controlSBML.constants as cn
from controlSBML.grid import Grid, Point
import helpers
from controlSBML.timeseries import Timeseries
from controlSBML.sbml_system import SBMLSystem
from controlSBML.staircase import Staircase
import controlSBML.util as util

import copy
import control
import numpy as np
import os
import pandas as pd
import sympy
import unittest

IGNORE_TEST = False
IS_PLOT = False
FIGSIZE = (5, 5)
helpers.setupPlotting(__file__)
MODEL = """
model *main_model()
S0 -> S1; k0*S0
S1 -> S2; k1*S1

k0 = 1
k1 = 1
S0 = 0
S1 = 0
S2 = 0
end
"""
MODEL2 = """
model *main2_model()
S0 -> S1; k0*S0
S1 -> S2; k1*S1
S2 -> ; k2*S2

k0 = 1
k1 = 1
k2 = 2
S0 = 10
S1 = 0
S2 = 0
end
"""
LINEAR_MDL = """
model *linear_model()
species S3

 -> S1; k0
S1 -> S2; k1*S1
S2 -> S3; k2*S2
S3 -> ; k3*S3

S1 = 0
S2 = 0
S3 = 0
k0 = 0
k1 = 1
k2 = 2
k3 = 3
end
"""
# Construct a transfer function for the model. This is a linear model, and so it should be accurate.
INPUT_NAME = "S0"
OUTPUT_NAME = "S2"
SYSTEM = SBMLSystem(MODEL2, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME], is_fixed_input_species=True)
TRANSFER_FUNCTION = control.TransferFunction(np.array([1.51083121, 2.01413339]), np.array([1.67214802, 1.24125478, 9.99999997]))
TIMES = np.linspace(0, 20, 200)
PARAMETER_DCT = {p: n+1 for n, p in enumerate(cn.CONTROL_PARAMETERS)}
SETPOINT = 3
SAVE_PATH = os.path.join(cn.TEST_DIR, "siso_closed_loop_designer.csv")
SAVE1_PATH = os.path.join(cn.TEST_DIR, "siso_closed_loop_designer1.csv")
REMOVE_FILES = [SAVE_PATH, SAVE1_PATH]
if False:
    # Required to construct the transfer function
    builder = SISOTransferFunctionBuilder(SYSTEM, input_name=INPUT_NAME, output_name=OUTPUT_NAME)
    staircase = ctl.Staircase(final_value=15, num_step=5)
    fitter_result = builder.fitTransferFunction(num_numerator=2, num_denominator=3, staircase=staircase)
    if False:
        builder.plotFitTransferFunction(fitter_result, figsize=(5,5))
    TRANSFER_FUNCTION = fitter_result.transfer_function


#############################
# Tests
#############################
class TestSISOClosedLoopDesigner(unittest.TestCase):

    def setUp(self):
        if IGNORE_TEST:
            return
        self.remove()
        self.init()

    def tearDown(self):
        self.remove()

    def remove(self):
        if IGNORE_TEST:
            return
        for ffile in REMOVE_FILES:
            if os.path.exists(ffile):
                os.remove(ffile)

    def init(self):
        if "sys_tf" in dir(self):
            # Already initialized
            return
        self.sys_tf = copy.deepcopy(TRANSFER_FUNCTION)
        self.system = copy.deepcopy(SYSTEM)
        self.designer = cld.SISOClosedLoopDesigner(self.system, self.sys_tf, times=TIMES, setpoint=SETPOINT,
                                                   save_path=SAVE_PATH)

    def testGetSet(self):
        if IGNORE_TEST:
            return
        self.designer.set(**PARAMETER_DCT)
        for name in cn.CONTROL_PARAMETERS:
            self.assertEqual(getattr(self.designer, name), PARAMETER_DCT[name])
        #
        dct = self.designer.get()
        for name in cn.CONTROL_PARAMETERS:
            self.assertEqual(dct[name], PARAMETER_DCT[name])

    def testCalculateClosedLoopTf(self):
        if IGNORE_TEST:
            return
        sys_tf = control.tf([1], [1, 2])
        closed_loop_tf_kp = cld._calculateClosedLoopTransferFunction(open_loop_transfer_function=sys_tf, kp=3)
        closed_loop_tf_ki = cld._calculateClosedLoopTransferFunction(open_loop_transfer_function=sys_tf, ki=3)
        _, ys_kp = control.step_response(closed_loop_tf_kp, TIMES)
        _, ys_ki = control.step_response(closed_loop_tf_ki, TIMES)
        self.assertTrue(ys_kp[-1] < ys_ki[-1])
        self.assertTrue(np.isclose(ys_ki[-1], 1, atol=0.01))

    def test_closed_loop_tf(self):
        if IGNORE_TEST:
            return
        sys_tf = control.tf([1], [1, 1])
        designer = cld.SISOClosedLoopDesigner(sys_tf)
        with self.assertRaises(ValueError):
            _ = designer.closed_loop_transfer_function()
        #
        designer.kf = 4
        closed_loop_tf = designer.closed_loop_transfer_function
        numr = np.array(closed_loop_tf.num[0][0])
        self.assertTrue(np.allclose(numr, [20000, 80000, 0]))
        denr = np.array(closed_loop_tf.den[0][0])
        self.assertTrue(np.allclose(denr, [1.00000e+00, 1.00050e+04, 1.30004e+05, 4.00000e+04]))

    def testDesign(self):
        if IGNORE_TEST:
            return
        self.init()
        def checkParams(names):
            for name in names:
                self.assertIsNotNone(getattr(designer, name))
            other_names = set(cn.CONTROL_PARAMETERS) - set(names)
            for name in other_names:
                self.assertIsNone(getattr(designer, name))
        designer = cld.SISOClosedLoopDesigner(SYSTEM, self.sys_tf, times=np.linspace(0, 200, 1000))
        designer.design(kp_spec=True, ki_spec=True, num_restart=1, max_value=10)
        param_dct = designer.get()
        designer.evaluate(is_plot=IS_PLOT)
        checkParams(["kp", "ki"])
        #
        designer = cld.SISOClosedLoopDesigner(SYSTEM, self.sys_tf, times=np.linspace(0, 200, 1000), setpoint=5)
        designer.set(**param_dct)
        designer.evaluate(is_plot=IS_PLOT)

    def testSimulate(self):
        if IGNORE_TEST:
            return
        def calcDiff(arr):
            return np.abs(arr[-1] - 1)
        #
        sys_tf = control.tf([1], [1, 1])
        designer = cld.SISOClosedLoopDesigner(SYSTEM, sys_tf)
        designer.set(kp=20)
        _, prediction1s = designer.simulateTransferFunction()
        designer.set(kp=20, ki=50)
        _, prediction2s = designer.simulateTransferFunction()
        self.assertLess(calcDiff(prediction2s), calcDiff(prediction1s))

    def testPlot(self):
        if IGNORE_TEST:
            return
        self.designer.set(**PARAMETER_DCT)
        self.designer.plot(is_plot=IS_PLOT, markers=["", ""])
        self.designer.set(kp=10, ki=5)
        self.designer.plot(is_plot=IS_PLOT)

    def makeDesigner(self, end_time=200):
        times = np.linspace(0, end_time, 10*end_time)
        system = copy.deepcopy(SYSTEM)
        transfer_function = copy.deepcopy(TRANSFER_FUNCTION)
        designer = cld.SISOClosedLoopDesigner(system, transfer_function, times=times, save_path=SAVE_PATH)
        return designer

    def testPlot2(self):
        if IGNORE_TEST:
            return
        designer = self.makeDesigner()
        designer.design(kp_spec=True, ki_spec=True, min_value=0.1, max_value=10, num_restart=1)
        designer.plot(is_plot=IS_PLOT, markers=["", ""])
        self.assertGreater(designer.kp, 0)
        self.assertGreater(designer.ki, 0)
        self.assertIsNone(designer.kf)

    def testPlot3(self):
        if IGNORE_TEST:
            return
        designer = self.makeDesigner()
        designer.design(kp_spec=True, ki_spec=True, kf_spec=True, min_value=0.1, max_value=10, 
                        num_coordinate=5, num_restart=1, is_report=IGNORE_TEST)
        designer.plot(is_plot=IS_PLOT, markers=["", ""])
        self.assertGreater(designer.kp, 0)
        self.assertGreater(designer.ki, 0)
        self.assertGreater(designer.kf, 0)

    def testDesignAlongGrid(self):
        if IGNORE_TEST:
            return
        designer = self.makeDesigner()
        grid = Grid(min_value=0.1, max_value=10, num_coordinate=5)
        grid.addAxis("kp")
        designer.designAlongGrid(grid)
        self.assertGreater(designer.kp, 0)
        self.assertIsNone(designer.ki)
        self.assertIsNone(designer.kf)
    
    def testDesignAlongGridParallel(self):
        if IGNORE_TEST:
            return
        designer = self.makeDesigner()
        grid = Grid(min_value=0.1, max_value=10, num_coordinate=11)
        grid = Grid(min_value=0.1, max_value=10, num_coordinate=5, is_random=False)
        grid.addAxis("kp")
        grid.addAxis("ki")
        designer.designAlongGrid(grid, is_report=IGNORE_TEST)
        self.assertGreater(designer.kp, 0)
        self.assertGreater(designer.ki, 0)
        self.assertIsNone(designer.kf)

    def testPlotDesignResult(self):
        if IGNORE_TEST:
            return
        designer = self.makeDesigner()
        def test(parameter_names):
            dct = {}
            for spec in cn.CONTROL_PARAMETER_SPECS:
                if spec in parameter_names:
                    dct[spec] = True
                else:
                    dct[spec] = False
            designer.design(min_value=0.1, max_value=10, 
                 num_coordinate=4, num_restart=1, is_report=IGNORE_TEST, **dct)
            designer.plotDesignResult(is_plot=IS_PLOT, figsize=(15,15))
        #
        test(["kp_spec", "ki_spec", "kf_spec"])
        test(["kp_spec", "ki_spec"])
        test(["kp_spec"])

    def test_closed_loop_tf(self):
        # Checks that the closed loop transfer function is calculated correctly
        if IGNORE_TEST:
            return
        # Setup
        s, kp, ki = sympy.symbols("s kp ki")
        # System transfer function
        sys_tf = control.tf([1], [1, 1])
        systf = 1/(s + 1)
        # Symbolic calculation of transfer function
        ctltf = kp + ki/s
        looptf = sympy.simplify(systf*ctltf)
        cltf = sympy.simplify(looptf/(1 + looptf))
        #
        designer = cld.SISOClosedLoopDesigner(self.system, sys_tf, setpoint=5)
        designer.set(kp=2, ki=3)
        closed_loop_tf = designer.closed_loop_transfer_function
        func1 = lambda x: np.real(closed_loop_tf(x))
        cltf_nums = cltf.subs({kp: 2, ki: 3})
        func2 = lambda x: sympy.N(cltf_nums.subs({s: x}))
        result = util.compareSingleArgumentFunctions(func1, func2, 0, 100)
        self.assertTrue(result)

    def testBug1(self):
        if IGNORE_TEST:
            return
        # Setup
        system = SBMLSystem(LINEAR_MDL, input_names=["k0"], output_names=["S3"], is_fixed_input_species=True,
                                is_steady_state=False)
        linear_bldr = SISOTransferFunctionBuilder(system)
        linear_staircase = Staircase(initial_value=0, final_value=10, num_step=5)
        fitter_result = linear_bldr.fitTransferFunction(num_numerator=2, num_denominator=3, 
                                                    staircase=linear_staircase, fit_start_time=20,
                                                start_time=0, end_time=200)
        linear_tf = fitter_result.transfer_function
        #
        times = np.linspace(0, 1000, 10000)
        designer = cld.SISOClosedLoopDesigner(system, linear_tf, times=times, setpoint=5)
        designer.design(kp_spec=True, ki_spec=True, num_restart=2, max_value=100)
        designer.evaluate(is_plot=IS_PLOT, figsize=FIGSIZE)
    
    def testBug3(self):
        if IGNORE_TEST:
            return
        # Setup
        url = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000823.2?filename=Varusai2018.xml"
        INPUT_NAME = "pIRS"
        OUTPUT_NAME = "pmTORC1"
        ctlsb = ControlSBML(url, save_path=SAVE1_PATH)
        ctlsb.setOptions(input_names=[INPUT_NAME], output_names=[OUTPUT_NAME])
        #
        grid = ctlsb.getGrid(kp_spec=True, ki_spec=False, num_coordinate=2, is_random=False)
        axis = grid.getAxis("kp")
        axis.setMinValue(0.1)
        axis.setMaxValue(1.1)
        ts, builder = ctlsb.plotGridDesign(grid, setpoint=120, num_restart=1, is_greedy=False, 
                                           is_plot=IS_PLOT, save_path=SAVE1_PATH)
        self.assertTrue(isinstance(ts, Timeseries))

    def testBug4(self):
        if IGNORE_TEST:
            return
        # Setup
        url = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000823.2?filename=Varusai2018.xml"
        INPUT_NAME = "pIRS"
        OUTPUT_NAME = "pmTORC1"
        ctlsb = ControlSBML(url, save_path=SAVE1_PATH)
        ctlsb.setOptions(input_names=[INPUT_NAME], output_names=[OUTPUT_NAME])
        #
        grid = ctlsb.getGrid(kp_spec=True, ki_spec=False, num_coordinate=40, is_random=False)
        axis = grid.getAxis("kp")
        axis.setMinValue(0.1)
        axis.setMaxValue(10)
        ts, builder = ctlsb.plotGridDesign(grid, setpoint=120, num_restart=1, is_plot=IS_PLOT,
                                           is_greedy=False, save_path=SAVE1_PATH)
        self.assertTrue(isinstance(ts, Timeseries))


if __name__ == '__main__':
    unittest.main()
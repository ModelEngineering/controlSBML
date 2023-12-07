import controlSBML.siso_closed_loop_designer as cld
from controlSBML.control_sbml import ControlSBML
from controlSBML.siso_transfer_function_builder import SISOTransferFunctionBuilder
import controlSBML as ctl
import controlSBML.constants as cn
from controlSBML.grid import Grid
import helpers

import copy
import control
import matplotlib.pyplot as plt
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
SYSTEM = ctl.SBMLSystem(MODEL2, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME], is_fixed_input_species=True)
TRANSFER_FUNCTION = control.TransferFunction(np.array([1.51083121, 2.01413339]), np.array([1.67214802, 1.24125478, 9.99999997]))
TIMES = np.linspace(0, 20, 200)
PARAMETER_DCT = {p: n+1 for n, p in enumerate(cn.CONTROL_PARAMETERS)}
SETPOINT = 3
SAVE_PATH = os.path.join(cn.TEST_DIR, "siso_closed_loop_designer.csv")
REMOVE_FILES = []
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
        designer.design(kp_spec=True, ki_spec=True, num_restart=3, max_value=10)
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
                        num_coordinate=5, num_restart=1, is_report=True)
        designer.plot(is_plot=IS_PLOT, markers=["", ""])
        self.assertGreater(designer.kp, 0)
        self.assertGreater(designer.ki, 0)
        self.assertGreater(designer.kf, 0)
    
    def testPlotDesignAlongGrid(self):
        if IGNORE_TEST:
            return
        designer = self.makeDesigner()
        grid = Grid(min_value=0.1, max_value=10, num_coordinate=5)
        grid.addAxis("kp")
        designer.designAlongGrid(grid, is_report=True)
        designer.plot(is_plot=IS_PLOT, markers=["", ""])
        self.assertGreater(designer.kp, 0)
        self.assertIsNone(designer.ki)
        self.assertIsNone(designer.kf)

    def testPlotDesignResult(self):
        if IGNORE_TEST:
            return
        designer = self.makeDesigner()
        grid = Grid(min_value=0.1, max_value=2, num_coordinate=8, is_random=False)
        grid.addAxis("kp", min_value=1, max_value=2, num_coordinate=8)
        grid.addAxis("ki", min_value=0.01, max_value=0.1, num_coordinate=8)
        grid.addAxis("kf", num_coordinate=8)
        if not os.path.isfile(SAVE_PATH):
            designer.designAlongGrid(grid)
        design_result_df = pd.read_csv(SAVE_PATH)
        def test(parameter_names):
            names = list(parameter_names)
            names.append(cn.MSE)
            designer._design_result_df = design_result_df[names]
            designer.plotDesignResult(is_plot=IS_PLOT, figsize=(5,5))
        #
        test(["kp", "ki", "kf"])
        test(["kp", "ki"])
        test(["kp"])

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
        func1 = lambda x: float(closed_loop_tf(x))
        cltf_nums = cltf.subs({kp: 2, ki: 3})
        func2 = lambda x: sympy.N(cltf_nums.subs({s: x}))
        result = ctl.util.compareSingleArgumentFunctions(func1, func2, 0, 100)
        self.assertTrue(result)

    def testBug1(self):
        if IGNORE_TEST:
            return
        # Setup
        system = ctl.SBMLSystem(LINEAR_MDL, input_names=["k0"], output_names=["S3"], is_fixed_input_species=True,
                                is_steady_state=False)
        linear_bldr = SISOTransferFunctionBuilder(system)
        linear_staircase = ctl.Staircase(initial_value=0, final_value=10, num_step=5)
        fitter_result = linear_bldr.fitTransferFunction(num_numerator=2, num_denominator=3, 
                                                    staircase=linear_staircase, fit_start_time=20,
                                                start_time=0, end_time=200)
        linear_tf = fitter_result.transfer_function
        #
        times = np.linspace(0, 1000, 10000)
        designer = cld.SISOClosedLoopDesigner(system, linear_tf, times=times, setpoint=5)
        designer.design(kp_spec=True, ki_spec=True, num_restart=2, max_value=100)
        designer.evaluate(is_plot=IS_PLOT, figsize=FIGSIZE)

    def testBug2(self):
        if IGNORE_TEST:
            return
        # Setup
        url = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000823.2?filename=Varusai2018.xml"
        INPUT_NAME = "pIRS"
        OUTPUT_NAME = "pmTORC1"
        ctlsb = ControlSBML(url)
        ctlsb.setOptions(input_names=[INPUT_NAME], output_names=[OUTPUT_NAME])
        _ = ctlsb.plotClosedLoop(setpoint=80, sign=-1, kp=1, ki=0.001, kf=None, figsize=FIGSIZE,
                                  times=np.linspace(0, 100, 1000), is_plot=False)
        _ = ctlsb.plotDesign(kp_spec=True, ki_spec=True, figsize=FIGSIZE, max_parameter_value=1, is_plot=IS_PLOT)
        self.assertIsNotNone(ctlsb.kp)

    def testFindFeasibleClosedLoopSystem(self):
        if IGNORE_TEST:
            return
        self.init()
        times = np.linspace(0, 100, 1000)
        designer = cld.SISOClosedLoopDesigner(self.system, self.sys_tf, times=times, setpoint=SETPOINT)
        initial_dct = {"kp": 1, "ki": 1}
        fixed_dct = {}
        value_dct = designer._searchForFeasibleClosedLoopSystem(initial_dct)
        designer.set(**value_dct)
        ts, _ = designer.evaluate(is_plot=IS_PLOT)
        self.assertTrue(np.isclose(SETPOINT, ts[OUTPUT_NAME].values[-1], atol=0.1))


#############################
class TestHistory(unittest.TestCase):

    def setUp(self):
        self.sys_tf = TRANSFER_FUNCTION
        self.designer = cld.SISOClosedLoopDesigner(SYSTEM, self.sys_tf)
        self.history = cld._History(self.designer)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertEqual(len(self.history), 0)
        diff = set(cn.CONTROL_PARAMETERS) - set(self.history._dct.keys())
        self.assertEqual(len(diff), 0)

    def addItems(self, count):
        for _ in range(count):
            self.history.add()
    
    def testAdd(self):
        if IGNORE_TEST:
            return
        self.addItems(1)
        self.assertEqual(len(self.history), 1)
        self.addItems(1)
        self.assertEqual(len(self.history), 2)

    def testReport(self):
        if IGNORE_TEST:
            return
        COUNT = 3
        self.addItems(COUNT)
        df = self.history.report()
        self.assertEqual(len(df), COUNT)
    
    def testGet(self):
        if IGNORE_TEST:
            return
        COUNT = 1
        self.addItems(COUNT)
        designer = self.history.get(0)
        self.assertTrue(isinstance(designer, cld.SISOClosedLoopDesigner))

    def testClear(self):
        if IGNORE_TEST:
            return
        self.addItems(3)
        self.history.clear()
        self.assertEqual(len(self.history), 0)


if __name__ == '__main__':
    unittest.main()

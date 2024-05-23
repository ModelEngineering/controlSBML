import controlSBML.siso_closed_loop_designer as cld # type: ignore
from controlSBML.siso_transfer_function_builder import SISOTransferFunctionBuilder # type: ignore  
import controlSBML.constants as cn  # type: ignore
from controlSBML.grid import Grid # type: ignore
import helpers
from controlSBML.timeseries import Timeseries # type: ignore
from controlSBML.sbml_system import SBMLSystem # type: ignore
from controlSBML.staircase import Staircase # type: ignore
import controlSBML.util as util # type: ignore

import copy
import control # type: ignore
import numpy as np
import os
import pandas as pd  # type: ignore
import sympy # type: ignore
from typing import List
import unittest

IGNORE_TEST = False
IS_PLOT = False
FIGSIZE = (5, 5)
#helpers.setupPlotting(__file__)
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
CONTROL_PARAMETERS = ["kP", "kI", "kD", "kF"]
INPUT_NAME = "S0"
OUTPUT_NAME = "S2"
SYSTEM = SBMLSystem(MODEL2, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME], is_fixed_input_species=True)
TRANSFER_FUNCTION = control.TransferFunction(np.array([1.51083121, 2.01413339]), np.array([1.67214802, 1.24125478, 9.99999997]))
TIMES = np.linspace(0, 20, 200)
PARAMETER_DCT = {p: n+1 for n, p in enumerate(CONTROL_PARAMETERS)}
SETPOINT = 3
REMOVE_FILES:List[str] = []
CONTROL_PARAMETER_SPECS = ["kP_spec", "kI_spec", "kD_spec", "kF_spec"]
#if False:
#    # Required to construct the transfer function
#    builder = SISOTransferFunctionBuilder(SYSTEM, input_name=INPUT_NAME, output_name=OUTPUT_NAME)
#    staircase = ctl.Staircase(final_value=15, num_step=5)
#    fitter_result = builder.fitTransferFunction(num_numerator=2, num_denominator=3, staircase=staircase)
#    if False:
#        builder.plotFitTransferFunction(fitter_result, figsize=(5,5))
#    TRANSFER_FUNCTION = fitter_result.transfer_function


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
        self.designer = cld.SISOClosedLoopDesigner(self.system, self.sys_tf, times=TIMES, setpoint=SETPOINT)
                                                   

    def testGetSet(self):
        if IGNORE_TEST:
            return
        self.designer.set(**PARAMETER_DCT)
        for name in CONTROL_PARAMETERS:
            self.assertEqual(getattr(self.designer, name), PARAMETER_DCT[name])
        #
        dct = self.designer.get()
        for name in CONTROL_PARAMETERS:
            self.assertEqual(dct[name], PARAMETER_DCT[name])

    def testCalculateClosedLoopTf(self):
        if IGNORE_TEST:
            return
        sys_tf = control.tf([1], [1, 2])
        closed_loop_tf_kP = cld._calculateClosedLoopTransferFunction(open_loop_transfer_function=sys_tf, kP=3)
        closed_loop_tf_kI = cld._calculateClosedLoopTransferFunction(open_loop_transfer_function=sys_tf, kI=3)
        closed_loop_tf_kP = cld._calculateClosedLoopTransferFunction(open_loop_transfer_function=sys_tf, kP=3, kD=2)
        _, ys_kP = control.step_response(closed_loop_tf_kP, TIMES)
        _, ys_kI = control.step_response(closed_loop_tf_kI, TIMES)
        self.assertTrue(ys_kP[-1] < ys_kI[-1])
        self.assertTrue(np.isclose(ys_kI[-1], 1, atol=0.01))

    def testDesign(self):
        if IGNORE_TEST:
            return
        self.init()
        def checkParams(names):
            for name in names:
                self.assertIsNotNone(getattr(designer, name))
            other_names = set(CONTROL_PARAMETERS) - set(names)
            for name in other_names:
                self.assertIsNone(getattr(designer, name))
        designer = cld.SISOClosedLoopDesigner(SYSTEM, self.sys_tf, times=np.linspace(0, 200, 1000))
        designer.design(kP_spec=True, kI_spec=True, kD_spec=0.1, num_restart=1, max_value=10)
        param_dct = designer.get()
        designer.evaluate(is_plot=IS_PLOT)
        checkParams(["kP", "kI", "kD"])
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
        designer.set(kP=20)
        _, prediction1s = designer.simulateTransferFunction()
        designer.set(kP=20, kI=50)
        designer.set(kP=20, kI=50, kD=2)
        _, prediction2s = designer.simulateTransferFunction()
        self.assertLess(calcDiff(prediction2s), calcDiff(prediction1s))

    def testPlot(self):
        if IGNORE_TEST:
            return
        self.init()
        self.designer.set(**PARAMETER_DCT)
        self.designer.plot(is_plot=IS_PLOT, markers=["", ""])
        self.designer.set(kP=10, kI=5, kD=2)
        self.designer.plot(is_plot=IS_PLOT)

    def makeDesigner(self, end_time=200):
        times = np.linspace(0, end_time, 10*end_time)
        system = copy.deepcopy(SYSTEM)
        transfer_function = copy.deepcopy(TRANSFER_FUNCTION)
        designer = cld.SISOClosedLoopDesigner(system, transfer_function, times=times)
        return designer

    def testPlot2(self):
        if IGNORE_TEST:
            return
        designer = self.makeDesigner()
        designer.design(kP_spec=True, kI_spec=True, min_value=0.1, max_value=10, num_restart=1)
        designer.plot(is_plot=IS_PLOT, markers=["", ""])
        self.assertGreater(designer.kP, 0)
        self.assertGreater(designer.kI, 0)
        self.assertIsNone(designer.kF)

    def testPlot3(self):
        if IGNORE_TEST:
            return
        designer = self.makeDesigner()
        designer.design(kP_spec=True, kI_spec=True, kF_spec=True, min_value=0.1, max_value=10, 
                        num_coordinate=5, num_restart=1, is_report=IGNORE_TEST)
        designer.plot(is_plot=IS_PLOT, markers=["", ""])
        self.assertGreater(designer.kP, 0)
        self.assertGreater(designer.kI, 0)

    def testDesignAlongGrid(self):
        if IGNORE_TEST:
            return
        designer = self.makeDesigner()
        grid = Grid(min_value=0.1, max_value=10, num_coordinate=5)
        grid.addAxis("kP")
        result = designer.designAlongGrid(grid, num_process=2)
        self.assertTrue(isinstance(result.dataframe, pd.DataFrame))
        self.assertGreater(designer.kP, 0)
        self.assertIsNone(designer.kI)
        self.assertIsNone(designer.kD)
        self.assertIsNone(designer.kF)
    
    def testDesignAlongGridParallel(self):
        if IGNORE_TEST:
            return
        designer = self.makeDesigner()
        grid = Grid(min_value=0.1, max_value=10, num_coordinate=11)
        grid = Grid(min_value=0.1, max_value=10, num_coordinate=5, is_random=False)
        grid.addAxis("kP")
        grid.addAxis("kD")
        designer.designAlongGrid(grid, is_report=IGNORE_TEST)
        self.assertGreater(designer.kP, 0)
        self.assertGreater(designer.kD, 0)
        self.assertIsNone(designer.kI)
        self.assertIsNone(designer.kF)

    def testPlotDesignResult(self):
        # Plots a previously computed result
        if IGNORE_TEST:
            return
        designer = self.makeDesigner()
        def test(parameter_names):
            dct = {}
            for spec in CONTROL_PARAMETER_SPECS:
                if spec in parameter_names:
                    dct[spec] = True
                else:
                    dct[spec] = False
            designer.design(min_value=0.1, max_value=10, 
                 num_coordinate=4, num_restart=1, is_report=IGNORE_TEST, **dct)
            designer.plotDesignResult(is_plot=IS_PLOT, figsize=(15,15))
        #
        test(["kP_spec", "kI_spec", "kD_spec"])
        test(["kP_spec", "kI_spec"])
        test(["kP_spec"])

    def testPlotDesignResult2(self):
        # Check for selection of parameters
        if IGNORE_TEST:
            return
        def test(parameter_names):
            designer = self.makeDesigner()
            dct = {}
            for spec in CONTROL_PARAMETER_SPECS:
                if spec in parameter_names:
                    dct[spec] = True
                else:
                    dct[spec] = False
            designer.design(min_value=0.1, max_value=10, 
                 num_coordinate=4, num_restart=1, is_report=IGNORE_TEST, **dct)
            designer.plotDesignResult(is_plot=IS_PLOT, figsize=(15,15))
        #
        test(["kP_spec", "kI_spec", "kD_spec"])
        with self.assertRaises(ValueError):
            test(["kP_spec", "kI_spec", "kD_spec", "kF_spec"])
        test(["kP_spec", "kI_spec"])
        test(["kP_spec"])

    def test_closed_loop_tf(self):
        # Checks that the closed loop transfer function is calculated correctly
        if IGNORE_TEST:
            return
        # Setup
        s, kP, kI = sympy.symbols("s kP kI")
        # System transfer function
        sys_tf = control.tf([1], [1, 1])
        systf = 1/(s + 1)
        # Symbolic calculation of transfer function
        ctltf = kP + kI/s
        looptf = sympy.simplify(systf*ctltf)
        cltf = sympy.simplify(looptf/(1 + looptf))
        #
        designer = cld.SISOClosedLoopDesigner(self.system, sys_tf, setpoint=5)
        designer.set(kP=2, kI=3)
        closed_loop_tf = designer.closed_loop_transfer_function
        func1 = lambda x: np.real(closed_loop_tf(x))
        cltf_nums = cltf.subs({kP: 2, kI: 3})
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
        fitter_result = linear_bldr.plotTransferFunctionFit(num_zero=1, num_pole=2, 
                                                    staircase=linear_staircase, fit_start_time=5,
                                                start_time=0, end_time=200, is_plot=IS_PLOT, figsize=FIGSIZE)
        linear_tf = fitter_result.transfer_function
        #
        times = np.linspace(0, 1000, 10000)
        designer = cld.SISOClosedLoopDesigner(system, linear_tf, times=times, setpoint=5)
        designer.design(kP_spec=True, kI_spec=True, num_restart=2, max_value=100)
        designer.evaluate(is_plot=IS_PLOT, figsize=FIGSIZE)
    

if __name__ == '__main__':
    unittest.main()
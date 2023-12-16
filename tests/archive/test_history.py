from controlSBML.archive.history import History
import controlSBML.constants as cn
from controlSBML.siso_closed_loop_designer import SISOClosedLoopDesigner
from controlSBML.sbml_system import SBMLSystem

import control
import numpy as np
import os
import unittest

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

IGNORE_TEST = False
IS_PLOT = False
# Construct a transfer function for the model. This is a linear model, and so it should be accurate.
INPUT_NAME = "S0"
OUTPUT_NAME = "S2"
SYSTEM = SBMLSystem(MODEL2, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME], is_fixed_input_species=True)
TRANSFER_FUNCTION = control.TransferFunction(np.array([1.51083121, 2.01413339]), np.array([1.67214802, 1.24125478, 9.99999997]))
TIMES = np.linspace(0, 20, 200)
PARAMETER_DCT = {p: n+1 for n, p in enumerate(cn.CONTROL_PARAMETERS)}
SETPOINT = 3



#############################
# Tests
#############################
class TestHistory(unittest.TestCase):

    def setUp(self):
        self.sys_tf = TRANSFER_FUNCTION
        self.designer = SISOClosedLoopDesigner(SYSTEM, self.sys_tf)
        self.history = History(self.designer)

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
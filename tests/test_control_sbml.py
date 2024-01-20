from controlSBML.control_sbml import ControlSBML # type: ignore
from controlSBML import constants as cn # type: ignore
from controlSBML.sbml_system import SBMLSystem # type: ignore
from controlSBML.siso_transfer_function_builder import SISOTransferFunctionBuilder # type: ignore
from controlSBML.timeseries import Timeseries # type: ignore
from controlSBML.antimony_builder import AntimonyBuilder # type: ignore
from controlSBML.grid import Grid # type: ignore

import control # type: ignore
import matplotlib.pyplot as plt
import numpy as np
import os
import unittest


IGNORE_TEST = False
IS_PLOT = False
FIGSIZE = (5, 5)
SAVE1_PATH = os.path.join(cn.TEST_DIR, "control_sbml_save_path.csv")
LINEAR_MDL = """
model *main_model()
// Illustrate Antimony File
species S1, S2, S3, S4

aa := S1
bb := S4

S1 -> S2; k1*S1
J1: S2 -> S3; k2*S2
J2: S3 -> S2; k3*S3
J3: S2 -> S4; k4*S2

k1 = 1
k2 = 1
k3 = 1
k4 = 1
S1 = 10
S2 = 0
S3 = 0
end
"""
SPECIES_NAMES = ["S1", "S2", "S3", "S4"]
CTLSB = ControlSBML(LINEAR_MDL, final_value=10)
CSV_FILE1 = os.path.join(cn.TEST_DIR, "test_control_sbml1.csv")
CSV_FILE2 = os.path.join(cn.TEST_DIR, "test_control_sbml2.csv")
REMOVE_FILES = [CSV_FILE1, CSV_FILE2]


#############################
# Tests
#############################
class TestControlSBML(unittest.TestCase):

    def setUp(self):
        self.ctlsb = CTLSB.copy()
        self.remove()

    def tearDown(self):
        self.remove()

    def remove(self):
        for file in REMOVE_FILES:
            if os.path.isfile(file):
                os.remove(file)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(LINEAR_MDL)
        self.assertTrue("RoadRunner" in str(type(ctlsb._roadrunner)))

    def testSetOptions(self):
        if IGNORE_TEST:
            return
        max_value = 10
        self.ctlsb.setOptions(times=np.linspace(0, max_value, 100))
        times = self.ctlsb.getTimes()
        self.assertEqual(times[-1], max_value)

    def testGetOptions(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S2",
                                              is_fixed_input_species=False, is_steady_state=False)
        system, builder = self.ctlsb.getSystem()
        self.assertTrue("S1" in system.input_names)
        self.assertTrue(isinstance(builder, SISOTransferFunctionBuilder))

    def testSetSystem(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S2",
                                              is_fixed_input_species=False, is_steady_state=False)
        system, builder = self.ctlsb.getSystem()
        self.assertTrue("S1" in system.input_names)
        self.assertTrue(isinstance(builder, SISOTransferFunctionBuilder))

    def testPlotModel(self):
        if IGNORE_TEST:
            return
        ts = self.ctlsb.plotModel(is_plot=IS_PLOT, figsize=FIGSIZE)
        ts = self.ctlsb.plotModel(is_plot=IS_PLOT, figsize=FIGSIZE, times=np.linspace(0, 100, 1000))
        ts = self.ctlsb.plotModel(is_plot=IS_PLOT, figsize=FIGSIZE, markers=False)
        if IS_PLOT:
            plt.savefig("test_plot_model.png")
        self.assertTrue(isinstance(ts, Timeseries))

    def testPlotStaircaseResponse(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S3")
        ts, builder = self.ctlsb.plotStaircaseResponse(is_plot=IS_PLOT, figsize=FIGSIZE,
                                                       times=np.linspace(0, 100, 1000))
        self.assertTrue(isinstance(ts, Timeseries))
        self.assertTrue(isinstance(builder, AntimonyBuilder))
        if IS_PLOT:
            plt.savefig("test_plot_staircase_response.png")

    def testPlotTransferFunctionFit(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S3")
        ts, builder = self.ctlsb.plotTransferFunctionFit(num_numerator=1, num_denominator=3, is_plot=IS_PLOT, 
                                                         figsize=FIGSIZE, times=np.linspace(0, 100, 1000))
        self.assertTrue(isinstance(ts, Timeseries))
        self.assertTrue(isinstance(builder, AntimonyBuilder))

    def testPlotSISOClosedLoop(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S3")
        ts, builder = self.ctlsb.plotClosedLoop(setpoint=3, is_plot=IS_PLOT, kP=1, figsize=FIGSIZE,
                                                          times=np.linspace(0, 100, 1000))
        self.assertTrue(isinstance(ts, Timeseries))
        self.assertTrue(isinstance(builder, AntimonyBuilder))
    
    def testPlotDesign(self):
        if IGNORE_TEST:
            return
        setpoint = 5
        self.ctlsb.setSystem(input_name="S1", output_name="S3")
        ts, builder = self.ctlsb.plotDesign(setpoint=setpoint, sign=-1, kP_spec=True, kI_spec=True, figsize=FIGSIZE, is_plot=IS_PLOT,
                                            min_parameter_value=0.001, max_parameter_value=10, num_restart=2,
                                            num_coordinate=5, num_process=1)
        # Show that kP, kI are now the defaults
        _ = self.ctlsb.plotClosedLoop(setpoint=setpoint, is_plot=IS_PLOT, kP=1, figsize=FIGSIZE,
                                                          times=np.linspace(0, 100, 1000))
        self.assertTrue(isinstance(ts, Timeseries))
        self.assertTrue(isinstance(builder, AntimonyBuilder))

    def testPlotGridDesign(self):
        if IGNORE_TEST:
            return
        setpoint = 5
        ctlsb = ControlSBML(LINEAR_MDL, setpoint=setpoint, final_value=10, input_names=["S1"], output_names=["S3"])
        grid = Grid(min_value=0.1, max_value=10, num_coordinate=5)
        grid.addAxis("kp")
        _ = ctlsb.plotGridDesign(grid, is_plot=IS_PLOT)

    def testPlotDesign1(self):
        if IGNORE_TEST:
            return
        setpoint = 5
        ctlsb = ControlSBML(LINEAR_MDL, final_value=10, input_names=["S1"], output_names=["S3"], save_path=CSV_FILE1)
        _ = ctlsb.plotDesign(setpoint=setpoint, sign=-1, kP_spec=True, kI_spec=False, is_plot=IS_PLOT,
                                            min_parameter_value=0.001, max_parameter_value=10, num_restart=1,
                                            num_coordinate=2)
        self.assertTrue(os.path.isfile(CSV_FILE1))

    def testEqualsCopy(self):
        if IGNORE_TEST:
            return
        ctlsb = self.ctlsb.copy()
        self.assertTrue(self.ctlsb.equals(ctlsb))
        ctlsb.setSystem(input_name="S1", output_name="S3")
        self.assertFalse(self.ctlsb.equals(ctlsb))

    def testGetters(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S3")
        self.ctlsb.setOptions(kP=5, sign=-1)
        self.ctlsb.plotTransferFunctionFit(is_plot=False)
        self.assertTrue(isinstance(self.ctlsb.getSystem()[0], SBMLSystem))
        self.assertTrue(isinstance(self.ctlsb.getSystem()[1], SISOTransferFunctionBuilder))
        self.assertTrue(isinstance(self.ctlsb.getAntimony(), str))
        self.assertTrue(isinstance(self.ctlsb.getInputName(), str))
        self.assertTrue(isinstance(self.ctlsb.getOutputName(), str))
        self.assertTrue(isinstance(self.ctlsb.getClosedLoopTransferFunction(), control.TransferFunction))
        self.assertTrue(isinstance(self.ctlsb.getOptions(), dict))
        self.assertTrue(isinstance(self.ctlsb.getFitterResult(), cn.FitterResult))
        self.assertTrue(isinstance(self.ctlsb.getTimes(), np.ndarray))

    def testFullAPI(self):
        if IGNORE_TEST:
            return
        TIMES = np.linspace(0, 10000, 10000)
        INPUT_NAME = "pIRS"
        OUTPUT_NAME = "pmTORC1"
        INPUT_NAME = "pIRS"
        URL = "https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1909250003/2/Varusai2018.xml"
        CTLSB = ControlSBML(URL, figsize=FIGSIZE, times=TIMES, markers=False)  # Specify default value of options
        ts = CTLSB.plotModel(ax2=0, is_plot=IS_PLOT)
        # Define the system and plot response over a controlled range
        CTLSB = ControlSBML(URL, figsize=FIGSIZE, input_names=[INPUT_NAME], output_names=[OUTPUT_NAME],
                        times=np.linspace(0, 10000, 10000), markers=False, sign=-1, is_plot=False)
        if True:
            _, builder = CTLSB.plotStaircaseResponse(is_plot=IS_PLOT)
            _, builder = CTLSB.plotStaircaseResponse(initial_value=20, final_value=25, is_plot=IS_PLOT)
            _ = CTLSB.plotTransferFunctionFit(figsize=FIGSIZE, num_numerator=2, num_denominator=3, initial_value=20, final_value=25,
                                            fit_start_time=2000, is_plot=IS_PLOT)
            _ = CTLSB.plotClosedLoop(setpoint=150, kP=1, kF=None, is_plot=IS_PLOT)
        ts, builder = CTLSB.plotDesign(setpoint=150, kP_spec=True, kI_spec=True, kF_spec=False, 
                                       num_restart=1, is_plot=IS_PLOT)
        _ = CTLSB.plotClosedLoop(setpoint=120, kP=0.002, kI=0.019, is_plot=IS_PLOT)
        _ = CTLSB.plotClosedLoop(setpoint=150, kP=1, is_plot=IS_PLOT)

    def testPlotDesignResult(self):
        if IGNORE_TEST:
            return
        setpoint = 5
        ctlsb = ControlSBML(LINEAR_MDL, final_value=10, input_names=["S1"], output_names=["S3"], save_path=CSV_FILE2)
        _ = ctlsb.plotDesign(setpoint=setpoint, sign=-1, kP_spec=True, kI_spec=True, is_plot=False,
                                            min_parameter_value=0.001, max_parameter_value=10, num_restart=1,
                                            num_coordinate=4)
        plt.close('all')
        _, ax = plt.subplots(1)
        ctlsb.plotDesignResult(is_plot=IS_PLOT, ax=ax)

    def testPerformance(self):
        if IGNORE_TEST:
            return
        URL = "https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1909250003/2/Varusai2018.xml"
        ctlsb = ControlSBML(URL, figsize=(5, 5), times=np.linspace(0, 2000, 20000), markers=False)
        INPUT_NAME = "pIRS"
        OUTPUT_NAME = "pmTORC1"
        ctlsb.setSystem(input_name=INPUT_NAME, output_name=OUTPUT_NAME)
        _, builder = ctlsb.plotTransferFunctionFit(num_numerator=2, num_denominator=3, initial_value=20, final_value=25,
                                  fit_start_time=1000, times=np.linspace(0, 10000, 100000),
                                  is_plot=IS_PLOT)
        
    def testBug1(self):
        if IGNORE_TEST:
            return
        URL = "https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1909250003/2/Varusai2018.xml"
        ctlsb = ControlSBML(URL, figsize=(5, 5), times=np.linspace(0, 2000, 20000), markers=False)
        INPUT_NAME = "pIRS"
        OUTPUT_NAME = "pmTORC1"
        ctlsb.setSystem(input_name=INPUT_NAME, output_name=OUTPUT_NAME)
        _ = ctlsb.plotTransferFunctionFit(num_numerator=1, num_denominator=2, initial_value=20, final_value=25,
                                  fit_start_time=2000, times=np.linspace(0, 10000, 100000), is_plot=IS_PLOT)
        
    def testGetPossibleInputs(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(LINEAR_MDL)
        input_str = ctlsb.getPossibleInputs()
        self.assertTrue(isinstance(input_str, str))
        self.assertTrue(len(input_str) > 0)
        
    def testGetPossibleOutputs(self):
        if IGNORE_TEST:
            return
        ctlsb = ControlSBML(LINEAR_MDL)
        output_str = ctlsb.getPossibleOutputs()
        self.assertTrue(isinstance(output_str, str))
        self.assertTrue(len(output_str) > 0)

    def testBug2(self):
        if IGNORE_TEST:
            return
        model = """
        S1 -> S2; k*S1
        S2 -> S1; k2*S2
        S2 -> ; k3*S2

        S1 = 10
        S2 = 0
        k = 1
        k2 = 0.5
        k3 = 0.5
        """
        ctlsb = ControlSBML(model)
        with self.assertRaises(ValueError):
            _ = ctlsb.plotTransferFunctionFit(num_numerator=1, num_denominator=2)
    
    def testBug3(self):
        #if IGNORE_TEST:
        #    return
        # Setup
        url = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000823.2?filename=Varusai2018.xml"
        INPUT_NAME = "pIRS"
        OUTPUT_NAME = "pmTORC1"
        ctlsb = ControlSBML(url, save_path=SAVE1_PATH)
        ctlsb.setOptions(input_names=[INPUT_NAME], output_names=[OUTPUT_NAME])
        #
        grid = ctlsb.getGrid(kP_spec=True, kI_spec=False, num_coordinate=2, is_random=False)
        axis = grid.getAxis("kP")
        axis.setMinValue(0.1)
        axis.setMaxValue(1.1)
        ts, _ = ctlsb.plotGridDesign(grid, setpoint=120, num_restart=1, is_greedy=False, 
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
        grid = ctlsb.getGrid(kP_spec=True, kI_spec=False, num_coordinate=40, is_random=False)
        axis = grid.getAxis("kP")
        axis.setMinValue(0.1)
        axis.setMaxValue(10)
        ts, builder = ctlsb.plotGridDesign(grid, setpoint=120, num_restart=1, is_plot=IS_PLOT,
                                           is_greedy=False, save_path=SAVE1_PATH)
        self.assertTrue(isinstance(ts, Timeseries))



if __name__ == '__main__':
  unittest.main()
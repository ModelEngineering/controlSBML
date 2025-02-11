from controlSBML.control_sbml import ControlSBML # type: ignore
from controlSBML import constants as cn # type: ignore
from controlSBML.sbml_system import SBMLSystem # type: ignore
from controlSBML.siso_transfer_function_builder import SISOTransferFunctionBuilder # type: ignore
from controlSBML.timeseries import Timeseries # type: ignore
from controlSBML.antimony_builder import AntimonyBuilder # type: ignore
from controlSBML.grid import Grid # type: ignore
import controlSBML.util as util  # type: ignore

import control # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import numpy as np
import os
from typing import List
import unittest


IGNORE_TEST = False
IS_PLOT = False
TIMES = cn.TIMES
FIGSIZE = (5, 5)
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
CTLSB = ControlSBML(LINEAR_MDL, input_name="S1", output_name="S3")
REMOVE_FILES:List[str] = []


#############################
# Tests
#############################
class TestControlSBML(unittest.TestCase):

    def setUp(self):
        if IGNORE_TEST:
            return
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

    def testGetOptions(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S2",
                                              is_fixed_input_species=False, is_steady_state=False)
        self.assertTrue("S1" in self.ctlsb._sbml_system.input_names)
        self.assertTrue(isinstance(self.ctlsb._transfer_function_builder, SISOTransferFunctionBuilder))

    def testSetSystem(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S2",
                                              is_fixed_input_species=False, is_steady_state=False)
        self.assertTrue("S1" in self.ctlsb._sbml_system.input_names)
        self.assertTrue(isinstance(self.ctlsb._transfer_function_builder, SISOTransferFunctionBuilder))

    def testPlotModel(self):
        if IGNORE_TEST:
            return
        self.ctlsb = CTLSB.copy()
        result = self.ctlsb.plotModel(is_plot=IS_PLOT, figsize=FIGSIZE, selections=["S2"],
                                      legend=["This is S2"], xlabel="TIME")
        result = self.ctlsb.plotModel(is_plot=IS_PLOT, figsize=FIGSIZE, selections=["S1", "S2", "S3"])
        result = self.ctlsb.plotModel(is_plot=IS_PLOT, figsize=FIGSIZE, times=np.linspace(0, 100, 1000))
        result = self.ctlsb.plotModel(is_plot=IS_PLOT, figsize=FIGSIZE, markers=False)
        self.assertTrue(isinstance(result.timeseries, Timeseries))

    def testDraw(self):
        if IGNORE_TEST:
            return
        ctlsb = CTLSB.copy()
        ctlsb.draw()

    def testPlotStaircaseResponse(self):
        if IGNORE_TEST:
            return
        self.ctlsb = CTLSB.copy()
        self.ctlsb.setSystem(input_name="S1", output_name="S3")
        result = self.ctlsb.plotStaircaseResponse(is_plot=IS_PLOT, figsize=FIGSIZE,
                                                       times=np.linspace(0, 100, 1000))
        self.assertTrue(isinstance(result.timeseries, Timeseries))
        self.assertTrue(isinstance(result.antimony_builder, AntimonyBuilder))

    def testPlotTransferFunctionFit(self):
        if IGNORE_TEST:
            return
        ctlsb = CTLSB.copy()
        result = ctlsb.plotTransferFunctionFit(num_zero=0, num_pole=2,
                                                         figsize=FIGSIZE, times=np.linspace(0, 100, 1000),
                                                         fitter_method="poly", is_plot=IS_PLOT)
        self.assertTrue(isinstance(result.timeseries, Timeseries))
        self.assertTrue(isinstance(result.antimony_builder, AntimonyBuilder))

    def testPlotSISOClosedLoop(self):
        if IGNORE_TEST:
            return
        ctlsb = CTLSB.copy()
        ctlsb.setSystem(input_name="S1", output_name="S3")
        ts, builder = ctlsb._plotClosedLoop(setpoint=3, is_plot=IS_PLOT, kP=1, figsize=FIGSIZE,
                                                          times=np.linspace(0, 100, 1000))
        self.assertTrue(isinstance(ts, Timeseries))
        self.assertTrue(isinstance(builder, AntimonyBuilder))
    
    def testPlotDesign(self):
        if IGNORE_TEST:
            return
        ctlsb = CTLSB.copy()
        setpoint = 5
        ctlsb.setSystem(input_name="S1", output_name="S3")
        design_result = ctlsb.plotDesign(setpoint=setpoint, kP_spec=True, kI_spec=True, figsize=FIGSIZE, is_plot=IS_PLOT,
                                            min_parameter_value=0.001, max_parameter_value=10, num_restart=2,
                                            num_coordinate=3, num_process=10)
        # Show that kP, kI are now the defaults
        _ = ctlsb._plotClosedLoop(setpoint=setpoint, is_plot=IS_PLOT, kP=1, figsize=FIGSIZE,
                                                          times=np.linspace(0, 100, 1000))
        self.assertTrue(isinstance(design_result.timeseries, Timeseries))
        self.assertTrue(isinstance(design_result.antimony_builder, AntimonyBuilder))

    def testPlotAllDesignResults(self):
        if IGNORE_TEST:
            return
        ctlsb = CTLSB.copy()
        setpoint = 5
        ctlsb.setSystem(input_name="S1", output_name="S3")
        _ = ctlsb.plotDesign(setpoint=setpoint, kP_spec=True, kI_spec=True, kD_spec=True,
                             figsize=FIGSIZE, is_plot=False,
                             min_parameter_value=0.001, max_parameter_value=10,
                             num_restart=1,
                             num_coordinate=3,
        )
        ctlsb.plotDesignResults(is_plot=IS_PLOT, columns=["kD", "kI", "kP"], num_top=15,
                                   round_digit=4)
     
    def testPlotDesignNoiseDisturbance(self):
        if IGNORE_TEST:
            return
        noise_spec = cn.NoiseSpec(sine_amp=0.1, sine_freq=10, random_mag=0.1, random_std=0.1)
        ctlsb = ControlSBML(LINEAR_MDL, input_name="S1", output_name="S3", noise_spec=noise_spec)
        ctlsb = CTLSB.copy()
        setpoint = 5
        ctlsb.setSystem(input_name="S1", output_name="S3")
        design_result = ctlsb.plotDesign(setpoint=setpoint, kP_spec=True, kI_spec=True, figsize=FIGSIZE, is_plot=IS_PLOT,
                                            min_parameter_value=0.001, max_parameter_value=10, num_restart=2,
                                            num_coordinate=5, num_process=10)
        # Show that kP, kI are now the defaults
        _ = ctlsb._plotClosedLoop(setpoint=setpoint, is_plot=IS_PLOT, kP=1, figsize=FIGSIZE,
                                                          times=np.linspace(0, 100, 1000))
        self.assertTrue(isinstance(design_result.timeseries, Timeseries))
        self.assertTrue(isinstance(design_result.antimony_builder, AntimonyBuilder))

    def testPlotDesignKwargs(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S3")
        for kwarg in ["kP", "kI", "kF"]:
            with self.assertRaises(ValueError):
                kwargs = {kwarg: 1}
                _ = self.ctlsb.plotDesign(**kwargs)

    def testPlotGridDesign(self):
        if IGNORE_TEST:
            return
        setpoint = 5
        ctlsb = ControlSBML(LINEAR_MDL, setpoint=setpoint, final_value=10, input_name="S1", output_name="S3")
        grid = Grid(min_value=0.1, max_value=10, num_coordinate=3)
        grid.addAxis("kP", num_coordinate=3)
        grid.addAxis("kI", num_coordinate=3)
        grid.addAxis("kD", num_coordinate=3)
        _ = ctlsb.plotGridDesign(grid, is_plot=IS_PLOT)

    def testEqualsCopy(self):
        if IGNORE_TEST:
            return
        ctlsb = self.ctlsb.copy()
        self.assertTrue(self.ctlsb.equals(ctlsb))
        ctlsb.setSystem(input_name="S1", output_name="S2")
        self.assertFalse(self.ctlsb.equals(ctlsb))

    def testGetters(self):
        if IGNORE_TEST:
            return
        self.ctlsb.setSystem(input_name="S1", output_name="S3")
        self.ctlsb.plotTransferFunctionFit(is_plot=False)
        self.assertTrue(isinstance(self.ctlsb._sbml_system, SBMLSystem))
        self.assertTrue(isinstance(self.ctlsb._transfer_function_builder, SISOTransferFunctionBuilder))
        self.assertTrue(isinstance(self.ctlsb.getAntimony(), str))
        self.assertTrue(isinstance(self.ctlsb.getClosedLoopTransferFunction(), control.TransferFunction))
        self.assertTrue(isinstance(self.ctlsb.getOptions(), dict))

    def testFullAPI(self):
        if IGNORE_TEST:
            return
        INPUT_NAME = "pIRS"
        OUTPUT_NAME = "pmTORC1"
        INPUT_NAME = "pIRS"
        path = util.getModelPath("Varusai2018")
        if False:
            CTLSB = ControlSBML(path, figsize=FIGSIZE, times=TIMES, markers=False,
                                input_name=INPUT_NAME, output_name=OUTPUT_NAME)  # Specify default value of options
            _ = CTLSB.plotModel(ax2=0, is_plot=IS_PLOT)
        # Define the system and plot response over a controlled range
        CTLSB = ControlSBML(path, figsize=FIGSIZE, input_name=INPUT_NAME, output_name=OUTPUT_NAME,
                         times=np.linspace(0, 1000, 10000),
                         markers=False, sign=-1, is_plot=False)
        _ = CTLSB.plotDesign(setpoint=150, kP_spec=True, kI_spec=True, kF_spec=False, 
                                       num_restart=1, is_plot=IS_PLOT, selections=[INPUT_NAME, OUTPUT_NAME],
                                       num_process=5)
        if False:
            _, builder = CTLSB.plotStaircaseResponse(is_plot=IS_PLOT)
            _, builder = CTLSB.plotStaircaseResponse(initial_value=20, final_value=25, is_plot=IS_PLOT)
            _ = CTLSB.plotTransferFunctionFit(figsize=FIGSIZE, num_zero=1, num_pole=2, initial_value=20, final_value=25,
                                            fit_start_time=200, is_plot=IS_PLOT)
            _ = CTLSB._plotClosedLoop(setpoint=150, kP=1, kF=None, is_plot=IS_PLOT)
        _ = CTLSB._plotClosedLoop(setpoint=120, kP=0.002, kI=0.019, is_plot=IS_PLOT)
        _ = CTLSB._plotClosedLoop(setpoint=150, kP=1, is_plot=IS_PLOT)

    def testPlotDesignResult(self):
        if IGNORE_TEST:
            return
        setpoint = 5
        ctlsb = ControlSBML(LINEAR_MDL, final_value=10, input_name="S1", output_name="S3")
        _ = ctlsb.plotDesign(setpoint=setpoint, kP_spec=True, kI_spec=True, kD_spec=True, is_plot=False,
                                            min_parameter_value=0.001, max_parameter_value=10, num_restart=1,
                                            num_coordinate=4)
        plt.close('all')
        _, ax = plt.subplots(1)
        ctlsb._plotDesignResult(is_plot=IS_PLOT, ax=ax)

    def testPerformance(self):
        if IGNORE_TEST:
            return
        path = util.getModelPath("Varusai2018")
        ctlsb = ControlSBML(path, figsize=(5, 5), times=np.linspace(0, 2000, 20000), markers=False)
        INPUT_NAME = "pIRS"
        OUTPUT_NAME = "pmTORC1"
        ctlsb.setSystem(input_name=INPUT_NAME, output_name=OUTPUT_NAME)
        _ = ctlsb.plotTransferFunctionFit(num_zero=1, num_pole=2, initial_value=20, final_value=25,
                                  fit_start_time=1000, times=np.linspace(0, 10000, 100000),
                                  is_plot=IS_PLOT)
        
    def testBug1(self):
        if IGNORE_TEST:
            return
        path = util.getModelPath("Varusai2018")
        ctlsb = ControlSBML(path, figsize=(5, 5), times=np.linspace(0, 2000, 20000), markers=False)
        INPUT_NAME = "pIRS"
        OUTPUT_NAME = "pmTORC1"
        ctlsb.setSystem(input_name=INPUT_NAME, output_name=OUTPUT_NAME)
        _ = ctlsb.plotTransferFunctionFit(num_zero=0, num_pole=1, initial_value=20, final_value=25,
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
            _ = ctlsb.plotTransferFunctionFit(num_zero=0, num_pole=1)
    
    def testBug3(self):
        if IGNORE_TEST:
            return
        # Setup
        path = util.getModelPath("Varusai2018")
        INPUT_NAME = "pIRS"
        OUTPUT_NAME = "pmTORC1"
        ctlsb = ControlSBML(path, input_name=INPUT_NAME, output_name=OUTPUT_NAME)
        #
        grid = ctlsb.getGrid(kP_spec=True, kI_spec=False, kD_spec=True, num_coordinate=10, is_random=False)
        axis = grid.getAxis("kP")
        axis.setMinValue(0.01)
        axis.setMaxValue(0.1)
        result = ctlsb.plotGridDesign(grid, setpoint=120, num_restart=1, is_plot=IS_PLOT)
        self.assertTrue(isinstance(result.timeseries, Timeseries))

    def testBug4(self):
        if IGNORE_TEST:
            return
        # Setup
        path = util.getModelPath("Varusai2018")
        INPUT_NAME = "pIRS"
        OUTPUT_NAME = "pmTORC1"
        ctlsb = ControlSBML(path, input_name=INPUT_NAME, output_name=OUTPUT_NAME)
        #
        grid = ctlsb.getGrid(kP_spec=True, kI_spec=False, num_coordinate=40, is_random=False)
        axis = grid.getAxis("kP")
        axis.setMinValue(0.1)
        axis.setMaxValue(10)
        result = ctlsb.plotGridDesign(grid, setpoint=120, num_restart=1, is_plot=IS_PLOT)
        self.assertTrue(isinstance(result.timeseries, Timeseries))

    def testBug5(self):
        if IGNORE_TEST:
            return
        # Setup
        times = np.linspace(0, 50, 500)
        WOLF_CTLSB = ControlSBML(cn.WOLF_PATH,
                        input_name="s1", output_name="s5", times=times)
        _ = WOLF_CTLSB.plotDesign(kP_spec=True, kI_spec=True, num_restart=1,
                                       is_plot=IS_PLOT, num_coordinate=20,
                                       min_parameter_value=0.001, max_parameter_value=5)

    def testBug6(self):
        # Bug with setting inputs that are fixed
        if IGNORE_TEST:
            return
        model = """
            model *main();
            species S0
            $S1 -> S2; k1*S1
            S2 -> S3; k2*S2
            S3 -> S4; k3*S3
            S4 ->; k4*S4
            S4 + S0 ->; k0*S0*S4

            $S0 = 4
            S1 = 10
            S2 = 0
            S3 = 0
            k1 =1
            k2 =1
            k3 = 1
            k4 = 1
            k0 = 1
            end
            """
        try:
            _ = ControlSBML(model, input_name="S0", output_name="S4") 
            self.assertTrue(True)
        except Exception as e:
            print(e)
            self.assertTrue(False)

    def testBug7(self):
        # Bug with setting inputs that are fixed
        if IGNORE_TEST:
            return
        # Global variables
        INPUT_NAME = "na"
        OUTPUT_NAME = "s5"
        _ = ControlSBML(cn.WOLF_PATH, input_name=INPUT_NAME, output_name=OUTPUT_NAME, times=TIMES)

    def testBug8(self):
        if IGNORE_TEST:
            return
        path = util.getModelPath("Alharbi2019_TNVM")
        ctlsb = ControlSBML(path, times=np.linspace(0, 30, 300), input_name="Vitamins",
                          output_name="Normal_cells", is_fixed_input_species=True)
        design_result = ctlsb.plotDesign(setpoint=2, kP_spec=0.2, kI_spec=0.1, figsize=FIGSIZE, times=np.linspace(0, 100, 1000), min_parameter_value=1,
                max_parameter_value=100,is_plot=IS_PLOT, num_coordinate=2) 
        self.assertTrue(isinstance(design_result.timeseries, Timeseries))
        
    def testBug9(self):
        if IGNORE_TEST:
            return
        path = util.getModelPath("Alharbi2019_TNVM")
        ctlsb = ControlSBML(path, input_name="Vitamins",
                          output_name="Normal_cells", is_fixed_input_species=True, 
                         times=np.linspace(0, 100, 1000), figsize=FIGSIZE)
        def test(kP_spec, kI_spec):
            design_result = ctlsb.plotDesign(setpoint=2, kP_spec=kP_spec, kI_spec=kI_spec, min_parameter_value=1,
                    max_parameter_value=100, num_restart=1, is_plot=IS_PLOT, num_coordinate=2)
            self.assertTrue("kP" in design_result.design_df.columns)
            self.assertTrue("kI" in design_result.design_df.columns)
            if design_result.timeseries is not None:
                self.assertGreater(len(design_result.timeseries.columns), 10)
        #
        test(0.1, 0.1)
        test(0.2, 0.1)

    def testBug10(self):
        if IGNORE_TEST:
            return
        TIMES = np.linspace(000, 15000, 50000)
        path = util.getModelPath("BIOMD0000000571")
        CTLSB = ControlSBML(path, xlabel="time (min)", times=TIMES)
        INPUT_NAME = "Mlc"
        OUTPUT_NAME = "EI_P"
        CTLSB.setSystem(input_name=INPUT_NAME, output_name=OUTPUT_NAME)
        grid = CTLSB.getGrid()
        grid.addAxis("kP", min_value=0.0, max_value=0.005, num_coordinate=3)
        grid.addAxis("kI", min_value=0.0, max_value=0.002, num_coordinate=3)
        design_result = CTLSB.plotGridDesign(grid, setpoint=0.0000003,
                                  is_plot=IS_PLOT)
        self.assertEqual(design_result.design_df.loc[0, cn.REASON], cn.DESIGN_RESULT_SUCCESS)

    def testBug11(self):
        # Bogus initial transient on fit
        if IGNORE_TEST:
            return
        path = util.getModelPath("Tsai2014")
        CTLSB = ControlSBML(path, times=np.linspace(0, 100, 1000))
        CTLSB.setSystem(input_name="Plx1_active", output_name="APC_C_active")
        _ = CTLSB.plotTransferFunctionFit(fit_start_time=20, final_value=1.0,
                                               is_plot=IS_PLOT)
        
    def testBug12(self):
        # Bad x-axis on plot
        if IGNORE_TEST:
            return
        path = util.getModelPath("Varusai2018")
        INPUT_NAME = "pIRS"
        OUTPUT_NAME = "pmTORC1"
        ctlsb= ControlSBML(path, figsize=(5, 5), times=np.linspace(0, 2000, 20000), markers=False,
                   input_name=INPUT_NAME, output_name=OUTPUT_NAME)
        _ = ctlsb.plotDesign(setpoint=80, kP_spec=1, kI_spec=0.01, times=np.linspace(0, 2000, 20000),
                             is_plot=IS_PLOT, num_coordinate=2)
        
    def testBug13(self):
        # Bad time axis
        if IGNORE_TEST:
            return
        path = util.getModelPath("FatehiChenar2018")
        INPUT_NAME = "E"
        OUTPUT_NAME = "R"
        CTLSB = ControlSBML(path, figsize=(5, 5), times=np.linspace(0, 10, 100), markers=False, is_fixed_input_species=True,
                   input_name=INPUT_NAME, output_name=OUTPUT_NAME)  # Specify default value of options
        TIMES = np.linspace(0, 10**5, 10**6)
        _ = CTLSB.plotDesign(setpoint=0.1, kP_spec=1, kI_spec=0.1, times=TIMES, num_restart=1, is_plot=IS_PLOT,
                             num_coordinate=2)

    def testBug14(self):
        # Not using the correct design parameters
        if IGNORE_TEST:
            return
        path = util.getModelPath("Smith2011_V1")
        INPUT_NAME = 'Pneumococci___P'
        OUTPUT_NAME = 'Neutrophils__N'
        CTLSB = ControlSBML(path, input_name=INPUT_NAME, output_name=OUTPUT_NAME, is_fixed_input_species=True,
                       figsize=FIGSIZE)
        TIMES = np.linspace(0, 200, 2000)
        SETPOINT = 1000
        _ = CTLSB.plotDesign(kP_spec=1.5555, kI_spec=0.018, times=TIMES, setpoint=SETPOINT, is_plot=IS_PLOT)
        grid = CTLSB.getGrid()
        grid.addAxis("kP", min_value=0.5, max_value=10, num_coordinate=3)
        grid.addAxis("kI", min_value=0.002, max_value=0.02, num_coordinate=3)
        _ = CTLSB.plotGridDesign(grid, times=TIMES, setpoint=SETPOINT, is_plot=IS_PLOT)

    def testBug15(self):
        # Not using the correct design parameters
        if IGNORE_TEST:
            return
        TIMES = np.linspace(0, 200, 2000)
        INPUT_NAME = 'Pneumococci___P'
        OUTPUT_NAME = 'Neutrophils__N'
        path = util.getModelPath("Smith2011_V1")
        noise_spec = cn.NoiseSpec(random_mag=10.0, random_std=1, offset=1)
        ctlsb = ControlSBML(path, times=TIMES, is_fixed_input_species=True, figsize=FIGSIZE,
                       input_name=INPUT_NAME, output_name=OUTPUT_NAME, noise_spec=noise_spec)
        grid = ctlsb.getGrid()
        grid.addAxis("kP", min_value=0.5, max_value=10, num_coordinate=5)
        grid.addAxis("kD", min_value=0.0, max_value=0.2, num_coordinate=5)
        grid.addAxis("kI", min_value=0.0, max_value=0.2, num_coordinate=5)
        #grid.addAxis("kF", min_value=0.01, max_value=0.02, num_coordinate=2)
        SETPOINT = 1000
        _ = ctlsb.plotGridDesign(grid, times=TIMES, setpoint=SETPOINT, num_restart=5,
                                             is_plot=IS_PLOT)




if __name__ == '__main__':
  unittest.main()
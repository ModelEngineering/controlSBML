"""
Top level API. 

There are 3 types of API calls:
* setters: change the internal state of ControlSBML
* getters: return a value
* plotters: plot a result; return Timeseries of data plotted and AntimonyBuilder used to generate the data; may update state

ControlSBML has "persistent options" that are maintained across calls to the plot* APIs. 
Persistent options relate to staircase definition (initial_value, final_value, num_step), times, closed loop parameters,
and plot options.

Typical Usage:

    # Examine an SBML model
    ctlsb = ControlSBML(model)
    _ = ctlsb.plotModel(times=np.linspace(0, 100, 1000))  # Plot the model

    # Define the system
    ctlsb.getPossibleInputs()
    ctlsb.getPossibleOutputs()
    ctlsb.setSystem(input_names=["S1"], output_names=["S3"])
    ctlsb.input_names, ctlsb.output_names  # Check the system definition

    # Evaluate controllability -- do the inputs affect the outputs?
    _ = ctlsb.plotStaircaseResponse(initial_value=0, final_value=10, num_step=5) # Displays the outputs associated with the inputs

    # Identify the system by fitting to the input staircase
    _ = ctlsb.plotTransferFunctionFit()    # Uses the same staircase function as above

    # Construct a SISO closed loop system and do a trail and error design
    _ = ctlsb.plotClosedLoop(kp=1, ki=0.1, kf=0.5, setpoint=3)
    _ = ctlsb.plotClosedLoop(kp=0.8, ki=0.2, kf=0.1, setpoint=5)

    # Automatically design
    _ = ctlsb.plotSISOClosedLoopDesign(kp=True, ki=True, min_value=0, max_value=10, num_iteration=20, setpoint=4)
    param_dct = {kp: ctlsb.kp, ki: ctlsb.ki, kf: ctlsb.kf}

    # Explore variations on the design
    new_param_dct = {n: v*1.1 for n, v param_dct.items() for n in ["kp", "ki", "kf"]}
    _ = ctlsb.plotClosedLoop(setpoint=5, **new_param_dct)  # Do further trail and error


Conventions for API names. The following prefixes are used for methods:
  get: returns a value
  plot: plots a result; returns Timeseries of data plotted and AntimonyBuilder used to generate the data; may update state
      is_plot: bool (if False, not plot is done. Default is True.)
      Follow matplotlib conventions for plot options: xlabel, ylabel, title, xlim, ylim, markers, etc.
  set: changes internal state of ControlSBML. No value is returned. Little computation is done.
"""

from controlSBML.sbml_system import SBMLSystem
from controlSBML.siso_transfer_function_builder import SISOTransferFunctionBuilder
from controlSBML.siso_closed_loop_designer import SISOClosedLoopDesigner
from controlSBML.timeseries import Timeseries
from controlSBML.staircase import Staircase
from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML import util
from controlSBML.option_management.option_manager import OptionManager
import controlSBML.constants as cn
import controlSBML.msgs as msgs
from controlSBML.option_set import OptionSet
from controlSBML.grid import Grid, Point

import os
import numpy as np
from typing import List, Dict, Tuple, Union

PLOT_KWARGS = list(set(cn.PLOT_KWARGS).union(cn.FIG_KWARGS))
SETPOINT = 1
C_SETPOINT = "setpoint"
TIMES = np.linspace(cn.START_TIME, cn.END_TIME, cn.POINTS_PER_TIME*(cn.END_TIME-cn.START_TIME))
C_INITIAL_VALUE = "initial_value"
C_FINAL_VALUE = "final_value"
C_NUM_STEP = "num_step"
C_TIMES = "times"
C_SIGN = "sign"
C_INPUT_NAMES = "input_names"
C_OUTPUT_NAMES = "output_names"
C_IS_FIXED_INPUT_SPECIES = "is_fixed_input_species"
C_IS_STEADY_STATE = "is_steady_state"
C_MARKERS = "markers"
STAIRCASE_OPTIONS = [C_INITIAL_VALUE, C_FINAL_VALUE, C_NUM_STEP]
TIMES_OPTIONS = [C_TIMES]
CLOSED_LOOP_PARAMETERS = [cn.CP_KP, cn.CP_KI, cn.CP_KF, C_SETPOINT, C_SIGN]
SYSTEM_SPECIFICATION = [C_INPUT_NAMES, C_OUTPUT_NAMES, C_IS_FIXED_INPUT_SPECIES, C_IS_STEADY_STATE]
PLOT_OPTIONS = list(cn.PLOT_KWARGS)
PLOT_OPTIONS.extend(cn.FIG_KWARGS)
PLOT_OPTIONS.append(C_MARKERS)
OPTIONS = STAIRCASE_OPTIONS + TIMES_OPTIONS + CLOSED_LOOP_PARAMETERS + PLOT_OPTIONS + SYSTEM_SPECIFICATION
CONTROL_PARAMETERS = [cn.CP_KP, cn.CP_KI, cn.CP_KF]
FIGSIZE = (5, 5)
INITIAL_OPTION_DCT = {cn.O_TITLE: "", cn.O_SUPTITLE: "", cn.O_WRITEFIG: False,
                        cn.O_XLABEL: "Time (sec)", cn.O_YLABEL: "Concentration",
                        cn.O_FIGSIZE: FIGSIZE,
                        cn.O_IS_PLOT: True,
                        C_IS_FIXED_INPUT_SPECIES: True,
                        C_IS_STEADY_STATE: False,
                        C_SIGN: -1,
                        C_SETPOINT: SETPOINT}
for key in OPTIONS:
    if not key in INITIAL_OPTION_DCT:
        INITIAL_OPTION_DCT[key] = None
SAVE_PATH = os.path.join(cn.DATA_DIR, "control_sbml.csv")


class ControlSBML(OptionSet):

    def __init__(self, model_reference:str, roadrunner=None, save_path:str=None,
                 **kwargs):
        """
        model_reference: str
            string, SBML file or Roadrunner object
        roadrunner: ExtendedRoadrunner
        save_path: str (path to file where results are saved after a design)
        Plot options:
            ax: axis for plot
            figure: figure object
            figsize: figure size (width, height)
            is_plot: bool (plot if True)
            markers: list-str (markers for plot lines; False, no markers)                                                                                       
            suptitle: str (subtitle)
            title: str (title)
            xlabel: str (x label)
            xlim: tupe-float (x lower limit, x upper limit)
            xticklabels: list-float (labels for x ticks)
            ylabel: str (y label)
            ylim: tupe-float (y lower limit, y upper limit)
            yticklabels: list-float (labels for y ticks)
        System options:
            input_names: list-str
            is_steady_state: bool (start system in steady state; default: False)
            is_fixed_input_species: bool (concentration of input species are controlled externally; default: True)
            output_names: list-str
        Closed loop options:
            kf: float (filter constant)
            ki: float (integral control)
            kp: float (proportional control)
            setpoint: float (regulation point)
            sign: -1/+1 (direction of feedback: default: -1) 
        Staircase options:
            final_value: float (last value of input in staircase; default: maximum input value in SBML model)
            initial_value: float (first value of input in staircase; default: minimum input value in SBML model)
            num_step: int (number of steps in staircase; default: 5)
        Miscellaneous options
            times: list-float (times of simulation; default: np.linspace(0, 5, 50))

        """
        # First initializations
        self.model_reference = model_reference
        if roadrunner is None:
            roadrunner = makeRoadrunner(model_reference)
        self._roadrunner = roadrunner
        self.save_path = save_path
        # Internal state
        self._fitter_result = cn.FitterResult()
        # Options
        super().__init__(**INITIAL_OPTION_DCT)
        self.setOptions(times=TIMES)
        self.setOptions(**kwargs)  # Override the base defaults
        initial_value, final_value, num_step = self.initializeStaircaseOptions()
        self.setOptionSet(initial_value=initial_value, final_value=final_value, num_step=num_step)
        self.setOptions(**kwargs)  # Override the base defaults
        
    def copy(self):
        ctlsb = ControlSBML(self.model_reference)
        ctlsb._fitter_result = self._fitter_result.copy()
        persistent_options = self.getOptions()
        ctlsb.setOptions(**persistent_options)
        return ctlsb
    
    def equals(self, other:object):
        if not isinstance(other, ControlSBML):
            return False
        if not self.model_reference == other.model_reference:
            return False
        if not self.getOptions() == other.getOptions():
            return False
        if not self._fitter_result.equals(other._fitter_result):
            return False
        return True

    ############ GETTERS ##############
    def getAntimony(self):
        return self._roadrunner.getAntimony()
    
    def getGrid(self, min_value:int=0, max_value:int=10, num_coordinate:int=10, kp_spec:bool=False,
                ki_spec:bool=False, kf_spec:bool=False, is_random=True)->Grid:
        """
        Creates a grid of values for the control parameters based on the specified control design.

        Args:
            min_value (int): _description_. Defaults to 0.
            max_value (int): _description_. Defaults to 10.
            num_coordinate (int): _description_. Defaults to 10.
            kp_spec (bool): Proportional control. Defaults to False.
            ki_spec (bool): Integral control. Defaults to False.
            kf_spec (bool): Filter. Defaults to False.
            is_random (bool): If True, the grid is randomly generated. Defaults to True.
        Returns:
            Grid
        """
        grid = Grid(min_value=min_value, max_value=max_value, num_coordinate=num_coordinate, is_random=is_random)
        def makeAxis(name, spec):
            if spec:
                grid.addAxis(name, min_value=min_value, max_value=max_value, num_coordinate=num_coordinate)
            return [spec]
        #
        makeAxis(cn.CP_KP, kp_spec)
        makeAxis(cn.CP_KI, ki_spec)
        makeAxis(cn.CP_KF, kf_spec)
        return grid

    def getPossibleInputs(self):
        sbml_system, _ = self.getSystem()
        return sbml_system.getValidInputs()
    
    def getPossibleOutputs(self):
        sbml_system, _ = self.getSystem()
        return sbml_system.getValidOutputs()

    def getInputName(self, option_set:OptionSet=None):
        """
        Args:
            option_set: OptionSet

        Returns:
           str (Name of input)
        """
        if option_set is None:
            input_names = self.input_names
        else:
            input_names = option_set.input_names
        return input_names[0]

    def getOutputName(self, option_set:OptionSet=None):
        """
        Args:
            option_set: OptionSet

        Returns:
           str (Name of output)
        """
        if option_set is None:
            output_names = self.output_names
        else:
            output_names = option_set.output_names
        return output_names[0]
    
    def getOptions(self, options:dict=None):
        """
        Gets current values of the persistent options

        Returns: dict. Keys are listed below by category.
            STAIRCASE_OPTIONS: initial_value, final_value, num_step
            TIMES_OPTIONS: times
            CLOSED_LOOP_PARAMETERS: kp, ki, kf, setpoint, sign
            PLOT_OPTIONS: ax ax2 end_time figure figsize is_plot
                            suptitle title writefig xlabel xlim xticklabels ylabel ylim yticklabels 
        """
        if options is None:
            options = OPTIONS
        return {k: getattr(self, k) for k in options}
    
    def getOpenLoopTransferFunction(self):
        return self._fitter_result.transfer_function
    
    def getFitterResult(self):
        return self._fitter_result

    def getTimes(self):
        return self.times

    def getStaircase(self, **kwargs):
        option_set = self.getOptionSet(**kwargs)
        return Staircase(initial_value=option_set.initial_value, final_value=option_set.final_value,
                         num_step=option_set.num_step, num_point=self._getNumpoint(option_set))
    
    def getSystem(self, **kwargs):
        """
        Returns the SBMLSystem and SISOTransferFunctionBuilder associated with this ControlSBML object.

        Args:
            kwargs: dict (options)

        Returns:
            SBMLSystem, TranferFunctionBuilder
        """
        option_set = self.getOptionSet(**kwargs)
        if option_set.input_names is None:
            input_names = []
        else:
            input_names = option_set.input_names
        if option_set.output_names is None:
            output_names = []
        else:
            output_names = option_set.output_names
        sbml_system = SBMLSystem(self.model_reference,
                                 input_names=input_names,
                output_names=output_names,
                is_fixed_input_species=option_set.is_fixed_input_species,
                is_steady_state=option_set.is_steady_state,
                roadrunner=self._roadrunner)
        if self.input_names is None:
            transfer_function_builder = None
        else:
            transfer_function_builder = SISOTransferFunctionBuilder(sbml_system, input_name=self.getInputName(option_set=option_set),
                output_name=self.getOutputName(option_set=option_set))
        return sbml_system, transfer_function_builder
    
    def getClosedLoopTransferFunction(self, **kwargs):
        """
        Calculates the closed loop transfer function

        Args:
            kwargs: dict (options)

        Returns: 
            control.TransferFunction
        """
        option_set = self.getOptionSet(**kwargs)
        sbml_system, _ = self.getSystem(**option_set)
        if len(sbml_system.input_names) != 1:
            raise ValueError("Must have exactly one input to use this method.")
        if len(sbml_system.output_names) != 1:
            raise ValueError("Must have exactly one output to use this method.")
        designer = SISOClosedLoopDesigner(sbml_system, self.getOpenLoopTransferFunction(),
                setpoint=option_set.setpoint, is_steady_state=sbml_system.is_steady_state,
                sign=option_set.sign)
        kp = 0
        ki = 0
        kf = 0
        for parameter_name in CONTROL_PARAMETERS:
            parameter = getattr(self, parameter_name)
            if parameter is not None:
                if not util.isNumber(parameter):
                    raise ValueError("Must assign float to kp, ki, and kf before using this method.")
                stmt = "%s = %f" % (parameter_name, parameter)
                exec(stmt)
        designer.set(kp=kp, ki=ki, kf=kf)
        return designer.closed_loop_transfer_function
    
    def getParameterStr(self, parameters:List[str], **kwargs):
        """
        Provides a string representation for parameter values.
        Args:
            parameters: list-str (list of parameters)
        Returns:
            str
        """
        stg = ""
        for name in parameters:
            if name in kwargs.keys():
                if util.isNumber(kwargs[name]):
                    stg += "{}={} ".format(name, np.round(kwargs[name], 4))
        return stg

    ############ SETTERS ##############
    def resetOptions(self):
        """
        Resets the persistent options to their default values.
        """
        self.setOptions(**INITIAL_OPTION_DCT)

    def setOptions(self, **kwargs):
        """
        Sets values of persistent options.

        Args:
            kwargs: dict of options
        """
        self.setOptionSet(**kwargs)

    def setSystem(self, input_name:str, output_name:str, is_fixed_input_species:bool=True, is_steady_state:bool=False):
        """
        Sets the options related to defining the system.

        Args:
            input_name: str
            output_name: str
            is_fixed_input_species: bool (if True, a species used as input is treated as fixed)
            is_fixed_steady_state: bool (if True, simulations are started in steady state)
        """
        input_names = [input_name]
        output_names = [output_name]
        self.setOptions(input_names=input_names, output_names=output_names,
                                   is_fixed_input_species=is_fixed_input_species, is_steady_state=is_steady_state)

    def initializeStaircaseOptions(self, initial_value=None, final_value=None, num_step=cn.DEFAULT_NUM_STEP):
        """
        Initializes the staircase options based on their simulation values.

        Args:
            initial_value: float (initial value of the staircase; default is the minimum value of the input if it's a species)
            final_value: float (final value of the staircase)
            num_step: int (number of steps in the staircase)
        """
        # Determine how defaults are assigned
        is_assign_from_simulation = False
        if (initial_value is None) or (final_value is None):
            if self.input_names is not None:
                input_name = self.getInputName()
                if input_name in self.getInputName():
                    is_assign_from_simulation = True
        # Calculate defaults if required
        if initial_value is None:
            if is_assign_from_simulation:
                ts = self.plotModel(is_plot=False)
                initial_value = ts[input_name].min()
            else:
                initial_value = cn.DEFAULT_INITIAL_VALUE
        if final_value is None:
            if is_assign_from_simulation:
                final_value = ts[input_name].max()
            else:
                final_value = cn.DEFAULT_INITIAL_VALUE
        return initial_value, final_value, num_step

    ############ PLOTTERS ##############
    def plotModel(self, **kwargs):
        """
        Plots the original model.

        Args:
            kwargs: dict (persistent options)
        
        Returns:
            Timeseries
        """
        kwargs[cn.O_AX2] = 0
        option_set = self.getOptionSet(**kwargs)
        self._roadrunner.reset()
        data = self._roadrunner.simulate(self._getStarttime(option_set),
                                         self._getEndtime(option_set), self._getNumpoint(option_set))
        ts = Timeseries(data)
        if option_set.is_plot:
            util.plotOneTS(ts, markers=option_set.markers, **self._getPlotOptions(**option_set))
        return ts
    
    def plotStaircaseResponse(self, **kwargs):
        """
        Simulates the staircase response and plots it.

        Args:
            initial_value: float (initial value of the staircase; default is the minimum value of the input if it's a species)
            final_value: float (final value of the staircase)
            num_step: int (number of steps in the staircase)
            kwargs: dict (persistent options)

        Returns:
            Timeseries
                index: time (ms)
                columns: <output_name>, staircase
            AntimonyBuilder
        """
        option_set = self.getOptionSet(**kwargs)
        _, siso_transfer_function_builder = self.getSystem(**option_set)
        staircase = self.getStaircase(**option_set)
        response_ts, builder = siso_transfer_function_builder.makeStaircaseResponse(
            staircase=staircase, is_steady_state=option_set.is_steady_state, times=option_set.times)
        siso_transfer_function_builder.plotStaircaseResponse(response_ts, **self._getPlotOptions(**option_set))
        return response_ts, builder
    
    def plotTransferFunctionFit(self, num_numerator:int=cn.DEFAULT_NUM_NUMERATOR,
                            num_denominator:int=cn.DEFAULT_NUM_DENOMINATOR, fit_start_time:float=None, fit_end_time:float=None, 
                            **kwargs):
        """
        Simulates the staircase response and plots it. Sets the fitter result.

        Args:
            num_numerator: int (number of numerator terms)
            num_denominator: int (number of denominator terms)
            fit_start_time: float (time at which fitting starts)
            fit_end_time: float (time at which fitting ends)
            kwargs: (Persistent options)

        Returns:
            Timeseries (predicted, staircase)
            AntimonyBuilder
        """
        option_set = self.getOptionSet(**kwargs)
        _, siso_transfer_function_builder = self.getSystem(**option_set)
        if fit_start_time is None:
            fit_start_time = option_set.times[0]
        if fit_end_time is None:
            fit_end_time = option_set.times[-1]
        staircase = self.getStaircase(**option_set)
        self._fitter_result = siso_transfer_function_builder.fitTransferFunction(num_numerator=num_numerator,
                            num_denominator=num_denominator, staircase=staircase,
                            fit_start_time=fit_start_time, fit_end_time=fit_end_time, times=option_set.times)
        if self.is_plot:
            siso_transfer_function_builder.plotFitterResult(self._fitter_result, **self._getPlotOptions(**option_set))
        return self._fitter_result.time_series, self._fitter_result.antimony_builder
    
    def plotClosedLoop(self, **kwargs):
        """
        Plots the closed loop response. Control parameters not explicity specified are None.

        Args:
            kwargs: dict (persistent options)
        Returns:
            Timeseries
            AntimonyBuilder
        """
        option_set = self.getOptionSet(**kwargs)
        # Only consider the explicitly specified parameters
        for parameter_name in CONTROL_PARAMETERS:
            if not parameter_name in kwargs.keys():
                option_set.setOptionSet(**{parameter_name: None})
        # Extract the plot options, but not axis 2
        plot_options = self._getPlotOptions(**option_set)
        # Construct the SBML system
        sbml_system, _ = self.getSystem(**option_set)
        response_ts, builder = sbml_system.simulateSISOClosedLoop(input_name=self.getInputName(option_set=option_set),
                                                                  output_name=self.getOutputName(option_set=option_set),
                kp=option_set.kp, ki=option_set.ki, kf=option_set.kf, setpoint=option_set.setpoint, sign=option_set.sign,
                times=option_set.times, is_steady_state=sbml_system.is_steady_state, inplace=False, initial_input_value=None,
                )
        if (not "title" in plot_options) or (len(plot_options["title"]) == 0):
            plot_options["title"] = self.getParameterStr([cn.CP_KP, cn.CP_KI, cn.CP_KF], **option_set)
        _ = plot_options.setdefault(C_MARKERS, False)
        sbml_system.plotSISOClosedLoop(response_ts, option_set.setpoint, **plot_options)
        return response_ts, builder

    def plotGridDesign(self, grid:Grid=None, num_restart:int=1, is_greedy:bool=False, is_plot_grid:bool=False, 
                       save_path:str=None, num_process:int=-1, **kwargs):
        """
        Plots the results of a closed loop design based a grid of values for the control parameters.
        Persists the closed loop design (kp, ki, kf) if a design is found.

        Args:
            grid: Grid (grid of values for the control parameters)
            is_greedy: bool (if True, use greedy algorithm)
            kwargs: dict (persistent options)
            is_plot_grid: bool (if True, plot the grid with mean squared error)
            num_process: int (number of processes to use; -1 means use all available)
        Returns:
            Timeseries
            AntimonyBuilder
        """
        option_set = self.getOptionSet(**kwargs)
        sbml_system, _ = self.getSystem(**option_set)
        if save_path is None:
            save_path = self.save_path
        if save_path is not None:
            if os.path.isfile(save_path):
                os.remove(save_path)
        designer = SISOClosedLoopDesigner(sbml_system, self.getOpenLoopTransferFunction(),
                setpoint=option_set.setpoint, is_steady_state=option_set.is_steady_state, times=option_set.times,
                sign=option_set.sign,
                input_name=self.getInputName(option_set=option_set),
                output_name=self.getOutputName(option_set=option_set), save_path=save_path)
        designer.designAlongGrid(grid, is_greedy=is_greedy, num_process=num_process)
        if designer.residual_mse is None:
            msgs.warn("No design found!")
            return None, None
        self.setOptions(kp=designer.kp, ki=designer.ki, kf=designer.kf)
        # Persist the design parameters
        self.setOptions(kp=designer.kp, ki=designer.ki, kf=designer.kf)
        _ = option_set.setdefault(cn.O_YLABEL, self.getOutputName(option_set=option_set))
        response_ts, antimony_builder = self.plotClosedLoop(
                times=option_set.times,
                setpoint=option_set.setpoint,
                sign=self.sign,
                kp=self.kp,
                ki=self.ki,
                kf=self.kf, 
                **self._getPlotOptions(**option_set))
        return response_ts, antimony_builder

    def plotDesign(self, kp_spec:bool=False, ki_spec:bool=False, kf_spec:bool=False, min_parameter_value:float=0,
                                 max_parameter_value:float=10, num_restart:int=3, 
                                 num_coordinate:int=3, is_greedy:bool=True, is_report:bool=False, 
                                 num_process:int=-1, **kwargs):
        """
        Plots the results of a closed loop design. The design is specified by the parameters kp_spec, ki_spec, and kf_spec.
           None or False: do not include the parameter
           True: include the parameter and find a value
           float: use this value for the parameter.
        Persists the closed loop design (kp, ki, kf) if a design is found.

        Args:
            kp_spec: float (specification of proportional gain)
            ki_spec: float (specification of integral gain)
            kf_spec: float (specification of filter gain)
            min_parameter_value: float (minimum value for kp, ki, kf; may be a dictionary)
            max_parameter_value: float (maximum value for kp, ki, kf; may be a dictionary)
            num_coordinate: int (number of coordinate descent iterations)
            is_greedy: bool (if True, use greedy algorithm)
            kwargs: dict (persistent options)
            is_report: bool (report progress on the design search)
        Returns:
            Timeseries
            AntimonyBuilder
        """
        option_set = self.getOptionSet(**kwargs)
        sbml_system, _ = self.getSystem(**option_set)
        designer = SISOClosedLoopDesigner(sbml_system, self.getOpenLoopTransferFunction(),
                setpoint=option_set.setpoint, is_steady_state=option_set.is_steady_state, times=option_set.times,
                sign=option_set.sign,
                input_name=self.getInputName(option_set=option_set),
                output_name=self.getOutputName(option_set=option_set), save_path=self.save_path)
        designer.design(kp_spec=kp_spec, ki_spec=ki_spec, kf_spec=kf_spec,
                num_restart=num_restart, min_value=min_parameter_value, max_value=max_parameter_value,
            num_coordinate=num_coordinate, is_greedy=is_greedy, is_report=is_report, num_process=num_process)
        if designer.residual_mse is None:
            msgs.warn("No design found!")
            return None, None
        self.setOptions(kp=designer.kp, ki=designer.ki, kf=designer.kf)
        # Persist the design parameters
        self.setOptions(kp=designer.kp, ki=designer.ki, kf=designer.kf)
        _ = option_set.setdefault(cn.O_YLABEL, self.getOutputName(option_set=option_set))
        response_ts, antimony_builder = self.plotClosedLoop(
                times=option_set.times,
                setpoint=option_set.setpoint,
                sign=self.sign,
                kp=self.kp,
                ki=self.ki,
                kf=self.kf, 
                **self._getPlotOptions(**option_set))
        return response_ts, antimony_builder
    
    def plotDesignResult(self, save_path:str=None, **kwargs):
        """
        Plots the mean squared error of all points searched in the last design.

        Args:
            save_path: str (path to CSV file where design results are saved)
            kwargs: dict (plot options)
            AntimonyBuilder
        """
        if save_path is None:
            save_path = self.save_path
        option_set = self.getOptionSet(**kwargs)
        sbml_system, _ = self.getSystem(**option_set)
        designer = SISOClosedLoopDesigner(sbml_system, self.getOpenLoopTransferFunction(), save_path=save_path)
        if designer.design_result_df is None:
            msgs.warn("No design found!")
            return None, None
        return designer.plotDesignResult(**kwargs)

    ############ PROPERTIES AND PRIVATE ##############
    @staticmethod 
    def _getOptions(options:dict)->Tuple[dict, dict]:
        """
        Gets simulation and plot/fig options.

        Args:
            options: dict

        Returns:
            dict (simulation options)
            dict (plot options)
        """
        mgr = OptionManager(options)
        plot_fig_opts = dict(mgr.plot_opts)
        plot_fig_opts.update(mgr.fig_opts)
        plot_fig_opts = {n: v for n, v in plot_fig_opts.items() if v is not None}
        sim_opts = {n: v for n, v in mgr.sim_opts.items() if v is not None}
        return sim_opts, plot_fig_opts
    
    def _getStarttime(self, option_set:OptionSet)->float:
        return option_set.times[0]
    
    def _getEndtime(self, option_set:OptionSet)->float:
        return option_set.times[-1]
    
    def _getNumpoint(self, option_set:OptionSet)->int:
        return len(option_set.times)
    
    def _getPlotOptions(self, **kwargs)->dict:
        option_set = self.getOptionSet(**kwargs)
        plot_options = {n: option_set[n] for n in PLOT_OPTIONS if n!=C_MARKERS}
        return plot_options
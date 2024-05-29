"""
Top level API. 

There are 3 types of API calls:
* setters: change the internal state of ControlSBML
* getters: return a value
* plotters: plot a result; return Timeseries of data plotted and AntimonyBuilder used to generate the data; may update state

Typical Usage:

    # Examine an SBML model
    ctlsb = ControlSBML(model)
    _ = ctlsb.plotModel(times=np.linspace(0, 100, 1000))  # Plot the model

    # Define the system
    ctlsb.getPossibleInputs()
    ctlsb.getPossibleOutputs()
    ctlsb.setSystem(input_name="S1", output_name="S3")

    # Evaluate controllability -- do the inputs affect the outputs?
    _ = ctlsb.plotStaircaseResponse(initial_value=0, final_value=10, num_step=5) # Displays the outputs associated with the inputs

    # Identify the system by fitting to the input staircase
    _ = ctlsb.plotTransferFunctionFit()    # Uses the same staircase function as above

    # Construct a SISO closed loop system and do a trail and error design
    _ = ctlsb.plotClosedLoop(kP=1, kI=0.1, kF=0.5, setpoint=3)
    _ = ctlsb.plotClosedLoop(kP=0.8, kI=0.2, kF=0.1, setpoint=5)

    # Automatically design
    _ = ctlsb.plotSISOClosedLoopDesign(kP_spec=True, kI_spec=True, min_value=0, max_value=10, num_iteration=20, setpoint=4)
    param_dct = {kP: ctlsb.kP, kI: ctlsb.kI, kD: ctlsb.kD, kF: ctlsb.kF}

    # Explore variations on the design
    new_param_dct = {n: v*1.1 for n, v param_dct.items() for n in ["kP", "kI", "kD", "kF"]}
    _ = ctlsb.plotClosedLoop(setpoint=5, **new_param_dct)  # Do further trail and error


Conventions for API names. The following prefixes are used for methods:
  get: returns a value
  plot: plots a result; returns Timeseries of data plotted and AntimonyBuilder used to generate the data; may update state
      is_plot: bool (if False, not plot is done. Default is True.)
      Follow matplotlib conventions for plot options: xlabel, ylabel, title, xlim, ylim, markers, etc.
  set: changes internal state of ControlSBML. No value is returned. Little computation is done.
"""

from controlSBML.sbml_system import SBMLSystem
from controlSBML.antimony_builder import AntimonyBuilder
from controlSBML.siso_transfer_function_builder import SISOTransferFunctionBuilder
from controlSBML.siso_closed_loop_designer import SISOClosedLoopDesigner
from controlSBML.timeseries import Timeseries
from controlSBML.staircase import Staircase
from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML import util
import controlSBML.constants as cn
import controlSBML.msgs as msgs
from controlSBML.grid import Grid
from controlSBML.parallel_coordinates import ParallelCoordinates

from collections import namedtuple
import os
import control  # type: ignore
import numpy as np
from typing import List, Tuple, Optional

SETPOINT = 1
FIGSIZE = (5, 5)
# Optional dictionaries
TIMES_DCT = {
    cn.O_TIMES: cn.TIMES,
}
STAIRCASE_DCT = {
    cn.O_INITIAL_VALUE: 0,
    cn.O_FINAL_VALUE: 10,
    cn.O_NUM_STEP: 5,
}
SIMULATION_DCT = {
    cn.O_END_TIME: cn.TIMES[-1],
    cn.O_START_TIME: cn.TIMES[0],
    cn.O_NUM_POINT: len(cn.TIMES),
    cn.O_SELECTIONS: None,
}
CLOSED_LOOP_DCT = {
    cn.CP_KP: 0, 
    cn.CP_KI: 0, 
    cn.CP_KD: 0, 
    cn.CP_KF: 0, 
    cn.O_SETPOINT: 1, 
    cn.O_SIGN: -1,
    cn.O_NUM_PROCESS: -1,
    cn.O_NUM_RESTART: 1,
    cn.O_KP_SPEC: False,
    cn.O_KI_SPEC: False,
    cn.O_KD_SPEC: False,
    cn.O_KF_SPEC: False,
    }
SYSTEM_OPTION_DCT = {
    cn.O_DISTURBANCE_SPEC: cn.DisturbanceSpec(),
    cn.O_FITTER_METHOD: cn.DEFAULT_FITTER_METHOD,
    cn.O_IS_FIXED_INPUT_SPECIES: True,
    cn.O_IS_STEADY_STATE: False,
    cn.O_IS_GREEDY: False,
    cn.O_INPUT_NAME: None,
    cn.O_NOISE_SPEC: cn.NoiseSpec(),
    cn.O_OUTPUT_NAME: None,
    }
PLOT_OPTION_DCT = { 
    cn.O_AX: None,
    cn.O_AX2: None,
    cn.O_FIGURE: None,
    cn.O_FIGSIZE: FIGSIZE,
    cn.O_IS_PLOT: True,
    cn.O_LEGEND: True,
    cn.O_LEGEND_CRD: None,
    cn.O_LEGEND_SPEC: None,
    cn.O_MARKERS: False,
    cn.O_SUPTITLE: "",
    cn.O_TITLE: "",
    cn.O_XLABEL: "time",
    cn.O_XLABEL_ANGLE: 0,
    cn.O_XTICKLABELS: None,
    cn.O_XLIM: None,
    cn.O_YLABEL: "concentration",
    cn.O_YLIM: None,
    cn.O_YTICKLABELS: None,
    cn.O_WRITEFIG: False,
    }
OPTION_DCT = {**STAIRCASE_DCT, **TIMES_DCT, **CLOSED_LOOP_DCT, **SYSTEM_OPTION_DCT, **PLOT_OPTION_DCT,
               **SIMULATION_DCT}
# Option keywords
PLOT_KEYS = list(PLOT_OPTION_DCT.keys())
STAIRCASE_KEYS = list(STAIRCASE_DCT.keys())
TIMES_KEYS = list(TIMES_DCT.keys())
OPTION_KEYS = list(OPTION_DCT.keys())
SIMULATION_KEYS = list(SIMULATION_DCT.keys())

# Returned results
DesignResult = namedtuple("DesignResult", ["timeseries", "antimony_builder", "designs"])
GridDesignResult = namedtuple("GridDesignResult", ["timeseries", "antimony_builder", "designs"])
ModelResult = namedtuple("ModelResult", ["timeseries"])
StaircaseResponseResult = namedtuple("StaircaseResponseResult", ["timeseries", "antimony_builder"])
TransferFunctionFitResult = namedtuple("TransferFunctionFitResult", ["timeseries", "antimony_builder"]) 


class ControlSBML(object):

    def __init__(self, model_reference:str, 
                 roadrunner=None,
                 input_name:Optional[str]=None,
                 output_name:Optional[str]=None,
                 is_fixed_input_species:Optional[bool]=OPTION_DCT[cn.O_IS_FIXED_INPUT_SPECIES],
                 is_steady_state:Optional[bool]=OPTION_DCT[cn.O_IS_STEADY_STATE],
                 noise_spec = OPTION_DCT[cn.O_NOISE_SPEC],
                 disturbance_spec = OPTION_DCT[cn.O_DISTURBANCE_SPEC],
                 fitter_method:Optional[str]=OPTION_DCT[cn.O_FITTER_METHOD],
                 setpoint:Optional[float]=OPTION_DCT[cn.O_SETPOINT],
                 sign:Optional[int]=OPTION_DCT[cn.O_SIGN],
                 times:Optional[np.ndarray[float]]=None,
                 **kwargs):
        """
        model_reference: str
            string, SBML file or Roadrunner object
        roadrunner: ExtendedRoadrunner
        input_name: str
        output_name: str
        is_fixed_input_species: bool
        is_steady_state: bool
        noise_spec: NoiseSpec
        disturbance_spec: DisturbanceSpec
        kwargs: dict with options listed below. These are the default options used unless overridden.
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
                xlabel_angle: int (angle of x label)
                ylim: tupe-float (y lower limit, y upper limit)
                yticklabels: list-float (labels for y ticks)
            System options:
                input_name: str
                is_steady_state: bool (start system in steady state; default: False)
                is_fixed_input_species: bool (concentration of input species are controlled externally; default: False)
                output_name: str
                noise_spec: NoiseSpec (specification of noise added to measured output)
                          sine_amp, sine_freq, rand_mag, rand_std, dc_gain, slope
                disturbance_spec: NoiseSpec (specification of disturbance added to control input)
                            sine_amp, sine_freq, rand_mag, rand_std, dc_gain, slope
            Closed loop options:
                kP: float (proportional control)
                kI: float (integral control)
                kD: float (differential constant)
                kF: float (filter constant)
                setpoint: float (regulation point)
                sign: -1/+1 (direction of feedback: default: -1) 
            Staircase options:
                final_value: float (last value of input in staircase; default: maximum input value in SBML model)
                initial_value: float (first value of input in staircase; default: minimum input value in SBML model)
                num_step: int (number of steps in staircase; default: 5)
            Miscellaneous options
                times: list-float (times of simulation; default: np.linspace(0, 5, 50))
                disturbance: NoiseSpec (specification of disturbance added to control input)
                       sine_amp, sine_freq, rand_mag, rand_std, dc_gain, slope
                noise: NoiseSpec (specification of noise added to measured output)
                       sine_amp, sine_freq, rand_mag, rand_std, dc_gain, slope

        """
        self._checkKwargs(OPTION_KEYS, **kwargs)
        # Set initial values of all options
        for key, value in OPTION_DCT.items():  # type: ignore
            new_value = kwargs.get(key, value)
            setattr(self, key, new_value)
        # Initializations
        self.model_reference = model_reference
        if roadrunner is None:
            roadrunner = makeRoadrunner(model_reference)
        self._roadrunner = roadrunner
        # Other assignments
        self.input_name = input_name
        self.output_name = output_name
        self.is_fixed_input_species = is_fixed_input_species
        self.is_steady_state = is_steady_state
        self.noise_spec = noise_spec
        self.disturbance_spec = disturbance_spec
        self.fitter_method = OPTION_DCT[cn.O_FITTER_METHOD] if fitter_method is None else fitter_method
        self.setpoint = OPTION_DCT[cn.O_SETPOINT] if setpoint is None else setpoint
        self.sign = OPTION_DCT[cn.O_SIGN] if sign is None else sign
        self.times = OPTION_DCT[cn.O_TIMES] if times is None else times
        # Internal state
        self._sbml_system, self._transfer_function_builder =  self.setSystem(input_name=input_name,
                       output_name=output_name, 
                       is_fixed_input_species=is_fixed_input_species, is_steady_state=is_steady_state)  # type: ignore
        self._fitter_result = cn.FitterResult()
        self.antimony_builder:Optional[AntimonyBuilder] = None  # Last antimony builder used
        self.design_result:Optional[DesignResult] = None  # Result of last design
        
    def copy(self):
        ctlsb = ControlSBML(self.model_reference,
                roadrunner=self._roadrunner,
                input_name=self.input_name,
                output_name=self.output_name,
                is_fixed_input_species=self.is_fixed_input_species,
                is_steady_state=self.is_steady_state,
                fitter_method=self.fitter_method,
                setpoint=self.setpoint,
                sign=self.sign,
                times=self.times)
        for key in cn.CONTROL_PARAMETERS:
            setattr(ctlsb, key, getattr(self, key))
        return ctlsb
    
    def _checkKwargs(self, valids:Optional[List[str]]=OPTION_KEYS, invalids:Optional[List[str]]=None,
                     **kwargs):
        """
        Checks that the kwargs are valid.
        Args:
            valids: valid keys
        """
        keys = set(kwargs.keys())
        not_valids = keys.difference(set(valids))   # type: ignore
        if len(not_valids) > 0:
            raise ValueError("Invalid keyword arguments: %s" % str(not_valids))
        if invalids is None:
            return
        is_invalids = set(keys).intersection(set(invalids))  # type: ignore
        if len(is_invalids) > 0:
            raise ValueError("Keyword arguments invalid for this method: %s" % str(is_invalids))
    
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
    def getRoadrunner(self):
        return self._roadrunner

    def getAntimony(self):
        return self._roadrunner.getAntimony()
    
    def getGrid(self, min_value:int=0, max_value:int=10, num_coordinate:int=10, kP_spec:bool=False,
                kI_spec:bool=False, kD_spec:bool=False, kF_spec:bool=False, is_random=True)->Grid:
        """
        Creates a grid of values for the control parameters based on the specified control design.

        Args:
            min_value (int): _description_. Defaults to 0.
            max_value (int): _description_. Defaults to 10.
            num_coordinate (int): _description_. Defaults to 10.
            kP_spec (bool): Proportional control. Defaults to False.
            kI_spec (bool): Integral control. Defaults to False.
            kD_spec (bool): Differential control. Defaults to False.
            kF_spec (bool): Filter. Defaults to False.
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
        makeAxis(cn.CP_KP, kP_spec)
        makeAxis(cn.CP_KI, kI_spec)
        makeAxis(cn.CP_KD, kD_spec)
        makeAxis(cn.CP_KF, kF_spec)
        return grid

    def getPossibleInputs(self):
        return self._sbml_system.getValidInputs()
    
    def getPossibleOutputs(self):
        return self._sbml_system.getValidOutputs()
    
    def _getTimes(self, **kwargs)->np.ndarray:
        times = kwargs.get(cn.O_TIMES, self.times)
        if times is None:
            times = self.times
        return times
    
    def getOpenLoopTransferFunction(self):
        return self._fitter_result.transfer_function

    def _getStaircase(self, initial_value:Optional[float]=None,
                      final_value:Optional[float]=None, num_step:Optional[int]=None):
        dct = self.getOptions(keys=STAIRCASE_KEYS, 
                              initial_value=initial_value, final_value=final_value, num_step=num_step)
        initial_value = dct[cn.O_INITIAL_VALUE]
        final_value = dct[cn.O_FINAL_VALUE]
        num_step = dct[cn.O_NUM_STEP]
        return Staircase(initial_value=initial_value, final_value=final_value,
                         num_step=num_step, num_point=len(self.times))
    
    def getClosedLoopTransferFunction(self, sign:Optional[int]=None)->control.TransferFunction:
        """
        Calculates the closed loop transfer function

        Args:


        Returns: 
            control.TransferFunction
        """
        sign = self.sign if sign is None else sign
        designer = SISOClosedLoopDesigner(self._sbml_system, self.getOpenLoopTransferFunction(),
                sign=self.sign)
        kP = 0
        kI = 0
        kD = 0
        kF = 0
        for parameter_name in cn.CONTROL_PARAMETERS:
            parameter = getattr(self, parameter_name)
            if parameter is not None:
                if not util.isNumber(parameter):
                    raise ValueError("Must assign float to kP, kI, kD, and kF before using this method.")
                setattr(self, parameter_name, parameter)
        designer.set(kP=kP, kI=kI, kD=kD, kF=kF)
        return designer.closed_loop_transfer_function
    
    def _getParameterStr(self, parameters:List[str], **kwargs):
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

    def setSystem(self, input_name:Optional[str]=None, output_name:Optional[str]=None,
                  is_fixed_input_species:bool=True,
                  is_steady_state:bool=False)->Tuple[SBMLSystem, SISOTransferFunctionBuilder]:
        """
        Sets the options related to defining the system.

        Args:
            input_name: str
            output_name: str
            is_fixed_input_species: bool (if True, a species used as input is treated as fixed)
            is_fixed_steady_state: bool (if True, simulations are started in steady state)
        """
        self.input_name = input_name
        self.output_name = output_name
        self.is_fixed_input_species = is_fixed_input_species if is_fixed_input_species is not None else self.is_fixed_input_species
        self.is_steady_state = is_steady_state if is_steady_state is not None else self.is_steady_state
        sbml_system = SBMLSystem(self.model_reference, input_names=[self.input_name], # type: ignore
                            output_names=[self.output_name], # type: ignore
                            noise_spec=self.noise_spec, disturbance_spec=self.disturbance_spec,
                            is_fixed_input_species= self.is_fixed_input_species,
                            is_steady_state=self.is_steady_state, roadrunner=self._roadrunner)
        builder = SISOTransferFunctionBuilder(
            sbml_system,
            input_name=self.input_name,  # type: ignore
            output_name=self.output_name,  # type: ignore
            fitter_method=self.fitter_method)
        self._sbml_system = sbml_system
        self._transfer_function_builder = builder
        return sbml_system, builder
    
    def setOption(self, key:str, value:object):
        """
        Sets an option.

        Args:
            key: str
            value: object
        """
        setattr(self, key, value)

    def _initializeStaircaseOptions(self, initial_value=None, final_value=None, num_step=cn.DEFAULT_NUM_STEP):
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
            if self.input_name is not None:
                input_name = self.input_name
                is_assign_from_simulation = True
        # Calculate defaults if required
        if initial_value is None:
            if is_assign_from_simulation:
                result = self.plotModel(is_plot=False, selections=[input_name])
                initial_value = result.timeseries[input_name].min()
            else:
                initial_value = cn.DEFAULT_INITIAL_VALUE
        if final_value is None:
            if is_assign_from_simulation:
                final_value = result.timeseries[input_name].max()
            else:
                final_value = cn.DEFAULT_INITIAL_VALUE
        if num_step is None:
            num_step = cn.DEFAULT_NUM_STEP
        return initial_value, final_value, num_step
  

    ############ PLOTTERS ##############
    def plotModel(self, 
                  times:Optional[np.ndarray[float]]=None,
                  selections:Optional[List[str]]=OPTION_DCT[cn.O_SELECTIONS],
                  **kwargs)->ModelResult:
        """
        Plots the SBML model without modification.

        Args:
            times: numpy array (times of simulation)
            selections: list-str (selections for simulation)
            kwargs: dict (plot options)
        
        Returns:
            ModelResult
                timeseries: Timeseries
        """
        self._checkKwargs(valids=PLOT_KEYS, **kwargs)
        dct = self.getOptions(times=times, selections=selections, **kwargs)
        plot_dct = util.subsetDct(dct, PLOT_KEYS)
        times = dct[cn.O_TIMES]
        #
        self._roadrunner.resetAll()
        self._roadrunner.resetSelectionLists()
        if selections is not None:
            selections.extend([cn.TIME])    # type: ignore
            selections = list(set(selections))
        data = self._roadrunner.simulate(times[0], times[-1], len(times), selections=selections)  # type: ignore
        ts = Timeseries(data)
        is_plot = plot_dct[cn.O_IS_PLOT]
        if is_plot:
            util.plotOneTS(ts, **plot_dct)
        return ModelResult(timeseries=ts)
    
    def plotStaircaseResponse(self,
                              initial_value:Optional[float]=OPTION_DCT[cn.O_INITIAL_VALUE],
                              final_value:Optional[float]=OPTION_DCT[cn.O_FINAL_VALUE],
                              num_step:Optional[int]=OPTION_DCT[cn.O_NUM_STEP],
                              times:Optional[np.ndarray]=None,
                              **kwargs)->StaircaseResponseResult:
        """
        Simulates the staircase response and plots it.

        Args:
            initial_value: float (initial value of the staircase; default is the minimum value of the input if it's a species)
            final_value: float (final value of the staircase)
            num_step: int (number of steps in the staircase)
            times: numpy array (times of simulation)
            kwargs: dict (plot options)

        Returns:
            StaircaseResponseResult
                timeseries: Timeseries
                antimony_builder: AntimonyBuilder
        """
        self._checkKwargs(valids=PLOT_KEYS, **kwargs)
        dct = self.getOptions(initial_value=initial_value, final_value=final_value,
                              num_step=num_step, times=times, **kwargs)
        times = dct[cn.O_TIMES]
        initial_value = dct[cn.O_INITIAL_VALUE]
        final_value = dct[cn.O_FINAL_VALUE]
        num_step = dct[cn.O_NUM_STEP]
        plot_dct = util.subsetDct(dct, PLOT_KEYS)
        staircase = self._getStaircase(initial_value=initial_value, final_value=final_value, num_step=num_step)
        response_ts, builder = self._transfer_function_builder.makeStaircaseResponse(times=times,
            staircase=staircase, is_steady_state=self.is_steady_state,
            )
        self.antimony_builder = builder
        self._transfer_function_builder.plotStaircaseResponse(response_ts, **plot_dct)
        return StaircaseResponseResult(timeseries=response_ts, antimony_builder=builder)
    
    def plotTransferFunctionFit(self,
            num_zero:int=cn.DEFAULT_NUM_ZERO,
            num_pole:int=cn.DEFAULT_NUM_POLE, 
            fit_start_time: Optional[float]=None,
            fit_end_time:Optional[float]=None,
            initial_value:Optional[float]=None,
            final_value:Optional[float]=None,
            num_step:Optional[int]=OPTION_DCT[cn.O_NUM_STEP],
            fitter_method:Optional[str]=OPTION_DCT[cn.O_FITTER_METHOD],
            times:Optional[np.ndarray]=None,
            **kwargs)->TransferFunctionFitResult:
        """
        Simulates the staircase response and plots it. Sets the fitter result.

        Args:
            num_zero: int (number of zeros)
            num_pole: int (number of poles)
            fit_start_time: float (time at which fitting starts)
            fit_end_time: float (time at which fitting ends)
            initial_value: float (initial value of staircase)
            final_value: float (final value of staircase)
            num_step: int (number of steps in staircase)
            times: numpy array (times of simulation)
            fitter:method: str
                'poly' uses a polynomial fit
                'gpz' fit gain, then poles, then zeros
            times: numpy array (times of simulation)
            kwargs: (plot options, noise, disturbance)

        Returns:
            TransferFunctionFitResult
                timeseries: Timeseries (predicted, staircase input, simulated output)
                antimony_builder: AntimonyBuilder
        """
        # Check the options
        self._checkKwargs(valids=PLOT_KEYS, **kwargs)
        # Setup values
        dct = self.getOptions(
            num_zero=num_zero,
            num_pole=num_pole,
            fit_start_time=fit_start_time,
            fit_end_time=fit_end_time,
            initial_value=initial_value,
            final_value=final_value,
            num_step=num_step,
            fitter_method=fitter_method,
            times=times,
            **kwargs)
        plot_dct = util.subsetDct(dct, PLOT_KEYS)
        times = dct[cn.O_TIMES]
        fitter_method = dct[cn.O_FITTER_METHOD]
        #
        if (self.input_name is None) or (self.output_name is None):
            raise ValueError("Must specify the input and output species to use this method.")
        if fit_start_time is None:
            fit_start_time = times[0]  # type: ignore
        if fit_end_time is None:
            fit_end_time = times[-1]   # type: ignore
        # Construct the staircase
        staircase = self._getStaircase(initial_value=initial_value, final_value=final_value, num_step=num_step)
        # Do the fit
        self._fitter_result = self._transfer_function_builder.plotTransferFunctionFit(
            num_zero=num_zero,
            num_pole=num_pole,
            staircase=staircase,
            fit_start_time=fit_start_time,
            fit_end_time=fit_end_time,
            fitter_method=fitter_method,
            times=times,
            **plot_dct)
        return TransferFunctionFitResult(
                timeseries=self._fitter_result.time_series, 
                antimony_builder=self._fitter_result.antimony_builder,
        )
    
    def _plotClosedLoop(self, 
                        kP:Optional[float]=OPTION_DCT[cn.CP_KP],
                        kI:Optional[float]=OPTION_DCT[cn.CP_KI],
                        kD:Optional[float]=OPTION_DCT[cn.CP_KD],
                        kF:Optional[float]=OPTION_DCT[cn.CP_KF],
                        setpoint:Optional[float]=OPTION_DCT[cn.O_SETPOINT],
                        sign:Optional[float]=OPTION_DCT[cn.O_SIGN],
                        selections:Optional[List[str]]=OPTION_DCT[cn.O_SELECTIONS],
                        times:Optional[np.ndarray]=None,
                        **kwargs)->Tuple[Timeseries, AntimonyBuilder]:
        """
        Plots the closed loop response. Control parameters not explicity specified are None.

        Args:
            kP: proportional control parameter
            kI: integral control parameter
            kD: differential control parameter
            kF: filter parameter
            sign: int (direction of feedback: -1 or 1)
            setpoint: float (regulation point)
            times: numpy array (times of simulation)
            selections: list-str (selections for simulation)
            kwargs: (plot options, noise, disturbance)
        Returns:
            Timeseries
            AntimonyBuilder
        """
        self._checkKwargs(PLOT_KEYS, **kwargs)
        #
        # Construct the SBML system
        dct = self.getOptions(sign=sign, setpoint=setpoint, kP=kP, kI=kI, kD=kD, kF=kF,
                              times=times, selections=selections, **kwargs)
        sign = dct[cn.O_SIGN]
        setpoint = dct[cn.O_SETPOINT]
        kP = dct[cn.CP_KP]
        kI = dct[cn.CP_KI]
        kD = dct[cn.CP_KD]
        kF = dct[cn.CP_KF]
        times = dct[cn.O_TIMES]
        plot_dct = util.subsetDct(dct, PLOT_KEYS)
        # Plot the response
        response_ts, builder = self._sbml_system.simulateSISOClosedLoop(input_name=self.input_name,
                output_name=self.output_name, sign=sign,
                kP=kP, kI=kI, kD=kD, kF=kF, setpoint=setpoint, selections=selections,
                times=times,
                )
        if response_ts is None:
            msg = "System is unstable for kP={}, kI={}, kD={}, kF={}".format(kP, kI, kD, kF)
            msgs.warn(msg)
            return response_ts, builder
        if plot_dct[cn.O_IS_PLOT]:
            self._sbml_system.plotSISOClosedLoop(response_ts, setpoint, times=times, **plot_dct)
        return response_ts, builder

    def plotGridDesign(self, 
                       grid:Grid,
                       setpoint:Optional[float]=OPTION_DCT[cn.O_SETPOINT],
                       sign:Optional[float]=OPTION_DCT[cn.O_SIGN],
                       times:Optional[np.ndarray]=None,
                       num_process:Optional[int]=OPTION_DCT[cn.O_NUM_PROCESS],
                       num_restart:Optional[int]=OPTION_DCT[cn.O_NUM_RESTART],
                       selections:Optional[List[str]]=OPTION_DCT[cn.O_SELECTIONS],
                       **kwargs)->GridDesignResult:
        """
        Plots the results of a closed loop design based a grid of values for the control parameters.
        Persists the closed loop design (kP, kI, kD, kF) if a design is found.

        Args:
            grid: Grid (grid of values for the control parameters)
            setpoint: float (regulation point)
            sign: float (direction of feedback: -1 or 1)
            times: numpy array (times of simulation)
            num_process: int (number of processes to use; -1 means use all available)
            selections: list-str (selections for the simulation)
            kwargs: dict (plot options)
        Returns:
            GridDesignResult
                timeseries: Timeseries
                antimony_builder: AntimonyBuilder
        """
        new_keys = list(PLOT_KEYS)
        self._checkKwargs(new_keys, **kwargs)
        option_dct = self.getOptions(sign=sign, setpoint=setpoint, times=times, selections=selections,
                                     num_process=num_process, num_restart=num_restart, **kwargs)
        # Initialize parameters
        setpoint = option_dct[cn.O_SETPOINT]
        sign = option_dct[cn.O_SIGN]
        times = option_dct[cn.O_TIMES]
        selections = option_dct[cn.O_SELECTIONS]
        num_process = option_dct[cn.O_NUM_PROCESS]
        num_restart = option_dct[cn.O_NUM_RESTART]
        plot_dct = util.subsetDct(option_dct, PLOT_KEYS)
        # Process the request
        designer = SISOClosedLoopDesigner(self._sbml_system, self.getOpenLoopTransferFunction(),
                is_steady_state=self.is_steady_state,
                times=times,
                sign=sign,
                setpoint=setpoint,
                input_name=self.input_name,
                output_name=self.output_name)
        # Translate axis names
        self.design_result = designer.designAlongGrid(grid,
                                 num_process=num_process,  # type: ignore
                                 num_restart=num_restart)    # type: ignore
        if designer.residual_mse is None:
            msgs.warn("No design found!")
            return GridDesignResult(timeseries=None, antimony_builder=None,
                                    designs=self.design_result)
        # Persist the design parameters
        self.setOption(cn.CP_KP, designer.kP)
        self.setOption(cn.CP_KI, designer.kI)
        self.setOption(cn.CP_KD, designer.kD)
        self.setOption(cn.CP_KF, designer.kF)
        # Plot the results 
        plot_dct[cn.O_YLABEL] = self.output_name if not cn.O_YLABEL in kwargs else kwargs[cn.O_YLABEL]
        if (not cn.O_TITLE in plot_dct) or (len(plot_dct[cn.O_TITLE]) == 0):
            plot_dct[cn.O_TITLE] = self._getParameterStr([cn.CP_KP, cn.CP_KI, cn.CP_KF, cn.CP_KD],
                                                     kP=getattr(self, cn.CP_KP),
                                                     kI=getattr(self, cn.CP_KI),
                                                     kD=getattr(self, cn.CP_KD),
                                                     kF=getattr(self, cn.CP_KF),
                                                     )
        response_ts, antimony_builder = self._plotClosedLoop(
                times=times,
                setpoint=setpoint,
                sign=self.sign,
                kP=designer.kP,
                kI=designer.kI,
                kD=designer.kD,
                kF=designer.kF,
                selections=selections,
                **plot_dct)
        return GridDesignResult(timeseries=response_ts, antimony_builder=antimony_builder,
                                designs=self.design_result)
    
    def plotAllDesignResults(self, design_result:Optional[DesignResult]=None,
                         title=None, figsize=None, columns=None, 
                         round_digit:int=4, is_plot=True):
        """
        Does a parallel coordinate plot of the design results.

        Args:
            design_result: DesignResult
            title: str
            figsize: tuple
            columns: list-str (order in whch columns appear on plot)
            round_digit: int (number of digits to round to)
        """
        if design_result is None:
            design_result = self.design_result
        df = design_result.dataframe.copy() # type: ignore
        del df[cn.REASON]
        ParallelCoordinates.plotParallelCoordinates(df,  # type: ignore
                                title=title, figsize=figsize,
                                columns=columns,
                                value_column=cn.SCORE,
                                round_digit=round_digit,
                                is_plot=is_plot)

    @staticmethod
    def setSpec(val):
        if isinstance(val, bool):
            return val
        if isinstance(val, int):
            return float(val)
        else:
            return val

    def plotDesign(self, 
                kP_spec:bool=OPTION_DCT[cn.O_KP_SPEC],
                kI_spec:bool=OPTION_DCT[cn.O_KI_SPEC],
                kD_spec:bool=OPTION_DCT[cn.O_KD_SPEC],
                kF_spec:bool=OPTION_DCT[cn.O_KF_SPEC],
                setpoint:Optional[float]=OPTION_DCT[cn.O_SETPOINT],
                sign:Optional[float]=OPTION_DCT[cn.O_SIGN],
                times:Optional[np.ndarray]=None,
                num_process:Optional[int]=OPTION_DCT[cn.O_NUM_PROCESS],
                num_restart:Optional[int]=OPTION_DCT[cn.O_NUM_RESTART],
                selections:Optional[List[str]]=OPTION_DCT[cn.O_SELECTIONS],
                min_parameter_value:float=0,
                max_parameter_value:float=10,
                num_coordinate:int=3,
                is_report:bool=False, 
                **kwargs)->DesignResult:
        """
        Plots the results of a closed loop design. The design is specified by the parameters
        kP_spec, kI_spec, kD_spec, and kF_spec.
           None or False: do not include the parameter
           True: include the parameter and find a value
           float: use this value for the parameter.
        Persists the closed loop design (kP, kI, kD, kF) if a design is found.

        Args:
            kP_spec: float (specification of proportional gain)
            kI_spec: float (specification of integral gain)
            kD_spec: float (specification of differential gain)
            kF_spec: float (specification of filter gain)
            setpoint: float (setpoint for the control)
            sign: float (direction of feedback: -1 or 1)
            min_parameter_value: float (minimum value for kP, kI, kD, kF; may be a dictionary)
            max_parameter_value: float (maximum value for kP, kI, kD, kF; may be a dictionary)
            num_coordinate: int (number of coordinate descent iterations)
            is_report: bool (report progress on the design search)
            times: numpy array (times of simulation)
            num_process: int (number of processes to use; -1 means use all available)
            selections: list-str (selections for the simulation)
            kwargs: dict (plot options)
        Returns:
            DesignResult
                timeseries: Timeseries of the response
                antimony_builder: AntimonyBuilder simulated with the design parameters
                designs
        """
        #
        def setValue(val):
            if util.isNumber(val):
                return val
            return 0.0
        def setSpec(key, val):
            if val is None:
                return self.__getattribute__(key)
            return val
        #
        kP_spec=self.setSpec(kP_spec)
        kI_spec=self.setSpec(kI_spec)
        kD_spec=self.setSpec(kD_spec)
        kF_spec=self.setSpec(kF_spec)
        new_keys = list(PLOT_KEYS)
        self._checkKwargs(new_keys, **kwargs)
        option_dct = self.getOptions(kP_spec=kP_spec,
                                     kI_spec=kI_spec,
                                     kD_spec=kD_spec,
                                     kF_spec=kF_spec,
                                     setpoint=setpoint,
                                     sign=sign,
                                     is_report=is_report,
                                     num_process=num_process,
                                     selections=selections,
                                     times=times,
                                     **kwargs)
        kP_spec = option_dct[cn.O_KP_SPEC]
        kI_spec = option_dct[cn.O_KI_SPEC]
        kD_spec = option_dct[cn.O_KD_SPEC]
        kF_spec = option_dct[cn.O_KF_SPEC]
        setpoint = option_dct[cn.O_SETPOINT]
        sign = option_dct[cn.O_SIGN]
        num_process = option_dct[cn.O_NUM_PROCESS]
        selections = option_dct[cn.O_SELECTIONS]
        times = option_dct[cn.O_TIMES]
        #
        is_greedy = option_dct[cn.O_IS_GREEDY]
        plot_dct = util.subsetDct(option_dct, PLOT_KEYS)
        #
        designer = SISOClosedLoopDesigner(self._sbml_system, self.getOpenLoopTransferFunction(),
                is_steady_state=self.is_steady_state,
                times=times,
                sign=sign,
                setpoint=setpoint,
                input_name=self.input_name,
                output_name=self.output_name)
        self.design_result = designer.design(
            kP_spec=kP_spec,
            kI_spec=kI_spec,
            kD_spec=kD_spec,
            kF_spec=kF_spec,
            num_restart=num_restart,
            min_value=min_parameter_value,
            max_value=max_parameter_value,
            num_coordinate=num_coordinate,
            is_greedy=is_greedy,
            is_report=is_report,
            num_process=num_process)   # type: ignore
        if designer.residual_mse is None:
            msgs.warn("No design found!")
            return DesignResult(timeseries=None, antimony_builder=None,
                                designs=self.design_result)
        # Persist the design parameters
        self.setOption(cn.CP_KP, designer.kP)
        self.setOption(cn.CP_KI, designer.kI)
        self.setOption(cn.CP_KD, designer.kD)
        self.setOption(cn.CP_KF, designer.kF)
        # Plot the results
        title = "" if not cn.O_TITLE in kwargs else kwargs[cn.O_TITLE]
        if len(title) == 0:
            if not kP_spec == False:
                title += "kP=%f " % setValue(designer.kP)
            if not kI_spec == False:
                title += "kI=%f " % setValue(designer.kI)
            if not kD_spec == False:
                title += "kD=%f " % setValue(designer.kD)
            if not kF_spec == False:
                title += "kF=%f " % setValue(designer.kF)
        plot_dct[cn.O_TITLE] = title
        response_ts, antimony_builder = self._plotClosedLoop(
                times=times,
                setpoint=setpoint,
                sign=sign,  # type: ignore
                kP=designer.kP,
                kI=designer.kI,
                kD=designer.kD,
                kF=designer.kF,
                selections=selections,
                **plot_dct)
        return DesignResult(timeseries=response_ts, antimony_builder=antimony_builder,
                            designs=self.design_result)
    
    def _plotDesignResult(self, **kwargs):
        """
        Plots the mean squared error of all points searched in the last design.

        Args:
            kwargs: dict (plot options)
            AntimonyBuilder
        """
        self._checkKwargs(**kwargs)
        designer = SISOClosedLoopDesigner(self._sbml_system, self.getOpenLoopTransferFunction())
        if designer.design_result_df is None:
            msgs.warn("No design found!")
            return None, None
        return designer.plotDesignResult(**kwargs)
    
    ############## MISC ################
    def getOptions(self, keys=None, **kwargs)->dict:
        """
        Optons the options in the object.
            STAIRCASE_OPTIONS: initial_value, final_value, num_step
            TIMES_OPTIONS: times
            CLOSED_LOOP_PARAMETERS: kP, kI, kD, kF, setpoint, sign
            PLOT_OPTIONS: ax ax2 end_time figure figsize is_plot
                            suptitle title writefig xlabel xlim xticklabels ylabel ylim yticklabels 
        Args:
            keys: list-str (keys in the dictionary returned)
            kwargs: dict (plot options)
        Returns:
            dict
        """
        keys = OPTION_KEYS if keys is None else keys
        new_kwargs = {}
        # FIXME: Doesn't allow user to change option to None
        for key in keys:
            new_kwargs[key] = getattr(self, key) if not key in kwargs.keys() else kwargs[key]
            new_kwargs[key] = getattr(self, key) if new_kwargs[key] is None else new_kwargs[key]
        return new_kwargs
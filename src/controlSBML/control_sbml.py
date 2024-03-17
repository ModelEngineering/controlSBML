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
    param_dct = {kP: ctlsb.kP, kI: ctlsb.kI, kF: ctlsb.kF}

    # Explore variations on the design
    new_param_dct = {n: v*1.1 for n, v param_dct.items() for n in ["kP", "kI", "kF"]}
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
from controlSBML.option_set import OptionSet
from controlSBML.grid import Grid, Point

import os
import control  # type: ignore
import numpy as np
from typing import List, Dict, Tuple, Optional

PLOT_KWARGS = list(set(cn.PLOT_KWARGS).union(cn.FIG_KWARGS))
SETPOINT = 1
STAIRCASE_OPTIONS = [cn.O_INITIAL_VALUE, cn.O_FINAL_VALUE, cn.O_NUM_STEP]
TIMES_OPTIONS = [cn.O_TIMES]
CLOSED_LOOP_PARAMETERS = [cn.CP_KP, cn.CP_KI, cn.CP_KF, cn.O_SETPOINT, cn.O_SIGN]
SYSTEM_SPECIFICATIONS = [cn.O_INPUT_NAME, cn.O_OUTPUT_NAME, cn.O_IS_FIXED_INPUT_SPECIES, cn.O_IS_STEADY_STATE,
                        cn.FITTER_METHOD]
PLOT_OPTIONS = list(cn.PLOT_KWARGS)
PLOT_OPTIONS.extend(cn.FIG_KWARGS)
PLOT_OPTIONS.append(cn.O_MARKERS)
OPTIONS = STAIRCASE_OPTIONS + TIMES_OPTIONS + CLOSED_LOOP_PARAMETERS + PLOT_OPTIONS + SYSTEM_SPECIFICATIONS
CONTROL_PARAMETERS = [cn.CP_KP, cn.CP_KI, cn.CP_KF]
FIGSIZE = (5, 5)
INITIAL_PLOT_OPTION_DCT = {cn.O_TITLE: "", cn.O_SUPTITLE: "", cn.O_WRITEFIG: False,
                        cn.O_XLABEL: "time", cn.O_YLABEL: "concentration",
                        cn.O_FIGSIZE: FIGSIZE,
                        cn.O_IS_PLOT: True,
                        cn.O_LEGEND_SPEC: None,
                        cn.O_LEGEND_CRD: None,
                        cn.O_YLIM: None,
                        cn.O_XLIM: None,
                        cn.O_XTICKLABELS: None,
                        cn.O_YTICKLABELS: None,
                        cn.O_AX: None,
                        cn.O_AX2: None,
                        cn.O_FIGURE: None,
                        cn.O_MARKERS: False,
                        }
SAVE_PATH = os.path.join(cn.DATA_DIR, "control_sbml.csv")


class ControlSBML(object):

    def __init__(self, model_reference:str, 
                 roadrunner=None,
                 input_name:Optional[str]=None,
                 output_name:Optional[str]=None,
                 is_fixed_input_species:Optional[bool]=False,
                 is_steady_state:Optional[bool]=False,
                 fitter_method:Optional[str]=cn.DEFAULT_FITTER_METHOD,
                 setpoint:Optional[float]=SETPOINT,
                 sign:Optional[int]=-1,
                 times:Optional[np.ndarray[float]]=None,
                 save_path:Optional[str]=None,
                 **kwargs):
        """
        model_reference: str
            string, SBML file or Roadrunner object
        roadrunner: ExtendedRoadrunner
        input_name: str
        output_name: str
        is_fixed_input_species: bool
        is_steady_state: bool
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
            input_name: str
            is_steady_state: bool (start system in steady state; default: False)
            is_fixed_input_species: bool (concentration of input species are controlled externally; default: False)
            output_name: str
        Closed loop options:
            kF: float (filter constant)
            kI: float (integral control)
            kP: float (proportional control)
            setpoint: float (regulation point)
            sign: -1/+1 (direction of feedback: default: -1) 
        Staircase options:
            final_value: float (last value of input in staircase; default: maximum input value in SBML model)
            initial_value: float (first value of input in staircase; default: minimum input value in SBML model)
            num_step: int (number of steps in staircase; default: 5)
        Miscellaneous options
            times: list-float (times of simulation; default: np.linspace(0, 5, 50))

        """
        self._checkKwargs(PLOT_OPTIONS, **kwargs)
        # Initializations
        self.model_reference = model_reference
        if roadrunner is None:
            roadrunner = makeRoadrunner(model_reference)
        self._roadrunner = roadrunner
        self.setpoint = setpoint
        self.times = cn.TIMES if times is None else times
        self.sign = sign
        self.fitter_method = fitter_method
        # Input and output names
        self.input_name = input_name
        self.output_name = output_name
        self.is_fixed_input_species = is_fixed_input_species
        self.is_steady_state = is_steady_state
        self.save_path = save_path
        # Internal state
        self._sbml_system, self._transfer_function_builder =  self.setSystem(input_name=input_name, output_name=output_name,
                       is_fixed_input_species=is_fixed_input_species, is_steady_state=is_steady_state)  # type: ignore
        self._fitter_result = cn.FitterResult()
        # Options
        for key, value in INITIAL_PLOT_OPTION_DCT.items():
            setattr(self, key, value)
        # Set initial values of options
        self.ax = None
        self.initial_value, self.final_value, self.num_step = self._initializeStaircaseOptions()
        self.figure = kwargs.get(cn.O_FIGURE, None)
        self.figsize = kwargs.get(cn.O_FIGSIZE, FIGSIZE)
        self.is_greedy = kwargs.get(cn.O_IS_GREEDY, False)
        self.is_plot = kwargs.get(cn.O_IS_PLOT, True)
        self.is_fixed_input_species = kwargs.get(cn.O_IS_FIXED_INPUT_SPECIES, True)    
        self.is_steady_state = kwargs.get(cn.O_IS_STEADY_STATE, False)
        self.markers = kwargs.get(cn.O_MARKERS, False)
        self.sign = kwargs.get(cn.O_SIGN, -1)
        self.setpoint = kwargs.get(cn.O_SETPOINT, SETPOINT)
        self.suptitle = kwargs.get(cn.O_SUPTITLE, "")
        self.title = kwargs.get(cn.O_TITLE, "")
        self.writefig = kwargs.get(cn.O_WRITEFIG, False)
        # Outputs
        self.kP = None
        self.kI = None
        self.kF = None
        
    def copy(self):
        options = list(PLOT_OPTIONS)
        options.extend(["kP", "kI", "kF"])
        kwargs = {}
        for key in options:
            kwargs[key] = getattr(self, key)
        return ControlSBML(self.model_reference,
                roadrunner=self._roadrunner,
                input_name=self.input_name,
                output_name=self.output_name,
                is_fixed_input_species=self.is_fixed_input_species,
                is_steady_state=self.is_steady_state,
                fitter_method=self.fitter_method,
                setpoint=self.setpoint,
                sign=self.sign,
                times=self.times,
                save_path=self.save_path,
                    **kwargs)
    
    def _checkKwargs(self, valids:Optional[List[str]]=OPTIONS, invalids:Optional[List[str]]=None,
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
        if not self._getOptions() == other._getOptions():
            return False
        if not self._fitter_result.equals(other._fitter_result):
            return False
        return True

    ############ GETTERS ##############
    def getAntimony(self):
        return self._roadrunner.getAntimony()
    
    def getGrid(self, min_value:int=0, max_value:int=10, num_coordinate:int=10, kP_spec:bool=False,
                kI_spec:bool=False, kF_spec:bool=False, is_random=True)->Grid:
        """
        Creates a grid of values for the control parameters based on the specified control design.

        Args:
            min_value (int): _description_. Defaults to 0.
            max_value (int): _description_. Defaults to 10.
            num_coordinate (int): _description_. Defaults to 10.
            kP_spec (bool): Proportional control. Defaults to False.
            kI_spec (bool): Integral control. Defaults to False.
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
        makeAxis(cn.CP_KF, kF_spec)
        return grid

    def getPossibleInputs(self):
        return self._sbml_system.getValidInputs()
    
    def getPossibleOutputs(self):
        return self._sbml_system.getValidOutputs()
    
    def _getOptions(self, options:Optional[dict]=None):
        """
        Gets current values of the options

        Returns: dict. Keys are listed below by category.
            STAIRCASE_OPTIONS: initial_value, final_value, num_step
            TIMES_OPTIONS: times
            CLOSED_LOOP_PARAMETERS: kP, kI, kF, setpoint, sign
            PLOT_OPTIONS: ax ax2 end_time figure figsize is_plot
                            suptitle title writefig xlabel xlim xticklabels ylabel ylim yticklabels 
        """
        if options is None:
            option_lst = OPTIONS
        return {k: getattr(self, k) for k in option_lst}
    
    def _getTimes(self, **kwargs)->np.ndarray:
        times = kwargs.get(cn.O_TIMES, self.times)
        if times is None:
            times = self.times
        return times
    
    def getOpenLoopTransferFunction(self):
        return self._fitter_result.transfer_function

    def _getStaircase(self, initial_value:Optional[float]=None,
                      final_value:Optional[float]=None, num_step:Optional[int]=None):
        initial_value = self.initial_value if initial_value is None else initial_value
        final_value = self.final_value if final_value is None else final_value
        num_step = self.num_step if num_step is None else num_step
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
        kF = 0
        for parameter_name in CONTROL_PARAMETERS:
            parameter = getattr(self, parameter_name)
            if parameter is not None:
                if not util.isNumber(parameter):
                    raise ValueError("Must assign float to kP, kI, and kF before using this method.")
                setattr(self, parameter_name, parameter)
        designer.set(kP=kP, kI=kI, kF=kF)
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
    def _setOptions(self, **kwargs):
        """
        Sets values of options.

        Args:
            kwargs: dict of options
        """
        self.setOptionSet(**kwargs)

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
                ts = self.plotModel(is_plot=False)
                initial_value = ts[input_name].min()
            else:
                initial_value = cn.DEFAULT_INITIAL_VALUE
        if final_value is None:
            if is_assign_from_simulation:
                final_value = ts[input_name].max()
            else:
                final_value = cn.DEFAULT_INITIAL_VALUE
        if num_step is None:
            num_step = cn.DEFAULT_NUM_STEP
        return initial_value, final_value, num_step

    ############ PLOTTERS ##############
    def plotModel(self, 
                  times:Optional[np.ndarray[float]]=None,
                  **kwargs)->Timeseries:
        """
        Plots the SBML model without modification.

        Args:
            times: numpy array (times of simulation)
            kwargs: dict (plot options)
        
        Returns:
            Timeseries
        """
        options = list(PLOT_OPTIONS)
        options.append(cn.O_IS_PLOT)
        options.extend(TIMES_OPTIONS)
        self._checkKwargs(options, **kwargs)
        is_plot = kwargs.get(cn.O_IS_PLOT, True)
        times = kwargs.get(cn.O_TIMES, self.times)
        self._roadrunner.reset()
        if (self.input_name is None) and (not "input_name" in kwargs.keys()):
            selections = None
        else:
            selections = [cn.TIME, self.input_name, self.output_name]
        data = self._roadrunner.simulate(times[0], times[-1], len(times), selections=selections)  # type: ignore
        ts = Timeseries(data)
        new_kwargs = dict(kwargs)
        if cn.O_TIMES in new_kwargs:
            new_kwargs.pop(cn.O_TIMES)
        if is_plot:
            util.plotOneTS(ts, markers=self.markers, **kwargs)
        return ts
    
    def plotStaircaseResponse(self,
                              initial_value:Optional[float]=None,
                              final_value:Optional[float]=None,
                              num_step:Optional[int]=cn.DEFAULT_NUM_STEP,
                              times:Optional[np.ndarray]=None,
                              **kwargs)->Tuple[Timeseries, AntimonyBuilder]:
        """
        Simulates the staircase response and plots it.

        Args:
            initial_value: float (initial value of the staircase; default is the minimum value of the input if it's a species)
            final_value: float (final value of the staircase)
            num_step: int (number of steps in the staircase)
            times: numpy array (times of simulation)
            kwargs: dict (plot options)

        Returns:
            Timeseries
                index: time (ms)
                columns: <output_name>, staircase
            AntimonyBuilder
        """
        self._checkKwargs(**kwargs, valids=PLOT_OPTIONS)
        times = self._getTimes(times=times)
        initial_value, final_value, num_step = self._initializeStaircaseOptions(initial_value=initial_value,
            final_value=final_value, num_step=num_step)
        staircase = self._getStaircase(initial_value=initial_value, final_value=final_value, num_step=num_step)
        response_ts, builder = self._transfer_function_builder.makeStaircaseResponse(times=times,
            staircase=staircase, is_steady_state=self.is_steady_state,
            )
        self._transfer_function_builder.plotStaircaseResponse(response_ts, **kwargs)
        return response_ts, builder
    
    def plotTransferFunctionFit(self,
            num_zero:int=cn.DEFAULT_NUM_ZERO,
            num_pole:int=cn.DEFAULT_NUM_POLE, 
            fit_start_time: Optional[float]=None,
            fit_end_time:Optional[float]=None,
            initial_value:Optional[float]=None,
            final_value:Optional[float]=None,
            num_step:Optional[int]=cn.DEFAULT_NUM_STEP,
            fitter_method:Optional[str]=cn.DEFAULT_FITTER_METHOD,
            times:Optional[np.ndarray]=None,
            **kwargs):
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
            fitter:method: str (method for fitting: 'poly' or 'gpz')
            times: numpy array (times of simulation)
            kwargs: (plot options)

        Returns:
            Timeseries (predicted, staircase)
            AntimonyBuilder
        """
        # Check the options
        self._checkKwargs(valids=PLOT_OPTIONS, **kwargs)
        # Setup values
        times = self.times if times is None else times
        fitter_method = self.fitter_method if fitter_method is None else fitter_method
        #
        if (self.input_name is None) or (self.output_name is None):
            raise ValueError("Must specify the input and output species to use this method.")
        if fit_start_time is None:
            fit_start_time = times[0]
        if fit_end_time is None:
            fit_end_time = times[-1]
        # Construct the staircase
        initial_value, final_value, num_step = self._initializeStaircaseOptions(initial_value=initial_value,
                                                                                final_value=final_value,
                                                                                num_step=num_step)
        staircase = self._getStaircase(initial_value=initial_value, final_value=final_value, num_step=num_step)
        new_kwargs = dict(kwargs)
        # Do the fit
        new_kwargs.update({cn.FITTER_METHOD: fitter_method,
                           cn.O_TIMES: times,})
        self._fitter_result = self._transfer_function_builder.plotTransferFunctionFit(num_zero=num_zero,
                            num_pole=num_pole, staircase=staircase,
                            fit_start_time=fit_start_time, fit_end_time=fit_end_time,
                            **new_kwargs)
        return self._fitter_result.time_series, self._fitter_result.antimony_builder
    
    def _plotClosedLoop(self, 
                        kP:Optional[float]=None, 
                        kI:Optional[float]=None, 
                        kF:Optional[float]=None, 
                        sign:Optional[int]=None,
                        setpoint:Optional[float]=None,
                        times:Optional[np.ndarray]=None,
                        **kwargs):
        """
        Plots the closed loop response. Control parameters not explicity specified are None.

        Args:
            kP: proportional control parameter
            kI: integral control parameter
            kF: filter parameter
            sign: int (direction of feedback: -1 or 1)
            setpoint: float (regulation point)
            times: numpy array (times of simulation)
            kwargs: plot options
        Returns:
            Timeseries
            AntimonyBuilder
        """
        self._checkKwargs(PLOT_OPTIONS, **kwargs)
        #
        # Construct the SBML system
        sign = self.sign if sign is None else sign
        kP = self.kP if kP is None else kP
        kI = self.kI if kI is None else kI
        kF = self.kF if kF is None else kF
        setpoint = self.setpoint if setpoint is None else setpoint
        response_ts, builder = self._sbml_system.simulateSISOClosedLoop(input_name=self.input_name,
                output_name=self.output_name, sign=sign,
                kP=kP, kI=kI, kF=kF, setpoint=setpoint,
                times=times,
                )
        if (not cn.O_TITLE in kwargs) or (len(kwargs[cn.O_TITLE]) == 0):
            kwargs["title"] = self._getParameterStr([cn.CP_KP, cn.CP_KI, cn.CP_KF], kP=kP, kI=kI, kF=kF)
        _ = kwargs.setdefault(cn.O_MARKERS, False)
        self._sbml_system.plotSISOClosedLoop(response_ts, setpoint, times=times, **kwargs)
        return response_ts, builder

    def plotGridDesign(self, 
                       grid:Grid,
                       setpoint:Optional[float]=None,
                       sign:Optional[float]=None,
                       times:Optional[np.ndarray]=None,
                       num_process:Optional[int]=-1,
                       num_restart:Optional[int]=3,
                       **kwargs):
        """
        Plots the results of a closed loop design based a grid of values for the control parameters.
        Persists the closed loop design (kP, kI, kF) if a design is found.

        Args:
            grid: Grid (grid of values for the control parameters)
            setpoint: float (regulation point)
            sign: float (direction of feedback: -1 or 1)
            times: numpy array (times of simulation)
            num_process: int (number of processes to use; -1 means use all available)
            kwargs: dict (plot options)
        Returns:
            Timeseries
            AntimonyBuilder
        """
        save_path = None   # Disable "save_path" feature
        #
        self._checkKwargs(PLOT_OPTIONS, **kwargs)
        # Initialize parameters
        setpoint = self.setpoint if setpoint is None else setpoint
        sign = self.sign if sign is None else sign
        times = self._getTimes(**kwargs)
        if save_path is not None:
            if len(save_path) == 0:
                save_path = self.save_path
            if os.path.isfile(save_path):
                os.remove(save_path)
        # Process the request
        designer = SISOClosedLoopDesigner(self._sbml_system, self.getOpenLoopTransferFunction(),
                is_steady_state=self.is_steady_state,
                times=times,
                sign=sign,
                setpoint=setpoint,
                input_name=self.input_name,
                output_name=self.output_name,
                save_path=self.save_path)
        # Translate axis names
        designer.designAlongGrid(grid, is_greedy=self.is_greedy, num_process=num_process,  # type: ignore
                                 num_restart=num_restart)    # type: ignore
        if designer.residual_mse is None:
            msgs.warn("No design found!")
            return None, None
        # Persist the design parameters
        self.kP = designer.kP
        self.kI = designer.kI
        self.kF = designer.kF
        options = dict(kwargs)
        options[cn.O_YLABEL] = self.output_name if not cn.O_YLABEL in kwargs else kwargs[cn.O_YLABEL]
        response_ts, antimony_builder = self._plotClosedLoop(
                times=times,
                setpoint=setpoint,
                sign=self.sign,
                kP=self.kP,
                kI=self.kI,
                kF=self.kF,
                **options)
        return response_ts, antimony_builder

    def plotDesign(self, 
                   kP_spec:bool=False, 
                   kI_spec:bool=False,
                   kF_spec:bool=False,
                   setpoint:Optional[float]=None,
                   sign:Optional[float]=None,
                   min_parameter_value:float=0,
                   max_parameter_value:float=10,
                   num_restart:int=3, 
                   num_coordinate:int=3,
                   is_report:bool=False, 
                   num_process:int=-1,
                   times:Optional[np.ndarray]=None, 
                   **kwargs)->Tuple[Timeseries, AntimonyBuilder]:
        """
        Plots the results of a closed loop design. The design is specified by the parameters kP_spec, kI_spec, and kF_spec.
           None or False: do not include the parameter
           True: include the parameter and find a value
           float: use this value for the parameter.
        Persists the closed loop design (kP, kI, kF) if a design is found.

        Args:
            kP_spec: float (specification of proportional gain)
            kI_spec: float (specification of integral gain)
            kF_spec: float (specification of filter gain)
            setpoint: float (setpoint for the control)
            sign: float (direction of feedback: -1 or 1)
            min_parameter_value: float (minimum value for kP, kI, kF; may be a dictionary)
            max_parameter_value: float (maximum value for kP, kI, kF; may be a dictionary)
            num_coordinate: int (number of coordinate descent iterations)
            is_report: bool (report progress on the design search)
            times: numpy array (times of simulation)
            kwargs: dict (plot options)
        Returns:
            Timeseries
            AntimonyBuilder
        """
        def setValue(val):
            if util.isNumber(val):
                return val
            return 0.0
        #
        times = self.times if times is None else times
        setpoint = self.setpoint if setpoint is None else setpoint
        sign = self.sign if sign is None else sign
        self._checkKwargs(PLOT_OPTIONS, **kwargs)
        #
        designer = SISOClosedLoopDesigner(self._sbml_system, self.getOpenLoopTransferFunction(),
                is_steady_state=self.is_steady_state,
                times=times,
                sign=sign,
                setpoint=setpoint,
                input_name=self.input_name,
                output_name=self.output_name,
                save_path=self.save_path)
        designer.design(kP_spec=kP_spec, kI_spec=kI_spec, kF_spec=kF_spec,
                num_restart=num_restart, min_value=min_parameter_value, max_value=max_parameter_value,
            num_coordinate=num_coordinate, is_greedy=self.is_greedy, is_report=is_report, num_process=num_process)
        if designer.residual_mse is None:
            msgs.warn("No design found!")
            return None, None  # type: ignore
        # Persist the design parameters
        self.kP = designer.kP
        self.kI = designer.kI
        self.kF = designer.kF
        # Plot the results
        title = "" if not cn.O_TITLE in kwargs else kwargs[cn.O_TITLE]
        if len(title) == 0:
            title = "kP=%f, kI=%f, kF=%f" % (setValue(designer.kP), setValue(designer.kI), setValue(designer.kF))
        new_kwargs = dict(kwargs)
        new_kwargs[cn.O_TITLE] = title
        response_ts, antimony_builder = self._plotClosedLoop(
                times=times,
                setpoint=setpoint,
                sign=sign,  # type: ignore
                kP=self.kP,
                kI=self.kI,
                kF=self.kF,
                **new_kwargs)
        return response_ts, antimony_builder
    
    def _plotDesignResult(self, save_path:Optional[str]=None, **kwargs):
        """
        Plots the mean squared error of all points searched in the last design.

        Args:
            save_path: str (path to CSV file where design results are saved)
            kwargs: dict (plot options)
            AntimonyBuilder
        """
        valids = ["save_path"]
        valids = list(set(valids).union(PLOT_OPTIONS))
        self._checkKwargs(**kwargs)
        if save_path is None:
            save_path = self.save_path
        designer = SISOClosedLoopDesigner(self._sbml_system, self.getOpenLoopTransferFunction(), save_path=save_path)
        if designer.design_result_df is None:
            msgs.warn("No design found!")
            return None, None
        return designer.plotDesignResult(**kwargs)
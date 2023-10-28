"""
Top level API

Usage:

    # Examine an SBML model
    ctlsb = ControlSBML(model)
    ctlsb.setTimes(np.linspace(0, 100, 1000))
    _ = ctlsb.plotModel()

    # Define the system
    ctlsb.getPossibleInputs()
    ctlsb.getPossibleOutputs()
    ctlsb.setSystem(input_names=["S1"], output_names=["S3"])
    ctlsb.getInputNames(), ctlsb.getOutputNames()

    # Evaluate controllability -- do the inputs affect the outputs?
    ctlsb.setStaircase(input_value=0, final_value=10, num_step=5)
    ctlsb.plotStaircaseResponse() # Displays the outputs associated with the inputs

    # Identify the system by fitting to the input staircase
    _ = ctlsb.plotTransferFunctionFit()
    open_loop_transfer_function = ctlsb.getOpenLoopTransferFunction()

    # Construct a SISO closed loop system and do a trail and error design
    _ = ctlsb.plotSISOClosedLoopSystem(kp=1, ki=0.1, kf=0.5, setpoint=3)
    _ = ctlsb.plotSISOClosedLoopSystem(kp=0.8, ki=0.2, kf=0.1, setpoint=5)

    # Automatically design
    _ = ctlsb.plotSISOClosedLoopDesign(kp=True, ki=True, min_value=0, max_value=10, num_iteration=20, setpoint=4)
    param_dct = ctlsb.getSISOClosedLoopParameters()

    # Explore deviations on the design
    new_param_dct = {n: v*1.1 for n, v param_dct.items()}
    _ = ctlsb.plotSISOClosedLoopSystem(setpoint=5, **new_param_dct)  # Do further trail and error


Conventions for API names. The following prefixes are used for methods:
  get: returns a value
  plot: plots a result; often returns an object. Can be significant computation. Keyword arguments are plot options.
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
from controlSBML.option_management.options import Options
import controlSBML.constants as cn

import numpy as np
from docstring_expander.expander import Expander

PLOT_KWARGS = list(set(cn.PLOT_KWARGS).union(cn.FIG_KWARGS))


class _StaircaseParameter(object):
    # Repository of parameters for staircase

    def __init__(self, initial_value=cn.DEFAULT_INITIAL_VALUE,
                 final_value=cn.DEFAULT_FINAL_VALUE, num_step=cn.DEFAULT_NUM_STEP):
        """
        initial_value: float
        final_value: float
        num_step: int
        """
        self.initial_value = initial_value
        self.final_value = final_value
        self.num_step = num_step

    def copy(self):
        return _StaircaseParameter(initial_value=self.initial_value,
                                   final_value=self.final_value,
                                   num_step=self.num_step)


class _SISODClosedLoopParameter(object):
    def __init__(self, kp=None, ki=None, kf=None, setpoint=1, sign=-1):
        """
        kp: float
        ki: float
        kf: float
        setpoint: float
        sign: int (-1: negative feedback; +1: positive feedback)
        """
        self.kp = kp
        self.ki = ki
        self.kf = kf
        self.setpoint = setpoint
        self.sign = sign

    def copy(self):
        return _SISODClosedLoopParameter(kp=self.kp, ki=self.ki, kf=self.kf, setpoint=self.setpoint, sign=self.sign)
    
    def get(self):
        return {"kp": self.kp, "ki": self.ki, "kf": self.kf, "setpoint": self.setpoint, "sign": self.sign}


class ControlSBML(object):

    def __init__(self, model_reference):
        """
        model_reference: str
            string, SBML file or Roadrunner object
        is_fixed_input_species: bool (input species are fixed)
        is_steady_state: bool (start the simulation at steady state)
        """
        # First initializations
        self.model_reference = model_reference
        self.roadrunner = makeRoadrunner(self.model_reference)
        self.setTimes()   # self.times
        self._sbml_system = None
        self.staircase_parameter = _StaircaseParameter()
        self.siso_transfer_function_builder = None
        self._fitter_result = cn.FitterResult()
        self.siso_closed_loop_parameter = _SISODClosedLoopParameter()

    def copy(self):
        ctlsb = ControlSBML(self.model_reference)
        ctlsb.setTimes(self.getTimes())
        if self._sbml_system is not None:
            ctlsb._sbml_system = self._sbml_system.copy()
        ctlsb.staircase_parameter = self.staircase_parameter.copy()
        ctlsb.siso_transfer_function_builder = self.siso_transfer_function_builder.copy()
        ctlsb._fitter_result = self._fitter_result.copy()
        ctlsb.siso_closed_loop_parameter = self.siso_closed_loop_parameter.copy()
        return ctlsb

    ############ GETTERS ##############
    def getSBMLSystem(self):
        if self._sbml_system is None:
            raise ValueError("Must no SBMLSystem. Use 'makeSystem'.")
        return self._sbml_system
        
    def getClosedLoopParameters(self):
        return self.siso_closed_loop_parameter.get()

    def getPossibleInputs(self):
        return self.getSBMLSystem().getValidInputs()
    
    def getPossibleOutputs(self):
        return self.getSBMLSystem().getValidOutputs()

    def getInputName(self):
        return self.getSBMLSystem().input_names[0]

    def getOutputName(self):
        return self.getSBMLSystem().output_names[0]
    
    def getOpenLoopTransferFunction(self):
        if self._fitter_result.transfer_function is None:
            raise ValueError("Must use 'plotTransferFunctionFit' before 'getOpenLoopTransferFunction.")
        return self._fitter_result.transfer_function
    
    def getFitterResult(self):
        return self._fitter_result
    
    def getSISOClosedLoopParameters(self):
        """
        Provides the parameters of the closed loop system
        """
        return self.getSBMLSystem().get()

    def getTimes(self):
        return self._times

    def getStaircase(self):
        return self._staircase
    
    def getClosedLoopTransferFunction(self):
        """
        Provides the closed loop transfer function
        """
        sbml_system = self.getSBMLSystem()
        parameters = self.siso_closed_loop_parameter
        if parameters is None:
            raise ValueError("Must use 'plotSISOClosedLoopSystem' before 'getClosedLoopTransferFunction'.") 
        designer = SISOClosedLoopDesigner(sbml_system, self.getOpenLoopTransferFunction(),
                setpoint=parameters.setpoint, is_steady_state=sbml_system.is_steady_state,
                sign=-parameters.sign)
        designer.set(kp=parameters.kp, ki=parameters.ki, kf=parameters.kf)
        return designer.closed_loop_tf

    ############ SETTERS ##############
    def setSystem(self, input_names=None, output_names=None, is_fixed_input_species=True, is_steady_state=False):
        """
        Creates an SBML system for the specified inputs and outputs.

        Args:
            input_names: list-str
            output_names: list-str
            is_fixed_input_species: bool (if True, a species used as input is treated as fixed)
            is_fixed_steady_state: bool (if True, simulations are started in steady state)
        """
        self._sbml_system = SBMLSystem(self.roadrunner, input_names=input_names, output_names=output_names,
                                      is_fixed_input_species=is_fixed_input_species, is_steady_state=is_steady_state)
    
    def setTimes(self, times=None):
        if times is None:
            times = np.linspace(cn.START_TIME, cn.END_TIME, cn.POINTS_PER_TIME*(cn.END_TIME-cn.START_TIME))
        self._times = times
    
    def setStaircase(self, initial_value=None, final_value=None, num_step=cn.DEFAULT_NUM_STEP):
        """
        Sets the staircase parameters.

        Args:
            initial_value: float (initial value of the staircase; default is the minimum value of the input if it's a species)
            final_value: float (final value of the staircase)
            num_step: int (number of steps in the staircase)
        """
        if (initial_value is None) or (final_value is None):
            if self._sbml_system is None:
                ts = self.plotModel(is_plot=False)
        # Determine how defaults are assigned
        is_assign_from_simulation = False
        if (initial_value is None) or (final_value is None):
            if self._sbml_system is not None:
                input_name = self.getInputName()
                if input_name in self.getSBMLSystem().antimony_builder.species_names:
                    is_assign_from_simulation = True
        # Calculate defaults if required
        if initial_value is None:
            if is_assign_from_simulation:
                initial_value = ts[input_name].min()
            else:
                initial_value = cn.DEFAULT_INITIAL_VALUE
        if final_value is None:
            if is_assign_from_simulation:
                final_value = ts[input_name].max()
            else:
                initial_value = cn.DEFAULT_FINAL_VALUE
        self.staircase_parameter = _StaircaseParameter(initial_value=initial_value, final_value=final_value,
                                                         num_step=num_step)

    ############ PLOTTERS ##############
    def plotModel(self, is_plot=True, **kwargs):
        """
        Plots the original model.

        Args:
            kwargs: plot options
        
        Returns:
            Timeseries
        """
        _, plot_fig_opts = self._getOptions(kwargs)
        data = self.roadrunner.simulate(self._times[0], self._times[-1], len(self._times))
        ts = Timeseries(data)
        if is_plot:
            util.plotOneTS(ts, **plot_fig_opts)
        return ts
    
    @Expander(cn.KWARGS, PLOT_KWARGS) 
    def plotStaircaseResponse(self, **kwargs):
        """
        Simulates the staircase response and plots it.

        Args:
            kwargs: plot options
            #@expand

        Returns:
            ctl.Timeseries
                index: time (ms)
                columns: <output_name>, staircase
            AntimonyBuilder
        """
        sim_opts, plot_fig_opts = self._getOptions(kwargs)
        sbml_system = self.getSBMLSystem()
        response_ts, builder = self.siso_transfer_function_builder.makeStaircaseResponse(
            staircase=self._staircase, is_steady_state=sbml_system.is_steady_state, **sim_opts)
        self.siso_transfer_function_builder.plotStaircaseResponse(response_ts, **plot_fig_opts)
        return response_ts, builder
    
    @Expander(cn.KWARGS, cn.ALL_KWARGS) 
    def plotTransferFunctionFit(self, num_numerator=cn.DEFAULT_NUM_NUMERATOR,
                            num_denominator=cn.DEFAULT_NUM_DENOMINATOR, fit_start_time=None, fit_end_time=None, **kwargs):
        """
        Simulates the staircase response and plots it. Sets the fitter result.

        Args:
            num_numerator: int (number of numerator terms)
            num_denominator: int (number of denominator terms)
            fit_start_time: float (time at which fitting starts)
            fit_end_time: float (time at which fitting ends)
            kwargs: dict (simulation options as described below)
            #@expand

        Returns:
            Timeseries (predicted, staircase)
            AntimonyBuilder
        """
        sim_opts, plot_fig_opts = self._getOptions(kwargs)
        self._fitter_result = self.siso_transfer_function_builder.fitTransferFunction(num_numerator=num_numerator,
                            num_denominator=num_denominator, staircase=self._staircase,
                            fit_start_time=fit_start_time, fit_end_time=fit_end_time, **sim_opts)
        self.siso_transfer_function_builder.plotFitterResult(self._fitter_result, **plot_fig_opts)
        return self._fitter_result.time_series, self._fitter_result.antimony_builder
    
    def plotSISOClosedLoopSystem(self, setpoint, sign=-1, kp=None, ki=None, kf=None, **kwargs):
        """
        Args:
            setpoint: float (setpoint for regulation)
            sign: int (-1: negative feedback; +1: positive feedback)
            kp: float (proportional gain)
            ki: float (integral gain)
            kf: float (filter gain)
            kwargs: dict (plot options)
            #@expand
        Returns:
            Timeseries
            AntimonyBuilder
        """
        _, plot_fig_opts = self._getOptions(kwargs)
        sbml_system = self.getSBMLSystem()
        self.siso_closed_loop_parameter = _SISODClosedLoopParameter(setpoint=setpoint, kp=kp, ki=ki, kf=kf, sign=sign)
        response_ts, builder = sbml_system.simulateSISOClosedLoop(input_name=self.getInputName(), output_name=self.getOutputNmae(),
                kp=kp, ki=ki, kf=kf, setpoint=setpoint,
                times=self._times, is_steady_state=sbml_system.is_steady_state, inplace=False, initial_input_value=None,
                sign=-1)
        sbml_system.plotSISOClosedLoop(response_ts, setpoint, markers=None, **plot_fig_opts)
        return response_ts, builder

    
    def plotSISOClosedLoopDesign(self, setpoint, sign=-1, kp=None, ki=None, kf=None, min_parameter_value=0,
                                 max_parameter_value=10, max_iteration=20, num_restart=3, **kwargs):
        """
        Args:
            setpoint: float (setpoint for regulation)
            sign: int (-1: negative feedback; +1: positive feedback)
            kp: float (proportional gain)
            ki: float (integral gain)
            kf: float (filter gain)
            min_parameter_value: float (minimum value for kp, ki, kf; may be a dictionary)
            max_parameter_value: float (maximum value for kp, ki, kf; may be a dictionary)
            kwargs: dict (plot options)
            #@expand
        Returns:
            Timeseries
            AntimonyBuilder
        """
        sbml_system = self.getSBMLSystem()
        designer = SISOClosedLoopDesigner(sbml_system, self.getOpenLoopTransferFunction(),
                setpoint=setpoint, is_steady_state=sbml_system.is_steady_state,
                sign=-sign)
        designer.design(self, input_name=self.getInputName(), output_name=self.getOutputName(),
                kp=kp, ki=ki, kf=kf, max_iteration=max_iteration,
                num_restart=num_restart, min_value=min_parameter_value, max_value=max_parameter_value)
        self.siso_closed_loop_parameter = _SISODClosedLoopParameter(setpoint=setpoint,
                kp=designer.kp, ki=designer.ki, kf=designer.kf, sign=sign)
        response_ts, antimony_builder = self.plotSISOClosedLoopSystem(
                self.siso_closed_loop_parameter.setpoint, 
                sign=self.siso_closed_loop_parameter.sign,
                kp=self.siso_closed_loop_parameter.kp, ki=None, kf=None, **kwargs)
        return response_ts, antimony_builder

   
    ############ PROPERTIES AND PRIVATE ##############
    @property
    def _staircase(self):
        """
        Creates a staircase to use for system evaluations.
        """
        return Staircase(initial_value=self.staircase_parameter.initial_value, final_value=self.staircase_parameter.final_value,
                         num_step=self.staircase_parameter.num_step, num_point=len(self._times))

    @staticmethod 
    def _getOptions(options):
        options = Options(**options)
        sim_opts, plot_opts, fig_opts = options.parse()
        plot_fig_opts = dict(plot_opts)
        plot_fig_opts.update(fig_opts)
        return sim_opts, plot_fig_opts
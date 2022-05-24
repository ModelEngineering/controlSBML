"""Builds and Evaluates SISO Systems"""

"""
TODO
1. makeNoiseAndDisturbance
2. addConnections
3. makeSISOFullStateClosedLoop
4. Consider a class hierarchy
   SISOCLosedLoop [common codes]
   SISOPIDClosedLoop, SISOFullStateClosedLoop
       make method constructs the system
"""

"""
The SISO closed loop system has consists of a set of elements that are combined
with a SISO system to form a closed loop system. The input to this system
is the reference signal with the name "cl_input.ref". The output
is the signal "cl_output.out".

A SISOClosedLoop may have the following elements:
  * ControlSBML SISO System
  * Controller
  * Filter (with DC gain = 1)
  * Noise
  * Disturbance
In addition to these elements, there may be summation nodes that add output
signals that are input to other elements.
Methods indicate which of these elements are present

Evaluations provided are:
  * Find dcgain of ControlSBML system with many inputs and one output
  * DC gain of the closed loop system
  * Evaluate the closed loop system for disturbances and noise

# TODO: Give better names to attachmentpoints
Below are the nodes and their attachment points. Connections are
made between attachment points.
  1 CL_INPUT: input to the network. A Passthru.
    * ref (in): input to the node (the reference signal)
    * out (out): output from the node
  2 SYSTEM: system under control. Inputs and outputs are species names.
  3 SUM_F_R: adds the reference signal to the filter output
    * ref (in): reference signal
    * fltr (in): filter signal
    * sum (out): output from the node
  4 CONTROLLER: has many variations. Some attachment points are:
    * err (in): error input
    * <species name> (in)
    * u (out): control output sent to the system input
  5 FLTR: filter
    * in (in): closed loop output
    * out (out): output of filter
  6 DISTURBANCE
    * out (out): disturbance generated
  7 SUM_D_U: Summation of controller output and disturbance
    * u (in): controller output 
    * d (in): disturbance
    * sum (out): summation
  8 NOISE
    * out (out): noise generated
  9 SUM_N_Y: Summation of system output and noise
    * n (in): noise
    * y (in): system output
    * sum (out): closed loop output
  3 SUM_D_U: adds the reference signal to the filter output
    * U (in): controller signal
    * D (in): disturbance
    * sum (out): output from the node
 10 CL_OUTPUT: output from closed loop system
    * in (in): measured output plus noise (if present)
    * out (out)

"""

import controlSBML as ctl
from controlSBML import constants as cn
from controlSBML.simulate_system import makeStateVector
from controlSBML import msgs
from controlSBML import util

import collections
import control
import numpy as np
import pandas as pd


CL_INPUT = "cl_input"
CL_INPUT_REF = "cl_input.ref"
CL_INPUT_OUT = "cl_input.out"
CL_OUTPUT = "cl_output"
CL_OUTPUT_IN = "cl_output.in"
CL_OUTPUT_OUT = "cl_output.out"
CONTROLLER = "controller"
CONTROLLER_IN = "controller.in"
CONTROLLER_REF = "controller.ref"
CONTROLLER_OUT = "controller.out"
DISTURBANCE = "disturbance"
DISTURBANCE_OUT = "disturbance.out"
FLTR = "fltr"
FLTR_IN = "fltr.in"
FLTR_OUT = "fltr.out"
NOISE = "noise"
NOISE_OUT = "noise.out"
SUM_D_U = "sum_D_U"
SUM_D_U_U = "sum_D_U.u"
SUM_D_U_D = "sum_D_U.d"
SUM_D_U_OUT = "sum_D_U.out"
SUM_F_R = "sum_F_R"
SUM_F_R_REF = "sum_F_R.ref"
SUM_F_R_FLTR = "sum_F_R.fltr"
SUM_F_R_OUT = "sum_F_R.out"
SUM_N_Y = "sum_N_Y"
SUM_N_Y_Y = "sum_N_Y.y"
SUM_N_Y_N = "sum_N_Y.n"
SUM_N_Y_OUT = "sum_N_Y.out"
SYSTEM = "system"


# sys: list-NonlinearIOSystem
# con: connections - list-str
# inp:inplist - list-str
# out:  outlist - list-str
ConnectionAssembly = collections.namedtuple("ConnectionAssembly",
      "sys con inp out")


class SISOClosedLoopSystem(object):
    # Constructs and evaluates a SISO closed loop system.

    def __init__(self, ctlsb):
        self.ctlsb = ctlsb
        self.factory = ctl.IOSystemFactory()
        self.closed_loop_system = None
        # Elements in the closed loop system 
        ## Subsysems
        self.system = None
        self.controller = None
        self.fltr = None
        ## External signals
        self.cl_input = None
        self.cl_output = None
        self.noise = None
        self.disturbance = None
        ## Summations
        self.sum_N_Y = None
        self.sum_F_R = None
        self.sum_D_U = None
    
    def report(self):
        """
        Report log from running the closed loop system.
        """
        return self.factory.report()

    def evaluateControllability(self, times, input_names=None,
           output_names=None):
        """
        Evaluates the controllability of the inputs on the outputs.
        If no input (output) is specified, then all inputs (outputs)
        are considered for self.ctlsb.

        Parameters
        ----------
        times: list-float
             times at which dc gain is evaluated
        input_names: list-str
        output_names: list-str
        
        Returns
        -------
        dict
            key: time
            value: pd.DataFrame
                 index: input_names
                 columns: output_names
                 value: dcgain of output for input
        """
        if input_names is None:
            input_names = self.ctlsb.input_names
        if output_names is None:
            output_names = self.ctlsb.output_names
        dct = {t: {} for t in times}
        for output_name in output_names:
            for input_name in input_names:
                ctlsb = ctl.ControlSBML(self.ctlsb.model_reference,
                      input_names=[input_name],
                      output_names=[output_name],
                      is_reduced=self.ctlsb.is_reduced)
                for time in times:
                    dct[time][output_name] = []
                for time in times:
                    tf = ctlsb.makeTransferFunction(time)
                    dct[time][output_name].append(tf.dcgain())
        # Construct the DataFrames
        result_dct = {}
        for time in dct.keys():
            result_dct[time] = pd.DataFrame(dct[time], index=input_names)
        #
        return result_dct

    def makePIDClosedLoopSystem(self,
          kp=1, ki=0, kd=0,                       # Controller parameters
          disturbance_amp=0, disturbance_frq=0,   # Disturbance
          noise_amp=0, noise_frq=0,               # Noise
          kf=None,                                # Filter
          system_input=None, system_name=None,    # SISO input, output
          system_output=None,
          closed_loop_outputs=None):              # list of outputs from closed loop system
        """
        Creates a closed loop system for a ControlSBML object.
        The closed loop system
        includes a PID controller and a filter.

        Parameters
        ----------
        kp: float
            Proportional constant for controller
        ki: float
            Integral constant for controller
        kd: float
            Differential constant for controller
        noise_amp: float
            amplitude of the sine wave for noise
        noise_frq: float
            frequency of the sine wave for noise
        noise_amp: float
            amplitude of the sine wave for noise
        noise_frq: float
            frequency of the sine wave for noise
        kf: float
            Constant for filter. If None, no filter.
        system_input: str
            name of the input used in the SISO system
            the name must be present in the ControlSBML object
        system_output: str
            name of the output used in the SISO system
            the name must be present in the ControlSBML object
        closed_loop_ouputs: list-str
            If None, CL_INPUT_REF, CL_OUTPUT_OUT
        
        Returns
        -------
        Timeseries
           CL_INPUT_REF: reference input to System
           CL_OUTPUT_OUT: output from the System
           noise.out: noise input to system
           disturbance.out: noise input to system
        """
        assembly = self._makeDisturbanceNoiseCLinputoutput(
              disturbance_amp=disturbance_amp,
              disturbance_frq=disturbance_frq,
              noise_amp=noise_amp,
              noise_frq=noise_frq,
              )
        sys_lst = list(assembly.sys)
        connections = list(assembly.con)
        if closed_loop_outputs is None:
            closed_loop_outputs = assembly.out
        # Make controller, system, filter
        self.sum_F_R = self.factory.makeAdder(SUM_F_R, input_names=["ref", "fltr"])
        # Choose first system input and last system output as defaults
        if system_input is None:
            system_input = self.ctlsb.input_names[0]
        if system_output is None:
            system_output = self.ctlsb.output_names[-1]
        ctlsb = ctl.ControlSBML(self.ctlsb.model_reference,
              input_names=[system_input], output_names=[system_output])
        self.system = ctlsb.makeNonlinearIOSystem(SYSTEM)
        system_in = "%s.%s" % (SYSTEM, system_input)
        system_out = "%s.%s" % (SYSTEM, system_output)
        #
        self.controller = self.factory.makePIDController(CONTROLLER,
              kp=kp, ki=ki, kd=kd)
        #
        if kf is None:
            self.fltr = self.factory.makePassthru(FLTR)
        else:
            self.fltr = self.factory.makeFilter(FLTR, kf)
        #
        sys_lst.extend([self.system, self.controller, self.fltr, self.sum_F_R])
        # Additional names
        filter_str = "-%s" % FLTR_OUT
        system_inp = "%s.%s" % (SYSTEM, system_input)
        system_out = "%s.%s" % (SYSTEM, system_output)
        connections.append([SUM_F_R_REF, CL_INPUT_OUT])
        connections.append([CONTROLLER_IN, SUM_F_R_OUT])
        connections.append([SUM_D_U_U, CONTROLLER_OUT])
        connections.append([system_in, SUM_D_U_OUT])
        connections.append([SUM_N_Y_Y, system_out])
        connections.append([FLTR_IN, CL_OUTPUT_OUT])
        connections.append([SUM_F_R_FLTR, filter_str])
        # Construct the interconnected system
        self.closed_loop_system = control.interconnect(sys_lst,
              connections=connections,
              inplist=[CL_INPUT_REF],
              outlist=closed_loop_outputs
            )

    def makeFullStateClosedLoopSystem(self,
          time=0,                                 # Time of linearization
          poles=-1,                               # Desired poles or dominant pole
          disturbance_amp=0, disturbance_frq=0,   # Disturbance
          noise_amp=0, noise_frq=0,               # Noise
          system_input=None, system_name=None,    # SISO input, output
          system_output=None,
          closed_loop_outputs=None):              # list of outputs from closed loop system
        """
        Creates a closed loop system for a ControlSBML object. The closed loop system
        includes a PID controller and a filter.

        Parameters
        ----------
        poles: float/list-flaot
            Desired poles (or just dominant pole)
        noise_amp: float
            amplitude of the sine wave for noise
        noise_frq: float
            frequency of the sine wave for noise
        noise_amp: float
            amplitude of the sine wave for noise
        noise_frq: float
            frequency of the sine wave for noise
        kf: float
            Constant for filter. If None, no filter.
        system_input: str
            name of the input used in the SISO system
            the name must be present in the ControlSBML object
        system_output: str
            name of the output used in the SISO system
            the name must be present in the ControlSBML object
        closed_loop_ouputs: list-str
            If None, CL_INPUT_REF, CL_OUTPUT_OUT
        
        Returns
        -------
        control.IOSystem.InterconnectedSystem
           CL_INPUT_IN: reference input to System
           CL_OUTPUT_OUT: output from the System
           noise.out: noise input to system
           disturbance.out: noise input to system
        """
        assembly = self._makeDisturbanceNoiseCLinputoutput(
              disturbance_amp=disturbance_amp,
              disturbance_frq=disturbance_frq,
              noise_amp=noise_amp,
              noise_frq=noise_frq,
              )
        sys_lst = list(assembly.sys)
        connections = list(assembly.con)
        if closed_loop_outputs is None:
            closed_loop_outputs = list(assembly.out)
        # Choose first system input and last system output as defaults
        if system_input is None:
            system_input = self.ctlsb.input_names[0]
        if system_output is None:
            system_output = self.ctlsb.output_names[-1]
        initial_ctlsb = ctl.ControlSBML(self.ctlsb.model_reference,
              input_names=[system_input], output_names=[system_output])
        output_names = list(initial_ctlsb.state_names)
        output_names.remove(system_input)
        ctlsb = ctl.ControlSBML(self.ctlsb.model_reference,
              input_names=[system_input], output_names=output_names)
        self.system = ctlsb.makeNonlinearIOSystem(SYSTEM)
        system_in = "%s.%s" % (SYSTEM, system_input)
        system_out = "%s.%s" % (SYSTEM, system_output)
        #
        self.controller = self.factory.makeFullStateController(CONTROLLER,
              ctlsb, poles=poles, time=time)
        #
        sys_lst.extend([self.system, self.controller])
        # Connections between system and controller
        for name in output_names:
            controller_input = "%s.%s" % (CONTROLLER, name)
            system_output = "%s.%s" % (SYSTEM, name)
            connections.append([controller_input, system_output])
        # Additional names
        connections.append([SUM_D_U_U, CONTROLLER_OUT])
        connections.append([system_in, SUM_D_U_OUT])
        connections.append([SUM_N_Y_Y, system_out])
        connections.append([CONTROLLER_REF, CL_INPUT_OUT])
        # Construct the interconnected system
        self.closed_loop_system = control.interconnect(sys_lst,
              connections=connections,
              inplist=[CL_INPUT_REF],
              outlist=closed_loop_outputs
            )

    def _makeDisturbanceNoiseCLinputoutput(self,
          disturbance_amp=0, disturbance_frq=0,   # Disturbance
          noise_amp=0, noise_frq=0):              # Noise
        """
        Creates the noise and disturbance and related elements.
        Constructs all internal connections.
            cl_input
            disturbance
            noise
            SUM_N_Y
            SUM_D_U
            cl_output
        connections
            DISTURBANCE_OUT -> SUM_D_U_D
            NOISE_OUT -> SUM_N_Y_N
            SUM_N_Y_OUT -> CL_OUTPUT_IN

        Parameters
        ----------
        noise_amp: float
            amplitude of the sine wave for noise
        noise_frq: float
            frequency of the sine wave for noise
        noise_amp: float
            amplitude of the sine wave for noise
        noise_frq: float
            frequency of the sine wave for noise
        
        Returns
        -------
        ConnectionAssembly
        """
        # Create the elements of the feedback loop
        self.cl_input = self.factory.makePassthru(CL_INPUT, input_name="ref")
        self.noise = self.factory.makeSinusoid(NOISE, noise_amp, noise_frq)
        self.disturbance = self.factory.makeSinusoid(DISTURBANCE,
              disturbance_amp, disturbance_frq)
        self.sum_N_Y = self.factory.makeAdder(SUM_N_Y, input_names=["y", "n"])
        self.sum_D_U = self.factory.makeAdder(SUM_D_U, input_names=["u", "d"])
        self.cl_output = self.factory.makePassthru(CL_OUTPUT)
        sys_lst = [self.noise, self.disturbance, self.cl_input,
             self.cl_output, self.sum_N_Y, self.sum_D_U]
        # Construct connections
        connections=[
                [SUM_D_U_D, DISTURBANCE_OUT],
                [SUM_N_Y_N, NOISE_OUT],
                [CL_OUTPUT_IN, SUM_N_Y_OUT],
              ]
        return ConnectionAssembly(
              con=connections,
              sys=sys_lst,
              inp=[CL_INPUT_REF],
              out=[CL_INPUT_REF, CL_OUTPUT_OUT],
              )

    def makeStepResponse(self, time=0, step_size=1, **time_opts):
        """
        Simulates the step response.

        Parameters
        ----------
        time: float
            time used for state in the ControlSBML system
        time_opts: Options
            start_time, end_time, points_per_time
        
        Returns
        -------
        Timeseries
        """
        self.factory.dt = time_opts.get("dt", 1.0/cn.POINTS_PER_TIME)
        self.factory.initializeLoggers()
        times = ctl.makeSimulationTimes(**time_opts)
        X0 = makeStateVector(self.closed_loop_system, start_time=time)
        timeresponse = control.input_output_response(self.closed_loop_system,
              times, U=step_size, X0=X0)
        return util.timeresponse2Timeseries(timeresponse)
  

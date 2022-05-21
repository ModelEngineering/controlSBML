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
  3 SUM_R_F: adds the reference signal to the filter output
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
  7 SUM_U_D: Summation of controller output and disturbance
    * u (in): controller output 
    * d (in): disturbance
    * sum (out): summation
  8 NOISE
    * out (out): noise generated
  9 SUM_Y_N: Summation of system output and noise
    * n (in): noise
    * y (in): system output
    * sum (out): closed loop output
 10 CL_OUTPUT: output from closed loop system
    * in (in): measured output plus noise (if present)
    * out (out)

"""

import controlSBML as ctl
from controlSBML.simulate_system import makeStateVector
from controlSBML import msgs
from controlSBML import util

import collections
import control
import numpy as np
import pandas as pd


CL_INPUT = "cl_input"
CL_INPUT_REF = "cl_input.ref"
CL_INPUT_IN = "cl_input.in"
CL_INPUT_OUT = "cl_input.out"
CL_OUTPUT = "cl_output"
CL_OUTPUT_IN = "cl_output.in"
CL_OUTPUT_OUT = "cl_output.out"
CONTROLLER = "controller"
CONTROLLER_ERR = "controller.err"
CONTROLLER_IN = "controller.in"
CONTROLLER_OUT = "controller.out"
CONTROLLER_U = "controller.u"
DISTURBANCE = "disturbance"
DISTURBANCE_OUT = "disturbance.out"
FLTR = "fltr"
FLTR_IN = "fltr.in"
FLTR_OUT = "fltr.out"
NOISE = "noise"
NOISE_OUT = "noise.out"
SUM_R_F = "sum_R_F"
SUM_R_F_REF = "sum_R_F.ref"
SUM_R_F_FLTR = "sum_R_F.fltr"
SUM_R_F_OUT = "sum_R_F.out"
SUM_U_D = "sum_U_D"
SUM_U_D_U = "sum_U_D.u"
SUM_U_D_D = "sum_U_D.d"
SUM_U_D_OUT = "sum_U_D.out"
SUM_Y_N = "sum_Y_N"
SUM_Y_N_Y = "sum_Y_N.y"
SUM_Y_N_N = "sum_Y_N.n"
SUM_Y_N_OUT = "sum_Y_N.out"
SYSTEM = "system"


# sys: dict with values sys_list  -list-NonlinearIOSystem
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
        self.sum_Y_N = None
        self.sum_R_F = None
        self.sum_U_D = None
    
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

    def makePIDClosedLoopSystem(self, name,
          kp=1, ki=0, kd=0,                       # Controller parameters
          disturbance_amp=0, disturbance_frq=0,   # Disturbance
          noise_amp=0, noise_frq=0,               # Noise
          kf=None,                                # Filter
          is_neg_feedback=True,                   # Negative feedback?
          system_input=None, system_name=None,    # SISO input, output
          system_output=None,
          closed_loop_outputs=None):              # list of outputs from closed loop system
        """
        Creates a closed loop system for a ControlSBML object. The closed loop system
        includes a PID controller and a filter.

        Parameters
        ----------
        name: str
            Name of the resulting system.
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
        is_neg_feedback: bool
            subtract filter output from reference
        system_input: str
            name of the input used in the SISO system
            the name must be present in the ControlSBML object
        system_output: str
            name of the output used in the SISO system
            the name must be present in the ControlSBML object
        closed_loop_ouputs: list-str
            If None, cl_input.in, cl_output.out
        
        Returns
        -------
        control.IOSystem.InterconnectedSystem
           cl_input.in: reference input to System
           exit.out: output from the System
           noise.out: noise input to system
           disturbance.out: noise input to system
        """
        # Can only do one build per incovation
        if self.closed_loop_system is not None:
            msg = "A closed loop system exists. "
            msg += "Create a new instance to build another one."
            raise ValueError(msg)
        # Initializations
        if closed_loop_outputs is None:
            closed_loop_outputs = [CL_INPUT_OUT, CL_OUTPUT_OUT]
        # Create the elements of the feedback loop
        self.cl_input = self.factory.makePassthru(CL_INPUT)
        self.cl_output = self.factory.makePassthru(CL_OUTPUT)
        self.noise = self.factory.makeSinusoid("noise", noise_amp, noise_frq)
        self.disturbance = self.factory.makeSinusoid("disturbance", 
              disturbance_amp, disturbance_frq)
        if system_input is None:
            system_input = self.ctlsb.input_names[0]
        if system_output is None:
            system_output = self.ctlsb.output_names[-1]
        ctlsb = ctl.ControlSBML(self.ctlsb.model_reference,
              input_names=[system_input], output_names=[system_output])
        self.system = ctlsb.makeNonlinearIOSystem(SYSTEM)
        self.controller = self.factory.makePIDController(CONTROLLER,
              kp=kp, ki=ki, kd=kd)
        if kf is None:
            self.fltr = self.factory.makePassthru(FLTR)
        else:
            self.fltr = self.factory.makeFilter(FLTR, kf)
        self.sum_Y_N = self.factory.makeAdder(SUM_Y_N, input_names=["y", "n"])
        self.sum_U_D = self.factory.makeAdder(SUM_U_D, input_names=["u", "d"])
        self.sum_R_F = self.factory.makeAdder(SUM_R_F, input_names=["ref", "fltr"])
        # Construct the interconnected system
        pfx = ""
        if is_neg_feedback:
            pfx = "-"
        filter_str = "%s%s" % (pfx, FLTR_OUT)
        system_inp = "%s.%s" % (SYSTEM, system_input)
        system_out = "%s.%s" % (SYSTEM, system_output)
        sys_lst = [self.noise, self.disturbance,
              self.cl_input, self.cl_output,
              self.sum_Y_N, self.sum_R_F, self.sum_U_D,
              self.system, self.fltr, self.controller]
        connections=[
          [SUM_R_F_REF, CL_INPUT_OUT],        # r(t)
          [SUM_R_F_FLTR, filter_str],
          [CONTROLLER_IN, SUM_R_F_OUT],    # e(t)
          [SUM_U_D_U, CONTROLLER_OUT],        # u(t)
          [SUM_U_D_D, DISTURBANCE_OUT],  # d(t)
          [system_inp,   SUM_U_D_OUT],
          [SUM_Y_N_Y, system_out],        # y(t)
          [SUM_Y_N_N, NOISE_OUT],        # n(t)
          [FLTR_IN,     SUM_Y_N_OUT],
          [CL_OUTPUT_IN, SUM_Y_N_OUT],
        ]
        self.closed_loop_system = control.interconnect(sys_lst,
              connections=connections,
              inplist=[CL_INPUT_IN],
              outlist=closed_loop_outputs,
            )

    def newmakePIDClosedLoopSystem(self, name,
          kp=1, ki=0, kd=0,                       # Controller parameters
          disturbance_amp=0, disturbance_frq=0,   # Disturbance
          noise_amp=0, noise_frq=0,               # Noise
          kf=None,                                # Filter
          is_neg_feedback=True,                   # Negative feedback?
          system_input=None, system_name=None,    # SISO input, output
          system_output=None,
          closed_loop_outputs=None):              # list of outputs from closed loop system
        """
        Creates a closed loop system for a ControlSBML object. The closed loop system
        includes a PID controller and a filter.

        Parameters
        ----------
        name: str
            Name of the resulting system.
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
        is_neg_feedback: bool
            subtract filter output from reference
        system_input: str
            name of the input used in the SISO system
            the name must be present in the ControlSBML object
        system_output: str
            name of the output used in the SISO system
            the name must be present in the ControlSBML object
        closed_loop_ouputs: list-str
            If None, cl_input.in, cl_output.out
        
        Returns
        -------
        control.IOSystem.InterconnectedSystem
           cl_input.in: reference input to System
           exit.out: output from the System
           noise.out: noise input to system
           disturbance.out: noise input to system
        """
        assembly = self._makeClosedLoopCommon(name,
              disturbance_amp=disturbance_amp,
              disturbance_frq=disturbance_frq,
              noise_amp=noise_amp,
              noise_frq=noise_frq,
              closed_loop_outputs=closed_loop_outputs,
              )
        # Make controller, system, filter
        self.sum_R_F = self.factory.makeAdder(SUM_R_F, input_names=["ref", "fltr"])
        #
        if system_input is None:
            system_input = self.ctlsb.input_names[0]
        if system_output is None:
            system_output = self.ctlsb.output_names[-1]
        ctlsb = ctl.ControlSBML(self.ctlsb.model_reference,
              input_names=[system_input], output_names=[system_output])
        self.system = ctlsb.makeNonlinearIOSystem(SYSTEM)
        #
        self.controller = self.factory.makePIDController(CONTROLLER,
              kp=kp, ki=ki, kd=kd)
        #
        if kf is None:
            self.fltr = self.factory.makePassthru(FLTR)
        else:
            self.fltr = self.factory.makeFilter(FLTR, kf)
        # Connect the pieces
        sys_lst = list(assembly.sys.values())
        sys_lst.extend([self.system, self.controller, self.fltr, self.sum_R_F])
        # Additional names
        pfx = ""
        if is_neg_feedback:
            pfx = "-"
        filter_str = "%s%s" % (pfx, FLTR_OUT)
        system_inp = "%s.%s" % (SYSTEM, system_input)
        system_out = "%s.%s" % (SYSTEM, system_output)
        # Add missing connections
        connections = list(assembly.con)
        connections.append([SUM_R_F_FLTR, CL_INPUT_OUT])
        connections.append([CONTROLLER_ERR, SUM_R_F_OUT])
        connections.append([FLTR_IN, CL_OUTPUT_OUT])
        connections.append([SUM_R_F_REF, filter_str])
        # Construct the interconnected system
        self.closed_loop_system = control.interconnect(sys_lst,
              connections=connections,
              inplist=[CL_INPUT_REF],
              outlist=assembly.out,
            )

    # TODO: _makeNoiseAndDisturbance
    def _makeClosedLoopCommon(self, name,
          disturbance_amp=0, disturbance_frq=0,   # Disturbance
          noise_amp=0, noise_frq=0,               # Noise
          system_input=None, system_name=None,    # SISO input, output
          system_output=None,
          closed_loop_outputs=None):              # list of outputs from closed loop system
        """
        Creates parts of the SISO closed loop system.
        Makes the following components and their internal connections:
            cl_input: closed loop input
            disturbance
            noise
            SUM_Y_N
            SUM_U_D
            cl_output: closed loop output
        The missing connections are (a) cl_input to system, (b) system to filter/controller.

        Parameters
        ----------
        name: str
            Name of the resulting system.
        noise_amp: float
            amplitude of the sine wave for noise
        noise_frq: float
            frequency of the sine wave for noise
        noise_amp: float
            amplitude of the sine wave for noise
        noise_frq: float
            frequency of the sine wave for noise
        closed_loop_ouputs: list-str
            If None, cl_output.out, exit.out
        
        Returns
        -------
        ConnectionAssembly
        """
        # Initializations
        if closed_loop_outputs is None:
            closed_loop_outputs = [CL_INPUT_OUT, CL_OUTPUT_OUT]
        # Create the elements of the feedback loop
        self.cl_input = self.factory.makePassthru(CL_INPUT)
        self.cl_output = self.factory.makePassthru(CL_OUTPUT)
        self.noise = self.factory.makeSinusoid(NOISE, noise_amp, noise_frq)
        self.disturbance = self.factory.makeSinusoid(DISTURBANCE,
              disturbance_amp, disturbance_frq)
        self.sum_Y_N = self.factory.makeAdder(SUM_Y_N, input_names=["y", "d"])
        self.sum_U_D = self.factory.makeAdder(SUM_U_D, input_names=["u", "d"])
        # Construct the interconnected system
        sys_dct = {NOISE: self.noise, 
              DISTURBANCE: self.disturbance,
              CL_INPUT: self.cl_input,
              CL_OUTPUT:  self.cl_output,
              SUM_Y_N: self.sum_Y_N,
              SUM_U_D: self.sum_U_D,
              }
        # Define system inputs and outputs
        if system_input is None:
            system_input = self.ctlsb.input_names[0]
        if system_output is None:
            system_output = self.ctlsb.output_names[-1]
        system_inp = "%s.%s" % (SYSTEM, system_input)
        system_out = "%s.%s" % (SYSTEM, system_output)
        # Construct connections
        connections=[
                [SUM_U_D_U, CONTROLLER_U],
                [SUM_U_D_D, DISTURBANCE_OUT],
                [system_inp, SUM_U_D_OUT],
                [SUM_Y_N_Y, system_out],
                [SUM_Y_N_N, NOISE_OUT],
                [CL_OUTPUT_IN, SUM_Y_N_OUT],
              ]
        return ConnectionAssembly(
              con=connections,
              sys=sys_dct,
              inp=[CL_INPUT_REF],
              out=closed_loop_outputs
              )

    def _makeDisturbanceNoise(self,
          disturbance_amp=0, disturbance_frq=0,   # Disturbance
          noise_amp=0, noise_frq=0):              # Noise
        """
        Creates the noise and disturbance elements.
            disturbance
            noise
            SUM_Y_N
            SUM_U_D

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
        self.noise = self.factory.makeSinusoid(NOISE, noise_amp, noise_frq)
        self.disturbance = self.factory.makeSinusoid(DISTURBANCE,
              disturbance_amp, disturbance_frq)
        self.sum_Y_N = self.factory.makeAdder(SUM_Y_N, input_names=["y", "d"])
        self.sum_U_D = self.factory.makeAdder(SUM_U_D, input_names=["u", "d"])
        # Construct the interconnected system
        sys_dct = {NOISE: self.noise, 
              DISTURBANCE: self.disturbance,
              SUM_Y_N: self.sum_Y_N,
              SUM_U_D: self.sum_U_D,
              }
        # Construct connections
        connections=[
                [SUM_U_D_D, DISTURBANCE_OUT],
                [SUM_Y_N_N, NOISE_OUT],
              ]
        return ConnectionAssembly(
              con=connections,
              sys=sys_dct,
              inp=None,
              out=None
              )

    def makeStepResponse(self, time=0, step_size=1, **sim_opts):
        """
        Simulates the step response.

        Parameters
        ----------
        time: float
            time used for state in the ControlSBML system
        sim_opts: Options
            start_time, end_time, points_per_time
        
        Returns
        -------
        Timeseries
        """
        self.factory.initializeLoggers()
        times = ctl.makeSimulationTimes(**sim_opts)
        X0 = makeStateVector(self.closed_loop_system, start_time=time)
        timeresponse = control.input_output_response(self.closed_loop_system,
              times, U=step_size, X0=X0)
        return util.timeresponse2Timeseries(timeresponse)
  

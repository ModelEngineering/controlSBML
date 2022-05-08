"""Builds and Evaluates SISO Systems"""

"""
The SISO closed loop system has:
  * A ControlSBML SISO System
  * PID Controller
  * Filter (with dcgain = 1)

The design parameters of the system are:
  * PID parameters: kp, ki, kd
  * Filter parameter: constant

Evaluations provided are:
  * Find dcgain of ControlSBML system with many inputs and one output
  * DC gain of the closed loop system
  * Evaluate the closed loop system for disturbances and noise
"""

import controlSBML as ctl
from controlSBML.simulate_system import makeStateVector
from controlSBML import msgs
from controlSBML import util

import control
import numpy as np
import pandas as pd


class SISOClosedLoopSystem(object):
    # Constructs and evaluates a SISO closed loop system.

    def __init__(self, ctlsb):
        self.ctlsb = ctlsb
        self.factory = ctl.IOSystemFactory()
        self.closed_loop_system = None
        self.component_names = ["system", "controller", "fltr", 
              "reference", "noise",
              "disturbance", "sum_Y_N", "sum_R_R", "sum_U_D"]
        # Elements in the closed loop system 
        ## Subsysems
        self.system = None
        self.controller = None
        self.fltr = None
        ## External signals
        self.entry = None
        self.exit = None
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

    def makeClosedLoopSystem(self, name,
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
            If None, entry.out, exit.out
        
        Returns
        -------
        control.IOSystem.InterconnectedSystem
           entry.in: reference input to System
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
            closed_loop_outputs = ["entry.out", "exit.out"]
        # Create the elements of the feedback loop
        self.entry = self.factory.makePassthru("entry")
        self.exit = self.factory.makePassthru("exit")
        self.noise = self.factory.makeSinusoid("noise", noise_amp, noise_frq)
        self.disturbance = self.factory.makeSinusoid("disturbance", 
              disturbance_amp, disturbance_frq)
        if system_input is None:
            system_input = self.ctlsb.input_names[0]
        if system_output is None:
            system_output = self.ctlsb.output_names[-1]
        ctlsb = ctl.ControlSBML(self.ctlsb.model_reference,
              input_names=[system_input], output_names=[system_output])
        self.system = ctlsb.makeNonlinearIOSystem("system")
        self.controller = self.factory.makePIDController("controller",
              kp=kp, ki=ki, kd=kd)
        if kf is None:
            self.fltr = self.factory.makePassthru("fltr")
        else:
            self.fltr = self.factory.makeFilter("fltr", kf)
        self.sum_Y_N = self.factory.makeAdder("sum_Y_N")
        self.sum_U_D = self.factory.makeAdder("sum_U_D")
        self.sum_R_F = self.factory.makeAdder("sum_R_F")
        # Construct the interconnected system
        if is_neg_feedback:
            filter_str = "-fltr.out"
        else:
            filter_str = "fltr.out"
        system_inp = "system.%s" % system_input
        system_out = "system.%s" % system_output
        sys_lst = [self.noise, self.disturbance,
              self.entry, self.exit,
              self.sum_Y_N, self.sum_R_F, self.sum_U_D,
              self.system, self.fltr, self.controller]
        self.closed_loop_system = control.interconnect(sys_lst,
              connections=[
                ["sum_R_F.in2", "entry.out"],        # r(t)
                ['controller.in', 'sum_R_F.out'],    # e(t)
                ['sum_U_D.in1', 'controller.out'],   # u(t)
                ['sum_U_D.in2', 'disturbance.out'],  # d(t)
                [system_inp,   'sum_U_D.out'],
                ['sum_Y_N.in1', system_out],        # y(t)
                ['sum_Y_N.in2', 'noise.out'],        # n(t)
                ['fltr.in',     'sum_Y_N.out'],
                ["exit.in", "sum_Y_N.out"],
                ['sum_R_F.in1', filter_str],
              ],
              inplist=["entry.in"],
              outlist=closed_loop_outputs,
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
  

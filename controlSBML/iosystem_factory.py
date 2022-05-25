"""Factory for making IOSystem objects"""

"""
Creates IOsystems to be used in constructing closed loop systems.
Inputs have the form in, in1, in2, ...
Outpus have the form out, out1, out2, ... 

IOSystems created:
    makeFilter: creates an exponential filter
    makeConstant: outputs a constant value
    makeAdder: creates an IOSystem that outputs the sum of the inputs
    makePassthur: creates an IOSystem that outputs the input
    makePIDController: creates a PID controller
    makeFullStateController: creates a PID controller
    makeSinusoid: creates an IOSystem that outputs a sinusoid
"""

from controlSBML import constants as cn
from controlSBML import logger as lg
from controlSBML import msgs

import control
import numpy as np
import pandas as pd

IN = "in"
IN1 = "in1"
IN2 = "in2"
OUT = "out"
REF = "ref"
STATE = "state"


class IOSystemFactory(object):

    def __init__(self, name="factory", dt=1/cn.POINTS_PER_TIME):
        """
        Parameters
        ----------
        name: str
            Name of the InterconnectedSystem
        dt: float
            delta between simulation times
        """
        self.name = name
        self.dt = dt
        self.registered_names = []  # List of names logged
        self.loggers = []  # Loggers for created objects

    def _registerLogger(self, logger_name, item_names):
        is_duplicate = any([logger_name == l.logger_name for l in self.loggers])
        if is_duplicate:
            raise ValueError("Duplicate logger name: %s" % logger_name)
        logger = lg.Logger(logger_name, item_names)
        self.loggers.append(logger)
        return logger

    @staticmethod
    def _array2scalar(vec):
        """
        Converts an array with a single element to a scalar.

        Parameters
        ----------
        vec: float/array-like
        
        Returns
        -------
        float
        """
        try:
            val = vec[0]
            if len(vec) > 1:
                raise RuntimeError("Array has more than 1 element.")
        except:
            val = vec
        return val

    def initializeLoggers(self):
        """
        Initialize loggers.
        """
        for logger in self.loggers:
            logger.initialize()

    def report(self):
        """
        Creates of report from the logs of objects produced by this factory.
        
        Returns
        -------
        pd.DataFrame
        """
        if len(self.loggers) == 0:
            return None
        if len(self.loggers) == 1:
            return self.loggers[0].report()
        # Handle >1 elements
        logger = self.loggers[0]
        others = self.loggers[1:]
        merged_logger = logger.merge(self.name, others)
        return merged_logger.report()

    def makePIDController(self, name, kp=2, ki=0, kd=0):
        """
        Creates a PID controller.
        NonlinearIOSystem with input IN and output OUT.
        States are:
            culmulative_err: cumulative control error
            last_err: last control error
        
        Parameters
        ----------
        name: str
            Name of the system
        kp: float
           proportional control constant
        ki: float
           integral control constant
        kd: float
           differential control constant.
        
        Returns
        -------
        control.NonlinearIOSystem
        """
        #
        logger = self._registerLogger(name, [IN, "last_err",
              "acc_err", OUT])
        def updfcn(time, x_vec, u_vec, __):
            # u: float (error signal)
            # x_vec[0] - last_u
            # x_vec[1] - accumulated_u
            u_val = self._array2scalar(u_vec)
            daccumulated_u = u_val/self.dt
            dlast_u_val = (u_val - x_vec[0])/self.dt
            return dlast_u_val, daccumulated_u
        #
        def outfcn(time, x_vec, u_vec, __):
            # u: float (error signal)
            u_val = self._array2scalar(u_vec)
            last_u_val = x_vec[0]
            accumulated_u_val = x_vec[1]
            control_out =  kp*u_val  \
                          + ki*accumulated_u_val \
                          + kd*(u_val - last_u_val)/self.dt
            logger.add(time, [u_val, x_vec[0], x_vec[1], control_out])
            return control_out
        #
        return control.NonlinearIOSystem(
            updfcn, outfcn, inputs=[IN], outputs=[OUT],
            states=["last_err", "accumulated_err"],
            name=name)

    def makeFullStateController(self, name, ctlsb, factor=1.0, poles=-2, time=0):
        """
        Creates a full state feedback controller for an SBML model
        where the system is linearized at the specified time.
        
        Parameters
        ----------
        name: str
            Name of the system
        ctlsb: ControlSBML
            SISO system
        factor: float
            factor for adjusting the reference input to get a closed loop
            transfer function of 1
        poles: list-float/float
            Desired poles; single float is the dominant pole
        time: float
            Time where system is lienarized
        
        Returns
        -------
        control.NonlinearIOSystem
           inputs:
               <state_variable> (excluding the input)
               ref: reference input
           outputs:
               out
        """
        # Validity Checks
        if len(ctlsb.input_names) != 1:
            msg = "SBML model must have a single input. Has: %s"  \
                  % str(ctlsb.input_names)
            raise ValueError(msg)
        # Initializations
        state_space = ctlsb.makeStateSpace(time=time)
        controller_input_names = [n for n in ctlsb.state_names
             if not n in ctlsb.input_names]
        controller_input_names.insert(0, REF)  # first input
        num_state_input = len(controller_input_names) - 1
        is_distinct_poles = True
        try:
            _ = len(poles)
            if len(set(poles)) < len(poles):
                is_distinct_poles = False
            poles = np.array(poles)
        except:
            # Insert that poles are distinct
            poles = np.array([poles + 0.1*n for n in range(num_state_input)])
        if not is_distinct_poles:
            raise ValueError("Poles must be distinct. Not: %s" % str(poles))
        # Calculate the gain matrix
        kp_vec = control.place(state_space.A, state_space.B, poles)
        def outfcn(time, _, u_vec, __):
            # u_vec: list-float - reference, state variables
            ref = factor*u_vec[0]
            arr = np.array(u_vec[1:])
            output = ref - kp_vec.dot(arr)
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, inputs=controller_input_names, outputs=['out'],
            name=name)

    def makeSinusoid(self, name, amp, frequency):
        """
        Construct a NonlinearIOSystem that outputs a sinusoid.
        The system has output OUT.
        
        Parameters
        ----------
        name: str
        amp: float
        frequency: float
        
        Returns
        -------
        NonlinearIOSystem
        """
        logger = self._registerLogger(name, [OUT])
        def outfcn(time, x_vec, u_vec, _):
            """
            Creates a sine wave at the frequency specified.
    
            Parameters
            ----------
            time: float
            """
            output = amp*np.sin(time*frequency)
            logger.add(time, [output])
            dct = {"amp": amp, "frequency": frequency}
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=[OUT], inputs=[],
            name=name)

    def makeAdder(self, name, num_input=2, input_names=None):
        """
        Inputs two or more elements. Outputs their sum. Name is "sum".
        The inputs are IN1, IN2, ... or the names in input_names.
        The output is OUT.
        
        Parameters
        ----------
        name: str
        num_input: int
        num_names: list-str
        
        Returns
        -------
        NonlinearIOSystem
        """
        if input_names is None:
            input_names = ["%s%d" % (IN, n) for n in range(1, num_input+1)]
        item_names = list(input_names)
        item_names.append(OUT)
        logger = self._registerLogger(name, item_names)
        dct = {}
        def outfcn(time, x_vec, u_vec, _):
            """
            Creates a sine wave at the frequency specified.
    
            Parameters
            ----------
            u: float, float
            """
            output = np.sum(u_vec)
            item_values = list(u_vec)
            item_values.append(output)
            logger.add(time, item_values)   
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, inputs=input_names, outputs=[OUT], name=name)

    def makeFilter(self, name, constant):
        """
        Construct a NonlinearIOSystem for x' = a*x + u, where u is the
        filter input.
        The system has input IN and output OUT. 
        The output is normalized so that the DC gain of the filter is 1.
        
        Parameters
        ----------
        name: str
        constant: float
            e**expo_constant*time
        
        Returns
        -------
        NonlinearIOSystem
        """
        logger = self._registerLogger(name, [IN, STATE, OUT])
        dct = {"constant": constant}
        def updfcn(time, x_vec, u_vec, _):
            """
            Returns the derivative of the state.
    
            Parameters
            ----------
            time: float
            x_vec: float
            u_vec: float
            """
            x_val = self._array2scalar(x_vec)
            u_val = self._array2scalar(u_vec)
            dx_val = constant*x_val + u_val
            dct = {}
            return dx_val
        #
        def outfcn(time, x_vec, u_vec, _):
            u_val = self._array2scalar(u_vec)
            x_val = self._array2scalar(x_vec)
            output = -constant*x_val
            logger.add(time, [u_val, x_val, output])
            return output
        #
        return control.NonlinearIOSystem(
            updfcn, outfcn, outputs=[OUT], inputs=[IN], states=[STATE],
            name=name)

    def makeConstant(self, name, constant):
        """
        Outputs a constant value.
        
        Parameters
        ----------
        name: str
        constant: float
        
        Returns
        -------
        NonlinearIOSystem
        """
        logger = self._registerLogger(name, [OUT])
        dct = {"constant": constant}
        def outfcn(time, x_vec, u_vec, _):
            output = constant
            logger.add(time, [output])
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=[OUT], inputs=[], name=name)

    def makePassthru(self, name, input_name=None, output_name=None):
        """
        Makes a pass through system that outputs its input.
        
        Parameters
        ----------
        name: str
        input_name: str
        output_name: str
        
        Returns
        -------
        NonlinearIOSystem
        """
        logger = self._registerLogger(name, [IN, OUT])
        dct = {}
        if input_name is None:
            input_name = IN
        if output_name is None:
            output_name = OUT
        def outfcn(time, x_vec, u_vec, ___):
            u_val = self._array2scalar(u_vec)
            output = u_val
            logger.add(time, [u_val, output])
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=[output_name], inputs=[input_name], name=name)

    def makeMultiplier(self, name, factor):
        """
        Outputs a multiple of the input signal.
        
        Parameters
        ----------
        name: str
        factor: float
        
        Returns
        -------
        NonlinearIOSystem
        """
        logger = self._registerLogger(name, [IN, OUT])
        dct = {"factor": factor}
        def outfcn(time, x_vec, u_vec, _):
            u_val = self._array2scalar(u_vec)
            output = factor*u_val
            logger.add(time, [u_val, output])
            return output
        
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=[OUT], inputs=[IN], name=name)

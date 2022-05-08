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
    makeSinusoid: creates an IOSystem that outputs a sinusoid
"""

from controlSBML import logger as lg
from controlSBML import msgs

import control
import numpy as np
import pandas as pd

IN = "in"
IN1 = "in1"
IN2 = "in2"
OUT = "out"
STATE = "state"


class IOSystemFactory(object):

    def __init__(self, name="factory"):
        self.name = name
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
        X_CUMU = 0  # cumulative control error
        X_LAST = 1 # last control error
        X_TIME = 2
        last_time = 0
        logger = self._registerLogger(name, [IN, "last_err",
              "acc_err", OUT])
        def updfcn(time, x_vec, u_vec, _):
            # Calculate the derivative of the state variables
            u_val = self._array2scalar(u_vec)
            last_time = x_vec[X_TIME]
            last_err = x_vec[X_LAST]
            last_cumu = x_vec[X_CUMU]
            time = max(time, 1e-8)
            # Calculate derivatives
            dtime = time - last_time
            dcumu = u_val/dtime
            derr = (u_val - last_err)/dtime
            return dcumu, derr, dtime
        #
        def outfcn(time, x_vec, u_vec, __):
            # u: float (error signal)
            u_val = self._array2scalar(u_vec)
            cumu_err = x_vec[X_CUMU]
            last_err = x_vec[X_LAST]
            control_out =  kp*u_val  \
                          + ki*cumu_err \
                          + kd*(u_val - last_err)
            logger.add(time, [u_val, x_vec[1], x_vec[0], control_out])
            return control_out
        #
        return control.NonlinearIOSystem(
            updfcn, outfcn, inputs=[IN], outputs=[OUT],
            states=["cumulative_error", "last_error", "last_time"],
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
        def outfcn(time, _, __, ___):
            """
            Creates a sine wave at the frequency specified.
    
            Parameters
            ----------
            time: float
            """
            output = amp*np.sin(time*frequency)
            logger.add(time, [output])
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=[OUT], inputs=[],
            name=name)

    def makeAdder(self, name, num_input=2):
        """
        Inputs two or more elements. Outputs their sum. Name is "sum".
        The inputs are IN1, IN2, ...
        The output is OUT.
        
        Parameters
        ----------
        name: str
        
        Returns
        -------
        NonlinearIOSystem
        """
        inputs = ["%s%d" % (IN, n) for n in range(1, num_input+1)]
        item_names = list(inputs)
        item_names.append(OUT)
        logger = self._registerLogger(name, item_names)
        def outfcn(time, __, u_vec, ___):
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
            None, outfcn, inputs=inputs, outputs=[OUT], name=name)

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
        def updfcn(time, x_vec, u_vec, ___):
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
        def outfcn(time, _, __, ___):
            output = constant
            logger.add(time, [output])
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=[OUT], inputs=[], name=name)

    def makePassthru(self, name):
        """
        Makes a pass through system that outputs its input.
        
        Parameters
        ----------
        name: str
        
        Returns
        -------
        NonlinearIOSystem
        """
        logger = self._registerLogger(name, [IN, OUT])
        def outfcn(time, _, u_vec, ___):
            u_val = self._array2scalar(u_vec)
            output = u_val
            logger.add(time, [u_val, output])
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=[OUT], inputs=[IN], name=name)

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
        def outfcn(time, _, u_vec, __):
            u_val = self._array2scalar(u_vec)
            output = factor*u_val
            logger.add(time, [u_val, output])
            return output
        
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=[OUT], inputs=[IN], name=name)

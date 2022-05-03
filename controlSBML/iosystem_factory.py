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


from controlSBML import msgs

import control
import numpy as np
import pandas as pd


class IOSystemFactory(object):

    def __init__(self):
        pass

    def makePIDController(self, name, kp=2, ki=0, kd=0):
        """
        Creates a PID controller.
        NonlinearIOSystem with input "in" and output "out".
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
        last_time = 0
        def updfcn(time, x_vec, u_vec, __):
            # Calculate the derivative of the state variables
            if not "len" in dir(u_vec):
                u_vec = [u_vec]
            last_err = x_vec[X_LAST]
            time = max(time, 1e-7)
            d_cumu = u_vec[0]/time
            d_last = (u_vec[0] - last_err)/time
            return d_cumu, d_last
        #
        def outfcn(_, x_vec, u_vec, __):
            # u: float (error signal)
            if not isinstance(x_vec, np.ndarray):
                x_vec = [u_vec]
            if not isinstance(u_vec, np.ndarray):
                u_vec = [u_vec]
            cumu_err = x_vec[X_CUMU]
            last_err = x_vec[X_LAST]
            control_out =  kp*u_vec[0]  \
                          + ki*cumu_err \
                          + kd*(u_vec[0] - last_err)
            return control_out
        #
        return control.NonlinearIOSystem(
            updfcn, outfcn, inputs=['in'], outputs=['out'],
            states=["cumulative_error", "last_error"],
            name=name)

    def makeSinusoid(self, name, amp, frequency):
        """
        Construct a NonlinearIOSystem that outputs a sinusoid.
        The system has output "out".
        
        Parameters
        ----------
        name: str
        amp: float
        frequency: float
        
        Returns
        -------
        NonlinearIOSystem
        """
        def outfcn(time, _, __, ___):
            """
            Creates a sine wave at the frequency specified.
    
            Parameters
            ----------
            time: float
            """
            return amp*np.sin(time*frequency)
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=['out'], inputs=[],
            name=name)

    def makeAdder(self, name, num_input=2):
        """
        Inputs two or more elements. Outputs their sum. Name is "sum".
        The inputs are "in1", "in2", ...
        The output is "out".
        
        Parameters
        ----------
        name: str
        
        Returns
        -------
        NonlinearIOSystem
        """
        inputs = ["in%d" % n for n in range(1, num_input+1)]
        def outfcn(_, __, u_vec, ___):
            """
            Creates a sine wave at the frequency specified.
    
            Parameters
            ----------
            u: float, float
            """
            return np.sum(u_vec)
        #
        return control.NonlinearIOSystem(
            None, outfcn, inputs=inputs, outputs=['out'], name=name)

    def makeFilter(self, name, constant):
        """
        Construct a NonlinearIOSystem for x' = a*x + u, where u is the
        filter input.
        The system has input "in" and output "out". 
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
        def updfcn(time, x_vec, u_vec, ___):
            """
            Returns the derivative of the state.
    
            Parameters
            ----------
            time: float
            x_vec: float
            u_vec: float
            """
            if not "len" in dir(x_vec):
                x_vec = [x_vec]
            if not "len" in dir(u_vec):
                u_vec = [u_vec]
            return constant*x_vec[0] + u_vec[0]
        #
        def outfcn(_, x_vec, __, ___):
            return -constant*x_vec[0]
        #
        return control.NonlinearIOSystem(
            updfcn, outfcn, outputs=['out'], inputs=["in"], states=["x"],
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
        def outfcn(_, __, ___, ____):
            return constant
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=['out'], inputs=[], name=name)

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
        def outfcn(_, __, u_vec, ____):
            return u_vec
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=['out'], inputs=['in'], name=name)

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
        def outfcn(_, __, u_vec, ____):
            if not "len" in dir(u_vec):
                u_vec = [u_vec]
            return [factor*u_vec[0]]
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=['out'], inputs=['in'], name=name)

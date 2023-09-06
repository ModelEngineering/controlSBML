"""Designs a closed loop SISO System with a PID controller and Filter."""

import controlSBML.constants as cn

import numpy as np
import control
import matplotlib.pyplot as plt


class SISOClosedLoopDesigner(object):

    def __init__(self, sys_tf):
        """
        Args:
            sys_tf: control.TransferFunction (open loop system)
        """
        self.sys_tf = sys_tf 
        self.kp = None
        self.ki = None
        self.kd = None
        self.kf = None

    def design(self,
          kp=None, ki=None, kd=None,              # Controller parameters
          kf=None):                               # Filter
        """
        Args:
            kp (float): proportional gain
            ki (float): integral gain
            kd (float): derivative gain
            kf (float): filter gain
        """
        def update(value, name):
            if value is not None:
                self.__setattr__(name, value)
        #
        update(kp, "kp")
        update(ki, "ki")
        update(kd, "kd")
        update(kf, "kf")
        # Construct the transfer functions
        controller_tf = control.tf([0], [1])
        if self.kp is not None:
            controller_tf += control.tf([self.kp], [1])
        if self.ki is not None:
            controller_tf += control.tf([self.ki], [1, 0])
        if self.kd is not None:
            controller_tf += control.tf([self.kd, 0], [1])
        # Filter
        filter_tf = 1
        if self.kf is not None:
            filter_tf = control.TransfrFunction([self.kf], [1, self.kf])
        # Closed loop --- FIX
        forward_tf = self.sys_tf*controller_tf
        loop_tf = filter_tf*forward_tf
        denom = 1 + np.array(loop_tf.den)
        numer = np.array(forward_tf.num) + np.array(loop_tf.den)
        closed_loop_tf = control.tf(numer, denom)
        #closed_loop_tf = control.feedback(self.sys_tf, controller_tf*filter_tf)
        import pdb; pdb.set_trace()
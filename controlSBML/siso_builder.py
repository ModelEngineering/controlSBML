"""Builds and Evaluates SISO Systems"""

"""
The SISO system has:
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
from controlSBML import msgs

import control
import numpy as np
import pandas as pd


class SISOBuilder(object):

    def __init__(self, ctlsb):
        pass

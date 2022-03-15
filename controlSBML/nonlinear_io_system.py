"""Extends control.NonlinearIOSystem for ControlSBML objects."""

import controlSBML.constants as cn
from controlSBML import msgs

import control
import numpy as np
import pandas as pd


ARG_LST = [cn.TIME, cn.STATE, cn.INPUT, cn.PARAMS]
MIN_ELAPSED_TIME = 1e-2

class NonlinearIOSystem(control.NonlinearIOSystem):

    def __init__(self, name, ctlsb, effector_dct=None,
         do_simulate_on_update=False,  **kwargs):
        """
        Parameters
        ----------
        name: Name of the non-linear system created
        ctlsb: ControlSBML
        effector_dct: dict
            key: input defined in ControlSBML object
            value: roadrunner muteable (e.g., floating species, parameter)
        do_simulate_on_update: bool
            simulate to the current time before each state update
        kwargs: dict (additional keyword arguments provided by caller)

        Usage:
            sys = NonlinearIOSystem(func, None,
                  inputs=["i1", i2"], outputs=["o1"],
                  state=["s1", "s2", "s3"])
            times, y_val1s = control.input_output_response(sys, input_times,
                  input1s, initial_state)
            times, y_val2s = control.input_output_response(sys, input_times,
                  input2s, initial_state)
        """
        self.ctlsb = ctlsb
        self.do_simulate_on_update = do_simulate_on_update
        # Useful properties
        self.state_names = ctlsb.state_names
        self.input_names = ctlsb.input_names
        self.output_names = ctlsb.output_names
        #
        self.effector_dct = effector_dct
        if self.effector_dct is None:
            self.effector_dct = {}
        for input_name in self.input_names:
            if not input_name in self.effector_dct.keys():
                # Make it its own effector
                self.effector_dct[input_name] = input_name
        self._ignored_inputs = []  # Note setable inputs
        # Initialize the controlNonlinearIOSystem object
        super().__init__(self.updfcn, self.outfcn,
              inputs=self.input_names, outputs=self.output_names,
              states=self.state_names, name=name)

    @property
    def outlist(self):
        return ["%s.%s" % (self.name, n) for n in self.output_names]

    @property
    def inplist(self):
        return ["%s.%s" % (self.name, self.effector_dct[n]) for n in self.input_names]

    def setTime(self, time):
        self.ctlsb.setTime(time)

    def getStateSer(self, time=0):
        """
        Gets the values of state at the specified time.

        Parameters
        ----------
        time: float
        
        Returns
        -------
        pd.Series
        """
        self.setTime(time)
        return self.ctlsb.state_ser

    def updfcn(self, time, x_vec, u_vec, _):
        """ 
        Computes the change in state. This is done by having roadrunner
        calculate fluxes. No simulation is run, and so this technique
        may not always work.
        
        Parameters
        ----------
        time: float: time
        x_vec: np.array(float): state vector
        u_vec: np.array(float): input vector (in log10 units)

        Returns
        -------
        np.array(float): change in state vector
        """
        if isinstance(u_vec, float):
            u_vec = [u_vec]
        if self.do_simulate_on_update:
            self.setTime(time)  # Consider time dependent functions
        # Adust the state
        state_dct = {n: x_vec[i] for i, n in enumerate(self.state_names)}
        self.ctlsb.set(state_dct)
        # Set values for input effector and check for errors
        input_dct = {self.effector_dct[n]: u_vec[i]
                     for i, n in enumerate(self.input_names)}
        for input_name, input_value in input_dct.items():
            try:
                self.ctlsb.set({input_name: input_value})
            except RuntimeError:
                if input_name not in self._ignored_inputs:
                    self._ignored_inputs.append(input_name)
                    text = "System %s: Input %s cannot be set. Ignored."  \
                          % (self.name, input_name)
                    msgs.warn(text)
        # Calculate the change in floating species in state
        dstate_names = ["%s'" % n for n in self.state_names]
        dstate_ser = pd.Series(self.ctlsb.get(dstate_names))
        return dstate_ser.values

    def outfcn(self, time, x_vec, u_vec, _): 
        """ 
        Calculates the values of outputs.
        
        Parameters
        ----------
        time: float: time
        x_vec: np.array(float): state vector
        u_vec: np.array(float): input vector (in log10 units)

        Returns
        -------
        np.array(float): change in state vector
        """
        out_vec = np.repeat(np.nan, self.ctlsb.num_output)
        for out_idx, name in enumerate(self.output_names):
            state_idx = self.state_names.index(name)
            out_vec[out_idx] = x_vec[state_idx]
        if np.isnan(np.sum(out_vec)):
            raise ValueError("Outputs could not be calculated.")
        return out_vec

"""Extends control.NonlinearIOSystem to consider ControlSBML objects."""

import controlSBML.constants as cn

import control
import numpy as np
import pandas as pd


ARG_LST = [cn.TIME, cn.STATE, cn.INPUT, cn.PARAMS]
MIN_ELAPSED_TIME = 1e-2

class NonlinearIOSystem(control.NonlinearIOSystem):

    def __init__(self, *pargs, ctlsb=None, effector_dct=None, **kwargs):
        """
        Parameters
        ----------
        pargs: list (positional arguments for control.NonlinearIOSystem)
        ctlsb: ControlSBML
        effector_dct: dict
            key: input defined in ControlSBML object
            value: roadrunner muteable (e.g., floating species, parameter)
        kwargs: dict (keyword arguments for control.NonlinearIOSystem)

        Usage:
            sys = NonlinearIOSystem(func, None,
                  inputs=["i1", i2"], outputs=["o1"],
                  state=["s1", "s2", "s3"])
            times, y_val1s = control.input_output_response(sys, input_times,
                  input1s, initial_state)
            sys.reset()
            times, y_val2s = control.input_output_response(sys, input_times,
                  input2s, initial_state)
        """
        self.ctlsb = ctlsb
        self.effector_dct = effector_dct
        # Set up the calls
        if ctlsb is not None:
            self._updfcn = self.ctlsbUpdfcn
            self._outfcn = self.ctlsbOutfcn
        else:
            self._updfcn = parms[0]
            self._outfcn = parms[1]
        # Initializations
        self.state_call_dct = None
        self.output_call_dct = None
        self._is_first_state_call = None
        self._last_time = None  # Last time at which state update was called
        self.reset()
        # Initialize the controlNonlinearIOSystem object
        newPargs = list(pargs)
        newPargs[0] = self._updfcnWrapper
        newPargs[0] = self._outfcnWrapper
        super.__init__(*pargs, **kwargs)

    def reset(self):
        """
        Sets to initial conditions.
        """
        self.state_call_dct = self._initializeCallDct()
        self.output_call_dct = self._initializeCallDct()
        self._is_first_state_call = True
        self._last_time = 0

    @staticmethod
    def _initializeCallDct():
        """
        Initializes the call dictionary history.
        """
        dct = {cn.EVENT: []}
        for key in ARG_LST:
            dct[key] = []
        dct[cn.RETURNS] = []
        return dct

    @staticmethod
    def _updateCallDct(dct, event, *pargs):
        """
        Updates information in the call dictionary history.
        """
        dct[cn.EVENT].append(event)
        for idx, key in enumerate(
            dct[ARG_LST[idx])  = pargs[idx]

    def _updfcnWrapper(self, *pargs, **kwargs):
        """
        Calls the state update function.

        Parameters
        ----------
        pargs: list (positional arguments passed to updfcn)
        kwargs: dict (keyword arguments passed to updfcn)
        
        Returns
        -------
        list (entries in updated state)
        """
        self._updateCallDct(self.state_call_dct, cn.STATE, *pargs)
        results = self._updfcn(*pargs)
        self._is_first_state_call = False
        dct[cn.RESULTS] = results
        return results

    def _outfcnWrapper(self, *pargs, **kwargs):
        """
        Calls the ouput update function.

        Parameters
        ----------
        pargs: list (positional arguments passed to updfcn)
        kwargs: dict (keyword arguments passed to updfcn)
        
        Returns
        -------
        list (entries in output)
        """
        self._updateCallDct(self.output_call_dct, cn.OUTPUT, *pargs)
        results = self._outfcn(*pargs)
        dct[cn.RESULTS] = results
        return results

    def ctlsbUpdfcn(self, time, x_vec, u_vec, _):
        """
        Computes updated states by running a short simulation.
        
        Parameters
        ----------
        time: float: time
        x_vec: np.array(float): state vector
        u_vec: np.array(float): input vector (in log10 units)
        """
        # Initializations
        if self._is_first_state_call:
            self.ctlsb.roadrunner.reset()
        if time - self._last_time < MIN_ELAPSED_TIME:
            # Verify that time is moving forward and sufficient time has passed
            pass
        else:
            # Update the state and inputs in the simulation
            if not "len" in dir(u_vec):
                u_vec = [u_vec]  # Ensure that this is vector-like
            input_dct = {self.effector_dct[n]: float(u_vec[i])
                         for i, n in enumerate(self.ctlsb.input_names)}
            self.ctlsb.set(input_dct)  # Update input values
            _ = self.ctlsb.roadrunner.simulate(self._last_time, time, 2)
            self._last_time = time
        return [self.ctlsb.get(n) for n in self.ctlsb.state_names]

    def ctlsbOutfcn(*_, **__):
        """
        Extracts the outputs from the simulation.

        Parameters
        ----------
         Arguments are passed but ignored.
        
        Returns
        -------
        np.array
        """
        return self.ctlsb.output_ser.values

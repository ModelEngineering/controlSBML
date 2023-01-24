"""Extends control.NonlinearIOSystem for ControlSBML objects."""

import controlSBML.constants as cn
from controlSBML import logger as lg
from controlSBML import msgs
from controlSBML import util

import control
import numpy as np
import pandas as pd


ARG_LST = [cn.TIME, cn.STATE, cn.INPUT, cn.PARAMS]
MIN_ELAPSED_TIME = 1e-2


class NonlinearIOSystem(control.NonlinearIOSystem):

    def __init__(self, name, ctlsb,
         do_simulate_on_update=False,  **kwargs):
        """
        Parameters
        ----------
        name: Name of the non-linear system created
        ctlsb: ControlSBML
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
        self.name = name
        self.ctlsb = ctlsb
        self.do_simulate_on_update = do_simulate_on_update
        # Useful properties
        self.input_names = ctlsb.input_names
        self.state_names = list(ctlsb.species_names)
        self.dstate_names = ["%s'" % n for n in self.state_names]
        self.output_names = ctlsb.output_names
        self.num_state = len(self.state_names)
        self.num_input = len(self.input_names)
        self.num_output = len(self.output_names)
        #
        self._ignored_inputs = []  # Note setable inputs
        self.logger = lg.Logger(self.name, self.state_names)
        # Initialize the controlNonlinearIOSystem object
        super().__init__(self._updfcn, self._outfcn,
              inputs=self.input_names, outputs=self.output_names,
              states=self.state_names, name=self.name)

    @property
    def outlist(self):
        return ["%s.%s" % (self.name, n) for n in self.output_names]

    def setTime(self, time):
        self.ctlsb.setTime(time)

    def makeStateSer(self, time=0):
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
        return util.makeRoadrunnerSer(self.ctlsb.roadrunner, self.state_names)

    def _updfcn(self, time, x_vec, u_vec, _):
        """
        Computes the change in state by adding the input to its current value.
        This is done by having roadrunner
        calculate fluxes. No simulation is run, and so this technique
        may not always work.
        Inputs that change floating species are viewed as additions to the species,
        not setting the level of the species.

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
        # Calculate the derivatives of floating species in state
        derivative_arr = np.array([v for v in self.ctlsb.get(self.dstate_names).values()])
        # Adjust the derivatives based on the input
        input_dct = {n: u_vec[i] for i, n in enumerate(self.input_names)}
        for idx, input_name in enumerate(input_dct.keys()):
            derivative_arr[idx] += input_dct[input_name]
        return derivative_arr

    def _outfcn(self, time, x_vec, _, __):
        """
        Calculates the values of outputs.

        Parameters
        ----------
        x_vec: np.array(float): state vector
        u_vec: np.array(float): input vector (in log10 units)

        Returns
        -------
        np.array(float): change in state vector
        """
        out_vec = np.repeat(np.nan, self.ctlsb.num_output)
        for out_idx, name in enumerate(self.output_names):
            if name in self.state_names:
                state_idx = self.state_names.index(name)
                out_vec[out_idx] = x_vec[state_idx]
        if np.isnan(np.sum(out_vec)):
            raise ValueError("Outputs could not be calculated.")
        #self.logger.add(time, self.ctlsb.get(self.state_names).values())
        return out_vec

    def _makeStairs(self, num_point, num_step, initial_value, final_value):
        """
        Constructs a monotone sequence of steps.
        
        Parameters
        ----------
        num_point: number of points in the staircase
        num_step: int (number of steps in the stair response.
        start_value: float (initial values of the inputs)
        final_value: float (ending values of the inputs)
        
        Returns
        -------
        np.ndarray
        """
        stairs = []  # Base value of stairs
        num_point_in_step = int(num_point/num_step)
        for num in range(num_step):
            stairs.extend(list(np.repeat(num, num_point_in_step)))
        num_added_point = num_point - len(stairs)
        if num_added_point < 0:
            raise RuntimeError("Negative residual count")
        elif num_added_point > 0:
            stairs.extend(list(np.repeat(num_step - 1, num_added_point)))
        stairs = np.array(stairs)
        # Rescale
        stairs = stairs*(final_value - initial_value)/(num_step - 1)
        stairs += initial_value
        #
        return stairs

    def plotStairResponse(self, times, num_step, initial_values, final_values, is_plot=True):
        """
        Plots the response to a monotonic sequence of step inputs.
        
        Parameters
        ----------
        times: list-float (times for the responses)
        num_step: int (number of steps in the stair response.
        start_values: list-float (initial values of the inputs)
        final_values: list-float (ending values of the inputs)
        is_plot: bool
        
        Returns
        -------
        Timeseries (response of outputs)
        """
        # Construct the staircase inputs
        num_point = len(times)
        stairs = []  # Base value of stairs
        for initial_value, final_value in zip(initial_values, final_values):
            new_stairs = self._makeStairs(num_point, num_step, initial_value, final_value)
            stairs.extend(new_stairs)
        # Restructure
        result_arr = np.array(stairs)
        result_arr = np.reshape(result_arr, (num_point, len(initial_values)))
        return result_arr


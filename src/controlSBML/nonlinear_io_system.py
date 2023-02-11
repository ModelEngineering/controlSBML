"""Extends control.NonlinearIOSystem for ControlSBML objects."""

import controlSBML.constants as cn
from controlSBML import logger as lg
from controlSBML import msgs
from controlSBML import util
from controlSBML.simulate_system import simulateSystem
import controlSBML.timeseries as ts
from controlSBML.option_management.option_manager import OptionManager
from controlSBML.option_management.options import Options

import control
import numpy as np
import pandas as pd
import warnings


ARG_LST = [cn.TIME, cn.STATE, cn.INPUT, cn.PARAMS]
MIN_ELAPSED_TIME = 1e-2


class NonlinearIOSystem(control.NonlinearIOSystem):

    def __init__(self, name, ctlsb,
         input_names=None, do_simulate_on_update=False):
        """
        Parameters
        ----------
        name: Name of the non-linear system created
        ctlsb: ControlSBML
        input_names: list-str (names of inputs to system)
        do_simulate_on_update: bool
            simulate to the current time before each state update

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
        if input_names is None:
            self.input_names = list(ctlsb.input_names)
        else:
            diff = set(input_names).difference(ctlsb.state_names)
            if len(diff) != 0:
                raise ValueError("Invalid input names: %s" % diff)
            self.input_names = [input_names[0]]
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

    def _makeStaircase(self, num_point, num_step, initial_value, final_value):
        """
        A staircase is a sequence of steps of the same magnitude and duration.

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
        steps = []  # Steps in the staircase
        num_point_in_step = int(num_point/num_step)
        for num in range(num_step):
            steps.extend(list(np.repeat(num, num_point_in_step)))
        num_added_point = num_point - len(steps)
        if num_added_point < 0:
            raise RuntimeError("Negative residual count")
        elif num_added_point > 0:
            steps.extend(list(np.repeat(num_step - 1, num_added_point)))
        staircase_arr = np.array(steps)
        # Rescale
        staircase_arr = staircase_arr*(final_value - initial_value)/(num_step - 1)
        staircase_arr += initial_value
        #
        return staircase_arr

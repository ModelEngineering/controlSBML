"""Extends control.NonlinearIOSystem for ControlSBML objects.
1. Inputs can be species, parameters, compartments.
2. Ouputs can be floating species, fluxes
3. If a floating species is an input, it can be manipulated in two ways:
   a. Set to a fixed value, which must be non-negative
   b. Specify an input rate, which can be positive or negative
   The choice here is indicated by the parameter is_fixed_input_species
"""

import controlSBML.constants as cn
from controlSBML import logger as lg
from controlSBML import util
from controlSBML import msgs

import control
import pandas as pd
import numpy as np


ARG_LST = [cn.TIME, cn.STATE, cn.INPUT, cn.PARAMS]
MIN_ELAPSED_TIME = 1e-2
DEFAULT_INPUT_SPECIES_DESCRIPTOR = False # True: fixed concentration, False: rate


class NonlinearIOSystem(control.NonlinearIOSystem):

    def __init__(self, name, ctlsb,
         input_names=None, output_names=None, do_simulate_on_update=False,
           is_fixed_input_species=DEFAULT_INPUT_SPECIES_DESCRIPTOR,
           is_log=False):
        """
        Parameters
        ----------
        name: Name of the non-linear system created
        ctlsb: ControlSBML
        input_names: list-str (names of inputs to system)
        do_simulate_on_update: bool
            simulate to the current time before each state update
        is_fixed_input_species: bool/dict (key: str (species name), value: bool (is fixed input species))

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
        self.is_fixed_input_species = is_fixed_input_species
        self.is_log = is_log
        if input_names is None:
            input_names = list(ctlsb.input_names)
        if output_names is None:
            output_names = list(ctlsb.output_names)
        self.input_names = input_names
        self.output_names = output_names
        self.state_names = list(ctlsb.floating_species_names)
        self.dstate_names = [self._makeDstateName(n) for n in self.state_names]
        self.num_state = len(self.state_names)
        self.num_input = len(self.input_names)
        self.num_output = len(self.output_names)
        # self.input_dct indicates if an input is held at a fixed at a value (True) or the input is a rate of change
        # FIXME: Only look at floating species
        species_names = set(ctlsb.floating_species_names).union(ctlsb.roadrunner.getBoundarySpeciesIds())
        input_species = set(self.input_names).intersection(species_names)
        self.input_dct = {n: DEFAULT_INPUT_SPECIES_DESCRIPTOR for n in self.input_names}
        if isinstance(is_fixed_input_species, bool):
            self.input_dct.update(
                  {n: is_fixed_input_species if n in self.ctlsb.floating_species_names
                    else True for n in self.input_names})
        else:
            if not isinstance(is_fixed_input_species, dict):
                msgs.error("Invalid type of is_fixed_input_species: %s" % type(is_fixed_input_species))
            self.input_dct.update((is_fixed_input_species))
        # Logging
        # Initialize the controlNonlinearIOSystem object
        item_names = list(self.input_names)
        item_names = list(self.state_names)
        out_names = list(item_names)
        self.outfcn_logger = lg.Logger(self.name, out_names)
        item_names.extend(self.dstate_names)
        self.updfcn_logger = lg.Logger(self.name, item_names)
        super().__init__(self._updfcn, self._outfcn,
              inputs=self.input_names, outputs=self.output_names,
              states=self.state_names, name=self.name)

    @property
    def outlist(self):
        return ["%s.%s" % (self.name, n) for n in self.output_names]

    def setTime(self, time):
        self.ctlsb.setTime(time)

    @staticmethod
    def _makeDstateName(name):
        return "%s'" % name

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
    
    def getSubsystem(self, name, input_names, output_names):
        """Creates a subsystem using a subset of the input and output names

        Parameters
        ----------
            input_names (_type_): _description_
            output_names (_type_): _description_
        """
        return NonlinearIOSystem(name, self.ctlsb,
              input_names=input_names, output_names=output_names,
              do_simulate_on_update=self.do_simulate_on_update,
              is_fixed_input_species=self.is_fixed_input_species 
              )
    
    def setSteadyState(self):
        """
        Sets the steady state of the system.
        
        Returns
        -------
            bool (success)
        """
        if not "roadrunner" in dir(self.ctlsb):
            return msgs.error("No roadrunner object.")
        return self.ctlsb.setSteadyState()

    def _updfcn(self, time, x_vec, u_vec, _):
        """
        Computes the change in state by controlling a species to a particular concentration.
        We don't do addition to or subtraction from a species since it's difficult to calculate
        the impact on rates.
        Roadrunner computes the derivatives.
        No simulation is run, and so this technique
        may not always work.

        Parameters
        ----------
        time: float: time
        x_vec: np.array(float): state vector (floating species)
        u_vec: np.array(float): input vector (in log10 units)

        Returns
        -------
        np.array(float): change in state vector
        """
        if isinstance(u_vec, float):
            u_vec = [u_vec]
        if self.do_simulate_on_update:
            self.setTime(time)  # Consider time dependent functions
        # Ensure that the state of the simulation is correct
        state_dct = {n: np.max(x_vec[i], 0)
              for i, n in enumerate(self.state_names)}
        self.ctlsb.set(state_dct)
        # Set the values of the inputs
        input_dct = {n: u_vec[i] for i, n in enumerate(self.input_names)}
        set_dct = dict(state_dct)
        # Handle fixed concentration species
        for name, is_fixed in self.input_dct.items():
            if is_fixed:
                set_dct[name] = input_dct[name]
        self.ctlsb.set(set_dct)
        # Calculate the derivatives of floating species in state
        derivative_dct = {n: v for n, v in self.ctlsb.get(self.dstate_names).items()}
        # Handle input species for which an additional rate is specified
        for name, is_fixed in self.input_dct.items():
            if not is_fixed:
                dname = self._makeDstateName(name)
                derivative_dct[dname] += input_dct[name]
        derivative_arr = np.array(list(derivative_dct.values()))
        # Update logger
        if self.is_log:
            arr = np.append(x_vec, u_vec)
            arr = np.append(arr, derivative_arr)
            self.updfcn_logger.add(time, arr)
        return derivative_arr

    def _outfcn(self, time, x_vec, u_vec, __):
        """
        Calculates the values of outputs.
        Ensures that all state variables are non-negative since they are
        rates or concentrations.

        Parameters
        ----------
        x_vec: np.array(float): state vector
        u_vec: np.array(float): input vector (in log10 units)

        Returns
        -------
        np.array(float): change in state vector
        """
        # Update the state of the simulation
        self.ctlsb.setTime(time)
        state_dct = {self.state_names[n]: x_vec[n] for n in range(len(self.state_names))}
        self.ctlsb.set(state_dct)
        # Calculate the values of the outputs
        out_vec = np.repeat(np.nan, self.num_output)
        for out_idx, name in enumerate(self.output_names):
            out_vec[out_idx] = self.ctlsb.get(name)
            if name in self.ctlsb.floating_species_names:
                outs = list(out_vec)
                outs[out_idx] = np.abs(out_vec[out_idx])
                out_vec = np.array(outs)
        if np.isnan(np.sum(out_vec)):
            raise ValueError("Outputs could not be calculated.")
        # Update logger
        if self.is_log:
            arr = np.append(x_vec, u_vec)
            self.outfcn_logger.add(time, arr)
        return out_vec
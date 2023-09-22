"""Creates a system from an SBML model: inputs, input type, outputs."""

import controlSBML.constants as cn
from controlSBML import msgs
from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML.timeseries import Timeseries
from controlSBML import util

import numpy as np
import pandas as pd


class SBMLSystem(object):

    def __init__(self, model_reference, input_names, output_names, is_fixed_input_species=False):
        """
        model_reference: str
            string, SBML file or Roadrunner object
        input_names: name (id) of reaction whose flux is being controlled
                     or the name of a chemical species
        output_names: list-str
            output species
        is_fixed_input_species: bool (input species are fixed)
        """
        # First initializations
        self.model_reference = model_reference
        self.is_fixed_input_species = is_fixed_input_species
        self.roadrunner = makeRoadrunner(self.model_reference)
        self.floating_species_names = list(set(self.roadrunner.getFloatingSpeciesIds()))
        # Handle inputs/outputs
        self.reaction_names = list(self.full_stoichiometry_df.columns)
        self.num_input = len(self.input_names)
        self.num_output = len(self.output_names)
        # Other calculations
        self.antimony = self.roadrunner.getAntimony()
        # Validation checks
        invalid_names = self._invalidInputNames()
        if len(invalid_names) > 0:
            text = "Inputs must be a species, parameter, or compartment."
            text += "   Invalid names are: %s" % str(invalid_names)
            raise ValueError(text)
        invalid_names = self._invalidOutputNames()
        if len(invalid_names) > 0:
            text = "Outputs must be species or fluxes."
            text += "The following outputs are invalid: %s" % str(invalid_names)
            raise ValueError(text)
        
    def _isExistingName(self, name, names):
        if name in self.roadrunner.keys():
            return True
        names.append(name)
        return False
        
    def _invalidInputNames(self):
        """
        Verifies that the input names are valid.

        Returns
        -------
        bool
        """
        invalid_names = []
        for name in self.input_names:
            if not self._isExistingName(name, invalid_names):
                continue
            if name in self.reaction_names:
                invalid_names.append(name)
                continue
            try:
                value = self.get(name)
                self.set({name: value})
            except Exception:
                invalid_names.append(name)
        return invalid_names
    
    def _invalidOutputNames(self):
        """
        Verifies that the output names are valid.

        Returns
        -------
        bool
        """
        invalid_names = []
        for name in self.output_names:
            if not self._isExistingName(name, invalid_names):
                continue
            try:
                # Must be able to retrieve a value
                _ = self.get(name)
            except RuntimeError:
                invalid_names.append(name)
        return invalid_names

    def get(self, names=None):
        """
        Provides the roadrunner values for a name. If no name,
        then all values of symbols are given.

        Parameters
        ----------
        name: str/list-str

        Returns
        -------
        object/dict
        """
        if names is None:
            names = [k for k in self.roadrunner.keys() if not ")" in k]
        if isinstance(names, str):
            return self.roadrunner[names]
        return {n: self.roadrunner[n] for n in names}

    def set(self, name_dct):
        """
        Sets the values of names and values.

        Parameters
        ----------
        name_dct: dict
            key: str
            value: value
        """
        util.setRoadrunnerValue(self.roadrunner, name_dct)
    
    def setSteadyState(self):
        """
        Sets the steady state of the system.
        
        Returns
        -------
            bool (success)
        """
        # Try to find the steady state
        for _ in range(3):
            try:
                self.roadrunner.steadyState()
                return True
            except RuntimeError:
                pass
        return False
    
    def simulate(self, start_time=cn.START_TIME, end_time=cn.END_TIME, num_point=None):
        """
        Simulates the system.

        Parameters
        ----------
        start_time: float
        end_time: float
        num_point: int

        Returns
        -------
        DataFrame
        """
        if num_point is None:
            num_point = int(cn.POINTS_PER_TIME*(end_time - start_time))
        self.roadrunner.reset()
        data = self.roadrunner.simulate(start_time, end_time, num_point)
        column_names = [c[1:-1] if c[0] == "[" else c for c in data.colnames]
        df = pd.DataFrame(data, columns=column_names)
        ts = Timeseries(df, times=df[cn.TIME])
        return ts
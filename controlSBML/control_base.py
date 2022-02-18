"""LTI Control for SBML models. Base handles initialization, get, set."""

"""
TO DO:

1. Plot difference between time jacoabian at reference vs. Current.
2. Plot TOTAL residual SSQ vs. jacobian difference
"""

from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML.options import Options

import control
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tellurium as te


TIME = "time"
START_TIME = 0  # Default start time
END_TIME = 5  # Default endtime
POINTS_PER_TIME = 10
TIME = "time"


class ControlBase(object):

    def __init__(self, model_reference):
        """
        Initializes instance variables
        :param str model_reference: string or SBML file or Roadrunner object
        """
        ##### PUBLIC #####
        self.model_reference = model_reference
        self.roadrunner = makeRoadrunner(model_reference)
        self.antimony = self.roadrunner.getAntimony()
        # Do the initializations
        self.boundary_species = self.roadrunner.getBoundarySpeciesIds()
        self._species_names = None
        self.roadrunner.reset()

    def _mkBoundarySpeciesFloating(self):
        for name in self.boundary_species:
            self.roadrunner.setBoundary(name, False)

    def _unmkBoundarySpeciesFloating(self):
        for name in self.boundary_species:
            self.roadrunner.setBoundary(name, True)

    @property
    def species_names(self):
        if self._species_names is None:
            self._mkBoundarySpeciesFloating()
            try:
                mat = self.roadrunner.getFullJacobian()
                self._unmkBoundarySpeciesFloating()
            except Exception:
                self._unmkBoundarySpeciesFloating()
                mat = self.roadrunner.getFullJacobian()
            self._species_names = list(mat.colnames)
        return self._species_names

    @property
    def jacobian(self):
        """
        Returns
        -------
        NamedArray with names for rows (rownames) and columns (colnames)
        Handle boundary species.
        """
        self._mkBoundarySpeciesFloating()
        try:
            mat = self.roadrunner.getFullJacobian()
            self._unmkBoundarySpeciesFloating()
        except Exception:
            self._unmkBoundarySpeciesFloating()
            mat = self.roadrunner.getFullJacobian()
        for idx, name in enumerate(mat.rownames):
            if name in self.boundary_species:
                mat[idx, :] = 0
        df = pd.DataFrame(mat, columns=mat.colnames, index=self.species_names)
        return df

    @property
    def current_state(self):
        """
        Contructs vector of current state values (floating and boundary species)

        Returns
        -------
        Series
            index: str (state names)
        """
        values = list(self.get(self.species_names).values())
        return pd.Series(values, index=self.species_names)

    @staticmethod
    def isRoadrunnerKey(key):
        return not ((key[0] == "_") or ("(" in key) or (key[-1] == "'"))

    def setTime(self, time):
        self.roadrunner.reset()
        _ = self.roadrunner.simulate(0, time)

    def copy(self):
        """
        Creates a copy of the object.

        Returns
        -------
        controlSBML
        """
        ctlsb = self.__class__(self.model_reference)
        # Update roadrunner
        for key, value in self.roadrunner.items():
            if self.isRoadrunnerKey(key):
                ctlsb.roadrunner[key] = value
        return ctlsb

    def equals(self, other):
        """
        Checks that they have the same information

        Parameters
        ----------
        other: ControlSBML

        Returns
        -------
        bool
        """
        bValue = self.antimony == other.antimony
        bValue = bValue and all([s1 == s2 for s1, s2
              in zip(self.species_names, self.jacobian.columns)])
        diff = set(self.roadrunner.keys()).symmetric_difference(
              other.roadrunner.keys())
        bValue = bValue and (len(diff) == 0)
        # Check the roadrunner state
        if bValue:
            for key, value in self.roadrunner.items():
                if self.isRoadrunnerKey(key):
                    bValue = bValue and (other.roadrunner[key] == value)
        return bValue

    def get(self, names=None):
        """
        Provides the roadrunner values for a name. If no name,
        then all values are given.

        Parameters
        ----------
        name: str/list-str

        Returns
        -------
        object/dict
        """
        if isinstance(names, str):
            return self.roadrunner[names]
        if names is None:
            names = self.roadrunner.keys()
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
        for name, value in name_dct.items():
            self.roadrunner[name] = value

    def makeStateSpace(self, A=None, B=None, C=None, D=None):
        """
        Creates a control system object for the n X n jacobian.

        Parameters
        ----------
        A: np.array(n X n)
        B: np.array(n X p)
        C: np.array(q X n)
        D: np.array(q X p)

        Returns
        -------
        control.StateSpace
        """
        # Construct the matrices
        if A is None:
            A = self.jacobian.values
        if B is None:
            B = np.repeat(0, A.shape[0])
            B = np.reshape(B, (A.shape[0], 1))
        if C is None:
            C = np.identity(A.shape[0])
        if D is None:
            D = B
        return control.StateSpace(A, B, C, D)

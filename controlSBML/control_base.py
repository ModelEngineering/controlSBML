"""LTI Control for SBML models. Base handles initialization, get, set."""

"""
TO DO:

1. Plot difference between time jacoabian at reference vs. Current.
2. Plot TOTAL residual SSQ vs. jacobian difference
"""

from controlSBML.make_roadrunner import makeRoadrunner

import control
import numpy as np
import pandas as pd


TIME = "time"
START_TIME = 0  # Default start time
END_TIME = 5  # Default endtime
POINTS_PER_TIME = 10
TIME = "time"


class ControlBase(object):

    def __init__(self, model_reference, include_boundary_species=True):
        """
        Initializes instance variables
        :param str model_reference: string or SBML file or Roadrunner object
        """
        ##### PUBLIC #####
        self.model_reference = model_reference
        self.include_boundary_species = include_boundary_species
        self.roadrunner = makeRoadrunner(model_reference)
        self.antimony = self.roadrunner.getAntimony()
        # Do the initializations
        self.boundary_species = self.roadrunner.getBoundarySpeciesIds()
        self.roadrunner.reset()

    def _mkBoundarySpeciesFloating(self):
        if self.include_boundary_species:
            for name in self.boundary_species:
                self.roadrunner.setBoundary(name, False)

    def _unmkBoundarySpeciesFloating(self):
        if self.include_boundary_species:
            for name in self.boundary_species:
                self.roadrunner.setBoundary(name, True)

    @property
    def jacobian(self):
        """
        Returns
        -------
        pd.DataFrame
        """
        return self._getJacobian(is_reduced=False)

    @property
    def reduced_jacobian(self):
        """
        Returns
        -------
        pd.DataFrame
        """
        return self._getJacobian(is_reduced=True)

    def _getJacobian(self, is_reduced=True):
        """
        Parameters
        ----------
        is_reduced: bool
            return a reduced Jacobian that eliminates
            conserved moeities and zero rows/columns

        Returns
        -------
        pd.DataFrame, species_names
        """
        def getMat():
            if is_reduced:
                mat = self.roadrunner.getReducedJacobian()
                all_names = list(mat.rownames)
                idxs = []
                for idx, name in enumerate(all_names):
                   sum_of_row_squares = np.sum(mat[idx, :]**2)
                   sum_of_column_squares = np.sum(mat[:, idx]**2)
                   if (not np.isclose(sum_of_row_squares, 0))  \
                       and (not np.isclose(sum_of_column_squares, 0)):
                       idxs.append(idx)
                names = [all_names[i] for i in idxs]
                jacobian_mat = np.array([[mat[i, j] for j in idxs] for i in idxs])
                jacobian_mat = np.reshape(jacobian_mat, (len(idxs), len(idxs)))
            else:
                jacobian_mat = self.roadrunner.getFullJacobian()
                names = list(jacobian_mat.rownames)
            return jacobian_mat, names
        #
        try:
            if self.include_boundary_species:
                self._mkBoundarySpeciesFloating()
                mat, names = getMat()
                self._unmkBoundarySpeciesFloating()
            else:
                mat, names = getMat()
        except Exception:
            if self.include_boundary_species:
                self._unmkBoundarySpeciesFloating()
            mat, names = getMat()
        for idx, name in enumerate(names):
            if name in self.boundary_species:
                mat[idx, :] = 0
        try:
            df = pd.DataFrame(mat, columns=names, index=names)
        except:
            import pdb; pdb.set_trace()
        return df

    def getSpeciesNames(self, is_reduced=False):
        """
        Gets the list of species names.

        Parameters
        ----------
        is_reduced: bool
            create a reduced model
        
        Returns
        -------
        list-str
        """
        if is_reduced:
            df = self.reduced_jacobian
        else:
            df = self.jacobian
        return list(df.columns)

    def getCurrentState(self, is_reduced=False, species_names=None):
        """
        Contructs vector of current state values (floating and boundary species)

        Parameters
        ----------
        is_reduced: bool
            Is a reduced model
        species_names: list-str
            List of species in state

        Returns
        -------
        Series
            index: str (state names)
        """
        if species_names is None:
            species_names = self.getSpeciesNames(is_reduced=is_reduced)
        values = list(self.get(species_names).values())
        return pd.Series(values, index=species_names)

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
              in zip(self.getSpeciesNames(), other.getSpeciesNames())])
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

    def makeStateSpace(self, A_mat=None, B_mat=None, C_mat=None, D_mat=None,
          is_reduced=True):
        """
        Creates a control system object for the n X n jacobian.

        Parameters
        ----------
        A_mat: np.array(n X n) or DataFrame
        B_mat: np.array(n X p) or DataFrame
        C_mat: np.array(q X n) or DataFrame
        D_mat: np.array(q X p) or DataFrame

        Returns
        -------
        control.StateSpace
        """
        def df2Mat(df):
            if isinstance(df, pd.DataFrame):
                return df.values
            else:
                return df
        #
        # Construct the matrices
        A_mat = df2Mat(A_mat)
        B_mat = df2Mat(B_mat)
        C_mat = df2Mat(C_mat)
        D_mat = df2Mat(D_mat)
        if A_mat is None:
            if is_reduced:
                A_mat = self.reduced_jacobian.values
            else:
                A_mat = self.jacobian.values
        if B_mat is None:
            B_mat = np.repeat(0, A_mat.shape[0])
            B_mat = np.reshape(B_mat, (A_mat.shape[0], 1))
        if C_mat is None:
            C_mat = np.identity(A_mat.shape[0])
        if D_mat is None:
            nrow = np.shape(C_mat)[0]
            ncol = np.shape(B_mat)[1]
            D_mat = np.repeat(0, nrow*ncol)
            D_mat = np.reshape(D_mat, (nrow, ncol))
        return control.StateSpace(A_mat, B_mat, C_mat, D_mat)

"""LTI Control for SBML models. Base handles initialization, get, set."""

"""
This module creates an LTI system for an SBML model.
States are chemical species. Inputs are unnormalized enzyme reaction elasticities.
Outputs are chemical species.

Notes:
1. Reaction enzymes are identified by the SBML reaction ID.
"""

"""
TO DO:

1. Plot difference between time jacoabian at reference vs. Current.
2. Plot TOTAL residual SSQ vs. jacobian difference
"""

from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML import msgs
import controlSBML as ctl

import control
import numpy as np
import pandas as pd


TIME = "time"
START_TIME = 0  # Default start time
END_TIME = 5  # Default endtime
POINTS_PER_TIME = 10
TIME = "time"
IS_DEBUG = False


cleanIt = lambda n: n if n[0] != "[" else n[1:-1]


class ControlBase(object):

    def __init__(self, model_reference, input_names=None, output_names=None):
        """
        Initializes instance variables
        model_reference: str
            string, SBML file or Roadrunner object
        input_names: name (id) of reaction whose flux is being controlled
        output_names: list-str
            output species
        """
        def calcNames(superset_names, subset_names, superset_name, subset_name):
            if subset_names is None:
                return list(superset_names)
            else:
                if not set(subset_names) <= set(superset_names):
                    text = "%s %s are not a subset of the %s %s" % \
                          (subset_name, subset_names, superset_name,
                          superset_names)
                    raise ValueError(text)
                return subset_names
        # Iinitial model calculations
        self.model_reference = model_reference
        self.roadrunner = makeRoadrunner(self.model_reference)
        # Set defaults
        self.state_names = list(self.roadrunner.getFullJacobian().rownames)
        self.species_names = list(
              self.roadrunner.getFullStoichiometryMatrix().rownames)
        self.depeendent_names = list(
              set(self.species_names).symmetric_difference(self.state_names))
        if not set(self.state_names) <= set(self.species_names):
            import pdb; pdb.set_trace()
            raise RuntimeError("State name is not a species name.")
        self.state_names = self._sortList(self.species_names, self.state_names)
        self.full_stoichiometry_df, self.reduced_stoichiometry_df  \
              = self._makeStoichiometryDF()
        # Check for consistency on the state specification
        self.reaction_names = list(self.reduced_stoichiometry_df.columns)
        self.num_state = len(self.state_names)
        # Handle defaults
        if input_names is None:
            self.input_names = []
        else:
            self.input_names = input_names
        self.input_names = self._sortList(self.reaction_names, self.input_names)
        self.num_input = len(self.input_names)
        self.output_names = calcNames(self.state_names, output_names,
              "states", "outputs")
        self.output_names = self._sortList(self.species_names, self.output_names)
        self.num_output = len(self.output_names)
        # Other calculations
        self.B_df = self._makeBDF()
        self.C_df = self._makeCDF()
        self.state_names = list(self.jacobian_df.columns)
        self.antimony = self.roadrunner.getAntimony()
        # Do the initializations
        self.roadrunner.reset()
        # Validation checks
        if set(self.state_names) != set(self.species_names):
            text = "State does not include some spaces.\n"
            text += "  Species are: %s" % str(self.species_names)
            text += "  States are: %s" % str(self.state_names)
            raise RuntimeError(text)
        if not set(self.output_names) <= set(self.species_names):
            diff = list(set(self.output_names).difference(self.species_names))
            text = "Outputs must be species. The following outputs are not species"
            text += "The following outputs are not species: %s" % str(diff)
            raise ValueError(text)

    @staticmethod
    def _makeDF(mat):
        df = pd.DataFrame(mat, columns=mat.colnames, index=mat.rownames)
        return df

    def _makeStoichiometryDF(self):
        """
        Creates the reduced stoichiometry matrix and the auxiliary matrix

        Returns
        -------
        DataFrame - full stoichiometry matrix
        DataFrame - reduced stoichiometry matrix
        """
        #
        reduced_stoichiometry_df = self._makeDF(
              self.roadrunner.getReducedStoichiometryMatrix())
        full_stoichiometry_df = self._makeDF(
              self.roadrunner.getFullStoichiometryMatrix())
        return full_stoichiometry_df, reduced_stoichiometry_df

    def _makeCDF(self):
        """
        Creates the output C dataframe based on the requested output_names.
        Columns should be the states. Rows are the outputs.

        Returns
        -------
        pd.DataFrame
        """
        mat = np.identity(len(self.state_names))
        C_df = pd.DataFrame(mat, columns=self.species_names,
              index=self.species_names)
        # Adjust to outputs
        for name in self.species_names:
            if not name in self.output_names:
                C_df = C_df.drop(name)
        # TODO: Address conserved moieties
        if False:
            # Handle computed outputs
            if self.roadrunner.getNumConservedMoieties() == 0:
                if not set(self.output_names) <= set(self.state_names):
                    import pdb; pdb.set_trace()
                    raise RuntimeError("Cannot compute output from states.")
            else:
                # Add the conservation matrix
                mat = self.roadrunner.getConservationMat()
                colnames = list(mat.colnames)
                if set(colnames) != set(self.state_names):
                    import pdb; pdb.set_trace()
                    raise RuntimeError(
                          "Not all states are columns of the conservation matrix.")
                rownames = list(mat.rownames)
                df = pd.DataFrame(self.roadrunner.getConservationMat(),
                       columns=colnames, index=rownames)
                # Remove overlapping rows
                for rowname in rownames:
                    if rowname in self.output_names:
                        C_df[rowname] = df[rowname]
        return C_df

    @property
    def roadrunner_namespace(self):
        """
        Constructs the roadrunner namespace and associated values.

        Parameters
        ----------

        Returns
        -------
        dict
        """
        dct = {}
        for id_lst in  [self.roadrunner.getCompartmentIds(),
               self.roadrunner.getBoundarySpeciesIds(),
               self.roadrunner.getFloatingSpeciesIds(),
               self.roadrunner.getBoundarySpeciesIds(),
               self.roadrunner.getGlobalParameterIds(),
               ]:
            for id in id_lst:
                 dct[id] = self.get(id)
        return dct

    @property
    def state_ser(self):
        """
        Contructs vector of current state values.

        Returns
        -------
        np.ndarray: N X 1
        """
        return self._makeSer(self.state_names)

    @property
    def output_ser(self):
        """
        Contructs vector of current state values.

        Returns
        -------
        np.ndarray: N X 1
        """
        return self._makeSer(self.output_names)

    def _makeSer(self, names):
        """
        Contructs a Series for the names.

        Returns
        -------
        pd.Series
        """
        return pd.Series(list(self.get(names).values()), index=names)

    @property
    def jacobian_df(self):
        """
        Returns
        -------
        pd.DataFrame, species_names
        """
        # FIXME: Use reduced stoichiometry matrix once it's clear how
        #        to calculate dependent spcies
        jacobian_mat = self.roadrunner.getFullJacobian()
        if len(jacobian_mat.rownames) != len(jacobian_mat.colnames):
            raise RuntimeError("Jacobian is not square!")
        names = list(jacobian_mat.colnames)
        jacobian_df = pd.DataFrame(jacobian_mat, columns=names, index=names)
        jacobian_df = jacobian_df.loc[self.state_names, self.state_names]
        return jacobian_df

    @property
    def A_df(self):
        return self.jacobian_df

    @staticmethod
    def isRoadrunnerKey(key):
        return not ((key[0] == "_") or ("(" in key) or (key[-1] == "'"))

    def setTime(self, time):
        self.roadrunner.reset()
        if time > 0.01:
            _ = self.roadrunner.simulate(0.0, time)

    # FIXME: Doesn't update "sets" done to roadrunner
    def copy(self):
        """
        Creates a copy of the object.

        Returns
        -------
        controlSBML
        """
        ctlsb = self.__class__(self.model_reference,
              input_names=self.input_names,
              output_names=self.output_names)
        ctlsb.setTime(self.getTime())
        return ctlsb

    def getTime(self):
        """
        Gets current simulation time.
        
        Returns
        -------
        float
        """
        return self.roadrunner.model.getTime()

    # TODO: More complete check of attributes?
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
        if IS_DEBUG:
             print("1: %d" % bValue)
        bValue = bValue and np.isclose(self.getTime(),  \
              other.getTime())
        if IS_DEBUG:
             print("2: %d" % bValue)
        bValue = bValue and all([s1 == s2 for s1, s2
              in zip(self.state_names, other.state_names)])
        if IS_DEBUG:
             print("3: %d" % bValue)
        diff = set(self.roadrunner.keys()).symmetric_difference(
              other.roadrunner.keys())
        bValue = bValue and (len(diff) == 0)
        if IS_DEBUG:
             print("4: %d" % bValue)
        for attr in ["state_names", "input_names", "output_names"]:
            expr1 = "self.%s" % attr
            expr2 = "other.%s" % attr
            try:
                np.array(eval(expr1)) == np.array(eval(expr2))
            except Exception:
                bValue = False
                break
        if IS_DEBUG:
             print("5: %d" % bValue)
        # Check the roadrunner state
        if bValue:
            for key, value in self.roadrunner.items():
                if self.isRoadrunnerKey(key):
                    bValue = bValue and (other.roadrunner[key] == value)
        if IS_DEBUG:
             print("6: %d" % bValue)
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
            if isinstance(value, int):
                value = float(value)
            self.roadrunner[name] = value

    @staticmethod
    def _sortList(super_lst, sub_lst):
        """
        Sorts the sub_lst in the same order as the super_lst.

        Parameters
        ----------
        super_lst: list
        sub_lst: list

        Returns
        -------
        list
        """
        new_super_lst = list(super_lst)
        return sorted(sub_lst, key=lambda v: new_super_lst.index(v))

    def _makeBDF(self):
        """
        Constructs a dataframe for the B matrix.

        Returns
        -------
        np.ndarray (n X p), where p = len(input_names)
        """
        # Select the columns needed from the stoichiometry matrix
        B_df = self.full_stoichiometry_df[self.input_names]
        # Subset to the states
        sub_names = list(set(B_df.index).intersection(self.state_names))
        sub_names = self._sortList(self.state_names, sub_names)
        B_df = B_df.loc[sub_names, :]
        #
        return B_df

    def makeStateSpace(self, A_mat=None, B_mat=None,
          C_mat=None, D_mat=None):
        """
        Creates a state space control object for
        the n X n jacobian. By default, the D matrix is always 0.

        Parameters
        ----------
        The default values of the matrices are calculated in the constructor.
        These can be overridden.
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
            A_mat = self.jacobian_df.values
        if B_mat is None:
            B_mat = self.B_df.values
        if C_mat is None:
            # Construct the output matrix
            C_mat = self.C_df.values
        if D_mat is None:
            nrow = len(self.output_names)
            ncol = np.shape(B_mat)[1]
            D_mat = np.repeat(0, nrow*ncol)
            D_mat = np.reshape(D_mat, (nrow, ncol))
        return control.StateSpace(A_mat, B_mat, C_mat, D_mat)

    def makeNonlinearIOSystem(self, name, effector_dct=None):
        """
        Creates an object that can be used in connections with the
        control package.
 
        Parameters
        ----------
        name: str (name of the system)
        effector_dct: dict (maps reaction inputs to roadrunner muteables)
            key: str (input name)
            value: str (name of roadrunner muteable)
        
        Returns
        -------
        controlSBML.NonelinearIOSystem
        """
        return ctl.NonlinearIOSystem(name, self, effector_dct=effector_dct)

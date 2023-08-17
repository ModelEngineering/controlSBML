"""LTI Control for SBML models. Base handles initialization, get, set."""

"""
This module creates an LTI system for an SBML model.
States are chemical species. Inputs are unnormalized enzyme reaction elasticities.
Outputs are chemical species.

Notes:
1. Reaction enzymes are identified by the SBML reaction ID.
2. The Jacobian is only recalculated if there is a change in time
"""


from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML import util
from controlSBML import msgs

import control
import numpy as np
import pandas as pd


ATOL = 1e-8  # Absolute tolerance for comparisons
TIME = "time"
START_TIME = 0  # Default start time
END_TIME = 5  # Default endtime
POINTS_PER_TIME = 10
TIME = "time"
IS_DEBUG = False
TIMEPOINT_NULL = 1


cleanIt = lambda n: n if n[0] != "[" else n[1:-1]


class ControlBase(object):

    def __init__(self, model_reference, input_names=None, output_names=None,
          is_reduced=False):
        """
        Initializes instance variables
        model_reference: str
            string, SBML file or Roadrunner object
        input_names: name (id) of reaction whose flux is being controlled
                     or the name of a chemical species
        output_names: list-str
            output species
        is_reduced: bool
            construct a reduced model so that the A matrix is nonsingular
        """
        # First initializations
        self.model_reference = model_reference
        self.roadrunner = makeRoadrunner(self.model_reference)
        self._default_state_names = list(
              self.roadrunner.getFullStoichiometryMatrix().rownames)
        # FIXME: Handle State Space models better so can allow BoundarySpecies
        if False:
            self.species_names = set(self.roadrunner.getFloatingSpeciesIds())
            self.species_names = list(
                  self.species_names.union(self.roadrunner.getBoundarySpeciesIds()))
        else:
            self.species_names = list(
                  self.roadrunner.getFullStoichiometryMatrix().rownames)
        if output_names is None:
            output_names = self.species_names
        self.output_names = output_names
        # Initializations with error checks
        self.is_reduced = is_reduced
        if len(set(self.roadrunner.getReactionIds()).intersection(self.output_names)) > 0:
            if self.is_reduced:
                self.is_reduced = False
                text = "Cannot have a flux output and is_reduced=True."
                text += "  Setting is_reduced=False."
                msgs.warn(text)
        # Iinitial model calculations
        self._jacobian_time = TIMEPOINT_NULL
        self._jacobian_df = None
        # Set defaults.
        self.depeendent_names = list(
              set(self.species_names).symmetric_difference(self.state_names))
        if not set(self.state_names) <= set(self.species_names):
            raise RuntimeError("State name is not a species name.")
        self.full_stoichiometry_df, self.reduced_stoichiometry_df  \
              = self._makeStoichiometryDF()
        # Check for consistency on the state specification
        if self.is_reduced:
            self.reaction_names = list(self.reduced_stoichiometry_df.columns)
        else:
            self.reaction_names = list(self.reduced_stoichiometry_df.columns)
        # Handle defaults
        if input_names is None:
            self.input_names = self.species_names
        else:
            self.input_names = input_names
        #self.input_names = self._sortList(self.reaction_names, self.input_names)
        self.num_input = len(self.input_names)
        self.num_output = len(self.output_names)
        # Other calculations
        self.antimony = self.roadrunner.getAntimony()
        # Do the initializations
        self.roadrunner.reset()
        # Validation checks
        if not set(self.state_names) <= set(self.species_names):
            text = "State does not include some spaces.\n"
            text += "  Species are: %s" % str(self.species_names)
            text += "  States are: %s" % str(self.state_names)
            raise RuntimeError(text)
        #
        possible_names = set(self.species_names).union(self.reaction_names)
        if not set(self.output_names) <= set(possible_names):
            diff = list(set(self.output_names).difference(self.species_names))
            text = "Outputs must be species or fluxes."
            text += "The following outputs are invalid: %s" % str(diff)
            raise ValueError(text)
        #
        possible_names = set(self.species_names).union(self.roadrunner.getBoundarySpeciesIds())
        if not set(self.input_names) <= set(possible_names):
            diff = list(set(self.input_names).difference(possible_names))
            text = "Inputs must be a species."
            text += "   Invalid names are: %s" % str(diff)
            raise ValueError(text)

    def _makeStoichiometryDF(self):
        """
        Creates the reduced stoichiometry matrix and the auxiliary matrix

        Returns
        -------
        DataFrame - full stoichiometry matrix
        DataFrame - reduced stoichiometry matrix
        """
        #
        reduced_stoichiometry_df = util.mat2DF(
              self.roadrunner.getReducedStoichiometryMatrix())
        full_stoichiometry_df = util.mat2DF(
              self.roadrunner.getFullStoichiometryMatrix())
        return full_stoichiometry_df, reduced_stoichiometry_df

    @property
    def roadrunner_namespace(self):
        """
        Constructs the roadrunner namespace and associated values.

        Parameters
        ----------
:
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
            for idx in id_lst:
                 dct[idx] = self.get(idx)
        return dct

    @property
    def state_ser(self):
        """
        Contructs Series of current state values.

        Returns
        -------
        pd.Series
        """
        ser = util.makeRoadrunnerSer(self.roadrunner, self.state_names)
        return ser

    @property
    def output_ser(self):
        """
        Contructs vector of current state values.

        Returns
        -------
        np.ndarray: N X 1
        """
        return util.makeSer(self.roadrunner, self.output_names)

    @property
    def jacobian_df(self):
        """
        Calculates the Jacobian, or a reduced Jacobian if the option is selected.
        Improves efficiency by using a previous calculation if the time has not changed.

        Returns
        -------
        pd.DataFrame, species_names
        """
        if np.isclose(self.getTime(), self._jacobian_time):
            if self._jacobian_df is not None:
                return self._jacobian_df
        if self.is_reduced:
            current_bool = self.roadrunner.conservedMoietyAnalysis
            self.roadrunner.conservedMoietyAnalysis = True
            jacobian_mat = self.roadrunner.getReducedJacobian()
            self.roadrunner.conservedMoietyAnalysis = current_bool
        else:
            try:
                jacobian_mat = self.roadrunner.getFullJacobian()
            except RuntimeError as excp:
                msgs.warn("Cannot calculate Jacobian because %s" % excp)
                return None
        if len(jacobian_mat.rownames) != len(jacobian_mat.colnames):
            raise RuntimeError("Jacobian is not square!")
        names = list(jacobian_mat.colnames)
        self._jacobian_df = pd.DataFrame(jacobian_mat, columns=names, index=names)
        self._jacobian_time = self.getTime()
        return self._jacobian_df

    @property
    def state_names(self):
        jacobian_df = self.jacobian_df
        if jacobian_df is not None:
            state_names = list(jacobian_df.columns)
            return self._sortList(self.species_names, state_names)
        else:
            return self._default_state_names

    @property
    def num_state(self):
        return len(self.state_names)

    @property
    def A_df(self):
        return self.jacobian_df

    @staticmethod
    def isRoadrunnerKey(key):
        return not ((key[0] == "_") or ("(" in key) or (key[-1] == "'"))

    def getJacobian(self, time=None):
        """
        Calculates the Jacobian at the specified time.

        Parameters
        ----------
        time: float

        Returns
        -------
        pd.DataFrame
        """
        # Calculate the Jacobian
        if time is not None:
            current_time = self.getTime()
            self.setTime(time)
        jacobian_df = self.jacobian_df.copy()
        if jacobian_df is None:
            raise ValueError("Cannot calculate Jacobian")
        if time is not None:
            self.setTime(current_time)
        return jacobian_df

    def setTime(self, time):
        self.roadrunner.reset()
        self._jacobian_time = TIMEPOINT_NULL
        if time > 0.01:
            _ = self.roadrunner.simulate(0.0, time)

    # FIXME: Doesn't update "sets" done to roadrunner. Can resolve by
    #        by assigning sets to new instance.
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

    def equals(self, other, is_quick_check=False):
        """
        Checks that they have the same information

        Parameters
        ----------
        other: ControlSBML
        is_quick_check: bool (don't do detailed comparisons)

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
        if not is_quick_check:
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


    def add(self, name_dct):
        """
        Adds the indicated value to the current value of the variable.

        Parameters
        ----------
        name_dct: dict
            key: str
            value: value
        """
        cur_dct = util.getRoadrunnerValue(self.roadrunner, name_dct.keys())
        new_dct = {n: cur_dct[n] + name_dct[n] for n in name_dct.keys()}
        util.setRoadrunnerValue(self.roadrunner, new_dct)

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

    def separateSpeciesReactionInputs(self):
        species_inputs = [n for n in self.input_names if n in self.species_names]
        reaction_inputs = [n for n in self.input_names if n in self.reaction_names]
        return species_inputs, reaction_inputs
    
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
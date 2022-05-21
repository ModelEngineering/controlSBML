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
import controlSBML as ctl
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
            self.input_names = []
        else:
            self.input_names = input_names
        #self.input_names = self._sortList(self.reaction_names, self.input_names)
        self.num_input = len(self.input_names)
        # FIXME: Is this necessary?
        # self.output_names = self._sortList(self.species_names, self.output_names)
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
        possible_names = set(self.species_names).union(self.reaction_names)
        if not set(self.output_names) <= set(possible_names):
            diff = list(set(self.output_names).difference(self.species_names))
            text = "Outputs must be species or fluxes."
            text += "The following outputs are invalid: %s" % str(diff)
            raise ValueError(text)
        possible_names = set(self.species_names).union(self.reaction_names)
        if not set(self.input_names) <= set(possible_names):
            diff = list(set(self.input_names).difference(possible_names))
            text = "Inputs must be a species or a reaction."
            text += "   Invalid names are: %s" % str(diff)
            raise ValueError(text)

    @property
    def B_df(self):
        return self._makeBDF()

    @property
    def C_df(self):
        return self._makeCDF()

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

    def _makeCDF(self):
        """
        Creates the output C dataframe based on the requested output_names.
        Columns should be the states. Rows are the outputs.
        There are 3 cases for outputs:
         1) The output is a state
         2) The output is a floating species whose concentration
            is a linear function of the other state variables
         3) The output is a flux and this
            is a linear function of the other state variables

        Returns
        -------
        pd.DataFrame
        """
        # FIXME: This may fail if is_reduced = True because
        #        and input_names includes state since these states are dropped
        is_state_input = len(set(self.input_names).intersection(
              self.state_names)) > 0
        if self.is_reduced and is_state_input:
            raise ValueError("Cannot handle input states for is_reduced=True")
        # Initializations 
        state_names = self.state_names
        # Delete states that are inputs
        state_names = list(set(self.state_names).difference(self.input_names))
        state_names = sorted(state_names,
              key=lambda n: self.state_names.index(n))
        num_state = len(state_names)
        num_output = len(self.output_names)
        if len(set(self.reaction_names).intersection(self.output_names)) > 0:
            flux_jacobian_df = self.makeFluxJacobian()
        else:
            flux_jacobian_df = None
        L0 = self.roadrunner.getL0Matrix()
        if len(L0) > 0:
            L0_df = pd.DataFrame(L0, columns=L0.colnames, index=L0.rownames)
        else:
            L0_df = None
        # Iterate across each output to construct the transpose of the C matrix
        C_T_dct = {}
        for name in self.output_names:
            if name in state_names:
                values = np.repeat(0, num_state)
                idx = state_names.index(name)
                values[idx] = 1
                C_T_dct[name] = values
            elif name in self.species_names:
                if L0_df is None:
                    raise RuntimeError("Species missing from L0: %s" % name)
                if not name in L0.rownames:
                    raise RuntimeError("Species missing from L0: %s" % name)
                C_T_dct[name] = L0_df.loc[name, :]
            elif name in self.reaction_names:
                C_T_dct[name] = flux_jacobian_df.loc[name, :]
        C_df_T = pd.DataFrame(C_T_dct, index=state_names)
        C_df = C_df_T.transpose()
        values = C_df.values.flatten()
        if any([np.isnan(v) for v in values]):
            raise RuntimeError("Nan value encountered.")
        return C_df

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
        Contructs vector of current state values.

        Returns
        -------
        np.ndarray: N X 1
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
            jacobian_mat = self.roadrunner.getFullJacobian()
        if len(jacobian_mat.rownames) != len(jacobian_mat.colnames):
            raise RuntimeError("Jacobian is not square!")
        names = list(jacobian_mat.colnames)
        self._jacobian_df = pd.DataFrame(jacobian_mat, columns=names, index=names)
        self._jacobian_time = self.getTime()
        return self._jacobian_df

    @property
    def state_names(self):
        state_names = list(self.jacobian_df.columns)
        return self._sortList(self.species_names, state_names)

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
        if time is not None:
            self.setTime(current_time)
        return jacobian_df

    def setTime(self, time):
        self.roadrunner.reset()
        self._jacobian_time = TIMEPOINT_NULL
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
        if names is None:
            names = self.roadrunner.keys()
        return util.getRoadrunnerValue(self.roadrunner, names)

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

    def _makeBDF(self, time=None):
        """
        Constructs a dataframe for the B matrix.
        The columns must be in the same order as the input_names.

        Parameters
        ---------
        time: float

        Returns
        -------
        np.ndarray (n X p), where p = len(input_names)
        """
        if len(self.input_names) > 0:
            # Determine which inputs are reactions and which inputs are species
            species_inputs, reaction_inputs = self.separateSpeciesReactionInputs()
            # Construct the matrix for species inputs
            if len(species_inputs) > 0:
                jacobian_df = self.getJacobian(time=time)
                # Don't include states that are input species
                df = jacobian_df.drop(species_inputs, axis=0)
                B_species_df = df[species_inputs]
                B_df = B_species_df
            if len(reaction_inputs) > 0:
                B_reaction_df = self.full_stoichiometry_df[reaction_inputs]
                B_df = B_reaction_df
            # Select the columns needed from the stoichiometry matrix
            # Merge the two
            if (len(species_inputs) > 0) and (len(reaction_inputs) > 0):
                B_df = pd.concat([B_species_df, B_reaction_df], axis=1)
            B_df = B_df[self.input_names]
        else:
            ncol = 1
            B_mat = np.repeat(0, self.num_state)
            B_mat = np.reshape(B_mat, (self.num_state, ncol))
            B_df = pd.DataFrame(B_mat, index=self.state_names)
        #
        return B_df

    def makeStateSpace(self, time=None, A_mat=None, B_mat=None,
          C_mat=None, D_mat=None):
        """
        Creates a state space control object for
        the n X n jacobian. By default, the D matrix is always 0.

        Parameters
        ----------
        The default values of the matrices are calculated in the constructor.
        These can be overridden.
        time: float (time at which Jacobian is obtained)
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
        if time is None:
            time = self.getTime()
        # The following linearization works if the following conditions hold:
        #  1. All floating species are in the state
        #  2. None of the matrices are specified
        #  3. States are used as input and output
        is_default_matrices = (A_mat is None) and (B_mat is None)  \
              and (C_mat is None)  and (D_mat is None)
        inputs_and_outputs = set(self.input_names).union(self.output_names)
        is_subset_state = inputs_and_outputs <= set(self.state_names)
        # FIXME: Delete?
        if False and (not self.is_reduced) and is_default_matrices and is_subset_state:
            cur_time = self.getTime()
            if time != self.getTime():
                self.setTime(time)
            nonlinear_sys = self.makeNonlinearIOSystem("system")
            state_values = self.state_ser.values
            ueq = np.repeat(1, nonlinear_sys.ninputs)
            state_space = control.linearize(nonlinear_sys, state_values,
                  ueq=ueq)
            if self.getTime() != cur_time:
                self.setTime(cur_time)
            return state_space
        # Do custom linearization
        # At least one matrix is specified
        # Construct the matrices
        A_mat = df2Mat(A_mat)
        B_mat = df2Mat(B_mat)
        C_mat = df2Mat(C_mat)
        D_mat = df2Mat(D_mat)
        # Calculate the Jacobian
        #
        if A_mat is None:
            A_df = self.getJacobian(time)
            columns = A_df.columns
            # Remove any state that's an input
            # Allow state to be generated internally
            for name in self.input_names:
                if name in columns:
                    A_df = A_df.drop(name, axis=0)
                    A_df = A_df.drop(name, axis=1)
            A_mat = A_df.values
        #
        if B_mat is None:
            B_df = self._makeBDF(time=time)
            if B_df is None:
                B_mat = None
            else:
                B_mat = B_df.values
        if C_mat is None:
            # Construct the output matrix
            C_mat = self.C_df.values
        if D_mat is None:
            nrow = len(self.output_names)
            if B_mat is None:
                ncol = 1
            else:
                ncol = np.shape(B_mat)[1]
            D_mat = np.repeat(0, nrow*ncol)
            D_mat = np.reshape(D_mat, (nrow, ncol))
        ss = control.StateSpace(A_mat, B_mat, C_mat, D_mat)
        return ss

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

    @staticmethod
    def reduceTransferFunction(tf, atol=ATOL):
        """
        Reduces the order of a transfer function if trailing zeroes.

        Parameters
        ----------
        tf: control.TransferFunction
        
        Returns
        -------
        tf: control.TransferFunction
        """
        def findOrderOfFirstNonzeroDigit(polynomial):
            for idx in range(len(polynomial)):
                pos = len(polynomial) - idx - 1
                if not np.isclose(polynomial[pos], 0, atol=atol):
                    return idx
            return len(polynomial) - 1
        def reduceOrder(polynomial, new_order):
            pos = len(polynomial) - new_order
            return polynomial[0:pos]
        #
        numerator = tf.num[0][0]
        denominator = tf.den[0][0]
        lowest_order = min(findOrderOfFirstNonzeroDigit(numerator),
              findOrderOfFirstNonzeroDigit(denominator))
        #
        new_numerator = reduceOrder(numerator, lowest_order)
        new_denominator = reduceOrder(denominator, lowest_order)
        new_tf = control.TransferFunction(new_numerator, new_denominator)
        return new_tf

    def makeTransferFunction(self, time=None, atol=ATOL):
        """
        Creates a transfer function for the system. Verifies that there
        is a single input and a single output. Reduces the order of the
        transfer function as needed.
        
        Parameters
        ----------
        time: float (time at which Jacobian is obtained)
        atol: absolute tolerance for comparison
        
        Returns
        -------
        control.TransferFunction
        """
        # Validity checks
        if len(self.input_names) != 1:
            raise ValueError("Must have exactly one input.")
        if len(self.output_names) != 1:
            raise ValueError("Must have exactly one output.")
        # Get initial transfer function
        state_space = self.makeStateSpace(time=time)
        tf = control.ss2tf(state_space)
        #
        return self.reduceTransferFunction(tf, atol=atol)

    def makeFluxJacobian(self, time=None):
        """
        Constructs the Jacobian of the reaction flux vector.

        Parameters
        ----------
        time: float
            Time at which this is calculated
        
        Returns
        -------
        pd.DataFrame
            index: reaction id
            column: species
        """
        # FIXME : Handle case where state_names < species_names
        # This calculation only works if state_names == species_names
        diff = set(self.state_names).symmetric_difference(self.species_names)
        if len(diff) > 0:
            raise RuntimeError(
                  "Code doesn't work if species_names != state_names")
        # Adjust time if necessary
        cur_time = self.getTime()
        if time is None:
            time = cur_time
        if not np.isclose(time, cur_time):
            self.setTime(time)
        #
        dct = {}
        for state_name in self.state_names:
            dct[state_name] = []
            for reaction_name in self.reaction_names:
                dct[state_name].append(self.roadrunner.getEE(
                      reaction_name, state_name))
        #
        df = pd.DataFrame(dct, index=self.reaction_names)
        if time != cur_time:
            self.setTime(cur_time)
        return df

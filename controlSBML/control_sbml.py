"""LTI Control for SBML models"""

"""
TO DO:

1. Plot difference between time jacoabian at reference vs. Current.
2. Plot TOTAL residual SSQ vs. jacobian difference
"""

import control
import numpy as np
import tellurium as te


ANT = "ant"
XML = "xml"


TYPE_MODEL = "type_model"  # libsbml model
TYPE_XML = "type_xml"  # XML string
TYPE_ANTIMONY = "type_xml"  # Antimony string
TYPE_FILE = "type_file" # File reference


class ControlSBML(object):

    def __init__(self, model_reference):
        """
        Initializes instance variables
        :param str model_reference: string or SBML file or Roadrunner object
        """
        ##### PUBLIC #####
        self.model_reference = model_reference
        self.roadrunner = None  # Roadrunner object
        # Read the model
        if "RoadRunner" in str(type(model_reference)):
            self.roadrunner = model_reference
            model_reference = model_reference.getSBML()
        elif isinstance(model_reference, str):
            parts = model_reference.split(".")
            if len(parts) == 2:
                if parts[1] == XML:
                    self.roadrunner = te.loadSBMLModel(model_reference)
                elif parts[1] == ANT:
                    self.roadrunner = te.loadAntimonyModel(model_reference)
                elif XML in model_reference.count:
                    self.roadrunner = te.loadSBMLModel(model_reference)
                else:
                    # Assume string for antimony model
                    self.roadrunner = te.loada(model_reference)
            else:
                self.roadrunner = te.loada(model_reference)
        else:
            raise ValueError("Invalid model reference")
        # Do the initializations
        self.antimony = self.roadrunner.getAntimony()
        self.roadrunner.reset()
        self.state_names = list(self.jacobian.colnames)

    @property
    def jacobian(self):
        """
        Returns
        -------
        NamedArray with names for rows (rownames) and columns (colnames)
        """
        mat = self.roadrunner.getFullJacobian()
        return mat

    @staticmethod
    def isRoadrunnerKey(key):
        return not ((key[0] == "_") or ("(" in key) or (key[-1] == "'"))

    def copy(self):
        """
        Creates a copy of the object.

        Returns
        -------
        controlSBML
        """
        ctlsb = ControlSBML(self.model_reference)
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
              in zip(self.state_names, self.jacobian.colnames)])
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

    def mkStateSpace(self, B=None, C=None, D=None):
        """
        Creates a control system object for the n X n jacobian.

        Parameters
        ----------
        B: np.array(n X p)
        C: np.array(q X n)
        D: np.array(q X p)

        Returns
        -------
        control.StateSpace
        """
        # Construct the matrices
        A = self.jacobian
        if B is None:
            B = np.repeat(0, A.shape[0])
            B = np.reshape(B, (A.shape[0], 1))
        if C is None:
            C = np.identity(A.shape[0])
        if D is None:
            D = B
        return control.StateSpace(A, B, C, D)

    def mkInitialState(self):
        """
        Contructs the initial state vector for StateSpace model.

        Returns
        -------
        np.array
        """
        values = list(self.get(self.state_names).values())
        return np.array(values)

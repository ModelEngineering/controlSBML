"""Constructs a sequential pathway in Antimony."""

import numpy as np


KINETICS_VALUE = 1
SPECIES_VALUE = 1
SPECIES_PFX = "S"
KINETICS_PFX = "k"

class SequentialModel(object):

    def __init__(self, num_reaction, kinetics_values=None,
          species_values=None, has_boundaries=True):
        """
        Parameters
        ----------
        num_reaction: int
            number of reactions in the sequential model
        kinetics_values: list-float
            len == num_reactions
        species_values: list-float
            len == num_reactions + 1
        has_boundaries: bool
            first and last species are fixed
        """
        self.num_reaction = num_reaction
        self.has_boundaries = has_boundaries
        self.num_species = self.num_reaction + 1
        if kinetics_values is None:
            kinetics_values = list(np.repeat(KINETICS_VALUE, self.num_reaction))
        self.kinetics_values = kinetics_values
        if species_values is None:
            species_values = list(np.repeat(SPECIES_VALUE, self.num_species))
        self.species_values = species_values

    @staticmethod
    def _makeKineticsName(idx):
        return "%s%d" % (KINETICS_PFX, idx)

    def _makeSpeciesName(self, idx):
        prefix = ""
        if self.has_boundaries:
            if (idx == 0) or (idx == self.num_species - 1):
                prefix = "$"
            else:
                prefix = ""
        return "%s%s%d" % (prefix, SPECIES_PFX, idx)

    def _makeKineticsAssignment(self, idx):
        return "%s = %f" % (self._makeKineticsName(idx),
              self.kinetics_values[idx])

    @staticmethod
    def _makeReactionLabel(idx):
        return "J%d: " % idx

    def _makeSpeciesInitialization(self, idx):
        return "%s = %f" % (self._makeSpeciesName(idx),
              self.species_values[idx])

    def _makeKineticsFormula(self, idx):
        """
        Creates a kinetics law.

        Parameters
        ----------
        idx: int
            index of the reaction
        
        Returns
        -------
        str
        """
        formula_str = "%s * %s" % (self._makeKineticsName(idx),
              self._makeSpeciesName(idx))
        return formula_str

    def _makeMassTransfer(self, idx1, idx2):
        return "%s -> %s" % (self._makeSpeciesName(idx1),
              self._makeSpeciesName(idx2))
        
    def generate(self):
        """
        Creates an antimony model.

        Returns
        -------
        str
        """
        stmts = []
        # Construct the reactions
        for idx in range(0, self.num_reaction):
            stmt = self._makeReactionLabel(idx)
            stmt = stmt + self._makeMassTransfer(idx, idx+1)
            stmt = stmt + "; " + self._makeKineticsFormula(idx)
            stmts.append(stmt)
        stmts.append("")
        # Construct the assignments
        for idx in range(0, self.num_reaction):
            stmts.append(self._makeKineticsAssignment(idx))
        for idx in range(self.num_species):
            stmts.append(self._makeSpeciesInitialization(idx))
        return "\n".join(stmts)
        

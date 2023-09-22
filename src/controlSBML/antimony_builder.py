"""Constructors Antimony to support control analysis and design."""

import controlSBML.constants as cn

import numpy as np

START_STR = "// ControlSBML: Start modifications VVVVVVVVVVVVVVV"
END_STR = "// ControlSBML: End modifications ^^^^^^^^^^^^^^^^"
COMMENT_STR = "//"
DEFAULT_NUM_STEP = 5
DEFAULT_INITIAL_VALUE = 0
DEFAULT_FINAL_VALUE = 10
DEFAULT_POINT_PER_STEP = 10

class AntimonyBuilder(object):

    def __init__(self, antimony, species_names):
        """
        Args:
            antimony: str (Antimony)
        """
        self.antimony = antimony
        self.species_names = species_names
        self.antimony_strs = antimony.split("\n")
        self.boundary_species = []
        for idx, stg in enumerate(self.antimony_strs):
            clean_stg = stg.strip()
            if len(clean_stg) == 0:
                continue
            if not clean_stg.startswith(COMMENT_STR):
                self.insert_pos = idx
                break

    def __repr__(self):
        return "\n".join(self.antimony_strs)

    def _insert(self, stg):
        """
        Args:
            stg: str
        """
        self.antimony_strs.insert(self.insert_pos, stg)
        self.insert_pos += 1

    def startModification(self):
        """
        Record the beginning of modifications to the antimony file.
        """
        self._insert(START_STR)

    def endModification(self):
        """
        Record the beginning of modifications to the antimony file.
        """
        self._insert(END_STR)

    def makeBoundarySpecies(self, species_name):
        """
        Args:
            species_name: str
        """
        self.boundary_species.append(species_name)
        self._insert("const %s" % species_name)

    def makeParameterNameForBoundaryReaction(self, species_name):
        """
        Args:
            species_name: str
        """
        return "_ControlSBML_k_%s" % species_name

    def makeBoundaryReaction(self, species_name):
        """
        Args:
            species_name: str
        """
        parameter_name = self.makeParameterNameForBoundaryReaction(species_name)
        reaction_str = " -> %s; %s" % (species_name, parameter_name)
        self._insert(reaction_str)
        initialization_str = "%s = 0" % parameter_name
        self._insert(initialization_str)

    def makeSISOClosedLoop(self, input_name, output_name, kp, ki, kd, kf):
        """
        Args:
            input_name: str
            output_name: str
            kp: float
            ki: float
            kd: float
            kf: float
        """
        raise NotImplementedError("makeSISOClosedLoop")
    
    def makeStaircase(self, input_name, times=cn.TIMES, initial_value=DEFAULT_INITIAL_VALUE,
                 num_step=DEFAULT_NUM_STEP, final_value=DEFAULT_FINAL_VALUE):
        """
        Adds events for the staircase.
        Args:
            species_name: str
            initial_value: float (value for first step)
            final_value: float (value for final step)
            num_step: int (number of steps in staircase)
            num_point_in_step: int (number of points in each step)
        Returns:
            array-float: values
        """
        # Find the name to use to effect the staircase
        if input_name in self.species_names:
            if not input_name in self.boundary_species:
                name = self.makeParameterNameForBoundaryReaction(input_name)
            else:
                name = input_name
        else:
            name = input_name
        #
        point_per_step = int(len(times)/num_step)
        step_size = (final_value - initial_value)/num_step
        values = []
        for nstep in range(num_step):
            values.extend([initial_value + nstep*step_size]*point_per_step)
            break_time = times[nstep*point_per_step]
            break_value = initial_value + nstep*step_size
            statement = "at (time>= %s): %s = %s" % (break_time, name, break_value)
            self._insert(statement)
        num_point_remaining = len(times) - len(values)
        values.extend([final_value for _ in range(num_point_remaining)])
        value_arr = np.array(values)
        return value_arr
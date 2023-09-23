"""Constructors Antimony to support control analysis and design."""

import controlSBML.constants as cn

import numpy as np


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
        self._initialized_output = False
        for idx, stg in enumerate(self.antimony_strs):
            clean_stg = stg.strip()
            if len(clean_stg) == 0:
                continue
            if not clean_stg.startswith(cn.COMMENT_STR):
                self.insert_pos = idx
                break

    def __repr__(self):
        return "\n".join(self.antimony_strs)

    def _insert(self, stg):
        """
        Args:
            stg: str
        """
        if not self._initialized_output:
            self._initialized_output = True
            self.antimony_strs.append("\n//--------------Aded by ControlSBML-----------------")
        self.antimony_strs.append(stg)
        self.insert_pos += 1

    def makeBoundarySpecies(self, species_name):
        """
        Args:
            species_name: str
        """
        self.boundary_species.append(species_name)
        self._insert("const %s" % species_name)

    def makeComment(self, comment):
        """
        Args:
            comment: str
        """
        self._insert("%s %s" % (cn.COMMENT_STR, comment))

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

    def makeClosedLoopSuffix(self, input_name, output_name):
        return "_%s_%s" % (input_name, output_name)

    def makeClosedLoopName(self, generic_name, input_name, output_name):
        """
        Creates the name of the symbol used in closed loop for the specified input_name and output_name

        Args:
            generic_name (_type_): _description_
            input_name (_type_): _description_
            output_name (_type_): _description_
        """
        suffix = self.makeClosedLoopSuffix(input_name, output_name)
        return "%s%s" % (generic_name, suffix)

    def makeSISOClosedLoop(self, input_name, output_name, kp=None, ki=None, kd=None, kf=None):
        """
        Args:
            input_name: str
            output_name: str
            kp: float
            ki: float
            kd: float
            kf: float
        """
        suffix = self.makeClosedLoopSuffix(input_name, output_name)
        # 
        if kf is not None:
            statement = "filter%s' = -%f*filter%s + %f*%s" % (
                suffix, kf, suffix, kf, output_name
            )
            self._insert(statement)
            statement = "filter%s = %s" % (suffix, output_name)
            self._insert(statement) 
        else:
            statement = "filter%s := %s" % (suffix, output_name)
            self._insert(statement) 
        statement = "control_error%s := reference%s - filter%s" % (
            suffix, suffix, suffix
        )
        self._insert(statement)
        statement = "integral_control_error%s' = control_error%s" % (
            suffix, suffix
        )
        self._insert(statement)
        if kp is not None:
            statement = "%s := %f*control_error%s" % (
                input_name, kp, suffix)
        else:
            statement = "%s = 0" % input_name
        if ki is not None:
            statement = statement + "+ %f*integral_control_error%s" % (ki, suffix)
        self._insert(statement)
        # Initialization statements
        statement = "integral_control_error%s = 0" % suffix
        self._insert(statement)
        statement = "reference%s = 0" % suffix
        self._insert(statement)
    
    def makeStaircase(self, input_name, times=cn.TIMES, initial_value=cn.DEFAULT_INITIAL_VALUE,
                 num_step=cn.DEFAULT_NUM_STEP, final_value=cn.DEFAULT_FINAL_VALUE):
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
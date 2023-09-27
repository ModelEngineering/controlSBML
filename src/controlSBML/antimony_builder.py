"""Constructors Antimony to support control analysis and design."""

import controlSBML.constants as cn

import numpy as np
import tellurium as te


class AntimonyBuilder(object):

    def __init__(self, antimony, symbol_dct, control_module_name=cn.DEFAULT_MODULE_NAME):
        """
        Args:
            antimony: str (Antimony)
            symbol_dct: dict
                key: str (symbol name)
                value: str (symbol type)
        """
        self.antimony = antimony
        self.symbol_dct = symbol_dct
        self.control_module_name = control_module_name
        self.antimony_strs = antimony.split("\n")
        self._initialized_output = False
        # Find the main module
        rr = te.loada(antimony)
        self.main_module_name = rr.model.getModelName()

    @property
    def floating_species_names(self):
        return [k for k, v in self.symbol_dct.items() if v == cn.TYPE_FLOATING_SPECIES]

    @property
    def boundary_species_names(self):
        return [k for k, v in self.symbol_dct.items() if v == cn.TYPE_BOUNDARY_SPECIES]

    def __repr__(self):
        outputs = list(self.antimony_strs)
        outputs.append("end")
        return "\n".join(outputs)

    def addStatement(self, stg):
        """
        Args:
            stg: str
        """
        if not self._initialized_output:
            self._initialized_output = True
            self.antimony_strs.append("\nmodule *%s()" % self.control_module_name)
            statement = "M: %s()" % self.main_module_name
            self.antimony_strs.append(statement)
            for name in self.symbol_dct.keys():
                statement = "M.%s is %s" % (name, name)
                self.antimony_strs.append(statement)
            self.antimony_strs.append("")
        self.antimony_strs.append(stg)

    def makeBoundarySpecies(self, species_name):
        """
        Args:
            species_name: str
        """
        self.symbol_dct[species_name] = cn.TYPE_BOUNDARY_SPECIES
        self.addStatement("const %s" % species_name)

    def makeComment(self, comment):
        """
        Args:
            comment: str
        """
        self.addStatement("%s %s" % (cn.COMMENT_STR, comment))

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
        self.addStatement(reaction_str)
        initialization_str = "%s = 0" % parameter_name
        self.addStatement(initialization_str)

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

    def makeSISOClosedLoop(self, input_name, output_name, kp=None, ki=None, kd=None, kf=None, reference=0):
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
            self.addStatement(statement)
            statement = "filter%s = %s" % (suffix, output_name)
            self.addStatement(statement) 
        else:
            statement = "filter%s := %s" % (suffix, output_name)
            self.addStatement(statement) 
        statement = "control_error%s := reference%s - filter%s" % (
            suffix, suffix, suffix
        )
        self.addStatement(statement)
        statement = "integral_control_error%s' = control_error%s" % (
            suffix, suffix
        )
        self.addStatement(statement)
        if kp is not None:
            statement = "%s := %f*control_error%s" % (
                input_name, kp, suffix)
        else:
            statement = "%s = 0" % input_name
        if ki is not None:
            statement = statement + "+ %f*integral_control_error%s" % (ki, suffix)
        self.addStatement(statement)
        # Initialization statements
        statement = "integral_control_error%s = 0" % suffix
        self.addStatement(statement)
        statement = "reference%s = %f" % (suffix, reference)
        self.addStatement(statement)

    def _makeInputManipulationName(self, input_name):
        """
        Constructs the name of the input that is being manipulated since floating species can be manipulated
        by a boundary reaction.

        Args:
            input_name: str
        """
        if input_name not in self.boundary_species_names:
            name = self.makeParameterNameForBoundaryReaction(input_name)
        else:
            name = input_name
        return name

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
        # Initialize the input
        name = self._makeInputManipulationName(input_name)
        statement = "%s = %f" % (name, initial_value)
        self.addStatement(statement)
        # Add events
        point_per_step = int(len(times)/num_step)
        step_size = (final_value - initial_value)/num_step
        values = []
        for nstep in range(num_step):
            values.extend([initial_value + nstep*step_size]*point_per_step)
            break_time = times[nstep*point_per_step]
            break_value = initial_value + nstep*step_size
            statement = "at (time>= %s): %s = %s" % (break_time, name, break_value)
            self.addStatement(statement)
        num_point_remaining = len(times) - len(values)
        values.extend([final_value for _ in range(num_point_remaining)])
        value_arr = np.array(values)
        return value_arr
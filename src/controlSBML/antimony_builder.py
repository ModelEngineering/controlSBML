"""
Constructors Antimony to support control analysis and design.

:Author: Joseph L. Hellerstein
:Date: 2023-09-01
:Email: joseph.hellerstein@gmail.com
:License: MIT

AntimonyBuilder creates Antimony code to construct control structures for design and analysis. Among
these are:
  PI controllers
  filters
  noise generators
Also provided is a mechanism for system identification by creating a staircase function.

Code is generated and placed in a new model module.

We use the following terms:
    signal: Something that generates data
    element: A component of a system with an input and an output
    system: A collection of elements

The following names are used for elements of a closed loop system:
  control_error: the difference between the setpoint and the output
  controller_in: the input to the controller
  controller_ot: the output of the controller
  filter_in: the input to the filter
  filter_ot: the output of the filter

Elements are connected by addition/subtraction. For example, we connect the filter output to the
control_error by the statement:
  control_error := filter_ot
This is an Antimony "assignment statement". It means that the values are equal for all time.
makeAdditionStatement(name1, name2, name3, ...) creates the statement:
  name1 := name2 + name3 + ...  
If a name is preceded by a "-" then it is subtracted.

The accommodate the presence of multiple controllers, the names of the elements are suffixed. A common
approach is the use the names of the input and output. For example, the control_error for a controller
that has an input named "S1" and an output named "S3" would be named "control_error_S1_S3".

"""

import controlSBML.constants as cn
from controlSBML import util

import numpy as np
import re
import tellurium as te

IN = "_in"
OT = "_ot"


class AntimonyBuilder(object):

    def __init__(self, antimony, symbol_dct=None):
        """
        Args:
            antimony: str (Antimony)
            symbol_dct: dict
                key: str (symbol name)
                value: str (symbol type)
        """
        self.antimony = antimony
        self.antimony_strs = antimony.split("\n")
        self.insert_pos = self._findMainModelEndPosition()  ## Add
        ##self.control_module_name = self._calculateControlModuleName(control_module_name)
        self._initialized_output = False
        # Find the main module
        rr = te.loada(antimony)
        if symbol_dct is None:
            symbol_dct = util.makeRoadrunnerSymbolDct(rr)
        self.symbol_dct = dict(symbol_dct)

    def copy(self):
        """
        Returns:
            AntimonyBuilder
        """
        builder = AntimonyBuilder(self.antimony, symbol_dct=self.symbol_dct.copy())
        builder.antimony_strs = list(self.antimony_strs)
        builder.insert_pos = self.insert_pos
        builder._initialized_output = self._initialized_output
        builder.symbol_dct = dict(self.symbol_dct)
        return builder
    
    def __eq__(self, other):
        is_equal = True
        is_debug = False
        is_equal &= self.antimony == other.antimony
        if is_debug and (not is_equal):
            print("Failed 1")
        is_equal &= util.allEqual(self.antimony_strs, other.antimony_strs)
        if is_debug and (not is_equal):
            print("Failed 2")
        is_equal &= self._initialized_output == other._initialized_output
        if is_debug and (not is_equal):
            print("Failed 3")
        is_equal = is_equal and (all([self.symbol_dct[k] == other.symbol_dct[k] for k in self.symbol_dct.keys()]))
        if is_debug and (not is_equal):
            print("Failed 4")
        is_equal = is_equal and (all([self.symbol_dct[k] == other.symbol_dct[k] for k in other.symbol_dct.keys()]))
        if is_debug and (not is_equal):
            print("Failed 5")
        is_equal &= self.insert_pos == other.insert_pos
        if is_debug and (not is_equal):
            print("Failed 6")
        return is_equal
    
    def _extractModelName(self, line):
        # Extracts the name of the model from the line
        start_pos = line.find("*") + 1
        end_pos = line.find("(")
        if (start_pos < 0) or (end_pos < 0) or (end_pos < start_pos):
            raise RuntimeError("Unable to extract model name from line: %s" % line)
        return line[start_pos:end_pos]

    def _findMainModelEndPosition(self):
        """
        Finds the position of the end statement for the main model.

        Returns:
            int (position in list of strings)
        """
        # Finds the name of the top level model
        model_names = []
        last_model_start = -1
        for pos, line in enumerate(self.antimony_strs):
            result = re.search("model .*[*].*()", line)
            if result:
                last_model_start = pos
        if last_model_start < 0:
            raise ValueError("Could not find a main model!")
        # Finds the end of the top level model
        for pos, line in enumerate(self.antimony_strs[last_model_start:]):
            new_line = line.strip()
            if new_line == "end":
                return pos + last_model_start
        raise ValueError("Could not find end of main model!")


    def _getSetOfType(self, antimony_type):
        return [k for k, v in self.symbol_dct.items() if v == antimony_type]
    
    @ property
    def species_names(self):
        return list(set(self.floating_species_names) + set(self.boundary_species_names))

    @property
    def floating_species_names(self):
        return self._getSetOfType(cn.TYPE_FLOATING_SPECIES)

    @property
    def boundary_species_names(self):
        names = self._getSetOfType(cn.TYPE_BOUNDARY_SPECIES)
        names = list(set(names))
        return names

    @property
    def reaction_names(self):
        return self._getSetOfType(cn.TYPE_REACTION)
    
    @property
    def parameter_names(self):
        return self._getSetOfType(cn.TYPE_PARAMETER)
    
    @property
    def assignment_names(self):
        return self._getSetOfType(cn.TYPE_ASSIGNMENT)

    def __repr__(self):
        return "\n".join([str(o) for o in self.antimony_strs])

    def addStatement(self, statement):
        """
        Args:
            statement: str
        """
        def insert(stg, increment=1):
            self.antimony_strs.insert(self.insert_pos, stg)
            self.insert_pos += increment
        #
        if not self._initialized_output:
            self._initialized_output = True
            insert("")
            insert("//vvvvvvvvvAdded by ControlSBMLvvvvvvvvvv")
            insert("//^^^^^^^^^Added by ControlSBML^^^^^^^^^^", increment=0)
        insert(statement)

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
        Makes a boundary reaction to regulate a species concentration.
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
    
    def makeAdditionStatement(self, *pargs, is_assignment=True):
        """
        Creates an addition statement that looks for a leanding "-".

        Args:
            *pargs: str
            is_assignment: bool (True means that the statement is an assignment statement)
        """
        statement = pargs[0]
        if is_assignment:
            statement += " := "
        else:
            statement += " = "
        for idx, argument in enumerate(pargs[1:]):
            argument_str = str(argument)
            if argument_str[0] == "-":
                if idx == 0:
                    statement += " -" + argument_str[1:]
                else:
                    statement += " - " + argument_str[1:]
            else:
                if idx > 0:
                    statement += " + " + argument_str
                else:
                    statement += argument_str
        self.addStatement(statement) 

    def makeSinusoidSignal(self, amplitude, frequency, prefix="sinusoid", suffix=""):
        """
        Makes a sinusoid. The created variable is prefix + suffix + "_ot"
        Prefix is used to scope within a control loop. Suffix is used to scope between control loops.

        Args:
            amplitude: float
            frequency: float
            prefix: str (beginning of the name)
            suffix: str (ending of the name)
        Returns:
            str: name of the sinusoid
        """
        self.addStatement("")
        self.makeComment("Make sinusoid: amplitude=%s, frequency=%s" % (str(amplitude), str(frequency)))
        name = prefix +  suffix +  OT
        statement = "%s := %f*sin(2*pi*%f*time)" % (name, amplitude, frequency)
        self.addStatement(statement) 
        return name
    
    def _makeInputOutputNames(self, prefix, suffix):
        base_name = prefix + suffix
        name_in = base_name + IN
        name_ot = base_name + OT
        return name_in, name_ot
    
    def makeFilterElement(self, kf, prefix="filter", suffix=""):
        """
        Makes a filter. prefix + suffix + IN is the input and prefix + suffix + OT is the output.
        Prefix is used to scope within a control loop. Suffix is used to scope between control loops.

        Args:
            kf: float
            prefix: str (beginning of the name)
            suffix: str (ending of the name)
        Returns:
            str: name of the filter input
            str: name of the filter output
        Usage:
            suffix = "_S1_S3"
            name_in, name_ot = self.makeFilter(0.5, suffix=suffix)
            self.makeAddition(name_in, "S3")   # S3 is the output of the system
            self.makeAddition("control_error", setpoint, "-"+name_ot)
        """
        self.addStatement("")
        self.makeComment("Filter: kf=%s" % (str(kf)))
        name_in, name_ot = self._makeInputOutputNames(prefix, suffix)
        if (kf is not None) and (kf > 0):
            statement = "%s' =  -%f*%s + %f*%s" % (name_ot, kf, name_ot, kf, name_in)
            self.addStatement(statement)
            self.makeAdditionStatement(name_ot, 0, is_assignment=False)   # Initialize the filter output
        else:
            self.makeAdditionStatement(name_ot, name_in)
        return name_in, name_ot
    
    def makeControlErrorSignal(self, setpoint, forward_output_name, sign, prefix="control_error", suffix=""):
        """
        Constructs the control error variable.

        Args:
            setpoint: float/str
            forward_output_name: str (Output for which the setpoint is compared)
            sign: float
            prefix: str (beginning of the name)
            suffix: str (ending of the name)
        Returns:
            str (name of the control error)
        Usage:
            suffix = "_S1_S3"
            control_error_name = self.makeControlError(setpoint, "S3", suffix=suffix)   # S3 is the output of the system
        """
        name_ot = prefix + suffix + OT
        if sign == 1:
            operator = "+"
        else:
            operator = "-"
        statement = "%s := %s %s %s" % (name_ot, str(setpoint), operator, forward_output_name)
        self.addStatement(statement)
        return name_ot
    
    def makePIControllerElement(self, kp=None, ki=None, prefix="controller", suffix=""):
        """
        Makes a PI controller. prefix + suffix + IN is the input and prefix + suffix + OT is the output.
        Prefix is used to scope within a control loop. Suffix is used to scope between control loops.

        Args:
            kp: float
            ki: float
            prefix: str (beginning of the name)
            suffix: str (ending of the name)
        Returns:
            str: name of the controller input
            str: name of the controller output
        Usage:
            suffix = "_S1_S3"
            name_in, name_ot = self.makeController(kp=1, suffix=suffix)
            control_error_name = self.makeControlError(setpoint, "S3", suffix=suffix)   # S3 is the output of the system
            self.makeAddition(name_in, control_error_name)  # control_error is input to the controller
            self.makeAddition("S1", name_ot)   # S1 is the input to the system and is set by the controller
        """
        self.addStatement("")
        self.makeComment("PI Controller: kp=%s, ki=%s" % (str(kp), str(ki)))
        name_in, name_ot = self._makeInputOutputNames(prefix, suffix)
        # Make the integral of the control error
        integral_error_name = prefix + "_integral_error" + suffix
        statement = "%s' = %s" % (integral_error_name, name_in)
        self.addStatement(statement)
        self.makeAdditionStatement(integral_error_name, 0, is_assignment=False)   #  integral_error_name = 0
        # Construct the control law
        if kp is not None:
            statement = "%s := %f*%s" % (name_ot, kp, name_in)
        else:
            statement = "%s = 0" % name_ot
        if ki is not None:
            statement = statement + "+ %f*%s" % (ki, integral_error_name)
        self.addStatement(statement)
        return name_in, name_ot

    def makeSISOClosedLoopSystem(self, input_name, output_name, kp=None, ki=None, kf=None, setpoint=0,
                           noise_amplitude=0, noise_frequency=20, disturbance_ampliude=0, disturbance_frequency=20,
                           initial_output_value=None, sign=-1):
        """
        Args:
            input_name: str (input to system)
            output_name: str (output from system)
            kp: float
            ki: float
            kd: float
            kf: float
            setpoint: float (setpoint)
            noise_amplitude: float (Amplitude of the additions to the output)
            noise_frequency: float (Frequency of the additions to the output)
            disturbance_amplitude: float (Amplitude of the disturbance)
            disturbance_frequency: float (Frequency of the disturbance)
            initial_input_value: float (initial value of the input)
        """
        suffix = self.makeClosedLoopSuffix(input_name, output_name)
        # Handle the initial value of the input
        if initial_output_value is not None:
            self.makeAdditionStatement(output_name, initial_output_value, is_assignment=False)
        # Make the elements of the closed loop
        noise_ot = self.makeSinusoidSignal(noise_amplitude, noise_frequency, prefix="noise", suffix=suffix)
        disturbance_ot = self.makeSinusoidSignal(disturbance_ampliude, disturbance_frequency, prefix="disturbance", suffix=suffix)
        filter_in, filter_ot = self.makeFilterElement(kf, prefix="filter", suffix=suffix)
        controller_in, controller_ot = self.makePIControllerElement(kp, ki, prefix="controller", suffix=suffix)
        control_error_name = self.makeControlErrorSignal(setpoint, filter_ot, sign, prefix="control_error", 
                                                         suffix=suffix)
        # Connect the pieces by specifying assignment statements
        self.addStatement("")
        self.makeComment("Connect the elements of the closed loop")
        self.makeAdditionStatement(controller_in, control_error_name)
        self.makeAdditionStatement(input_name, controller_ot, disturbance_ot)
        self.makeAdditionStatement(filter_in, output_name, noise_ot)
    
    def getInputManipulationName(self, input_name):
        """
        Constructs the name of the input that is being manipulated since floating species can be manipulated
        by a boundary reaction rate.

        Args:
            input_name: str
        """
        if (input_name in self.floating_species_names):
            # The input is a floating species that is being controlled by a boundary reaction
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
        name = self.getInputManipulationName(input_name)
        statement = "%s = %f" % (name, initial_value)
        self.addStatement(statement)
        # Add events
        point_per_step = int(len(times)/(num_step+1))
        step_size = (final_value - initial_value)/(num_step)
        values = []
        for nstep in range(num_step + 1):
            values.extend([initial_value + nstep*step_size]*point_per_step)
            break_time = times[nstep*point_per_step]
            break_value = initial_value + nstep*step_size
            statement = "at (time>= %s): %s = %s" % (break_time, name, break_value)
            self.addStatement(statement)
        num_point_remaining = len(times) - len(values)
        values.extend([final_value for _ in range(num_point_remaining)])
        value_arr = np.array(values)
        return value_arr
"""Creates a system from an SBML model: inputs, input type, outputs."""

import controlSBML.constants as cn
from controlSBML.antimony_builder import AntimonyBuilder
from controlSBML import msgs
from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML.timeseries import Timeseries
from controlSBML import util

import re
import tellurium as te


REFERENCE = "reference"


class SBMLSystem(object):

    def __init__(self, model_reference, input_names=None, output_names=None, is_fixed_input_species=False,
                 control_module_name=cn.DEFAULT_MODULE_NAME, model_id="model"):
        """
        model_reference: str
            string, SBML file or Roadrunner object
        input_names: name (id) of reaction whose flux is being controlled
                     or the name of a chemical species
        output_names: list-str
            output species
        is_fixed_input_species: bool (input species are fixed)
        control_module_name: str
            name of the control module
        model_id: str (identifier of the model)
        """
        if input_names is None:
            input_names = []
        if output_names is None:
            output_names = []
        # First initializations
        self.model_reference = model_reference
        self.model_id = model_id
        self.is_fixed_input_species = is_fixed_input_species
        self.roadrunner = makeRoadrunner(self.model_reference)
        # Validate the input and output names
        self.input_names = [self.makeInputName(n, self.roadrunner) for n in input_names]
        self.output_names = [self.makeOutputName(n, self.roadrunner) for n in output_names]
        # Verify that the main model is a module
        self.antimony = self._getAntimony()
        self.symbol_dct = self._makeSymbolDct()
        try:
            self.antimony_builder = AntimonyBuilder(self.antimony, self.symbol_dct)
        except Exception as exp:
            msgs.error("Cannot create AntimonyBuilder: %s" % exp)
        if self.antimony_builder.parent_model_name in [None, cn.DEFAULT_ROADRUNNER_MODULE]:
            msgs.error("Cannot find the name of the parent model. Is it a modular model?")
        # Add boundary information depending on the type of input
        for name in self.input_names:
            if name in self.antimony_builder.floating_species_names:
                if is_fixed_input_species:
                    self.antimony_builder.makeBoundarySpecies(name)
                else:
                    self.antimony_builder.makeBoundaryReaction(name)
        
    def _getAntimony(self):
        """
        Extracts antimony from the roadrunner object. Cleans it as necessary.

        Returns:
            str
        """
        antimony = self.roadrunner.getAntimony()
        # FIXME: Handle this better?
        if "unknown_model_qual" in antimony:
            msgs.warn("Antimony contains 'unknown_model_qual' which is not supported by tellurium. Replaced with 'description'.")
            antimony = antimony.replace(" unknown_model_qual ", " description ")
        return antimony
        
    def model_name(self):
        model_name = self.roadrunner.model.getModelName()
        if not " " in model_name:
            return model_name
        else:
            return self.antimony_builder.model_name
        
    def _makeSymbolDct(self):
        """
        Returns
        -------
        dict
            key: str (symbol name)
            value: str (symbol type)
        """
        symbol_dct = util.makeRoadrunnerSymbolDct(self.roadrunner)
        names = list(self.input_names)
        names.extend(self.output_names)
        final_symbol_dct = {k: v for k, v in symbol_dct.items() if k in names}
        return final_symbol_dct
        
    def _isExistingName(self, name, names):
        if name in self.roadrunner.keys():
            return True
        names.append(name)
        return False

    @staticmethod
    def _getValidNames(names):
        return [n for n in names if not n in ["at", "in"]]

    @staticmethod
    def _getName(names, name, name_type):
        """"
        Finds a valid name for input or output.

        names: list-str (permitted names)
        """
        valid_names = SBMLSystem._getValidNames(names)
        if len(valid_names) == 0:
            msgs.error("No %s name." % (name_type))
        if name is None:
            name = valid_names[0]
        if not name in valid_names:
            msgs.error("name %s not in %s" % (name, valid_names))
        return name

    @classmethod
    def makeInputName(cls, name, roadrunner):
        """
        Finds valid name to use for model input. Checks validity of input names.

        Args:
            name: str
            roadrunner: ExtendedRoadrunner 
        Returns:
            str (valid input name)
        """
        # Find the valid input and output names. Input and output should be different.
        input_names = roadrunner.getFloatingSpeciesIds()
        input_names.extend(roadrunner.getBoundarySpeciesIds())
        input_names.extend(roadrunner.getAssignmentRuleIds())
        input_names.extend(roadrunner.getGlobalParameterIds())
        new_input_name = cls._getName(input_names, name, name_type="input")
        return new_input_name
    
    @classmethod
    def makeOutputName(cls, name, roadrunner, input_names=None):
        """
        Finds valid name to use for model output. Checks validity of output names.

        Args:
            name: str
            roadrunner: ExtendedRoadrunner
            input_names: list-str (name of input. If None, then not checked to see if output is differs from input)
        Returns:
            str (valid input name)
        """
        if input_names is None:
            input_names = []
        output_names = roadrunner.getFloatingSpeciesIds()
        output_names.extend(roadrunner.getAssignmentRuleIds())
        output_names.extend(roadrunner.getReactionIds())
        new_name = cls._getName(output_names, name, name_type="output")
        for _ in input_names:
            if new_name in input_names:
                output_names.remove(new_name)
                new_name = cls._getName(output_names, name, name_type="output")
            else:
                break
        return new_name
    
    def get(self, name):
        """
        Provides the roadrunner values for a name.

        Parameters
        ----------
        name: str

        Returns
        -------
        object
        """
        return self.roadrunner[name]

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
    
    def simulate(self, start_time=cn.START_TIME, end_time=cn.END_TIME, num_point=None, is_steady_state=False):
        """
        Simulates the system.

        Parameters
        ----------
        start_time: float
        end_time: float
        num_point: int
        is_steady_state: bool (start the simulation at steady state)

        Returns
        -------
        DataFrame
        """
        if num_point is None:
            num_point = int(cn.POINTS_PER_TIME*(end_time - start_time))
        self.roadrunner.reset()
        if is_steady_state:
            self.setSteadyState()
        return self._simulate(start_time, end_time, num_point, is_steady_state, is_reload=False)
    
    def _simulate(self, start_time, end_time, num_point, is_steady_state=False, is_reload=True):
        """
        Simulates the system the roadrunner object.

        Parameters
        ----------
        start_time: float
        end_time: float
        num_point: int
        is_steady_state: bool (start the simulation at steady state)

        Returns
        -------
        DataFrame
        """
        if is_reload:
            self.roadrunner = te.loada(str(self.antimony_builder))
        if is_steady_state:
            self.setSteadyState()
        selections = list(self.input_names)
        selections.extend(self.output_names)
        selections.insert(0, cn.TIME)
        data = self.roadrunner.simulate(start_time, end_time, num_point, selections=selections)
        ts = Timeseries(data)
        return ts
    
    def simulateSISOClosedLoop(self, input_name=None, output_name=None, kp=None, ki=None, kf=None, reference=1,
                               start_time=cn.START_TIME, end_time=cn.END_TIME, num_point=None, is_steady_state=False):
        """
        Simulates a closed loop system.

        Args:
            input_name: str
            output_name: str
            kp: float
            ki float
            kf: float
            reference: float (setpoint)
        Returns:
            Timeseries
        """
        if input_name is None:
            input_name = self.input_names[0]
        if output_name is None:
            output_name = self.output_names[0]
        self.antimony_builder.addStatement("")
        comment = "Closed loop: %s -> %s" % (input_name, output_name)
        self.antimony_builder.makeComment(comment)
        #
        if self.is_fixed_input_species and (input_name in self.antimony_builder.floating_species_names):
            new_input_name = input_name
        else:
            new_input_name = self.antimony_builder.makeParameterNameForBoundaryReaction(input_name)
        self.antimony_builder.makeSISOClosedLoopSystem(new_input_name, output_name, kp=kp, ki=ki, kf=kf, reference=reference)
        # Run the simulation
        return self._simulate(start_time, end_time, num_point, is_steady_state, is_reload=True)
    
    def simulateStaircase(self, input_name, output_name, times=cn.TIMES, initial_value=cn.DEFAULT_INITIAL_VALUE,
                 num_step=cn.DEFAULT_NUM_STEP, final_value=cn.DEFAULT_FINAL_VALUE, is_steady_state=True):
        """
        Adds events for the staircase.
        Args:
            input_name: str
            output_name: str
            initial_value: float (value for first step)
            final_value: float (value for final step)
            num_step: int (number of steps in staircase)
            num_point_in_step: int (number of points in each step)
        Returns:
            Timeseries
        """
        self.antimony_builder.addStatement("")
        self.antimony_builder.makeComment("Staircase: %s->%s" % (input_name, output_name))
        self.antimony_builder.makeStaircase(input_name, times=times, initial_value=initial_value,
                                            num_step=num_step, final_value=final_value)
        ts = self._simulate(start_time=times[0], end_time=times[-1], num_point=len(times),
                            is_steady_state=is_steady_state)
        return ts
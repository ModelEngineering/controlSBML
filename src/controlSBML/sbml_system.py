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

    def __init__(self, model_reference, input_names, output_names, is_fixed_input_species=False,
                 control_module_name=cn.DEFAULT_MODULE_NAME):
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
        """
        # First initializations
        self.model_reference = model_reference
        self.input_names = input_names
        self.output_names = output_names
        self.is_fixed_input_species = is_fixed_input_species
        self.roadrunner = makeRoadrunner(self.model_reference)
        # Verify that the main model is a module
        self.antimony = self._getAntimony()
        self.symbol_dct = self._makeSymbolDct()
        try:
            self.antimony_builder = AntimonyBuilder(self.antimony, self.symbol_dct)
        except Exception as exp:
            import pdb; pdb.set_trace()
            pass
        if self.antimony_builder.parent_model_name in [None, cn.DEFAULT_ROADRUNNER_MODULE]:
            msgs.error("Cannot find the name of the parent model. Is it a modular model?")
        # Add boundary information depending on the type of input
        for name in self.input_names:
            if name in self.antimony_builder.floating_species_names:
                if is_fixed_input_species:
                    self.antimony_builder.makeBoundarySpecies(name)
                else:
                    self.antimony_builder.makeBoundaryReaction(name)
        # Validation checks
        invalid_names = self._verifyInputNames()
        if len(invalid_names) > 0:
            text = "Inputs must be a species, parameter, or compartment."
            text += "   Invalid names are: %s" % str(invalid_names)
            raise ValueError(text)
        invalid_names = self._verifyOutputNames()
        if len(invalid_names) > 0:
            text = "Outputs must be species or fluxes."
            text += "The following outputs are invalid: %s" % str(invalid_names)
            raise ValueError(text)
        
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
        
    def _verifyInputNames(self):
        """
        Verifies that the input names are valid.

        Returns
        -------
        bool
        """
        invalid_names = []
        for name in self.input_names:
            if not self._isExistingName(name, invalid_names):
                continue
            if name in self.antimony_builder.reaction_names:
                invalid_names.append(name)
                continue
            try:
                value = self.get(name)
                self.set({name: value})
            except Exception:
                invalid_names.append(name)
        return invalid_names
    
    def _verifyOutputNames(self):
        """
        Verifies that the output names are valid.

        Returns
        -------
        bool
        """
        invalid_names = []
        for name in self.output_names:
            if not self._isExistingName(name, invalid_names):
                continue
            if name in self.antimony_builder.reaction_names:
                invalid_names.append(name)
                continue
            if name in self.antimony_builder.boundary_species_names:
                invalid_names.append(name)
                continue
            try:
                # Must be able to retrieve a value
                _ = self.get(name)
            except RuntimeError:
                invalid_names.append(name)
        return invalid_names

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
        data = self.roadrunner.simulate(start_time, end_time, num_point)
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
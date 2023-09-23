"""Creates a system from an SBML model: inputs, input type, outputs."""

import controlSBML.constants as cn
from controlSBML.antimony_builder import AntimonyBuilder
from controlSBML import msgs
from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML.timeseries import Timeseries
from controlSBML import util

import tellurium as te


REFERENCE = "reference"


class SBMLSystem(object):

    def __init__(self, model_reference, input_names, output_names, is_fixed_input_species=False):
        """
        model_reference: str
            string, SBML file or Roadrunner object
        input_names: name (id) of reaction whose flux is being controlled
                     or the name of a chemical species
        output_names: list-str
            output species
        is_fixed_input_species: bool (input species are fixed)
        """
        # First initializations
        self.model_reference = model_reference
        self.input_names = input_names
        self.output_names = output_names
        self.is_fixed_input_species = is_fixed_input_species
        self.roadrunner = makeRoadrunner(self.model_reference)
        self.floating_species_names = list(set(self.roadrunner.getFloatingSpeciesIds()))
        self.antimony = self.roadrunner.getAntimony()
        self.antimony_builder = AntimonyBuilder(self.antimony, self.floating_species_names)
        for name in self.input_names:
            if name in self.floating_species_names:
                if is_fixed_input_species:
                    self.antimony_builder.makeBoundarySpecies(name)
                else:
                    self.antimony_builder.makeBoundaryReaction(name)
        # Handle inputs/outputs
        stoichiometry_mat = self.roadrunner.getFullStoichiometryMatrix()
        self.floating_species_names = list(self.roadrunner.getFloatingSpeciesIds())
        self.reaction_names = list(stoichiometry_mat.colnames)
        self.parameter_names = list(self.roadrunner.getGlobalParameterIds())
        self.num_input = len(self.input_names)
        self.num_output = len(self.output_names)
        # Validation checks
        invalid_names = self._invalidInputNames()
        if len(invalid_names) > 0:
            text = "Inputs must be a species, parameter, or compartment."
            text += "   Invalid names are: %s" % str(invalid_names)
            raise ValueError(text)
        invalid_names = self._invalidOutputNames()
        if len(invalid_names) > 0:
            text = "Outputs must be species or fluxes."
            text += "The following outputs are invalid: %s" % str(invalid_names)
            raise ValueError(text)
        
    def _isExistingName(self, name, names):
        if name in self.roadrunner.keys():
            return True
        names.append(name)
        return False
        
    def _invalidInputNames(self):
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
            if name in self.reaction_names:
                invalid_names.append(name)
                continue
            try:
                value = self.get(name)
                self.set({name: value})
            except Exception:
                invalid_names.append(name)
        return invalid_names
    
    def _invalidOutputNames(self):
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
        comment = "Closed loop: %s -> %s" % (input_name, output_name)
        self.antimony_builder.makeComment(comment)
        #
        if self.is_fixed_input_species and (input_name in self.floating_species_names):
            self.antimony_builder.makeBoundarySpecies(input_name)
            new_input_name = input_name
        else:
            self.antimony_builder.makeBoundaryReaction(input_name)
            new_input_name = self.antimony_builder.makeParameterNameForBoundaryReaction(input_name)
        self.antimony_builder.makeSISOClosedLoop(new_input_name, output_name, kp=kp, ki=ki, kf=kf)
        reference_name = self.antimony_builder.makeClosedLoopName(REFERENCE, input_name, output_name)
        self.roadrunner = te.loada(str(self.antimony_builder))
        self.set({reference_name: reference})
        return self._simulate(start_time, end_time, num_point, is_steady_state)
    
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
        self.antimony_builder.makeComment("Staircase: %s->%s" % (input_name, output_name))
        self.antimony_builder.makeStaircase(input_name, times=times, initial_value=initial_value,
                                            num_step=num_step, final_value=final_value)
        ts = self._simulate(start_time=times[0], end_time=times[-1], num_point=len(times),
                            is_steady_state=is_steady_state)
        return ts
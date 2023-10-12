"""Creates a system from an SBML model: inputs, input type, outputs."""

import controlSBML.constants as cn
from controlSBML.antimony_builder import AntimonyBuilder
from controlSBML import msgs
from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML.timeseries import Timeseries
from controlSBML import util

import matplotlib.pyplot as plt
import tellurium as te


REFERENCE = "reference"


class SBMLSystem(object):

    def __init__(self, model_reference, input_names=None, output_names=None, is_fixed_input_species=False,
                 model_id="model"):
        """
        model_reference: str
            string, SBML file or Roadrunner object
        input_names: name (id) of reaction whose flux is being controlled
                     or the name of a chemical species
        output_names: list-str
            output species
        is_fixed_input_species: bool (input species are fixed)
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
            raise ValueError("Cannot create AntimonyBuilder: %s" % exp)
        # Add boundary information depending on the type of input
        for name in self.input_names:
            if name in self.antimony_builder.floating_species_names:
                if is_fixed_input_species:
                    self.antimony_builder.makeBoundarySpecies(name)
                else:
                    self.antimony_builder.makeBoundaryReaction(name)

    def copy(self):
        model_reference = str(self)
        system = SBMLSystem(self.model_reference, input_names=self.input_names, output_names=self.output_names,
                            is_fixed_input_species=self.is_fixed_input_species, model_id=self.model_id)
        system.antimony_builder = self.antimony_builder.copy()
        system.symbol_dct = dict(self.symbol_dct)
        return system
    
    def __eq__(self, other):
        is_debug = False
        is_equal = True
        is_equal = is_equal and (self.model_reference == other.model_reference)
        if is_debug:
            print("Failed 1")
        is_equal = is_equal and (self.model_id == other.model_id)
        if is_debug:
            print("Failed 2")
        is_equal = is_equal and (self.input_names == other.input_names)
        if is_debug:
            print("Failed 3")
        is_equal = is_equal and (self.is_fixed_input_species == other.is_fixed_input_species)
        if is_debug:
            print("Failed 4")
        is_equal = is_equal and (util.allEqual(self.input_names, other.input_names))
        if is_debug:
            print("Failed 3")
        is_equal = is_equal and (util.allEqual(self.output_names, other.output_names))
        if is_debug:
            print("Failed 4")
        is_equal = is_equal and (self.antimony == other.antimony)
        if is_debug:
            print("Failed 5")
        is_equal = is_equal and (self.antimony_builder == other.antimony_builder)
        if is_debug:
            print("Failed 6")
        is_equal = is_equal and (all([self.symbol_dct[k] == other.symbol_dct[k] for k in self.symbol_dct.keys()]))
        if is_debug:
            print("Failed 7")
        is_equal = is_equal and (all([self.symbol_dct[k] == other.symbol_dct[k] for k in other.symbol_dct.keys()]))
        if is_debug:
            print("Failed 8")
        return is_equal
        
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
            raise ValueError("No %s name." % (name_type))
        if name is None:
            name = valid_names[0]
        if not name in valid_names:
            raise ValueError("name %s not in %s" % (name, valid_names))
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
        output_names = []
        output_names.extend(roadrunner.getFloatingSpeciesIds())
        output_names.extend(roadrunner.getReactionIds())
        output_names.extend(roadrunner.getAssignmentRuleIds())
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
    
    def _simulate(self, start_time, end_time, num_point, antimony_builder=None, is_steady_state=False, is_reload=True):
        """
        Simulates the system the roadrunner object.

        Parameters
        ----------
        start_time: float
        end_time: float
        num_point: int
        antimoney_builder: AntimonyBuilder
        is_steady_state: bool (start the simulation at steady state)

        Returns
        -------
        DataFrame
        """
        if antimony_builder is None:
            antimony_builder = self.antimony_builder
        if is_reload:
            self.roadrunner = te.loada(str(antimony_builder))
        if is_steady_state:
            self.setSteadyState()
        selections = list(self.input_names)
        selections.extend(self.output_names)
        selections.insert(0, cn.TIME)
        data = self.roadrunner.simulate(start_time, end_time, num_point, selections=selections)
        ts = Timeseries(data)
        return ts
    
    def simulateSISOClosedLoop(self, input_name=None, output_name=None, kp=None, ki=None, kf=None, setpoint=1,
                               start_time=cn.START_TIME, end_time=cn.END_TIME, num_point=None,
                               is_steady_state=False, inplace=False):
        """
        Simulates a closed loop system.

        Args:
            input_name: str
            output_name: str
            kp: float
            ki float
            kf: float
            setpoint: float (setpoint)
            inplace: bool (update the existing model with the closed loop statements)
        Returns:
            Timeseries
            AntimonyBuilder
        """
        if input_name is None:
            input_name = self.input_names[0]
        if output_name is None:
            output_name = self.output_names[0]
        if inplace:
            builder = self.antimony_builder
        else:
            builder = self.antimony_builder.copy()
        builder.addStatement("")
        comment = "Closed loop: %s -> %s" % (input_name, output_name)
        builder.makeComment(comment)
        #
        if input_name in builder.boundary_species_names:
            new_input_name = input_name
        else:
            new_input_name = builder.makeParameterNameForBoundaryReaction(input_name)
        builder.makeSISOClosedLoopSystem(new_input_name, output_name, kp=kp, ki=ki, kf=kf, setpoint=setpoint)
        # Run the simulation
        return self._simulate(start_time, end_time, num_point, is_steady_state=is_steady_state,
                              antimony_builder=builder, is_reload=True), builder
    
    def simulateStaircase(self, input_name, output_name, times=cn.TIMES, initial_value=cn.DEFAULT_INITIAL_VALUE,
                 num_step=cn.DEFAULT_NUM_STEP, final_value=cn.DEFAULT_FINAL_VALUE, is_steady_state=True, inplace=True):
        """
        Adds events for the staircase.
        Args:
            input_name: str
            output_name: str
            initial_value: float (value for first step)
            final_value: float (value for final step)
            num_step: int (number of steps in staircase)
            num_point_in_step: int (number of points in each step)
            inplace: bool (update the existing model with the Staircase statements)
        Returns:
            Timeseries
            AntimonyBuilder
        """
        if inplace:
            builder = self.antimony_builder
        else:
            builder = self.antimony_builder.copy()
        builder.addStatement("")
        builder.makeComment("Staircase: %s->%s" % (input_name, output_name))
        builder.makeStaircase(input_name, times=times, initial_value=initial_value,
                                            num_step=num_step, final_value=final_value)
        ts = self._simulate(start_time=times[0], antimony_builder=builder, end_time=times[-1], num_point=len(times),
                            is_steady_state=is_steady_state)
        return ts, builder
    
    def plotSISOClosedLoop(self, timeseries, setpoint, **kwargs):
        """
        Plots the results of a closed lop simulation. Input and output are defined in the SBMLSystem constructor.

        Args:
            timeseries: Timeseries
            setpoint: float
            kwargs: dict (kwargs for plotOneTS)
        """
        is_plot, new_kwargs = util.setNoPlot(kwargs)
        new_kwargs["is_plot"] = False
        plot_result = util.plotOneTS(timeseries, ax2=0, **new_kwargs)
        ax = plot_result.ax
        times = timeseries.index/cn.MS_IN_SEC
        ax.plot(times, [setpoint]*len(times), color="red", linestyle="--")
        legends = list(timeseries.columns)
        legends.append("setpoint")
        ax.legend(legends)
        if is_plot:
            plt.show()
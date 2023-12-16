"""Creates a system from an SBML model: inputs, input type, outputs."""

import controlSBML.constants as cn
from controlSBML.antimony_builder import AntimonyBuilder
from controlSBML import msgs
from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML.timeseries import Timeseries
from controlSBML import util
from controlSBML.option_management.option_manager import OptionManager

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tellurium as te


class SBMLSystem(object):

    def __init__(self, model_reference, input_names=None, output_names=None, is_fixed_input_species=False,
                 model_id="model", is_steady_state=False, roadrunner=None):
        """
        model_reference: str
            string, SBML file or Roadrunner object
        input_names: name (id) of reaction whose flux is being controlled
                     or the name of a chemical species
        output_names: list-str
            output species
        is_fixed_input_species: bool (input species are fixed)
        model_id: str (identifier of the model)
        is_steady_state: bool (start the simulation at steady state)
        """
        if input_names is None:
            input_names = []
        if output_names is None:
            output_names = []
        self._specified_input_names = input_names
        self._specified_output_names = output_names
        # First initializations
        self.model_reference = model_reference
        self.model_id = model_id
        self.is_fixed_input_species = is_fixed_input_species
        self.is_steady_state = is_steady_state
        self.is_fixed_input_species = is_fixed_input_species
        # The following are calculated on first reference
        self._antimony_builder = None
        self._roadrunner = roadrunner
        self._input_names = None
        self._output_names = None
        self._original_antimony = None
        self._symbol_dct = None
        

    #################### PROPERTIES ####################
    @property
    def original_antimony(self):
        if self._original_antimony is None:
            self._original_antimony = self._getAntimony()
        return self._original_antimony
    
    @property
    def symbol_dct(self):
        if self._symbol_dct is None:
            self._symbol_dct = self._makeSymbolDct()
        return self._symbol_dct

    @property
    def roadrunner(self):
        if self._roadrunner is None:
            self._roadrunner = makeRoadrunner(self.model_reference)
        return self._roadrunner

    # FIXME: Does not copy updates to the antimony?
    @property
    def antimony_builder(self):
        if self._antimony_builder is None:
            try:
                self._antimony_builder = AntimonyBuilder(self.original_antimony, self.symbol_dct)
            except Exception as exp:
                raise ValueError("Cannot create AntimonyBuilder: %s" % exp)
            # Create boundary information depending on the type of input
            for name in self.input_names:
                if name in self._antimony_builder.floating_species_names:
                    if self.is_fixed_input_species:
                        self._antimony_builder.makeBoundarySpecies(name)
                    else:
                        self._antimony_builder.makeBoundaryReaction(name)
        return self._antimony_builder
    
    @property
    def input_names(self):
        if self._input_names is None:
            self._input_names = [self.makeInputName(n, self.roadrunner) for n in self._specified_input_names]
        return self._input_names
    
    @property
    def output_names(self):
        if self._output_names is None:
            self._output_names = [self.makeOutputName(n, self.roadrunner) for n in self._specified_output_names]
        return self._output_names
            
    #################### METHODS ####################
    def copy(self):
        """
        Copies the SBML object.

        Returns:
            _type_: _description_
        """
        system = SBMLSystem(self.model_reference, input_names=self._specified_input_names,
                            output_names=self._specified_output_names,
                            is_fixed_input_species=self.is_fixed_input_species, model_id=self.model_id,
                            is_steady_state=self.is_steady_state)
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
        is_equal = is_equal and (self.original_antimony == other.original_antimony)
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
    
    def simulate(self, start_time=cn.START_TIME, end_time=cn.END_TIME, num_point=None, is_steady_state=None):
        """
        Simulates the system.

        Parameters
        ----------
        start_time: float
        end_time: float
        num_point: int
        is_steady_state: bool (start the simulation at steady state, overrides constructor setting)

        Returns
        -------
        DataFrame
        """
        if is_steady_state is None:
            is_steady_state = self.is_steady_state
        if num_point is None:
            num_point = int(cn.POINTS_PER_TIME*(end_time - start_time))
        if is_steady_state:
            self.setSteadyState()
        return self._simulate(start_time, end_time, num_point, is_steady_state, is_reload=False)
    
    def _simulate(self, start_time, end_time, num_point, antimony_builder=None, is_steady_state=False, is_reload=False):
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
        antimony = str(antimony_builder)
        if antimony[-1] == "\n":
            antimony = antimony[:-1]
        if is_reload:
            try:
                self._roadrunner = te.loada(str(antimony_builder))
            except Exception as exp:
                raise ValueError("Cannot reload antimony: %s" % exp)
        else:
            self.roadrunner.reset()
        if is_steady_state:
            self.setSteadyState()
        selections = list(self.input_names)
        selections.extend(self.output_names)
        selections.insert(0, cn.TIME)
        data = self.roadrunner.simulate(start_time, end_time, num_point, selections=selections)
        ts = Timeseries(data)
        return ts
    
    def isInitialized(self)->bool:
        """
        Determines if the system has been initialized.

        Returns
        -------
        bool
        """
        return self._antimony_builder is not None

    # FIXME: Delete since duplicated in controlSBML 
    def plotModel(self, start_time=cn.START_TIME, end_time=cn.END_TIME, num_point=None, **kwargs):
        """
        Plots the original model.

        Args:
            start_time: float
            end_time: float
            num_point: int
            kwargs: dict (kwargs for plotOneTS)
        Returns:
            Timeseries
        """
        if num_point is None:
            num_point = int(cn.POINTS_PER_TIME*(end_time - start_time))
        roadrunner = makeRoadrunner(self.model_reference)
        data = roadrunner.simulate(start_time, end_time, num_point)
        ts = Timeseries(data)
        util.plotOneTS(ts, **kwargs)
        return ts

    def simulateSISOClosedLoop(self, input_name=None, output_name=None, kp=None, ki=None, kf=None, setpoint=1,
                               start_time=cn.START_TIME, end_time=cn.END_TIME, times=None, num_point=None,
                               is_steady_state=False, inplace=False, initial_input_value=None,
                               sign=-1):
        """
        Simulates a closed loop system.

        Args:
            input_name: str
            output_name: str
            kp: float
            ki float
            kf: float
            setpoint: float (setpoint)
            times: np.ndarray (times for the simulation)
            start_time: float (overridden by times)
            end_time: float (overridden by times)
            num_point: int (overridden by times)
            inplace: bool (update the existing model with the closed loop statements)
            initial_input_value: float (initial value of the input)
            sign: float (sign of the feedback)
        Returns:
            Timeseries
            AntimonyBuilder
        """
        if times is not None:
            start_time = times[0]
            end_time = times[-1]
            num_point = len(times)
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
        elif input_name in builder.floating_species_names:
            new_input_name = builder.makeParameterNameForBoundaryReaction(input_name)
        else:
            new_input_name = input_name
        builder.makeSISOClosedLoopSystem(new_input_name, output_name, kp=kp, ki=ki, kf=kf, setpoint=setpoint,
                                         initial_output_value=initial_input_value, sign=sign)
        # Run the simulation
        result = self._simulate(start_time, end_time, num_point, is_steady_state=is_steady_state,
                            antimony_builder=builder, is_reload=True), builder
        return result
    
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
                            is_steady_state=is_steady_state, is_reload=True)
        return ts, builder
    
    @staticmethod 
    def setYAxColor(ax, position, color):
        # Set the colors of the labels, axes, and spines
        ax.tick_params(axis='y', labelcolor=color)
        ax.spines[position].set_color(color)
        ax.yaxis.label.set_color(color)

    def getName(self, is_input=True):
        if is_input:
            names = self.input_names
            name_type ="input"
        else:
            names = self.output_names
            name_type ="output"
        if len(names) < 1:
            raise ValueError("No %s name is specified." % name_type)
        return names[0]
    
    def plotSISOClosedLoop(self, timeseries, setpoint, mgr=None, markers=None, **kwargs):
        """
        Plots the results of a closed lop simulation. Input and output are defined in the SBMLSystem constructor.

        Args:
            timeseries: Timeseries
            setpoint: float
            kwargs: dict (kwargs for plotOneTS)
        """
        input_name = self.getName(is_input=True)
        output_name = self.getName(is_input=False)
        is_plot, new_kwargs = util.setNoPlot(kwargs)
        if mgr is None:
            mgr = OptionManager(new_kwargs)
        new_kwargs["is_plot"] = False
        df = pd.DataFrame(timeseries[output_name], columns=[output_name])
        new_kwargs.setdefault(cn.O_AX2, 0)
        plot_result = util.plotOneTS(df, colors=[cn.SIMULATED_COLOR], markers=markers, **new_kwargs)
        ax = plot_result.ax
        ax.set_ylabel(output_name)
        # Plot the setpoint
        setpoint_arr = np.ones(len(timeseries))*setpoint
        times = np.array(timeseries.index)/cn.MS_IN_SEC
        ax.plot(times, setpoint_arr, color=cn.SIMULATED_COLOR,
            linestyle="--")
        # Do the plots
        mgr.plot_opts.set(cn.O_AX, ax)
        if mgr.plot_opts[cn.O_AX2] is None:
            ax2 = ax.twinx()
            mgr.plot_opts[cn.O_AX2] = ax2
        else:
            ax2 = mgr.plot_opts[cn.O_AX2]
        # Plot the staircase
        ax2.plot(times, timeseries[input_name], color=cn.INPUT_COLOR)
        self.setYAxColor(ax, "left", cn.SIMULATED_COLOR)
        self.setYAxColor(ax2, "right", cn.INPUT_COLOR)
        ax2.set_ylabel(input_name)
        #mgr.doPlotOpts()
        ax.legend([])
        #mgr.doFigOpts()
        if is_plot:
            plt.show()

    def printModel(self, is_updated_model=True):
        """
        Prints the antimony for the model.

        Args:
            is_updated_model: bool (include updates generated for system identification, closed loop design, ...)
        """
        if is_updated_model:
            print(self.antimony_builder)
        else:
            print(self.original_antimony)

    def getValidSymbols(self, is_input=True, is_str=True, is_updated_model=False):
        """
        Provides the names of valid input or ouput for the model.

        Args:
            is_input: bool (True: input, False: output)
            is_str: bool (True: return as a string, False: return as a Series)
            is_updated_model: bool (True: include updates generated for system identification, closed loop design, ...)

        Returns:
            pd.Series or str
                index: type
                value: list of symbol names
        """
        if is_updated_model:
            roadrunner = self.roadrunner
        else:
            roadrunner = makeRoadrunner(self.model_reference)
        #
        symbol_dct = util.makeRoadrunnerSymbolDct(roadrunner)
        if is_input:
            valid_types = cn.VALID_INPUT_TYPES
        else:
            valid_types = cn.VALID_OUTPUT_TYPES
        dct = {}
        for rr_type in valid_types:
            dct[rr_type] = ", ".join([k for k, v in symbol_dct.items() if v == rr_type])
        ser = pd.Series(dct)
        #
        if not is_str:
            return ser
        outs = []
        for rr_type in ser.index:
            if len(ser[rr_type]) > 0:
                out = rr_type + ":" + "\t" + ser[rr_type]
                outs.append(out)
        return "\n\n".join(outs)

    def getValidInputs(self):
        """
        Provides the names of valid model inputs.

        Returns:
            str
        """
        return self.getValidSymbols(is_input=True, is_str=True)
        

    def getValidOutputs(self):
        """
        Provides the names of valid model inputs.

        Returns:
            str
        """
        return self.getValidSymbols(is_input=False, is_str=True)
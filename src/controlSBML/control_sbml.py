"""Control engineering for SBML models"""

"""
The features of ControlSBML are organized as follows:
* control_base: intializations, get, set, copy, make
* control_analysis: creates analysis dataframes and series
* control_plot: plot data

This module exists to avoid circular imports so that objects that depend on ControlSBML can be constructed.
"""

from controlSBML import constants as cn
from controlSBML.control_plot import ControlPlot
import controlSBML.siso_transfer_function_builder as stb
from controlSBML.mimo_transfer_function_builder import MIMOTransferFunctionBuilder
from controlSBML.staircase import Staircase
from controlSBML.option_management.option_manager import OptionManager

from docstring_expander.expander import Expander


class ControlSBML(ControlPlot):

    def __init__(self, model_reference, input_names=None, output_names=None,
          is_reduced=False):
        """
        Initializes instance variables
        model_reference: str
            string, SBML file or Roadrunner object
        input_names: name (id) of reaction whose flux is being controlled
        output_names: list-str
            output species
        is_reduced: bool
            construct a reduced model so that the A matrix is nonsingular
        """
        super().__init__(model_reference,
              input_names=input_names, output_names=output_names,
              is_reduced=is_reduced)
        
    def makeSISOTransferFunctionBuilder(self, system_name="sys",
          input_name=None, output_name=None, **kwargs):
        """
        Creates a SISOTransferFunctionBuilder to construct a SISO transfer function.
        The default input and output names are the first input and output names.

        Parameters
        ----------
        system_name: str
        input_name: str
        output_name: str
        kwargs: dict (additional arguments for NonlinearIOSystem)

        Returns
        -------
        SISOTransferFunctionBuilder
        """
        def getName(specified_name, names):
            if specified_name is None:
                if len(names) < 1:
                    raise ValueError("Must specify at least one name")
                name = names[0]
            else:
                name = specified_name
            return name
        #
        input_name = getName(input_name, self.input_names)
        output_name = getName(output_name, self.output_names)
        sys = self.makeNonlinearIOSystem(system_name, **kwargs)
        return stb.SISOTransferFunctionBuilder(sys, input_name=input_name,
              output_name=output_name)
    
    def makeMIMOTransferFunctionBuilder(self, **kwargs):
        """
        Constructs transfer functions for SBML systems.

        Parameters
        ----------
        kwargs: dict (additional arguments for TransferFunctionBuilder)
        """
        return MIMOTransferFunctionBuilder(self, **kwargs)
    
    @Expander(cn.KWARGS, cn.SIM_KWARGS)
    def makeMIMOTransferFunctionDF(self, num_numerator=cn.DEFAULT_NUM_NUMERATOR, num_denominator=cn.DEFAULT_NUM_DENOMINATOR,
                 staircase=None, 
                 is_fixed_input_species=False, do_simulate_on_update=False,
                 is_steady_state=True, **kwargs):
        """
        Constructs transfer functions between all inputs and outputs.

        Parameters
        ----------
        num_numerator: int (number of numerator coefficients)
        num_denominator: int (number of denominator coefficients)
        staircase:
            None: use default relative staircase (steadystate +/- 0.5 fractional devition)
            float: use relative staircase (steadystate +/- staircase)
            Staircase: use staircase
            dct: dictionary of staircases
        is_steady_state: bool (operating point is steady state)
        is_fixed_input_species: bool (input species are fixed)
        do_simulate_on_update: bool (simulate on update to state to handle functions of time)
        #@expand

        Returns
        ------
        pd.DataFrame
            index: input_name
            columns: output_name
            values: control.TransferFunction
        """
        mgr = OptionManager(kwargs)
        if (staircase is None) or (isinstance(staircase, float)):
            # Find the centers for the staircase functions used in system identification
            self.setSteadyState()
            center_dct = {n: self.get(n) for n in self.state_names if n in self.input_names}
            if isinstance(staircase, float):
                deviation = staircase
            else:
                deviation = cn.DEFAULT_DEVIATION
            staircase_dct = {n: Staircase.makeRelativeStaircase(
                center=center_dct[n], fractional_deviation=deviation) for n in self.state_names if n in self.input_names}
        else:
            staircase_dct = staircase
        builder = MIMOTransferFunctionBuilder(self, is_fixed_input_species=is_fixed_input_species,
                                              do_simulate_on_update=do_simulate_on_update)
        fitter_result_df = builder.fitTransferFunction(num_numerator=num_numerator, num_denominator=num_denominator,
                                    staircase=staircase_dct, is_steady_state=is_steady_state, **mgr.sim_opts)
        transfer_function_df = fitter_result_df.copy()
        for input_name in transfer_function_df.index:
            for output_name in transfer_function_df.columns:
                transfer_function_df.loc[input_name, output_name] = fitter_result_df.loc[input_name, output_name].transfer_function
        return transfer_function_df
    
    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotBode(self, num_numerator=cn.DEFAULT_NUM_NUMERATOR, num_denominator=cn.DEFAULT_NUM_DENOMINATOR,
                 staircase=None, is_magnitude=True, is_phase=True, **kwargs):
        """
        Constructs bode plots for a MIMO system whose tansfer functions are obtained by system identification.

        Parameters
        ----------
        staircase:
            None: use default relative staircase (steadystate +/ 0.5)
            float: use relative staircase (steadystate +/- staircase)
            Staircase: use staircase
            dct: dictionary of staircases
        is_magnitude: bool
            Do magnitude plots
        is_phase: bool
            Do phase plots
        is_plot: bool
            Display plots
        #@expand

        Returns
        ------
        PlotResult
        """
        transfer_function_df = self.makeMIMOTransferFunctionDF(num_numerator=num_numerator,
                                                               num_denominator=num_denominator,
                                                               staircase=staircase, **kwargs)
        self.plotBodeTF(transfer_function_df, is_magnitude=is_magnitude, is_phase=is_phase, **kwargs)
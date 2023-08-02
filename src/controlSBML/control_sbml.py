"""Control engineering for SBML models"""

"""
The features of ControlSBML are organized as follows:
* control_base: intializations, get, set, copy, make
* control_analysis: creates analysis dataframes and series
* control_plot: plot data

This module exists to avoid circular imports so that objects that depend on ControlSBML can be constructed.
"""

from controlSBML.control_plot import ControlPlot
import controlSBML.siso_transfer_function_builder as stb
from controlSBML.nonlinear_io_system import NonlinearIOSystem


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
    
    def makeNonlinearIOSystem(self, name, **kwargs):
        """
        Creates an object that can be used in connections with the
        control package.

        Parameters
        ----------
        name: str (name of the system)
        kwargs: dict (additional arguments for NonlinearIOSystem)

        Returns
        -------
        controlSBML.NonelinearIOSystem
        """
        return NonlinearIOSystem(name, self, input_names=self.input_names, output_names=self.output_names,
                                     **kwargs)
        

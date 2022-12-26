"""LTI Control for SBML models"""

"""
The features of ControlSBML are organized as follows:
* control_base: intializations, get, set, copy, make
* control_analysis: creates analysis dataframes and series
* control_plot: plot data
"""

from controlSBML.control_plot import ControlPlot


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
        

"""LTI Control for SBML models"""

"""
The features of ControlSBML are organized as follows:
* control_base: intializations, get, set, copy, make
* control_analysis: creates analysis dataframes and series
* control_plot: plot data
"""

from controlSBML.control_plot import ControlPlot


class ControlSBML(ControlPlot):

    def __init__(self, model_reference, include_boundary_species=True):
        """
        Initializes instance variables
        :param str model_reference: string or SBML file or Roadrunner object
        """
        super().__init__(model_reference,
              include_boundary_species=include_boundary_species)

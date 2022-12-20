.. controlSBML documentation master file, created by
   sphinx-quickstart on Sun Dec 18 20:06:56 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

``controlSBML``: A package for control analysis and design for SBML models
==========================================================================

``controlSBML`` provides a bridge between the CalTech
control systems library
(https://sourceforge.net/projects/python-control/)
and
models represented in the System Biology Markup Language (SBML)
community standard.
The approach is:

1. Construct a ``ControlSBML`` object for the SBML model.

2. From this object, obtain objects used in the CalTech package:

   a. ``NonlinearIOSystem`` object that wraps the entire SBML model; and
   
   b. ``StateSpace`` object that is a linear approximation to the SBML model

3. Do control analysis and design using the CalTech package.

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

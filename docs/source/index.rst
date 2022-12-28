.. controlSBML documentation master file, created by
   sphinx-quickstart on Sun Dec 18 20:06:56 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

``controlSBML``: Control engineering with SBML models
======================================================

``controlSBML`` provides a bridge between SBML models and the
control systems library.
At the lowest level, this bridge wraps an SBML model as an object known
to the controls package so that SBML models can be analyzed.
The approach is:

1. Construct a ``ControlSBML`` object for the SBML model.

2. From this object, obtain objects used in the CalTech package:

   a. ``NonlinearIOSystem`` object that wraps the entire SBML model; and
   
   b. ``StateSpace`` object that is a linear approximation to the SBML model

3. Do control analysis and design using the CalTech package.

Below, we provide an overview of the features of ``controlSBML``.
Readers who are less familiar with control engineering
are encouraged to read the "Concepts" section.

Beyond this, ``controlSBML`` facilitates the construction of
models of controllers and filters to facilitate the evaluation of candidate
designs.
This is done through the use of *element factories* that create
``NonlinearIOSystem`` objects that can be used in the evaluation of
closed loop designs.

A third capability of ``controlSBML`` is to automate model construction
for an entire closed loop system.
At present, the package only supports models that have a single control input
and a single signal used in feedback (SISO systems).

There is a fourth category of features provided by ``controlSBML``.
These are convenience methods for control analysis and design.
Typically, these are a simplified way to access capabilities that
are provided by the ``control`` pacakge, possibly with some new features.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   concepts
   system_models
   closed_loop_models
   element_factories
   template_closed_loop_systems
   convenience_methods



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

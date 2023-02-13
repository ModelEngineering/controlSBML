.. controlSBML documentation master file, created by
   sphinx-quickstart on Sun Dec 18 20:06:56 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

``controlSBML``: Control engineering with SBML models
=============================================================

``controlSBML`` provides a bridge between SBML models and the
`CalTech control systems library <https://python-control.readthedocs.io/en/latest/intro.html>`_.
At the lowest level, this bridge wraps an SBML model as an object known
to the ``control`` package so that its tools can be
applied to SBML models.
The approach is:

1. Construct a ``controlSBML.ControlSBML`` object for the SBML model.

2. From this object, obtain objects used in the ``control`` package.

   a. ``control.NonlinearIOSystem`` objects wrap the entire SBML model to construct computational testbeds for closed loop systems.
   
   b. ``control.StateSpace`` objects are linear approximations to the SBML model for control analysis and design using transfer functions.

3. Do control analysis and design using the ``control`` package.

Below, we provide an overview of the features of ``controlSBML``.
Readers who are less familiar with control engineering
are encouraged to read the
:ref:`Concepts section <concepts>`:.

In addition to providing a bridge between SBML and the ``control`` package,
``controlSBML`` provides the following.

1. Automating the construction of models of elements in a closed loop system; in particular, controllers and filters.  This is done through the use of **element factories** that create ``NonlinearIOSystem`` objects that can be used in the evaluation of closed loop designs.

2. Methods to aid in system identification of SBML models as transfer functions.

3. Automating model construction for an entire closed loop system.  These are referred to as  **system factories**.  At present, the package only supports models that have a single control input and a single signal used in feedback (SISO systems).

4. Convenience methods for control analysis and design.  Typically, these are a simplified way to access capabilities that are provided by the ``control`` pacakge, possibly with some new features.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   concepts
   system_models
   siso_transfer_function_builder
   modeling_closed_loops
   element_factories
   system_factories
   convenience_methods
   detailed_example



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

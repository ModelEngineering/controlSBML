.. controlSBML documentation master file, created by
   sphinx-quickstart on Sun Dec 18 20:06:56 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

``controlSBML``: Control engineering with SBML models
=============================================================

``controlSBML`` is a python package for doing control analysis and design in
systems biology
where the open loop system is a model described by the
community standard Systems Biology Markup Language
(SBML).
SBML models may be specified as strings, files, or URLs. They may be in the
XML representation or in Antimony.

At present, ``controlSBML`` only supports the development of SISO (single input, single output)
control loops. The input can be a chemical species or an SMBL parameter (e.g.,
a kinetic constant). The output can be a chemical species or a reaction flux.
The package supports the following workflow:

* access the SBML model
* analyze controllability (the ability of an input to control the output)
* system identification (find a transfer function that fits the input-output response)
* control design for proportional-integral-differential (PID) control with a filter
* analyze the effect of perturbations (either disturbances to the control input to the open loop system or noise in the output measured from the open loop system)

``controlSBML`` leverages the
`CalTech control systems library <https://python-control.readthedocs.io/en/latest/intro.html>`_.
The control package is used
to do classical control analysis such as bode plots and root locus plots.


Below, we provide an overview of the features of ``controlSBML``.
Readers who are less familiar with control engineering
are encouraged to read the
:ref:`Concepts section <concepts>`:.

In addition to providing a bridge between SBML and the ``control`` package,
``controlSBML`` provides the following.

1. Automating the construction of models of elements in a closed loop system; in particular, controllers and filters.  This is done through the use of **element factories** that create ``NonlinearIOSystem`` objects that can be used in the evaluation of closed loop designs.

2. Methods to aid in system identification of SBML models as transfer functions.

3. **System factories** for constructing commonly used elements of closed loop control.

4. **Closed loop factories** for building an entire closed loop systems.

5. Convenience methods for control analysis and design.  Typically, these are a simplified way to access capabilities that are provided by the ``control`` pacakge, possibly with some new features.

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   installation
   concepts
   system_models
   siso_transfer_function_builder
   modeling_closed_loops
   element_factories
   closed_loop_factories
   convenience_methods
   detailed_example



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

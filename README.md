# controlSBML
``controlSBML`` provides: (1) a bridge between SBML models and the python ``control`` package and (2) capabilities for the analysis and design of control systems for SBML models.

## Installation
``pip install controlSBML``

To find the current version:
```
import controlSBML as ctl
ctl.__version__
```

## Constructing ControlSBML
``ctlsb = ControlSBML(model_reference)`` where
``model_reference``can be any of the following:
* path to a local file with the extension ``.ant`` or ``.xml``
* a python string for an SBML or Antimony model representation
* a URL to an XML file

## Key properties
* ``species_names`` is a list of species in the model
* ``roadrunner`` is the roadrunner object for the model
* ``jacobian`` is the Jacobian of the model at the current simulation time. Changes to ``roadrunner`` can change this property.
* ``antimony`` is the antimony representation of the model

## Key methods
* ``setTime(new_time)`` resets ``roadrunner`` and runs a simulation from time zero to the time specified.
* ``makeStateSpace()`` creates a ``control`` state space object using the Jacobian of the model at the current simulation time.
* ``plotTrueModel`` plots a roadrunner simulation.
* ``plotLinearApproximation`` plots the linear approximation provided by an ${\bf A}$ matrix. The default ${\bf A}$ matrix
is the Jacobian of the current ``RoadRunner`` instance. You can either specify a different ${\bf A}$ or set the ``RoadRunner`` instance to a different time.

## Example

[This notebook](https://github.com/ModelEngineering/controlSBML/blob/main/notebooks/UsingControlSBML.ipynb) provides a simple example of using ``controlSBML``.

## Version History
* 0.1.4,
  * More options for plotting and simulations
* 0.1.3, 2/13/2022. First release.

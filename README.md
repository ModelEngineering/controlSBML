# controlSBML
``controlSBML`` provides a bridge between SBML models and the python ``control`` library. Some analysis facilities are provided as well.

## Installation
``pip install controlSBML``

## Usage
``ctlsb = ControlSBML(model_reference)``
``model_reference``can be any of the following:
* path to a local file with the extension ``.ant`` or ``.xml``
* a python string for an SBML or Antimony model representation
* a URL to an XML file

## Key properties
* ``species_names`` is a list of species in the model
* ``roadrunner`` is the roadrunner object for the model
* ``jacobian`` is the Jacobian of the model at the current simulation time. Changes to ``roadrunner`` may change this property.
* ``antimony`` is the antimony representation of the model

## Key methods
* ``setTime(new_time)`` resets ``roadrunner`` and runs a simulation from time zero to the time specified.
* ``makeStateSpace()`` creates a ``control`` state space object using the Jacobian of the model at the current simulation time.

## Example



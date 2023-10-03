"""Analyzes BioModels for Verification."""

from controlSBML import constants as cn
from controlSBML.sbml_system import SBMLSystem
from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML.iterate_biomodels import iterateBiomodels

import matplotlib.pyplot as plt
import numpy as np
import os
import tellurium as te

INPUT_COLOR = "red"
OUTPUT_COLOR = "blue"

def checkFile(filename, contents):
    """Check the file for errors.
    Args:
        filename: str
        contents: str
    Returns:
        str
    """
    plot_name = filename + ".png"
    if os.path.isfile(plot_name):
        return ("Already processed %s" % filename)
    else:
        return ""

def getValidNames(names):
    return [n for n in names if not n in ["at", "in"]]

def doStaircase(filename, model, input_name=None, output_name=None, initial_value=0,
                end_time=20, final_value=10):
    # Initializations
    plot_name = filename.replace(".xml", ".png")
    try:
        roadrunner = te.loads(model)
    except Exception as error:
        print("Error for %s: %s" % (filename, error))
        return
    if (input_name is None) or (output_name is None):
        floating_species = roadrunner.getFloatingSpeciesIds()
        floating_species = getValidNames(floating_species)
    if input_name is None:
        if len(floating_species) < 1:
            print("No floating species for %s" % filename)
            return
        input_name = floating_species[0]
    if output_name is None:
        if len(floating_species) < 2:
            print("Only one floating species for %s" % filename)
            return
        output_name = floating_species[1]  
    # Construct the staircase
    try:
        system = SBMLSystem(roadrunner, [input_name], [output_name], is_fixed_input_species=False)
        times = np.linspace(0, end_time, 10*end_time)
        ts = system.simulateStaircase(input_name, output_name, times=times, final_value=final_value, num_step=5, is_steady_state=False)
    except Exception as error:
        print("Error for %s: %s" % (filename, error))
        return
    # Make the plot
    inputs = ts[input_name].values
    outputs = ts[output_name].values
    times = ts.index/1000
    _, ax = plt.subplots(1)
    ax.scatter(times, outputs, color=OUTPUT_COLOR)
    ax2 = ax.twinx()
    ax2.plot(times, inputs, color=INPUT_COLOR)
    ax2.set_ylabel(input_name, color=INPUT_COLOR)
    ax.set_ylabel(output_name, color=OUTPUT_COLOR)
    ax.set_xlabel("time")
    ax.set_title(filename)
    plt.savefig(plot_name)

def main(start=0, end=1e4):
    iterator = iterateBiomodels(start=start, end=end, is_report=True, checkerFunctions=[checkFile])
    for filename, model in iterator:
        doStaircase(filename, model)

if __name__ == '__main__':
    main(start=1)
"""Analyzes BioModels for Verification."""

from controlSBML import constants as cn
from controlSBML.sbml_system import SBMLSystem
from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML import util

import matplotlib.pyplot as plt
import numpy as np
import os

MODEL_DIR = os.path.abspath(__file__)
for _ in range(3):
    MODEL_DIR = os.path.dirname(MODEL_DIR)
MODEL_DIR = os.path.join(MODEL_DIR, "SBMLModel")
MODEL_DIR = os.path.join(MODEL_DIR, "data")
INPUT_COLOR = "red"
OUTPUT_COLOR = "blue"
IGNORE_FILES = ["BIOMD0000000075.xml", "BIOMD0000000081.xml", "BIOMD0000000353.xml",
                  "BIOMD0000000573.xml",
                  "BIOMD0000000627.xml"]


def getValidNames(names):
    return [n for n in names if not n in ["at", "in"]]

def doStaircase(filename, input_name=None, output_name=None, initial_value=0,
                end_time=20, final_value=10):
    # Initializations
    path = os.path.join(MODEL_DIR, filename)
    try:
        roadrunner = makeRoadrunner(path)
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
    plot_name = filename.replace(".xml", ".png")
    plt.savefig(plot_name)

def main(start=0, end=None):
    filenames = [f for f in os.listdir(MODEL_DIR) if f.endswith(".xml")]
    filenames.sort()
    for idx, filename in enumerate(filenames):
        if filename in IGNORE_FILES:
            print("%d: Ignoring %s" % (idx, filename))
            continue
        if idx < start:
            print("%d: Skipping %s" % (idx, filename))
            continue
        if (end is not None) and (idx >= end):
            break
        print("%d: Processing %s" % (idx, filename))
        doStaircase(filename)

if __name__ == '__main__':
    main(start=620, end=None)
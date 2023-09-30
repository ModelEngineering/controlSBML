"""Analyzes BioModels for Verification."""

from controlSBML import constants as cn
from controlSBML.sbml_system import SBMLSystem
from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML import util

import matplotlib.pyplot as plt
import numpy as np
import os

DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_DIR = os.path.join(DIR, "SBMLModel")
MODEL_DIR = os.path.join(MODEL_DIR, "data")


def doStaircase(filename, input_name=None, output_name=None, initial_value=0,
                end_time=20, final_value=10, is_plot=True, **plotKwargs):
    path = os.path.join(MODEL_DIR, filename)
    roadrunner = makeRoadrunner(path)
    if (input_name is None) or (output_name is None):
        floating_species = roadrunner.getFloatingSpeciesIds()
    if input_name is None:
        input_name = floating_species[0]
    if output_name is None:
        output_name = floating_species[1]  
    names = [input_name, output_name]
    system = SBMLSystem(roadrunner, [input_name], [output_name], is_fixed_input_species=False)
    times = np.linspace(0, end_time, 10*end_time)
    ts = system.simulateStaircase(input_name, output_name, times=times, final_value=final_value, num_step=5, is_steady_state=False)
    inputs = ts[input_name].values
    outputs = ts[output_name].values
    times = ts.index/1000
    #for column in ts.columns:
    #   if not column in names:
    #       del ts[column]
    plt.scatter(times, outputs)
    ax = plt.gca()
    ax2 = ax.twinx()
    ax2.plot(times, inputs, color="red")
    if not is_plot:
        plt.close()

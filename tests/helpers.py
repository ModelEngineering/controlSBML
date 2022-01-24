import os
from controlSBML.control_sbml import ControlSBML


DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PATH_1 = os.path.join(DIR, "Model_antimony.ant")
TEST_PATH_2 = os.path.join(DIR, "BIOMD0000000006.xml")
#TEST_PATH_3 = os.path.join(DIR, "BIOMD0000000056.xml")
TEST_PATH_3 = os.path.join(DIR, "BIOMD56.ant")

def getControl():
    return control_sbml.ControlSBML(TEST_PATH_1)

def getControl_BIOMD6():
    return control_sbml.ControlSBML(TEST_PATH_2)

def getControl_BIOMD56():
    return control_sbml.ControlSBML(TEST_PATH_3)

"""Constants for Project."""
import collections
import os


################ DIRECTORIES #################
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
for _ in range(1):
  PROJECT_DIR = os.path.dirname(PROJECT_DIR)
CODE_DIR = os.path.join(PROJECT_DIR, "controlSBML")
TEST_DIR = os.path.join(PROJECT_DIR, "tests")
DATA_DIR = os.path.join(PROJECT_DIR, "data")
BIOMODELS_ZIP_FILENAME = "biomodels.zip"

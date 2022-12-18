"""Finds the version in pyproject.toml"""

import tomli
import os

DIR = os.path.dirname(os.path.abspath(__file__))
DIR = os.path.dirname(DIR)
PATH = os.path.join(DIR, "pyproject.toml")

with open(PATH, "rb") as f:
    data = tomli.load(f)

__version__ = data["project"]["version"]

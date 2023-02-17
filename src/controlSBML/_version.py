"""Finds the version in pyproject.toml"""

import tomli
import os

DIR = os.path.dirname(os.path.abspath(__file__))
PATH = os.path.join(DIR, "pyproject.toml")

try:
    with open(PATH, "rb") as f:
        data = tomli.load(f)
    __version__ = data["project"]["version"]
except FileNotFoundError:
    __version__ = "1.0.10"


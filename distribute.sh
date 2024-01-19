#!/bin/bash
# $HOME/.pypirc contains API token information
# To update the API token, log into PyPI>Account settings, look for
#   API tokens
python3 -m build
python3 -m twine upload dist/*

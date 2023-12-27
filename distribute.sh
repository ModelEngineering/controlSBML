#!/bin/bash
# $HOME/.pypirc contains API token information
# To update the API token, log into PyPI>Account settings, look for
#   API tokens
python -m build
python -m twine upload dist/*

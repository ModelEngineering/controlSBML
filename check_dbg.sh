#!/bin/bash
# Checks if there are debug codes present
for f in tests/*.py
  do
    echo "**$f"
    grep "IGNORE_TEST = T" $f
    grep "IS_PLOT = T" $f
  done

#!/bin/bash
# Ensures that notebooks can run

# Create the script files
cd notebooks
for f in *.ipynb
  do
    jupyter nbconvert --to script $f
  done
cd ..
# Execute the script files
for f in `ls notebooks/*.py`
  do
    echo "Running $f"
    python $f
  done

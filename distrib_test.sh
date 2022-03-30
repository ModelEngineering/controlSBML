#!/bin/bash
# Setup the virtual environment to test a PyPI distribution
cd $HOME
DIR=testing_control_sbml
if [ -d $DIR ] 
then
    rm -rf $DIR
    echo "Deleting existing $DIR"
fi
python3 -m venv ${DIR}
source ${DIR}/bin/activate
pip install --upgrade pip
pip install controlSBML
echo "Testing the install"
nosetests controlSBML/tests

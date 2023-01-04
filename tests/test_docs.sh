#!/bin/bash
# Tests the code in documents
# Usage: bash test_docs.sh
REPOPATH=`get_repo_path.sh`
TESTFILE=test_docs.py
cd ${REPOPATH}/controlSBML
source activate.sh
cd docs/source
if [ -f "$TESTFILE" ]; then
    rm ${TESTFILE}
fi
for f in installation.rst system_models.rst modeling_closed_loops.rst
    do
    grep-line.sh $f >> ${TESTFILE}
done
#
python ${TESTFILE}
cp test_docs.py /tmp
rm test_docs.py

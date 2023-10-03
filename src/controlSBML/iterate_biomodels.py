"""Iterator for Biomdels SBML files."""

import controlSBML.constants as cn

import os
import zipfile

ARCHIVE_PATH = os.path.dirname(os.path.abspath(__file__))
ARCHIVE_PATH = os.path.join(ARCHIVE_PATH, "temp_biomodels.zip")
# The following models have defects and should be ignored.
IGNORE_FILES = ["BIOMD0000000075.xml", "BIOMD0000000081.xml", "BIOMD0000000353.xml",
                  "BIOMD0000000573.xml",
                  "BIOMD0000000627.xml"]

def iterateBiomodels(start=0, end=1000000, is_report=False):
    """Iterate over Biomodels SBML files.
    Args:
        start_num: int
        end_num: int
        is_report: bool (report progress)

    Yields
    ------
    str - name of SBML file
    str - contents of SBML file
    """
    with zipfile.ZipFile(ARCHIVE_PATH, "r") as archive:
        for name in archive.namelist():
            if not name.endswith(".xml"):
                continue
            if name in IGNORE_FILES:
                if is_report:
                    print("Ignoring %s" % name)
                continue
            model_num = name.split(".")[0]
            model_num = int(model_num[5:])
            if model_num < start:
                if is_report:
                    print("Skipping %s" % name)
                continue
            if model_num > end:
                if is_report:
                    print("Skipping %s" % name)
                break
            print("Acquiring %s" % name)
            contents = archive.read(name)
            yield name, contents.decode()
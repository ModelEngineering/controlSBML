"""Iterator for Biomdels SBML files."""

import controlSBML.constants as cn

import os
import zipfile

ARCHIVE_PATH = os.path.dirname(os.path.abspath(__file__))
ARCHIVE_PATH = os.path.join(ARCHIVE_PATH, "temp_biomodels.zip")
# The following models have defects and should be ignored.
IGNORE_FILES = ["BIOMD0000000056",
                "BIOMD0000000075",
                "BIOMD0000000081",
                "BIOMD0000000255",
                "BIOMD0000000353",
                "BIOMD0000000437",
                "BIOMD0000001061",
                "BIOMD0000001064",
                 ]

def iterateBiomodels(start=0, end=1000000, is_report=False, checkerFunctions=None):
    """Iterate over Biomodels SBML files.
    Args:
        start_num: int
        end_num: int
        is_report: bool (report progress)
        checkerFunctions: list-Function
            input: filename, contents
            output: str (if len > 0, then skip the file)

    Yields
    ------
    str - name of SBML file
    str - contents of SBML file
    """
    if checkerFunctions is None:
        checkerFunctions = []
    with zipfile.ZipFile(ARCHIVE_PATH, "r") as archive:
        for name in archive.namelist():
            filename, extension = os.path.splitext(name)
            if not extension == ".xml":
                continue
            if filename in IGNORE_FILES:
                if is_report:
                    print("Ignoring %s" % name)
                continue
            model_num = int(filename[5:])
            if model_num < start:
                if is_report:
                    print("Skipping %s" % name)
                continue
            if model_num > end:
                if is_report:
                    print("Skipping %s" % name)
                break
            # Get the contents
            contents = archive.read(name)
            contents = contents.decode()
            # Handle checker functions
            failed_check = False
            for func in checkerFunctions:
                result = func(name, contents)
                if len(result) > 0:
                    failed_check = True
                    if is_report:
                        print("%s: %s" % (result, name))
                    break 
            if failed_check:
                continue
            print("Processing %s" % name)
            yield name, contents
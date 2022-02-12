import numpy as np
import pandas as pd
import urllib.request


REPO_URL =  "https://github.com/ModelEngineering/controlSBML/raw/main"
MODULE_URLPAT =  REPO_URL + "/controlSBML/%s.py"  # Pattern for module URL
MODEL_URLPAT =  REPO_URL + "/models/%s"
MODEL_823_FILE = "biomodels_823.ant"
LOCAL_FILE = "local.txt"


def calculateMatrixDistance(mat1, mat2):
    """
    Calculates the distance between two matrices with the same shape.

    Parameters
    ----------
    mat1 - np.array
    mat2 - np.array
    
    Returns
    -------
    float
    """
    if np.shape(mat1) != np.shape(mat2):
        raise ValueError("Matrices must be the same shape.")
    return np.sum( (mat1 - mat2)**2)

def getSharedCodes(module_name="util"):
    """
    Obtains common codes from the github repository.

    Parameters
    ----------
    module_name: str
        name of the python module in the src directory
    """
    url = "https://github.com/ModelEngineering/controlSBML/raw/main/controlSBML/%s.py" % module_name
    local_python = "python.py"
    _, _ = urllib.request.urlretrieve(url=url, filename=local_python)
    with open(local_python, "r") as fd:
        code_str = "".join(fd.readlines())
    print(code_str)
    exec(code_str, globals())

def getModel(file_name=MODEL_823_FILE):
    """
    
    Parameters
    ----------
    file_name: str (filename with extension)
    
    Returns
    -------
    str
    """
    url = MODEL_URLPAT % file_name
    _, _ = urllib.request.urlretrieve(url=url, filename=LOCAL_FILE)
    with open(LOCAL_FILE, "r") as fd:
        model_str = "".join(fd.readlines())
    return model_str

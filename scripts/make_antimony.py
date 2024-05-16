""" Iterates across all models in BioModels to do staircase and make a SISOClosedLoop"""
import controlSBML.constants as cn  # type: ignore
import os
import tellurium as te # type: ignore

def makeAntimony(url:str, file_name:str):
    """Writes the antimony model to a file.

    Args:
        url (str)
        path (str)
    """
    rr = te.loadSBMLModel(url)
    path = os.path.join(cn.MODEL_DIR, file_name)
    antimony_str = rr.getAntimony()
    with open(path, 'w') as f:
        f.write(antimony_str)

if __name__ == '__main__':
    if False:
        makeAntimony( "https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1808280007/6/Smith2011_V1.xml",
                    "Smith2011_V1.ant")
    makeAntimony("https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1911110001/4/FatehiChenar2018.xml",
                    "FatehiChenar2018.ant")
    makeAntimony("https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1909250003/2/Varusai2018.xml",
                    "Varusai2018.ant")
    makeAntimony("https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1809060006/5/Tsai2014.xml",
                    "Tsai2014.ant")
    makeAntimony("https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1501300000/3/BIOMD0000000571_url.xml",
                    "BIOMD0000000571_url.ant")
    makeAntimony("https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL2108260003/3/Alharbi2019%20TNVM.xml",
                    "Alharbi2019_TNVM.ant")
from controlSBML.make_roadrunner import makeRoadrunner

import os
import unittest
import tellurium as te


IGNORE_TEST = False
IS_PLOT = False
HTTP_FILE = "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000206.2?filename=BIOMD0000000206_url.xml"
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
ANTIMONY_FILE = os.path.join(TEST_DIR, "Model_antimony.ant")
XML_FILE_56 = os.path.join(TEST_DIR, "BIOMD0000000056.xml")
LINEAR_MDL = """
$S0 -> S1; $S0
S1 -> S2; S1
S2 -> S3; S2

S0 = 1
S1 = 10
S2 = 0
S3 = 0
"""


#############################
# Tests
#############################
class TestMakeRoadrunner(unittest.TestCase):

    def setUp(self):
      pass

    def testAntimonyFile(self):
        if IGNORE_TEST:
            return
        rr = makeRoadrunner(ANTIMONY_FILE)
        self.assertTrue("RoadRunner" in str(type(rr)))

    def testHTTPFile(self):
        if IGNORE_TEST:
            return
        rr = makeRoadrunner(HTTP_FILE)
        self.assertTrue("RoadRunner" in str(type(rr)))

    def testAntimonyString(self):
        if IGNORE_TEST:
            return
        rr = makeRoadrunner(LINEAR_MDL)
        self.assertTrue("RoadRunner" in str(type(rr)))

    def testXMLFile(self):
        if IGNORE_TEST:
            return
        rr = makeRoadrunner(XML_FILE_56)
        self.assertTrue("RoadRunner" in str(type(rr)))

    def testXMLString(self):
        if IGNORE_TEST:
            return
        with open(XML_FILE_56, "r") as fd:
            lines = fd.readlines()
        model_str = "".join(lines)
        rr = makeRoadrunner(model_str)
        self.assertTrue("RoadRunner" in str(type(rr)))



if __name__ == '__main__':
  unittest.main()

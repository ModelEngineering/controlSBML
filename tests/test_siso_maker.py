from controlSBML.siso_maker import SISOMaker

import unittest


IGNORE_TEST = False
IS_PLOT = False
IMPROPER_LINEAR_MDL = """
// Illustrate Antimony File


S1 -> S2; k1*S1
J1: S2 -> S3; k2*S2
J2: S3 -> S2; k3*S3
J3: S2 -> S4; k4*S2

k1 = 1
k2 = 1
k3 = 1
k4 = 1
S1 = 10
S2 = 0
S3 = 0
"""
LINEAR_MDL = "model *main_model()\n" + IMPROPER_LINEAR_MDL + "\nend"
SPECIES_NAMES = ["S1", "S2", "S3"]
FILENAME = "linear.xml"


#############################
# Tests
#############################
class TestModelMaker(unittest.TestCase):

    def setUp(self):
        self.maker = SISOMaker(LINEAR_MDL, model_id="linear")

    def testConstructor(self):
        #if IGNORE_TEST:
        #    return
        self.assertTrue("RoadRunner" in str(type(self.maker.roadrunner)))
        self.assertEqual(self.maker.input_name, "S1")
        self.assertEqual(self.maker.output_name, "S2")
        self.assertTrue("SBMLSystem" in str(type(self.maker.system)))

    def testMakeStaircase(self):
        if IGNORE_TEST:
            return
        self.maker.makeStaircase(is_plot=IS_PLOT)

    def testMakeClosedLoop(self):
        if IGNORE_TEST:
            return
        self.maker.makeClosedLoop(is_plot=IS_PLOT)

    def testRunModel(self):
        if IGNORE_TEST:
            return
        self.maker.runModel(LINEAR_MDL, is_plot=IS_PLOT)

    def testRunBiomodels(self):
        #if IGNORE_TEST:
        #    return
        self.maker.runBiomodels(start=0, end=5, is_report=True, end_time=20)


if __name__ == '__main__':
  unittest.main()
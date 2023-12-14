from controlSBML.dict_array import DictArray

import unittest


IGNORE_TEST = False
DICTS = [{"kp": 1}, {"ki": 2, "kp": 3}, {"kf": 4, "kp": 5}]

#############################
# Tests
#############################
class TestDictList(unittest.TestCase):

    def setUp(self):
        self.dl = DictArray.makeFromDicts(DICTS)

    def check(self, dl=None):
        if dl is None:
            dl = self.dl
        keys = list(dl.keys())
        first_key = keys[0]
        length = len(dl[first_key])
        for value in dl.values():
            self.assertTrue(len(value) == length)

    def testConstructor(self):
        if IGNORE_TEST:
            return
        self.assertTrue(isinstance(self.dl, DictArray))
        self.check()
        
    def testMakeFromDicts(self):
        if IGNORE_TEST:
            return
        self.check()

    def testCopy(self):
        if IGNORE_TEST:
            return
        dl = self.dl.copy()
        self.assertTrue(isinstance(dl, DictArray))
        self.check(dl=dl)
        self.assertTrue(dl == self.dl)

    def testAppend(self):
        if IGNORE_TEST:
            return
        dl = self.dl.copy()
        dl.append(kp=6)
        self.check(dl=dl)

    def testMakeDicts(self):
        if IGNORE_TEST:
            return
        dcts = self.dl.makeDicts()
        dl = self.dl.makeFromDicts(dcts)
        self.assertEqual(dl, self.dl)

    def testMakeFromDataframe(self):
        if IGNORE_TEST:
            return
        df = self.dl.getDataframe()
        dl = DictArray.makeFromDataframe(df)
        self.assertEqual(dl, self.dl)

    def testEquals(self):
        if IGNORE_TEST:
            return
        dl = self.dl.copy()
        self.assertTrue(dl == self.dl)
        dl.append(kp=6)
        self.assertFalse(dl == self.dl)



if __name__ == '__main__':
    unittest.main()
from controlSBML.option_management.options import Options

import unittest

IGNORE_TEST = False
IS_PLOT = False
PLOT_DCT = dict(
      ylim=None,           # maximum and minimum value of y
      xlabel="",           
      ylabel="",           
      title="",             # plot title
      legend_spec=None,     # LegendSpec
      ax=None,              # axis to plot
      xticklabels=None,
      yticklabels=None,
      )
OTHER_DCT = dict(
      first=1,
      last=2,
      )


class TestOptions(unittest.TestCase):

    def setUp(self):
        self.options = Options(PLOT_DCT, [PLOT_DCT, OTHER_DCT])

    def testConstructor(self):
        if IGNORE_TEST:
          return
        self.assertTrue("xlabel" in self.options.keys())
        #
        options = Options(self.options, [PLOT_DCT, OTHER_DCT])
        self.assertTrue("xlabel" in options.keys())

    def testSet(self):
        if IGNORE_TEST:
          return
        DEFAULT = (1, 2)
        OVERRIDE = (1, 3)
        # Set a default value
        self.options.set("ylim", default=DEFAULT)
        self.assertEqual(self.options["ylim"], DEFAULT)
        # Override the default
        self.options.set("ylim", override=OVERRIDE)
        self.assertEqual(self.options["ylim"], OVERRIDE)
        # Setting a default value has no effect
        self.options.set("ylim", default=DEFAULT)
        self.assertEqual(self.options["ylim"], OVERRIDE)

    def testParse(self):
        if IGNORE_TEST:
          return
        def isSameDct(dct1, dct2):
            diff = set(dct1.keys()).symmetric_difference(dct2.keys())
            self.assertEqual(len(diff), 0)
            diff = set(dct1.values()).symmetric_difference(dct2.values())
            self.assertEqual(len(diff), 0)
        #
        default_dcts = [PLOT_DCT, OTHER_DCT]
        options = Options({}, default_dcts)
        plot_opts, other_opts = options.parse()
        isSameDct(plot_opts, PLOT_DCT)
        isSameDct(other_opts, OTHER_DCT)


if __name__ == '__main__':
  unittest.main()

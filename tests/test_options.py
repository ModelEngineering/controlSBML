from controlSBML.options import Options

import unittest

IGNORE_TEST = False
IS_PLOT = False
PLOT_OPTS = dict(
      ylim=None,           # maximum and minimum value of y
      xlabel="",           
      ylabel="",           
      title="",             # plot title
      legend_spec=None,     # LegendSpec
      ax=None,              # axis to plot
      xticklabels=None,
      yticklabels=None,
      )
OTHER_OPTS = dict(
      first=1,
      last=2,
      )


class TestOptions(unittest.TestCase):

    def setUp(self):
        self.options = Options(PLOT_OPTS)

    def testConstructor(self):
        if IGNORE_TEST:
          return
        self.assertTrue("xlabel" in self.options.keys())
        #
        options = Options(self.options)
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
        kwargs = Options(PLOT_OPTS)
        kwargs.update(OTHER_OPTS)
        plot_opts, other_opts = self.options.parse([PLOT_OPTS, OTHER_OPTS])
        isSameDct(plot_opts, PLOT_OPTS)
        isSameDct(other_opts, OTHER_OPTS)


if __name__ == '__main__':
  unittest.main()

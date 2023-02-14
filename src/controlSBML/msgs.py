"""Generate messages for the application."""

import sys
import warnings

def error(text):
  print("***Error. Reason follows.")
  print("   %s" % text)
  sys.exit()

def warn(text):
  new_text = "\n\n***Warning*** %s" % text
  warnings.warn(new_text)

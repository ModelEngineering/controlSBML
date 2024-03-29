"""Generate messages for the application."""

import sys
import warnings

def error(text):
  print("***Error. Reason follows.")
  print("   %s" % text)
  raise ValueError("Cannot continue")

def warn(text):
  new_text = "\n\n***Warning*** %s" % text
  warnings.warn(new_text)

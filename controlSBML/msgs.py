"""Generate messages for the application."""

import sys

def error(text):
  print("***Error. Reason follows.")
  print("   %s" % text)
  sys.exit()

def warn(text):
  print("\n\n***Warning***\n")
  print("   %s" % text)

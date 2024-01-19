# NOTES FOR BUILDING
1. Use distribute.sh to create a new release. Update versions in _version.py and pyproject.toml
2. Testing the release
  1. Deactivate existing virtual environment
  1. Create a new virtual environment
  1. Eliminate PYTHONPATH
  1. Copy tests outside of controlSBML
  1. nose2 tests

from setuptools import setup, find_packages
import os.path
import codecs

with open("README.md", "r") as fh:
    long_description = fh.read()

# The following two methods were copied from
# https://packaging.python.org/guides/single-sourcing-package-version/#single-sourcing-the-version
def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

# parse_requirements() returns generator of pip.req.InstallRequirement objects
INSTALL_REQUIREMENTS = [
    "coverage",
    "control",
    "docstring-expander",
    "jupyterlab",
    "matplotlib",
    "nose",
    "numpy",
    "pandas",
    "pip",
    "pylint",
    "python-libsbml",
    "seaborn",
    "symmathSBML",
    "sympy",
    "tellurium",
    "tk",
    "urllib3",
]


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            parts = line.split(" ")
            version = parts[2]
            version = version.replace('"', '')
            version = version.replace("'", '')
            print ("Processing version %s" % version)
            return version
    else:
        raise RuntimeError("Unable to find version string.")

def doSetup(install_requires):
  try:
    setup(
        name='controlSBML',
        version=get_version("controlSBML/_version.py"),
        author='Joseph L. Hellerstein',
        author_email='jlheller@uw.edu',
        packages=find_packages(exclude=['tests', 'sym', 'images', 'data',
            'notebook', 'docs']),
        url='https://github.com/ModelEngineering/controlSBML',
        description='Control analysis and design for SBML Models.',
        long_description=open('README.md').read(),
        long_description_content_type='text/markdown',
        package_dir={'controlSBML': 'controlSBML'},
        python_requires='>=3.6',
        install_requires=install_requires,
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Developers',      # Define that your audience are developers
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'License :: OSI Approved :: MIT License',   # Again, pick a license
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
          ],
        )
  except Exception as exp:
    import pdb; pdb.set_trace()
    pass


if __name__ == '__main__':
  doSetup(INSTALL_REQUIREMENTS)

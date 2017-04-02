from codecs import open as codecs_open
from setuptools import setup, find_packages
from warnings import warn


# Get the long description from the relevant file
with codecs_open('README.md', encoding='utf-8') as f:
    long_description = f.read()

HAVE_SIMUPOP = False
HAVE_MSPRIME = False
try:
    import msprime
except ImportError:
    warn("`msprime` not present and must be installed")


setup(name='ftprime',
      version='0.0.1',
      description=u"fill msprime data structure in forward time",
      long_description=long_description,
      classifiers=[],
      keywords='',
      author=u"Jaime Ashander",
      author_email='jashander@ucdavis.edu',
      url='https://github.com/ashander/ftprime',
      license='GPL3',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          'click'
      ],
      extras_require={
          'dev': [
              'pytest',
              'pytest-cov',
              'sphinx',
              'recommonmark',
              'sphinx_rtd_theme',
              'simuPOP>=1.1.7'
          ],
      },
      )

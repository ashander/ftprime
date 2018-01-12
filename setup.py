from codecs import open as codecs_open
from setuptools import setup, find_packages
from warnings import warn


# Get the long description from the relevant file
with codecs_open('README.md', encoding='utf-8') as f:
    long_description = f.read()

try:
    import msprime  #NOQA
except ImportError:
    warn("`msprime` not present and must be installed")
try:
    import simuPOP  #NOQA
except ImportError:
    warn("`simuPOP` not present and must be installed")


setup(name='ftprime',
      version='0.0.6-rc',
      description=u"Simulate a `msprime` tree sequence in forward time",
      long_description=long_description,
      keywords=['tree sequence'],
      author=u"Jaime Ashander",
      author_email='jashander@ucla.edu',
      url='https://github.com/ashander/ftprime',
      license='GPL3',
      packages=find_packages(exclude=['writeups', 'ez_setup', 'examples',
                                      'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          'msprime',
          'numpy',
      ],
      extras_require={
          'dev': [
              'pytest',
              'pytest-cov',
              'sphinx',
              'recommonmark',
              'sphinx_rtd_theme',
              'six',
          ],
          'recomb_collector': [
              'simuPOP>=1.1.8.3',
          ]
      },
      platforms=["POSIX"],
      classifiers=['Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU General Public License v3 ' +
                   'or later (GPLv3+)'
                   'Operating System :: POSIX',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Programming Language :: Python :: 3.5',
                   'Programming Language :: Python :: 3.6',
                   'Development Status :: 3 - Alpha',
                   'Environment :: Other Environment',
                   'Intended Audience :: Science/Research',
                   ],
      )

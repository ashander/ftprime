from codecs import open as codecs_open
from setuptools import setup, find_packages


# Get the long description from the relevant file
with codecs_open('README.rst', encoding='utf-8') as f:
    long_description = f.read()


setup(name='ftprime',
      version='0.0.1',
      description=u"msprime in forward time",
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
          'test': ['pytest'],
      },
      )

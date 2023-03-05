#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/5/23 4:52 PM

@author: alejandrobertolet
"""

import io, os, sys
from shutil import rmtree

from distutils.core import setup
from setuptools import find_packages, setup, Command

# Package meta-data
NAME = 'mgm'
DESCRIPTION = 'Microdosimetric Gamma Model'
URL = 'https://github.com/mghro/mgm'
EMAIL = 'abertoletreina@mgh.harvard.edu'
AUTHOR = 'Alejandro Bertolet'
REQUIRES_PYTHON = '>=3.6.0'
VERSION = '1.0.0'

# What packages are required for this module to be executed?
REQUIRED = ['numpy', 'scipy', 'matplotlib']

# What packages are optional?
EXTRAS = { }

# Import the README and use it as the long-description
try:
    with io.open('README.md', encoding='utf-8') as f:
        LONG_DESCRIPTION = 'n' + f.read()
except:
    LONG_DESCRIPTION = DESCRIPTION

# Load the package's __version__.py module as a dictionary
about = {}
if not VERSION:
    project_slug = NAME.lower().replace("-", "_").replace(" ", "_")
    with open(os.path.join(project_slug, '__version__.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = VERSION

# Support setup.py upload
class UploadCommand(Command):
    description = 'Build and publish the package.'
    user_options = []

    @staticmethod
    def status(s):
        print('n', s)

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status('Removing previous builds...')
            rmtree(os.path.join(os.path.dirname(__file__), 'dist'))
        except OSError:
            pass

        self.status('Building Source and Wheel (universal) distribution...')
        os.system('{0} setup.py sdist bdist_wheel --universal'.format(sys.executable))

        self.status('Uploading the package to PyPI via Twine...')
        os.system('twine upload dist/*')

        self.status('Pushing git tags...')
        os.system('git tag v{0}'.format(about['__version__']))
        os.system('git push --tags')

# Execute setup
setup(
    name=NAME,
    version=about['__version__'],
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    #include_package_data=True,
    license='MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    # $ setup.py publish support.
    cmdclass={
        'upload': UploadCommand,
    },
)

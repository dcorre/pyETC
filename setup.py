#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from __future__ import print_function

import os, sys
from distutils.version import LooseVersion

from setuptools import (setup, find_packages,
                        __version__ as setuptools_version)


#import versioneer
##from setup_utils import (CMDCLASS, get_setup_requires, get_scripts)
#__version__ = versioneer.get_version()
#CMDCLASS=versioneer.get_cmdclass()



with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('docs/pyETC/history.rst') as history_file:
    history = history_file.read()


# Cython is required by scikit-image
install_requires = [
        'Click>=6.0',
        'cython',
        'matplotlib',
        'scipy',
        'astropy',
        'jupyter',
        'scikit-image',
        'hjson'
        ]


setup_requirements = ['pytest-runner', 'numpy']

test_requirements = ['pytest', ]

setup(
    author="David Corre, Alain Klotz",
    author_email='david.corre.fr@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    description="Exposure Time Calculator for optical/NIR telescope",
    entry_points={
        'console_scripts': [
            'pyETC=pyETC.cli:main',
        ],
    },
    install_requires=install_requires,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='pyETC',
    name='pyETC',
    packages=find_packages(),
    #scripts=get_scripts(),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/dcorre/pyETC',
    version='0.1.0',
    zip_safe=False,
)

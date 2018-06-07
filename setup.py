#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('docs/pyETC/history.rst') as history_file:
    history = history_file.read()

requirements = ['Click>=6.0', 'numpy', 'matplotlib', 'scipy', 'astropy',
                'jupyter', 'scikit-image', 'hjson']

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author="David Corre",
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
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='pyETC',
    name='pyETC',
    packages=find_packages(include=['pyETC']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/dcorre/pyETC',
    version='0.1.0',
    zip_safe=False,
)

# Exposure Time Calculator for optical/NIR telescope

* Free software: MIT license
* Documentation: https://pyETC.readthedocs.io.

Release status
--------------

[![PyPI version](https://badge.fury.io/py/pyETC.svg)](http://badge.fury.io/py/pyETC)
![Supported Python versions](https://img.shields.io/pypi/pyversions/pyETC.svg)


Development status
--------------------

[![Build Status](https://travis-ci.com/dcorre/pyETC.svg?branch=master)](https://travis-ci.com/dcorre/pyETC)
[![codecov](https://codecov.io/gh/dcorre/pyETC/branch/master/graphs/badge.svg)](https://codecov.io/gh/dcorre/pyETC/branch/master)
[![Documentation Status](https://readthedocs.org/projects/pyetc/badge/?version=latest)](https://pyetc.readthedocs.io/en/latest/?badge=latest)

Features
--------
* Exposure Time Calculator configurable for any telescope in the optical/NIR domain.
* Telescope configuration and response defined by transmission curve of each optical element (mirrors, lenses, coating, filters, dichroics) and camera quntum efficiency. Also possible to provide only a global transmission curve.
* Object can be a givn magnitude (AB or Vega), a file containing a spectrum or a light curve. Can also be a simulated GRB light curve through the connection with [pyGRBaglow](https://github.com/dcorre/pyGRBaglow).
* Output is provided through a python dictionary.
* Have a look to the notebook(s) to see how it works.


Installation
------------
See the doc: https://pyETC.readthedocs.io.

Credits
-------

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [audreyr/cookiecutter-pypackage](https://github.com/audreyr/cookiecutter-pypackage) project template.


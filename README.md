<!-- # Release status

[![PyPI version](https://img.shields.io/pypi/v/pyETC.svg)](https://pypi.python.org/pypi/pyETC)
-->
# Development status

[![Build Status](https://travis-ci.com/dcorre/pyETC.svg?branch=master)](https://travis-ci.com/dcorre/pyETC)
[![codecov](https://codecov.io/gh/dcorre/pyETC/branch/master/graphs/badge.svg)](https://codecov.io/gh/dcorre/pyETC/branch/master)
[![Documentation Status](https://readthedocs.org/projects/pyetc/badge/?version=latest)](https://pyetc.readthedocs.io/en/latest/?badge=latest)
<!-- [![Linux](https://img.shields.io/travis/dcorre/pyETC/master.svg?label=Linux)](https://travis-ci.com/pyETC/pyETC)
[![OSX](https://img.shields.io/travis/dcorre/pyETC/master.svg?label=OSX)](https://travis-ci.com/pyETC/pyETC)
[![Windows](https://img.shields.io/travis/dcorre/pyETC/master.svg?label=Windows)](https://travis-ci.com/pyETC/pyETC)
-->
Exposure Time Calculator for optical/NIR telescope


* Free software: MIT license
* Documentation: https://pyETC.readthedocs.io.


Features
--------

* Telescope configuration and response defined by transmission curve of each optical element (mirrors, lenses, coating, filters, dichroics) and camera quntum efficiency. Also possible to provide only a global transmission curve.
* Object can be a givn magnitude (AB or Vega), a file containing a spectrum or a light curve. Can also be a simulated GRB light curve through the connection with [pyGRBaglow](https://github.com/dcorre/pyGRBaglow)].
* Output is provided through a python dictionary.

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

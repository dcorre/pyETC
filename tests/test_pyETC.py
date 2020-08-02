#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `pyETC` package."""

import pytest
import numpy as np

from click.testing import CliRunner

from pyETC import pyETC
from pyETC import cli
from pyETC.pyETC import etc

@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert result.exit_code == 0
    assert 'pyETC.cli.main' in result.output
    help_result = runner.invoke(cli.main, ['--help'])
    assert help_result.exit_code == 0
    assert '--help  Show this message and exit.' in help_result.output

def test_simple_run():
    """Test that one can make a run"""
    # load ETC with a config file and the COLIBRI caracteristics
    COLIBRI_ETC=etc(configFile='example.hjson',name_telescope='colibri')

    # Execute
    COLIBRI_ETC.sim()

    assert np.round(COLIBRI_ETC.information['mag'],3) == 9.812


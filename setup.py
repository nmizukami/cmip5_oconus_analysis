#!/usr/bin/env python

"""
The setup script
https://stackoverflow.com/questions/1471994/what-is-setup-py
https://packaging.python.org/tutorials/packaging-projects/
"""

import sys
from setuptools import setup, find_packages

# with open("requirements.txt") as f:
#     INSTALL_REQUIRES = f.read().strip().split("\n")

with open("README.md") as f:
    LONG_DESCRIPTION = f.read()

PYTHON_REQUIRES = '>=3.6'

description = ("CMIP5 OCONUS hydrologic projection analysis")
setup(
    name="cmip5_oconus",
    version='0.0.1',
    description=description,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    maintainer="Naoki Mizukami",
    maintainer_email="mizukami@ucar.edu",
    url="https://github.com/mizukami/cmip5_oconus_analysis",
    py_modules=['cmip5_oconus'],
    packages=find_packages(),
    python_requires=PYTHON_REQUIRES,
    license="Apache",
    keywords="oconus",
)

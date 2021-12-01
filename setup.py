#!/usr/bin/env python
from setuptools import setup, find_packages

DISTNAME = "lib_position_transition_zone"
VERSION = "0.0.1"
LICENSE = "MIT"
AUTHOR = "Romain Caneill"
AUTHOR_EMAIL = "romain.caneill@gu.se"
DESCRIPTION = "Library used for the analyze of the NEMO runs, regarding the position of the transition zone."

setup(
    name=DISTNAME,
    version=VERSION,
    license=LICENSE,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    long_description=DESCRIPTION,
    packages=find_packages(exclude=["docs", "tests"]),
)

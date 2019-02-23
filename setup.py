#!/usr/bin/env python

import setuptools
import pathlib

name = "mutational_signature"
version = "0.4"
release = "0.4.0"
here = pathlib.Path(__file__).parent.resolve()

setuptools.setup(
    name=name,
    version=version,
    packages=setuptools.find_packages()
)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2021  Elija Feigl
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see https://www.gnu.org/licenses/gpl-3.0.html.

from ifs_data.version import get_version
from setuptools import setup, find_packages

description = """\
Initial folding screen DNAOrigami utility dieztlab\
    processes database design files. matlabscript on IFS has to be run first.
    creates database.csv in specified folder"
"""

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("LICENSE", "r") as fh:
    license = fh.read()

setup(
    name="ifs_data",
    version=get_version(),
    author="Elija Feigl",
    author_email="elija.feigl@tum.de",
    description=description,
    license=license,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/elija-feigl/ifs_data",
    packages=find_packages(),
    include_package_data=True,
    install_requires=(
        'Bio',
        'numpy',
        'attrs',
        'click',
        'pandas',
        'scipy',
    ),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: GNU General Public License Version 3",
        "Operating System :: OS Independent",
    ),
    entry_points='''
        [console_scripts]
        ifs=ifs_data.scripts.ifs:cli
    ''',
)

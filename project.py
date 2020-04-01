#!/usr/bin/env python
# -*- coding: utf-8 -*-3

import attr
from pathlib import Path
import contextlib


@attr.s(slots=True)
class Project(object):
    input: Path = attr.ib()
    output: Path = attr.ib()
    filename: str = attr.ib()


@contextlib.contextmanager
def ignored(*exceptions):
    try:
        yield
    except exceptions:
        pass


@contextlib.contextmanager
def get_file(logger, folder: Path, rex: str, *exceptions):
    try:
        yield
    except exceptions as e:
        e_ = "{} missing    | folder: {}".format(rex[-4:], folder)
        logger.error(e_)

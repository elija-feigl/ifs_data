#!/usr/bin/env python
# -*- coding: utf-8 -*-3

import numpy as np
import attr
from pathlib import Path
import contextlib


EXC_TXT = "______________| Folder: {}\n     | Exception: {}"
GEL_PROPERTIES = ["user", "project", "design_name", "date",
                  "scaffold_type", "lattice_type", "scaffold_concentration",
                  "staple_concentration", "gelsize", "agarose_concentration",
                  "staining", "mg_concentration", "voltage", "running_time",
                  "cooling"]
FOLD_PROPERTIES = ["qualityMetric", "bestTscrn", "bestMgscrn", "qualityMetric",
                   "fractionMonomer", "bandWidthNormalized", "migrationDistanceNormalized", "fractionPocket", "fractionSmear"]


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
    except exceptions:
        e_ = ("{} missing  " + EXC_TXT[14:26]).format(rex[-4:], folder)
        logger.error(e_)


def get_statistics(data_list, data_name):
    """[summary]

    Arguments:
        data_list {[type]} -- [description]
        data_name {[type]} -- [description]

    Returns:
        [type] -- [description]
    """
    return {data_name + "_avg": np.average(data_list),
            data_name + "_std": np.std(data_list),
            data_name + "_max": np.max(data_list),
            data_name + "_min": np.min(data_list),
            }

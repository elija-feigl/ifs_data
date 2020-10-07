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
    if len(data_list) != 0:
        return {"avg_" + data_name: np.average(data_list),
                "std_" + data_name: np.std(data_list),
                "max_" + data_name: np.max(data_list),
                "min_" + data_name: np.min(data_list),
                }
    else:
        return {"avg_" + data_name: 0,
                "std_" + data_name: 0,
                "max_" + data_name: 0,
                "min_" + data_name: 0,
                }


def get_full_scaff_co_typ_stat(design):
    """[numbers of ]

    Returns:
        [dict]: [data for database; writung to csv]
    """
    # TODO: the designprocess data are not consistant
    data = {
        'full_scaf_co_type_1': 0,
        'full_scaf_co_type_2': 0,
        'full_scaf_co_type_3': 0
    }
    for full in design.full_crossovers:
        if full.strand_typ == 'scaffold':
            if full.scaff_full_type == 1:
                data['full_scaf_co_type_1'] += 1
            elif full.scaff_full_type == 2:
                data['full_scaf_co_type_2'] += 1
            elif full.scaff_full_type == 3:
                data['full_scaf_co_type_3'] += 1
    return data

import numpy as np
from pathlib import Path
import contextlib

from typing import Tuple

from nanodesign.data.base import DnaBase as Base
from nanodesign.data.strand import DnaStrand as Strand

EXC_TXT = "______________| Folder: {}\n     | Exception: {}"
GEL_PROPERTIES = ["user", "project", "design_name", "date", "tem_verified",
                  "scaffold_type", "lattice_type", "scaffold_concentration",
                  "staple_concentration", "gelsize", "agarose_concentration",
                  "staining", "mg_concentration", "voltage", "running_time",
                  "cooling"]
FOLD_PROPERTIES = ["qualityMetric", "bestTscrn", "bestMgscrn", "qualityMetric",
                   "fractionMonomer", "bandWidthNormalized", "migrationDistanceNormalized", "fractionPocket", "fractionSmear"]
scaffold_dict_len = {"8064": 8064, "7560": 7560, "cs11": 7560, "cs16": 7560,
                     "cs15": 7560, "2873": 2873, "1317": 1317, "2048": 2057,
                     "7249": 7249, "9072": 9072, "cs12": 7560,
                     "7704": 7704, "cs17": 8039, "cs13": 7560, "4536": 4536}
scaffold_dict_name = {"8064": "8064", "7560": "7560", "cs11": "cs11", "cs16": "cs16",
                      "cs15": "cs15", "2873": "CS3_XS", "1317": "RFP", "2048": "Pippin",
                      "7249": "7249", "9072": "CS3_XL", "cs12": "cs12",
                      "7704": "7704", "cs17": "cs17", "cs13": "cs13", "4536": "CS3_S"}
scaffold_dict_circ = {"8064": True, "7560": True, "cs11": True, "cs16": False,
                      "cs15": False, "2873": True, "1317": True, "2048": True,
                      "7249": True, "9072": True, "cs12": False,
                      "7704": True, "cs17": False, "cs13": False, "4536": True}
scaffold_dict_gc = {"8064": 0.44, "7560": 0.42, "cs11": 0.435, "cs16": 0.6,
                    "cs15": 0.58, "2873": 0.5, "1317": 0.51, "2048": 0.51,
                    "7249": 0.42, "9072": 0.49, "cs12": 0.3,
                    "7704": 0.427, "cs17": 0.448, "cs13": 0.34, "4536": 0.49}



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
    if data_list:
        return {"avg_" + data_name: np.average(data_list),
                "std_" + data_name: np.std(data_list),
                "max_" + data_name: np.max(data_list),
                "min_" + data_name: np.min(data_list),
                }
    else:
        return {"avg_" + data_name: 0.,
                "std_" + data_name: 0.,
                "max_" + data_name: 0.,
                "min_" + data_name: 0.,
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


def _hps(base: Base) -> Tuple[int, int, bool]:
    return (base.h, base.p, base.is_scaf)


def _close_strand(strand: Strand) -> None:
    """ closes a given strand, making it a loop."""
    start = strand.tour[0]
    end = strand.tour[-1]
    start.up, end.down = end, start
    strand.is_circular = True

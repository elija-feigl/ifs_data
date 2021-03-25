import numpy as np
from pathlib import Path

from typing import Tuple, Union, Any

from nanodesign.data.base import DnaBase as Base
from nanodesign.data.strand import DnaStrand as Strand

EXC_TXT = "______________| Folder: {}\n     | Exception: {}"
GEL_PROPERTIES = ["user", "project", "design_name", "date", "tem_verified",
                  "scaffold_type", "lattice_type", "scaffold_concentration",
                  "staple_concentration", "gelsize", "agarose_concentration",
                  "staining", "mg_concentration", "voltage", "running_time",
                  "cooling"]
FOLD_PROPERTIES = ["qualityMetric", "bestTscrn", "bestMgscrn", "qualityMetric",
                   "fractionMonomer", "bandWidthNormalized", "migrationDistanceNormalized",
                   "fractionPocket", "fractionSmear"]
scaffold_dict_len = {"8064": 8064, "7560": 7560, "CS11": 7560, "CS16": 7560,
                     "CS15": 7560, "2873": 2873, "1317": 1317, "2048": 2057,
                     "7249": 7249, "9072": 9072, "CS12": 7560,
                     "7704": 7704, "CS17": 8039, "CS13": 7560, "4536": 4536}
scaffold_dict_name = {"8064": "8064", "7560": "7560", "CS11": "CS11", "CS16": "CS16",
                      "CS15": "CS15", "2873": "CS3_XS", "1317": "RFP", "2048": "Pippin",
                      "7249": "7249", "9072": "CS3_XL", "CS12": "CS12",
                      "7704": "7704", "CS17": "CS17", "CS13": "CS13", "4536": "CS3_S"}
scaffold_dict_circ = {"8064": True, "7560": True, "CS11": True, "CS16": False,
                      "CS15": False, "2873": True, "1317": True, "2048": True,
                      "7249": True, "9072": True, "CS12": False,
                      "7704": True, "CS17": False, "CS13": False, "4536": True}
scaffold_dict_gc = {"8064": 0.44, "7560": 0.42, "CS11": 0.435, "CS16": 0.6,
                    "CS15": 0.58, "2873": 0.5, "1317": 0.51, "2048": 0.51,
                    "7249": 0.42, "9072": 0.49, "CS12": 0.3,
                    "7704": 0.427, "CS17": 0.448, "CS13": 0.34, "4536": 0.49}
T_screen = {f"T{x}": list(range(47 + 2 * (x-1), 51 + 2 * (x-1)))
            for x in range(1, 9)}
Mg_screen = {f"M{x}": x for x in range(5, 35, 5)}


def get_file(logger, folder: Path, suffix: str) -> Path:
    files = list(folder.glob(f"*.{suffix}"))
    if not files:
        logger.error(f"Could not find .{suffix} file. Abort!")
        raise IOError
    elif len(files) > 1:
        logger.error(f"Found more than one .{suffix} file. Abort!")
        raise IOError
    else:
        return files.pop()


def save_division(dividend,  divisor) -> float:
    return 0. if not divisor else dividend / divisor


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


def _hps(base: Base) -> Tuple[int, int, bool]:
    return (base.h, base.p, base.is_scaf)


def _close_strand(strand: Strand) -> None:
    """ closes a given strand, making it a loop."""
    start = strand.tour[0]
    end = strand.tour[-1]
    start.up, end.down = end, start
    strand.is_circular = True


def _check_base_crossover(base: Base) -> Union[Base, bool]:
    """ Check if there is a crossover to a different helix at the given
        base.
        returns base if True else returns False
    """
    if base.down is None:
        return False
    if base.down.h != base.h:
        return base.down

    if base.up is None:
        return False
    if base.up.h != base.h:
        return base.up
    return False


def _change_strand_type(strand: Strand):
    """switch strand type. Propagates to bases of strand."""
    typ = strand.is_scaffold
    strand.is_scaffold = not typ
    for base in strand.tour:
        base.is_scaf = not typ


def _str(value: Any):
    if isinstance(value, float):
        return f'{value:.5g}'
    else:
        return str(value)

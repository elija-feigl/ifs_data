#!/usr/bin/env python
# -*- coding: utf-8 -*-3

import os
import argparse
import logging
import scipy.io as sio

from pathlib import Path
from typing import IO
import datetime

from utils import Project, ignored, get_file, EXC_TXT, GEL_PROPERTIES, FOLD_PROPERTIES
from designData import DesignData
from compute_data import Compute

__authors__ = ["Elija Feigl", "Kourosh Zargari"]
__version__ = "0.2"
__descr__ = "processes database design files. matlabscript on IFS has to be\
             run first. creates database.csv in specified folder"


def export_data(data: dict, fdb_file: IO) -> None:
    if not fdb_file.tell():  # true if empty
        header = ",".join(str(k) for k in data.keys())
        fdb_file.write(header + "\n")
    export = ",".join(str(v) for v in data.values())
    fdb_file.write(export + "\n")


def process_mat_file(mat_file: IO) -> dict:

    mat = sio.loadmat(mat_file, squeeze_me=True)
    try:
        gel_info = mat["gelInfo"]
    except KeyError:
        raise

    try:
        fold_info = mat["foldingAnalysis"]
    except KeyError:
        raise

    best_idx = int(fold_info["bestFoldingIndex"])
    data = dict()

    for prop in GEL_PROPERTIES + FOLD_PROPERTIES:
        info = gel_info if prop in GEL_PROPERTIES else fold_info
        is_set = (prop in info.dtype.names)
        # NOTE: why is this in here twice?
        if prop in ["qualityMetric", "fractionMonomer", "bandWidthNormalized",
                    "migrationDistanceNormalized", "fractionPocket", "fractionSmear"] and is_set:
            prop_str = str(info[prop].item()[best_idx])
        elif is_set:
            prop_str = str(info[prop]).replace(",", ".")
        else:
            prop_str = " "
        data.update({prop: prop_str.lower()})

    for prop in ["qualityMetric", "fractionMonomer", "bandWidthNormalized",
                 "migrationDistanceNormalized", "fractionPocket", "fractionSmear"]:
        if prop in fold_info.dtype.names:
            prop_float = fold_info[prop].item()[best_idx]
        else:
            prop_float = 0
        data.update({prop: prop_float})
    return data


def proc_input() -> Project:
    date_str = str(datetime.date.today().strftime("%y-%b-%d"))

    def get_description() -> str:
        return "{}\n {}\n {}".format(__descr__, __version__, __authors__)

    parser = argparse.ArgumentParser(
        description=get_description(),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input",
                        help="input folder",
                        type=str,
                        default="D:/work/Dropbox (DIETZ LAB)/FOLDINGSCREENS",
                        )
    parser.add_argument("-o", "--output",
                        help="output folder",
                        type=str,
                        default=f"./AAA__database__/{date_str}"
                        )
    parser.add_argument("-d", "--datafile",
                        help="database-file name",
                        type=str,
                        default="fdb"
                        )

    args = parser.parse_args()
    project = Project(input=Path(args.input),
                      output=Path(args.output),
                      filename=args.datafile
                      )
    with ignored(FileExistsError):
        os.mkdir(project.output)
    return project


def main():
    project = proc_input()
    date_str = str(datetime.date.today().strftime("%y-%b-%d"))
    log_time = str(datetime.datetime.now().strftime("%H-%M"))

    # NOTE: the log file name is the hour and minute that it has been created.
    logname = f"fdb-{log_time}.log"
    logging.basicConfig(filename=project.output / logname)

    handle = "folding-DB"
    logger = logging.getLogger(handle)

    filename = "{}-{}.csv".format(project.filename, date_str)
    fdb_filepath = project.output / filename

    with open(fdb_filepath, mode="w+", encoding="utf-8") as fdb_file:
        for child in project.input.iterdir():
            if child.name.startswith(".") or child.name[-2:] == "__":
                continue
            with get_file(logger, child, "*.json", IndexError):
                json = list(child.glob("*.json")).pop()
            with get_file(logger, child, "*.mat", IndexError):
                mat = list(child.glob("*.mat")).pop()

            try:
                mat_data = process_mat_file(mat)
            except Exception as e:
                e_ = "info.mat file " + EXC_TXT[14:].format(child.name, e)
                logger.error(e_)
                continue

            design_name = mat_data["design_name"]
            design_seq = mat_data["scaffold_type"].upper()

            try:
                designdata = DesignData(
                    json=json, name=design_name, seq=design_seq)
            except Exception as e:
                e_ = "nanodesign    " + EXC_TXT[14:].format(child.name, e)
                logger.error(e_)
                continue
            try:
                Compute.compute_data(designdata)
            except Exception as e:
                e_ = "designdata    " + EXC_TXT[14:].format(child.name, e)
                logger.error(e_)
                continue

            json_data = Compute.prep_data_for_export(designdata)
            data = {**mat_data, **json_data}

            try:
                export_data(data=data, fdb_file=fdb_file)
            except Exception as e:
                e_ = "data export   " + EXC_TXT[14:].format(child.name, e)
                logger.error(e_)
                continue

    return


if __name__ == "__main__":
    main()

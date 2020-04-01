#!/usr/bin/env python
# -*- coding: utf-8 -*-3

import os
import argparse
import logging
import csv
import scipy.io as sio

from pathlib import Path
from typing import IO
from datetime import date

from utils import Project, ignored, get_file, EXC_TXT
from designData import DesignData


__authors__ = ["Elija Feigl", "Kuorosh Zargari"]
__version__ = "0.1"
__descr__ = "processes database design files. matlabscript on IFS has to be\
             run first. creates database.csv in specified folder"


def export_data(data: dict, fdb_file: IO) -> None:
    if not fdb_file.tell():  # true if empty
        header = ", ".join(str(k) for k in data.keys())
        fdb_file.write(header + "\n")
    export = ", ".join(str(v) for v in data.values())
    fdb_file.write(export + "\n")


def process_mat_file(mat_file: IO) -> dict:
    data = dict()
    mat = sio.loadmat(mat_file, squeeze_me=True)

    types = ['user', 'project', 'design_name', 'date',
             'scaffold_type', 'lattice_type', 'scaffold_concentration',
             'staple_concentration', 'gelsize', 'agarose_concentration',
             'staining', 'mg_concentration', 'voltage', 'running_time',
             'cooling']

    for typ in types:
        try:
            if ',' in str(mat['gelInfo'][typ]):
                new = str(mat['gelInfo'][typ]).replace(',', '.')
                info = {typ: new}
            else:
                info = {typ: str(mat['gelInfo'][typ])}

            data.update(info)
        except ValueError:
            data.update({typ: " "})

    types2 = ["qualityMetric", "bestTscrn", "bestMgscrn"]

    for typ in types2:
        try:
            if ',' in str(mat['foldingAnalysis'][typ]):
                new = str(mat['foldingAnalysis'][typ]).replace(',', '.')
                info = {typ: new}
            else:
                info = {typ: str(mat['foldingAnalysis'][typ])}

            data.update(info)
        except ValueError:
            data.update({typ: " "})

    best_idx = int(mat['foldingAnalysis']["bestFoldingIndex"])

    for typ in ["qualityMetric", "fractionMonomer"]:
        try:
            new = float(mat['foldingAnalysis'][typ].tolist()[best_idx])
            info = {typ: new}
            data.update(info)
        except ValueError:
            data.update({typ: " "})

    return data


def proc_input() -> Project:
    def get_description() -> str:
        return "{}\n {}\n {}".format(__descr__, __version__, __authors__)

    parser = argparse.ArgumentParser(
        description=get_description(),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input",
                        help="input folder",
                        type=str,
                        default="./",
                        )
    parser.add_argument("-o", "--output",
                        help="output folder",
                        type=str,
                        default="./AAA__database__"
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
    logging.basicConfig()
    handle = "folding-DB"
    logger = logging.getLogger(handle)
    project = proc_input()

    date_str = str(date.today().strftime("%y-%b-%d"))
    filename = "{}-{}.csv".format(project.filename, date_str)
    fdb_filepath = project.output / filename

    with open(fdb_filepath, mode="w+") as fdb_file:
        for child in project.input.iterdir():
            if child.name.startswith('.') or child.name[-2:] == "__": continue
            with get_file(logger, child, "*.json", IndexError):
                json = list(child.glob("*.json")).pop()
            with get_file(logger, child, "*.mat", IndexError):
                mat = list(child.glob("*.mat")).pop()

            try: mat_data = process_mat_file(mat)
            except Exception as e:
                e_ = "info.mat file " + EXC_TXT[14:].format(child.name, e)
                logger.error(e_)

            design_name = mat_data["design_name"]

            try: designdata = DesignData(json=json, name=design_name)
            except Exception as e:
                e_ = "nanodesign    " + EXC_TXT[14:].format(child.name, e)
                logger.error(e_)
            try: json_data = designdata.compute_data()
            except Exception as e:
                e_ = "designdata    " + EXC_TXT[14:].format(child.name, e)
                logger.error(e_)

            json_data = designdata.prep_data_for_export()
            data = {**mat_data, **json_data}

            try: export_data(data=data, fdb_file=fdb_file)
            except Exception as e:
                e_ = "data export   " + EXC_TXT[14:].format(child.name, e)
                logger.error(e_)
    return


if __name__ == "__main__":
    main()

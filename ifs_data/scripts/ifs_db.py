#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import click
import logging
import scipy.io as sio
import pandas as pd
import numpy as np

from pathlib import Path
from typing import IO
from datetime import date

from ..version import get_version
from ..core.utils import (get_file, EXC_TXT, GEL_PROPERTIES, FOLD_PROPERTIES,
                          scaffold_dict_len, scaffold_dict_name, scaffold_dict_circ, scaffold_dict_gc)
from ..core.designData import Design
from ..database.compute_data import DesignStats


logger = logging.getLogger(__name__)


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(get_version())
    ctx.exit()


def export_data(data: dict, fdb_file: IO) -> None:
    if not fdb_file.tell():  # true if empty
        header = ",".join(str(k) for k in data.keys())
        fdb_file.write(header + "\n")
    export = ",".join(str(v) for v in data.values())
    fdb_file.write(export + "\n")


def process_txt_file(txt_file: str) -> dict:
    data = dict()
    with open(txt_file) as f:
        for line in f:
            line = line.lower().replace(":", "=")
            for prop in ["tem_verified", "comment"]:
                if line.startswith(prop):
                    line.replace
                    data[prop] = line.split("=")[1].split(
                        "#")[0].strip().replace(",", "-")
    return data


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

    try:
        profile_data = mat["profileData"]
    except KeyError:
        raise

    sc_idx = 1 if bool(profile_data["has_ladder"]) else 0
    best_idx = int(fold_info["bestFoldingIndex"]) - 1  # NOTE: matlab indexing!
    data = dict()

    for prop in GEL_PROPERTIES + FOLD_PROPERTIES:
        info = gel_info if prop in GEL_PROPERTIES else fold_info
        is_set = (prop in info.dtype.names)
        # NOTE: why is this in here twice?

        if prop in ["qualityMetric", "fractionMonomer", "bandWidthNormalized",
                    "migrationDistanceNormalized", "fractionPocket", "fractionSmear"] and is_set:
            prop_str = str(info[prop].item()[best_idx])
            prop_str_sc = str(info[prop].item()[sc_idx]) + "_scaffold"
            data.update({prop + "_scaffold": prop_str_sc.lower()})
        elif is_set:
            prop_str = str(info[prop]).replace(",", ".")
        else:
            prop_str = " "
        data.update({prop: prop_str.lower()})

    for prop in ["qualityMetric", "fractionMonomer", "bandWidthNormalized",
                 "migrationDistanceNormalized", "fractionPocket", "fractionSmear"]:
        if prop in fold_info.dtype.names:
            prop_float = fold_info[prop].item()[best_idx]
            prop_float_sc = fold_info[prop].item()[sc_idx]
            data.update({prop + "_scaffold": prop_float_sc})
        else:
            prop_float = 0
        data.update({prop: prop_float})
    return data


def complete_dataframe(f):

    df = pd.read_csv(f)
    df['date'] = pd.to_datetime(df['date'])
    df['year'] = pd.to_datetime(df['date']).dt.year
    df['month'] = pd.to_datetime(df['date']).dt.month

    df["co_density_scaffold"] = df["co_set_scaffold_co"] / \
        df["co_possible_scaffold_co"]

    df["co_density_staple"] = df["co_set_staple_co"] / \
        df["co_possible_staple_co"]
    df["co_density"] = (df["co_set_staple_co"] + df["co_set_scaffold_co"]) / \
        (df["co_possible_staple_co"]+df["co_possible_scaffold_co"])

    dfall = (
        df["full_scaf_co_type_1"] +
        df["full_scaf_co_type_2"] +
        df["full_scaf_co_type_3"]
    )
    for x in [1, 2, 3]:
        df[f'full_scaf_co_type_{x}_%'] = df[f"full_scaf_co_type_{x}"] / dfall

    df["tem_verified"].replace(np.nan, "no", inplace=True)
    df["tem_verified"].replace(' ', "no", inplace=True)
    df["tem_verified"].replace('tem_verified', "yes", inplace=True)
    df["tem_verified"] = df["tem_verified"].astype('category')

    df["scaffold"] = df["scaffold_type"]
    df["lattice"] = df["lattice_type"]

    df["scaffold_length"] = df["scaffold"].map(scaffold_dict_len)
    df["scaffold_gc"] = df["scaffold"].map(scaffold_dict_gc)
    df["scaffold_name"] = df["scaffold"].map(scaffold_dict_name)
    df["scaffold_circ"] = df["scaffold"].map(scaffold_dict_circ)
    df["rel_avg_loops_length"] = df["avg_loops_length"] / df["scaffold_length"]

    df["yield"] = df["fractionMonomer"] * \
        df["tem_verified"].map(dict(yes=1, no=0)).astype(int)
    df["purity"] = df["bandWidthNormalized"] * \
        df["tem_verified"].map(dict(yes=1, no=0)).astype(int)
    df["quality"] = df["yield"] * df["purity"]

    return df


@click.group()
@click.option('--version', is_flag=True, callback=print_version,
              expose_value=False, is_eager=True)
def cli_db():
    pass


# @click.argument("-i", "--db_folder", "db_folder",   default=".", help="input folder")
# @click.option("-o", "--output",  default="./__database", help="output folder")
# @click.option("-d", "--datafile",  default="fdb", help="database-file name")

@cli_db.command()
def create_database():
    """ processes database design files. matlabscript on IFS has to be run first.
        creates database.csv in specified folder
    """

    db_folder = Path(".")
    output = Path("__database")
    datafile = "fdb"
    output.mkdir(parents=True, exist_ok=True)

    date_str = str(date.today().strftime("%y-%b-%d"))
    fdb_filepath = output / f"{datafile}-{date_str}.csv"

    exclude_count = 0
    with open(fdb_filepath, mode="w+", encoding="utf-8") as fdb_file:
        for child in db_folder.iterdir():
            # try:
            exclude = (
                child.name.startswith(".")
                or "__" in child.name[-3:]
                or child.name.startswith("__")
                or not os.path.isdir(child.name)
            )
            if exclude:
                exclude_count += 1
                continue
            # TODO: error if more than one .mat file
            print(child.name)
            with get_file(logger, child, "*.json", IndexError):
                json = list(child.glob("*.json")).pop()
            with get_file(logger, child, "*.mat", IndexError):
                mat = list(child.glob("*.mat")).pop()
            with get_file(logger, child, "*.txt", IndexError):
                txt = list(child.glob("*.txt")).pop()

            try:
                mat_data = process_mat_file(mat)
            except Exception as e:
                e_ = "info.mat file " + EXC_TXT[14:].format(child.name, e)
                logger.error(e_)
                continue

            try:
                # NOTE: TEM_verified and comment not in yet in .mat 2021.02.19
                # trying to retrieve from text file
                txt_data = process_txt_file(txt)
                mat_data.update(txt_data)
            except Exception as e:
                e_ = "info.txt file " + EXC_TXT[14:].format(child.name, e)
                logger.error(e_)
                continue

            design_name = mat_data["design_name"]
            design_seq = mat_data["scaffold_type"].upper()

            # try:
            designdata = Design(
                json=json, name=design_name, seq=design_seq)
            # except Exception as e:
            #    e_ = "nanodesign    " + EXC_TXT[14:].format(child.name, e)
            #    logger.error(e_)
            #   raise e

            # try:
            compute = DesignStats(design=designdata)
            compute.compute_data()
            # except Exception as e:
            #    e_ = "designdata    " + EXC_TXT[14:].format(child.name, e)
            #    logger.error(e_)
            #    continue

            # try:
            json_data = compute.prep_data_for_export()
            data = {**mat_data, **json_data}
            # except Exception as e:
            #    e_ = "exportdata    " + EXC_TXT[14:].format(child.name, e)
            #    logger.error(e_)
            #    continue

            # try:
            export_data(data=data, fdb_file=fdb_file)
            # except Exception as e:
            #    e_ = "data export   " + EXC_TXT[14:].format(child.name, e)
            #    logger.error(e_)
            #    continue

            # except Exception as e:
            #    e_ = "critical eror   " + EXC_TXT[14:].format(child.name, e)
            #    logger.error(e_)
            #    continue
    print(exclude_count, " folders deletionped")

    df = complete_dataframe(fdb_filepath)
    df.to_csv(fdb_filepath)

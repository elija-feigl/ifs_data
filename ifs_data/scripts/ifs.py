#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import click
import logging
import scipy.io as sio
import pandas as pd
import numpy as np

from typing import IO
from datetime import date
from pathlib import Path
from shutil import copyfile

from ifs_data import _init_logging
from ifs_data.version import get_version
from ifs_data.core.utils import (get_file, GEL_PROPERTIES, FOLD_PROPERTIES,
                                 scaffold_dict_len, scaffold_dict_name, scaffold_dict_circ,
                                 scaffold_dict_gc, T_screen, Mg_screen)
from ifs_data.core.designData import Design
from ifs_data.database.compute_data import DesignStats


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
                    data[prop] = line.split("=")[1].split(
                        "#")[0].strip().replace(",", "-")
    return data


def process_txt_file_publication(txt_file: str) -> dict:
    props = ["tem_verified", "scaffold_type", "published", "lattice_type",
             "date", "user", "design_name", "project"]
    data = dict()
    with open(txt_file) as f:
        for line in f:
            line = line.replace(":", "=", 1)
            for prop in props:
                if line.lower().startswith(prop):
                    try:
                        data[prop] = line.split("=")[1].split("#")[0].strip()
                        continue
                    except Exception:
                        logger.exception(
                            f"faulty data in property ({prop}) for file: {txt_file}")

        missing = set(props) - set(data.keys())
        if missing:
            logger.error(f"missing data ({missing}) for file: {txt_file}")

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
    best_idx = int(fold_info["bestFoldingIndex"]) - 1  # matlab indexing!
    data = dict()

    for prop in GEL_PROPERTIES + FOLD_PROPERTIES:
        info = gel_info if prop in GEL_PROPERTIES else fold_info
        is_set = (prop in info.dtype.names)

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


# TODO: move to compute_data
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

    df["bestTs"] = df["bestTscrn"].map(T_screen)
    df["bestMs"] = df["bestMgscrn"].map(Mg_screen)

    df.loc[df["comment"].str.contains(
        "temperature screen only over 3", na=False), 'bestTs'] = df.loc[df["comment"].str.contains(
            "temperature screen only over 3", na=False), 'bestTs'].apply(lambda x: x[1:])

    df["bestTs_max"] = df["bestTs"].apply(lambda x: max(x))
    df["bestTs_min"] = df["bestTs"].apply(lambda x: min(x))
    return df


@click.group()
@click.option('--version', is_flag=True, callback=print_version,
              expose_value=False, is_eager=True)
def cli():
    pass


@cli.command()
@click.option("-i", "--db_folder",  type=click.Path(exists=True), default=Path("."), help="input folder")
@click.option("-o", "--output",  type=click.Path(), default=Path("./__database"), help="output folder")
@click.option("-d", "--datafile",  default="fdb", help="database-file name")
def create_database(db_folder, output, datafile):
    """ processes database design files. matlabscript on IFS has to be run first.
        creates database.csv in specified folder
    """
    logger = logging.getLogger(__name__)
    output.mkdir(parents=True, exist_ok=True)

    date_str = str(date.today().strftime("%y-%b-%d"))
    fdb_filepath = output / f"{datafile}-{date_str}.csv"

    exclude_count = 0
    with open(fdb_filepath, mode="w+", encoding="utf-8") as fdb_file:
        for child in db_folder.iterdir():
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
            logger.info(child.name)
            with get_file(logger, child, "*.json", IndexError):
                json = list(child.glob("*.json")).pop()
            with get_file(logger, child, "*.mat", IndexError):
                mat = list(child.glob("*.mat")).pop()
            with get_file(logger, child, "*.txt", IndexError):
                txt = list(child.glob("*.txt")).pop()

            mat_data = process_mat_file(mat)

            # NOTE: TEM_verified and comment not in yet in .mat 2021.02.19
            #    trying to retrieve from text file
            txt_data = process_txt_file(txt)
            mat_data.update(txt_data)

            design_name = mat_data["design_name"]
            design_seq = mat_data["scaffold_type"].upper()

            designdata = Design(
                json=json, name=design_name, seq=design_seq, circ_scaffold=True)

            compute = DesignStats(design=designdata)
            compute.compute_data()

            json_data = compute.prep_data_for_export()
            data = {**mat_data, **json_data}

            export_data(data=data, fdb_file=fdb_file)

    logger.debug(f"{exclude_count} folders skiped")

    df = complete_dataframe(fdb_filepath)
    df.to_csv(fdb_filepath)


@cli.command()
@click.argument('json', type=click.Path(exists=True))
@click.argument('sequence', type=str)
@click.option("-d", "--datafile",  default="fdb", help="database-file name")
def analyse_design(json, sequence, datafile):
    """ analyse a single design file and link it to the ifs data as a csv

        JSON is the name of the cadnano design file [.json]\n
        SEQUENCE is the name of scaffold strand\n
    """
    logger = logging.getLogger(__name__)

    outname = f"{json.name}-stat.csv"
    design = Design(json=json, name=json.name,
                    seq=sequence, circ_scaffold=True)
    logger.debut(f"Successfully read design {json}")

    compute = DesignStats(design=design)
    compute.compute_data()
    data = compute.prep_data_for_export()

    with open(outname, mode="w") as outfile:
        header = ",".join(str(k) for k in data.keys())
        outfile.write(header + "\n")
        export = ",".join(str(v) for v in data.values())
        outfile.write(export + "\n")


@cli.command()
@click.argument('json', type=click.Path(exists=True))
@click.argument('sequence', type=str)
@click.argument('database', type=click.Path(exists=True))
@click.option("-o", "--output",  type=click.Path(), default=None, help="output file name")
@click.option('--to-best', is_flag=True,
              help='compare to best designs')
@click.option('--to-same-scaffold', is_flag=True,
              help='compare to same scaffold')
@click.option('--to-same-lattice', is_flag=True,
              help='compare to same lattice')
def compare_design(json, sequence, database, output, to_best, to_same_scaffold, to_same_lattice):
    """ analyse a single design file and compare it to an existing database

        JSON is the cadnano design file [.json]\n
        DATABASE is the database file [.csv]\n
        SEQUENCE is the name of scaffold strand\n
    """
    logger = logging.getLogger(__name__)

    if output is None:
        output = f"{json.name}-comparison.csv"
    if output.suffix != ".csv":
        logger.error(f"{output} is not a database file [.csv].")

    design = Design(json=json, name=json.name,
                    seq=sequence, circ_scaffold=True)
    logger.debut(f"Successfully read design {json}")

    compute = DesignStats(design=design)
    compute.compute_data()
    data = compute.prep_data_for_export()

    # TODO: get average data from database
    #   add same scaffold if selcted
    #   add same lattice if selected
    #   add same lattice and scaffold if both
    # TODO: print "prediction" for T and Mg

    # TODO: print comparative csv
    with open(output, mode="w") as outfile:
        header = ",".join(str(k) for k in data.keys())
        outfile.write(header + "\n")
        export = ",".join(str(v) for v in data.values())
        outfile.write(export + "\n")


@cli.command()
@click.option("-i", "--db_folder",  type=click.Path(exists=True), default=Path("."), help="input folder")
@click.option("-o", "--output",  type=click.Path(), default=Path("./__publications"), help="output folder")
def create_publication_db(db_folder, output, ):
    """ parse all folders in the database and extract publication data.
        creates a new folder "__publications" wich contains a folder for each publication:
            these folders contain cadnano design file and a short info file with scaffold, user, date, etc.
    """
    logger = logging.getLogger(__name__)
    output.mkdir(parents=True, exist_ok=True)

    exclude_count = 0
    unpublished_count = 0
    for child in db_folder.iterdir():
        exclude = (
            child.name.startswith(".")
            or "__M" in child.name[-3:]
            or child.name.startswith("__")
            or not os.path.isdir(child.name)
        )
        if exclude:
            logger.debug(f"skipping folder {child.name}")
            exclude_count += 1
            continue

        logger.info(child.name)
        with get_file(logger, child, "*.json", IndexError):
            json = list(child.glob("*.json")).pop()
        with get_file(logger, child, "*.txt", IndexError):
            txt = list(child.glob("*.txt")).pop()

        txt_data = process_txt_file_publication(txt)
        publication = txt_data["published"]

        if publication.lower() not in ["", "no", "none"]:
            pub_path = output / publication
            if not (pub_path.exists() and pub_path.is_dir()):
                pub_path.mkdir()
            project = txt_data["project"]
            design_name = txt_data["design_name"]
            name = f"{project}-{design_name}"
            copyfile(json, pub_path / f"{name}.json")

            txt_info = pub_path / f"{name}.txt"
            with txt_info.open(mode='w') as f:
                for key, value in txt_data.items():
                    f.write(f"{key} = {value}\n")
        else:
            unpublished_count += 1

    logger.debug(f"{exclude_count} folders skiped")
    logger.info(f"{unpublished_count} unpublished designs ")


if __name__ == "__main__":
    _init_logging()
    create_database()

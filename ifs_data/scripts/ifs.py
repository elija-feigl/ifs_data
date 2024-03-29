#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
import os
import click
import logging
import pandas as pd
import scipy.io as sio

from typing import IO, Any, Dict
from datetime import date
from pathlib import Path
from shutil import copyfile

from ifs_data import _init_logging
from ifs_data import get_version
from ifs_data.core.utils import (get_file, GEL_PROPERTIES, FOLD_PROPERTIES,
                                 scaffold_dict_len, scaffold_dict_name, scaffold_dict_circ,
                                 scaffold_dict_gc, T_screen, Mg_screen, _str)
from ifs_data.core.designData import Design
from ifs_data.database.compute_data import DesignStats

TOP_X = 0.2
logger = logging.getLogger(__name__)


def print_version(ctx, _, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(get_version())
    ctx.exit()


def export_data(data: dict, fdb_file: IO) -> None:
    if not fdb_file.tell():  # true if empty
        header = ",".join(str(k) for k in data.keys())
        fdb_file.write(header + "\n")
    export = ",".join(_str(v) for v in data.values())
    fdb_file.write(export + "\n")


def _process_txt_file(txt_file: Path) -> dict:
    # NOTE: remove once matlab script includes these values
    data = dict()
    with txt_file.open() as f:
        for line in f:
            line = line.replace(":", "=", 1)
            for prop in ["tem_verified", "comment"]:
                if line.lower().startswith(prop):
                    data[prop] = line.split("=")[1].split(
                        "#")[0].strip().replace(",", "-")

    missing = set(["tem_verified", "comment"]) - set(data.keys())
    if missing:
        logger.error(f"missing data ({missing}) for file: {txt_file}")
    return data


def process_txt_file_publication(txt_file: Path) -> dict:
    props = ["published", "date", "user", "design_name", "project"]
    data = dict()
    with txt_file.open() as f:
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


def process_mat_file(mat_file: Path, txt_file: Path) -> dict:

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
        # legacy version
        has_ladder = mat["profileData"]["has_ladder"]
        sc_idx = 1 if bool(has_ladder) else 0
    except KeyError:
        # new version
        try:
            sc_idx = gel_info["species"].item()["scaffold"].item()[
                "indices"].item()[0] - 1
        except KeyError:
            raise

    best_idx = int(fold_info["bestFoldingIndex"]) - 1  # matlab indexing!
    data: Dict[str, Any] = dict()

    for prop in GEL_PROPERTIES + FOLD_PROPERTIES:
        info = gel_info if prop in GEL_PROPERTIES else fold_info
        is_set = (prop in info.dtype.names)

        if prop in ["qualityMetric", "fractionMonomer", "bandWidthNormalized",
                    "migrationDistanceNormalized", "fractionPocket", "fractionSmear"] and is_set:
            prop_str = str(info[prop].item()[best_idx])
            prop_str_sc = str(info[prop].item()[sc_idx]) + "_scaffold"
            data.update({prop + "_scaffold": prop_str_sc})
        elif is_set:
            prop_str = str(info[prop]).replace(",", ".")
        else:
            prop_str = " "
        if prop == "name":
            prop_str = prop_str.upper()
        if prop == "lattice_type":
            # correct some erroneous input in database
            prop_str = prop_str.split("_")[0].lower()
        data.update({prop: prop_str})

    for prop in ["qualityMetric", "fractionMonomer", "bandWidthNormalized",
                 "migrationDistanceNormalized", "fractionPocket", "fractionSmear"]:
        if prop in fold_info.dtype.names:
            prop_float = fold_info[prop].item()[best_idx]
            prop_float_sc = fold_info[prop].item()[sc_idx]
            data.update({prop + "_scaffold": prop_float_sc})
        else:
            prop_float = 0
        data.update({prop: prop_float})

    # NOTE: TEM_verified and comment not in yet in .mat 2021.02.19
    #    trying to retrieve from text file
    txt_data = _process_txt_file(txt_file)
    data.update(txt_data)

    data["lattice"] = data["lattice_type"]
    data["scaffold"] = data["scaffold_type"]
    try:
        data["scaffold_name"] = scaffold_dict_name[data["scaffold"]]
        data["scaffold_length"] = scaffold_dict_len[data["scaffold_name"]]
        data["scaffold_gc"] = scaffold_dict_gc[data["scaffold_name"]]
        data["scaffold_circ"] = scaffold_dict_circ[data["scaffold_name"]]
    except KeyError:
        d_sc = data["scaffold"]
        logger.error(f"Scaffold type unknown {d_sc}")
        raise

    tem_verified = True if data["tem_verified"] == "yes" else False
    if not tem_verified:
        logger_str = data["tem_verified"]
        logger.debug(f"Not tem_verified: {logger_str}")

    data["yield"] = float(data["fractionMonomer"]) * tem_verified
    data["purity"] = float(data["bandWidthNormalized"]) * tem_verified
    data["quality"] = data["yield"] * data["purity"]

    data["bestTs"] = T_screen[data["bestTscrn"]]
    data["bestMs"] = Mg_screen[data["bestMgscrn"]]

    if "temperature screen only over 3" in data["comment"]:
        # TODO: check copy
        data["bestTs"] = data["bestTs"][1:]
        logger.debug("adapting bestTs for shortened T-screen")

    data["bestTs_max"] = max(data["bestTs"])
    data["bestTs_min"] = min(data["bestTs"])
    data.pop("bestTs", None)
    return data


@ click.group()
@ click.option('--version', is_flag=True, callback=print_version,
               expose_value=False, is_eager=True)
def cli():
    pass


@ cli.command()
@ click.option("-i", "--db_folder", default=Path("."), help="input folder", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@ click.option("-o", "--output",  default=Path("./__database"), help="output folder", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@ click.option("-d", "--datafile",  default="fdb", help="database-file name")
def create_database(db_folder, output, datafile):
    """ combine stats from from database design files & IFS_matlab data.
            creates .csv in specified folder
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

            try:
                json = get_file(logger, child, "*json")
                mat = get_file(logger, child, "*mat")
                txt = get_file(logger, child, "*txt")
            except IOError:
                logger.error(f"Folder {child.name} incomplete.")
                continue

            try:
                mat_data = process_mat_file(mat, txt)
            except KeyError:
                logger.error(
                    f"Folder {child.name} has incompatible .mat file.")
                continue
            logger.info(f"Folder {child.name}")

            design_name = mat_data["design_name"]
            design_seq = mat_data["scaffold_type"].upper()
            design_data = Design(
                json=json, name=design_name, seq=design_seq, circ_scaffold=True)

            compute = DesignStats(design=design_data)
            compute.compute_data()

            json_data = compute.prep_data_for_export()
            data = {**mat_data, **json_data}
            export_data(data=data, fdb_file=fdb_file)
    logger.debug(f"{exclude_count} folders skipped.")


@ cli.command()
@ click.argument('json', type=click.Path(exists=True, resolve_path=True, path_type=Path))
@ click.argument('sequence', type=str)
@ click.option("-d", "--datafile",  default="fdb", help="database-file name")
def analyze_design(json, sequence, datafile):
    """ analyze a single design file with output as a csv

        JSON is the name of the cadnano design file [.json]\n
        SEQUENCE is the name of scaffold strand\n
    """
    logger = logging.getLogger(__name__)
    out_name = f"{json.stem}-stat.csv"
    design = Design(json=json, name=json.stem,
                    seq=sequence, circ_scaffold=True)
    logger.debug(f"Successfully read design {json}")

    compute = DesignStats(design=design)
    compute.compute_data()
    data = compute.prep_data_for_export()

    with open(out_name, mode="w") as outfile:
        header = ",".join(str(k) for k in data.keys())
        outfile.write(header + "\n")
        export = ",".join(_str(v) for v in data.values())
        outfile.write(export + "\n")


@ cli.command()
@ click.argument('json',  type=click.Path(exists=True, resolve_path=True, path_type=Path))
@ click.argument('sequence', type=str)
@ click.argument('database', type=click.Path(exists=True, resolve_path=True, path_type=Path))
@ click.option("-o", "--output",  default=None, help="output file name", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@ click.option('--to-best', is_flag=True,
               help='compare to best designs')
@ click.option('--to-scaffold', is_flag=True,
               help='compare to same scaffold')
@ click.option('--to-lattice', is_flag=True,
               help='compare to same lattice')
def compare_design(json, sequence, database, output, to_best, to_scaffold, to_lattice):
    """ analyze a single design file and compare it to an existing database

        JSON is the cadnano design file [.json]\n
        SEQUENCE is the name of scaffold strand\n
        DATABASE is the database file [.csv]\n
    """
    logger = logging.getLogger(__name__)

    if output is None:
        output = Path(f"{json.stem}-comparison.csv")
    else:
        output = Path(output)

    design = Design(json=json, name=json.stem,
                    seq=sequence, circ_scaffold=True)
    logger.debug(f"Successfully read design {json}")
    compute = DesignStats(design=design)
    compute.compute_data()
    data = compute.prep_data_for_export()

    df_data = pd.DataFrame(data, index=[f"{json.stem}"])

    if database.suffix != ".csv":
        logger.error(f"{output} is not a database file [.csv].")
    db = pd.read_csv(database)
    df_data.loc["mean"] = db.mean()

    if to_scaffold and to_lattice:
        db_select = db[(db["scaffold"] == sequence)
                       & (db["lattice"] == design.lattice)]
        select = "_scaffold-lattice"
    elif to_scaffold:
        db_select = db[db["scaffold"] == sequence]
        select = "_scaffold"
    elif to_lattice:
        db_select = db[db["lattice"] == design.lattice]
        select = "_lattice"

    if to_best:
        db = db[db["quality"] > (1.-TOP_X) * db["quality"].max()]
        db_select = db_select[db_select["quality"]
                              > (1.-TOP_X) * db_select["quality"].max()]
        top = f"_top{int(TOP_X*100)}%"
    else:
        top = ""

    df_data.loc[f"mean{top}"] = db.mean()
    df_data.loc[f"comp_{top}{select}"] = db_select.mean()

    # TODO: print "prediction" for T and Mg
    # TODO: print reduced selection to std?
    df_data.to_csv(output, float_format='%.2f')


@ cli.command()
@ click.option("-i", "--db_folder",  default=Path("."), help="input folder", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@ click.option("-o", "--output", default=Path("./__publications"), help="output folder", type=click.Path(exists=True, resolve_path=True, path_type=Path))
def create_publication_db(db_folder, output):
    """ parse all folders in the database and extract publication data.
        creates a new folder "__publications" which contains a folder for each publication:
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

        try:
            json = get_file(logger, child, "*json")
            txt = get_file(logger, child, "*txt")
        except IOError:
            logger.error(f"Folder {child.name} incomplete.")
            continue
        logger.info(f"Folder {child.name}")

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

    logger.debug(f"{exclude_count} folders skipped")
    logger.info(f"{unpublished_count} unpublished designs ")


if __name__ == "__main__":
    _init_logging()
    create_database(Path("."), Path("./__publications"), "fdb")

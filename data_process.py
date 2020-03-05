import os
import glob
import designData
from pathlib import Path
import logging
from datetime import date
import csv
import scipy.io as sio


FOLDER_EXCEPTION = [
    "Foldingscreen_analysis-master", "AAA_TEMPLATE",
    "ZZZ_output", "ZZZfolder_of_shame_aka_missing_data",
    "ZZZnon_standard_folding_screens", "Icon", ".DS_Store",
    "JF_MP_ExcessFiles"
]


def export_data(data: dict, name: str, outputname: str, filename: str, mat_data: dict) -> None:

    export = designData.prep_data_for_export(data)
    header = ", ".join([str(i) for i in mat_data.keys()]) + \
        ", " + ", ".join([str(i) for i in export.keys()])
    export_str = ", ".join([str(i) for i in mat_data.values()]) + \
        ", " + ", ".join([str(i) for i in export.values()])

    try:
        with open(filename) as out:
            has_header = csv.Sniffer().has_header(out.read(1024))
            if not has_header:
                out.write(header + "\n")

        with open(filename, mode="a") as out:
            out.write(export_str + "\n")

    except FileNotFoundError:
        with open(filename, mode="w+") as out:
            out.write(header + "\n")
        with open(filename, mode="a") as out:
            out.write(export_str + "\n")

    return


def read_mat(mat_file: dict):
    data = dict()
    mat = sio.loadmat(mat_file[0], squeeze_me=True)

    types = ['user', 'project', 'design_name', 'date',
             'scaffold_type', 'lattice_type', 'scaffold_concentration',
             'staple_concentration', 'gelsize', 'agarose_concentration', 'staining',
             'mg_concentration', 'voltage', 'running_time', 'cooling']

    """
    types = list(mat['gelInfo'].dtype.names)
    exception = ['lanes', 'log_file', 'lanes_unparsed', 'comment', 'filename']
    """
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
    return data


def main():

    logging.basicConfig()
    handle = "folding-DB"
    logger = logging.getLogger(handle)

    outputname = "database"
    try:
        os.mkdir("./" + outputname)
    except FileExistsError:
        pass

    filename = "./" + outputname + "/" + "foldingdatabase-" + \
        str(date.today().strftime("%y-%b-%d")) + ".csv"

    try:
        os.remove(filename)

    except FileNotFoundError:
        pass

    parent_folder = "../Foldingscreens_202002/"
    folders = Path("../Foldingscreens_202002")

    for folder in folders.iterdir():
        if folder.name in FOLDER_EXCEPTION + [outputname]:
            continue
        txt_file = glob.glob(parent_folder + folder.name + "/*.txt")
        mat_file = glob.glob(parent_folder + folder.name + "/*.mat")

        if txt_file:
            with open(txt_file[0], 'r', encoding="utf8") as gel_info:
                jsons = glob.glob(parent_folder + folder.name + "/*.json")

                if jsons:
                    json = jsons[0]

                    try:
                        for line in gel_info:
                            if line.startswith("Design_name"):
                                try:

                                    mat_data = read_mat(mat_file)

                                    name = mat_data["design_name"]

                                    designdata = designData.DesignData(
                                        json=json, name=name)

                                except Exception as e:
                                    e_ = "nanodesign:  {} | Error: {}".format(
                                        name, e)
                                    logger.error(e_)
                                try:
                                    data = designdata.compute_data()
                                    export_data(data=data, name=name,
                                                outputname=outputname, filename=filename, mat_data=mat_data)
                                except Exception as e:
                                    e_ = "stats:       {} | Error: {}".format(
                                        name, e)
                                    logger.error(e_)

                    except UnicodeDecodeError as e:
                        e_ = "stats:       {} | Error: {}".format(
                            folder.name, e)
                        logger.error(e_)
                else:
                    e_ = "json missing:   | Folder: {}".format(folder.name)
                    logger.warning(e_)
        else:
            e_ = "???  :          | Folder: {}".format(folder.name)
            logger.warning(e_)

    return


if __name__ == "__main__":
    main()

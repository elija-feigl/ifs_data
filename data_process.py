import os
import glob
import designData
from pathlib import Path
import logging
from datetime import date


FOLDER_EXCEPTION = [
    "Foldingscreen_analysis-master", "AAA_TEMPLATE",
    "ZZZ_output", "ZZZfolder_of_shame_aka_missing_data",
    "ZZZnon_standard_folding_screens", "Icon", ".DS_Store",
    "JF_MP_ExcessFiles"
]


def export_data(data: dict, name: str) -> None:

    export = designData.prep_data_for_export(data)
    header = ", ".join([str(i) for i in export.keys()])
    export_str = ", ".join([str(i) for i in export.values()])

    filename = "foldingdatabase-" + \
        str(date.today.strftime("%y-%b-%d")) + ".csv"

    try:
        with open("./database/" + filename, mode="r+") as out:
            if header != out.readline(0):
                out.write(header + "\n")

        with open("./database/" + filename, mode="a") as out:
            out.write(export_str + "\n")

    except FileNotFoundError:
        with open("./database/" + filename, mode="w+") as out:
            out.write(header + "\n")

        with open("./database/" + filename, mode="a") as out:
            out.write(export_str + "\n")

    return


def main():

    logging.basicConfig()
    handle = "folding-DB"
    logger = logging.getLogger(handle)

    outputname = "database"
    try:
        os.mkdir("./" + outputname)
    except FileExistsError:
        pass

    folders = Path("./")
    for folder in folders.iterdir():
        if folder.name in FOLDER_EXCEPTION + [outputname]:
            continue
        folder_content = glob.glob(folder.name + "/*.txt")
        if folder_content:
            with open(folder_content[0], 'r', encoding="utf8") as gel_info:

                jsons = glob.glob(folder.name + "/*.json")
                if jsons:
                    json = jsons[0]

                    for line in gel_info:
                        if line.startswith("Design_name ="):
                            try:
                                name = line[13:-1].strip()
                                designdata = designData.DesignData(
                                    json=json, name=name
                                )
                            except Exception as e:
                                e_ = "nanodesign:  {} | Error: {}".format(
                                    name, e)
                                logger.error(e_)
                            try:
                                data = designdata.compute_data()
                                export_data(data=data, name=name)
                            except Exception as e:
                                e_ = "stats:       {} | Error: {}".format(
                                    name, e)
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

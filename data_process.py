import os
import glob
import designData
from pathlib import Path


def export_data(data: dict, name: str) -> None:

    export = designData.prep_data_for_export(data)
    header = ", ".join([str(i) for i in export.keys()])
    export_str = ", ".join([str(i) for i in export.values()])

    try:
        with open("./database/designdata.csv", mode="r+") as out:
            if header == out.readline(0):
                pass
    except FileNotFoundError:
        with open("./database/designdata.csv", mode="w+") as out:
            out.write(header + "\n")

    with open("./database/designdata.csv", mode="a") as out:
        out.write(export_str + "\n")
        # out.write("\nEND")
    return


def main():

    try:
        os.mkdir("./database")
    except FileExistsError:
        pass

    folders = Path("./")
    for folder in folders.iterdir():
        #import ipdb; ipdb.set_trace()
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
                                data = designdata.compute_data()
                                export_data(data=data, name=name)
                            except BaseException:
                                print("issue with ", json, "of design ", name)
               else:
                   print("it was probably floris who forgot to upload a designfile... ")
        else:
            print("issue in folder:", folder.name)

    return


if __name__ == "__main__":
    main()

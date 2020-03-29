"""
[this code is to get a single stucture data from the designData.py file.
It will store structer data as a csv file in a folder called "out" in the repository.
]
"""

import designData
import os


def export_data(data: dict, name: str) -> None:

    export = designData.prep_data_for_export(data)
    header = ", ".join([str(i) for i in export.keys()])
    export_str = ", ".join([str(i) for i in export.values()])

    try:
        os.mkdir("out")
    except FileExistsError:
        pass
    with open("./out/" + name + "-designdata.csv", mode="w+") as out:

        out.write(header + "\n")
        out.write(export_str + "\n")
        # out.write("\nEND")
    return


def main():

    print("master, I am awaiting the name of your design")
    json = "TTcorr"  # enter file name here. it has to be in the same location of this file.'
    name = json
    print("Thank you Sir")

    designdata = designData.DesignData(name=name, json=json + ".json")
    data = designdata.compute_data()
    export_data(data=data, name=name)
    return


if __name__ == "__main__":
    main()

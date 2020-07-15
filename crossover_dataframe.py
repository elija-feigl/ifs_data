from designData import DesignData
import pandas as pd
import argparse
from pathlib import Path


def crossover_dataframe(designdata):
    data = {
        'type': [], 'strand_type': [], 'helices': [],
        'orientation': [], 'positions': [], 'Cadnano_type(integers)': []
    }
    for co in designdata.all_crossovers:
        data['type'].append(co.typ)
        data['strand_type'].append(co.strand_typ)
        data['helices'].append(tuple(co.h))
        data['orientation'].append(co.orientation)
        data['positions'].append(co.coordinate)
        data['Cadnano_type(integers)'].append(co.scaff_full_type)

    df_crossover = pd.DataFrame(data, columns=list(data.keys()))
    df_crossover.to_csv(r'crossover_df.csv')


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input",
                        help="input file",
                        type=str,
                        default="TTcorr.json",
                        )
    args = parser.parse_args()
    json = Path(args.input)
    designdata = DesignData(json=json, name=json.name, seq='8064')
    crossover_dataframe(designdata)


if __name__ == "__main__":
    main()

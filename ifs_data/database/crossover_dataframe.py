import pandas as pd
import argparse
from pathlib import Path

from ..core.designData import Design
from ..data.crossover import Crossover


def crossover_dataframe(designdata: Design):
    data = {
        'type': [], 'strand_type': [], 'helices': [],
        'orientation': [], 'positions': [], 'Cadnano_type(integers)': []
    }
    for co in designdata.crossovers:
        bases = co.get_bases()
        data['type'].append(co.typ)
        data['strand_type'].append(co.is_scaffold)
        data['helices'].append({base.h for base in bases})
        data['orientation'].append(co.orientation)
        data['positions'].append({(base.h, base.p) for base in bases})
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
    designdata = Design(json=json, name=json.name,
                        seq='8064', circ_scaffold=False)
    crossover_dataframe(designdata)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
import pandas as pd
import argparse
from pathlib import Path

from ..core.designData import Design


def crossover_dataframe(design_data: Design):
    data: dict = {
        'type': [], 'strand_type': [], 'helices': [],
        'orientation': [], 'positions': [], 'Cadnano_type(integers)': []
    }
    for co in design_data.crossovers:
        bases = co.bases
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
    design_data = Design(json=json, name=json.name,
                        seq='8064', circ_scaffold=False)
    crossover_dataframe(design_data)


if __name__ == "__main__":
    main()

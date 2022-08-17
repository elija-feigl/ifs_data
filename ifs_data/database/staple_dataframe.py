#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
import pandas as pd
import argparse
from pathlib import Path

from ..core.designData import DesignData


def staple_dataframe(design_data):
    data = {
        'ID': [], 'length': [], 'n_helices': [], 'helices': [], "n_domains": [], "domain_lengths": [],
        "domain_melt_t": [], "n_long_domains": [], 'position_5prime': []
    }
    for staple in design_data.staples:
        data['ID'].append(staple.id)
        data['length'].append(len(staple.tour))
        data['n_helices'].append(
            len(design_data.staple_helix_dict[staple]))
        data['helices'].append(
            list(design_data.staple_helix_dict[staple]))
        data["n_domains"].append(len(design_data.domain_data[staple]))
        data['domain_lengths'].append(design_data.domain_lengths_data[staple])
        data["domain_melt_t"].append(design_data.staple_domains_melt_t[staple])
        data['n_long_domains'].append(design_data.long_domains[staple])
        first = (staple.tour[0].h, staple.tour[0].p)
        data['position_5prime'].append(first)

    df_staple = pd.DataFrame(data, columns=list(data.keys()))
    df_staple.to_csv(r'staple.csv')


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
    design_data = DesignData(json=json, name=json.name, seq='8064')
    staple_dataframe(design_data)


if __name__ == "__main__":
    main()

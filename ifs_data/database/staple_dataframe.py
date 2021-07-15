#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2021  Elija Feigl
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see https://www.gnu.org/licenses/gpl-3.0.html.

import pandas as pd
import argparse
from pathlib import Path

from ..core.designData import DesignData


def staple_dataframe(designdata):
    data = {
        'ID': [], 'length': [], 'n_helices': [], 'helices': [], "n_domains": [], "domain_lengths": [],
        "domain_melt_t": [], "n_long_domains": [], 'position_5prime': []
    }
    for staple in designdata.staples:
        data['ID'].append(staple.id)
        data['length'].append(len(staple.tour))
        data['n_helices'].append(
            len(designdata.staple_helix_dict[staple]))
        data['helices'].append(
            list(designdata.staple_helix_dict[staple]))
        data["n_domains"].append(len(designdata.domain_data[staple]))
        data['domain_lengths'].append(designdata.domain_lengths_data[staple])
        data["domain_melt_t"].append(designdata.staple_domains_melt_t[staple])
        data['n_long_domains'].append(designdata.long_domains[staple])
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
    designdata = DesignData(json=json, name=json.name, seq='8064')
    staple_dataframe(designdata)


if __name__ == "__main__":
    main()

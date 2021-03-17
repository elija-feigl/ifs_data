#!/usr/bin/env python
# -*- coding: utf-8 -*-

import click
import logging
from pathlib import Path

import scipy.io as sio

from ..core.compute_data import Compute
from ..core.designData import DesignData


logger = logging.getLogger(__name__)


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(get_version())
    ctx.exit()


@click.group()
@click.option('--version', is_flag=True, callback=print_version,
              expose_value=False, is_eager=True)
def cli():
    pass


@cli.command()
@click.argument('json', type=click.Path(exists=True))
@click.argument('sequence', type=str)
def analyse_design():
    """ analyse a single design file and link it to the ifs data as a csv

        JSON is the name of the cadnano design file [.json]\n
        SEQUENCE is the name of scaffold strand\n
    """

    outname = f"{json.name}-stat.csv"
    designdata = DesignData(json=json, name=json.name, seq=sequence)
    compute = Compute()
    compute.compute_data(designdata)
    data = compute.prep_data_for_export(designdata)

    with open(outname, mode="w") as outfile:
        header = ",".join(str(k) for k in data.keys())
        outfile.write(header + "\n")
        export = ",".join(str(v) for v in data.values())
        outfile.write(export + "\n")

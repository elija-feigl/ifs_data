#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
from typing import Optional, Set

import attr
import numpy as np
from nanodesign.data.base import DnaBase as Base
from ..core.utils import _hps


@ attr.s(slots=True, frozen=True)
class Connection(object):
    """ interhelical connection of two bases """
    base1: Base = attr.ib()
    base2: Base = attr.ib()

    def __str__(self):
        return f"Connection: {[_hps(base) for base in self.bases]}"

    @ property
    def bases(self) -> Set[Base]:
        return {self.base1, self.base2}


@ attr.s(slots=True)
class Crossover(object):
    """ Dna Origami Crossover object
        composed of up to two connections, hence up to 4 bases
    """
    connection1: Connection = attr.ib()
    connection2: Optional[Connection] = attr.ib(default=None)
    is_end: bool = attr.ib(default=False)

    is_scaffold: bool = attr.ib(default=False)
    typ: str = attr.ib(default="")

    scaff_full_type: int = attr.ib(default=-1)
    orientation: str = attr.ib(default="")

    def __attrs_post_init__(self):
        self.is_scaffold = self.connection1.base1.is_scaf
        self.typ = self._get_typ()

    def __str__(self):
        strand_type = "scafffold" if self.is_scaffold else "staple"
        return f"Crossover: {strand_type}-{self.typ}-{self.orientation}: {[_hps(base) for base in self.bases]}"

    def _get_typ(self):
        if self.is_end:
            return "end"
        elif self.connection2 is None:
            return "half"
        else:
            return "full"

    @ property
    def bases(self) -> Set[Base]:
        bases = self.connection1.bases
        if self.connection2 is not None:
            bases |= self.connection1.bases
        return bases

    def set_orientation(self, helices):
        con = self.connection1
        is_vertical = (
            helices[con.base1.h].lattice_row
            == helices[con.base2.h].lattice_row
        )
        self.orientation = "vertical" if is_vertical else "horizontal"

    def set_scaffold_subtyp(self, helices, lattice):
        if self.is_scaffold and self.typ == "full":
            self.scaff_full_type = self._get_scaffold_subtyp(helices, lattice)

    def _get_scaffold_subtyp(self, helices, lattice):
        """[assign type to full crossover: position of the full scaffold crossover depending on
        the position suggested by cadnano]
        """

        helix1 = helices[self.connection1.base1.h]
        helix2 = helices[self.connection1.base2.h]
        co_p = (self.connection1.base1.p + self.connection2.base1.p) // 2

        # NOTE: nanodesign.crossover is equivalent to this.Connection
        possible_ps = [
            p for (helix, p, _) in helix1.possible_staple_crossovers if helix is helix2]

        # find closest possible Connection
        sub = np.Infinity
        for possible in possible_ps:
            sub_new = co_p - possible
            if abs(sub_new) <= abs(sub):
                sub = sub_new

        if sub is np.Infinity:
            return 0

        # calculate lattice depended type
        if lattice == "square":
            mod = sub % 32
            if 0 <= mod <= 11:
                return 1
            elif 21 <= mod < 32:
                return 3
            else:
                return 2
        else:
            mod = sub % 21
            if 0 <= mod <= 11:
                return 1
            else:
                return 3

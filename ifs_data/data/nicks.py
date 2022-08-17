#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
from dataclasses import dataclass
from typing import Set

from nanodesign.data.base import DnaBase as Base


@dataclass
class Nick(object):
    base1: Base
    base2: Base

    def __post_init__(self):
        if self.base1.h != self.base2.h:
            raise ValueError

    @property
    def bases(self) -> Set[Base]:
        return {self.base1, self.base2}

    @property
    def p(self) -> Set[int]:
        """positions of bases"""
        return {self.base1.p, self.base2.p}

    @property
    def h(self) -> Set[int]:
        """helix of the nick"""
        return self.base1.h

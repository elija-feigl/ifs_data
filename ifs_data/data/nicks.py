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

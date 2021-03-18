import attr
from typing import Tuple, Set

from nanodesign.data.base import DnaBase as Base


@attr.s
class Nick(object):
    base1: Base = attr.ib()
    base2: Base = attr.ib()

    def __attrs_post_init__(self):
        self.bases: Tuple[Base, Base] = (self.base1, self.base2)
        self.p: Set[int] = set([(self.base1.p, self.base2.p)])
        self.h: Set[int] = set([(self.base1.h, self.base2.h)])

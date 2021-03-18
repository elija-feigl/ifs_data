import attr
from typing import Tuple, Set

from nanodesign.data.base import DnaBase as Base

# NOTE: do we actually need this? waht additional benefit does it gives


@attr.s
class Nick(object):
    bases: Tuple[Base, Base] = attr.ib()
    ps: Set[int] = attr.ib()
    hs: Set[int] = attr.ib()

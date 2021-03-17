from typing import List

__version__ = "0.4"
__authors__ = ["Elija Feigl", "Kourosh Zargari"]
__copyright__ = "Copyright 2021, Dietzlab (TUM)"
__license__ = "GPL-3.0"
__maintainer__ = "Elija Feigl"
__email__ = "elija.feigl@tum.de"
__status__ = "Development"


def get_version() -> str:
    return __version__


def get_authors() -> List[str]:
    return __authors__

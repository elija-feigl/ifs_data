
import nanodesign


class Crossover(object):

    def __init__(self, typ, is_scaf, is_vertical, coordinate, bases):
        """
            Initialize a Crossover object.

            Arguments:
                typ (str): type of crossover (full, half, endloop).
                is_scaf (bool): is a crossover for a (scaffold = True) or a (staple = False).
                is_vertical (bool): indicates the orientation of the crossover
                                    in the DNA origami structure.[vertical = True] & [horizontal = False]
                coordinates (list): the position and helix (base.p, base.h)
                                    number of each base that is connected in the crossover.
                bases (tuple): the bases that are connected togather via the crossover.

        """

        self.typ = typ
        self.is_scaf = is_scaf
        self.is_vertical = is_vertical
        self.coordinate = coordinate
        self.bases = bases

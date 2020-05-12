
import nanodesign


class Crossover(object):

    def __init__(self, typ, strand_typ, p, h, is_vertical, coordinate, bases, scaff_full_type):
        """
            Initialize a Crossover object.

            Arguments:
                typ (str): type of crossover (full, half, endloop).
                strand_typ (bool): is a crossover for a (scaffold = True) or a (staple = False).
                is_vertical (bool): indicates the orientation of the crossover
                                    in the DNA origami structure.[vertical = True] & [horizontal = False]
                coordinates (list): the position and helix (base.p, base.h)
                                    number of each base that is connected in the crossover.
                bases (tuple): the bases that are connected togather via the crossover.

        """

        self.typ = typ
        self.strand_typ = strand_typ
        self.p = p
        self.h = h
        self.is_vertical = is_vertical
        self.coordinate = coordinate
        self.bases = bases
        self.scaff_full_type = None

    def set_full_scaf_type(self, typ):
        self.scaff_full_type = typ

    def get_typ(self, typ):
        self.typ = typ

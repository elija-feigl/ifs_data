from typing import Optional


#@attr.s(slots=True)
class Crossover(object):
    __slots__ = ["typ", "scaff_full_type", "coordinate",
                 "h", "p", "orientation", "bases", "strand_typ"]

    def __init__(self, typ, con_tuple, helices):
        """
            Initialize a Crossover object.

            Arguments:
                typ (str): type of crossover (full, half, endloop).
                con_tuple tuple(tuple(base,base), optinal(tuple(base,base)))

                strand_typ (bool): is a crossover for a (scaffold = True) or a (staple = False).
                is_vertical (bool): indicates the orientation of the crossover
                                    in the DNA origami structure.[vertical = True] & [horizontal = False]
                coordinates (list): the position and helix (base.h, base.p)
                                    number of each base that is connected in the crossover.
                bases (tuple): the bases that are connected togather via the crossover.

        """
        self.typ: str = typ
        self.scaff_full_type: Optional[int] = None
        self.create_crossover(con_tuple, helices)

    def create_crossover(self, con_tuple, helices):
        """[this function creats crossover objects with attributes which are available in crossover class]

        Arguments:
            typ {[str]} -- [crossover type: full, half, endloop]
            con_tuple  tuple(tuple(base,base), optinal(tuple(base,base))) -- [a resemblance of a crossover made of a
                            tuple consisting two basis which indicates a crossover]

        Returns:
            crossover Object -- [description]
        """
        tmp_bases = list([None, None])
        for i, con in enumerate(con_tuple):
            if con is None:
                continue
            tmp_bases[i] = con if con[0].h < con[1].h else con[::-1]
        if con_tuple[1] is None:
            self.bases = tuple(tmp_bases)
        else:
            self.bases = tuple(
                tmp_bases) if tmp_bases[0][0].p < tmp_bases[1][0].p else tuple(tmp_bases[::-1])

        first_base = self.bases[0][0]
        self.strand_typ = 'scaffold' if first_base.is_scaf else 'staple'

        self.h = set()
        self.p = set()
        self.coordinate = list()

        for con in self.bases:
            if con is None:
                continue
            for base in con:
                self.coordinate.append((base.h, base.p))
                self.p.add(base.p)
                self.h.add(base.h)

        is_vertical = (
            helices[self.bases[0][0].h].lattice_row
            == helices[self.bases[0][1].h].lattice_row
        )
        self.orientation = "vertical" if is_vertical else "horizontal"

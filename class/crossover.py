import nanodesign
import designData


class crossover(object):

    def __init__(self):

    def _get_all_co_tuple(self) -> list:
        all_co_tuples = set()
        all_co_tuples_list = []
        for strand in self.all_strands:
            if strand.is_scaffold:

                new_strand = self._close_strand(strand)
                self.all_strands[strand.id] = new_strand

        for strand in self.all_strands:
            for base in strand.tour:
                if self.dna_structure._check_base_crossover(base):
                    co_tuple = tuple()
                    if base.up.h != base.h:
                        co_tuple = (base, base.up)
                        all_co_tuples.add(tuple(set(co_tuple)))
                    elif base.down.h != base.h:
                        co_tuple = (base.down, base)
                        all_co_tuples.add(tuple(set(co_tuple)))
        for co in all_co_tuples:
            all_co_tuples_list.append(co)

        return all_co_tuples_list

    def _get_horozantal_vertical_co(self):
        all_co_tuples_h = set()
        all_co_tuples_v = set()
        helix_row = []
        map_id_helices = self.dna_structure.structure_helices_map
        for co_tuple in self.all_co_tuples_list:
            for base in co_tuple:
                helix_row.append(map_id_helices[base.h].lattice_row)
            # helix2_row = map_id_helices[co_tuple[1].helix].lattice_row

            if helix_row[0] == helix_row[1]:
                all_co_tuples_h.add(co_tuple)
            else:
                all_co_tuples_v.add(co_tuple)
            helix_row = []

        return all_co_tuples_h, all_co_tuples_v

    def _get_full_co_list(self) -> list:
        """[gets the full crossovers]

        Returns:

            list -- [full_co_list_packed: get full co as a pack of two crossover (representation: [Co[B,B],Co[B,B], all in lists)
                     full_co_list_seperate: also seperately as individual crossovers (every Co in frozenset of two bases)]
        """

        full_co_list = []
        for co in self.all_co_tuples_list:
            co_neighbours = {"co_tuples_plus": tuple(),
                             "co_tuples_minus": tuple()
                             }
            co_tuple_plus = set()
            co_tuple_minus = set()

            for base in co:
                base_plus = self.get_base_from_hps(
                    base.h, base.p + 1, base.is_scaf)
                base_minus = self.get_base_from_hps(
                    base.h, base.p - 1, base.is_scaf, dir=-1)
                co_tuple_plus.add(base_plus)
                co_tuple_minus.add(base_minus)
            co_neighbours["co_tuples_plus"] = tuple(co_tuple_plus)
            co_neighbours["co_tuples_minus"] = tuple(co_tuple_minus)

            fullco = set()
            for typ in ["co_tuples_plus", "co_tuples_minus"]:
                if co_neighbours[typ] in self.all_co_tuples_list:
                    fullco.add(co)
                    fullco.add(co_neighbours[typ])
                    full_co_list.append(frozenset(fullco))

                co_neighbours[typ] = []

        # putting all full_co in a list configuration as [(Co(B,B),Co(B,B))] two parallel Co in a tuple and two bases also in a tuple

        full_co_set = set(full_co_list)
        full_co_list_seperate = []
        full_co_list_packed = []

        for full in full_co_set:
            full_co_list_packed.append(tuple(full))
            for co in full:
                full_co_list_seperate.append(co)

        return full_co_list_seperate, full_co_list_packed

    def _get_endloop_co_list(self) -> list:
        end_co_list = []
        for co in self.all_co_tuples_list:
            for base in co:
                base_plus = self.get_base_from_hps(
                    base.h, base.p + 1, base.is_scaf)
                base_minus = self.get_base_from_hps(
                    base.h, base.p - 1, base.is_scaf, dir=-1)

                if (base_plus is None) or (base_minus is None):
                    if co not in end_co_list:
                        end_co_list.append(co)

        return end_co_list

    def _get_half_co_list(self) -> list:
        half_co_list = []

        for co in self.all_co_tuples_list:
            if (co not in self.end_co_list) and (co not in self.full_co_list_seperate):
                half_co_list.append(co)

        return half_co_list

    def get_n_scaf_staple_co_types(self):
        """ dependes on:
                "full": self.full_co_list,
                "half": self.half_co_list,
                "end":  self.end_co_set
            """
        data = {"scaffold": dict(), "staple": dict()}
        types = {"full": self.full_co_list_seperate,
                 "half": self.half_co_list,
                 "end": self.end_co_list,
                 }

        for typ, crossovers in types.items():
            co_subsets = {"scaffold": {"": set(), "-h": set(), "-v": set()},
                          "staple": {"": set(), "-h": set(), "-v": set()}}
            for co in crossovers:
                for base in co:
                    strand = "scaffold" if base.is_scaf else "staple"
                    co_subsets[strand][""].add(co)
                    if co in self.all_co_tuples_h:
                        co_subsets[strand]["-h"].add(co)
                    else:
                        co_subsets[strand]["-v"].add(co)

            for s, direction_sets in co_subsets.items():  # scaffold, staple
                for dir in direction_sets:  # h, v
                    len_subset = len(co_subsets[s][dir])
                    n_co = len_subset / 2 if typ == "full" else len_subset
                    data[s][typ + dir] = n_co

        for strand in ["scaffold", "staple"]:
            data[strand]["co"] = data[strand]["half"] + data[strand]["full"]
            for typ in ["v", "h"]:
                data[strand]["co-" + typ] = (data[strand]["half-" + typ]
                                             + data[strand]["full-" + typ])
        return data

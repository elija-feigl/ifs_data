#!/usr/bin/env python
# -*- coding: utf-8 -*-3

from utils import get_statistics
import nanodesign as nd
import numpy as np
import os
import operator
from pathlib import Path
import argparse

from nanodesign.converters import Converter
from classes.crossover import Crossover
pass


class DesignData(object):

    def __init__(self, json: str, name: str):
        self.json: str = json
        self.name: str = name
        self.dna_structure, self.dna_structure_skips = self.init_design()
        self.all_strands: list = self.dna_structure.strands
        self.all_bases: list = self.get_all_bases()
        self.all_staples: list = self.get_all_staple()
        self.hps_base = self.init_hps()
        self.data: dict = {}
        self.n_st_domains = self.get_staple_domain()

        # crossover
        self.all_co = self._get_all_co()
        self.full_co_list_seperate, self.full_co_list_packed = self._get_full_co_list()
        self.end_co_list = self._get_endloop()
        self.half_co_list = self._get_half_co()
        self.all_crossovers, self.full_crossovers, self.half_crossovers, self.endloops = self.creat_crossover_lists()
        self.stacks = self.get_stacks()

        self.st_helix_dict: dict = self.init_helix_dict()
        self.first_bases, self.last_bases = self._get_first_last_bases_of_strands()
        self.helices = self.dna_structure.structure_helices_map
        self.nicks: list = self.get_nicks()
        self.long_domains = self.get_staples_with_long_domains()
        self.blunt_ends = self.get_blunt_ends()

    def compute_data(self) -> None:
        data = {}
        data["name"] = self.name
        data["lattice_type"] = self.get_lattice_type()
        data["n_helices"] = len(self.dna_structure.structure_helices_map)
        data["n_skips"] = self.get_n_skips()
        data["n_nicks"] = len(self.nicks)
        data["stacks"] = len(self.get_stacks())
        data["n_stacks"] = self.get_n_stacks()
        data.update(self.get_insertion_deletion_density())
        data["n_blunt_ends"] = len(self.get_blunt_ends())

        # staple stats
        data["n_staples"] = len(self.get_all_staple())
        data["staple_length"] = self.get_staples_length()
        data["helices_staples_pass"] = list(self.init_helix_dict().values())

        # domains
        data["n_staple_domain"] = self.get_staple_domain()
        data["long_domains"] = self.get_staples_with_long_domains()
        data.update(self.divide_domain_lengths())

        # crossovers
        data["co_set"] = self.classify_crossovers()
        data["co_possible"], data["co_density"] = self.get_co_density()

        self.data = data

    def init_design(self):
        # seq_file = self.name + ".seq"
        seq_name = None
        converter = Converter(modify=True, logg=False)
        converter.read_cadnano_file(self.json, None, seq_name)
        converter.dna_structure.compute_aux_data()
        converter_skip = Converter(modify=False, logg=False)
        converter_skip.read_cadnano_file(self.json, None, seq_name)
        converter_skip.dna_structure.compute_aux_data()
        return converter.dna_structure, converter_skip.dna_structure

    def _close_strand(self, strand):
        """[closes the scaffold andmaking it a loop]

        Returns:
            [strand] -- [closed scaffold strand]
        """
        start = strand.tour[0]
        end = strand.tour[-1]
        start.up, end.down = end, start
        self.dna_structure.strands[start.strand].is_circular = True
        return strand

    def init_hps(self) -> dict:
        hps_base = {}
        for strand in self.all_strands:
            for base in strand.tour:
                position = (base.h, base.p, base.is_scaf)
                hps_base[position] = base
        return hps_base

    def get_base_from_hps(self, h, p, is_scaffold, dir=1):
        """[get the base object from its coordination: (h, p is_scaffold)]

        Arguments:
            h {[int]} -- [helix]
            p {[int]} -- [position]
            is_scaffold {bool} -- [scaffold = True or staple = False]

        Keyword Arguments:
            dir {int} -- [direction] (default: {1})

        Returns:
            [base] -- [base of the giver coordination]
        """

        if (h, p) in self.dna_structure.Dhp_skips:
            p += np.sign(dir)

        return self.hps_base.get((h, p, is_scaffold), None)

    def get_base_plus_minus(self, base):
        """[given a base, it returns a neighbour base]

        Returns:
            [base_plus] -- [base with one position up along the helix]
            [base_minus] -- [base with on position down along the helix]
        """
        base_plus = self.get_base_from_hps(base.h, base.p + 1, base.is_scaf)
        base_minus = self.get_base_from_hps(
            base.h, base.p - 1, base.is_scaf, dir=-1)

        return base_plus, base_minus

    def get_lattice_type(self):
        if type(self.dna_structure.lattice) == nd.data.lattice.SquareLattice:
            return "Square"
        else:
            return "Honeycomb"

    def get_all_bases(self) -> list:
        all_bases = list()
        for strand in self.all_strands:
            for base in strand.tour:
                all_bases.append(base)
        return all_bases

    def get_staples_length(self) -> list:
        len_strands = list()
        for strand in self.all_strands:
            if not strand.is_scaffold:
                tour_clean = [base for base in strand.tour]
                len_strands.append(len(tour_clean))
        return len_strands

    def get_n_skips(self) -> int:
        return len(self.dna_structure.Dhp_skips)

    def get_staple_domain(self) -> list:
        data = {}
        domain_data = {}
        n_st_domains = list()  # number of domains for each staple
        for strand in self.all_strands:
            if not strand.is_scaffold:
                data = {str(strand.id): strand.domain_list}
                domain_data.update(data)
                data = {}
        for strand in domain_data.keys():
            n_st_domains.append(len(domain_data[strand]))
        # domin_data : is a dict that map the strand ID to the domains it has
        return n_st_domains

    def get_staples_with_long_domains(self) -> list:
        """[long domain are domains with 14 and more bases]

        Returns:
            list -- [numbers of long_domains for the list of all strands]
        """
        n_long_domains = list()
        for strand in self.all_strands:
            if not strand.is_scaffold:
                dummy = list()
                for domain in strand.domain_list:
                    if len(domain.base_list) >= 14:
                        dummy.append(strand)
                n_long_domains.append(len(dummy))

        return n_long_domains

    def divide_domain_lengths(self) -> dict:
        long_st_domain = list()
        domain_unpaired = list()

        data = {
            "2_long_domains": list(),
            "1_long_domains": list(),
            "0_long_domains": list(),
            "co_rule_violation": list()
        }

        for strand in self.all_staples:
            for domain in strand.domain_list:
                if len(domain.base_list) >= 14:
                    long_st_domain.append(strand)

                for base in domain.base_list:
                    if base.across is None:
                        domain_unpaired.append(domain)
                        break

            if len(long_st_domain) >= 2:
                data["2_long_domains"].append(strand)
            elif len(long_st_domain) == 1:
                data["1_long_domains"].append(strand)
            elif len(long_st_domain) == 0:
                data["0_long_domains"].append(strand)
            long_st_domain = list()

        for strand in self.all_staples:
            for domain in strand.domain_list:
                if domain not in domain_unpaired:
                    if len(domain.base_list) < 5:
                        data["co_rule_violation"].append(domain)

        return data

    def get_all_staple(self) -> int:
        all_staples = list()
        for strand in self.all_strands:
            if not strand.is_scaffold:
                all_staples.append(strand)
        return all_staples

    def init_helix_dict(self) -> dict:
        """
        [creates a dict with staple ID as key and the number of helices that it passes through as values]
        """
        st_helix_dict = {}
        helices_list = list()
        for strand in self.all_staples:
            helices = set()
            helices.add(strand.tour[0].h)
            for base in strand.tour:
                if self.dna_structure._check_base_crossover(base):
                    if base.h not in helices:
                        helices.add(base.h)
            helices_list.append(len(helices))

            if not strand.is_scaffold:
                st_helix_dict.update({str(strand.id): helices_list[-1]})

        return st_helix_dict

    def _get_first_last_bases_of_strands(self) -> list:
        first_bases = set()
        last_bases = set()
        for strand in self.all_strands:
            if not strand.is_scaffold:
                first_bases.add(strand.tour[0])
                last_bases.add(strand.tour[-1])

        return first_bases, last_bases

    def get_nicks(self) -> int:
        nicks = list()
        for base in self.first_bases:
            base_plus = self.get_base_from_hps(
                base.h, base.p + 1, base.is_scaf)
            base_minus = self.get_base_from_hps(
                base.h, base.p - 1, base.is_scaf, dir=-1)
            if base_plus in self.last_bases:
                # order of nick is always (first,last)
                nicks.append((base, base_plus))
            elif base_minus in self.last_bases:
                nicks.append((base, base_minus))

        return nicks

    def creat_crossover_lists(self):
        """[this method creates crossover objects in the dna origami structure]

        Returns:
            [list] -- [all_crossovers: list of all crossover object]
            [list] -- [full_crossovers: list of all full crossover object]
            [list] -- [half_crossovers: list of all half crossover object]
            [list] -- [endloops: list of all endloop object]
        """
        all_crossovers = list()
        endloops = list()
        full_crossovers = list()
        half_crossovers = list()

        for full_co in self.full_co_list_packed:
            typ = 'full'
            full_crossover = self.creat_crossover(typ, full_co)
            full_crossovers.append(full_crossover)
            all_crossovers.append(full_crossover)

        for ends in self.end_co_list:
            typ = 'end'
            end = self.creat_crossover(typ, ends)
            endloops.append(end)
            all_crossovers.append(end)

        for half in self.half_co_list:
            typ = 'half'
            half = self.creat_crossover(typ, half)
            half_crossovers.append(half)
            all_crossovers.append(half)

        return all_crossovers, full_crossovers, half_crossovers, endloops

    def creat_crossover(self, typ, co) -> list:
        """[this function creats crossover objects with attributes which are available in crossover class]

        Arguments:
            typ {[str]} -- [crossover type: full, half, endloop]
            co {[tuple]} -- [a resemblance of a crossover made of a tuple consisting two basis which indicates a crossover]

        Returns:
            crossover Object -- [description]
        """
        crossovers = list()
        coordinate = list()
        helix_row = list()

        if typ == 'full':
            base_typ = co[0][0]
            co_typ = co[0]
            for bases in co:
                for base in bases:
                    coordinate.append([base.p, base.h])
        else:
            base_typ = co[0]
            co_typ = co

            for base in co:
                coordinate.append((base.p, base.h))

        if base_typ.is_scaf:
            strand_typ = 'scaffold'
        else:
            strand_typ = 'staple'

        for base in co_typ:
            map_id_helices = self.dna_structure.structure_helices_map
            helix_row.append(
                map_id_helices[base.h].lattice_row)

        if helix_row[0] == helix_row[1]:
            is_vertical = False
        else:
            is_vertical = True

        crossover = Crossover(typ,
                              strand_typ, is_vertical, coordinate, co)

        return crossover

    def _get_all_co(self) -> list:
        """[get a list of crossovers but not as objects but as a tuple of the two bases connected via a crossover]

        Returns:
            list -- [list of all crossovers]
        """
        all_co_tuples = set()
        all_co = list()
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
            all_co.append(co)

        return all_co

    def _get_full_co_list(self) -> list:
        """[gets the full crossovers as tuples of bases]

        Returns:

            list -- [full_co_list_packed: get full co as a pack of two crossover (representation: [Co[B,B],Co[B,B], all in lists)
                     full_co_list_seperate: also seperately as individual crossovers (every Co in frozenset of two bases)]
        """

        full_co_list = list()
        for co in self.all_co:
            co_neighbours = {"co_tuples_plus": tuple(),
                             "co_tuples_minus": tuple()
                             }
            co_tuple_plus = set()
            co_tuple_minus = set()

            for base in co:
                base_plus, base_minus = self.get_base_plus_minus(base)
                co_tuple_plus.add(base_plus)
                co_tuple_minus.add(base_minus)
            co_neighbours["co_tuples_plus"] = tuple(co_tuple_plus)
            co_neighbours["co_tuples_minus"] = tuple(co_tuple_minus)

            fullco = set()
            for typ in ["co_tuples_plus", "co_tuples_minus"]:
                if co_neighbours[typ] in self.all_co:
                    fullco.add(co)
                    fullco.add(co_neighbours[typ])
                    full_co_list.append(frozenset(fullco))

                co_neighbours[typ] = list()

        # putting all full_co in a list configuration as [(Co(B,B),Co(B,B))] two parallel Co in a tuple and two bases also in a tuple

        full_co_set = set(full_co_list)
        full_co_list_seperate = list()
        full_co_list_packed = list()

        for full in full_co_set:
            full_co_list_packed.append(tuple(full))
            for co in full:
                full_co_list_seperate.append(co)

        return full_co_list_seperate, full_co_list_packed

    def _get_endloop(self) -> list:
        end_co_list = list()
        for co in self.all_co:
            for base in co:
                base_plus = self.get_base_from_hps(
                    base.h, base.p + 1, base.is_scaf)
                base_minus = self.get_base_from_hps(
                    base.h, base.p - 1, base.is_scaf, dir=-1)

                if (base_plus is None) or (base_minus is None):
                    if co not in end_co_list:
                        end_co_list.append(co)

        return end_co_list

    def _get_half_co(self) -> list:
        half_co_list = list()
        for co in self.all_co:
            if (co not in self.end_co_list) and (co not in self.full_co_list_seperate):
                half_co_list.append(co)

        return half_co_list

    def classify_crossovers(self):
        data = {"scaffold": dict(), "staple": dict()}
        types = {"full": self.full_crossovers,
                 "half": self.half_crossovers,
                 "end": self.endloops,
                 }

        for typ, crossovers in types.items():
            co_subsets = {"scaffold": {"": set(), "-h": set(), "-v": set()},
                          "staple": {"": set(), "-h": set(), "-v": set()}}

            for co in crossovers:
                strand = "scaffold" if co.strand_typ == 'scaffold' else "staple"
                co_subsets[strand][""].add(co)
                if not co.is_vertical:
                    co_subsets[strand]["-h"].add(co)
                else:
                    co_subsets[strand]["-v"].add(co)

                for s, direction_sets in co_subsets.items():  # scaffold, staple
                    for dir in direction_sets:  # h, v
                        len_subset = len(co_subsets[s][dir])
                        data[s][typ + dir] = len_subset

        for strand in ["scaffold", "staple"]:
            data[strand]["co"] = data[strand]["half"] + data[strand]["full"]
            for typ in ["v", "h"]:
                data[strand][
                    "co-" + typ] = (data[strand]["half-" + typ] + data[strand]["full-" + typ])

        return data

    def get_insertion_deletion_density(self):
        data = {"del_density": 0,
                "ins_density": 0}
        base_ins = 0
        for strand in self.dna_structure_skips.strands:
            for base in strand.tour:
                if base.num_insertions != 0:
                    base_ins += 1

        data["del_density"] = len(
            self.dna_structure.Dhp_skips) / len(self.all_bases)
        data["ins_density"] = base_ins / len(self.all_bases)

        return data

    def get_stacks(self):
        """[get the list of all stacks in the structure]

        return  [stacks: a list of stacks
                n_stacks: a list of length of each stack ]
        """

        stacks = list()
        added = list()
        J = 0
        K = 0

        full_packed_co = list()
        same_pos = list()
        dummy = list()
        full_packed_co.extend(self.full_co_list_packed)
        for f in full_packed_co:
            added.append(False)

        for full in full_packed_co:
            if added[J] is False:
                K = 0
                for full_1 in full_packed_co:
                    if (full != full_1) and (added[K] is False) and (np.abs(full[0][0].p - full_1[0][0].p) <= 3):
                        dummy.append(full_1)
                        added[K] = True
                    K = K + 1
                if len(dummy) >= 1:
                    dummy.append(full)
                    added[J] = True
                # dummy = tuple(dummy)
                same_pos.append(dummy)
                dummy = list()
            J = J + 1

        def common(lst1, lst2):
            return list(set(lst1) & set(lst2))

        def checker(full, group, dummy):
            h = list()
            h_1 = list()
            dummy = list()

            for co in full:
                for base in co:
                    if base.h not in h:
                        h.append(base.h)

                for full_1 in group:
                    for co in full_1:
                        for base in co:
                            if base.h not in h_1:
                                h_1.append(base.h)

                    if len(common(h, h_1)) == 1:
                        if full_1 not in dummy:
                            dummy.append(full_1)
                        else:
                            pass
                        if full not in dummy:
                            dummy.append(full)
                        else:
                            pass

                    h_1 = list()
            h = list()

            return dummy

        dummy = list()
        checked = set()

        for group in same_pos:
            for f in group:
                if f not in checked:
                    dummy.extend(checker(f, group, dummy))
                    checked.add(f)

                    if len(dummy) >= 1:
                        n = 0
                        while n < len(group):
                            for ff in dummy:
                                if ff not in checked:
                                    dummy.extend(checker(ff, group, dummy))
                                    checked.add(ff)
                            n = n + 1
                        n = 0

                        dummy = tuple(set(tuple(dummy)))
                        if dummy not in stacks:
                            stacks.append(dummy)

                        dummy = list()

                    else:
                        checked.add(f)
                        continue

            checked = set()

        return stacks

    def get_n_stacks(self):
        n_stacks = list()
        for stack in self.stacks:
            n_stacks.append(len(stack))
        return n_stacks

    def get_co_density(self):
        def is_ds(pos, hid):
            is_sc = (hid, pos, True) in self.hps_base
            is_st = (hid, pos, False) in self.hps_base
            # (hid, pos) in self.skips (note: list of (h,p) for all skips)
            is_skip = False
            return ((is_sc and is_st) or is_skip)

        def cleanup_co(co_list):
            n_ends = 0
            if not co_list:
                return 0, 0
            if len(co_list) == 1:
                return 1, 0
            if len(co_list) == 2 and co_list[0] != co_list[1]:
                return 2, 0

            if co_list[0] + 1 != co_list[1]:
                n_ends += 1
                co_list = co_list[1:]
            if co_list[-1] - 1 != co_list[-2]:
                n_ends += 1
                co_list = co_list[:-1]
            return n_ends, len(co_list) / 2

        possible_crossovers = {"scaffold": {"co": 0, "co-h": 0, "co-v": 0, "end": 0},
                               "staple": {"co": 0, "co-h": 0, "co-v": 0, "end": 0}
                               }
        # part 1: number of possible crossovers
        helices = self.dna_structure.structure_helices_map.values()

        for helix in helices:
            helix_row = helix.lattice_row

            for strand in ["scaffold", "staple"]:
                for typ in ["v", "h"]:
                    if strand == "scaffold":
                        p_co = helix.possible_scaffold_crossovers
                    else:
                        p_co = helix.possible_staple_crossovers

                    if typ == "h":
                        x = [co[1] for co in p_co
                             if (is_ds(pos=co[1], hid=helix.id) and (helix_row == co[0].lattice_row))
                             ]
                    else:
                        x = [co[1] for co in p_co
                             if (is_ds(pos=co[1], hid=helix.id) and (helix_row != co[0].lattice_row))
                             ]
                    end, co = cleanup_co(sorted(x))
                    possible_crossovers[strand]["co"] += co
                    possible_crossovers[strand]["co-" + typ] += co
                    possible_crossovers[strand]["end"] += end

        # part2 get actual crossovers
        set_crossovers = self.classify_crossovers()

        co_density = dict()
        for strand in ["scaffold", "staple"]:
            co_density[strand] = dict()
            for typ, n_possible in possible_crossovers[strand].items():
                n_set = set_crossovers[strand][typ]
                if n_possible == 0:
                    co_density[strand][typ] = 0
                else:
                    co_density[strand][typ] = n_set / n_possible

        return possible_crossovers, co_density

    def get_blunt_ends(self):
        blunt_ends = set()
        for co in self.end_co_list:
            if co[0].is_scaf:
                if (co[0].across is True) and (co[1].across is True):
                    for base in co:
                        if (base.across in self.first_bases) or (base.across in self.last_bases):
                            blunt_ends.add(co)
        return blunt_ends

    def prep_data_for_export(self) -> dict:
        export = dict()
        for name, value in self.data.items():
            if name in ["co_set", "co_possible", "co_density"]:
                for strand_name, subtypes in value.items():
                    for typ, n_co in subtypes.items():
                        export["{}-{}-{}".format(name,
                                                 strand_name, typ)] = n_co

            elif name in ["staple_length", "helices_staples_pass", "n_staple_domain", "long_domains", "n_stacks"]:
                stats = get_statistics(value, name)
                for stat_name, stat in stats.items():
                    export[stat_name] = stat
            elif name in ["2_long_domains", "1_long_domains", "0_long_domains", "co_rule_violation"]:
                export[name] = len(value)
            else:
                export[name] = value
        return export


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input",
                        help="input file",
                        type=str,
                        default="TTcorr.json",
                        )
    args = parser.parse_args()
    json = Path(args.input)
    outname = "{}-stat.csv".format(json.name)
    designdata = DesignData(json=json, name=json.name)
    designdata.compute_data()
    data = designdata.prep_data_for_export()
    with open(outname, mode="w+") as outfile:
        header = ",".join(str(k) for k in data.keys())
        outfile.write(header + "\n")
        export = ",".join(str(v) for v in data.values())
        outfile.write(export + "\n")


if __name__ == "__main__":
    main()

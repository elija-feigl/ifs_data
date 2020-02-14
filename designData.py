#!/usr/bin/env python
from nanodesign.converters import Converter
import nanodesign as nd
import numpy as np
import re
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt


class DesignData(object):

    def __init__(self, json: str, name: str):
        self.json: str = json
        self.name: str = name
        self.dna_structure, self.dna_structure_skips = self.init_design()
        self.all_strands: list = self.dna_structure.strands
        self.all_staples: list = self.get_all_staple()
        self.data: dict = {}
        self.all_bases: list = self.get_all_bases()
        self.hps_base: dict = self.init_hps()
        self.n_st_domains = self.get_staple_domain()
        # crossover
        self.all_co_tuples_list: list = self._get_all_co_tuple()
        self.all_co_tuples_h, self.all_co_tuples_v = self._get_horozantal_vertical_co()
        self.full_co_list_seperate, self.full_co_list_packed, self.full_co_set = self._get_full_co_list()
        self.end_co_set = self._get_endloop_co_list()
        self.half_co_list = self._get_half_co_list()
        self.stacks = self.get_stacks()

        self.st_helix_dict: dict = self.init_helix_dict()
        self.first_bases, self.last_bases = self._get_first_last_bases_of_strands()
        self.helices = self.dna_structure.structure_helices_map
        self.nicks: list = self.get_nicks()
        self.long_domains = self.get_staples_with_long_domains()
        self.blunt_ends = self.get_blunt_ends()

    def compute_data(self) -> dict:
        data = {}
        data["name"] = self.name
        data["Lattice type"] = self.get_lattice_type()
        data["n_helices"] = len(self.dna_structure.structure_helices_map)
        data["n_skips"] = self.get_n_skips()
        data["n_nicks"] = len(self.nicks)
        # domains
        data["n_staple_domain"] = self.get_staple_domain()
        data["long_domains"] = self.get_staples_with_long_domains()
        data.update(self.divide_domain_lengths())
        # staple stats
        data["n_staples"] = len(self.get_all_staple())
        data["staple_length"] = self.get_staples_length()
        data["helices_staples_pass"] = list(self.init_helix_dict().values())
        # crossovers
        data["co_set"] = self.get_n_scaf_staple_co_types()
        data["co_possible"], data["co_density"] = self.get_co_density()
        data["stacks"] = len(self.get_stacks())
        data["n_stacks"] = self.get_n_stacks()

        data.update(self.get_insertion_deletion_density())

        data["n_blunt_ends"] = len(self.get_blunt_ends())

        self.data = data
        return self.data

    def init_design(self):
        # seq_file = self.name + ".seq"
        seq_name = None
        converter = Converter(modify=True)
        converter.read_cadnano_file(self.json, None, seq_name)
        converter.dna_structure.compute_aux_data()
        converter_skip = Converter(modify=False)
        converter_skip.read_cadnano_file(self.json, None, seq_name)
        converter_skip.dna_structure.compute_aux_data()
        return converter.dna_structure, converter_skip.dna_structure

    def init_hps(self) -> dict:
        hps_base = {}
        for strand in self.all_strands:
            for base in strand.tour:
                position = (base.h, base.p, base.is_scaf)
                hps_base[position] = base
        return hps_base

    def _close_strand(self, strand):
        """[closes the scaffold andmaking it a loop]

        Returns:
            [strand] -- [description]
        """
        start = strand.tour[0]
        end = strand.tour[-1]
        start.up, end.down = end, start
        self.dna_structure.strands[start.strand].is_circular = True
        return strand

    def get_base_from_hps(self, h, p, is_scaffold, dir=1):
        if (h, p) in self.dna_structure.Dhp_skips:
            p += np.sign(dir)
        return self.hps_base.get((h, p, is_scaffold), None)

    def get_lattice_type(self):
        if type(self.dna_structure.lattice) == nd.data.lattice.SquareLattice:
            return "Square"
        else:
            return "Honeycomb"

    def get_all_bases(self) -> list:
        all_bases = []
        for strand in self.all_strands:
            for base in strand.tour:
                all_bases.append(base)
        return all_bases

    def get_staples_length(self):
        len_strands = []
        for strand in self.all_strands:
            if not strand.is_scaffold:
                tour_clean = [base for base in strand.tour]
                len_strands.append(len(tour_clean))
        return len_strands

    def get_structure_size(self):
        # TODO:...
        return

    def get_n_skips(self) -> int:
        return len(self.dna_structure.Dhp_skips)

    def get_staple_domain(self) -> list:
        data = {}
        domain_data = {}
        n_st_domains = []  # number of domains for each staple
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
        n_long_domains = []
        for strand in self.all_strands:
            if not strand.is_scaffold:
                dummy = []
                for domain in strand.domain_list:
                    if len(domain.base_list) >= 14:
                        dummy.append(strand)
                n_long_domains.append(len(dummy))

        return n_long_domains

    def divide_domain_lengths(self) -> dict:
        long_st_domain = []
        domain_unpaired = []

        data = {
            "2_long_domains": [],
            "1_long_domains": [],
            "0_long_domains": [],
            "co_rule_violation": []
        }

        for strand in self.all_staples:
            for domain in strand.domain_list:
                if len(domain.base_list) >= 14:
                    long_st_domain.append(strand)
            for domain in strand.domain_list:
                for base in domain.base_list:
                    if base.across == None:
                        domain_unpaired.append(domain)
                        break

#TODO###
            if len(long_st_domain) >= 2:
                data["2_long_domains"].append(strand)
            elif len(long_st_domain) == 1:
                data["1_long_domains"].append(strand)
            elif len(long_st_domain) == 0:
                data["0_long_domains"].append(strand)
            long_st_domain = []

        for strand in self.all_staples:
            for domain in strand.domain_list:
                if not domain in domain_unpaired:
                    if len(domain.base_list) < 5:
                        data["co_rule_violation"].append(domain)

        return data

    def get_all_staple(self) -> int:
        all_staples = []
        for strand in self.all_strands:
            if not strand.is_scaffold:
                all_staples.append(strand)
        return all_staples

    def init_helix_dict(self) -> dict:
        """
        [creates a dict with staple ID as key and the number of helices that it passes through as values]
        """
        st_helix_dict = {}
        helices_list = []
        for strand in self.all_staples:
            helices = set()
            helices.add(strand.tour[0].h)
            for base in strand.tour:
                if self.dna_structure._check_base_crossover(base):
                    if base.h not in helices:
                        helices.add(base.h)
            helices_list.append(len(helices))

        for strand, helices in zip(self.all_strands, helices_list):
            if not strand.is_scaffold:
                dic = {str(strand.id): helices}
                st_helix_dict.update(dic)

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
        nicks = []
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

    def _get_all_co_tuple(self) -> list:
        all_co_tuples = set()
        all_co_tuple_list = []
        for strand in self.all_strands:
            if strand.is_scaffold:
                self.all_strands.remove(strand)
                strand = self._close_strand(strand)
                self.all_strands.append(strand)

        for strand in self.all_strands:
            for base in strand.tour:
                if self.dna_structure._check_base_crossover(base):
                    co_tuple = set()
                    if base.up.h != base.h:
                        co_tuple.add(base.up)
                        co_tuple.add(base)
                        all_co_tuples.add(frozenset(co_tuple))
                        co_tuple = set()
                    elif base.down.h != base.h:
                        co_tuple.add(base.down)
                        co_tuple.add(base)
                        all_co_tuples.add(frozenset(co_tuple))
                        co_tuple = set()
        for co in all_co_tuples:
            all_co_tuple_list.append(co)

        return all_co_tuple_list

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
        co_plus_tuples = []
        co_minus_tuples = []
        full_co_list = []
        fullco = set()
        for co in self.all_co_tuples_list:
            co_tuple_plus = set()
            co_tuple_minus = set()
            for base in co:
                base_plus = self.get_base_from_hps(
                    base.h, base.p + 1, base.is_scaf)
                base_minus = self.get_base_from_hps(
                    base.h, base.p - 1, base.is_scaf, dir=-1)
                co_tuple_plus.add(base_plus)
                co_tuple_minus.add(base_minus)
            co_plus_tuples.append(frozenset(co_tuple_plus))
            co_minus_tuples.append(frozenset(co_tuple_minus))

            if co_plus_tuples[0] in self.all_co_tuples_list:
                fullco.add(co)
                fullco.add(co_plus_tuples[0])
                full_co_list.append(frozenset(fullco))
                fullco = set()
            elif co_minus_tuples[0] in self.all_co_tuples_list:
                fullco.add(co)
                fullco.add(co_minus_tuples[0])
                full_co_list.append(frozenset(fullco))
                fullco = set()

            co_plus_tuples = []
            co_minus_tuples = []

        # putting all full_co in a list configuration as [(Co(B,B),Co(B,B))] two parallel Co in a tuple and two bases also in a tuple

        full_co_set = set(full_co_list)
        full_co_list_seperate = []
        for full in full_co_set:
            for co in full:
                full_co_list_seperate.append(co)

        full_co_list_packed = []
        full_co_tuple = []
        for full_co in full_co_set:
            for co in full_co:
                full_co_tuple.append(tuple(co))
            full_co_list_packed.append(tuple(full_co_tuple))
            full_co_tuple = []

        dummy = []
        dummy_m = []
        for full in full_co_list_packed:
            for co in full:
                for base in co:
                    dummy.append(base)
            dummy_m.append(dummy)
            dummy = []

        return full_co_list_seperate, full_co_list_packed, full_co_set

    def _get_endloop_co_list(self) -> list:
        end_co_set = set()
        for co in self.all_co_tuples_list:
            for base in co:
                base_plus = self.get_base_from_hps(
                    base.h, base.p + 1, base.is_scaf)
                base_minus = self.get_base_from_hps(
                    base.h, base.p - 1, base.is_scaf, dir=-1)

                if (base_plus is None) or (base_minus is None):
                    end_co_set.add(co)

        return end_co_set

    def _get_half_co_list(self) -> list:
        half_co_list = []

        for co in self.all_co_tuples_list:
            if (co not in self.end_co_set) and (co not in self.full_co_list_seperate):
                half_co_list.append(co)

        return half_co_list

    # def get_n_co_types(self) -> int:
    #    return len(self.full_co_list)/2., len(self.half_co_list), len(self.end_co_set)

    def get_n_scaf_staple_co_types(self):
        """ dependes on:
                self.full_co_list,
                 "half": self.half_co_list,
                 "end": self.end_co_set,
                  self.all_co_tuples_h:

            """
        data = {"scaffold": dict(), "staple": dict()}
        types = {"full": self.full_co_list_seperate,
                 "half": self.half_co_list,
                 "end": self.end_co_set,
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
                    n_co = len_subset/2 if typ == "full" else len_subset
                    data[s][typ + dir] = n_co

        for strand in ["scaffold", "staple"]:
            data[strand]["co"] = data[strand]["half"] + data[strand]["full"]
            for typ in ["v", "h"]:
                data[strand]["co-"+typ] = (data[strand]["half-" + typ]
                                           + data[strand]["full-"+typ])
        return data

    def get_insertion_deletion_density(self):
        data = {"del_density": 0,
                "ins_density": 0}
        base_ins = 0
        for strand in self.dna_structure_skips.strands:
            for base in strand.tour:
                if base.num_insertions == -1:
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

        stacks = []
        added = []
        J = 0
        K = 0

        full_packed_co = []
        same_pos = []
        dummy = []
        full_packed_co.extend(self.full_co_list_packed)
        for f in full_packed_co:
            added.append(False)

        for full in full_packed_co:
            if added[J] == False:
                K = 0
                for full_1 in full_packed_co:
                    if (full != full_1) and (added[K] == False) and (np.abs(full[0][0].p - full_1[0][0].p) <= 3):
                        dummy.append(full_1)
                        added[K] = True
                    K = K + 1
                if len(dummy) >= 1:
                    dummy.append(full)
                    added[J] = True
                # dummy = tuple(dummy)
                same_pos.append(dummy)
                dummy = []
            J = J + 1

        def common(lst1, lst2):
            return list(set(lst1) & set(lst2))

        def checker(full, group, dummy):
            h = []
            h_1 = []
            dummy = []

            for co in full:
                for base in co:
                    if not base.h in h:
                        h.append(base.h)

                for full_1 in group:
                    for co in full_1:
                        for base in co:
                            if not base.h in h_1:
                                h_1.append(base.h)

                    if len(common(h, h_1)) == 1:
                        if not full_1 in dummy:
                            dummy.append(full_1)
                        else:
                            pass
                        if not full in dummy:
                            dummy.append(full)
                        else:
                            pass

                    h_1 = []
            h = []

            return dummy

        dummy = []
        checked = set()

        for group in same_pos:
            for f in group:
                if not f in checked:
                    dummy.extend(checker(f, group, dummy))
                    checked.add(f)

                    if len(dummy) >= 1:
                        n = 0
                        while n < len(group):
                            for ff in dummy:
                                if not ff in checked:
                                    dummy.extend(checker(ff, group, dummy))
                                    checked.add(ff)
                            n = n+1
                        n = 0

                        dummy = tuple(set(tuple(dummy)))
                        if not dummy in stacks:
                            stacks.append(dummy)

                        dummy = []

                    else:
                        checked.add(f)
                        continue

            checked = set()

        return stacks

    def get_n_stacks(self):
        n_stacks = []
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

            if co_list[0]+1 != co_list[1]:
                n_ends += 1
                co_list = co_list[1:]
            if co_list[-1]-1 != co_list[-2]:
                n_ends += 1
                co_list = co_list[:-1]
            return n_ends, len(co_list)/2

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
                             if (is_ds(pos=co[1], hid=helix.id)
                                 and (helix_row == co[0].lattice_row)
                                 )
                             ]
                    else:
                        x = [co[1] for co in p_co
                             if (is_ds(pos=co[1], hid=helix.id)
                                 and (helix_row != co[0].lattice_row)
                                 )
                             ]
                    end, co = cleanup_co(sorted(x))
                    possible_crossovers[strand]["co"] += co
                    possible_crossovers[strand]["co-" + typ] += co
                    possible_crossovers[strand]["end"] += end

        # part2 get actual crossovers
        set_crossovers = self.get_n_scaf_staple_co_types()

        co_density = dict()
        for strand in ["scaffold", "staple"]:
            co_density[strand] = dict()
            for typ, n_possible in possible_crossovers[strand].items():
                n_set = set_crossovers[strand][typ]
                co_density[strand][typ] = n_set / n_possible

        return possible_crossovers, co_density

    def get_blunt_ends(self):
        blunt_ends = set()
        for co in self.end_co_set:
            for base in co:
                if base.is_scaf:
                    if (base.across in self.first_bases) or (base.across in self.last_bases):
                        blunt_ends.add(co)
        return blunt_ends


def get_statistics(data_list, data_name):
    """[summary]

    Arguments:
        data_list {[type]} -- [description]
        data_name {[type]} -- [description]

    Returns:
        [type] -- [description]
    """
    return {data_name + "_avg": np.average(data_list),
            data_name + "_std": np.std(data_list),
            data_name + "_max": np.max(data_list),
            data_name + "_min": np.min(data_list),
            }


def prep_data_for_export(data):
    # TODO: split lists and sets in statistics (avg, std, min, max)
    export = dict()
    for name, value in data.items():
        if name in ["co_set", "co_possible", "co_density"]:
            for strand_name, subtypes in value.items():
                for typ, n_co in subtypes.items():
                    export["{}-{}-{}".format(name, strand_name, typ)] = n_co

        elif name in ["staple_length", "helices_staples_pass", "n_staple_domain", "long_domains", "n_stacks"]:
            stats = get_statistics(value, name)
            for stat_name, stat in stats.items():
                export[stat_name] = stat
        elif name in ["2_long_domains", "1_long_domains", "0_long_domains", "co_rule_violation"]:
            export[name] = len(value)
        else:
            export[name] = value

    return export


"""
def export_data(data: dict, name: str) -> None:

    export = prep_data_for_export(data)
    header = ", ".join([str(i) for i in export.keys()])
    export_str = ", ".join([str(i) for i in export.values()])

    try:
        os.mkdir("out")
    except FileExistsError:
        pass
    with open("./out/" + name + "-designdata.csv", mode="w+") as out:

        out.write(header + "\n")
        out.write(export_str + "\n")
        # out.write("\nEND")
    return


def main():
    f = open(Path("./txt_file.txt"), 'rt', encoding="utf8")
    for line in f:
        if line.startswith('Project ='):
            name = line[9:-1].strip()
            break
    # print("master, I am awaiting the name of your design")
    # name = input()
    # print("Thank you Sir")
    designData = DesignData(name=name)
    data = designData.compute_data()
    export_data(data=data, name=name)
    return


if __name__ == "__main__":
    main()
"""

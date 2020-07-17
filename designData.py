#!/usr/bin/env python
# -*- coding: utf-8 -*-3

import nanodesign as nd
import numpy as np
from Bio.SeqUtils import MeltingTemp  # compute melting temperatures
from Bio.Seq import Seq
from classes.crossover import Crossover
from classes.nicks import Nick
from statistics import mean


class DesignData(object):

    def __init__(self, json: str, name: str, seq: str):

        self.json: str = json
        self.name: str = name
        self.seq: str = seq
        self.data = dict()
        self.dna_structure, self.dna_structure_skips = self.init_design()

        # strand, base
        self.strands: list = self.dna_structure.strands
        self._get_staple()
        self.scaf_bases: list = self.get_all_scaf_bases()
        self.hps_base = self._init_hps()
        self.hps_base_skips = self._init_hps_skips()
        self._get_first_last_bases_of_strands()

        # helix
        self.helices = self.dna_structure.structure_helices_map
        self._init_helix_data()

        # domain
        self._init_domain_data()
        self.long_domains = self.get_staples_with_long_domains()
        self.staple_domains_melt_t: dict = self.staple_domains_melt_t()
        # NOTE: dictionary of staples and maximum melting temp of domains for each staple
        self.max_staple_melt_t = {key: max(value) for (
            key, value) in self.staple_domains_melt_t.items()}

        # crossover
        self.all_co_sets, self.all_co_lists = self._get_all_connections()
        self.full_co_sets_seperate, self.full_co_tuples = self._get_full_co_list()
        self.end_co_sets = list()
        self.end_co_tuples = self._get_endloop()
        self.half_co_tuples = self._get_half_co()
        self.all_crossovers = list()
        self.full_crossovers = list()
        self.half_crossovers = list()
        self.endloops = list()
        self._create_crossover_lists()
        self.pos = self.full_scaff_type()

        # structure
        self.alpha_value = self.alpha_value()
        self.stacks = self.get_stacks()
        self.stacks_lengths = self.get_stacks_lengths()
        self.loops_length_list = self.get_loops()
        self.nicks: list = self.get_nicks()
        self.blunt_ends = self.get_blunt_ends()

    def init_design(self):
        converter = nd.converters.Converter(modify=True, logg=False)
        converter.read_cadnano_file(self.json, None, self.seq)
        converter.dna_structure.compute_aux_data()
        converter_skip = nd.converters.Converter(modify=False, logg=False)
        converter_skip.read_cadnano_file(self.json, None, self.seq)
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

    def _get_staple(self) -> list:
        self.staples = [
            strand for strand in self.strands if not strand.is_scaffold]

    def _init_hps(self) -> dict:
        """[create a dictionary of bases positions and the bases itself for all bases in the structure excluding the skips]

        Returns:
            dict: [keys: position = (base.h, base.p, base.is_scaf), values: base object]
        """
        hps_base = dict()
        for strand in self.dna_structure.strands:
            hps_base.update({(base.h, base.p, base.is_scaf)                             : base for base in strand.tour})

        return hps_base

    def _init_hps_skips(self) -> dict:
        """[create a dictionary of bases positions and the bases itself for all bases in the structure including the skips]

        Returns:
            dict: [keys: position = (base.h, base.p, base.is_scaf), values: base object]
        """
        hps_base_skips = dict()
        for strand in self.dna_structure_skips.strands:
            hps_base_skips.update(
                {(base.h, base.p, base.is_scaf): base for base in strand.tour})

        return hps_base_skips

    def _get_base_from_hps(self, h, p, is_scaffold, dir=1):
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
        """[given a base, it returns the neighbour bases]

        Returns:
            [base_plus] -- [base with one position up along the helix]
            [base_minus] -- [base with on position down along the helix]
        """
        base_plus = self._get_base_from_hps(base.h, base.p + 1, base.is_scaf)
        base_minus = self._get_base_from_hps(
            base.h, base.p - 1, base.is_scaf, dir=-1)

        return base_plus, base_minus

    def get_lattice_type(self):
        if type(self.dna_structure.lattice) == nd.data.lattice.SquareLattice:
            return "Square"
        else:
            return "Honeycomb"

    def get_dimension(self):
        """[getting the dimension of the structure.]

        Returns:
            [tuple]: [
                a: max length of the columns in CadNano
                b: max length of the rows in Cadnano
                c: diffrence of min base postion in the struction and the max base position in Cadnano (structure's depth)
            ]
        NOTE: the values could be wrong for exotic designs
        """

        lattice_cols = [helix.lattice_col for helix in self.helices.values()]
        lattice_rows = [helix.lattice_row for helix in self.helices.values()]
        base_pos = [base.p for base in self.scaf_bases]

        a = max(lattice_cols) - min(lattice_cols) + 1
        b = max(lattice_rows) - min(lattice_rows) + 1
        c = max(base_pos) - min(base_pos) + 1

        dimension = (a, b, c)

        return dimension

    def get_all_scaf_bases(self) -> list:

        all_scaf_bases = list()
        for strand in self.strands:
            all_scaf_bases.extend(
                [base for base in strand.tour if strand.is_scaffold])

        return all_scaf_bases

    def get_staples_length(self) -> list:
        """[creates a list of staples length]

        Returns:
            list: [staples length]
        """
        staples_length = list()
        for staple in self.staples:
            staples_length.append(len([base for base in staple.tour]))

        return staples_length

    def _init_domain_data(self) -> dict:
        """[creates a three dictionary of staples and their domains informations]

        Returns:
            dict -- domain_data:  [dictionary of staples and list of their domains]
            dict -- domain_lengths_data:  [dictionary of staples and list of their domains lengths]
            dict -- n_staples_domains:  [dictionary of staples and number of their domains]
        """
        self.domain_data = dict()
        self.domain_lengths_data = dict()
        self.n_staples_domains = dict()

        self.domain_data.update(
            {staple: staple.domain_list for staple in self.staples})

        for staple, domains in self.domain_data.items():
            self.domain_lengths_data.update(
                {staple: [len(domain.base_list) for domain in domains]})

        self.n_staples_domains = {staple: len(
            self.domain_data[staple]) for staple in self.domain_data.keys()}

    def staple_domains_melt_t(self) -> dict:
        """[staples domain melting temperature.]

        Returns:
            dict: [key: staple. values: list of domains' melting temperature.]
        """
        staple_domains_melt_t = {}

        for staple, domains in self.domain_data.items():
            for domain in domains:
                if "N" not in domain.sequence:
                    # NOTE: using nearest neighbor for domain with length higher than 14
                    if len(domain.base_list) > 14:
                        staple_domains_melt_t.setdefault(staple, []).append(MeltingTemp.Tm_NN(
                            Seq(domain.sequence), Na=0, Mg=17.5))

                    # NOTE: using 'Wallace rule' for domain with length less than 14
                    else:
                        staple_domains_melt_t.setdefault(staple, []).append(MeltingTemp.Tm_Wallace(
                            Seq(domain.sequence)))

        return staple_domains_melt_t

    def alpha_value(self):
        """[alpha value : The ratio of number of staples having doamins with melting
        temperature higher than critical temperature to the number of all staples in the structure]

        Returns:
            [dict]: [alpha values of the structure for different critical temperatures]
        """

        T_crit = {40: int, 55: int, 70: int}

        def calculate(self, T_crit):
            """
            [calculates alpha value for a given critical temperature]
            """
            return sum([True for T in list(self.max_staple_melt_t.values()) if T >= T_crit]) / len(self.max_staple_melt_t)

        alpha_values = {T: calculate(self, T) for T in T_crit}

        return alpha_values

    def get_staples_with_long_domains(self) -> dict:
        """[long domain are domains with 14 and more bases]

        Returns:
            dict -- [staples : numbers of long_domains for the each staple]
        """
        n_long_domains = dict()
        for staple, domains in self.domain_data.items():
            n_long_domains.setdefault(staple, 0)
            for domain in domains:
                if len(domain.base_list) > 14:
                    n_long_domains[staple] += 1

        return n_long_domains

    def divide_domain_lengths(self) -> dict:
        """[divide staples having 0, 1, 2 or more long domains ]

        Returns:
            dict: [
                2_long_domains: two or more long domains,
                1_long_domains: having only one long domain,
                0_long_domains: having no long domain,
                co_rule_violation: unpaired domains with less than 5 bases]
        """

        data = dict()
        domain_unpaired = list()

        for staple, n_longs in self.long_domains.items():
            if n_longs >= 2:
                data.setdefault("2_long_domains", []).append(staple)
            elif n_longs == 1:
                data.setdefault("1_long_domains", []).append(staple)
            elif n_longs == 0:
                data.setdefault("0_long_domains", []).append(staple)

        for staple, domains in self.domain_data.items():
            domain_unpaired.extend([
                domain for domain in domains if domain.base_list[0].across is None])

        for staple, domains in self.domain_data.items():
            data.setdefault("co_rule_violation", []).extend(
                [domain for domain in domains if domain not in domain_unpaired and len(domain.base_list) < 5])

        return data

    def _init_helix_data(self) -> dict:
        """
        [creates a two dictionary with staple and helices that it passes through]

        Returns:
            dict: [
                staple_helix_dict: [staple: helices ID the staple passes through],
                num_staple_helix_dict: [staple: number of helices the staple passes through],
        """
        self.staple_helix_dict = dict()
        self.num_staple_helix_dict = dict()
        for staple in self.staples:
            self.staple_helix_dict.update(
                {staple: tuple({base.h for base in staple.tour})})
            self.num_staple_helix_dict.update(
                {staple: len(tuple({base.h for base in staple.tour}))})

    def _get_first_last_bases_of_strands(self) -> list:

        self.first_bases = {staple.tour[0] for staple in self.staples}
        self.last_bases = {staple.tour[-1] for staple in self.staples}

    def get_nicks(self) -> int:
        # order of nick is always (first,last)
        nicks = list()

        def create_nick(base, neighbor_base):
            bases = (base, neighbor_base)
            p = (base.p, neighbor_base.p)
            h = (base.h, neighbor_base.h)
            nick = Nick(bases, set(p), set(h))
            return nick

        for base in self.first_bases:
            base_plus, base_minus = self.get_base_plus_minus(base)
            if base_plus in self.last_bases:
                nicks.append(create_nick(base, base_plus))
            if base_minus in self.last_bases:
                nicks.append(create_nick(base, base_minus))

        return nicks

    def _create_crossover_lists(self):
        """[this method creates crossover objects in the dna origami structure]

        Returns:
            [list] -- [all_crossovers: list of all crossover object]
            [list] -- [full_crossovers: list of all full crossover object]
            [list] -- [half_crossovers: list of all half crossover object]
            [list] -- [endloops: list of all endloop object]
        """

        crossovers_dict = {
            "full": self.full_co_tuples,
            "half": self.half_co_tuples,
            "end": self.end_co_tuples
        }

        def create_list(typ, list_name, all_crossovers, co_list):
            for co in co_list:
                crossover = Crossover(typ, co, self.helices)
                list_name.append(crossover)
                self.all_crossovers.append(crossover)

        for typ, crossovers in crossovers_dict.items():

            if typ == 'full':
                list_name = self.full_crossovers

            elif typ == 'end':
                list_name = self.endloops
            elif typ == 'half':
                list_name = self.half_crossovers

            create_list(typ, list_name, self.all_crossovers, crossovers)

    def full_scaff_type(self):
        """[type of the full scaffold crossover depending on the position suggested by cadnano]
        """
        for full in self.full_crossovers:
            if full.strand_typ == 'scaffold':
                sub_new = np.Infinity
                # find closest crossover
                # TODO: change to possibe staple co!!! or even just check p
                for co in self.all_crossovers:
                    if co.strand_typ == 'staple':

                        if co.h == full.h:
                            sub = (mean(full.p) - mean(co.p))

                            if abs(sub) <= abs(sub_new):
                                sub_new = sub

                if sub_new is np.Infinity:
                    full.scaff_full_type = 0
                    return

                # calculate type
                if self.get_lattice_type() == 'Square':
                    mod = sub_new % 32

                    if 0 <= mod <= 11:
                        typ = 1
                    elif 21 <= mod < 32:
                        typ = 3
                    else:
                        typ = 2

                else:
                    mod = sub_new % 21
                    if 0 <= mod <= 11:
                        typ = 1
                    elif 11 < mod < 21:
                        typ = 3

                full.scaff_full_type = typ

    def get_all_full_scaff_crossover_types(self):
        data = {
            'full_scaf_co_type_1': 0,
            'full_scaf_co_type_2': 0,
            'full_scaf_co_type_3': 0
        }
        for full in self.full_crossovers:
            if full.strand_typ == 'scaffold':
                if full.scaff_full_type == 1:
                    data['full_scaf_co_type_1'] += 1
                elif full.scaff_full_type == 2:
                    data['full_scaf_co_type_2'] += 1
                elif full.scaff_full_type == 3:
                    data['full_scaf_co_type_3'] += 1
        return data

    def _get_all_connections(self) -> list:
        """[get a list of connections but not as objects but as a tuple of the two bases connected via a crossover]

        Returns:
            list -- [list of all crossovers]
        """
        all_co_tuples = set()
        all_co_lists = list()
        all_co_sets = list()
        for strand in self.strands:
            if strand.is_scaffold:

                new_strand = self._close_strand(strand)
                self.strands[strand.id] = new_strand

        for strand in self.strands:
            for base in strand.tour:
                if self.dna_structure._check_base_crossover(base):
                    co_tuple = set()
                    if base.up.h != base.h:
                        co_tuple = (base, base.up)
                        all_co_tuples.add(tuple(set(co_tuple)))
                    elif base.down.h != base.h:
                        co_tuple = (base.down, base)
                        all_co_tuples.add(tuple(set(co_tuple)))

        for co in all_co_tuples:
            all_co_sets.append(set(co))
            all_co_lists.append(co)

        return all_co_sets, all_co_lists

    def _get_full_co_list(self) -> list:
        """[gets the full crossovers as tuples of bases]

        Returns:

            list -- [
                fulfull_packed_co.extend(self.full_co_tuples: get full co as a pack of two crossover
                (representation: [Co[B,B],Co[B,B], all in lists)
                full_co_list_seperate: also seperately as individual crossovers
                (every Co in frozenset of two bases)
                    ]
        """

        full_co_list = list()
        # lista = list()
        for co in self.all_co_sets:
            co_neighbours = dict()
            co_tuple_plus = set()
            co_tuple_minus = set()

            for base in co:
                base_plus, base_minus = self.get_base_plus_minus(base)
                co_tuple_plus.add(base_plus)
                co_tuple_minus.add(base_minus)

            co_neighbours["plus"] = set(co_tuple_plus)
            co_neighbours["minus"] = set(co_tuple_minus)

            fullco = set()
            for typ in ["plus", "minus"]:
                if co_neighbours[typ] in self.all_co_sets:
                    fullco.add(frozenset(co))
                    fullco.add(frozenset(co_neighbours[typ]))
                    full_co_list.append(frozenset(fullco))

                co_neighbours[typ] = list()

        """
        putting all full_co in a list configuration as [(Co(B,B),Co(B,B))]
        two parallel Co in a tuple and two bases also in a tuple
        """
        full_co_set = set(full_co_list)
        full_co_sets_seperate = list()
        full_co_tuples = list()

        for full_set in full_co_set:
            full = []
            for co in full_set:
                full.append(tuple(co))
            full_co_tuples.append(tuple(full))
            for co in full_set:
                full_co_sets_seperate.append(co)

        return full_co_sets_seperate, full_co_tuples

    def _get_endloop(self) -> list:
        for co in self.all_co_sets:
            for base in co:
                base_plus, base_minus = self.get_base_plus_minus(base)

                if (base_plus is None) or (base_minus is None):
                    if co not in self.end_co_sets:
                        self.end_co_sets.append(frozenset(co))

        end_tuples = list()
        for end in self.end_co_sets:
            # NOTE: we want to ensure bases has consistent type regardless of type
            end_tuples.append(tuple([tuple(end), None]))

        return end_tuples

    def _get_half_co(self) -> list:
        half_co_sets = list()
        half_co_tuples = list()
        for co in self.all_co_sets:
            if (co not in self.end_co_sets) and (co not in self.full_co_sets_seperate):
                half_co_sets.append(co)

        for co in half_co_sets:
            # NOTE: we want to ensure bases has consistent type regardless of type
            half_co_tuples.append(tuple([tuple(co), None]))

        return half_co_tuples

    def classify_crossovers(self):
        data = {"scaffold": dict(), "staple": dict()}
        types = {"full": self.full_crossovers,
                 "half": self.half_crossovers,
                 "end": self.endloops,
                 }

        for typ, crossovers in types.items():
            co_subsets = {"scaffold": {"": set(), "_h": set(), "_v": set()},
                          "staple": {"": set(), "_h": set(), "_v": set()}}

            for co in crossovers:
                strand = "scaffold" if co.strand_typ == 'scaffold' else "staple"
                co_subsets[strand][""].add(co)
                if co.orientation == "horizontal":
                    co_subsets[strand]["_h"].add(co)
                else:
                    co_subsets[strand]["_v"].add(co)

                for s, direction_sets in co_subsets.items():  # scaffold, staple
                    for dir in direction_sets:  # h, v
                        len_subset = len(co_subsets[s][dir])
                        data[s][typ + dir] = len_subset

        for strand in ["scaffold", "staple"]:
            data[strand]["co"] = data[strand]["half"] + data[strand]["full"]
            for typ in ["v", "h"]:
                data[strand][
                    "co_" + typ] = (data[strand]["half_" + typ] + data[strand]["full_" + typ])

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
            self.dna_structure.Dhp_skips) / len(self.scaf_bases)
        data["ins_density"] = base_ins / len(self.scaf_bases)

        return data

    def get_stacks(self):
        """[get the list of all stacks in the structure]

        return  [stacks: a list of stacks]
        """

        stacks = list()
        added = list()
        J = 0
        K = 0

        full_packed_co = list()
        same_pos = list()
        dummy = list()
        full_packed_co.extend(self.full_co_tuples)
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

    def get_stacks_lengths(self):
        stacks_lengths = list()
        for stack in self.stacks:
            stacks_lengths.append(len(stack))
        return stacks_lengths

    def get_co_density(self):
        """[calculate crossover density (number of crossovers is the structure divided by possible crossovers CadNano)]
        NOTE: the values for number of possible_co and co_desity are not exact but close to the true value
        """
        # TODO: possible_ co numbers are not exactly correct
        def is_ds(pos, hid):
            is_sc = (hid, pos, True) in self.hps_base_skips
            is_st = (hid, pos, False) in self.hps_base_skips
            # (hid, pos) in self.skips (note: list of (h,p) for all skips)
            is_skip = False

            return ((is_sc or is_st) or is_skip)

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
            # TODO: devision by two is assumed for possible_full_co(two connections)
            return n_ends, len(co_list) // 2

        def neighbour_bases(strand_typ, helix):
            """[gives a list of all bases in the neighbouring helix to ckeck if there
            exist a base in the neighbouring helix to connect to the base that is a possible_co in the main helix]

            Args:
                strand_typ ([str]): [Scaffold or Staple]
                helix ([DnaStructureHelix]): [neighbouring helix]

            Returns:
                [list]: [list of all the bases in the neighbouring helix]

            """
            if strand_typ == 'scaffold':
                neighbour_bases = [base.p for base in helix.scaffold_bases]
            else:
                neighbour_bases = [base.p for base in helix.staple_bases]

            return neighbour_bases

        def orientation(typ, helix, helix_row):
            if typ == 'h':
                return helix_row == helix.lattice_row
            else:
                return helix_row != helix.lattice_row

        possible_crossovers = {"scaffold": {"co": 0, "co_h": 0, "co_v": 0},
                               "staple": {"co": 0, "co_h": 0, "co_v": 0}
                               }
        # part 1: number of possible crossovers
        helices = self.dna_structure_skips.structure_helices_map.values()

        for helix in helices:
            helix_row = helix.lattice_row

            for strand in ["scaffold", "staple"]:
                for typ in ["v", "h"]:

                    # NOTE: nanodesign crossoevers are actually connections
                    if strand == "scaffold":
                        p_co = helix.possible_scaffold_crossovers

                    else:
                        p_co = helix.possible_staple_crossovers

                    x = [co[1] for co in p_co if (is_ds(pos=co[1], hid=helix.id)
                                                  and orientation(typ, co[0], helix_row)
                                                  and (co[1] in neighbour_bases(strand, co[0])))]

                    end, co = cleanup_co(sorted(x))
                    # TODO: devision by two is assumed for counting each possible_co two times for a helix and its neighbour

                    possible_crossovers[strand]["co"] += co // 2
                    possible_crossovers[strand]["co_" + typ] += co // 2
                    # possible_crossovers[strand]["end"] += end

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

        for end in self.endloops:
            has_across = (end.bases[0][0].across is True) and (
                end.bases[0][1].across is True)
            if end.strand_typ == 'scaffold' and has_across:
                blunt_ends = {end.bases[0] for base in end.bases[0] if (
                    base.across in self.first_bases) or (base.across in self.last_bases)}

        return blunt_ends

    def get_loops(self):
        loops = list()
        for co in self.full_crossovers + self.half_crossovers:
            sub = np.inf
            if co.strand_typ == 'scaffold':
                stacks = tuple([
                    tuple([co.bases[0][0], co.bases[1][0]]),
                    tuple([co.bases[0][1], co.bases[1][1]])
                ])
                for stack in stacks:
                    if stack[0].across is None or stack[1].across is None:
                        continue
                    same_staple = (stack[0].across.strand
                                   == stack[1].across.strand)
                    sc = self.strands[stack[0].strand]
                    same_scaffold = (sc.id == stack[1].strand)
                    if same_staple and same_scaffold:
                        # NOTE: potentially (stack[0].residue -1) istead of tour(index)
                        sub_new = abs(sc.tour.index(stack[0])
                                      - sc.tour.index(stack[1])
                                      )
                        if sub_new > len(sc.tour) / 2:
                            sub_new = len(sc.tour) - sub_new
                        if sub_new < sub:
                            sub = sub_new

            else:  # staple
                for connection in co.bases:
                    if connection is None:  # NOTE: only 1 for half_co
                        continue
                    if connection[0].across is None or connection[1].across is None:
                        continue
                    sc = self.strands[connection[0].across.strand]
                    same_scaffold = (sc.id == connection[1].across.strand)
                    if same_scaffold:
                        base_1 = sc.tour.index(connection[0].across)
                        base_2 = sc.tour.index(connection[1].across)
                        sub_new = abs(base_1 - base_2)

                        if sub_new > len(sc.tour) / 2:
                            sub_new = len(sc.tour) - sub_new
                        if sub_new < sub:
                            sub = sub_new
            if not np.isinf(sub):
                loops.append(sub)

        return loops

import attr
import logging
import numpy as np

from typing import List, Dict, Tuple, Set, Any

from nanodesign.data.base import DnaBase as Base
from nanodesign.data.strand import DnaStrand as Strand

from ..core.utils import get_statistics
from ..core.designData import Design


@attr.s
class DesignStats(object):
    design: Design = attr.ib()

    logger = logging.getLogger(__name__)

    def compute_data(self) -> None:
        design = self.design
        data: Dict[str, Any] = dict()
        data["name"] = design.name
        data["lattice_type"] = design.lattice
        data["dim_x"], data["dim_y"], data["dim_z"] = self.get_dimension()
        data['alpha_value'] = self.get_alpha_value()
        data["n_helices"] = self.get_n_helices()
        data["n_deletions"] = len(design.hps_deletions) // 2
        data["n_insertions"] = len(design.hps_insertions) // 2
        data["insertions_denisty"], data["deletion_denisty"] = self.get_insertion_deletion_density()

        data["n_nicks"] = self.get_number_nicks()
        data["n_stacks"] = len(self.get_stacks_lengths())
        data["stacks_length"] = self.get_stacks_lengths()
        data["loops_length"] = self.get_loops()

        data["n_bluntends"] = len(self.get_blunt_ends())

        # staple stats
        data["n_staples"] = len(design.staples)
        data["staples_length"] = self.get_staples_length()
        data["helices_staples_pass"] = self.get_num_staple_helix()

        # domains
        data["n_staples_domains"] = self.get_n_staples_domains()
        data["long_domains"] = list(
            self.get_staples_with_long_domains().values())
        data.update(self.divide_domain_lengths())
        data["staple_domain_melt_T"] = list(
            design.max_staple_melt_t.values())

        # crossovers
        data["co_set"] = self.classify_crossovers()
        data.update(self.get_full_scaff_co_typ_stat())
        data["co_possible"], data["co_density"] = self.get_co_density()

        self.data = data

    def get_full_scaff_co_typ_stat(self):
        """[numbers of ]

        Returns:
            [dict]: [data for database; writung to csv]
        """
        # TODO: the designprocess data are not consistant
        data = {
            'full_scaf_co_type_1': 0,
            'full_scaf_co_type_2': 0,
            'full_scaf_co_type_3': 0
        }
        for full in (co for co in self.design.crossovers if co.typ == "full"):
            if full.is_scaffold == 'scaffold':
                if full.scaff_full_type == 1:
                    data['full_scaf_co_type_1'] += 1
                elif full.scaff_full_type == 2:
                    data['full_scaf_co_type_2'] += 1
                elif full.scaff_full_type == 3:
                    data['full_scaf_co_type_3'] += 1
        return data

    def prep_data_for_export(self) -> dict:
        export = dict()
        for name, value in self.data.items():
            if name in ["co_set", "co_possible", "co_density"]:
                for strand_name, subtypes in value.items():
                    for typ, n_co in subtypes.items():
                        if (strand_name == 'scaffold') and (typ in ['half', 'half_v', 'half_h']):
                            pass
                        else:
                            export["{}_{}_{}".format(
                                name, strand_name, typ)] = n_co

            elif name in ["staples_length", "helices_staples_pass", "n_staples_domains",
                          "long_domains", "staple_domain_melt_T", "stacks_length", "loops_length"]:

                stats = get_statistics(value, name)
                for stat_name, stat in stats.items():
                    export[stat_name] = stat

            elif name == "alpha_value":
                for temp, alpha_value in value.items():
                    export[f"{name}_{temp}"] = alpha_value

            elif name in ["2_long_domains", "1_long_domains", "0_long_domains", "co_rule_violation"]:
                export[name] = len(value)
            else:
                export[name] = value

        return export

    def get_alpha_value(self) -> Dict[int, float]:
        """ alpha value : The ratio of number of staples having doamins with melting
            temperature higher than critical temperature to the number of all staples in the structure]
        """
        def calculate(max_Ts, T_crit):
            """ calculates alpha value for a given critical temperature
            """
            maxT_staple_domains = list(max_Ts.values())
            n_domains_total = len(maxT_staple_domains)
            gen_domains_critical = (
                1 for T in maxT_staple_domains if T >= T_crit)
            return sum(gen_domains_critical) / n_domains_total

        T_crit = {40: int, 55: int, 70: int}

        alpha_values = {T: calculate(self.design.max_staple_melt_t, T)
                        for T in T_crit}
        return alpha_values

    def get_staples_length(self) -> List[int]:
        """creates a list of the length of each staple"""
        return [len(staple.tour) for staple in self.design.staples]

    def get_num_staple_helix(self) -> List[int]:
        # TODO cleanup
        """ creates a dictionary with a set of helices that it passes through for each staple"""
        staple_helix_dict = {staple: {base.h for base in staple.tour}
                             for staple in self.design.staples}
        return [len(helix_ids) for helix_ids in staple_helix_dict.values()]

    def get_number_nicks(self) -> int:
        return len(self.design.nicks)

    def get_n_helices(self, min_bases=20) -> int:
        helices = [h for h in self.design.helices.values() if (len(
            h.scaffold_bases) + len(h.scaffold_bases)) > min_bases]
        return len(helices)

    def get_dimension(self, min_bases=20) -> Tuple[int, int, int]:
        """gets the dimension of the structure as columns, rows, and position (min-max)
            NOTE: these values are wrong for multidomain designs
        """
        helices = [h for h in self.design.helices.values() if (len(
            h.scaffold_bases) + len(h.scaffold_bases)) > min_bases]

        min_cols = min(helices, key=lambda h: h.lattice_col)
        max_cols = max(helices, key=lambda h: h.lattice_col)
        columns = max_cols.lattice_col - min_cols.lattice_col + 1

        min_rows = min(helices, key=lambda h: h.lattice_row)
        max_rows = max(helices, key=lambda h: h.lattice_row)
        rows = max_rows.lattice_row - min_rows.lattice_row + 1

        min_base = min(helices, key=lambda h: h.min_scaffold_pos)
        max_base = max(helices, key=lambda h: h.max_scaffold_pos)
        bases = max_base.max_scaffold_pos - min_base.min_scaffold_pos + 1

        return (columns, rows, bases)

    def get_insertion_deletion_density(self):
        scaffolds_length = sum(len(s.tour) for s in self.design.scaffolds)
        ins_density = len(self.design.hps_insertions) / scaffolds_length
        del_density = len(self.design.hps_deletions) / scaffolds_length
        return ins_density, del_density

    def get_stacks_lengths(self):
        return [(len(s) - 1) for stacks in self.design.stacks.values() for s in stacks]

    def get_n_domains(self):
        # TODO: cleaunp
        domain_data = self.design.domain_data
        n_staples_domains = {staple: len(
            domain_data[staple]) for staple in domain_data.keys()}
        return list(n_staples_domains.values())

    def get_n_staples_domains(self):
        # TODO: cleaunp (with previous)
        domain_data = self.design.domain_data
        n_staples_domains = {staple: len(
            domain_data[staple]) for staple in domain_data.keys()}
        return list(n_staples_domains.values())

    def get_staples_with_long_domains(self) -> Dict[Strand, int]:

        # TODO: cleanup and move to compute
        """ computes numbers of long_domains for the each staple
                long domain are domains with 14 or more bases
        """
        domain_data = self.design.domain_data
        lengths_data = {staple: [len(domain.base_list) for domain in domain]
                        for staple, domain in domain_data.items()}
        return {staple: len(list(filter(lambda x: x >= 14, domain_length)))
                for staple, domain_length in lengths_data.items()}

    def divide_domain_lengths(self) -> dict:

        # TODO: cleanup and move to compute
        """[divide staples having 0, 1, 2 or more long domains ]

        Returns:
            dict: [
                2_long_domains: two or more long domains,
                1_long_domains: having only one long domain,
                0_long_domains: having no long domain,
                co_rule_violation: unpaired domains with less than 5 bases]
        """
        # TODO: cleanup
        data = dict()
        domain_unpaired = list()
        domain_data = self.design.domain_data
        long_domains = self.get_staples_with_long_domains()

        for staple, n_longs in long_domains.items():
            if n_longs >= 2:
                data.setdefault("2_long_domains", []).append(staple)
            elif n_longs == 1:
                data.setdefault("1_long_domains", []).append(staple)
            elif n_longs == 0:
                data.setdefault("0_long_domains", []).append(staple)

        for staple, domains in domain_data.items():
            domain_unpaired.extend(
                [domain for domain in domains if domain.base_list[0].across is None])

            data.setdefault("co_rule_violation", []).extend(
                [domain for domain in domains if (domain not in domain_unpaired) and (len(domain.base_list) < 5)])
        return data

    def get_co_density(self):
        """ calculate crossover density (number of crossovers is the structure divided by possible crossovers CadNano)
            NOTE: the values for number of possible_co and co_desity are not exact but close to the true value
        """
        # TODO: possible_ co numbers are not exactly correct
        # TODO: what does it do?

        def is_ds(pos, hid):
            # NOTE: called super often
            is_sc = (hid, pos, True) in self.design.hps_deletions
            is_sc = (hid, pos, True) in self.design.hps_deletions
            is_st = (hid, pos, False) in self.design.hps_deletions
            # (hid, pos) in self.deletions (note: list of (h,p) for all deletions)
            is_deletion = False

            return ((is_sc or is_st) or is_deletion)

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
        helices = self.design.helices.values()

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
                    # TODO: devision by two is assumed for counting each possible_co
                    # two times for a helix and its neighbour

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

    def classify_crossovers(self):
        data = {"scaffold": dict(), "staple": dict()}
        types = {"full": [co for co in self.design.crossovers if co.typ == "full"],
                 "half": [co for co in self.design.crossovers if co.typ == "half"],
                 "end": [co for co in self.design.crossovers if co.typ == "end"],
                 }

        for typ, crossovers in types.items():
            co_subsets = {"scaffold": {"": list(), "_h": list(), "_v": list()},
                          "staple": {"": list(), "_h": list(), "_v": list()}}

            for co in crossovers:
                strand = "scaffold" if co.is_scaffold == 'scaffold' else "staple"
                co_subsets[strand][""].append(co)
                if co.orientation == "horizontal":
                    co_subsets[strand]["_h"].append(co)
                else:
                    co_subsets[strand]["_v"].append(co)

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

    def get_blunt_ends(self):
        blunt_ends = set()
        first_bases = {staple.tour[0] for staple in self.design.staples}
        last_bases = {staple.tour[-1] for staple in self.design.staples}

        for end in (co for co in self.design.crossovers if co.typ == "end"):
            has_across = (end.connection1.base1.across is True) and (
                end.connection1.base2.across is True)
            if end.is_scaffold == 'scaffold' and has_across:
                blunt_ends = {end.connection1 for base in end.connection1 if (
                    base.across in first_bases) or (base.across in last_bases)}

        return blunt_ends

    def get_loops(self):
        loops = list()
        real_crossover = (
            co for co in self.design.crossovers if co.typ != "end")
        for co in real_crossover:
            sub = np.inf
            if co.is_scaffold == 'scaffold':
                stacks = tuple([
                    tuple([co.connection1.base1, co.connection2.base1]),
                    tuple([co.connection1.base2, co.connection2.base2])
                ])
                for stack in stacks:
                    if stack[0].across is None or stack[1].across is None:
                        continue
                    same_staple = (stack[0].across.strand
                                   == stack[1].across.strand)
                    sc = self.design.strands[stack[0].strand]
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
                for connection in [co.connection1, co.connection2]:
                    if connection is None:  # NOTE: only 1 for half_co
                        continue
                    if connection.base1.across is None or connection.base2.across is None:
                        continue
                    sc = self.design.strands[connection.base1.across.strand]
                    same_scaffold = (sc.id == connection.base2.across.strand)
                    if same_scaffold:
                        base_1 = sc.tour.index(connection.base1.across)
                        base_2 = sc.tour.index(connection.base2.across)
                        sub_new = abs(base_1 - base_2)

                        if sub_new > len(sc.tour) / 2:
                            sub_new = len(sc.tour) - sub_new
                        if sub_new < sub:
                            sub = sub_new
            if not np.isinf(sub):
                loops.append(sub)

        return loops

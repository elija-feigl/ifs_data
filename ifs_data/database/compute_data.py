import attr
import logging
import numpy as np

from typing import List, Dict, Tuple, Any
from operator import attrgetter

from ..core.utils import get_statistics, n_bases, save_division
from ..core.designData import Design

# TODO: classmethod etc
STRAND_TYPES = ["scaffold", "staple"]
CO_TYPES = ["full", "half", "end"]
ORIENTS = ["v", "h"]


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
        data['alpha_value'] = self.get_alpha_values()
        data["n_helices"] = self.get_n_helices()

        data["n_nicks"] = self.get_number_nicks()
        data["n_stacks"] = len(self.get_stacks_lengths())
        data["stacks_length"] = self.get_stacks_lengths()
        data["loops_length"] = self.get_loops()

        data["n_insertions"] = len(design.hps_insertions) // 2
        data["n_deletions"] = len(design.hps_deletions) // 2
        data["insertions_denisty"], data["deletion_denisty"] = self.get_insertion_deletion_density()

        data["n_bluntends"] = len(self.get_blunt_ends())

        # staple stats
        data["n_staples"] = len(design.staples)
        data["staples_length"] = self.get_staples_length()
        data["helices_staples_pass"] = self.get_num_staple_helix()

        # domains
        data["n_staples_domains"] = self.get_n_staples_domains()
        data.update(self.get_long_domains())
        data["co_rule_violation"] = self.get_co_rule_violation()
        data["staple_domain_melt_T"] = self.get_staple_domain_melt_T()

        # crossovers
        data.update(self.get_full_scaff_co_typ_stat())
        data["co_set"], data["co_possible"], data["co_density"] = self.get_co_density()

        self.data = data

    def prep_data_for_export(self) -> dict:
        export = dict()
        for name, value in self.data.items():
            if name in ["co_set", "co_possible", "co_density"]:
                for strand_name, subtypes in value.items():
                    for typ, n_co in subtypes.items():
                        if (strand_name == 'scaffold') and (typ in ['half', 'half_v', 'half_h']):
                            self.logger.debug(
                                "statistics of design reveal scaffold half crossover.")
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

            elif name in ["co_rule_violation"]:
                export[name] = len(value)
            else:
                export[name] = value

        return export

    def get_alpha_values(self, T_crit={40: int, 55: int, 70: int}) -> Dict[int, float]:
        """ alpha value : The ratio of number of staples having doamins with melting
            temperature higher than critical temperature to the number of all staples in the structure]
        """
        def alpha(max_Ts, T_thres):
            """ calculates alpha value for a given critical temperature."""
            domains_critical = sum(1 for T in max_Ts if T >= T_thres)
            return domains_critical / len(max_Ts)

        domain_max_T = list(self.design.max_staple_melt_t.values())
        alpha_values = {T: alpha(domain_max_T, T) for T in T_crit}
        return alpha_values

    def get_staples_length(self) -> List[int]:
        """ creates a list of the length of each staple"""
        return [len(staple.tour) for staple in self.design.staples]

    def get_num_staple_helix(self) -> List[int]:
        """ creates a dictionary with a set of helices that it passes through for each staple"""
        return [len(staple.helix_list) for staple in self.design.staples]

    def get_number_nicks(self) -> int:
        return len(self.design.nicks)

    def get_n_helices(self, min_bases=20) -> int:
        return len([h for h in self.design.helices if n_bases(h) > min_bases])

    def get_dimension(self, min_bases=20) -> Tuple[int, int, int]:
        """ gets the dimension of the structure as columns, rows, and position (min-max)
            NOTE: these values are wrong for multidomain designs
        """
        def dim_attr(helices, attr: str):
            min_attr = min(helices, key=attrgetter(attr))
            max_attr = max(helices, key=attrgetter(attr))
            return getattr(max_attr, attr) - getattr(min_attr, attr) + 1

        helices = [h for h in self.design.helices if n_bases(h) > min_bases]
        return (dim_attr(helices, "lattice_col"),
                dim_attr(helices, "lattice_row"),
                dim_attr(helices, "min_scaffold_pos"))

    def get_insertion_deletion_density(self):
        scaffolds_length = sum(len(s.tour) for s in self.design.scaffolds)
        ins_density = len(self.design.hps_insertions) / scaffolds_length
        del_density = len(self.design.hps_deletions) / scaffolds_length
        return ins_density, del_density

    def get_stacks_lengths(self):
        """ creates a list of the length of each stack."""
        return [(len(s) - 1) for stacks in self.design.stacks.values() for s in stacks]

    def get_n_staples_domains(self):
        """ creates a list of the length of each domain of all staples."""
        return [len(staple.domain_list) for staple in self.design.staples]

    def get_long_domains(self) -> dict:
        """ list all long domains per staple
            count how many staples have either 0,1 or 2 long segments."""
        def _is_long_(domain):
            return True if (len(domain.base_list) >= 14) else False

        long_domains_staple = list()
        for staple in self.design.staples:
            long_domains = [d.id for d in staple.domain_list if _is_long_(d)]
            long_domains_staple.append(len(long_domains))

        n_long_types = [0, 0, 0]
        for n_longs in long_domains_staple:
            n_typ = 2 if n_longs > 2 else n_longs
            n_long_types[n_typ] += 1

        return {
            "long_domains": long_domains_staple,
            '0_long_domains': n_long_types[0],
            '1_long_domains': n_long_types[1],
            '2_long_domains': n_long_types[2],
        }

    def get_co_rule_violation(self) -> List:
        """ co_rule_violation: unpaired domains with less than min_length base."""
        def _violates_co_rules(domain, min_length=5):
            is_short = True if (len(domain.base_list) < min_length) else False
            is_ds = (domain.base_list[0].across is None)
            return is_short and is_ds

        return [d for staple in self.design.staples
                for d in staple.domain_list if _violates_co_rules(d)]

    def get_staple_domain_melt_T(self) -> List:
        return list(self.design.max_staple_melt_t.values())

    def get_co_density(self):
        """ calculate crossover density (number of crossovers is the structure
            divided by possible crossovers CadNano)."""
        def init_co_dict(seperate_types: bool = True):
            types = CO_TYPES if seperate_types else[]
            co_dict: Dict[str, Dict[str, int]] = dict()
            for strand_type in STRAND_TYPES:
                for typ in types:
                    for orient in ORIENTS:
                        keys = [f"co_{orient}", "co"]
                        if seperate_types:
                            keys += [f"{typ}_{orient}", f"{typ}"]
                        for key in keys:
                            co_dict[strand_type][key] = 0
            return co_dict

        def count_co(position_list):
            # NOTE: innaccurate for helices that are not continuous
            ends = 0
            n = len(position_list)
            if n == 0:
                return 0, 0
            if n == 1:
                return 1, 0
            if n == 2:
                if position_list[0] + 1 != position_list[1]:
                    return 1, 0
                else:
                    return 0, 1

            if position_list[0] + 1 != position_list[1]:
                position_list.pop(0)
                ends += 1
            if position_list[-1] - 1 != position_list[-2]:
                position_list.pop()
                ends += 1
            # NOTE: divison by 2 as two connections per crossover (not 100% accurate)
            return ends, len(position_list) // 2

        def is_possible(typ, p, helix):
            """ A base in the neighbouring helix to connect to exists."""
            bases = getattr(helix, f"{typ}_bases")
            return any(b.p == p for b in bases)

        def is_oriented(typ, helix, helix2):
            """ Horizontal connections are in the same row, verticals in same column"""
            orient = "lattice_row" if typ == 'h' else "lattice_col"
            return getattr(helix, orient) == getattr(helix2, orient)

        possible_connections = init_co_dict(seperate_types=False)
        for helix in self.design.helices:
            for strand_type in STRAND_TYPES:
                for orient in ORIENTS:
                    connections = getattr(
                        helix, f"possible_{strand_type}_crossovers")

                    positions = set()
                    for helix2, p, _ in connections:
                        is_valid = (is_oriented(orient, helix, helix2)
                                    and is_possible(strand_type, p, helix2))
                        if is_valid:
                            positions.add(p)

                    n_end, n_co = count_co(sorted(list(positions)))

                    # NOTE: 2021-03-19: ends will be excluded from co density estimates
                    n_connections = n_co  # + n_end
                    possible_connections[strand_type]["co"] += n_connections
                    possible_connections[strand_type][f"co_{orient}"] += n_connections

        # NOTE: divison by 2 for counting each connection twice (for each helix)
        # TODO: is there prettier way?
        for strand_type in STRAND_TYPES:
            possible_connections[strand_type]["co"] //= 2
            for orient in ORIENTS:
                possible_connections[strand_type][f"co_{orient}"] //= 2

        # part2 get actual crossovers
        set_crossovers = init_co_dict()
        for co in self.design.crossovers:
            strand_type = "scaffold" if co.is_scaffold else "staple"
            orient = "h" if co.orientation == "horizontal" else "v"
            set_crossovers[strand_type][f"{co.typ}_{orient}"] += 1
            set_crossovers[strand_type][f"{co.typ}"] += 1
            # NOTE: 2021-03-19: ends will be excluded from co density estimates
            if co.typ != "end":
                set_crossovers[strand_type][f"co_{orient}"] += 1
                set_crossovers[strand_type]["co"] += 1

        # part3 get crossover denisity
        co_density = dict()
        for strand_type in STRAND_TYPES:
            co_density[strand_type] = dict()
            for key, n_possible in possible_connections[strand_type].items():
                if "co" in key:
                    n_set = set_crossovers[strand_type][key]
                    co_density[strand_type][key] = save_division(
                        n_set, n_possible)

        return set_crossovers, possible_connections, co_density

    def get_full_scaff_co_typ_stat(self):
        type_count = [0, 0, 0]
        for full in (co for co in self.design.crossovers if co.typ == "full" and co.is_scaffold):
            type_count[full.scaff_full_type-1] += 1

        return {
            'full_scaf_co_type_1': type_count[0],
            'full_scaf_co_type_2': type_count[1],
            'full_scaf_co_type_3': type_count[2],
        }

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

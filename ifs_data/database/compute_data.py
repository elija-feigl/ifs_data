import logging
from dataclasses import dataclass
from itertools import chain
from operator import attrgetter
from typing import Any, Dict, List, Tuple

from nanodesign.data.base import DnaBase as Base
from nanodesign.data.strand import DnaStrand as Strand

from ..core.designData import Design
from ..core.utils import get_statistics, save_division
from ..data.crossover import Crossover

STRAND_TYPES = ["scaffold", "staple"]
CO_TYPES = ["full", "half", "end"]
ORIENTS = ["v", "h"]


@dataclass
class DesignStats(object):
    design: Design

    def __post_init__(self):
        self.logger = logging.getLogger(__name__)

    def compute_data(self) -> None:
        design = self.design
        data: Dict[str, Any] = dict()
        data["name"] = design.name
        data["lattice_type"] = design.lattice
        data.update(self.get_dimension())
        data["n_helices"] = self.get_n_helices()

        data["n_nicks"] = self.get_number_nicks()
        stacks = self.get_stacks_lengths()
        data["n_stacks"] = len(stacks)
        data["stacks_length"] = stacks
        data["loops_length"] = self.get_loops()

        # NOTE: skips and insertions are only counted once per bp
        data["n_insertions"] = len(design.hps_insertions) // 2
        data["n_deletions"] = len(design.hps_deletions) // 2
        data.update(self.get_insertion_deletion_density())

        # staple stats
        data["n_staples"] = len(design.staples)
        data["staples_length"] = self.get_staples_length()
        data["helices_staples_pass"] = self.get_num_staple_helix()
        data['alpha_value'] = self.get_alpha_values()

        # domains
        data["n_staples_domains"] = self.get_n_staples_domains()
        data.update(self.get_long_domains())
        data["co_rule_violation"] = self.get_co_rule_violation()
        data["staple_domain_melt_T"] = self.get_staple_domain_melt_T()

        # crossovers
        data.update(self.get_full_scaff_co_typ_stat())
        data.update(self.get_co_density())
        data["n_bluntends"] = self.get_n_blunt_ends()

        self.data = data

    def prep_data_for_export(self) -> dict:
        export = dict()
        for name, value in self.data.items():
            if name in ["co_set", "co_possible", "co_density"]:
                for strand_name, subtypes in value.items():
                    for typ, n_co in subtypes.items():
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
        """ creates a list of the length of each staple."""
        return [len(staple.tour) for staple in self.design.staples]

    def get_num_staple_helix(self) -> List[int]:
        """ creates a list with a set of helices that it passes through for each staple"""
        return [len(staple.helix_list) for staple in self.design.staples]

    def get_number_nicks(self) -> int:
        return len(self.design.nicks)

    def get_n_helices(self, min_bases=20) -> int:
        return len([h for h in self.design.helices if len(h.scaffold_bases) > min_bases])

    def get_dimension(self, min_bases=20) -> Dict[str, int]:
        """ gets the dimension of the structure as columns, rows, and position (min-max)
            NOTE: these values are wrong for multidomain designs
        """
        def dim_attr(helices, attr: str):
            min_attr = min(helices, key=attrgetter(attr))
            max_attr = max(helices, key=attrgetter(attr))
            return getattr(max_attr, attr) - getattr(min_attr, attr) + 1

        helices = [h for h in self.design.helices if len(
            h.scaffold_bases) > min_bases]
        return {
            "dim_x": dim_attr(helices, "lattice_col"),
            "dim_y": dim_attr(helices, "lattice_row"),
            "dim_z": dim_attr(helices, "min_scaffold_pos"),
        }

    def get_insertion_deletion_density(self):
        # NOTE: skips and insertions are only counted once per bp
        scaffolds_length = sum(len(s.tour) for s in self.design.scaffolds)
        ins_density = len(self.design.hps_insertions) / 2 / scaffolds_length
        del_density = len(self.design.hps_deletions) / 2 / scaffolds_length
        return {"insertions_denisty": ins_density, "deletion_denisty": del_density}

    def get_stacks_lengths(self):
        """ creates a list of the length of each stack."""
        return [(len(s) - 1) for stacks in self.design.stacks.values() for s in stacks]

    def get_n_staples_domains(self):
        """ creates a list of the length of each domain of all staples."""
        return [len(staple.domain_list) for staple in self.design.staples]

    def get_long_domains(self) -> dict:
        """ list all long domains per staple
            count how many staples have either 0, 1 or 2 long segments."""
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

    def get_co_density(self) -> dict:
        """ calculate crossover density (number of crossovers is the structure
            divided by possible crossovers CadNano)."""
        def init_co_dict(seperate_types: bool = True):
            types = CO_TYPES if seperate_types else []
            co_dict: Dict[str, Dict[str, int]] = dict()
            for strand_type in STRAND_TYPES:
                co_dict[strand_type] = dict()
                for orient in ORIENTS:
                    keys = [f"co_{orient}", "co"]
                    if seperate_types:
                        for typ in types:
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
            """ Horizontal connections are in the same row, verticals in same column."""
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
        co_density: dict = dict()
        for strand_type in STRAND_TYPES:
            co_density[strand_type] = dict()
            for key, n_possible in possible_connections[strand_type].items():
                if "co" in key:
                    n_set = set_crossovers[strand_type][key]
                    co_density[strand_type][key] = save_division(
                        n_set, n_possible)

        return {
            "co_set": set_crossovers,
            'co_possible': possible_connections,
            'co_density': co_density,
        }

    def get_full_scaff_co_typ_stat(self) -> Dict[str, int]:
        type_count = [0, 0, 0]
        for full in (co for co in self.design.crossovers if co.typ == "full" and co.is_scaffold):
            type_count[full.scaff_full_type-1] += 1

        return {
            'full_scaf_co_type_1': type_count[0],
            'full_scaf_co_type_2': type_count[1],
            'full_scaf_co_type_3': type_count[2],
        }

    def get_n_blunt_ends(self) -> int:
        """ Blunt ends are scaffold end crossovers where the staple ends at the crossover."""
        def is_scaff_end(co: Crossover) -> bool:
            return co.typ == "end" and co.is_scaffold

        def termini(staple: Strand) -> Tuple[Base, Base]:
            return (staple.tour[0], staple.tour[-1])

        def is_blunt(co: Crossover) -> bool:
            return all((b.across in staple_termini) for b in co.bases)

        ends = [co for co in self.design.crossovers if is_scaff_end(co)]

        staple_termini_tuple = (termini(s) for s in self.design.staples)
        staple_termini = set(chain.from_iterable(staple_termini_tuple))
        return sum(1 for co in ends if is_blunt(co))

    def get_loops(self):
        """ A scaffold loop is the distance along the scaffold of two scaffold bases
                conected via a crossover. There are two types:
                 * connected by staple crossover
                 * connected via staple at full scaffold crossover
        """
        def correct_diff(co, diff):
            if co.is_scaffold:
                sc_id = list(co.bases)[0].strand
            else:
                sc_id = list(co.bases)[0].across.strand

            sc = next((s for s in self.design.scaffolds if s.id == sc_id), None)
            if sc is None:
                import ipdb
                ipdb.set_trace()
            sc_length = len(sc.tour)

            return diff if diff > (sc_length / 2) else (sc_length - diff)

        loops = list()
        # = [co for co in self.design.crossovers if co.typ != "end"]
        for co in (co for co in self.design.crossovers if co.typ != "end"):
            if co.is_scaffold:
                if co.typ == "half":
                    self.logger.debug(f"Ignoring scaffold half {co}.")
                    continue
                p_diff = abs(co.connection1.scaffold_pos[0] -
                             co.connection2.scaffold_pos[0])

            else:
                if not any(b.across.is_scaf for b in co.bases):
                    self.logger.debug(f"Ignoring staple-staple {co}.")
                    continue
                sc_pos = co.connection1.scaffold_pos
                p_diff = sc_pos[1] - sc_pos[0]

            distance = correct_diff(co, p_diff)
            loops.append(distance)
        return loops

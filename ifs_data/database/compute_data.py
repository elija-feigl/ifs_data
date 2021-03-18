import attr
import logging
from typing import List, Dict, Tuple, Set, Any

from ..core.utils import get_statistics, get_full_scaff_co_typ_stat
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
        data["n_helices"] = len(design.helices.keys())
        data["n_deletions"] = len(design.hps_deletions) // 2
        data["n_insertions"] = len(design.hps_insertions) // 2
        data["insertions_denisty"], data["deletion_denisty"] = self.get_insertion_deletion_density()

        data["n_nicks"] = self.get_number_nicks()
        data["n_stacks"] = len(self.get_stacks_lengths())
        data["stacks_length"] = self.get_stacks_lengths()
        data["loops_length"] = design.get_loops()

        data["n_bluntends"] = len(design.get_blunt_ends())

        # staple stats
        data["n_staples"] = len(design.staples)
        data["staples_length"] = self.get_staples_length()
        data["helices_staples_pass"] = self.get_num_staple_helix()

        # domains
        data["n_staples_domains"] = self.get_n_staples_domains()
        data["long_domains"] = list(self.design.get_staples_with_long_domains().values())
        data.update(design.divide_domain_lengths())
        data["staple_domain_melt_T"] = list(
            design.max_staple_melt_t.values())

        # crossovers
        data["co_set"] = design.classify_crossovers()
        data.update(get_full_scaff_co_typ_stat(design))
        data["co_possible"], data["co_density"] = design.get_co_density()

        self.data = data

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

    def get_dimension(self) -> Tuple[int, int, int]:
        """gets the dimension of the structure as columns, rows, and position (min-max)
            NOTE: these values are wrong for multidomain designs
        """
        helices = self.design.helices.values()

        min_cols = min(helices, key=lambda h: h.lattice_col)
        max_cols = max(helices, key=lambda h: h.lattice_col)
        columns = max_cols.lattice_col - min_cols.lattice_col + 1

        min_rows = min(helices, key=lambda h: h.lattice_row)
        max_rows = max(helices, key=lambda h: h.lattice_row)
        rows = max_rows.lattice_row - min_rows.lattice_row + 1

        min_base = min(helices, key=lambda h: h.min_scaffold_pos)
        max_base = max(helices, key=lambda h: h.max_scaffold_pos)
        bases = max_base.min_scaffold_pos - min_base.min_scaffold_pos + 1

        return (columns, rows, bases)

    def get_insertion_deletion_density(self):
        scaffolds_length = sum(len(s.tour) for s in self.design.scaffolds)
        ins_density = len(self.design.hps_insertions) / scaffolds_length
        del_density = len(self.design.hps_deletions) / scaffolds_length
        return ins_density, del_density

    def get_stacks_lengths(self):
        return [(len(s) - 1) for stacks in self.design.stacks.values() for s in stacks]

    def get_n_domains(self):
        #TODO: cleaunp
        domain_data = self.design.domain_data
        n_staples_domains = {staple: len(
            domain_data[staple]) for staple in domain_data.keys()}
        return list(n_staples_domains.values())

    def get_n_staples_domains(self):
        #TODO: cleaunp (with previous)
        domain_data = self.design.domain_data
        n_staples_domains = {staple: len(
            domain_data[staple]) for staple in domain_data.keys()}
        return list(n_staples_domains.values())

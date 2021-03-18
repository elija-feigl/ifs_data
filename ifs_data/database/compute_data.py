import attr
import logging

from ..core.utils import get_statistics, get_full_scaff_co_typ_stat
from ..core.designData import Design


@attr.s
class DesignStats(object):
    design: Design = attr.ib()

    logger = logging.getLogger(__name__)

    def compute_data(self) -> None:
        design = self.design
        data = {}
        data["name"] = design.name
        data["lattice_type"] = design.lattice
        data["dim_x"], data["dim_y"], data["dim_z"] = design.get_dimension()
        data['alpha_value'] = design.get_alpha_value()
        data["n_helices"] = len(design.helices.keys())
        data["n_deletions"] = len(design.hps_deletions)
        data["n_insertions"] = len(design.hps_insertions)
        data.update(design.get_insertion_deletion_density())

        data["n_nicks"] = len(design.get_nicks())
        data["n_stacks"] = len(design.get_stacks_lengths())
        data["stacks_length"] = design.get_stacks_lengths()
        data["loops_length"] = design.get_loops()

        data["n_bluntends"] = len(design.get_blunt_ends())

        # staple stats
        data["n_staples"] = len(design.staples)
        data["staples_length"] = design.get_staples_length()
        data["helices_staples_pass"] = design.get_num_staple_helix()

        # domains
        data["n_staples_domains"] = list(design.n_staples_domains.values())
        data["long_domains"] = list(design.long_domains.values())
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

        

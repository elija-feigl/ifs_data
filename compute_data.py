from designData import DesignData
import argparse
from pathlib import Path
from utils import get_statistics


class Compute(object):

    def compute_data(designdata) -> None:
        data = {}
        data["name"] = designdata.name
        data["lattice_type"] = designdata.get_lattice_type()
        data["dim_x"], data["dim_y"], data["dim_z"] = designdata.get_dimension()
        data['alpha_value'] = designdata.alpha_value
        data["n_helices"] = len(designdata.dna_structure.structure_helices_map)
        data["n_skips"] = len(designdata.dna_structure.Dhp_skips)
        data["n_nicks"] = len(designdata.nicks)
        data["n_stacks"] = len(designdata.stacks)
        data["stacks_length"] = designdata.stacks_lengths
        data["loops_length"] = designdata.loops_length_list
        data.update(designdata.get_insertion_deletion_density())
        data["n_bluntends"] = len(designdata.get_blunt_ends())

        # staple stats
        data["n_staples"] = len(designdata.staples)
        data["staples_length"] = designdata.get_staples_length()
        data["helices_staples_pass"] = list(
            designdata.num_staple_helix_dict.values())

        # domains
        data["n_staples_domains"] = list(designdata.n_staples_domains.values())
        data["long_domains"] = list(designdata.long_domains.values())
        data.update(designdata.divide_domain_lengths())
        data["staple_domain_melt_T"] = list(
            designdata.max_staple_melt_t.values())

        # crossovers
        data["co_set"] = designdata.classify_crossovers()
        data.update(designdata.get_all_full_scaff_crossover_types())
        data["co_possible"], data["co_density"] = designdata.get_co_density()

        designdata.data = data

    def prep_data_for_export(designdata) -> dict:
        export = dict()
        for name, value in designdata.data.items():
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
    designdata = DesignData(json=json, name=json.name, seq='8064')
    Compute.compute_data(designdata)
    data = Compute.prep_data_for_export(designdata)
    with open(outname, mode="w+") as outfile:
        header = ",".join(str(k) for k in data.keys())
        outfile.write(header + "\n")
        export = ",".join(str(v) for v in data.values())
        outfile.write(export + "\n")


if __name__ == "__main__":
    main()

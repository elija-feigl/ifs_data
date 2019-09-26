#!/usr/bin/env python
from nanodesign.converters import Converter
import numpy as np
import ipdb  # use this for debugging instead of print() (ipdb.set_trace())

class DesignData(object):

    def __init__(self, name: str):
        self.name: str = name
        self.dna_structure = self.init_design()
        self.all_strands: list = self.dna_structure.strands
        self.pairedbases: list = self.init_pairedbases()
        self.data: dict = {}
        self.all_bases: list = []
        for strand in self.all_strands:
            for base in strand.tour:
                self.all_bases.append(base)

    def get_base_hps(self, h, p, is_scaffold):
        base=None 
        # TODO: not yet implemented
        return base
    
    def compute_data(self) -> dict:
        data = {}
        data["total_co"] = self.get_total_co()
        data["half_co"] = self.get_half_co()
        data["end_loop"] = self.get_end_loop()
        data["full_co"] = self.get_full_co()
        data["staples"] = self.get_number_of_staple()
        scaf_co, staple_co = self.get_staple_scaffold_co()
        data["scaf_co"] = scaf_co
        data["staple_co"] = staple_co
        staple_length_avg, staple_length_std = self.get_staple_length_parameter()
        data["staple_length_avg"] = staple_length_avg
        data["staple_length_std"] = staple_length_std
        data["num_skips"] = self.get_number_skips()

        self.data = data
        return self.data

    def init_design(self):
        file_name = self.name + ".json"
        seq_file = self.name + ".seq"
        seq_name = None
        converter = Converter()
        converter.read_cadnano_file(file_name, seq_file, seq_name)
        converter.dna_structure.get_domains()
        return converter.dna_structure

    def init_pairedbases(self) -> list:
        pairedbases = []
        for strand in self.all_strands:
            for base in strand.tour:
                if base.across is not None:
                    pairedbases.append(base)

        return pairedbases
    
    def get_all_bases(self) -> list:
        all_bases = []
        for strand in self.all_strands:
            for base in strand.tour:
              all_bases.append(base)
        return all_bases
        
    def get_staple_length_parameter(self):
        len_strands = []
        for strand in self.all_strands:
            if not strand.is_scaffold:
                tour_clean = [base for base in strand.tour if base.num_deletions == 0]
                len_strands.append(len(tour_clean))
        return np.average(len_strands), np.std(len_strands)

    def get_number_skips(self) -> int:
        num_skips = 0
        for strand in self.all_strands:
            skips = [base for base in strand.tour if base.num_deletions == -1]
            num_skips += len(skips)
        return num_skips

    def get_number_of_staple(self) -> int:
        num_staples = []
        for strand in self.all_strands:
            if not strand.is_scaffold:
                num_staples.append(strand)
        return len(num_staples)

    def get_helix_per_staple(self) -> list:
        helices = set()
        helix_list = []
        for strand in self.all_strands:
            if not strand.is_scaffold:
                for base in strand.tour:
                    if base.h not in helices:
                        helices.add(base.h)
                helix_list.append(len(helices))
        return helix_list

    def init_helix_dict(self) -> dict:
        helix_dic = {}
        helix_list = self.get_helix_per_staple()
        for strand, helices in zip(self.all_strands, helix_list):
            dic = {str(strand): helices}
            helix_dic.update(dic)
        return helix_dic

    def get_staple_scaffold_co(self):
        N_scaf_co = 0
        N_staple_co = 0

        for strand in self.all_strands:
                for base in strand.tour:
                    if self.dna_structure._check_base_crossover(base):
                        if base.is_scaf:
                            N_scaf_co += 1
                        else:
                            N_staple_co += 1
        return N_scaf_co/2, N_staple_co/2

    def get_total_co(self) -> int:
        number = 0
        for strand in self.all_strands:
                for base in strand.tour:
                    if self.dna_structure._check_base_crossover(base):
                        number += 1
        number = number/2
        return number

    def get_end_loop(self) -> int:
        number = 0
        for strand in self.all_strands:
            for base in strand.tour:
                if base.is_scaf:
                    if base.h + 1 not in self.all_bases or base.h - 1 not in self.all_bases:
                        number += 1
        return number

    def get_full_co(self) -> int:
        number = 0
        for base in self.pairedbases:
            if self.dna_structure._check_base_crossover(base) and self.dna_structure._check_base_crossover(base.across):
                number += 1
        number = number/2
        return number

    def get_half_co(self) -> int:
        number = 0
        for base in self.pairedbases:
            if not base.is_scaf:
                if self.dna_structure._check_base_crossover(base) and not self.dna_structure._check_base_crossover(base.across):
                    number += 1

        number = number/2
        return number

def tester(all_strands, dna_structure) -> None:
    for strand in all_strands:
        if not strand.is_scaffold:
            for base in strand.tour:
                if base.h == 0 and base.p == 71:
                    sin_base = base
                    if not dna_structure._check_base_crossover(sin_base.across):
                        print('bingo')
                    else:
                        print('no')
    return

def export_data(data: dict, name: str) -> None:
    header = ", ".join(data.keys())
    data_str = ", ".join([str(i) for i in data.values()])
    with open(name+"-designdata.csv", mode="w+") as out:
        out.write(header + "\n")
        out.write(data_str)
        out.write("\nEND")
    return

def main():
    print("master, I am awaiting the name of your design")
    name = input()
    print("Thank you Sir")
    designData = DesignData(name=name)
    data = designData.compute_data()
    export_data(data=data, name=name)
    return

if __name__ == "__main__":
    main()


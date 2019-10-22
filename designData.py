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
        self.all_bases: list = self.get_all_bases()
        self.hps_base: dict = self.init_hps()
        #self.position_of_base: dict = self.get_hps_from_base()
        self.all_co_skips: list = self.get_all_co_skips()
        self.all_co: list = self.get_all_co()
        self.scaf_co_bases, self.staple_co_bases = self.get_staple_scaffold_co_bases()
        self.helices_n: list = self.get_helix_per_staple()
        self.helix_dic: dict = self.init_helix_dict()
        self.first_bases, self.last_bases = self.get_first_last_bases_of_strands()
        self.nicks: list = self.get_nicks()
         
    def compute_data(self) -> dict:
        data = {}
        data["n_staples"] = self.get_number_of_staple()
        staple_length_avg, staple_length_std = self.get_staple_length_parameter()
        data["staple_length_avg"] = staple_length_avg
        data["staple_length_std"] = staple_length_std
        data["avg_helices_staples_pass"] = self.avg_helices_staples_pass()
        data["num_skips"] = self.get_n_skips()
        data["num_nicks"] = self.get_n_nicks()
        data["total_co"] = self.get_total_n_co()
        n_half, n_full, n_end = self.get_co_type()
        data["avg_n_staples_co"] = self.avg_staples_crossovers()
        data["half_co"] = n_half
        data["full_co"] = n_full
        data["endloops"] = n_end
        scaf_co, staple_co = self.get_staple_scaffold_co_n()
        data["scaf_co"] = scaf_co
        data["staple_co"] = staple_co
        n_staple_half, n_staple_full, n_staple_end = self.get_staple_co_types()
        data["half_staple_co"] = n_staple_half
        data["full_staple_co"] = n_staple_full
        data["staple_endloops"] = n_staple_end
        n_scaf_half, n_scaf_full, n_scaf_end = self.get_scaf_co_types()
        data["half_scaff_co"] = n_scaf_half
        data["full_scaff_co"] = n_scaf_full
        data["scaff_endloops"] = n_scaf_end

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
    
    def init_hps(self) -> dict:
        hps_base = {}
        for base in self.all_bases :
            position = (base.h, base.p, base.is_scaf)
            hps_base[position] = base
        return hps_base
            
    def get_base_from_hps(self, h, p, is_scaffold):
        try:
            base = self.hps_base[(h, p, is_scaffold)]
        except KeyError:
            return None
        if base.num_deletions != 0:
            return 'skip'
        else:
            return base
        
    def get_hps_from_base(self, base):
        return (base.h, base.p, base.is_scaf)

        
    def get_all_bases(self) -> list:
        all_bases = []
        for strand in self.all_strands:
            for base in strand.tour: 
                if base.num_deletions == 0:
                    all_bases.append(base)
        return all_bases
        
    def get_staple_length_parameter(self):
        len_strands = []
        for strand in self.all_strands:
            if not strand.is_scaffold:
                tour_clean = [base for base in strand.tour if base.num_deletions == 0]
                len_strands.append(len(tour_clean))
        return np.average(len_strands), np.std(len_strands)

    def get_n_skips(self) -> int:
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
        helices_n = []
        for strand in self.all_strands:
            helices = set()
            if not strand.is_scaffold:
                for base in strand.tour:
                    if self.dna_structure._check_base_crossover(base):
                        if base.h not in helices:
                            helices.add(base.h)
                helices_n.append(len(helices))
        #ipdb.set_trace()
        return helices_n
            
    def avg_helices_staples_pass(self) -> int:
        return np.average(self.helices_n)
    
    def init_helix_dict(self) -> dict:
        helix_dic = {}
        helix_list = self.get_helix_per_staple()
        for strand, helices in zip(self.all_strands, helix_list):
            dic = {str(strand): helices}
            helix_dic.update(dic)
        #ipdb.set_trace()
        return helix_dic
    
    def get_first_last_bases_of_strands(self) -> list:
        first_bases = []
        last_bases = []
        for strand in self.all_strands:
            if not strand.is_scaffold:
                first_bases.append(strand.tour[0])
                last_bases.append(strand.tour[-1])
        #ipdb.set_trace()
        return first_bases, last_bases
    
    def get_nicks(self) ->int:
        nicks = []
        for base in self.first_bases:
            base_plus = self.get_base_from_hps(base.h, base.p + 1, base.is_scaf)
            base_minus = self.get_base_from_hps(base.h, base.p - 1, base.is_scaf)
            if base_plus in self.last_bases:
                nicks.append((base,base_plus))#order of nick is always (first,last)
            elif base_minus in self.last_bases:
                nicks.append((base,base_minus))
        #ipdb.set_trace()       
        return nicks
    
    def get_n_nicks(self) -> int:
        return len(self.nicks)
    
    def get_all_co_skips(self) -> list:
        all_co_skips = []
        for strand in self.all_strands:
                for base in strand.tour:
                    if self.dna_structure._check_base_crossover(base):
                        if base.num_deletions == -1:
                            if base.up.h != base.h:
                                all_co_skips.append(base.down)
                            elif base.up.h == base.h:
                                all_co_skips.append(base.up)
                        else:
                            all_co_skips.append(base)
        #ipdb.set_trace()
        return all_co_skips
    
    def get_all_co(self):
        all_co = []
        for strand in self.all_strands:
                for base in strand.tour:
                    if self.dna_structure._check_base_crossover(base):
                            all_co.append(base)
        #ipdb.set_trace()
        return all_co
            
    def get_total_n_co(self) -> int:
        return len(self.all_co)/2.
    
    def avg_staples_crossovers(self) -> int:
        co_staple = []
        n_co_staples = []
        for strand in self.all_strands:
            if not strand.is_scaffold:
                for base in strand.tour:
                    if self.dna_structure._check_base_crossover(base):
                        co_staple.append(base)
                n_co_staples.append(len(co_staple)/2.)
                co_staple = []
        return np.average(n_co_staples)       
        
    def get_co_type(self) -> int:
        n_end = 0
        n_half = 0
        n_full = 0

        for co in self.all_co_skips:
            base_plus = self.get_base_from_hps(co.h, co.p + 1, co.is_scaf)
            base_minus = self.get_base_from_hps(co.h, co.p - 1, co.is_scaf)
            is_end = base_plus is None or base_minus is None
            is_full = base_plus in self.all_co_skips or base_minus in self.all_co_skips
            is_half = not is_end and not is_full
            if is_end:
                n_end += 1
            elif is_full:
                n_full += 1
            elif is_half:
                n_half += 1
            #ipdb.set_trace()
        return n_half/2., n_full/4., n_end/2.
 
    def get_staple_scaffold_co_bases(self):
        scaf_co_bases = []
        staple_co_bases = []
        
        for base in self.all_co_skips:
        #considering the skips in crossovers
            if base.is_scaf:
                scaf_co_bases.append(base)
            else:
                staple_co_bases.append(base)
                            
        return scaf_co_bases, staple_co_bases
    
    def get_staple_scaffold_co_n(self):
        return len(self.scaf_co_bases)/2., len(self.staple_co_bases)/2.

    def get_scaf_co_types(self) -> int:
        n_scaf_end = 0
        n_scaf_half = 0
        n_scaf_full = 0
        for co in self.scaf_co_bases:
            base_plus = self.get_base_from_hps(co.h, co.p + 1, co.is_scaf)
            base_minus = self.get_base_from_hps(co.h, co.p - 1, co.is_scaf)
            is_end = base_plus is None or base_minus is None
            is_full = base_plus in self.scaf_co_bases or base_minus in self.scaf_co_bases
            is_half = not is_end and not is_full
            if is_end:
                n_scaf_end += 1
            elif is_full:
                n_scaf_full += 1
            elif is_half:
                n_scaf_half += 1

        return n_scaf_half/2., n_scaf_full/4., n_scaf_end/2.
            
    def get_staple_co_types(self) -> int:
        n_staple_end = 0
        n_staple_half = 0
        n_staple_full = 0
        for co in self.staple_co_bases:
            base_plus = self.get_base_from_hps(co.h, co.p + 1, co.is_scaf)
            base_minus = self.get_base_from_hps(co.h, co.p - 1, co.is_scaf)
            is_end = base_plus is None or base_minus is None
            is_full = base_plus in self.staple_co_bases or base_minus in self.staple_co_bases
            is_half = not is_end and not is_full
            if is_end:
                n_staple_end += 1
            elif is_full:
                n_staple_full += 1
            elif is_half:
                n_staple_half += 1
        #ipdb.set_trace()
        return n_staple_half/2., n_staple_full/4., n_staple_end/2.
  
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
    #ipdb.set_trace()
    export_data(data=data, name=name)
    return

if __name__ == "__main__":
    main()


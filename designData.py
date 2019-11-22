#!/usr/bin/env python
from nanodesign.converters import Converter
import nanodesign as nd
import numpy as np
import ipdb  # use this for debugging instead of print() (ipdb.set_trace())
import re

class DesignData(object):

    def __init__(self, name: str):
        self.name: str = name
        self.dna_structure = self.init_design()
        self.all_strands: list = self.dna_structure.strands
        self.pairedbases: list = self.init_pairedbases()
        self.data: dict = {}
        self.all_bases: list = self.get_all_bases()
        self.hps_base: dict = self.init_hps()
        self.all_co_tuples_list: list = self.get_all_co_tuple()
        self.all_co_bases: list = self.get_all_co_bases()
        self.full_co_list: list = self.get_full_co_list()
        self.end_co_set: list = self.get_endloop_co_list()
        self.half_co_list: list = self.get_half_co_list()
        self.scaf_co_bases, self.staple_co_bases = self.get_staple_scaffold_co_bases()
        self.helices_n: list = self.get_helix_per_staple()
        self.helix_dic: dict = self.init_helix_dict()
        self.first_bases, self.last_bases = self.get_first_last_bases_of_strands()
        self.helices = self.get_helices()
        self.nicks: list = self.get_nicks()
        #self.possible_co: list = self.co_density()

         
    def compute_data(self) -> dict:
        data = {}
        
        data["Lattice type"] = self.get_lattice_type()
        data["n_staples"] = self.get_number_of_staple()
        
        st_l_avg, s_l_std, max_l, min_l = self.get_staple_length_statistics()
        data["staple_length_avg"] = st_l_avg
        data["staple_length_std"] = s_l_std
        data["staple_max_length"] = max_l
        data["staple_min_length"] = min_l
        
        avg_h_st_p, std_h_st_p, max_h_st_p, min_h_s_p = self.helices_staples_pass_statistics()
        data["avg_helices_staples_pass"] = avg_h_st_p 
        data["std_helices_staples_pass"] = std_h_st_p
        data["max_helices_a_staples_pass"] = max_h_st_p
        data["min_helices_staples_pass"] = min_h_s_p
        
        data["n_helices"] = len(self.get_helices())
        data["n_skips"] = self.get_n_skips()
        data["n_nicks"] = self.get_n_nicks()
        data["n_co"] = self.get_total_n_co()
        
        avg_n_st_do, std_n_st_do, max_n_st_do, min_n_st_do = self.get_staple_domain_statistics()
        data["avg_n_staple_domain"] = avg_n_st_do 
        data["std_n_staple_domain"] = std_n_st_do
        data["max_n_staple_domain"] = max_n_st_do
        data["min_n_staple_domain"] = min_n_st_do
        
        n_st_l_do, n_st_s_do = self.get_n_staples_with_and_without_long_domains()
        data["n_staple_with_long_domains"] = n_st_l_do
        data["n_staple_with_no_long_domains"] = n_st_s_do
        
        avg_n_l_do, std_n_l_do, max_n_l_do, min_n_l_do = self.get_long_domains_statistic()
        data["avg_n_long_domain_in_staples"] = avg_n_l_do
        data["std_n_long_domain_in_staples"] = std_n_l_do
        data["max_n_long_domain_in_staples"] = max_n_l_do
        data["min_n_long_domain_in_staples"] = min_n_l_do
        
        scaf_co, staple_co = self.get_n_staple_scaffold_co()
        data["scaf_co"] = scaf_co
        data["staple_co"] = staple_co
        
        avg_n_st_co, std_n_st_co, max_n_st_co, min_n_st_co = self.staples_crossovers_statistics()
        data["avg_n_staples_co"] = avg_n_st_co
        data["std_n_staples_co"] = std_n_st_co
        data["max_n_staples_co"] = max_n_st_co
        data["min_n_staples_co"] = min_n_st_co
        
        full_co, half_co, endloop_co = self.get_n_co_types()
        data["n_all_full_co"] = full_co
        data["n_all_half_co"] = half_co
        data["n_all_endloops"] = endloop_co
        
        full_scaf_co, full_staple_co, half_scaf_co, half_staple_co, end_scaf_co, end_staple_co = self.get_n_scaf_staple_co_types()
        data["full_scaff_co"] = full_scaf_co
        data["full_staple_co"] = full_staple_co
        data["half_scaff_co"] = half_scaf_co
        data["half_staple_co"] = half_staple_co
        data["scaff_endloops"] = end_scaf_co
        data["staple_endloops"] = end_staple_co
        
        #bluntends = self.get_blunt_ends()
        #data["n_bluntends"] = len(bluntends)
        
        scaf_co_density, st_co_density, all_co_density, all_co_density_noends = self.get_co_density()
        data["scaf_co_density"] = scaf_co_density
        data["st_co_density"] = st_co_density
        data["all_co_density"] = all_co_density
        data["all_co_density_noends"] = all_co_density_noends
        
        self.data = data
        return self.data

    def init_design(self):
        file_name = self.name + ".json"
        seq_file = self.name + ".seq"
        seq_name = None
        converter = Converter(modify=True)
        converter.read_cadnano_file(file_name, None, "p8064")
        converter.dna_structure.compute_aux_data()
        #ipdb.set_trace()
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
        for strand in self.all_strands:  
            for base in strand.tour:
                position = (base.h, base.p, base.is_scaf)
                hps_base[position] = base
        return hps_base
            
    def get_lattice_type(self):
        if type(self.dna_structure.lattice) == nd.data.lattice.SquareLattice:
            return "Square"
        else:
            return "Honeycomb"
    
    def get_base_from_hps(self, h, p, is_scaffold, dir=1):
        if (h, p) in self.dna_structure.Dhp_skips:
            p += np.sign(dir) 
        return self.hps_base.get((h, p, is_scaffold), None)

        
    def get_hps_from_base(self, base):
        return (base.h, base.p, base.is_scaf)

        
    def get_all_bases(self) -> list:
        all_bases = []
        for strand in self.all_strands:
            for base in strand.tour: 
                all_bases.append(base)
        return all_bases
        
    def get_staple_length_statistics(self):
        len_strands = []
        for strand in self.all_strands:
            if not strand.is_scaffold:
                tour_clean = [base for base in strand.tour]
                len_strands.append(len(tour_clean))
        return (np.average(len_strands), np.std(len_strands),
                np.max(len_strands), np.min(len_strands))
        
    def get_structure_size(self):
        
        return
    
    def get_helices(self) -> int:
        helices = set()
        for strand in self.all_strands:
            if strand.is_scaffold:
                for base in strand.tour:
                    if self.dna_structure._check_base_crossover(base):
                        if base.h not in helices:
                            helices.add(base.h)
        #ipdb.set_trace()
        return helices
                    
    def get_n_skips(self) -> int:
        return len(self.dna_structure.Dhp_skips)

    def get_staple_domain_statistics(self) -> list:
        n_st_domains = []
        for strand in self.all_strands:
            if not strand.is_scaffold:
                n_st_domains.append(len(strand.domain_list))
        #ipdb.set_trace()
        return (np.average(n_st_domains), np.std(n_st_domains), 
                np.max(n_st_domains), np.min(n_st_domains) )     
    
    def get_n_staples_with_and_without_long_domains(self) -> int:
        long_st_domain = []
        short_st_domain = []
        for strand in self.all_strands:
            if not strand.is_scaffold:
                for domain in strand.domain_list:
                    if len(domain.base_list) >= 14:
                        long_st_domain.append(strand)
                if strand not in long_st_domain:
                    short_st_domain.append(strand)
        #ipdb.set_trace()
        return len(long_st_domain), len(short_st_domain)
    
    def get_long_domains_statistic(self) -> int:
        long_st_domain = []
        n_long_domains = []
        for strand in self.all_strands:
            n_long_do = 0
            if not strand.is_scaffold:
                for domain in strand.domain_list:
                    if len(domain.base_list) >= 14:
                        long_st_domain.append(strand)
                        n_long_do +=1 
                n_long_domains.append(n_long_do)
                
        return (np.average(n_long_domains), np.std(n_long_domains),
                np.max(n_long_domains), np.min(n_long_domains))
    
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
            
    def helices_staples_pass_statistics(self) -> int:
        return (np.average(self.helices_n), np.std(self.helices_n),
                np.max(self.helices_n), np.min(self.helices_n))
    
    def init_helix_dict(self) -> dict:
        helix_dic = {}
        helix_list = self.get_helix_per_staple()
        for strand, helices in zip(self.all_strands, helix_list):
            dic = {str(strand): helices}
            helix_dic.update(dic)
        #ipdb.set_trace()
        return helix_dic
    
    def get_first_last_bases_of_strands(self) -> list:
        first_bases = set()
        last_bases = set()
        for strand in self.all_strands:
            if not strand.is_scaffold:
                first_bases.add(strand.tour[0])
                last_bases.add(strand.tour[-1])
        #ipdb.set_trace()
        return first_bases, last_bases
    
    def get_nicks(self) ->int:
        nicks = []
        for base in self.first_bases:
            base_plus = self.get_base_from_hps(base.h, base.p + 1, base.is_scaf)
            base_minus = self.get_base_from_hps(base.h, base.p - 1, base.is_scaf, dir=-1)
            if base_plus in self.last_bases:
                nicks.append((base,base_plus))#order of nick is always (first,last)
            elif base_minus in self.last_bases:
                nicks.append((base,base_minus))
        #ipdb.set_trace()       
        return nicks
    
    def get_n_nicks(self) -> int:
        return len(self.nicks)
    
    def get_all_co_tuple(self) -> list:
        all_co_tuples = set()
        all_co_tuple_list = []
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
        #ipdb.set_trace()
        return all_co_tuple_list
    
    def get_all_co_bases(self):
        all_co_bases = []
        for strand in self.all_strands:
                for base in strand.tour:
                    if self.dna_structure._check_base_crossover(base):
                            all_co_bases.append(base)
        #ipdb.set_trace()
        return all_co_bases
            
    def get_total_n_co(self) -> int:
        return len(self.all_co_bases)/2.
    
    def staples_crossovers_statistics(self) -> int:
        co_staple = []
        n_co_staples = []
        for strand in self.all_strands:
            if not strand.is_scaffold:
                for base in strand.tour:
                    if self.dna_structure._check_base_crossover(base):
                        co_staple.append(base)
                n_co_staples.append(len(co_staple)/2.)
                co_staple = []
        return (np.average(n_co_staples), np.std(n_co_staples),
                np.max(n_co_staples), np.min(n_co_staples))     
        
    def get_full_co_list(self) -> list:
        co_plus_tuples = []
        co_minus_tuples = []
        full_co_list = []
        for co in self.all_co_tuples_list:
            co_tuple_plus = set()
            co_tuple_minus = set()
            for base in co:
                base_plus = self.get_base_from_hps(base.h, base.p + 1, base.is_scaf)
                base_minus = self.get_base_from_hps(base.h, base.p - 1, base.is_scaf, dir=-1)
                co_tuple_plus.add(base_plus)
                co_tuple_minus.add(base_minus)
            #ipdb.set_trace()
            co_plus_tuples.append(frozenset(co_tuple_plus))    
            co_minus_tuples.append(frozenset(co_tuple_minus))
            
            if co_plus_tuples[0] in self.all_co_tuples_list:
                full_co_list.append(co)
            elif co_minus_tuples[0] in self.all_co_tuples_list:
                full_co_list.append(co)
                
            co_plus_tuples = []
            co_minus_tuples = []

        #ipdb.set_trace()    
        return full_co_list
    
    def get_endloop_co_list(self) -> list:    
        end_co_set = set()
        #base_plus_list = []
        for co in self.all_co_tuples_list:   
            for base in co:
                base_plus = self.get_base_from_hps(base.h, base.p + 1, base.is_scaf)
                base_minus = self.get_base_from_hps(base.h, base.p - 1, base.is_scaf, dir=-1)
                #base_plus_list.append(base_plus)
                if (base_plus is None) or (base_minus is None):
                    end_co_set.add(co)
        #ipdb.set_trace()
        return end_co_set
 
    def get_half_co_list(self) -> list:
        half_co_list = []
        for co in self.all_co_tuples_list:
            if (co not in self.end_co_set) and (co not in self.full_co_list):
                half_co_list.append(co)
        #ipdb.set_trace()        
        return half_co_list
    
    def get_n_co_types(self) -> int:
        return len(self.full_co_list)/2., len(self.half_co_list), len(self.end_co_set)
    
    def get_staple_scaffold_co_bases(self):
        scaf_co_bases = []
        staple_co_bases = []
        
        for base in self.all_co_bases:
        #considering the skips in crossovers
            if base.is_scaf:
                scaf_co_bases.append(base)
            else:
                staple_co_bases.append(base)
                            
        return scaf_co_bases, staple_co_bases
    
    def get_n_staple_scaffold_co(self):
        return len(self.scaf_co_bases)/2., len(self.staple_co_bases)/2.
    
    def get_n_scaf_staple_co_types(self):
        full_scaf_co_set = set()
        full_staple_co_set = set()
        half_scaf_co_set = set()
        half_staple_co_set = set()
        end_scaf_co_set = set()
        end_staple_co_set = set()
        for co in self.full_co_list:
            for base in co:
                if base.is_scaf:
                    full_scaf_co_set.add(co)
                else:
                    full_staple_co_set.add(co)
                    
        for co in self.half_co_list:
            for base in co:
                if base.is_scaf:
                    half_scaf_co_set.add(co)
                    #we put it to zero by hand
                    half_scaf_co_set = set()
                else:
                    half_staple_co_set.add(co)
                    
        for co in self.end_co_set:
            for base in co:
                if base.is_scaf:
                    end_scaf_co_set.add(co)
                else:
                    end_staple_co_set.add(co)
                
                
       # ipdb.set_trace()         
        return (len(full_scaf_co_set)/2., len(full_staple_co_set)/2., len(half_scaf_co_set), 
                len(half_staple_co_set), len(end_scaf_co_set), len(end_staple_co_set))

    def get_co_density(self):
        def is_ds(pos, hid):
            is_sc = (hid, pos , True) in self.hps_base
            is_st = (hid, pos, False) in self.hps_base
            is_skip = False # (hid, pos) in self.skips (note: list of (h,p) for all skips)
            return ((is_sc and is_st) or is_skip)
        
        def cleanup_co(co_list):
            n_ends = 0
            if not co_list:
               return 0, 0
            if co_list[0]+1 != co_list[1]:
                n_ends += 1
                co_list = co_list[1:]
            if co_list[-1]-1 != co_list[-2]:
                n_ends += 1
                co_list = co_list[:-1]
            return n_ends, len(co_list)/2
        
        #possible_co = []
        #possible_co.extend(self.dna_structure.connector_points)
        helices = self.dna_structure.structure_helices_map.values()
        n_possible_co_st = 0
        n_possible_co_sc = 0
        n_possible_end_co_st = 0
        n_possible_end_co_sc = 0
        
        for helix in helices:
            stple_co = [co[1] for co in helix.possible_staple_crossovers if is_ds(pos=co[1], hid=helix.id)]
            stple_co = sorted(stple_co)
            scaf_co = [co[1] for co in helix.possible_scaffold_crossovers if is_ds(pos=co[1], hid=helix.id)]
            scaf_co = sorted(scaf_co)
            n_end_st, n_co_st = cleanup_co(stple_co)
            n_possible_co_st += (n_co_st)
            n_possible_end_co_st += n_end_st
            n_end_sc, n_co_sc = cleanup_co(scaf_co)
            n_possible_co_sc += (n_co_sc)
            n_possible_end_co_sc += n_end_sc
        n_possible_co_sc = (n_possible_co_sc)/2  
        n_possible_co_st = (n_possible_co_st)/2
        n_possible_end_co_sc = (n_possible_end_co_sc)/2
        n_possible_end_co_st = (n_possible_end_co_st)/2
         
        full_scaf_co, full_staple_co, half_scaf_co, half_staple_co, end_scaf_co, end_staple_co = self.get_n_scaf_staple_co_types()   
        
        scaf_co_density = (full_scaf_co + half_scaf_co)/ (n_possible_co_sc)
        
        st_co_density = (full_staple_co + half_staple_co) / (n_possible_co_st)
        
        all_co_density = ((full_scaf_co + full_staple_co + half_scaf_co + 
                          half_staple_co +  end_scaf_co + end_staple_co) /
                        (n_possible_co_sc + n_possible_co_st + 
                         n_possible_end_co_sc + n_possible_end_co_st ))
                           
        all_co_density_noends = (full_scaf_co + full_staple_co + half_scaf_co + 
                                 half_staple_co) / (n_possible_co_sc + n_possible_co_st)

        #ipdb.set_trace()    
        return scaf_co_density, st_co_density, all_co_density, all_co_density_noends
    
    def get_blunt_ends(self):
        blunt_ends = set()
        blunt_end = set()
        for co in self.end_co_set:
            for base in co:
                    for co_s in self.all_co_tuples_list:
                        for base_s in co_s:
                            if not base_s.is_scaf:
                                if (base.h == base_s.h) and (base.p == base_s.p):
                                    blunt_end.add(frozenset(co))
                                    blunt_end.add(frozenset(co_s))
                                    blunt_ends.add(frozenset(blunt_end)) 
        return blunt_ends
    
def export_data(data: dict, name: str) -> None:
    header = ", ".join(data.keys())
    data_str = ", ".join([str(i) for i in data.values()])
    with open(name+"-designdata.csv", mode="w+") as out:
        out.write(header + "\n")
        out.write(data_str)
        out.write("\nEND")
    return

def main():
    file  = open("txt_file.txt", 'rt', encoding="utf8")
    for i, line in enumerate(file):
        if i == 3:
            if line.startswith('Project =  '):
                name = line[11:-1]
    #print("master, I am awaiting the name of your design")
    #name = input()
    #print("Thank you Sir")
    designData = DesignData(name=name)
    data = designData.compute_data()
    export_data(data=data, name=name)
    return

if __name__ == "__main__":
    main()


from nanodesign.converters import Converter

from nanodesign.data import dna_structure_helix

import numpy as np

import scipy.ndimage as nd

import math





def get_design(path):

    file_name = path + ".json"

    seq_file = path + ".seq"

    seq_name = None



    converter = Converter()

    converter.read_cadnano_file(file_name, seq_file, seq_name)

    converter.dna_structure.get_domains()

    converter.dna_structure._compute_helix_design_crossovers()



    return converter.dna_structure

raw_input



dna_structure = get_design(raw_input())
#'./TTcorr'
# list_wo_skips = [base for base in list_of_bases if base.num_deletions == 0]

all_strands = dna_structure.strands

len_strands = []

# for i in range(0,len(all_strands)):

# m    len_strands.append(len(all_strands[i].tour))

# len_strands.remove(np.max(len_strands))

for strand in all_strands:

    if not strand.is_scaffold:

        tour_clean = [base for base in strand.tour if base.num_deletions == 0]

        len_strands.append(len(tour_clean))


num_skips = 0

for strand in all_strands:

    # if not strand.is_scaffold:

    skips = [base for base in strand.tour if base.num_deletions == -1]

    num_skips += len(skips)

# print("number of Skips = " , num_bases)


# for i in range(len(skips)):

#   print(skips[i].h, skips[i].p)

num_strands = []

for strand in all_strands:

    if not strand.is_scaffold:

        num_strands.append(strand)

# print('Number of Strands = ',len(num_strands))

N_co = 0
N_d_co = 0
for strand in all_strands:

    for base in strand.tour:

        if dna_structure._check_base_crossover(base):

            N_co += 1
        if dna_structure._check_base_double_crossover(base):
         
            N_d_co += 1


header = "AvgStapleLength, StdStapleLength, N_skips, N_staples, N_crossovers, N_double_crossovers"

data = [np.average(len_strands), np.std(len_strands),

        num_skips, len(num_strands), N_co, N_d_co]



with open("data.txt", mode="w+") as out:

    out.write(header + "\n")

    for x in data:

        out.write(str(x) + ", ")

    out.write("\nEND")


def main():

    return

if __name__ == "__main__":

    main()
    
    
    
    
    
    
    
    
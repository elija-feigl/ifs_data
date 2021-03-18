# ifs_data: dietzlab initial folding screen database

A Python3 package that contains a collection of scripts related to the analysis of initial folding screens of DNA-Origami structures.


# Usage
...

# Dependencies

* Python >= 3.7
  * numpy >= 1.14
  * scipy >= 1.1
  * attrs >= ??
  * nanodesign3  >= ?? # not publised yet


# Installation

install all required Python Packages. the following are not available via the PyPI
## nanodesign3
https://github.com/elija-feigl/nanodesign_dietz


git clone https://github.com/elija-feigl/ifs_data
cd ifs_data
python setup.py install
cd ..


# Citation

When using ifs_data Python Package in published work, please cite the following paper:

...


# documentation of design_statistics

"name": name of the structure.

"Lattice type": lattice type if the stucture.
	possible values: Honeycomb, square.

"n_helices": number of helices in the structure.

"n_deletions": number of deletions in the structure.

"n_nicks": number of the nicks in the structure.

"n_staple_domain": number of domains for each staple. 
	possible values: Avg, Std, Max, Min numbr of staples domains.

"2_long_domains": number of staples with more than 2 long domains.

"1_long_domains": number of staples with 1 long domains.

"0_long_domains": number of staples with no long domains.

"co_rule_violation": number of domains with less than 5 bases which violates crossover rule.

"n_staples": number of staples for the structure.
	
"staple_length": lengh of each staple.
	possible values: Avg, Std, Max, Min number of staple lengths.

"helices_staples_pass":numbers of helices each staple passes through.
	possible values: Avg, Std, Max, Min number of helices a staple passes.

##################################################################################################
"co_set-scaffold-full": number of scaffold full crossovers.

"co_set-scaffold-full-h": number of scaffold horizontal full crossovers.

"co_set-scaffold-full-v": number of scaffold vertical full crossovers.

"co_set-scaffold-half": number of scaffold half crossovers.

"co_set-scaffold-half-h": number of scaffold horizontal half crossovers.

"co_set-scaffold-half-v": number of scaffold vertical half crossovers.

"co_set-scaffold-end": number of scaffold endloops.

"co_set-scaffold-end-h": number of scaffold horizontal endloops.

"co_set-scaffold-end-v": number of scaffold vertical endloops.

"co_set-scaffold-co": number of scaffold crossovers.

"co_set-scaffold-co-h": number of scaffold horizontal crossovers.

"co_set-scaffold-co-v": number of scaffold vertical crossovers.

##################################################################################################

"co_set-staple-full": number of staple full crossovers.

"co_set-staple-full-h": number of staple horizontal full crossovers.

"co_set-staple-full-v": number of staple vertical full crossovers.

"co_set-staple-half": number of staple half crossovers.

"co_set-staple-half-h": number of staple horizontal half crossovers.

"co_set-staple-half-v": number of staple vertical half crossovers.

"co_set-staple-end": number of staple endloops.

"co_set-staple-end-h": number of staple horizontal endloops.

"co_set-staple-end-v": number of staple vertical endloops.

"co_set-staple-co" : number of staple crossovers.

"co_set-staple-co-h" : number of staple horizontal crossovers.

"co_set-staple-co-v" : number of scaffold vertical crossovers.

##################################################################################################

"co_possible-scaffold-co": number of possible scaffold crossovers.

"co_possible-scaffold-co-h": number of possible scaffold horizontal crossovers.

"co_possible-scaffold-co-v": number of possible scaffold vertical crossovers.

"co_possible-scaffold-end": number of possible scaffold enloops.

##################################################################################################

"co_possible-staple-co": number of possible staple crossovers.

"co_possible-staple-co-h": number of possible staple horizontal crossovers.

"co_possible-staple-co-v": number of possible staple vertical crossovers.

"co_possible-staple-end": number of possible staple enloops.

##################################################################################################

"co_density-scaffold-co": number of scaffold crossovers divided by number of possible scaffold crossovers.

"co_density-scaffold-co-h": number of scaffold horizontal crossovers divided by number of possible scaffold horizontal crossovers.

"co_density-scaffold-co-v": number of scaffold vertical crossovers divided by number of possible scaffold crossovers.

"co_density-scaffold-end": number of scaffold endloops divided by number of possible scaffold endloops.

"co_density-staple-co": number of staple crossovers divided by number of possible staple crossovers.

"co_density-staple-co-h": number of staple horizontal crossovers divided by number of possible staple horizontal crossovers.

"co_density-staple-co-v": number of staple vertical crossovers divided by number of possible staple crossovers.

"co_density-staple-end": number of staple endloops divided by number of possible staple endloops.

"n_stacks":length of different stacks in the structure.(how many full crossovers are aligned subsequently)
	possible values: Avg, Std, Max, Min stack lengths of the structure.

"stacks" : number of stacks.

"del_density": number of deletions(deletions) divided by number of all bases including the deletions.

"ins_density": number of insertions divided by number of all bases.

"n_blunt_ends": number of blunt ends.




























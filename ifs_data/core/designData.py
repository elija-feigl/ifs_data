import numpy as np
import logging
import itertools
import attr
import networkx as nx

from typing import List, Dict, Tuple, Set, Optional
from collections import defaultdict

import nanodesign as nd
from nanodesign.converters import Converter
from nanodesign.data.base import DnaBase as Base
from nanodesign.data.strand import DnaStrand as Strand
from nanodesign.data.dna_structure import DnaStructure as Structure

from Bio.SeqUtils import MeltingTemp
from Bio.Seq import Seq

from ..data.crossover import Crossover, Connection
from ..data.nicks import Nick
from .utils import _hps, _close_strand, _check_base_crossover, _change_strand_type


@attr.s
class Design(object):
    json: str = attr.ib()
    name: str = attr.ib()
    seq: str = attr.ib()
    circ_scaffold: bool = attr.ib(default=True)

    def __attrs_post_init__(self):

        self.logger = logging.getLogger(__name__)
        logging.getLogger("nanodesign").setLevel(logging.WARNING)

        self.dna_structure = self._init_design()

        # NOTE: methods called by _init... are named _create...
        self._init_configuration()
        self._init_topology()
        self._init_domains()
        self._init_connectivity()

    ###########################################################################
    # initialisation
    def _init_design(self) -> Structure:
        """ convert into nanodesign structure"""
        converter = Converter()
        converter.modify = True
        converter.read_cadnano_file(self.json, None, self.seq)
        converter.dna_structure.compute_aux_data()
        return converter.dna_structure

    def _init_configuration(self) -> None:
        """ retrive Origami specific information"""
        self.strands = self.dna_structure.strands
        self._fix_strand_classification()
        self.scaffolds = self._create_scaffolds()
        self.staples = self._create_staples()
        self.lattice = self._create_lattice()
        self.helix_maps = self.dna_structure.structure_helices_map
        self.helices = list(self.dna_structure.structure_helices_map.values())

    def _init_topology(self) -> None:
        """ add topological data required for direct access of bases"""
        self.hps_base = self._create_hps()
        self.hps_deletions, self.hps_insertions = self._create_hps_modifications()
        self.nicks = self._create_nicks()

    def _init_domains(self) -> None:
        """ everything related to domains (distribution & metling T)"""
        self.domain_data = self._create_domain_data()
        self.max_staple_melt_t = self._create_staple_max_melt_T()

    def _init_connectivity(self) -> None:
        """ crossover data
                full crossovers are composed of two connections
                half crossovers and endloops of one connection
            stacks
        """
        self.crossovers = self._create_crossovers()
        self.stacks = self._create_stacks()

    ###########################################################################
    # attribute creators
    def _fix_strand_classification(self, scaffold_min_length: int = 500) -> None:
        """scaffolds-strands smaller than scaffold_min_length bases are to be consedered as staples."""
        for strand in self.strands:
            base = strand.tour[0]
            is_false_scaf = (strand.is_scaffold and len(
                strand.tour) < scaffold_min_length)
            is_false_staple = (not strand.is_scaffold and len(
                strand.tour) > scaffold_min_length)
            if is_false_scaf or is_false_staple:
                logger_string = "staple" if is_false_scaf else "scaffold"
                self.logger.debug(
                    f"strand starting at{(base.h, base.p)} considered as {logger_string}")
                _change_strand_type(strand)

    def _create_scaffolds(self) -> List[Strand]:
        """ regular cadnano designs have a linear scaffold to indicate sequence start even for circular scaffolds."""
        scaffolds = [strand for strand in self.strands if strand.is_scaffold]
        if self.circ_scaffold:
            for strand in scaffolds:
                _close_strand(strand)
        return scaffolds

    def _create_staples(self) -> List[Strand]:
        """ NOTE: call afer _create_scaffolds to catch converted strands"""
        return [strand for strand in self.strands if not strand.is_scaffold]

    def _create_hps(self) -> Dict[Tuple[int, int, bool], Base]:
        """create a dictionary from cadnano bases positions to bases."""
        return {_hps(base): base for strand in self.strands for base in strand.tour}

    def _create_hps_modifications(self) -> Tuple[List, List]:
        """ create a list of all modified cadnano bases positions."""
        # TODO: implement without reusing converter
        converter = Converter()
        converter.modify = False
        converter.read_cadnano_file(self.json, None, self.seq)
        converter.dna_structure.compute_aux_data()
        dna_structure_del_ins = converter.dna_structure

        hps_deletions = list()
        hps_insertions = list()
        for strand in dna_structure_del_ins.strands:
            for base in strand.tour:
                if base.num_insertions != 0:
                    hps_insertions.append(_hps(base))
                elif base.num_deletions != 0:
                    hps_deletions.append(_hps(base))
        return hps_deletions, hps_insertions

    def _create_nicks(self) -> List[Nick]:
        """ order of nick is always (first base,last base)"""
        def _2(base):
            base_plus, base_minus = self._get_base_plus_minus(base)
            if base_plus in last_bases:
                return base_plus
            elif base_minus in last_bases:
                return base_minus
            else:
                return False

        first_bases = {staple.tour[0] for staple in self.staples}
        last_bases = {staple.tour[-1] for staple in self.staples}
        return [Nick(base, _2(base)) for base in first_bases if _2(base)]

    def _create_lattice(self) -> str:
        is_square = (type(self.dna_structure.lattice)
                     == nd.data.lattice.SquareLattice)
        return "square" if is_square else "honeycomb"

    def _create_domain_data(self) -> Dict[Strand, List]:
        """ dict linking every staple to a list of its domains."""
        return {staple: staple.domain_list for staple in self.staples}

    def _create_staple_max_melt_T(self) -> Dict[Strand, float]:
        """ max_melt_T is the staple domain with the highest metling temperature"""
        staple_domains_melt_t: Dict[Strand, List[float]] = dict()
        for staple in self.staples:
            domains = staple.domain_list
            for domain in domains:
                if "N" not in domain.sequence:
                    # NOTE: using nearest neighbor for domain with length higher
                    #   than 14 using 'Wallace rule' else
                    if len(domain.base_list) > 14:
                        staple_domains_melt_t.setdefault(staple, []).append(MeltingTemp.Tm_NN(
                            Seq(domain.sequence), Na=0, Mg=17.5))
                    else:
                        staple_domains_melt_t.setdefault(staple, []).append(MeltingTemp.Tm_Wallace(
                            Seq(domain.sequence)))
        max_staple_melt_t = {key: max(value) for (
            key, value) in staple_domains_melt_t.items()}
        return max_staple_melt_t

    def _create_all_connections(self) -> Set[Connection]:
        """get a set of all interhelical connection of two bases."""
        all_con = set()
        passed_bases = set()
        for strand in self.strands:
            for base1 in (b for b in strand.tour if b not in passed_bases):
                base2 = _check_base_crossover(base1)
                if base2:
                    con = Connection(base1, base2)
                    passed_bases.add(base1)
                    passed_bases.add(base2)

                    all_con.add(con)
        return all_con

    def _create_crossovers(self) -> List[Crossover]:
        """organize connections into crossover objects with full syntax and typ annotation"""
        def find_second_connections(base: Base) -> Optional[Connection]:
            """ check for possible 2nd connection for full crossover."""
            if not base:
                return None
            for con2 in _connections:
                con2_bases = con2.get_bases()
                is_base_in_con = base in con2_bases
                is_same_hs = {b.h for b in con2_bases} == {
                    b.h for b in con2_bases}
                if is_base_in_con and is_same_hs:
                    _connections.remove(con2)
                    connected_bases.difference_update(con2.get_bases())
                    return con2
            else:
                self.logger.debug(
                    f"Connection {(con2_base.h, con2_base.p, con2_base.is_scaf)}  ignored for Crossovers. Likely odd use of Crossover in Cadnano")
                return None

        _connections = self._create_all_connections()
        connected_bases = {b for con in _connections for b in con.get_bases()}

        crossovers = list()
        while _connections:
            con1 = _connections.pop()
            con1_bases = con1.get_bases()
            connected_bases.difference_update(con1_bases)

            is_ssDNA_loop = (None in [b.across for b in con1_bases])
            if is_ssDNA_loop:
                self.logger.debug(
                    f"connection {(con1.base1.h, con1.base1.p)}  removed as ss loop")
                continue

            base_plus, base_minus = self._get_base_plus_minus(con1.base1)
            is_end_connection = (base_plus is None) or (base_minus is None)
            if is_end_connection:
                is_end, con2 = True, None
            else:
                is_end = False
                if base_plus in connected_bases:
                    con2_base = base_plus
                elif base_minus in connected_bases:
                    con2_base = base_minus
                else:
                    con2_base = False

                con2 = find_second_connections(con2_base)
            co = Crossover(con1, con2, is_end)
            co.set_orientation(self.helix_maps)
            co.set_scaffold_subtyp(self.helix_maps, self.lattice)
            crossovers.append(co)
        return crossovers

    def _create_stacks(self) -> Dict[int, List[Set[int]]]:
        """ A stack is a Sets of connected helices (helixID)
            Stacks are sorted to Lists according to their position (basePosition)"""
        def _sort_by_position():
            position_crossovers = defaultdict(list)
            for co in full_crossovers:
                # NOTE: assuming "perfect crossovers": position is marked by left base in left connection
                co_p = min(base.p for base in co.get_bases())
                position_crossovers[co_p].append(
                    (co.connection1.base1.h, co.connection1.base2.h))
            return position_crossovers

        def _sort_to_stacks(position_crossovers):
            stacks = defaultdict(list)
            for stack_p, edges in position_crossovers.items():
                G = nx.Graph()
                G.add_edges_from(edges)
                connected = list(
                    comp for comp in nx.connected_components(G) if len(comp) > 2)
                if connected:
                    self.logger.debug(f"Found stack at {stack_p}: {connected}")
                    stacks[stack_p] = connected
            return stacks

        full_crossovers = [co for co in self.crossovers if co.typ == "full"]
        same_pos = _sort_by_position()
        stacks = _sort_to_stacks(same_pos)
        return stacks

    ###########################################################################
    # internal getters
    def _get_neighbor_from_hps(self, h: int, p: int, s: bool, dir: int = 1) -> Base:
        """ get the base object from its coordination: (h, p is_scaffold)
            jumps deletions into given direction dir (default: {1})"""
        if (h, p, s) in self.hps_deletions:
            p += np.sign(dir)
        return self.hps_base.get((h, p, s), None)

    def _get_base_plus_minus(self, base: Base) -> Tuple[Base, Base]:
        """ given a base, it returns the neighbour bases."""
        base_plus = self._get_neighbor_from_hps(
            base.h, base.p + 1, base.is_scaf)
        base_minus = self._get_neighbor_from_hps(
            base.h, base.p - 1, base.is_scaf, dir=-1)
        return base_plus, base_minus

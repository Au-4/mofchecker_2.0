# -*- coding: utf-8 -*-
"""Check if there are fused rings with 5-ring N which is possible charged and have to be checked manaully."""
from pymatgen.analysis.graphs import StructureGraph

from mofchecker.types import StructureIStructureType
from collections import Counter
from .base_coordination_check import BaseCoordinationCheck
from ..utils.get_indices import get_n_indices, get_c_indices,num_neighbor_metal, non_metal_neighbors,num_neighbor_H, non_H_neighbors
from .geometry import _guess_underbound_nitrogen_cn2


class Fusedring_Check(BaseCoordinationCheck):


    def __init__(self, structure: StructureIStructureType, structure_graph: StructureGraph):
        """Initialize the fusedring check.

        Args:
            structure (StructureIStructureType): The structure to check.
            structure_graph (StructureGraph): The structure graph to use for the check.
        """
        self.structure = structure
        self.c_indices = get_c_indices(self.structure)
        self.n_indices = get_n_indices(self.structure)
        self.structure_graph = structure_graph

    @property
    def name(self):
        """Return the name of the check."""
        return "Possible charged fused ring with N"

    @property
    def description(self):
        """Return a description of the check."""
        return "Checks if there are fused ring with 5-ring N."

    def _run_check(self):
        fused_ring = self._get_fused_ring()
        return len(fused_ring) == 0, fused_ring
    def _get_fused_ring(self):
        N_sum = []
        N_jump = []
        connection_pool = []
        for site_index in self.n_indices:
            if site_index not in N_jump:
                cn = self.get_cn(site_index)
                neighbors = self.get_connected_sites(site_index)                
                cm = num_neighbor_metal(neighbors)
                ch = num_neighbor_H(neighbors)
                if cn-cm-ch==2:
                    non_H = non_H_neighbors(neighbors)
                    non_metals = non_metal_neighbors(non_H)
                    undercoordinated_nitrogen = _guess_underbound_nitrogen_cn2(
                    self.structure,
                    site_index,
                    non_metals,
                    self.get_connected_sites(non_metals[0].index),
                    self.get_connected_sites(non_metals[1].index),
                    25,
                )
                    if not undercoordinated_nitrogen:
                        ring_2 = non_metals[0]
                        ring_3 = non_metals[1]
                        neighbor_0_neighbors = self.get_connected_sites(non_metals[0].index)
                        neighbor_1_neighbors = self.get_connected_sites(non_metals[1].index)
                        for neighbor_0_neighbor in neighbor_0_neighbors:
                            for neighbor_1_neighbor in neighbor_1_neighbors:
                                dis4_5 = self.structure.get_distance(neighbor_0_neighbor.index,neighbor_1_neighbor.index)
                                if dis4_5 < 1.85 and dis4_5 > 1:
                                    ring_4 = neighbor_0_neighbor
                                    ring_5 = neighbor_1_neighbor
                                    rings = [ring_2,ring_3,ring_4,ring_5]
                                    rings_index = [site_index, ring_2.index, ring_3.index, ring_4.index, ring_5.index]
                                    N_sum.append(site_index)
                                    for ring in rings:
                                        connection_pool.append(ring.index)
                                        ring_neighbors = non_metal_neighbors(self.get_connected_sites(ring.index))
                                        for ring_neighbor in ring_neighbors:
                                            if ring_neighbor.index not in rings_index:
                                                connection_pool.append(ring_neighbor.index)
                                        if str(ring.site.specie) == 'N':
                                            N_jump.append(ring.index)
        connection_count=Counter(connection_pool)
        i = 0 
        for element, count in connection_count.items(): 
            if count > 1:
                i = i+1
        if i >= len(N_sum):
            return N_sum
        else :
            N_sum = []
            return N_sum
                                   
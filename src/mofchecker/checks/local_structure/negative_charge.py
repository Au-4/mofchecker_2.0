# -*- coding: utf-8 -*-
"""Check the negative charge from linkers like N, O, S....."""
from pymatgen.analysis.graphs import StructureGraph
from mofchecker.checker_types import StructureIStructureType
from structuregraph_helpers.create import construct_clean_graph
from .base_coordination_check import BaseCoordinationCheck
from ..utils.get_indices import get_n_indices, get_o_indices, get_halogen_indices,\
    num_neighbor_metal, num_neighbor_halogen, get_s_indices, non_metal_neighbor, non_metal_neighbors, is_metal
from .geometry import _guess_underbound_nitrogen_cn2, get_angle_between_site_and_neighbors
import networkx as nx
import numpy as np
class Negative_charge_Check(BaseCoordinationCheck):

    def __init__(self, structure: StructureIStructureType, structure_graph: StructureGraph):
        """Initialize a new PositivechargeCheck.

        Args:
            structure (StructureIStructureType): Structure to check.
            structure_graph (StructureGraph): StructureGraph of the structure.
        """
        self.structure = structure
        self.n_indices = get_n_indices(self.structure)
        self.o_indices = get_o_indices(self.structure)
        self.s_indices = get_s_indices(self.structure)
        self.halogen_indices = get_halogen_indices(self.structure)
        self.structure_graph = structure_graph

    @property
    def name(self):
        """Return the name of the check."""
        return "Negative charge from linkers"

    @property
    def description(self):
        """Return a description of the check."""
        return "Check the negative charge from the linkers."

    def _run_check(self):
        number = self._get_halogen() + self._get_O_charge() + self._get_S_charge() + self._get_N_charge()
        return len(number) == 0, number

    def _get_halogen(self):
        """Check for all negative charge from halogen."""
        halogen_sum = []
        F_jump = []
        for site_index in self.halogen_indices:
            if site_index not in F_jump:
                cn = self.get_cn(site_index)
                neighbors = self.get_connected_sites(site_index)                
                cm = num_neighbor_metal(neighbors)
                non_metal = non_metal_neighbor(neighbors)
                cf = num_neighbor_halogen(self.get_connected_sites(non_metal))
                if cn-cm==0:
                    halogen_sum.append(site_index)
                    """Check for all negative charge from halogen-non_metal_element cluster 
                    like (BF4)-, (SiF5)-, (SiF6)2-, (PF6)-..."""
                elif str(non_metal.site.specie) == 'B' and cf == 4:
                    halogen_sum.append(site_index)
                    for neighbor in self.get_connected_sites(non_metal):
                        F_jump.append(neighbor.index)
                elif str(non_metal.site.specie) == 'Si' and cf == 5:
                    halogen_sum.append(site_index)
                    for neighbor in self.get_connected_sites(non_metal):
                        F_jump.append(neighbor.index)
                elif str(non_metal.site.specie) == 'Si' and cf == 6:
                    halogen_sum.append(site_index)
                    halogen_sum.append(site_index)
                    for neighbor in self.get_connected_sites(non_metal):
                        F_jump.append(neighbor.index)
                elif str(non_metal.site.specie) == 'P' and cf == 6:
                    halogen_sum.append(site_index)
                    for neighbor in self.get_connected_sites(non_metal):
                        F_jump.append(neighbor.index)                
                else:
                    continue
        return halogen_sum
    
    def _get_O_charge(self):
        """Check for all negative charge from O."""
        O_sum = []
        O_jump = []
        for site_index in self.o_indices:
            if site_index not in O_jump:
                cn = self.get_cn(site_index)
                neighbors = self.get_connected_sites(site_index)                
                cm = num_neighbor_metal(neighbors)
                """check for single O atom"""
                if cn-cm==0:
                    O_sum.append(site_index)
                    O_sum.append(site_index)
                    """check for netural O"""
                elif cn-cm==2:
                    continue

                elif cn-cm==1:
                    non_metal = non_metal_neighbor(neighbors)
                    """check for single OH group""" 
                    if str(non_metal.site.specie) == 'H':
                        O_sum.append(site_index)
                        """check for C-O connection"""
                    elif str(non_metal.site.specie) == 'C':
                        neighbor_C = non_metal.index
                        disC_O = self.structure.get_distance(neighbor_C,site_index)
                        C_neighbors = self.get_connected_sites(neighbor_C)
                        num_O = 0
                        n_O = 0
                        num_n = 0
                        for C_neighbor in C_neighbors:
                                C_neighbor_site = C_neighbor.site
                                neighbor_cn = self.get_cn(C_neighbor.index)
                                neighbor_cm = num_neighbor_metal(self.get_connected_sites(C_neighbor.index))
                                if str(C_neighbor_site.specie) == 'N':
                                    num_n = num_n+1
                                if str(C_neighbor_site.specie) == 'O':
                                    n_O = n_O+1
                                if str(C_neighbor_site.specie) == 'O' and neighbor_cn-neighbor_cm ==1:
                                    num_O = num_O+1
                                    """only have to check one atom out of two O from carboxylic group"""
                                    if C_neighbor.index != site_index:
                                        O_possible_jump = C_neighbor.index
                        """check for charged carboxylic group by finding one C atom coordinated with two O atoms and without H or R on O site"""
                        if num_O == 2:
                            O_sum.append(site_index)
                            O_jump.append(O_possible_jump)
                        """check for alcohol or phenol by analysing C-O bond distance to predict if it is a single bond"""
                        if n_O ==1:
                            if disC_O > 1.315:
                                O_sum.append(site_index)
                            elif disC_O <= 1.315 and num_n == 1:
                                O_sum.append(site_index)

                        """check for N-O connection"""
                    elif str(non_metal.site.specie) == 'N':

                        neighbor_N = non_metal.index
                        N_neighbors = self.get_connected_sites(neighbor_N)
                        neighbor_cn = self.get_cn(neighbor_N)
                        num_O = 0
                        O_NOx = []
                        num_metal = num_neighbor_metal(self.get_connected_sites(neighbor_N))
                        for N_neighbor in N_neighbors:
                            if str(N_neighbor.site.specie) == 'O':
                                num_O = num_O+1
                                O_NOx.append(N_neighbor)

                        """check for charged nitrate group rather than nitro group"""

                        negative_O = 1-(neighbor_cn-num_metal-num_O)
                        if negative_O > 0:
                            for O_atom in O_NOx:
                                O_NOx_cn = self.get_cn(O_atom.index)
                                O_NOx_neighbors = self.get_connected_sites(O_atom.index)                
                                O_NOx_cm = num_neighbor_metal(O_NOx_neighbors)
                                if O_atom.index != site_index:
                                    O_jump.append(O_atom.index)
                                if O_NOx_cn-O_NOx_cm != 1:
                                    negative_O = negative_O-1
                            if negative_O==1:
                                O_sum.append(site_index)
                            elif negative_O==0:
                                continue

                        """check for S-O connection"""
                    elif str(non_metal.site.specie) == 'S':                   
                        neighbor_S = non_metal.index
                        S_neighbors = self.get_connected_sites(neighbor_S)
                        neighbor_cn = self.get_cn(neighbor_S)
                        num_O = 0
                        O_SOx = []
                        num_metal = num_neighbor_metal(self.get_connected_sites(neighbor_S))
                        for S_neighbor in S_neighbors:
                            if str(S_neighbor.site.specie) == 'O':
                                num_O = num_O+1
                                O_SOx.append(S_neighbor)

                        """check for sulfate,sulfite,sulfonic group and sulfinic group SOx"""

                        negative_O = 2-(neighbor_cn-num_metal-num_O)
                        if negative_O > 0:
                            for O_atom in O_SOx:
                                O_SOx_cn = self.get_cn(O_atom.index)
                                O_SOx_neighbors = self.get_connected_sites(O_atom.index)                
                                O_SOx_cm = num_neighbor_metal(O_SOx_neighbors)
                                if O_atom.index != site_index:
                                    O_jump.append(O_atom.index)
                                if O_SOx_cn-O_SOx_cm != 1:
                                    negative_O = negative_O-1
                            if negative_O==2:
                                O_sum.append(site_index)
                                O_sum.append(site_index)
                            elif negative_O==1:
                                O_sum.append(site_index)
                            elif negative_O==0:
                                continue

                        """check for P-O connection"""
                    elif str(non_metal.site.specie) == 'P':
                        neighbor_P = non_metal.index
                        P_neighbors = self.get_connected_sites(neighbor_P)
                        neighbor_cn = self.get_cn(neighbor_P)
                        num_O = 0
                        O_POx = []
                        num_metal = num_neighbor_metal(self.get_connected_sites(neighbor_P))
                        for P_neighbor in P_neighbors:
                            if str(P_neighbor.site.specie) == 'O':
                                num_O = num_O+1
                                O_POx.append(P_neighbor)
                        """check for POx"""

                        negative_O = 3-(neighbor_cn-num_metal-num_O)

                        if negative_O > 0:
                            for O_atom in O_POx:
                                O_POx_cn = self.get_cn(O_atom.index)
                                O_POx_neighbors = self.get_connected_sites(O_atom.index)                
                                O_POx_cm = num_neighbor_metal(O_POx_neighbors)
                                if O_atom.index != site_index:
                                    O_jump.append(O_atom.index)
                                if O_POx_cn-O_POx_cm != 1:
                                    negative_O = negative_O-1

                            if negative_O==2:
                                O_sum.append(site_index)
                                O_sum.append(site_index)
                            elif negative_O==1:
                                O_sum.append(site_index)
                            elif negative_O==0:
                                continue
                            elif negative_O == 3:
                                O_sum.append(site_index)
                                O_sum.append(site_index)
                                O_sum.append(site_index)                                

                        """check for halogen-O connection like ClO4-, IO3-..."""
                    elif str(non_metal.site.specie) == 'Cl' or str(non_metal.site.specie) == 'Br' or str(non_metal.site.specie) == 'I':
                        neighbor_F = non_metal.index
                        F_neighbors = self.get_connected_sites(neighbor_F)
                        neighbor_cn = self.get_cn(neighbor_F)
                        num_O = 0
                        O_ClOx = []
                        num_metal = num_neighbor_metal(self.get_connected_sites(neighbor_F))
                        for F_neighbor in F_neighbors:
                            if str(F_neighbor.site.specie) == 'O':
                                num_O = num_O+1
                                O_ClOx.append(F_neighbor)

                        """check for charged group ClOx"""

                        negative_O = 1-(neighbor_cn-num_metal-num_O)
                        if negative_O > 0:
                            for O_atom in O_ClOx:
                                O_ClOx_cn = self.get_cn(O_atom.index)
                                O_ClOx_neighbors = self.get_connected_sites(O_atom.index)                
                                O_ClOx_cm = num_neighbor_metal(O_ClOx_neighbors)
                                if O_atom.index != site_index:
                                    O_jump.append(O_atom.index)
                                if O_ClOx_cn-O_ClOx_cm != 1:
                                    negative_O = negative_O-1
                            if negative_O==1:
                                O_sum.append(site_index)
                            elif negative_O==0:
                                continue


                        """other non_metal elements connection"""
                    else:
                        O_sum.append(site_index)
                else:
                    continue

        return O_sum
    
    def _get_S_charge(self):
        """Check for all negative charge from S."""
        S_sum = []
        S_jump = []
        for site_index in self.s_indices:
            if site_index not in S_jump:
                
                cn = self.get_cn(site_index)
                neighbors = self.get_connected_sites(site_index)                
                cm = num_neighbor_metal(neighbors)

                """check for single S atom"""
                if cn-cm==0:
                    S_sum.append(site_index)
                    S_sum.append(site_index)
                    """check for netural S"""
                elif cn-cm==2:
                    continue
                    """check for thiol group"""
                elif cn-cm==1: 
 
                    non_metal = non_metal_neighbor(neighbors)
                    """check for single SH group""" 
                    if str(non_metal.site.specie) == 'H':
                        S_sum.append(site_index)
                        """check for C-S connection"""
                    elif str(non_metal.site.specie) == 'C':
                        neighbor_C = non_metal.index
                        disC_S = self.structure.get_distance(neighbor_C,site_index)
                        """check for one C atom coordinated with two S atoms and no H or R on S"""
                        C_neighbors = self.get_connected_sites(neighbor_C)
                        num_S = 1
                        negative_S = 0
                        num_n = 0
                        for C_neighbor in C_neighbors:
                            C_neighbor_site = C_neighbor.site
                            neighbor_cn = self.get_cn(C_neighbor.index)
                            neighbor_cm = num_neighbor_metal(self.get_connected_sites(C_neighbor.index))
                            if str(C_neighbor_site.specie) == 'S' and C_neighbor.index != site_index:
                                num_S = num_S+1
                                S_possible_jump = C_neighbor.index
                            if str(C_neighbor_site.specie) == 'N':
                                num_n = num_n+1
                            if neighbor_cn-neighbor_cm ==1:
                                negative_S = negative_S + 1
                            """S-C single bonds usually contribute coornation with metal but S-C double bonds don't"""
                        if num_S == 1:
                            if disC_S > 1.66:
                                S_sum.append(site_index)
                            elif disC_S <= 1.66 and num_n == 1:
                                S_sum.append(site_index)
                        if num_S == 2 and negative_S >= 1:
                            S_sum.append(site_index)
                            S_jump.append(S_possible_jump)
                    else :
                        S_sum.append(site_index)

                    '''already check S in SOx'''
                else:
                    continue
        return S_sum
    
    def _get_N_charge(self):
        """Check for all negative charge from N. In some rare case there is amide anion and will also be checked as undercoordinated N .
           Besides of this, we only focus on heterocycle anion with N like pyrrole and Pyrazole and pseudohalogen like CN- or SCN-"""
        N_sum = []
        N_jump = []
        cycles = []
        nx_graph = construct_clean_graph(self.structure_graph)
        simple_cycles = nx.simple_cycles(nx_graph, length_bound=16)
        for cycle in simple_cycles:
            cycles.append(cycle)
        for site_index in self.n_indices:
            if site_index not in N_jump:
                cn = self.get_cn(site_index)
                neighbors = self.get_connected_sites(site_index)             
                cm = num_neighbor_metal(neighbors)
                """check for CN- and N3-; SCN- and OCN- has already been checked by previous S_check and O_check"""
                if cn-cm==1:
                    non_metal = non_metal_neighbor(neighbors)
                    if str(non_metal.site.specie) == 'C':
                        neighbor_C = non_metal.index
                        C_neighbors = self.get_connected_sites(neighbor_C)
                        cn_neighbor_C = self.get_cn(neighbor_C)
                        cm_neighbor_C = num_neighbor_metal(C_neighbors)
                        if cn_neighbor_C-cm_neighbor_C == 1:
                            N_sum.append(site_index)
                    if str(non_metal.site.specie) == 'N':
                        neighbor_N = non_metal.index
                        N_neighbors = self.get_connected_sites(neighbor_N)
                        cn_neighbor_N = self.get_cn(neighbor_N)
                        if cn_neighbor_N == 2 and (str(N_neighbors[0].site.specie) == 'N' and str(N_neighbors[1].site.specie) == 'N'):
                            angle = get_angle_between_site_and_neighbors(non_metal.site, N_neighbors)
                            if np.abs(180 - angle) < 25 or np.abs(0 - angle) < 25 :
                                N_sum.append(site_index)
                                N_jump.append(N_neighbors[0].index)
                                N_jump.append(N_neighbors[1].index)
                """check for 5 aromatic rings with N and Porphyrin"""
                if cn-cm==2:
                    non_metals = non_metal_neighbors(neighbors)
                    undercoordinated_nitrogen = _guess_underbound_nitrogen_cn2(
                    self.structure,
                    site_index,
                    non_metals,
                    self.get_connected_sites(non_metals[0].index),
                    self.get_connected_sites(non_metals[1].index),
                    30,
                )
                    if undercoordinated_nitrogen:
                        N_sum.append(site_index)
                    elif not undercoordinated_nitrogen:
                        for cycle in cycles:
                            N_possible_jump = []
                            if len(cycle) == 16 and site_index in cycle and all( not is_metal(self.structure[ring]) for ring in cycle):
                                num_n = 0
                                num_c = 0
                                num_n3 = 0
                                for ring in cycle:
                                    if str(self.structure[ring].specie) == 'N':
                                        N_possible_jump.append(ring)
                                        num_n += 1
                                        n_cn = self.get_cn(ring)
                                        n_cm = num_neighbor_metal(self.get_connected_sites(ring))
                                        if n_cn-n_cm == 3:
                                            num_n3 += 1
                                    if str(self.structure[ring].specie) == 'C':
                                        num_c += 1
                                if (num_n == 4) and (num_n3 == 0) and (num_c == 12):
                                    N_sum.append(site_index)
                                    N_sum.append(site_index)
                                if (num_n == 4) and (num_n3 == 1) and (num_c == 12):
                                    N_sum.append(site_index)
                                if (num_n == 4) and (num_c == 12):
                                    for n in N_possible_jump:
                                        N_jump.append(n)
                                    break
                        """"check Porphyrin first, and if there is a Porphyrin will skip this 5 ring check."""
                        CN3_N = 0
                        if site_index not in N_jump:
                            for cycle in cycles:
                                if len(cycle) == 5 and site_index in cycle and all( not is_metal(self.structure[ring]) for ring in cycle) :
                                    for ring in cycle:
                                        cn = len(non_metal_neighbors(self.get_connected_sites(ring)))
                                        if str(self.structure[ring].specie) == 'N':
                                            N_jump.append(ring)
                                        if str(self.structure[ring].specie) == 'N' and cn == 3:
                                            CN3_N = CN3_N+1
                                        if str(self.structure[ring].specie) == 'O' or str(self.structure[ring].specie) == 'S':
                                            CN3_N = CN3_N+1
                                    if CN3_N == 0 or CN3_N == 2 or CN3_N == 4:
                                        N_sum.append(site_index)                         
        return N_sum

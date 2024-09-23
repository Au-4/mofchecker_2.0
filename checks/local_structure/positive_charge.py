# -*- coding: utf-8 -*-
"""Check the positive charge from linkers like N with CN4 and O with CN3."""
from pymatgen.analysis.graphs import StructureGraph

from mofchecker.types import StructureIStructureType

from .base_coordination_check import BaseCoordinationCheck
from ..utils.get_indices import _is_any_neighbor_metal, get_n_indices, get_o_indices


class Positive_charge_Check(BaseCoordinationCheck):

    def __init__(self, structure: StructureIStructureType, structure_graph: StructureGraph):
        """Initialize a new PositivechargeCheck.

        Args:
            structure (StructureIStructureType): Structure to check.
            structure_graph (StructureGraph): StructureGraph of the structure.
        """
        self.structure = structure
        self.n_indices = get_n_indices(self.structure)
        self.o_indices = get_o_indices(self.structure)
        self.structure_graph = structure_graph

    @property
    def name(self):
        """Return the name of the check."""
        return "Positive charge from linkers"

    @property
    def description(self):
        """Return a description of the check."""
        return "Check the positive charge from the linkers."

    def _run_check(self):
        nitrogen_oxygen = self._get_overcoordinated_nitrogen()
        return len(nitrogen_oxygen) == 0, nitrogen_oxygen

    def _get_overcoordinated_nitrogen(self):
        """Check for all N if CN=4 and all O if CN = 3, ignore metal bonds."""
        overcoordinated_nitrogen_oxygen = []

        for site_index in self.n_indices:
            cn = self.get_cn(site_index)  # pylint:disable=invalid-name
            if cn == 4 and not _is_any_neighbor_metal(self.get_connected_sites(site_index)):
                overcoordinated_nitrogen_oxygen.append(site_index)
        for site_index in self.o_indices:
            cn = self.get_cn(site_index)  # pylint:disable=invalid-name
            if cn == 3 and not _is_any_neighbor_metal(self.get_connected_sites(site_index)):
                overcoordinated_nitrogen_oxygen.append(site_index)
        return overcoordinated_nitrogen_oxygen
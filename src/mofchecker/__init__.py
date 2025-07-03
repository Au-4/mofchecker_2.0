# -*- coding: utf-8 -*-
"""Basic sanity checks for MOFs."""
import os
import warnings
import numpy as np
from collections import OrderedDict
from pathlib import Path
from typing import Iterable, List, Union
from scipy.spatial.transform import Rotation as R
import psutil
import sys
import networkx as nx
from ase import Atoms
from backports.cached_property import cached_property
from pymatgen.analysis.graphs import ConnectedSite, StructureGraph
from pymatgen.core import IStructure, Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from structuregraph_helpers.analysis import get_cn
from structuregraph_helpers.create import construct_clean_graph, get_structure_graph
from structuregraph_helpers.hash import (
    decorated_graph_hash,
    decorated_scaffold_hash,
    undecorated_graph_hash,
    undecorated_scaffold_hash,
)

from mofchecker.checks.local_structure.geometrically_exposed_metal import GeometricallyExposedMetal
from mofchecker.checks.local_structure.undercoordinated_alkaline import (
    UnderCoordinatedAlkaliAlkaline,
)
from mofchecker.checks.local_structure.undercoordinated_rare_earth import (
    UnderCoordinatedRareEarthCheck,
)

from .checks.charge_check import ChargeCheck
#from .checks.oxi_check import OxiCheck
from .checks.floating_solvent import FloatingSolventCheck
from .checks.global_structure import HasCarbon, HasHydrogen, HasMetal, HasNitrogen
from .checks.global_structure.graphcheck import IsThreeDimensional
from .checks.local_structure import (
    AtomicOverlapCheck,
    FalseOxoCheck,
    OverCoordinatedCarbonCheck,
    OverCoordinatedHydrogenCheck,
    OverCoordinatedNitrogenCheck,
    UnderCoordinatedCarbonCheck,
    UnderCoordinatedNitrogenCheck,
    Positive_charge_Check,
    Negative_charge_Check,
    Fusedring_Check,
    O_site_adding_hydrogen,
    adding_linker,
)
from .checks.oms import MOFOMS
from .checks.utils.get_indices import get_c_indices, get_h_indices, get_metal_indices, get_n_indices, get_o_indices,get_ge_indices, get_sb_indices
from .checks.zeopp import PorosityCheck
from .symmetry import get_spacegroup_symbol_and_number, get_symmetry_hash
from .utils import _check_if_ordered
from .version import get_version
from .errors import Too_many_adding_linkers

__version__ = get_version()

__all__ = ["__version__", "MOFChecker", "DESCRIPTORS"]


DESCRIPTORS = [
    "name",
    "linker_name",
    "graph_hash",
    "spacegroup_symbol",
    "spacegroup_number",
    "undecorated_graph_hash",
    "scaffold_hash",
    "undecorated_scaffold_hash",
    "symmetry_hash",
    "formula",
    "path",
    "linker_path",
    "density",
    "has_carbon",
    "has_hydrogen",
    "has_atomic_overlaps",
    "has_overcoordinated_c",
    "has_overcoordinated_n",
    "has_overcoordinated_h",
    "has_undercoordinated_c",
    "has_undercoordinated_n",
    "has_undercoordinated_rare_earth",
    "has_metal",
    "has_lone_molecule",
    "has_high_charges",
    "metal_number",
    "positive_charge_from_linkers",
    "negative_charge_from_linkers",
    "adding_hydrogen",
    "adding_linker",
    "possible_charged_fused_ring",
    "is_porous",
    "has_suspicicious_terminal_oxo",
    "has_undercoordinated_alkali_alkaline",
    "has_geometrically_exposed_metal",
    "has_3d_connected_graph",
]
def check_min_ram_gb(min_gb=4):
    total_gb = psutil.virtual_memory().total / (1024**3)
    if total_gb < min_gb:
        raise MemoryError(f"Need at least {min_gb} GB RAM，now only {total_gb:.2f} GB")
        sys.exit(1)

class MOFChecker:
    """MOFChecker performs basic sanity checks for MOFs."""

    def __init__(
        self,
        structure: Union[Structure, IStructure],
        linker_structure: Union[Structure, IStructure],
        symprec: float = 0.5,
        angle_tolerance: float = 5,
        h_number: int = 0,
        linker_number: int = 0,
        primitive: bool = True,
    ):
        """Construct a MOFChecker instance.

        Args:
            structure (Structure): pymatgen Structure object
            symprec (float): Symmetry precision
            angle_tolerance (float): Angle tolerance
            primitive (bool): If True,
                use primitive cell for structure

        Raises:
            NotImplementedError in the case of partial occupancies
        """
        check_min_ram_gb(4)
        _check_if_ordered(structure)

        if (symprec is not None) or (angle_tolerance is not None):
            try:
                structure = SpacegroupAnalyzer(
                    structure, symprec=symprec, angle_tolerance=angle_tolerance
                ).get_symmetrized_structure()
                
            except TypeError:
                # If symmetrization fails
                pass

        if primitive:
            structure = structure.get_primitive_structure()

        if isinstance(structure, Structure):
            structure = IStructure.from_sites(structure)

        self.structure = structure
        self.h_number = h_number
        self.linker_number = linker_number
        self.linker_structure = linker_structure

        self.metal_indices = get_metal_indices(self.structure)

        self.charges = None
        self._porous = ""
        self.metal_features = None
        self._cnn_method = "vesta"
        self._filename = None
        self._name = None
        self._linker_filename = None
        self._linker_name = None
        self.c_indices = get_c_indices(self.structure)
        self.h_indices = get_h_indices(self.structure)
        self.n_indices = get_n_indices(self.structure)
        self.o_indices = get_o_indices(self.structure)
        self._overvalent_c = None
        self._overvalent_n = None
        self._overvalent_h = None
        self._graph = None
        self._nx_graph = None

        self._connected_sites = {}
        self._cns = {}
        self._checks = {
            "has_c": HasCarbon(self.structure),
            "has_h": HasHydrogen(self.structure),
            "has_metal": HasMetal(self.structure),
            "has_nitrogen": HasNitrogen(self.structure),
            "no_atomic_overlaps": AtomicOverlapCheck(self.structure),
            "no_undercoordinated_carbon": UnderCoordinatedCarbonCheck.from_mofchecker(self),
            "no_overcoordinated_carbon": OverCoordinatedCarbonCheck.from_mofchecker(self),
            "no_overcoordinated_hydrogen": OverCoordinatedHydrogenCheck.from_mofchecker(self),
            "no_overcoordinated_nitrogen": OverCoordinatedNitrogenCheck.from_mofchecker(self),
            "no_undercoordinated_nitrogen": UnderCoordinatedNitrogenCheck.from_mofchecker(self),
            "no_undercoordinated_rare_earth": UnderCoordinatedRareEarthCheck.from_mofchecker(self),
            "no_undercoordinated_alkali_alkaline": UnderCoordinatedAlkaliAlkaline.from_mofchecker(
                self
            ),
            "no_geometrically_exposed_metal": GeometricallyExposedMetal.from_mofchecker(self),
            "no_floating_molecule": FloatingSolventCheck.from_mofchecker(self),
            "no_high_charges": ChargeCheck(self.structure),
#            "metal_oxi": OxiCheck(self.structure),
            "positive_charge": Positive_charge_Check.from_mofchecker(self),
            "negative_charge": Negative_charge_Check.from_mofchecker(self),
            "add_h": O_site_adding_hydrogen.from_mofchecker(self),
            "add_linker": adding_linker.from_mofchecker(self),
            "fused_ring": Fusedring_Check.from_mofchecker(self),
            "is_porous": PorosityCheck(self.structure),
            "no_oms": MOFOMS.from_mofchecker(self),
            "no_false_terminal_oxo": FalseOxoCheck.from_mofchecker(self),
            "has_3d_connected_graph": IsThreeDimensional.from_mofchecker(self),
        }

    @property
    def checks(self):
        """Get a dictionary of all check classes."""
        return self._checks

    def _set_filename(self, path):
        self._filename = os.path.abspath(path)
        self._name = Path(path).stem
    def _set_linkername(self, linker_path):
        self._linker_filename = os.path.abspath(linker_path)
        self._linker_name = Path(linker_path).stem
    def get_overlapping_indices(self):
        """Return the indices of overlapping atoms."""
        return self.checks["no_atomic_overlaps"].flagged_indices

    @property
    def graph_hash(self) -> str:
        """Return the Weisfeiler-Lehman graph hash.

        Hashes are identical for isomorphic graphs
        (taking the atomic kinds into account)
        and there are guarantees that non-isomorphic graphs will get different hashes.

        Returns:
            str: Graph hash
        """
        return decorated_graph_hash(self._graph, lqg=False)

    @property
    def spacegroup_symbol(self) -> str:
        """Return the international spacegroup symbol."""
        return get_spacegroup_symbol_and_number(self.structure)["symbol"]

    @property
    def spacegroup_number(self) -> int:
        """Return the international spacegroup number."""
        return get_spacegroup_symbol_and_number(self.structure)["number"]

    @cached_property
    def symmetry_hash(self) -> str:
        """Hash the structure based on its symmetrized versions.

        That is, the spacegroup and Wyckoff letters.

        Returns:
            str: Symmetry hash
        """
        return get_symmetry_hash(self.structure)

    @property
    def undercoordinated_c_candidate_positions(self) -> Iterable[Iterable[float]]:
        """Candidate positions for addition H on undercoordinated C.

        Returns:
            Iterable[Iterable[float]]: Candidate positions
        """
        return self.checks["no_undercoordinated_carbon"].candidate_positions

    @property
    def undercoordinated_n_candidate_positions(self) -> Iterable[Iterable[float]]:
        """Candidate positions for addition H on undercoordinated N.

        Returns:
            Iterable[Iterable[float]]: Candidate positions
        """
        return self.checks["no_undercoordinated_nitrogen"].candidate_positions

    @property
    def undecorated_graph_hash(self) -> str:
        """Return the Weisfeiler-Lehman graph hash.

        Undecorated means that the atomic kinds are not taken into account.
        Hashes are identical for isomorphic graphs and there are
        guarantees that non-isomorphic graphs will get different hashes.

        Returns:
            str: Graph hash without atomic kinds
        """
        return undecorated_graph_hash(self._graph, lqg=False)

    @property
    def scaffold_hash(self) -> str:
        """Return the Weisfeiler-Lehman graph hash for the scaffold.

        The scaffold is the graph with the all terminal groups and
        atoms removed (i.e., formally, bridges are broken).
        Hashes are identical for isomorphic graphs and there are
        guarantees that non-isomorphic graphs will get different hashes.

        Returns:
            str: Graph hash for the scaffold
        """
        return decorated_scaffold_hash(self._graph, lqg=False)

    @property
    def undecorated_scaffold_hash(self) -> str:
        """Return the Weisfeiler-Lehman graph hash for the undecorated scaffold.

        The scaffold is the graph with the all terminal groups and
        atoms removed (i.e., formally, bridges are broken).
        Undecorated means that the atomic numbers are not taken into account.

        Hashes are identical for isomorphic graphs and there are
        guarantees that non-isomorphic graphs will get different hashes.

        Returns:
            str: Graph hash for the undecorated scaffold
        """
        return undecorated_scaffold_hash(self._graph, lqg=False)

    @property
    def has_atomic_overlaps(self) -> bool:
        """Check if there are any overlaps in the structure."""
        return not self.checks["no_atomic_overlaps"].is_ok

    @property
    def name(self) -> str:
        """Return filename stem if the MOFChecker instance was created based on a file."""
        return self._name
    @property
    def linker_name(self) -> str:
        """Return filename stem if the MOFChecker instance was created based on a file."""
        return self._linker_name
    @property
    def path(self) -> str:
        """Return filepath if created from a file."""
        return self._filename
    @property
    def linker_path(self) -> str:
        """Return filepath if created from a file."""
        return self._linker_filename
    @property
    def has_carbon(self) -> bool:
        """Check if there is any carbon atom in the structure."""
        return self.checks["has_c"].is_ok

    @property
    def has_hydrogen(self) -> bool:
        """Check if there is any hydrogen atom in the structure."""
        return self.checks["has_h"].is_ok

    @property
    def density(self) -> float:
        """Density of structure."""
        return self.structure.density

    @property
    def volume(self) -> float:
        """Volume of structure in A^3."""
        return self.structure.volume

    @property
    def formula(self) -> str:
        """Return the chemical formula of the structure."""
        return self.structure.formula

    @property
    def has_3d_connected_graph(self) -> bool:
        """Check if the graph is 3D connected."""
        return self.checks["has_3d_connected_graph"].is_ok

    # ToDo: deprecate one of overvalent/overcoordinated
    @property
    def has_overvalent_c(self) -> bool:
        """Return true ifsome carbon in the structure has more than 4 neighbors.

        Returns:
            bool: True if carbon with CN > 4 in structure.
        """
        return not self.checks["no_overcoordinated_carbon"].is_ok

    @property
    def has_overcoordinated_c(self) -> bool:
        """Return true ifsome carbon in the structure has more than 4 neighbors.

        Alias for has_overvalent_c.

        Returns:
            bool: True if carbon with CN > 4 in structure.
        """
        return self.has_overvalent_c

    @property
    def overvalent_c_indices(self) -> List[int]:
        """Return indices of carbons with more than 4 neighbors.

        Returns:
            List[int]: Indices of carbons with CN > 4.
        """
        return self.checks["no_overcoordinated_carbon"].flagged_indices

    @property
    def has_overvalent_h(self) -> bool:
        """Return true if some hydrogen has more than 1 neighbor.

        Returns:
            bool: True if hydrogen with CN > 1 in structure.
        """
        return not self.checks["no_overcoordinated_hydrogen"].is_ok

    @property
    def has_overcoordinated_h(self) -> bool:
        """See has_overvalent_h."""
        return self.has_overvalent_h

    @property
    def overvalent_h_indices(self) -> List[int]:
        """Return indices of hydrogens with more than 1 neighbors.

        Returns:
            List[int]: Indices of hydrogens with CN > 1.
        """
        return self.checks["no_overcoordinated_hydrogen"].flagged_indices

    @property
    def has_undercoordinated_c(self) -> bool:
        """Check if there is a carbon that likely misses hydrogen."""
        return not self.checks["no_undercoordinated_carbon"].is_ok

    @property
    def undercoordinated_c_indices(self) -> List[int]:
        """Return indices of carbon in the structure that likely miss some neighbors.

        Returns:
            List[int]: Indices of carbons with CN < 4.
        """
        return self.checks["no_undercoordinated_carbon"].flagged_indices

    @property
    def has_undercoordinated_n(self) -> bool:
        """Return False if there is a nitrogen that likely misses hydrogen."""
        return not self.checks["no_undercoordinated_nitrogen"].is_ok

    @property
    def undercoordinated_n_indices(self) -> List[int]:
        """Return indices of nitrogen that likely miss some neighbors.

        Returns:
            List[int]: Indices of nitrogens with CN < 4.
        """
        return self.checks["no_undercoordinated_nitrogen"].flagged_indices

    @property
    def has_undercoordinated_rare_earth(self) -> bool:
        """Return True if there is a rare earth metal that likely misses hydrogen."""
        return not self.checks["no_undercoordinated_rare_earth"].is_ok

    @property
    def undercoordinated_rare_earth_indices(self) -> List[int]:
        """Return indices of rare earth metals in the structure that likely miss some neighbors."""
        return self.checks["no_undercoordinated_rare_earth"].flagged_indices

    @property
    def has_undercoordinated_alkali_alkaline(self) -> bool:
        """Return True if there is a alkali or alkaline earth metal that likely misses some neighbors."""
        return not self.checks["no_undercoordinated_alkali_alkaline"].is_ok

    @property
    def has_geometrically_exposed_metal(self) -> bool:
        """Check if there is a metal that is geometrically exposed."""
        return not self.checks["no_geometrically_exposed_metal"].is_ok

    @property
    def geometrically_exposed_metal_indice(self) -> List[int]:
        """Check if there is a metal that is geometrically exposed."""
        return self.checks["no_geometrically_exposed_metal"].flagged_indices

    @property
    def is_porous(self) -> Union[bool, None]:
        """Return True if the MOF is porous according to the CoRE-MOF definition.

        Returns None if the check could not be run successfully.

        Returns:
            Union[bool, None]: True if porous, False if not porous, None if check could not be run.
        """
        return self.checks["is_porous"].is_ok

    @property
    def has_high_charges(self) -> Union[bool, None]:
        """Check if the structure has unreasonably high EqEq charges.

        Returns None if the check could not be run successfully.

        Returns:
            Union[bool, None]: True if charges are too high,
                False if charges are ok, None if check could not be run.
        """
        return not self.checks["no_high_charges"].is_ok
    
 #   @property
 #   def metal_oxidation_states(self) -> int:
 #       """Check the oxidationstate.
#
#        Returns:
#            the sum of metal oxidation state.
#        """
#        return self.checks["metal_oxi"]

    @property
    def has_suspicicious_terminal_oxo(self) -> bool:
        """Flag metals with a potentially wrong terminal oxo group."""
        return not self.checks["no_false_terminal_oxo"].is_ok

    @property
    def suspicicious_terminal_oxo_indices(self) -> List[int]:
        """Return indices of metals with a potentially wrong terminal oxo group."""
        return self.checks["no_false_terminal_oxo"].flagged_indices

    @property
    def nx_graph(self) -> nx.Graph:
        """Return a networkx graph with atom numbers as node labels."""
        if self._nx_graph is None:
            _ = self.graph
        return self._nx_graph

    @property
    def graph(self) -> StructureGraph:
        """Return a pymatgen structure graph."""
        if self._graph is None:
            self._graph = get_structure_graph(self.structure, self._cnn_method)
            self._nx_graph = construct_clean_graph(self._graph)
        return self._graph

    def get_connected_sites(self, site_index: int) -> List[ConnectedSite]:
        """Get connected sites for given index.

        Uses internal cache for speedup.

        Args:
            site_index (int): Index of the site to get connected sites for.

        Returns:
            List[ConnectedSite]: List of connected sites.
        """
        if site_index not in self._connected_sites:
            self._connected_sites[site_index] = self.graph.get_connected_sites(site_index)
        return self._connected_sites[site_index]

    def get_cn(self, site_index: int) -> int:
        """Get coordination number for site.

        Uses internal cache for speedup.

        Args:
            site_index (int): index of site in pymatgen Structure

        Returns:
            int: Coordination number
        """
        if site_index not in self._cns:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self._cns[site_index] = get_cn(self.graph, site_index)
        return self._cns[site_index]

    @property
    def has_overvalent_n(self) -> bool:
        """Return True if some nitrogen has more than 4 neighbors.

        Returns:
            bool: True if nitrogen with CN > 4 in structure.
        """
        return not self.checks["no_overcoordinated_nitrogen"].is_ok

    @property
    def positive_charge_from_linkers(self) -> int:
        """Return positive charge value"""
        return len(self.checks["positive_charge"].flagged_indices)
    @property
    def negative_charge_from_linkers(self) -> int:
        """Return negative charge value"""
        return len(self.checks["negative_charge"].flagged_indices)
    @property
    def adding_hydrogen(self) -> Iterable[Iterable[float]]:
        """Return adding H coordinates cif """
        h_coord = []
        new_structure = Structure.from_sites(self.structure.sites)
        atom_type = "H"
        n = 0
        h_number = self.h_number
  #      print(self.checks["no_undercoordinated_carbon"].is_ok, self.checks["no_undercoordinated_nitrogen"].is_ok)
        if not (len(self.checks["no_undercoordinated_carbon"].candidate_positions) == 0):
            h_coord.append(self.checks["no_undercoordinated_carbon"].candidate_positions)
    #        for group in h_coord:
    #            for vec in group:
    #                fra_vec = new_structure.lattice.get_fractional_coords(vec)
    #                new_structure.append(atom_type, fra_vec)
    #    h_coord = []
        if not (len(self.checks["no_undercoordinated_nitrogen"].candidate_positions) == 0):
            h_coord.append(self.checks["no_undercoordinated_nitrogen"].candidate_positions)
        h_coord.append(self.checks["add_h"].candidate_positions)
        for group in h_coord:
            for vec in group:
                if n < h_number:
                    fra_vec = new_structure.lattice.get_fractional_coords(vec)
                    new_structure.append(atom_type, fra_vec)
                    n += 1
                else:
                    break
        new_name = str(self.name)+'_add_H.cif'
        new_structure.to(new_name,"cif")
        return new_structure

    @property
    def adding_linker(self) -> Iterable[Iterable[float]]:
        """Return adding linker coordinates cif """
        def normal_vector(a,b,c):
            ab = b-a
            ac = c-a
            normal = np.cross(ab, ac)
            return normal
        
        def distance_between_groups(group1, structure):
            min_dist = np.inf
            for atom1 in group1:
                fra_coords = structure.lattice.get_fractional_coords(atom1)
                for site in structure:
                    site_fra_coords = structure.lattice.get_fractional_coords(site.coords)
                    dist = structure.lattice.get_distance_and_image(fra_coords, site_fra_coords)[0]
                    if dist < min_dist:
                        min_dist = dist
            return min_dist
        
        def rotate_group(coords, axis, angle):
            rotation = R.from_rotvec(angle * np.array(axis))
            return rotation.apply(coords)
        
        i=0
        linker_number = self.linker_number
        mof_coords = []
        new_structure = Structure.from_sites(self.structure.sites)
        for site in self.structure.sites:
            mof_coords.append(site.coords)
        z_sites = self.checks["add_linker"].candidate_positions
        x_indices = self.checks["add_linker"].flagged_indices
        new_name = str(self.name)+'_add_'+str(self.linker_name)+'.cif'
        if linker_number > len(z_sites):
            raise Too_many_adding_linkers('Too_many_adding_linkers, should check manully')
        for z_site in z_sites:
            if i < linker_number:
                x_index = x_indices[i]
                x_coords=new_structure[x_index].coords
                i+=1
                mindis = 5
                mindis_coords = [0,0,0]
                linker_coords = [z_site]
                angle = 0
                for linker_site in self.linker_structure.sites:
                    if str(linker_site.label) == 'X':
                        vec = z_site - linker_site.coords
                        linker_site.coords = linker_site.coords + vec
                if len(self.linker_structure.sites)>1:
                    for linker_site in self.linker_structure.sites:
                        if str(linker_site.label) != 'X':
                            linker_site.coords = linker_site.coords + vec 
                            dis = np.linalg.norm(linker_site.coords-z_site)
                            if dis < mindis:
                                mindis = dis
                                mindis_coords = linker_site.coords
                            linker_coords.append(linker_site.coords)
                    axis = normal_vector(z_site, x_coords,mindis_coords)
                    step_angle = np.radians(5)
                    origin = z_site
                    linker_centered = linker_coords - origin
                    dist_data={}
                    while angle < np.radians(360):
                        linker_centered = rotate_group(linker_centered, axis, step_angle) 
                        linker_rotated = linker_centered + origin
                        min_dist = distance_between_groups(linker_rotated, new_structure)
                        dist_data[min_dist] = angle
                        angle += step_angle
                        if min_dist > 2.4:
                            break
                    max_dist = max(dist_data.keys())
                    angle = dist_data[max_dist]
                    linker_centered = linker_coords - origin
                    linker_centered = rotate_group(linker_centered, axis, angle) 
                    linker_rotated = linker_centered + origin

                    j=1
                    for linker_site in self.linker_structure.sites:
                        if str(linker_site.label) != 'X':
                            linker_site.coords = linker_rotated[j]
                            j+=1
            
                for linker_site in self.linker_structure.sites:
                    atom_type = linker_site.species_string
                    fra_coords = new_structure.lattice.get_fractional_coords(linker_site.coords)
                    new_structure.append(atom_type, fra_coords)
        new_structure.to(new_name,"cif")
        return new_structure

    @property
    def metal_number(self) -> int:
        """Return metal number in primitive cell, just a simple test to compare with oximachine"""

        metal_num = len(get_metal_indices(self.structure))-len(get_sb_indices(self.structure))
        return metal_num
    @property
    def possible_charged_fused_ring(self) -> bool:
        """Return if there is a complicated possible charged fused ring, have to check mannully"""
        return not self.checks["fused_ring"].is_ok
    @property
    def has_overcoordinated_n(self) -> bool:
        """Return True if some nitrogen has more than 4 neighbors.

        Alias for has_overvalent_n.

        Returns:
            bool: True if nitrogen with CN > 4 in structure.
        """
        return self.has_overvalent_n

    @property
    def has_lone_molecule(self) -> bool:
        """Return true if there is a isolated floating atom or molecule."""
        return not self.checks["no_floating_molecule"].is_ok

    @property
    def lone_molecule_indices(self) -> List[int]:
        """Return indices of non-periodic connected component in the structure."""
        return self.checks["no_floating_molecule"].flagged_indices

    @classmethod
    def _from_file(cls, path: str, linker_path:str, **kwargs):
        structure = Structure.from_file(path)
        mofchecker = cls(structure, **kwargs)
        mofchecker._set_filename(path)
        mofchecker._set_linkername(linker_path)  # pylint:disable=protected-access
        return mofchecker

    @classmethod
    def from_cif(
        cls,
        path: Union[str, Path],
        linker_path: Union[str, Path]= None,
        symprec: float = 0.5,
        angle_tolerance: float = 5,
        h_number: int = 0,
        linker_number: int = 0,
        primitive: bool = True,
    ) -> "MOFChecker":
        """Create a MOFChecker instance from a CIF file.

        Args:
            path (Union[str, Path]): Path to string file
            symprec (float): Symmetry tolerance
            angle_tolerance (float): Angle tolerance
            primitive (bool): Whether to use primitive cell

        Returns:
            MOFChecker: Instance of MOFChecker
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cifparser = CifParser(path)
            structure = cifparser.get_structures()[0]
            linker_structure = None
            if linker_path != None:
                linker_cifparser = CifParser(linker_path)
                linker_structure = linker_cifparser.get_structures()[0]
                
            omscls = cls(
                structure, linker_structure, symprec=symprec, angle_tolerance=angle_tolerance, h_number=h_number, linker_number = linker_number, primitive=primitive
            )
            omscls._set_filename(path)
            if linker_path != None:
                omscls._set_linkername(linker_path) 
            return omscls

    @classmethod
    def from_ase(
        cls, atoms: Atoms, linker_atoms: Optional[Atoms] = None, symprec: float = 0.5, angle_tolerance: float = 5, primitive: bool = True
    ) -> "MOFChecker":
        """Create a MOFChecker instance from an ASE atoms object.

        Args:
            atoms (Atoms): ase atoms object
            linker_atoms (Atoms): ase atoms object for linker
            symprec (float): Symmetry tolerance
            angle_tolerance (float): Angle tolerance
            primitive (bool): Whether to use primitive cell

        Returns:
            MOFChecker: Instance of MOFChecker
        """
        adaptor = AseAtomsAdaptor()
        structure = adaptor.get_structure(atoms)
        linker_structure = None
        if linker_atoms is not None:
            linker_structure = adaptor.get_structure(linker_atoms)
        omscls = cls(
            structure, linker_structure, symprec=symprec, angle_tolerance=angle_tolerance, primitive=primitive
        )
        return omscls

    @property
    def has_metal(self) -> bool:
        """Return True if the structure has a metal."""
        return self.checks["has_metal"].is_ok

    @property
    def has_oms(self) -> bool:
        """Return true if open metal sites are detected."""
        return not self.checks["no_oms"].is_ok
    @property
    def oms_indice(self) -> List[int]:
        """Return true if open metal sites are detected."""
        return self.checks["no_oms"].flagged_indices
    
    def _set_cnn(self, method="vesta"):
        if self._cnn_method == method.lower():
            return
        self._cnn_method = method.lower()

    def get_mof_descriptors(self, descriptors=None) -> OrderedDict:
        """Run sanity checks and get a dictionary with the result.

        Args:
            descriptors (List): If provided, compute only the passed descriptors

        Returns:
            OrderedDict: result of overall checks
        """
        if descriptors is None:
            descriptors = DESCRIPTORS

        result_dict = OrderedDict(
            ((descriptor, getattr(self, descriptor)) for descriptor in descriptors)
        )
        return result_dict

# -*- coding: utf-8 -*-
"""Utility function for getting the indices for certain atoms in the structure."""
import functools
from typing import Union

import pymatgen
from pymatgen.core import IStructure, Structure

from ..data import _get_vdw_radius
from ...definitions import METALS


def _vdw_radius_neighbors(structure, site_index, tolerance: float = 1.5):
    elem = str(structure[site_index].specie)
    radius = _get_vdw_radius(elem)
    return structure.get_neighbors(structure[site_index], tolerance * radius)


def is_metal(site: pymatgen.core.Site) -> bool:
    """Return True if the site is a metal.

    Considers transition metal, lanthanide, actinide,
    or Al, Ga, In, Tl, Ge, Sn, Pb, Sb, Bi, Po

    Args:
        site: pymatgen.core.Site

    Returns:
        bool: True if the site is a metal
    """
    if str(site.specie) in METALS:
        return True
    return False
def is_halogen(site: pymatgen.core.Site) -> bool:
    if str(site.specie) == 'F' or str(site.specie) == 'Cl' or str(site.specie) == 'Br' or str(site.specie) == 'I':
        return True
    return False


def is_boron(site: pymatgen.core.Site) -> bool:
    if str(site.specie) == 'B':
        return True
    return False
def is_H(site: pymatgen.core.Site) -> bool:
    if str(site.specie) == 'H':
        return True
    return False
def is_C(site: pymatgen.core.Site) -> bool:
    if str(site.specie) == 'C':
        return True
    return False
def is_N(site: pymatgen.core.Site) -> bool:
    if str(site.specie) == 'N':
        return True
    return False
def is_S(site: pymatgen.core.Site) -> bool:
    if str(site.specie) == 'S':
        return True
    return False
def is_P(site: pymatgen.core.Site) -> bool:
    if str(site.specie) == 'P':
        return True
    return False
def is_O(site: pymatgen.core.Site) -> bool:
    if str(site.specie) == 'O':
        return True
    return False
@functools.lru_cache(maxsize=2, typed=False)
def _get_indices(immutable_structure: IStructure) -> dict:
    return {
        "c": _get_c_indices(immutable_structure),
        "h": _get_h_indices(immutable_structure),
        "n": _get_n_indices(immutable_structure),
        "s": _get_s_indices(immutable_structure),
        "o": _get_o_indices(immutable_structure),
        "metal": _get_metal_indices(immutable_structure),
        "rare_earth": _get_rare_earth_indices(immutable_structure),
        "alkali_alkaline": _get_alkali_alkaline_indices(immutable_structure),
        "halogen": _get_halogen_indices(immutable_structure),
        "all": _get_all_indices(immutable_structure),
        "x":_get_x_indices(immutable_structure),
        "x2":_get_x2_indices(immutable_structure),
        "ge":_get_ge_indices(immutable_structure),
        "sb":_get_sb_indices(immutable_structure),
    }


def _get_c_indices(structure):
    return [i for i, species in enumerate(structure.species) if str(species) == "C"]
def _get_x_indices(structure):
    return [i for i, site in enumerate(structure.sites) if str(site.label) == "X"]
def _get_x2_indices(structure):
    return [i for i, site in enumerate(structure.sites) if str(site.label) == "X2"]
def _get_all_indices(structure):
    return [i for i, species in enumerate(structure.species) if str(species) != "M"]

def _get_h_indices(structure):
    return [i for i, species in enumerate(structure.species) if str(species) == "H"]


def _get_n_indices(structure):
    return [i for i, species in enumerate(structure.species) if str(species) == "N"]

def _get_s_indices(structure):
    return [i for i, species in enumerate(structure.species) if str(species) == "S"]

def _get_o_indices(structure):
    return [i for i, site in enumerate(structure) if str(site.specie) == "O"]


def _get_metal_indices(structure):
    return [i for i, site in enumerate(structure) if is_metal(site)]
def _get_ge_indices(structure):
    return [i for i, species in enumerate(structure.species) if str(species) == "Ge"]
def _get_sb_indices(structure):
    return [i for i, species in enumerate(structure.species) if str(species) == "Sb"]
def _get_rare_earth_indices(structure):
    return [i for i, site in enumerate(structure) if site.specie.is_rare_earth]


def _get_alkali_alkaline_indices(structure):
    return [
        i for i, site in enumerate(structure) if site.specie.is_alkali or site.specie.is_alkaline
    ]


def _get_halogen_indices(structure):
    return [i for i, site in enumerate(structure) if site.specie.is_halogen]


def get_h_indices(structure):
    """Get the indices of all H."""
    return get_indices(structure)["h"]
def get_x_indices(structure):
    """Get the indices of all X."""
    return get_indices(structure)["x"]
def get_x2_indices(structure):
    """Get the indices of all X2."""
    return get_indices(structure)["x2"]

def get_c_indices(structure):
    """Get the indices of all C."""
    return get_indices(structure)["c"]
def get_all_indices(structure):
    """Get the indices of all C."""
    return get_indices(structure)["all"]

def get_s_indices(structure):
    """Get the indices of all S."""
    return get_indices(structure)["s"]

def get_n_indices(structure):
    """Get the indices of all N."""
    return get_indices(structure)["n"]
def get_ge_indices(structure):
    """Get the indices of all Ge."""
    return get_indices(structure)["ge"]
def get_sb_indices(structure):
    """Get the indices of all Sb."""
    return get_indices(structure)["sb"]
def get_o_indices(structure):
    """Get the indices of all O."""
    return get_indices(structure)["o"]


def get_metal_indices(structure):
    """Get the indices of all metals."""
    return get_indices(structure)["metal"]


def get_rare_earth_indices(structure):
    """Get the indices of all rare-earth metals."""
    return get_indices(structure)["rare_earth"]


def get_halogen_indices(structure):
    """Get the indices of all halogens."""
    return get_indices(structure)["halogen"]


def get_alkali_alkaline_indices(structure):
    """Get the indices of all alkali and alkaline earth metals."""
    return get_indices(structure)["alkali_alkaline"]


def get_indices(structure: Union[Structure, IStructure]) -> dict:
    """Get all the relevant indices."""
    if isinstance(structure, Structure):
        structure = IStructure.from_sites(structure)
    return _get_indices(structure)


def _is_any_neighbor_metal(neighbors):
    return any(is_metal(neighbor.site) for neighbor in neighbors)

def num_neighbor_metal(neighbors):
    num = 0
    for neighbor in neighbors:
        if is_metal(neighbor.site):
            num = num+1
    return int(num)

def num_neighbor_H(neighbors):
    num = 0
    for neighbor in neighbors:
        if is_H(neighbor.site):
            num = num+1
    return int(num)
def num_neighbor_O(neighbors):
    num = 0
    for neighbor in neighbors:
        if is_O(neighbor.site):
            num = num+1
    return int(num)
def num_neighbor_S(neighbors):
    num = 0
    for neighbor in neighbors:
        if is_S(neighbor.site):
            num = num+1
    return int(num)
def num_neighbor_N(neighbors):
    num = 0
    for neighbor in neighbors:
        if is_N(neighbor.site):
            num = num+1
    return int(num)

def num_neighbor_halogen(neighbors):
    num = 0
    for neighbor in neighbors:
        if is_halogen(neighbor.site):
            num = num+1
    return int(num)

def non_metal_neighbor(neighbors):
    for neighbor in neighbors:
        if not is_metal(neighbor.site):
            return neighbor
def non_metal_neighbors(neighbors):
    non_metals = []
    for neighbor in neighbors:
        if not is_metal(neighbor.site):
            non_metals.append(neighbor)
    return non_metals

def non_H_neighbors(neighbors):
    non_H = []
    for neighbor in neighbors:
        if not is_H(neighbor.site):
            non_H.append(neighbor)
    return non_H

def _is_any_neighbor_boron(neighbors):
    return any(is_boron(neighbor.site) for neighbor in neighbors)

def _is_any_neighbor_H(neighbors):
    return any(is_H(neighbor.site) for neighbor in neighbors)

def _is_any_neighbor_C(neighbors):
    return any(is_C(neighbor.site) for neighbor in neighbors)

def _is_any_neighbor_N(neighbors):
    return any(is_N(neighbor.site) for neighbor in neighbors)

def _is_any_neighbor_S(neighbors):
    return any(is_S(neighbor.site) for neighbor in neighbors)

def _is_any_neighbor_P(neighbors):
    return any(is_P(neighbor.site) for neighbor in neighbors)

def _is_all_neighbor_O(neighbors):
    return all(is_O(neighbor.site) for neighbor in neighbors)
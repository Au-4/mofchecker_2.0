# -*- coding: utf-8 -*-
"""Check the oxidation states of the structure"""
import warnings
from tempfile import NamedTemporaryFile

import numpy as np
from oximachinerunner import OximachineRunner
from .check_base import AbstractCheck
from ..types import StructureIStructureType


class OxiCheck(AbstractCheck):
    """Check the oxidation states of the structure"""

    def __init__(self, structure: StructureIStructureType):
        """Create an oxidation state check instance.

        Args:
            structure (StructureIStructureType): The structure to check.
        """
        self.structure = structure


    @property
    def name(self) -> str:
        """Return the name of the check."""
        return "Sum of Oxidation state"

    @property
    def description(self) -> str:
        """Return a description of the check."""
        return f"Check the oxidation states of the structure"

    def _run_check(self):
        try:

            with NamedTemporaryFile("w", suffix=".cif") as file:
                runner = OximachineRunner()
                self.structure.to(fmt="cif", filename=file.name)
                oxidation = runner.run_oximachine(file.name, verbose=False)
                oxidationlist = oxidation['prediction']
                metal_charges = np.sum(oxidationlist)

            return metal_charges
        except ImportError:
            warnings.warn("Install the oximachine to run the oxidation check")
            return None
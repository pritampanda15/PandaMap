"""
PandaMap: Protein AND ligAnd interaction MAPper

A Python package for visualizing protein-ligand interactions
with 2D ligand structure representation and minimal external dependencies.

Importing the package prints nothing.  Progress messages go to the
``pandamap`` logger; to see them from Python:

    import logging
    logging.basicConfig(level=logging.INFO, format="%(message)s")

At most once a day a background thread asks PyPI whether a newer PandaMap
exists and, if so, prints a one-line notice asking you to upgrade.  Set
``PANDAMAP_NO_UPDATE_CHECK=1`` to turn that off.
"""

import logging
import warnings
from importlib.metadata import version

from .core import (SimpleLigandStructure, HybridProtLigMapper,
                   analyze_trajectory, list_ligands)
from .render_3d import render_contacts_3d
from .batch import analyze_batch

from .utils.versioning import check_for_updates_in_background

logging.getLogger("pandamap").addHandler(logging.NullHandler())

# Tell people on an old release that a newer one exists (non-blocking)
check_for_updates_in_background()

__all__ = ["SimpleLigandStructure", "HybridProtLigMapper", "analyze_trajectory",
           "analyze_batch", "list_ligands", "render_contacts_3d"]
__version__ = "4.4.0"  # fallback when the package metadata is unavailable
try:
    __version__ = version("pandamap")
except Exception:
    warnings.warn(
        "Failed to retrieve version from importlib.metadata. "
        "Using fallback version: {}".format(__version__),
        UserWarning,
        stacklevel=2
    )

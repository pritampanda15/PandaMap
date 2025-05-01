"""
PandaMap: Protein AND ligAnd interaction MAPper

A Python package for visualizing protein-ligand interactions 
with 2D ligand structure representation and minimal external dependencies.
"""

import warnings
from importlib.metadata import version
import threading
from .utils.versioning import check_for_updates

# Run the check in a background thread
threading.Thread(target=check_for_updates, daemon=True).start()


# Core imports
from .core import SimpleLigandStructure, HybridProtLigMapper

__all__ = ["SimpleLigandStructure", "HybridProtLigMapper"]
__version__ = "4.1.0"  # Keep this as fallback if importlib fails
try:
    __version__ = version("pandamap")
except Exception:
    warnings.warn(
        "Failed to retrieve version from importlib.metadata. "
        "Using fallback version: {}".format(__version__),
        UserWarning,
        stacklevel=2
    )
# Ensure the version is set correctly
# Check for updates in a background thread
# This is to avoid blocking the main thread during import
"""
PandaMap: Protein AND ligAnd interaction MAPper

A Python package for visualizing protein-ligand interactions 
with 2D ligand structure representation and minimal external dependencies.
"""

import warnings
from importlib.metadata import version
import threading

# Core imports
from .core import SimpleLigandStructure, HybridProtLigMapper

__all__ = ["SimpleLigandStructure", "HybridProtLigMapper"]
__version__ = "3.6"  # Keep this as fallback if importlib fails

# --- Auto-update checker (non-blocking) ---
def _check_for_updates():
    """Check PyPI for newer versions without blocking imports."""
    try:
        import requests  # Lazy import to avoid adding dependency unless needed
        
        package_name = "pandamap"  # PyPI package name (adjust if different)
        current_version = version(package_name)  # Gets installed version
        
        # Fetch latest version from PyPI
        response = requests.get(
            f"https://pypi.org/pypi/{package_name}/json",
            timeout=2  # Fail quickly if PyPI is slow/unreachable
        )
        latest_version = response.json()["info"]["version"]
        
        if current_version != latest_version:
            warnings.warn(
                f"PandaMap {latest_version} is available (you have {current_version}). "
                f"Run `pip install --upgrade {package_name}` to update.\n"
                "To disable this check, set env var PANDAMAP_NO_UPDATE_CHECK=1.",
                UserWarning,
                stacklevel=2
            )
    except Exception:
        pass  # Silently fail on any error (network, PyPI down, etc.)

# Run check only if not disabled
if not __import__('os').getenv("PANDAMAP_NO_UPDATE_CHECK"):
    threading.Thread(target=_check_for_updates, daemon=True).start()
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
__version__ = "4.0.0"  # Keep this as fallback if importlib fails

# --- Auto-update checker (non-blocking) ---
def _check_for_updates():
    """Check PyPI for newer versions and notify user with red box."""
    try:
        import requests  # Lazy import
        import importlib.util
        import shutil

        package_name = "pandamap"
        current_version = version(package_name)

        response = requests.get(
            f"https://pypi.org/pypi/{package_name}/json",
            timeout=2
        )
        latest_version = response.json()["info"]["version"]

        if current_version != latest_version:
            message = (
                f"\n🚨 [bold red]PandaMap {latest_version} is available![/bold red] "
                f"[dim](you have {current_version})[/dim]\n\n"
                f"[yellow]Update with:[/yellow] [green]pip install --upgrade {package_name}[/green]\n"
                f"[dim]To disable update checks, set: PANDAMAP_NO_UPDATE_CHECK=1[/dim]\n"
            )

            # Use rich if available and stdout is a terminal
            if shutil.which("rich"):
                try:
                    from rich.console import Console
                    console = Console()
                    console.print(message)
                except ImportError:
                    warnings.warn(message, UserWarning, stacklevel=2)
            else:
                warnings.warn(message, UserWarning, stacklevel=2)

    except Exception:
        pass  # Don't crash anything if this fails

# Run check only if not disabled
if not __import__('os').getenv("PANDAMAP_NO_UPDATE_CHECK"):
    threading.Thread(target=_check_for_updates, daemon=True).start()

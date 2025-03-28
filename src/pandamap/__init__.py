"""
PandaMap: Protein AND ligAnd interaction MAPper

A Python package for visualizing protein-ligand interactions 
with 2D ligand structure representation and minimal external dependencies.
"""

__version__ = "0.1.0"

from .core import SimpleLigandStructure, HybridProtLigMapper

__all__ = ["SimpleLigandStructure", "HybridProtLigMapper"]

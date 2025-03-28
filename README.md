# PandaMap: A Python package for visualizing protein-ligand interactions with 2D ligand structure representation. 

**P**rotein **AND** lig**A**nd interaction **MAP**per: A Python package for visualizing protein-ligand interactions with 2D ligand structure representation

## Overview

PandaMap is a lightweight tool for visualizing protein-ligand interactions from PDB files. It generates intuitive 2D interaction diagrams that display both the ligand structure and its interactions with protein residues.

Key features:
- Visualization of protein-ligand interactions with minimal dependencies
- 2D representation of ligand structure without requiring RDKit
- Detection of multiple interaction types (hydrogen bonds, π-stacking, hydrophobic)
- Command-line interface for quick analysis
- Python API for integration into computational workflows

## Installation

```bash
pip install pandamap
```

## Dependencies

- NumPy
- Matplotlib
- BioPython

## Basic Usage

### Command Line Interface

```bash
# Basic usage
pandamap protein_ligand.pdb --output interactions.png

# Specify a particular ligand by residue name
pandamap protein_ligand.pdb --ligand LIG
```

### Python API

```python
from pandamap import HybridProtLigMapper

# Initialize with PDB file
mapper = HybridProtLigMapper("protein_ligand.pdb", ligand_resname="LIG")

# Run analysis and generate visualization
output_file = mapper.run_analysis(output_file="interactions.png")

# Or run steps separately
mapper.detect_interactions()
mapper.estimate_solvent_accessibility()
mapper.visualize(output_file="interactions.png")
```

## Example Output

![PandaMap](test/complex_interactions.png)

## Citation

If you use PandaMap in your research, please cite:

```
Pritam Kumar Panda. (2025). Protein AND ligAnd interaction MAPper: A Python package for visualizing protein-ligand interactions with 2D ligand structure representation. GitHub repository. https://github.com/pritampanda15/PandaMap
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

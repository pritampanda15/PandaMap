# PandaMap: A Python package for visualizing protein-ligand interactions with 2D ligand structure representation. 

**P**rotein **AND** lig**A**nd interaction **MAP**per: A Python package for visualizing protein-ligand interactions with 2D ligand structure representation

<p align="center">
  <img src="logo/pandamap-logo.svg" alt="PandaMap Logo" width="400">
</p>

## Overview

PandaMap is a lightweight tool for visualizing protein-ligand interactions from PDB files. It generates intuitive 2D interaction diagrams that display both the ligand structure and its interactions with protein residues.

Key features:
- Multiple structure format support e.g., pdb, mmcif, cif, pdbqt
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
- dssp #It can be installed externally
```bash
brew install dssp #mac users
sudo apt-get install dssp #linux users
Windows: Download from https://swift.cmbi.umcn.nl/gv/dssp/ 
```
- NumPy
- Matplotlib
- BioPython

## Basic Usage
```bash
usage: pandamap [-h] [--output OUTPUT] [--ligand LIGAND] [--dpi DPI]
                [--title TITLE] [--version]
                structure_file

PandaMap: Visualize protein-ligand interactions from structure files

positional arguments:
  structure_file        Path to structure file (PDB, mmCIF/CIF, or PDBQT
                        format)

options:
  -h, --help            show this help message and exit
  --output OUTPUT, -o OUTPUT
                        Output image file path
  --ligand LIGAND, -l LIGAND
                        Specific ligand residue name to analyze
  --dpi DPI             Image resolution (default: 300 dpi)
  --title TITLE, -t TITLE
                        Custom title for the visualization
  --version, -v         Show version information

```

### Command Line Interface

```bash
# Basic usage
pandamap protein_ligand.pdb --output interactions.png
pandamap complex.cif --output cif_interaction.png


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
# Using external DSSP (recommended)
mapper = HybridProtLigMapper("protein_ligand.pdb")
mapper.run_analysis(use_dssp=True)

# Using pure Python implementation
mapper.run_analysis(use_dssp=False)
## Example Output

![PandaMap](test/PAH_1.png)
![PandaMap](test/complex_interactions.png)
![PandaMap](test/HEM.png)


## Citation

If you use PandaMap in your research, please cite:

```
Pritam Kumar Panda. (2025). Protein AND ligAnd interaction MAPper: A Python package for visualizing protein-ligand interactions with 2D ligand structure representation. GitHub repository. https://github.com/pritampanda15/PandaMap
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

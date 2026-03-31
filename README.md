# PandaMap: A Python Package for Visualizing Protein–Ligand Interactions

**P**rotein **AND** lig**A**nd interaction **MAP**per — comprehensive detection, visualization, and empirical binding affinity estimation for protein–ligand complexes.

<p align="center">
  <img src="https://raw.githubusercontent.com/pritampanda15/PandaMap/main/logo/pandamap-logo.svg" alt="PandaMap Logo" width="400">
</p>
<p align="center">
  <a href="https://pypi.org/project/pandamap/">
    <img src="https://img.shields.io/pypi/v/pandamap.svg" alt="PyPI Version">
  </a>
  <a href="https://github.com/pritampanda15/PandaMap/blob/main/LICENSE">
    <img src="https://img.shields.io/github/license/pritampanda15/PandaMap" alt="License">
  </a>
  <a href="https://github.com/pritampanda15/PandaMap/stargazers">
    <img src="https://img.shields.io/github/stars/pritampanda15/PandaMap?style=social" alt="GitHub Stars">
  </a>
  <a href="https://github.com/pritampanda15/PandaMap/issues">
    <img src="https://img.shields.io/github/issues/pritampanda15/PandaMap" alt="GitHub Issues">
  </a>
  <a href="https://github.com/pritampanda15/PandaMap/network/members">
    <img src="https://img.shields.io/github/forks/pritampanda15/PandaMap?style=social" alt="GitHub Forks">
  </a>
  <a href="https://pepy.tech/project/PandaMap">
    <img src="https://static.pepy.tech/badge/PandaMap" alt="Downloads">
  </a>
</p>

---

## What's New in v4.2

| Feature | Description |
|---|---|
| **Scientifically validated cutoffs** | H-bond lower bound 2.5 Å, ionic/salt bridge 5.5 Å, halogen lower bound 2.5 Å — aligned with PLIP and crystallographic surveys |
| **RDKit 2D coordinates** | Chemically accurate 2D ligand layout when RDKit is available; PCA projection fallback requires no extra dependencies |
| **Topology-based ring detection** | Iterative leaf-node pruning on the bond graph — identifies rings of any size without SMILES |
| **Exact aromatic atom filtering** | Per-residue ring-atom name sets (PHE/TYR/TRP/HIS) eliminate false π-system contacts from β-carbons |
| **Improved ΔG estimation** | Per-residue deduplication + distance decay + rotatable bond entropy penalty; thermodynamically calibrated Kd labels at 298 K |
| **Trajectory analysis** | Multi-frame PDB/NMR ensemble analysis with per-interaction occupancy statistics and CSV export |
| **Arc stagger in 2D diagram** | Multiple interactions to the same residue fan out with alternating curvature — no marker overlap |
| **`--deltaG` CLI flag** | Estimate binding free energy directly from the command line |

---

## Features

- **16 interaction classes** detected with crystallographically validated distance thresholds:
  - Hydrogen bonds (2.5–3.5 Å, N/O pairs)
  - π–π stacking, cation–π, pi–cation, carbon–π, donor–π, amide–π, alkyl–π
  - Hydrophobic contacts (4.0 Å)
  - Ionic interactions and salt bridges (5.5 Å)
  - Halogen bonds (2.5–3.5 Å)
  - Metal coordination (2.8 Å)
  - Covalent bonds (≤2.1 Å, CYS/SER/LYS/HIS)
  - Attractive and repulsive charge interactions
- **Three output formats**: 2D PNG diagram, interactive 3D HTML (3Dmol.js), plain-text report
- **Empirical ΔG scoring** with per-residue deduplication, distance decay, and rotor penalty
- **Multi-frame trajectory analysis** with occupancy statistics and CSV export
- **Solvent accessibility**: DSSP (preferred) → Shrake–Rupley fallback → geometric fallback
- **Multi-format input**: PDB, mmCIF/CIF, PDBQT (AutoDock Vina)

👉 [Live interactive 3D example](https://github.com/pritampanda15/PandaMap/blob/main/test/complex_3d_visualization.html)

![PandaMap](https://raw.githubusercontent.com/pritampanda15/PandaMap/main/test/pandamap_3d_visualization.png)

---

## Installation

```bash
pip install pandamap
```

### Optional dependencies

```bash
pip install pandamap[fancy]   # coloured CLI output (rich)
pip install pandamap[viz]     # programmatic 3D viewer (py3Dmol)
pip install pandamap[full]    # all extras
```

### External optional: DSSP (accurate solvent accessibility)

```bash
brew install dssp                          # macOS
sudo apt-get install dssp                  # Linux
# Windows: https://swift.cmbi.umcn.nl/gv/dssp/
```

RDKit (for chemically accurate 2D ligand coordinates):
```bash
conda install -c conda-forge rdkit        # recommended
pip install rdkit                          # pip alternative
```
PandaMap works without RDKit — it falls back to PCA-based 2D projection automatically.

---

## Quick Start

```bash
# 2D interaction diagram
pandamap structure.pdb

# Specify ligand, generate report and 3D viewer
pandamap complex.pdb --ligand PFL --report --3d

# Estimate binding free energy
pandamap complex.pdb --ligand PFL --deltaG

# Full analysis
pandamap complex.pdb --ligand LIG --report --3d --deltaG --dpi 300
```

---

## Command-Line Reference

```
pandamap <structure_file> [options]

Positional arguments:
  structure_file          Path to PDB, mmCIF/CIF, or PDBQT file

Options:
  -l, --ligand NAME       Three-letter residue code of the ligand (default: auto-detect)
  -o, --output FILE       Output PNG file path
  -r, --report            Generate plain-text interaction report
  --report-file FILE      Path for text report
  --3d                    Generate interactive 3D HTML visualization
  --3d-output FILE        Path for 3D HTML file
  --deltaG                Estimate binding free energy (ΔG, kcal/mol)
  --dpi DPI               Output PNG resolution (default: 300)
  -t, --title TEXT        Custom diagram title
  --width PX              3D viewer width in pixels (default: 800)
  --height PX             3D viewer height in pixels (default: 600)
  --no-surface            Hide protein surface in 3D viewer
  --no-3d-cues            Disable depth cues in 2D diagram
  -v, --version           Show version
  -h, --help              Show help
```

---

## Python API

### Single-structure analysis

```python
from pandamap import HybridProtLigMapper

mapper = HybridProtLigMapper("complex.pdb", ligand_resname="LIG")
mapper.detect_interactions()

# 2D diagram
mapper.visualize(output_file="interactions.png")

# Text report
from pandamap.improved_interaction_detection import ImprovedInteractionDetection
detector = ImprovedInteractionDetection()
detector.generate_report(
    ligand_metadata={
        'hetid': mapper.ligand_residue.resname,
        'chain': mapper.ligand_residue.parent.id,
        'position': mapper.ligand_residue.id[1],
        'longname': mapper.ligand_residue.resname,
        'type': 'LIGAND',
    },
    interaction_data=mapper.interactions,
    output_file="report.txt"
)

# Inspect raw interactions
for itype, contacts in mapper.interactions.items():
    if contacts:
        print(f"{itype}: {len(contacts)} contacts")
```

### Empirical ΔG estimation

```python
result = mapper.estimate_binding_affinity()

print(f"ΔG ≈ {result['dG_estimated']:.2f} kcal/mol")
print(result['interpretation'])
print(result['note'])

for itype, info in result['breakdown'].items():
    print(f"  {itype} (n={info['unique_residues']}): {info['contribution_kcal_mol']:+.2f} kcal/mol")
```

**ΔG interpretation (298 K, ΔG = −RT·ln Kd):**

| ΔG (kcal/mol) | Kd range | Label |
|---|---|---|
| < −12 | ~nM or better | Very strong binder |
| −9 to −12 | nM–µM | Strong binder |
| −6 to −9 | µM | Moderate binder |
| −3 to −6 | mM | Weak binder |
| ≥ −3 | — | Very weak / no binding |

> Note: Empirical estimate ±2–3 kcal/mol. Not a substitute for FEP or MM-GBSA.

### 3D visualization

```python
from pandamap.create_3d_view import create_pandamap_3d_viz

create_pandamap_3d_viz(
    mapper=mapper,
    output_file="interactions_3d.html",
    width=1024,
    height=768,
    show_surface=True
)
```

### Multi-frame trajectory analysis

```python
from pandamap import analyze_trajectory

summary = analyze_trajectory(
    trajectory_file="simulation.pdb",    # multi-MODEL PDB
    ligand_resname="LIG",
    output_dir="./trajectory_output",
    visualize_frames=False               # set True to generate per-frame PNGs
)

print(f"Frames analysed: {summary['n_frames']}")
print(f"Mean ΔG: {summary['mean_dG']:.2f} ± {summary['std_dG']:.2f} kcal/mol")
# Per-residue occupancy CSV written to ./trajectory_output/trajectory_analysis.csv
```

---

## Distance Cutoffs

All cutoffs are validated against PLIP and published crystallographic surveys:

| Interaction | Cutoff | Reference |
|---|---|---|
| Hydrogen bond | 2.5–3.5 Å | PLIP; Auffinger 2004 |
| π–π stacking | 5.5 Å (atom–atom) | McGaughey 1998 |
| Hydrophobic | 4.0 Å | Bissantz 2010 |
| Ionic / salt bridge | 5.5 Å | Kumar & Nussinov 1999 |
| Halogen bond | 2.5–3.5 Å | Auffinger 2004 |
| Metal coordination | 2.8 Å | CSD surveys |
| Covalent | ≤2.1 Å | — |
| Repulsion | 4.0 Å | — |

---

## Example Outputs

### 2D Interaction Diagram

![PandaMap](https://raw.githubusercontent.com/pritampanda15/PandaMap/main/test/complex_interactions.png)
![PandaMap](https://raw.githubusercontent.com/pritampanda15/PandaMap/main/test/1m17_interactions.png)
![PandaMap](https://raw.githubusercontent.com/pritampanda15/PandaMap/main/test/4jmz_interactions.png)

### Text Report

```
=============================================================================
PandaMap Interaction Report
=============================================================================

Ligand: PAH:A:439
Name: PAH
Type: LIGAND

Interacting Chains: A
Interacting Residues: 13

Interaction Summary:
  Hydrogen Bonds: 10
  Carbon-π Interactions: 1
  Metal Coordination: 4
  Ionic Interactions: 2
  Salt Bridges: 2
  Alkyl-π Interactions: 1
  Attractive Charge: 2
  Repulsion: 5

Hydrogen Bonds:
  1. GLU168A  -- 2.66Å -- PAH
  2. ASP246A  -- 2.60Å -- PAH
  3. GLN167A  -- 3.10Å -- PAH
  4. ASP320A  -- 3.46Å -- PAH
  5. LYS396A  -- 3.05Å -- PAH
  ...

=============================================================================
```

### ΔG Estimation Output

```
--- Estimated Binding Affinity ---
  ΔG ≈ -7.42 kcal/mol
  Strong binder (Kd ~nM–µM range)
  Breakdown:
    hydrogen_bonds (n=10): -8.63 kcal/mol
    metal_coordination (n=4): -6.80 kcal/mol
    ionic (n=2): -2.91 kcal/mol
    hydrophobic (n=3): -0.72 kcal/mol
    rotatable_bond_penalty (n=2): +1.00 kcal/mol
----------------------------------
```

---

## Citation

If you use PandaMap in your research, please cite:

```
Pritam Kumar Panda. (2025). PandaMap: A Python Package for Comprehensive
Visualization of Protein–Ligand Interaction Networks and Empirical Binding
Affinity Estimation. Boring Science LLC.
https://github.com/pritampanda15/PandaMap
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.

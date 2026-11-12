# PandaMap examples

Five structures and everything PandaMap 4.4 writes for each of them. Use this
page to see what the outputs look like before running anything, or to check
that your installation gives the same results.

Each structure was run with one command:

```bash
pandamap example/1m17.pdb --ligand AQ4 --outdir example/outputs/1m17 \
    --report --plots --csv --json --score --3d
```

To regenerate every file here with your installed version, run
`bash example/regenerate.sh` from the repository root.

| Structure | Ligand | What it is | Contacts | Residues |
|---|---|---|---:|---:|
| [`1m17.pdb`](#1m17) | AQ4 | EGFR kinase with erlotinib | 19 | 15 |
| [`1els.pdb`](#1els) | PAH | Enolase with a phosphonate inhibitor | 28 | 14 |
| [`4jmz.pdb`](#4jmz) | HEM | Heme in a peroxidase | 51 | 23 |
| [`complex.pdb`](#complex) | PFL | Propofol at a subunit interface | 13 | 10 |
| [`A0152a.pdb`](#a0152a) | LIG | Ligand in a protein with a structural zinc site | 6 | 5 |

## What is in each folder

| File | Contents |
|---|---|
| `<name>_interactions.png` | Modern 2D diagram |
| `<name>_interactions_classic.png` | Classic diagram |
| `<name>_interactions_3d.png` | 3D contact view |
| `<name>_interactions_report.png` | Five-panel summary figure |
| `<name>_interactions_report.txt` | Text report |
| `<name>_interactions.csv` | One row per contact |
| `<name>_interactions.json` | The whole analysis, including the interaction score |
| `<name>_interaction_score.txt` | What `--score` prints |
| `<name>_3d_visualization.html` | Interactive 3D viewer; download it and open it in a browser |

## 1m17

**EGFR kinase with erlotinib** (ligand `AQ4`, chain A, residue 999). A drug-like ligand in a kinase pocket: hinge hydrogen bond to MET 769 and a water bridge to THR 766.

Detected: 8 alkyl pi, 1 amide pi, 2 donor pi, 4 hydrogen bonds, 1 hydrophobic, 1 ionic, 1 pi cation, 1 water bridges.

| Modern 2D | Classic | 3D contacts |
|:---:|:---:|:---:|
| [<img src="outputs/1m17/1m17_interactions.png" width="260">](outputs/1m17/1m17_interactions.png) | [<img src="outputs/1m17/1m17_interactions_classic.png" width="260">](outputs/1m17/1m17_interactions_classic.png) | [<img src="outputs/1m17/1m17_interactions_3d.png" width="260">](outputs/1m17/1m17_interactions_3d.png) |

[<img src="outputs/1m17/1m17_interactions_report.png" width="640">](outputs/1m17/1m17_interactions_report.png)

Data: [text report](outputs/1m17/1m17_interactions_report.txt) · [CSV](outputs/1m17/1m17_interactions.csv) · [JSON](outputs/1m17/1m17_interactions.json) · [interaction score](outputs/1m17/1m17_interaction_score.txt) · [3D viewer (HTML)](outputs/1m17/1m17_3d_visualization.html)

## 1els

**Enolase with a phosphonate inhibitor** (ligand `PAH`, chain A, residue 439). A small, highly charged ligand held by two metal ions, many hydrogen bonds and water bridges.

Detected: 1 alkyl pi, 13 hydrogen bonds, 2 ionic, 4 metal coordination, 8 water bridges.

| Modern 2D | Classic | 3D contacts |
|:---:|:---:|:---:|
| [<img src="outputs/1els/1els_interactions.png" width="260">](outputs/1els/1els_interactions.png) | [<img src="outputs/1els/1els_interactions_classic.png" width="260">](outputs/1els/1els_interactions_classic.png) | [<img src="outputs/1els/1els_interactions_3d.png" width="260">](outputs/1els/1els_interactions_3d.png) |

[<img src="outputs/1els/1els_interactions_report.png" width="640">](outputs/1els/1els_interactions_report.png)

Data: [text report](outputs/1els/1els_interactions_report.txt) · [CSV](outputs/1els/1els_interactions.csv) · [JSON](outputs/1els/1els_interactions.json) · [interaction score](outputs/1els/1els_interaction_score.txt) · [3D viewer (HTML)](outputs/1els/1els_3d_visualization.html)

## 4jmz

**Heme in a peroxidase** (ligand `HEM`, chain A, residue 301). A large cofactor with 23 contacting residues; shows how the diagrams scale and how a metal complex is drawn.

Detected: 16 alkyl pi, 2 carbon pi, 1 donor pi, 13 hydrogen bonds, 5 hydrophobic, 3 ionic, 3 pi cation, 4 pi pi stacking, 4 water bridges.

| Modern 2D | Classic | 3D contacts |
|:---:|:---:|:---:|
| [<img src="outputs/4jmz/4jmz_interactions.png" width="260">](outputs/4jmz/4jmz_interactions.png) | [<img src="outputs/4jmz/4jmz_interactions_classic.png" width="260">](outputs/4jmz/4jmz_interactions_classic.png) | [<img src="outputs/4jmz/4jmz_interactions_3d.png" width="260">](outputs/4jmz/4jmz_interactions_3d.png) |

[<img src="outputs/4jmz/4jmz_interactions_report.png" width="640">](outputs/4jmz/4jmz_interactions_report.png)

Data: [text report](outputs/4jmz/4jmz_interactions_report.txt) · [CSV](outputs/4jmz/4jmz_interactions.csv) · [JSON](outputs/4jmz/4jmz_interactions.json) · [interaction score](outputs/4jmz/4jmz_interaction_score.txt) · [3D viewer (HTML)](outputs/4jmz/4jmz_3d_visualization.html)

## complex

**Propofol at a subunit interface** (ligand `PFL`, chain D, residue 601). A binding site made by two chains (C and D); residue labels carry the chain.

Detected: 7 alkyl pi, 1 amide pi, 1 donor pi, 1 hydrogen bonds, 2 hydrophobic, 1 pi pi stacking.

| Modern 2D | Classic | 3D contacts |
|:---:|:---:|:---:|
| [<img src="outputs/complex/complex_interactions.png" width="260">](outputs/complex/complex_interactions.png) | [<img src="outputs/complex/complex_interactions_classic.png" width="260">](outputs/complex/complex_interactions_classic.png) | [<img src="outputs/complex/complex_interactions_3d.png" width="260">](outputs/complex/complex_interactions_3d.png) |

[<img src="outputs/complex/complex_interactions_report.png" width="640">](outputs/complex/complex_interactions_report.png)

Data: [text report](outputs/complex/complex_interactions_report.txt) · [CSV](outputs/complex/complex_interactions.csv) · [JSON](outputs/complex/complex_interactions.json) · [interaction score](outputs/complex/complex_interaction_score.txt) · [3D viewer (HTML)](outputs/complex/complex_3d_visualization.html)

## A0152a

**Ligand in a protein with a structural zinc site** (ligand `LIG`, chain A, residue 147). A small, mostly hydrophobic site. The zinc is 32 Å from the ligand and is correctly left out.

Detected: 2 alkyl pi, 1 donor pi, 2 hydrophobic, 1 ionic.

| Modern 2D | Classic | 3D contacts |
|:---:|:---:|:---:|
| [<img src="outputs/A0152a/A0152a_interactions.png" width="260">](outputs/A0152a/A0152a_interactions.png) | [<img src="outputs/A0152a/A0152a_interactions_classic.png" width="260">](outputs/A0152a/A0152a_interactions_classic.png) | [<img src="outputs/A0152a/A0152a_interactions_3d.png" width="260">](outputs/A0152a/A0152a_interactions_3d.png) |

[<img src="outputs/A0152a/A0152a_interactions_report.png" width="640">](outputs/A0152a/A0152a_interactions_report.png)

Data: [text report](outputs/A0152a/A0152a_interactions_report.txt) · [CSV](outputs/A0152a/A0152a_interactions.csv) · [JSON](outputs/A0152a/A0152a_interactions.json) · [interaction score](outputs/A0152a/A0152a_interaction_score.txt) · [3D viewer (HTML)](outputs/A0152a/A0152a_3d_visualization.html)

## Batch mode

All five structures in one run, without images:

```bash
pandamap example/ --outdir example/outputs/batch --diagrams none
```

- [`batch_summary.csv`](outputs/batch/batch_summary.csv): one row per structure.
- [`batch_interactions.csv`](outputs/batch/batch_interactions.csv): one row per contact, with a `structure` column.

Without `--ligand`, batch mode picks the ligand in each file automatically.

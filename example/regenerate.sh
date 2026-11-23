#!/usr/bin/env bash
# Regenerate every file under example/outputs/ with the installed PandaMap.
# Run from the repository root:  bash example/regenerate.sh
set -euo pipefail
cd "$(dirname "$0")/.."

run() {  # run <structure> <ligand>
  pandamap "example/$1.pdb" --ligand "$2" --outdir "example/outputs/$1" \
    --report --plots --csv --json --score --3d --dpi 150 --quiet \
    > "example/outputs/$1/$1_interaction_score.txt"
}

rm -rf example/outputs
for name in 1m17 1els 4jmz complex A0152a; do mkdir -p "example/outputs/$name"; done

run 1m17    AQ4   # EGFR kinase with erlotinib
run 1els    PAH   # enolase with a phosphonate inhibitor and two metal ions
run 4jmz    HEM   # heme in a peroxidase
run complex PFL   # propofol at a subunit interface (two chains)
run A0152a  LIG   # ligand with a structural zinc site elsewhere in the protein

# Batch mode over the whole folder: two combined tables, no images
pandamap example/ --outdir example/outputs/batch --diagrams none --quiet

# Water bridges with 3D-RISM waters: PLIP angle test, then distances only
wb() {  # wb <folder> [extra options]
  local out="example/water_bridge/$1"; shift
  rm -rf "$out"; mkdir -p "$out"
  pandamap example/water_bridge/complex_rism_waters.pdb --ligand UNK --outdir "$out" \
    --report --plots --csv --json --score --3d --dpi 150 --quiet "$@" \
    > "$out/complex_rism_waters_interaction_score.txt"
}
wb angle_test
wb distance_only --water-distance-only

# The modern diagram of 1M17 in each colour scheme (from PandaMap-Color)
rm -rf example/color_schemes; mkdir -p example/color_schemes
for scheme in default colorblind monochrome dark publication; do
  pandamap example/1m17.pdb --ligand AQ4 --diagrams modern --color-scheme "$scheme" \
    --output "example/color_schemes/1m17_$scheme.png" --dpi 150 --quiet
done

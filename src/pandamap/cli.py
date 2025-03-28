#!/usr/bin/env python
"""
Command-line interface for PandaMap.
"""

import sys
import argparse
from pandamap.core import HybridProtLigMapper

def main():
    """Command-line interface for PandaMap."""
    parser = argparse.ArgumentParser(
        description='PandaMap: Visualize protein-ligand interactions from PDB files')
    
    parser.add_argument('pdb_file', help='Path to PDB file with protein-ligand complex')
    parser.add_argument('--output', '-o', help='Output image file path')
    parser.add_argument('--ligand', '-l', help='Specific ligand residue name to analyze')
    parser.add_argument('--dpi', type=int, default=300, help='Image resolution (default: 300 dpi)')
    parser.add_argument('--title', '-t', help='Custom title for the visualization')
    parser.add_argument('--version', '-v', action='store_true', help='Show version information')
    
    args = parser.parse_args()
    
    # Handle version flag
    if args.version:
        from pandamap import __version__
        print(f"PandaMap version {__version__}")
        return 0
    
    # Check for required pdb_file
    if not hasattr(args, 'pdb_file') or args.pdb_file is None:
        parser.print_help()
        return 1
    
    try:
        # Initialize and run the analysis
        mapper = HybridProtLigMapper(args.pdb_file, ligand_resname=args.ligand)
        output_file = mapper.run_analysis(output_file=args.output)
        
        print(f"Analysis complete. Visualization saved to: {output_file}")
        return 0
        
    except Exception as e:
        print(f"Error: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python
"""
Command-line interface for PandaMap.
"""

import sys
import os
import argparse
from pandamap.core import HybridProtLigMapper

def main():
    """Command-line interface for PandaMap."""
    parser = argparse.ArgumentParser(
        description='PandaMap: Visualize protein-ligand interactions from structure files')
    
    parser.add_argument('structure_file', 
                      help='Path to structure file (PDB, mmCIF/CIF, or PDBQT format)')
    parser.add_argument('--output', '-o', help='Output image file path')
    parser.add_argument('--ligand', '-l', help='Specific ligand residue name to analyze')
    parser.add_argument('--dpi', type=int, default=300, help='Image resolution (default: 300 dpi)')
    parser.add_argument('--title', '-t', help='Custom title for the visualization')
    parser.add_argument('--version', '-v', action='store_true', help='Show version information')
    parser.add_argument('--report', '-r', action='store_true', 
                       help='Generate text report')
    parser.add_argument('--report-file', 
                       help='Output file for the text report (default: based on structure filename)')
    
    
    args = parser.parse_args()
    
    # Handle version flag
    if args.version:
        from pandamap import __version__
        print(f"PandaMap version {__version__}")
        return 0
    
    # Check file existence
    if not os.path.exists(args.structure_file):
        print(f"Error: File not found: {args.structure_file}")
        return 1
    
    # Check file extension
    file_ext = os.path.splitext(args.structure_file)[1].lower()
    if file_ext not in ['.pdb', '.cif', '.mmcif', '.pdbqt']:
        print(f"Warning: Unrecognized file extension: {file_ext}")
        print("Supported formats: .pdb, .cif, .mmcif, .pdbqt")
        choice = input("Attempt to parse anyway? (y/n): ")
        if choice.lower() != 'y':
            return 1
    
    try:
        mapper = HybridProtLigMapper(args.structure_file, ligand_resname=args.ligand)
        output_file = mapper.run_analysis(
            output_file=args.output,
            generate_report=args.report,
            report_file=args.report_file
        )
        
        print(f"Analysis complete. Visualization saved to: {output_file}")
        if args.report:
            report_file = args.report_file or f"{os.path.splitext(output_file)[0]}_report.txt"
            print(f"Interaction report saved to: {report_file}")
        return 0
        
    except Exception as e:
        print(f"Error: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
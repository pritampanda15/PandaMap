#!/usr/bin/env python
"""
Command-line interface for PandaMap.
"""
import sys
import os
import argparse
from pandamap.core import HybridProtLigMapper
from pandamap.create_3d_view import create_pandamap_3d_viz


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
    
    # Add 3D visualization options
    parser.add_argument('--3d', dest='use_3d', action='store_true',
        help='Generate an interactive 3D visualization of protein-ligand interactions')
    parser.add_argument('--3d-output', dest='output_3d',
        help='Output file path for 3D visualization (default: based on input filename)')
    parser.add_argument('--no-surface', action='store_true',
        help='Do not show protein surface in 3D visualization')
    parser.add_argument('--width', type=int, default=800,
        help='Width of 3D visualization in pixels (default: 800)')
    parser.add_argument('--height', type=int, default=600,
        help='Height of 3D visualization in pixels (default: 600)')
    parser.add_argument('--no-3d-cues', action='store_true',
        help='Disable 3D cues in 2D visualization')
    
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
    
    # Set up 3D visualization output file if needed
    if args.use_3d and args.output_3d is None:
        base_name = os.path.splitext(os.path.basename(args.structure_file))[0]
        args.output_3d = f"{base_name}_3d_visualization.html"
    
    try:
        # Create the mapper and run the analysis
        mapper = HybridProtLigMapper(args.structure_file, ligand_resname=args.ligand)
        
        # Run the standard analysis
        output_file = mapper.run_analysis(
            output_file=args.output,
            generate_report=args.report,
            report_file=args.report_file
        )
        
        print(f"Analysis complete. Visualization saved to: {output_file}")
        
        if args.report:
            report_file = args.report_file or f"{os.path.splitext(output_file)[0]}_report.txt"
            print(f"Interaction report saved to: {report_file}")
        
        # Generate 3D visualization if requested
        if args.use_3d:
            try:
                print("Generating 3D visualization...")
                create_pandamap_3d_viz(
                    mapper=mapper,
                    output_file=args.output_3d,
                    width=args.width,
                    height=args.height,
                    show_surface=not args.no_surface
                )
            except Exception as e:
                print(f"\nError generating 3D visualization: {str(e)}")
                import traceback
                traceback.print_exc()
                return 1
        
        return 0
        
    except Exception as e:
        print(f"Error: {e}")
        # Uncomment the next line for debugging
        # import traceback; traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
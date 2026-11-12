#!/usr/bin/env python
"""
Command-line interface for PandaMap.
"""
import sys
import os
import argparse
from pandamap.core import HybridProtLigMapper, list_ligands
from pandamap.create_3d_view import create_pandamap_3d_viz


def _make_console_safe():
    """Never let a non-UTF-8 console (e.g. Windows cp1252, or redirected
    output) crash the run on characters such as the Greek letters in the
    interaction names: unencodable characters are replaced instead."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors='replace')
        except (AttributeError, ValueError):
            pass


def main():
    """Command-line interface for PandaMap."""
    _make_console_safe()
    from pandamap import __version__
    parser = argparse.ArgumentParser(
        description='PandaMap: Visualize protein-ligand interactions from structure files')
    parser.add_argument('structure_file',
        help='Path to structure file (PDB, mmCIF/CIF, or PDBQT format)')
    parser.add_argument('--output', '-o',
        help='Output image file path; the extension selects the format '
             '(.png, .svg, .pdf)')
    parser.add_argument('--ligand', '-l',
        help='Ligand residue name to analyze (default: auto-detect the largest '
             'non-solvent HETATM group)')
    parser.add_argument('--chain', help='Chain ID of the ligand, when it occurs more than once')
    parser.add_argument('--resnum', type=int,
        help='Residue number of the ligand, when it occurs more than once')
    parser.add_argument('--list-ligands', action='store_true',
        help='List the candidate ligands in the structure file and exit')
    parser.add_argument('--debug', action='store_true',
        help='Show the full Python traceback when an error occurs')
    parser.add_argument('--dpi', type=int, default=300, help='Image resolution (default: 300 dpi)')
    parser.add_argument('--title', '-t', help='Custom title for the visualization')
    parser.add_argument('--version', '-v', action='version',
        version=f'PandaMap version {__version__}', help='Show version information')
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
        help='Disable depth cues in the classic diagram')
    parser.add_argument('--diagrams', '-d', nargs='+', default=['all'],
        choices=['modern', 'classic', '3d', 'all'], metavar='STYLE',
        help="Which diagrams to draw: 'modern' (2D, residues placed next to the "
             "ligand atoms they contact), 'classic' (the original circular "
             "PandaMap diagram), '3d' (transparent residues around the ligand "
             "with each contacting atom highlighted), or 'all' (default). "
             "Several can be given, e.g. --diagrams modern 3d. The first is "
             "written to --output; the others get a _classic / _3d suffix")
    parser.add_argument('--show-nearby', action='store_true',
        help='Also show binding-site residues without a specific interaction '
             '(modern diagram)')
    parser.add_argument('--deltaG', action='store_true',
        help='Estimate binding free energy (ΔG) using empirical scoring')
    parser.add_argument('--csv', nargs='?', const=True, default=False,
        metavar='FILE',
        help='Export all detected interactions as a machine-readable CSV')
    parser.add_argument('--json', nargs='?', const=True, default=False,
        metavar='FILE',
        help='Export the analysis as JSON: ligand, summary, every interaction '
             'with its atoms and distance, and a per-residue roll-up '
             '(includes the ΔG estimate when --deltaG is given)')
    parser.add_argument('--plots', nargs='?', const=True, default=False,
        metavar='FILE',
        help='Generate a four-panel graphical interaction report (PNG)')
    parser.add_argument('--strict-hbond-geometry', action='store_true',
        dest='strict_hbond_geometry',
        help='Also apply the heavy-atom angular proxy when the structure has '
             'no explicit hydrogens (the exact D-H...A test always applies '
             'when hydrogens are present)')
    parser.add_argument('--strict-halogens', action='store_true',
        dest='strict_halogens',
        help='Restrict halogen bonds to Cl/Br/I (Auffinger 2004); '
             'by default C-F...O contacts are also reported, as in PLIP')

    args = parser.parse_args()
    
    # Check file existence
    if not os.path.exists(args.structure_file):
        print(f"Error: File not found: {args.structure_file}")
        return 1
    
    # Unknown extensions are parsed by content (PDB unless it looks like mmCIF)
    file_ext = os.path.splitext(args.structure_file)[1].lower()
    if file_ext not in ['.pdb', '.cif', '.mmcif', '.pdbqt']:
        print(f"Note: unrecognized file extension '{file_ext}'; "
              "detecting the format from the file content.")

    if args.list_ligands:
        try:
            ligands = list_ligands(args.structure_file)
        except Exception as e:
            print(f"Error: {e}")
            return 1
        if not ligands:
            print("No ligand (HETATM) found.")
            return 1
        print(f"{'Ligand':<8}{'Chain':<7}{'Residue':<9}Atoms")
        for lig in ligands:
            print(f"{lig['resname']:<8}{lig['chain']:<7}{lig['resnum']:<9}{lig['n_atoms']}")
        return 0

    # Set up 3D visualization output file if needed
    if args.use_3d and args.output_3d is None:
        base_name = os.path.splitext(os.path.basename(args.structure_file))[0]
        args.output_3d = f"{base_name}_3d_visualization.html"
    
    try:
        diagrams = ['modern', 'classic', '3d'] if 'all' in args.diagrams \
            else list(dict.fromkeys(args.diagrams))
        # Create the mapper and run the analysis
        mapper = HybridProtLigMapper(args.structure_file, ligand_resname=args.ligand,
                                     ligand_chain=args.chain,
                                     ligand_resnum=args.resnum,
                                     strict_halogens=args.strict_halogens,
                                     strict_hbond_geometry=args.strict_hbond_geometry)
        
        # Run the standard analysis
        output_file = mapper.run_analysis(
            output_file=args.output,
            generate_report=args.report,
            report_file=args.report_file,
            dpi=args.dpi,
            title=args.title,
            show_3d_cues=not args.no_3d_cues,
            diagrams=diagrams,
            show_nearby_residues=args.show_nearby
        )
        
        print("Analysis complete. Diagrams saved:")
        for style, path in mapper.diagram_files.items():
            print(f"  {style:<8}{path}")

        if args.deltaG:
            result = mapper.estimate_binding_affinity()
            print("\n--- Estimated Binding Affinity ---")
            print(f"  ΔG ≈ {result['dG_estimated']:.2f} kcal/mol")
            print(f"  {result['interpretation']}")
            print("  Breakdown:")
            for itype, info in result['breakdown'].items():
                n = info.get('unique_residues', info.get('count', '?'))
                print(f"    {itype} (n={n}): {info['contribution_kcal_mol']:+.2f} kcal/mol")
            print("----------------------------------")
        
        if args.csv:
            csv_path = args.csv if isinstance(args.csv, str) else \
                f"{os.path.splitext(output_file)[0]}.csv"
            mapper.export_interactions_csv(csv_path)

        if args.json:
            json_path = args.json if isinstance(args.json, str) else \
                f"{os.path.splitext(output_file)[0]}.json"
            mapper.export_interactions_json(json_path, include_affinity=args.deltaG)

        if args.plots:
            plot_path = args.plots if isinstance(args.plots, str) else \
                f"{os.path.splitext(output_file)[0]}_report.png"
            mapper.generate_interaction_plots(plot_path, dpi=args.dpi)

        if args.report:
            report_file = args.report_file or (
                os.path.splitext(os.path.basename(args.structure_file))[0]
                + "_interactions_report.txt")
            if os.path.exists(report_file):
                print(f"Interaction report saved to: {report_file}")
            else:
                print(f"Warning: the interaction report could not be written to {report_file}")
        
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
        if args.debug:
            import traceback
            traceback.print_exc()
        elif not isinstance(e, (ValueError, OSError)):
            print("Run again with --debug for the full traceback.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
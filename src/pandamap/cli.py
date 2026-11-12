#!/usr/bin/env python
"""
Command-line interface for PandaMap.
"""
import argparse
import logging
import os
import sys
import threading

from pandamap.core import HybridProtLigMapper, list_ligands
from pandamap.create_3d_view import create_pandamap_3d_viz
from pandamap.batch import analyze_batch, expand_structure_paths


def _make_console_safe():
    """Never let a non-UTF-8 console (e.g. Windows cp1252, or redirected
    output) crash the run on characters such as the Greek letters in the
    interaction names: unencodable characters are replaced instead."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors='replace')
        except (AttributeError, ValueError):
            pass


def _configure_logging(quiet):
    """Send the library's progress messages to the terminal (or, with
    --quiet, only its warnings).  Used from Python, PandaMap stays silent
    unless the caller configures the 'pandamap' logger."""
    logger = logging.getLogger("pandamap")
    for handler in [h for h in logger.handlers if getattr(h, '_pandamap_cli', False)]:
        logger.removeHandler(handler)
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter("%(message)s"))
    handler._pandamap_cli = True
    logger.addHandler(handler)
    logger.setLevel(logging.WARNING if quiet else logging.INFO)
    logger.propagate = False


def build_parser(version):
    parser = argparse.ArgumentParser(
        description='PandaMap: visualize protein-ligand interactions from structure files')
    parser.add_argument('structures', nargs='+', metavar='structure_file',
        help='Structure file (PDB, mmCIF/CIF or PDBQT). Give several files or a '
             'folder to analyse them all in batch mode')
    parser.add_argument('--output', '-o',
        help='Output image file path; the extension selects the format '
             '(.png, .svg, .pdf)')
    parser.add_argument('--outdir', '-O', metavar='DIR',
        help='Folder for all output files (created if needed; default: the '
             'current folder, or pandamap_batch in batch mode)')
    parser.add_argument('--ligand', '-l',
        help='Ligand residue name to analyze (default: auto-detect the largest '
             'non-solvent HETATM group)')
    parser.add_argument('--chain', help='Chain ID of the ligand, when it occurs more than once')
    parser.add_argument('--resnum', type=int,
        help='Residue number of the ligand, when it occurs more than once')
    parser.add_argument('--list-ligands', action='store_true',
        help='List the candidate ligands in the structure file and exit')
    parser.add_argument('--quiet', '-q', action='store_true',
        help='Print only warnings and errors')
    parser.add_argument('--debug', action='store_true',
        help='Show the full Python traceback when an error occurs')
    parser.add_argument('--dpi', type=int, default=300, help='Image resolution (default: 300 dpi)')
    parser.add_argument('--title', '-t', help='Custom title for the visualization')
    parser.add_argument('--version', '-v', action='version',
        version=f'PandaMap version {version}', help='Show version information')
    parser.add_argument('--report', '-r', action='store_true',
        help='Generate text report')
    parser.add_argument('--report-file',
        help='Output file for the text report (default: based on structure filename)')

    # 3D visualization options
    parser.add_argument('--3d', dest='use_3d', action='store_true',
        help='Generate an interactive 3D visualization of protein-ligand interactions')
    parser.add_argument('--3d-output', dest='output_3d',
        help='Output file path for 3D visualization (default: based on input filename)')
    parser.add_argument('--3d-offline', dest='offline_3d', action='store_true',
        help='Embed the 3Dmol.js library in the 3D page so it works without an '
             'internet connection (implies --3d; the library is downloaded once '
             'and cached)')
    parser.add_argument('--no-surface', action='store_true',
        help='Do not show protein surface in 3D visualization')
    parser.add_argument('--width', type=int, default=800,
        help='Width of 3D visualization in pixels (default: 800)')
    parser.add_argument('--height', type=int, default=600,
        help='Height of 3D visualization in pixels (default: 600)')
    parser.add_argument('--no-3d-cues', action='store_true',
        help='Disable depth cues in the classic diagram')
    parser.add_argument('--diagrams', '-d', nargs='+', default=None,
        choices=['modern', 'classic', '3d', 'all', 'none'], metavar='STYLE',
        help="Which diagrams to draw: 'modern' (2D, residues placed next to the "
             "ligand atoms they contact), 'classic' (the original circular "
             "PandaMap diagram), '3d' (transparent residues around the ligand "
             "with each contacting atom highlighted), 'all' or 'none'. "
             "Default: all for one structure, modern in batch mode. Several "
             "can be given, e.g. --diagrams modern 3d. The first is written to "
             "--output; the others get a _classic / _3d suffix")
    parser.add_argument('--show-nearby', action='store_true',
        help='Also show binding-site residues without a specific interaction '
             '(modern diagram)')
    parser.add_argument('--score', '--deltaG', dest='deltaG', action='store_true',
        help='Print the empirical interaction score and its breakdown. It is on a '
             'kcal/mol scale but is not calibrated against measured affinities '
             '(--deltaG is the older name for this option)')
    parser.add_argument('--csv', nargs='?', const=True, default=False,
        metavar='FILE',
        help='Export all detected interactions as a machine-readable CSV')
    parser.add_argument('--json', nargs='?', const=True, default=False,
        metavar='FILE',
        help='Export the analysis as JSON: ligand, summary, every interaction '
             'with its atoms and distance, and a per-residue roll-up '
             '(includes the interaction score when --score is given)')
    parser.add_argument('--plots', nargs='?', const=True, default=False,
        metavar='FILE',
        help='Draw the five-panel summary figure, built around the 2D interaction '
             'map (.png, .pdf or .svg)')
    parser.add_argument('--strict-hbond-geometry', action='store_true',
        dest='strict_hbond_geometry',
        help='Also apply the heavy-atom angular proxy when the structure has '
             'no explicit hydrogens (the exact D-H...A test always applies '
             'when hydrogens are present)')
    parser.add_argument('--strict-halogens', action='store_true',
        dest='strict_halogens',
        help='Restrict halogen bonds to Cl/Br/I (Auffinger 2004); '
             'by default C-F...O contacts are also reported, as in PLIP')
    return parser


def _selected_diagrams(choice, default):
    if choice is None:
        return list(default)
    if 'none' in choice:
        return []
    if 'all' in choice:
        return ['modern', 'classic', '3d']
    return list(dict.fromkeys(choice))


def _run_batch(args, files, say):
    outdir = args.outdir or 'pandamap_batch'
    result = analyze_batch(
        files, output_dir=outdir, ligand_resname=args.ligand,
        diagrams=_selected_diagrams(args.diagrams, ['modern']),
        report=args.report, plots=bool(args.plots), json_export=bool(args.json),
        dpi=args.dpi, strict_halogens=args.strict_halogens,
        strict_hbond_geometry=args.strict_hbond_geometry)
    failed = [r for r in result['results'] if r['status'] != 'ok']
    say(f"Batch summary saved to: {result['summary_csv']}")
    say(f"All contacts saved to: {result['interactions_csv']}")
    for row in failed:
        print(f"Failed: {row['file']}: {row['error']}")
    return 1 if len(failed) == len(result['results']) else 0


def _run_single(args, structure_file, say):
    stem = os.path.splitext(os.path.basename(structure_file))[0]
    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)
    base = os.path.join(args.outdir or '', stem)

    file_ext = os.path.splitext(structure_file)[1].lower()
    if file_ext not in ['.pdb', '.cif', '.mmcif', '.pdbqt']:
        say(f"Note: unrecognized file extension '{file_ext}'; "
            "detecting the format from the file content.")

    if args.list_ligands:
        ligands = list_ligands(structure_file)
        if not ligands:
            print("No ligand (HETATM) found.")
            return 1
        print(f"{'Ligand':<8}{'Chain':<7}{'Residue':<9}Atoms")
        for lig in ligands:
            print(f"{lig['resname']:<8}{lig['chain']:<7}{lig['resnum']:<9}{lig['n_atoms']}")
        return 0

    diagrams = _selected_diagrams(args.diagrams, ['modern', 'classic', '3d'])
    output_file = args.output or f"{base}_interactions.png"
    report_file = args.report_file or f"{base}_interactions_report.txt"
    out_stem = os.path.splitext(output_file)[0]

    mapper = HybridProtLigMapper(structure_file, ligand_resname=args.ligand,
                                 ligand_chain=args.chain,
                                 ligand_resnum=args.resnum,
                                 strict_halogens=args.strict_halogens,
                                 strict_hbond_geometry=args.strict_hbond_geometry)
    mapper.run_analysis(
        output_file=output_file,
        generate_report=args.report,
        report_file=report_file,
        dpi=args.dpi,
        title=args.title,
        show_3d_cues=not args.no_3d_cues,
        diagrams=diagrams,
        show_nearby_residues=args.show_nearby
    )

    if mapper.diagram_files:
        say("Analysis complete. Diagrams saved:")
        for style, path in mapper.diagram_files.items():
            say(f"  {style:<8}{path}")
    else:
        say("Analysis complete (no diagrams requested).")

    if args.deltaG:
        result = mapper.interaction_score()
        print("\n--- Empirical interaction score ---")
        print(f"  Score: {result['score']:.2f} (kcal/mol scale, not calibrated)")
        print(f"  {result['interpretation']}")
        print(f"  {result['note']}")
        print("  Breakdown:")
        for itype, info in result['breakdown'].items():
            n = info.get('unique_residues', info.get('count', '?'))
            print(f"    {itype} (n={n}): {info['contribution_kcal_mol']:+.2f}")
        print("----------------------------------")

    if args.csv:
        mapper.export_interactions_csv(
            args.csv if isinstance(args.csv, str) else f"{out_stem}.csv")
    if args.json:
        mapper.export_interactions_json(
            args.json if isinstance(args.json, str) else f"{out_stem}.json",
            include_score=args.deltaG)
    if args.plots:
        mapper.generate_interaction_plots(
            args.plots if isinstance(args.plots, str) else f"{out_stem}_report.png",
            dpi=args.dpi)

    if args.report:
        if os.path.exists(report_file):
            say(f"Interaction report saved to: {report_file}")
        else:
            print(f"Warning: the interaction report could not be written to {report_file}")

    if args.use_3d or args.offline_3d:
        create_pandamap_3d_viz(
            mapper=mapper,
            output_file=args.output_3d or f"{base}_3d_visualization.html",
            width=args.width,
            height=args.height,
            show_surface=not args.no_surface,
            self_contained=args.offline_3d
        )
    return 0


def main():
    """Command-line interface for PandaMap."""
    _make_console_safe()
    from pandamap import __version__
    args = build_parser(__version__).parse_args()
    _configure_logging(args.quiet)
    say = (lambda *a, **k: None) if args.quiet else print

    # Tell people on an old release that a newer one exists.  Importing
    # pandamap has already started this check; here we only make sure a
    # short command still gets to print the notice before it exits.
    from pandamap.utils.versioning import check_for_updates_in_background
    checker = check_for_updates_in_background()

    missing = [p for p in args.structures if not os.path.exists(p)]
    if missing:
        print(f"Error: File not found: {', '.join(missing)}")
        return 1

    try:
        files = expand_structure_paths(args.structures)
        batch = len(args.structures) > 1 or os.path.isdir(args.structures[0])
        if batch and not files:
            print("Error: no structure files (.pdb, .cif, .mmcif, .pdbqt) found in: "
                  + ", ".join(args.structures))
            return 1
        if batch and args.list_ligands:
            for path in files:
                print(f"\n{path}")
                _run_single(args, path, say)
            code = 0
        elif batch:
            code = _run_batch(args, files, say)
        else:
            code = _run_single(args, args.structures[0], say)
    except Exception as e:
        print(f"Error: {e}")
        if args.debug:
            import traceback
            traceback.print_exc()
        elif not isinstance(e, (ValueError, OSError)):
            print("Run again with --debug for the full traceback.")
        return 1

    if checker is not None:
        checker.join(timeout=0.3)
    return code


if __name__ == "__main__":
    sys.exit(main())

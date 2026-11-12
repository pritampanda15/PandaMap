"""
Batch analysis: run PandaMap over many structures (a folder of complexes or
docking poses) and collect the results in two tables.

    batch_summary.csv        one row per structure
    batch_interactions.csv   one row per contact, with a ``structure`` column
"""
import csv
import logging
import os

logger = logging.getLogger("pandamap")

STRUCTURE_EXTENSIONS = ('.pdb', '.cif', '.mmcif', '.pdbqt', '.ent')


def expand_structure_paths(paths):
    """Files as given; folders expanded to the structure files they contain."""
    files = []
    for path in paths:
        if os.path.isdir(path):
            files += sorted(os.path.join(path, name) for name in os.listdir(path)
                            if name.lower().endswith(STRUCTURE_EXTENSIONS))
        else:
            files.append(path)
    return files


def analyze_batch(paths, output_dir='pandamap_batch', ligand_resname=None,
                  diagrams=('modern',), report=False, plots=False, json_export=False,
                  dpi=150, **mapper_options):
    """
    Analyse several structures and write combined tables.

    Parameters
    ----------
    paths : list of str
        Structure files and/or folders containing them.
    output_dir : str
        Where the two combined CSV files go; each structure's own outputs go
        into ``output_dir/<structure name>/``.
    ligand_resname : str, optional
        Ligand residue name used for every structure (default: auto-detect
        per structure).
    diagrams : sequence of str
        Diagrams to draw per structure ('modern', 'classic', '3d'); pass an
        empty tuple to skip images, which is much faster for screening.
    report, plots, json_export : bool
        Also write the text report, the summary figure and the JSON export
        for each structure.
    mapper_options :
        Passed to ``HybridProtLigMapper`` (for example ``strict_halogens``).

    A structure that fails is recorded with its error and the batch goes on.

    Returns
    -------
    dict with ``results`` (one dict per structure), ``summary_csv`` and
    ``interactions_csv``.
    """
    from .core import HybridProtLigMapper

    files = expand_structure_paths(paths)
    if not files:
        raise ValueError("No structure files found in: " + ", ".join(map(str, paths)))
    os.makedirs(output_dir, exist_ok=True)

    results, contact_rows, all_types = [], [], []
    seen = {}
    for number, path in enumerate(files, 1):
        stem = os.path.splitext(os.path.basename(path))[0]
        seen[stem] = seen.get(stem, 0) + 1
        name = stem if seen[stem] == 1 else f"{stem}_{seen[stem]}"   # same name twice
        row = {'structure': name, 'file': path, 'status': 'ok', 'error': ''}
        logger.info("[%d/%d] %s", number, len(files), path)
        try:
            mapper = HybridProtLigMapper(path, ligand_resname=ligand_resname, **mapper_options)
            folder = os.path.join(output_dir, name)
            wants_files = bool(diagrams) or report or plots or json_export
            if wants_files:
                os.makedirs(folder, exist_ok=True)
            base = os.path.join(folder, f"{name}_interactions")
            mapper.run_analysis(output_file=base + ".png", diagrams=tuple(diagrams or ()),
                                generate_report=report, report_file=base + "_report.txt",
                                dpi=dpi)
            if plots:
                mapper.generate_interaction_plots(base + "_report.png", dpi=dpi)
            if json_export:
                mapper.export_interactions_json(base + ".json", include_score=True)

            data = mapper.interactions_as_dict()
            score = mapper.interaction_score()
            row.update({
                'ligand': data['ligand']['name'], 'ligand_chain': data['ligand']['chain'],
                'ligand_resnum': data['ligand']['number'],
                'n_interactions': data['summary']['n_interactions'],
                'n_residues': data['summary']['n_residues'],
                'interaction_score': score['score'],
            })
            for itype, count in data['summary']['by_type'].items():
                row[itype] = count
                if itype not in all_types:
                    all_types.append(itype)
            for c in data['interactions']:
                contact_rows.append({
                    'structure': name, 'ligand': data['ligand']['name'],
                    'interaction_type': c['type'],
                    'protein_residue': c['residue']['name'],
                    'protein_resnum': c['residue']['number'],
                    'protein_chain': c['residue']['chain'],
                    'protein_atom': c['protein_atom'], 'ligand_atom': c['ligand_atom'],
                    'distance_A': round(c['distance'], 2),
                    'solvent_accessible': c['solvent_accessible'],
                    'note': c.get('note', ''),
                })
        except Exception as exc:
            row.update({'status': 'failed', 'error': str(exc)})
            logger.warning("  failed: %s", exc)
        results.append(row)

    summary_csv = os.path.join(output_dir, 'batch_summary.csv')
    fields = ['structure', 'file', 'status', 'ligand', 'ligand_chain', 'ligand_resnum',
              'n_interactions', 'n_residues', 'interaction_score'] + sorted(all_types) + ['error']
    with open(summary_csv, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, restval='')
        writer.writeheader()
        for row in results:
            # a structure without a given contact type has zero of them
            writer.writerow({**{t: 0 for t in all_types if row['status'] == 'ok'}, **row})

    interactions_csv = os.path.join(output_dir, 'batch_interactions.csv')
    with open(interactions_csv, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            'structure', 'ligand', 'interaction_type', 'protein_residue', 'protein_resnum',
            'protein_chain', 'protein_atom', 'ligand_atom', 'distance_A',
            'solvent_accessible', 'note'])
        writer.writeheader()
        writer.writerows(contact_rows)

    failed = sum(1 for r in results if r['status'] != 'ok')
    logger.info("Batch complete: %d structure%s analysed, %d failed.",
                len(results) - failed, '' if len(results) - failed == 1 else 's', failed)
    logger.info("Summary: %s", summary_csv)
    logger.info("All contacts: %s", interactions_csv)
    return {'results': results, 'summary_csv': summary_csv,
            'interactions_csv': interactions_csv}

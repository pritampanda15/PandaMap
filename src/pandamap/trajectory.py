"""
Multi-model analysis: run PandaMap on every model of an NMR ensemble or a
trajectory saved as a multi-MODEL PDB file, and aggregate the results.
"""
import logging
import os
import tempfile
from collections import defaultdict

from .structure_io import MultiFormatParser

logger = logging.getLogger("pandamap")


def analyze_trajectory(structure_file, ligand_resname=None,
                       output_dir='.', dpi=300,
                       generate_reports=False,
                       frame_step=1):
    """
    Analyse every structural model (frame) in a multi-model PDB/mmCIF file and
    aggregate interaction statistics across the ensemble.

    Multi-model PDB files are produced by NMR refinement, conformational
    sampling, or converted from MD trajectories (e.g., via MDAnalysis or
    GROMACS trjconv).  For true MD trajectories in DCD/XTC format, convert
    to multi-model PDB first:
        ``gmx trjconv -f traj.xtc -s topol.tpr -o traj_frames.pdb``
        ``mdconvert traj.dcd -o traj_frames.pdb``

    Parameters
    ----------
    structure_file : str
        Path to a multi-model PDB or mmCIF file.
    ligand_resname : str, optional
        Three-letter residue name of the ligand (auto-detected if None).
    output_dir : str
        Directory for per-frame PNGs and an aggregated summary CSV.
    dpi : int
        Resolution for per-frame PNG images.
    generate_reports : bool
        If True, generate a text report for each frame.
    frame_step : int
        Analyse every Nth frame (1 = every frame, 2 = every other, etc.).

    Returns
    -------
    dict
        'per_frame'   : list of per-frame result dicts (frame_idx, interactions, dG)
        'mean_counts' : dict of interaction_type → mean count across frames
        'occupancy'   : dict of residue_id → fraction of frames it appears
        'summary_csv' : path to the written CSV file
    """
    import csv
    from .core import HybridProtLigMapper

    os.makedirs(output_dir, exist_ok=True)

    # Count models in the file
    parser = MultiFormatParser()
    structure = parser.parse_structure(structure_file)
    models = list(structure.get_models())

    if len(models) == 0:
        raise ValueError(f"No structural models found in {structure_file}")

    if len(models) == 1:
        logger.warning("Warning: Only one model found. "
              "For trajectory analysis supply a multi-model PDB file.")

    # Per-frame results
    per_frame = []
    residue_frame_counts = defaultdict(int)
    interaction_sums = defaultdict(int)
    n_analysed = 0

    # Write a single-model PDB from each model for HybridProtLigMapper
    tmp_pdb = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
    tmp_path = tmp_pdb.name
    tmp_pdb.close()

    for model_idx, model in enumerate(models):
        if model_idx % frame_step != 0:
            continue

        # Write this model to a temporary single-model PDB
        from Bio.PDB import PDBIO, Select

        class SingleModelSelect(Select):
            def __init__(self, target_model_id):
                self._mid = target_model_id
            def accept_model(self, m):
                return 1 if m.id == self._mid else 0

        io = PDBIO()
        io.set_structure(structure)
        io.save(tmp_path, select=SingleModelSelect(model.id))

        # Run PandaMap on this frame
        try:
            mapper = HybridProtLigMapper(tmp_path, ligand_resname=ligand_resname)
            mapper.detect_interactions()
            affinity = mapper.estimate_binding_affinity()

            frame_counts = {k: len(v) for k, v in mapper.interactions.items()}
            for itype, cnt in frame_counts.items():
                interaction_sums[itype] += cnt

            for res_id in mapper.interacting_residues:
                residue_frame_counts[res_id] += 1

            frame_result = {
                'frame_idx': model_idx,
                'model_id': model.id,
                'interaction_counts': frame_counts,
                'dG_estimated': affinity['dG_estimated'],
                'interacting_residues': list(mapper.interacting_residues),
            }

            # Optionally generate visualisation
            png_path = os.path.join(
                output_dir,
                f"{os.path.splitext(os.path.basename(structure_file))[0]}_frame{model_idx:04d}.png"
            )
            try:
                mapper.calculate_realistic_solvent_accessibility()
                mapper.visualize(output_file=png_path, dpi=dpi,
                                 title=f"Frame {model_idx}")
                frame_result['png'] = png_path
            except Exception as viz_err:
                frame_result['png'] = None
                logger.info(f"  Frame {model_idx} visualisation skipped: {viz_err}")

            if generate_reports:
                rpt_path = os.path.join(
                    output_dir,
                    f"{os.path.splitext(os.path.basename(structure_file))[0]}_frame{model_idx:04d}.txt"
                )
                try:
                    if hasattr(mapper, 'generate_interaction_report'):
                        mapper.generate_interaction_report(output_file=rpt_path)
                        frame_result['report'] = rpt_path
                except Exception:
                    pass

            per_frame.append(frame_result)
            n_analysed += 1
            logger.info(f"Frame {model_idx:4d}/{len(models)-1}  "
                  f"score {affinity['dG_estimated']:+.1f}  "
                  f"residues={len(mapper.interacting_residues)}")

        except Exception as frame_err:
            logger.info(f"Frame {model_idx} failed: {frame_err}")

    # Clean up temporary PDB
    try:
        os.unlink(tmp_path)
    except Exception:
        pass

    if n_analysed == 0:
        raise RuntimeError("No frames could be analysed successfully.")

    # Compute aggregate statistics
    mean_counts = {k: round(v / n_analysed, 2)
                   for k, v in interaction_sums.items()}
    total_frames = n_analysed
    occupancy = {res_id: round(cnt / total_frames, 3)
                 for res_id, cnt in residue_frame_counts.items()}

    # Write summary CSV
    csv_path = os.path.join(
        output_dir,
        f"{os.path.splitext(os.path.basename(structure_file))[0]}_trajectory_summary.csv"
    )
    fieldnames = ['frame_idx', 'model_id', 'dG_estimated'] + \
                 list(interaction_sums.keys())
    with open(csv_path, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames,
                                extrasaction='ignore')
        writer.writeheader()
        for fr in per_frame:
            row = {
                'frame_idx': fr['frame_idx'],
                'model_id': fr['model_id'],
                'dG_estimated': fr['dG_estimated'],
            }
            row.update(fr['interaction_counts'])
            writer.writerow(row)

    logger.info(f"\nTrajectory analysis complete: {n_analysed} frames analysed.")
    logger.info(f"Summary CSV written to: {csv_path}")

    return {
        'per_frame': per_frame,
        'mean_counts': mean_counts,
        'occupancy': occupancy,
        'summary_csv': csv_path,
        'n_frames_analysed': n_analysed,
    }

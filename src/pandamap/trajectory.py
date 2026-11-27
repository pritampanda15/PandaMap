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
                       frame_step=1,
                       ligand_chain=None,
                       ligand_resnum=None,
                       diagrams=('modern',),
                       color_scheme=None,
                       strict_halogens=False,
                       strict_hbond_geometry=False,
                       water_angle_check=True):
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
            mapper = HybridProtLigMapper(
                tmp_path, ligand_resname=ligand_resname,
                ligand_chain=ligand_chain, ligand_resnum=ligand_resnum,
                strict_halogens=strict_halogens,
                strict_hbond_geometry=strict_hbond_geometry,
                water_angle_check=water_angle_check,
            )
            mapper.detect_interactions()
            affinity = mapper.estimate_binding_affinity()

            frame_counts = {k: len(v) for k, v in mapper.interactions.items()}
            for itype, cnt in frame_counts.items():
                interaction_sums[itype] += cnt

            # Count actual contacting residues (not mere proximity neighbors)
            contacting = set()
            for hits in mapper.interactions.values():
                for hit in hits:
                    res = hit.get('protein_residue')
                    if res is not None:
                        contacting.add((res.resname, res.id[1]))
            for res_id in contacting:
                residue_frame_counts[res_id] += 1

            frame_result = {
                'frame_idx': model_idx,
                'model_id': model.id,
                'interaction_counts': frame_counts,
                'dG_estimated': affinity['dG_estimated'],
                'interacting_residues': sorted(contacting),
            }

            # Optionally generate visualisation
            if diagrams:
                png_path = os.path.join(
                    output_dir,
                    f"{os.path.splitext(os.path.basename(structure_file))[0]}_frame{model_idx:04d}.png"
                )
                try:
                    mapper.calculate_realistic_solvent_accessibility()
                    mapper.run_analysis(
                        output_file=png_path, dpi=dpi, use_dssp=False,
                        title=f"Frame {model_idx}", diagrams=diagrams,
                        color_scheme=color_scheme,
                    )
                    frame_result['png'] = png_path
                except Exception as viz_err:
                    frame_result['png'] = None
                    logger.info(f"  Frame {model_idx} visualisation skipped: {viz_err}")
            else:
                frame_result['png'] = None

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


def _draw_trajectory_timeline(summary_csv, occupancy_csv, snapshots, output_file,
                              title=None, dpi=200):
    """
    Render a publication-style three-panel MD trajectory overview figure:
      a) Top residue contact occupancy (%) across all frames
      b) Total interactions & hydrogen bonds vs simulation time (ns)
         (vertical markers indicate rendered snapshot timestamps)
      c) Empirical interaction score vs simulation time (ns)
    """
    import csv
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    times_ns, total_ints, hbonds, water_b, scores = [], [], [], [], []
    with open(summary_csv, encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            times_ns.append(float(row['time_ns']))
            total_ints.append(int(row['n_interactions']))
            hbonds.append(int(row.get('hydrogen_bonds', 0)))
            water_b.append(int(row.get('water_bridges', 0)))
            scores.append(float(row['score']))

    if not times_ns:
        return None

    occ_rows = []
    if os.path.isfile(occupancy_csv):
        with open(occupancy_csv, encoding='utf-8') as fh:
            for row in csv.DictReader(fh):
                occ_rows.append((row['label'], float(row['occupancy_pct']), row.get('dominant_types', '')))
    occ_rows = occ_rows[:18]

    t = np.array(times_ns, dtype=float)
    use_frame_axis = float(np.ptp(t)) < 1e-6
    x = np.arange(len(t), dtype=float) if use_frame_axis else t
    xlabel = 'Frame / Pose Index' if use_frame_axis else 'Simulation Time (ns)'
    snap_xs = [float(s['frame_idx']) if use_frame_axis else float(s['time_ns']) for s in snapshots]

    fig = plt.figure(figsize=(11.5, 7.6), facecolor='white')
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1.35], height_ratios=[1.0, 1.0],
                          wspace=0.30, hspace=0.32, left=0.11, right=0.97, top=0.90, bottom=0.09)

    ax_occ = fig.add_subplot(gs[:, 0])
    ax_cnt = fig.add_subplot(gs[0, 1])
    ax_scr = fig.add_subplot(gs[1, 1], sharex=ax_cnt)

    # Panel a: Residue contact occupancy
    if occ_rows:
        labels = [r[0] for r in occ_rows]
        vals = [r[1] for r in occ_rows]
        ys = np.arange(len(labels))
        colors = ['#00A087' if v >= 70 else ('#3C5488' if v >= 35 else '#8491B4') for v in vals]
        ax_occ.barh(ys, vals, height=0.68, color=colors, edgecolor='none')
        for y, v in zip(ys, vals):
            ax_occ.text(min(v + 1.5, 102), y, f"{v:.1f}%", va='center', fontsize=8, color='#333333')
        ax_occ.set_yticks(ys)
        ax_occ.set_yticklabels(labels, fontsize=8.5)
        ax_occ.invert_yaxis()
        ax_occ.set_xlim(0, 114)
        ax_occ.set_xlabel('Contact Occupancy (%)', fontsize=9.5)
    else:
        ax_occ.text(0.5, 0.5, 'No contacts detected', ha='center', va='center', transform=ax_occ.transAxes)
    ax_occ.set_title('a  Top Residue Occupancy (All Frames)', loc='left', fontsize=10.5, fontweight='bold')
    ax_occ.spines['top'].set_visible(False)
    ax_occ.spines['right'].set_visible(False)
    ax_occ.grid(axis='x', color='#EAEAEA', lw=0.6)
    ax_occ.set_axisbelow(True)

    # Panel b: Contact counts over time
    ax_cnt.plot(x, total_ints, color='#3C5488', lw=1.2, alpha=0.85, label='Total contacts')
    ax_cnt.plot(x, hbonds, color='#00A087', lw=1.4, label='Hydrogen bonds')
    if any(water_b):
        ax_cnt.plot(x, water_b, color='#4DBBD5', lw=1.2, ls='--', label='Water bridges')
    for sx in snap_xs:
        ax_cnt.axvline(sx, color='#E64B35', lw=0.7, ls=':', alpha=0.55)
    ax_cnt.set_ylabel('Contacts (n)', fontsize=9.5)
    ax_cnt.set_title('b  Contacts Along Trajectory', loc='left', fontsize=10.5, fontweight='bold')
    ax_cnt.legend(frameon=False, fontsize=8.5, loc='upper right')
    ax_cnt.spines['top'].set_visible(False)
    ax_cnt.spines['right'].set_visible(False)
    ax_cnt.grid(color='#EFEFEF', lw=0.6)

    # Panel c: Interaction score over time
    sc = np.array(scores, dtype=float)
    ax_scr.plot(x, sc, color='#B09C85', lw=0.9, alpha=0.75, label='Per-frame score')
    if len(sc) >= 15:
        win = max(5, len(sc) // 25)
        kernel = np.ones(win) / win
        smooth = np.convolve(sc, kernel, mode='valid')
        offset = win // 2
        ax_scr.plot(x[offset:offset + len(smooth)], smooth, color='#E64B35', lw=1.8,
                    label=f'Rolling mean ({win} frames)')
    for sx in snap_xs:
        ax_scr.axvline(sx, color='#E64B35', lw=0.7, ls=':', alpha=0.55)
    ax_scr.set_xlabel(xlabel, fontsize=9.5)
    ax_scr.set_ylabel('Interaction Score\n(kcal/mol scale)', fontsize=9.5)
    ax_scr.set_title('c  Empirical Interaction Score', loc='left', fontsize=10.5, fontweight='bold')
    ax_scr.legend(frameon=False, fontsize=8.5, loc='upper right')
    ax_scr.spines['top'].set_visible(False)
    ax_scr.spines['right'].set_visible(False)
    ax_scr.grid(color='#EFEFEF', lw=0.6)

    suptitle = title or f"PandaMap Trajectory Analysis ({len(times_ns):,} frames)"
    if snap_xs:
        suptitle += f" · {len(snap_xs)} rendered snapshots (dotted red)"
    fig.suptitle(suptitle, fontsize=12.5, fontweight='bold', y=0.97)
    fig.savefig(output_file, dpi=dpi, facecolor='white')
    plt.close(fig)
    logger.info("Trajectory timeline saved to: %s", output_file)
    return output_file


def analyze_md_trajectory(
    topology_file,
    traj_file=None,
    output_dir='pandamap_traj',
    ligand_resname=None,
    ligand_chain=None,
    ligand_resnum=None,
    every=None,
    stride=None,
    start=None,
    end=None,
    dt_ps=0.0,
    diagrams=('modern',),
    report=False,
    plots=False,
    csv_export=False,
    json_export=False,
    use_3d=False,
    offline_3d=False,
    dpi=150,
    color_scheme=None,
    show_nearby=False,
    show_3d_cues=True,
    strict_halogens=False,
    strict_hbond_geometry=False,
    water_angle_check=True,
    write_fingerprints=True,
    write_timeline=True,
    threads=0,
    water_shell_a=10.0,
):
    """
    Fast Rust-accelerated MD trajectory (PDB/GRO + XTC) or multi-pose ensemble
    (multi-MODEL PDB / multi-pose PDBQT) analysis.

    Stage 1 (Rust engine):
        Uses 3D spatial hashing + multi-core frame parallelism to fingerprint
        every frame in the trajectory and write:
          - ``trajectory_summary.csv``
          - ``residue_occupancy.csv``
          - ``trajectory_fingerprints.csv`` (ML-ready feature matrix)
          - ``snapshots/<label>/<label>.pdb`` at the requested interval
            (e.g. ``every="10ns"`` or ``every="1ns"``)

    Stage 2 (PandaMap Python renderers):
        Runs full PandaMap visualization/reporting on each selected snapshot
        inside ``<output_dir>/snapshots/<label>/`` and draws
        ``<output_dir>/trajectory_timeline.png``.
    """
    from .core import HybridProtLigMapper
    from .create_3d_view import create_pandamap_3d_viz
    from .rust_engine import run_rust_trajectory

    os.makedirs(output_dir, exist_ok=True)
    diagrams = tuple(diagrams or ())
    wants_snapshot_outputs = bool(
        diagrams or report or plots or csv_export or json_export or use_3d or offline_3d
    )

    manifest = run_rust_trajectory(
        topology_file=topology_file,
        traj_file=traj_file,
        output_dir=output_dir,
        ligand_resname=ligand_resname,
        ligand_chain=ligand_chain,
        ligand_resnum=ligand_resnum,
        every=every,
        stride=stride,
        start=start,
        end=end,
        dt_ps=dt_ps,
        water_shell_a=water_shell_a,
        threads=threads,
        strict_halogens=strict_halogens,
        strict_hbond_geometry=strict_hbond_geometry,
        water_angle_check=water_angle_check,
        write_fingerprints=write_fingerprints,
        export_snapshots=wants_snapshot_outputs,
    )

    lig_name = manifest['ligand']['name']
    lig_chain = manifest['ligand']['chain'] or None
    lig_resnum = manifest['ligand']['resnum']
    rendered_snapshots = []

    if wants_snapshot_outputs and manifest.get('snapshots'):
        total_snaps = len(manifest['snapshots'])
        for idx, snap in enumerate(manifest['snapshots'], 1):
            label = snap['label']
            pdb_path = snap['pdb_path']
            snap_dir = os.path.dirname(pdb_path)
            base = os.path.join(snap_dir, f"{label}_interactions")
            time_ns = float(snap['time_ns'])
            frame_idx = int(snap['frame_idx'])
            snap_title = (
                f"{os.path.splitext(os.path.basename(topology_file))[0]} · {lig_name} · "
                f"t = {time_ns:.2f} ns (frame {frame_idx})"
                if traj_file or dt_ps > 0
                else f"{os.path.splitext(os.path.basename(topology_file))[0]} · {lig_name} · frame {frame_idx}"
            )
            logger.info("[%d/%d] Rendering snapshot %s (t=%.2f ns)...", idx, total_snaps, label, time_ns)
            snap_out = dict(snap)
            snap_out['outputs'] = {}
            try:
                mapper = HybridProtLigMapper(
                    pdb_path,
                    ligand_resname=lig_name,
                    ligand_chain=lig_chain,
                    ligand_resnum=lig_resnum,
                    strict_halogens=strict_halogens,
                    strict_hbond_geometry=strict_hbond_geometry,
                    water_angle_check=water_angle_check,
                )
                mapper.run_analysis(
                    output_file=base + ".png",
                    use_dssp=False,
                    generate_report=report,
                    report_file=base + "_report.txt",
                    dpi=dpi,
                    title=snap_title,
                    show_3d_cues=show_3d_cues,
                    diagrams=diagrams,
                    show_nearby_residues=show_nearby,
                    color_scheme=color_scheme,
                )
                snap_out['outputs']['diagrams'] = dict(mapper.diagram_files)
                if report:
                    snap_out['outputs']['report'] = base + "_report.txt"
                if csv_export:
                    snap_out['outputs']['csv'] = mapper.export_interactions_csv(base + ".csv")
                if json_export:
                    snap_out['outputs']['json'] = mapper.export_interactions_json(
                        base + ".json", include_score=True
                    )
                if plots:
                    p = mapper.generate_interaction_plots(base + "_report.png", dpi=dpi)
                    if p:
                        snap_out['outputs']['plots'] = p
                if use_3d or offline_3d:
                    html_out = os.path.join(snap_dir, f"{label}_3d_visualization.html")
                    snap_out['outputs']['html_3d'] = create_pandamap_3d_viz(
                        mapper=mapper,
                        output_file=html_out,
                        self_contained=offline_3d,
                    )
                snap_out['status'] = 'ok'
            except Exception as exc:
                snap_out['status'] = 'failed'
                snap_out['error'] = str(exc)
                logger.warning("  Snapshot %s failed: %s", label, exc)
            rendered_snapshots.append(snap_out)

    timeline_png = None
    if write_timeline:
        try:
            timeline_png = _draw_trajectory_timeline(
                manifest['summary_csv'],
                manifest['occupancy_csv'],
                manifest.get('snapshots', []),
                os.path.join(output_dir, 'trajectory_timeline.png'),
                title=f"PandaMap MD Analysis: {os.path.basename(topology_file)} · {lig_name}",
                dpi=max(150, min(dpi, 300)),
            )
        except Exception as exc:
            logger.warning("Could not draw trajectory timeline: %s", exc)

    return {
        'engine': manifest['engine'],
        'n_frames': manifest['n_frames'],
        'elapsed_sec': manifest['elapsed_sec'],
        'frames_per_sec': manifest['frames_per_sec'],
        'ligand': manifest['ligand'],
        'summary_csv': manifest['summary_csv'],
        'occupancy_csv': manifest['occupancy_csv'],
        'fingerprints_csv': manifest.get('fingerprints_csv'),
        'manifest_json': manifest['manifest_json'],
        'timeline_png': timeline_png,
        'snapshots': rendered_snapshots,
        'output_dir': output_dir,
    }

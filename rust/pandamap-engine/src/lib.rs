//! PandaMap Rust interaction fingerprinting and trajectory engine.

pub mod geometry;
pub mod interactions;
pub mod spatial_hash;
pub mod topology;
pub mod xtc;

use std::collections::HashSet;
use std::io;
use std::path::Path;
use std::time::Instant;

use interactions::{
    analyze_models_parallel, analyze_xtc_batch_parallel, export_snapshot_pdbs,
    select_snapshot_indices, write_analysis_outputs, DetectConfig, ExportedSnapshot,
    FrameAnalysis, SnapshotSpec,
};
use topology::{load_topology, prepare_system};
use xtc::{RawXtcFrame, XtcReader};

#[derive(Debug, Clone)]
pub struct RunOptions {
    pub topology_path: String,
    pub traj_path: Option<String>,
    pub out_dir: String,
    pub ligand_resname: Option<String>,
    pub ligand_chain: Option<String>,
    pub ligand_resnum: Option<i32>,
    pub every_ps: Option<f64>,
    pub stride: Option<usize>,
    pub start_ps: Option<f64>,
    pub end_ps: Option<f64>,
    pub dt_ps: f32,
    pub water_shell_a: f64,
    pub write_fingerprints: bool,
    pub export_snapshots: bool,
    pub detect: DetectConfig,
}

pub struct RunResult {
    pub n_frames: usize,
    pub elapsed_sec: f64,
    pub snapshots: Vec<ExportedSnapshot>,
}

/// Execute full trajectory (PDB/GRO + XTC) or multi-model (PDB/PDBQT) analysis.
pub fn run_trajectory_analysis(opts: &RunOptions) -> io::Result<RunResult> {
    let t0 = Instant::now();
    let topo = load_topology(&opts.topology_path)?;
    let natoms_topo = topo.atoms.len();
    let sys = prepare_system(
        topo.atoms,
        opts.ligand_resname.as_deref(),
        opts.ligand_chain.as_deref(),
        opts.ligand_resnum,
    )?;

    let out_dir = Path::new(&opts.out_dir);
    let snap_spec = SnapshotSpec {
        every_ps: opts.every_ps,
        stride: opts.stride,
        start_ps: opts.start_ps,
        end_ps: opts.end_ps,
        water_shell_a: opts.water_shell_a,
    };

    if let Some(ref xtc_path) = opts.traj_path {
        // Streaming chunked parallel XTC pass
        let mut reader = XtcReader::open(xtc_path)?;
        let batch_size = 512usize;
        let mut all_results: Vec<FrameAnalysis> = Vec::new();
        let mut raw_batch: Vec<RawXtcFrame> = Vec::with_capacity(batch_size);

        loop {
            raw_batch.clear();
            for _ in 0..batch_size {
                match reader.read_next_raw()? {
                    Some(raw) => {
                        if raw.natoms != natoms_topo {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "Atom count mismatch: topology has {} atoms but XTC frame {} has {}",
                                    natoms_topo, raw.frame_idx, raw.natoms
                                ),
                            ));
                        }
                        raw_batch.push(raw);
                    }
                    None => break,
                }
            }
            if raw_batch.is_empty() {
                break;
            }
            let analyzed = analyze_xtc_batch_parallel(&sys, &raw_batch, &opts.detect)?;
            all_results.extend(analyzed);
        }

        if all_results.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("No frames found in XTC trajectory: {}", xtc_path),
            ));
        }

        let snapshots = if opts.export_snapshots {
            let chosen_indices = select_snapshot_indices(&all_results, &snap_spec);
            let chosen_set: HashSet<usize> = chosen_indices.iter().copied().collect();
            let mut decoded_selected: Vec<(usize, &FrameAnalysis, Vec<[f32; 3]>)> =
                Vec::with_capacity(chosen_set.len());

            if !chosen_set.is_empty() {
                // Second fast sequential pass skipping non-selected frames without decompression!
                let mut r2 = XtcReader::open(xtc_path)?;
                let max_needed = *chosen_indices.iter().max().unwrap();
                for idx in 0..=max_needed {
                    if chosen_set.contains(&idx) {
                        if let Some(raw) = r2.read_next_raw()? {
                            let dec = raw.decode()?;
                            decoded_selected.push((idx, &all_results[idx], dec.coords_angstrom));
                        }
                    } else if r2.skip_next()?.is_none() {
                        break;
                    }
                }
            }
            export_snapshot_pdbs(out_dir, &sys, &decoded_selected, opts.water_shell_a, true)?
        } else {
            Vec::new()
        };

        let elapsed = t0.elapsed().as_secs_f64();
        write_analysis_outputs(
            out_dir,
            &sys,
            &all_results,
            &snapshots,
            elapsed,
            opts.write_fingerprints,
        )?;
        Ok(RunResult {
            n_frames: all_results.len(),
            elapsed_sec: elapsed,
            snapshots,
        })
    } else {
        // Multi-MODEL PDB or multi-pose PDBQT
        let models = topo.models;
        let all_results = analyze_models_parallel(&sys, &models, opts.dt_ps, &opts.detect);
        let use_time = opts.dt_ps > 0.0 && opts.every_ps.is_some();
        let snapshots = if opts.export_snapshots {
            let chosen = select_snapshot_indices(&all_results, &snap_spec);
            let selected: Vec<(usize, &FrameAnalysis, Vec<[f32; 3]>)> = chosen
                .into_iter()
                .map(|i| (i, &all_results[i], models[i].clone()))
                .collect();
            export_snapshot_pdbs(out_dir, &sys, &selected, opts.water_shell_a, use_time)?
        } else {
            Vec::new()
        };
        let elapsed = t0.elapsed().as_secs_f64();
        write_analysis_outputs(
            out_dir,
            &sys,
            &all_results,
            &snapshots,
            elapsed,
            opts.write_fingerprints,
        )?;
        Ok(RunResult {
            n_frames: all_results.len(),
            elapsed_sec: elapsed,
            snapshots,
        })
    }
}

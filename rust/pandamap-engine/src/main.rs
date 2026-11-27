//! CLI binary for PandaMap's Rust trajectory and docking interaction engine.

use std::env;
use std::io;
use std::process;

use pandamap_engine::interactions::DetectConfig;
use pandamap_engine::topology::load_topology;
use pandamap_engine::xtc::XtcWriter;
use pandamap_engine::{run_trajectory_analysis, RunOptions};

fn print_usage() {
    eprintln!(
        "pandamap-engine: high-performance Rust trajectory & pose engine for PandaMap\n\
         \n\
         Usage:\n\
           pandamap-engine analyze --topology <PDB|GRO|PDBQT> [options]\n\
           pandamap-engine make-xtc --topology <PDB|GRO> --output <OUT.xtc> [--frames N] [--dt-ps PS]\n\
         \n\
         Options for analyze:\n\
           --traj <FILE.xtc>        GROMACS XTC trajectory file\n\
           --outdir <DIR>           Output directory (default: pandamap_traj)\n\
           --ligand <NAME>          Ligand residue name\n\
           --chain <ID>             Ligand chain ID\n\
           --resnum <N>             Ligand residue number\n\
           --every-ps <PS>          Snapshot interval in picoseconds (e.g. 1000 = 1ns, 10000 = 10ns)\n\
           --stride <N>             Snapshot frame stride (used when --every-ps is not set)\n\
           --start-ps <PS>          Start time window in picoseconds\n\
           --end-ps <PS>            End time window in picoseconds\n\
           --dt-ps <PS>             Timestep for multi-MODEL PDB/PDBQT without timestamps (default: 0)\n\
           --water-shell <A>        Water shell radius around ligand in snapshot PDBs (default: 10.0)\n\
           --threads <N>            Worker threads (default: 0 = auto)\n\
           --strict-halogens        Restrict halogens to Cl/Br/I\n\
           --strict-hbond-geometry  Apply heavy-atom H-bond proxy when H absent\n\
           --water-distance-only    Skip 71-140 deg water-bridge angle filter\n\
           --no-fingerprints        Skip writing trajectory_fingerprints.csv\n\
           --no-snapshots           Skip writing snapshot PDBs\n"
    );
}

fn run_make_xtc(args: &[String]) -> io::Result<()> {
    let mut topo_path = String::new();
    let mut out_path = String::new();
    let mut n_frames: usize = 100;
    let mut dt_ps: f32 = 100.0;
    let mut jitter: f32 = 0.25;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--topology" | "-s" => {
                i += 1;
                topo_path = args.get(i).cloned().unwrap_or_default();
            }
            "--output" | "-o" => {
                i += 1;
                out_path = args.get(i).cloned().unwrap_or_default();
            }
            "--frames" => {
                i += 1;
                n_frames = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(100);
            }
            "--dt-ps" => {
                i += 1;
                dt_ps = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(100.0);
            }
            "--jitter" => {
                i += 1;
                jitter = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(0.25);
            }
            other => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("Unknown make-xtc option: {}", other),
                ));
            }
        }
        i += 1;
    }

    if topo_path.is_empty() || out_path.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "make-xtc requires --topology and --output",
        ));
    }

    let topo = load_topology(&topo_path)?;
    let base: Vec<[f32; 3]> = topo
        .atoms
        .iter()
        .map(|a| [a.coord.x as f32, a.coord.y as f32, a.coord.z as f32])
        .collect();
    let box_a = [[80.0, 0.0, 0.0], [0.0, 80.0, 0.0], [0.0, 0.0, 80.0]];

    let mut writer = XtcWriter::create(&out_path)?;
    let mut frame_coords = base.clone();
    for f in 0..n_frames {
        let phase = (f as f32) * 0.15;
        for (idx, (dst, src)) in frame_coords.iter_mut().zip(base.iter()).enumerate() {
            let a = &topo.atoms[idx];
            // Gently oscillate ligand atoms relative to protein to simulate MD fluctuations
            let amp = if a.is_hetatm { jitter } else { jitter * 0.25 };
            let fi = (idx as f32) * 0.37;
            dst[0] = src[0] + amp * (phase + fi).sin();
            dst[1] = src[1] + amp * (phase * 0.8 + fi * 1.3).cos();
            dst[2] = src[2] + amp * (phase * 0.6 + fi * 0.7).sin();
        }
        writer.write_frame((f * 500) as i32, f as f32 * dt_ps, box_a, &frame_coords, 1000.0)?;
    }
    writer.flush()?;
    println!(
        "Wrote {} frames ({} atoms, dt={} ps) to {}",
        n_frames,
        base.len(),
        dt_ps,
        out_path
    );
    Ok(())
}

fn run_analyze(args: &[String]) -> io::Result<()> {
    let mut opts = RunOptions {
        topology_path: String::new(),
        traj_path: None,
        out_dir: "pandamap_traj".to_string(),
        ligand_resname: None,
        ligand_chain: None,
        ligand_resnum: None,
        every_ps: None,
        stride: None,
        start_ps: None,
        end_ps: None,
        dt_ps: 0.0,
        water_shell_a: 10.0,
        write_fingerprints: true,
        export_snapshots: true,
        detect: DetectConfig::default(),
    };

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--topology" | "-s" => {
                i += 1;
                opts.topology_path = args.get(i).cloned().unwrap_or_default();
            }
            "--traj" | "-x" => {
                i += 1;
                opts.traj_path = args.get(i).cloned();
            }
            "--outdir" | "-O" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    opts.out_dir = v.clone();
                }
            }
            "--ligand" | "-l" => {
                i += 1;
                opts.ligand_resname = args.get(i).cloned();
            }
            "--chain" => {
                i += 1;
                opts.ligand_chain = args.get(i).cloned();
            }
            "--resnum" => {
                i += 1;
                opts.ligand_resnum = args.get(i).and_then(|s| s.parse().ok());
            }
            "--every-ps" => {
                i += 1;
                opts.every_ps = args.get(i).and_then(|s| s.parse().ok());
            }
            "--stride" => {
                i += 1;
                opts.stride = args.get(i).and_then(|s| s.parse().ok());
            }
            "--start-ps" => {
                i += 1;
                opts.start_ps = args.get(i).and_then(|s| s.parse().ok());
            }
            "--end-ps" => {
                i += 1;
                opts.end_ps = args.get(i).and_then(|s| s.parse().ok());
            }
            "--dt-ps" => {
                i += 1;
                opts.dt_ps = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            }
            "--water-shell" => {
                i += 1;
                opts.water_shell_a = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(10.0);
            }
            "--threads" => {
                i += 1;
                opts.detect.threads = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(0);
            }
            "--strict-halogens" => opts.detect.strict_halogens = true,
            "--strict-hbond-geometry" => opts.detect.strict_hbond_geometry = true,
            "--water-distance-only" => opts.detect.water_angle_check = false,
            "--no-fingerprints" => opts.write_fingerprints = false,
            "--no-snapshots" => opts.export_snapshots = false,
            other => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("Unknown analyze option: {}", other),
                ));
            }
        }
        i += 1;
    }

    if opts.topology_path.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "analyze requires --topology <PDB|GRO|PDBQT>",
        ));
    }

    let res = run_trajectory_analysis(&opts)?;
    let fps = if res.elapsed_sec > 0.0 {
        res.n_frames as f64 / res.elapsed_sec
    } else {
        0.0
    };
    println!(
        "Rust engine analyzed {} frames in {:.3}s ({:.0} frames/s); exported {} snapshot(s) -> {}",
        res.n_frames,
        res.elapsed_sec,
        fps,
        res.snapshots.len(),
        opts.out_dir
    );
    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() || matches!(args[0].as_str(), "-h" | "--help") {
        print_usage();
        process::exit(0);
    }

    let result = match args[0].as_str() {
        "analyze" => run_analyze(&args[1..]),
        "make-xtc" => run_make_xtc(&args[1..]),
        other => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Unknown subcommand '{}'. Use 'analyze' or 'make-xtc'.", other),
        )),
    };

    if let Err(err) = result {
        eprintln!("pandamap-engine error: {}", err);
        process::exit(1);
    }
}

"""
Python bridge to PandaMap's compiled Rust trajectory & pose engine
(``pandamap-engine``).

Provides:
- Time string parsing (``"10ns"`` -> ``10000.0`` ps, ``"500ps"`` -> ``500.0`` ps)
- Discovery (and optional on-demand Cargo build) of the ``pandamap-engine`` binary
- Clean subprocess wrapper returning the parsed ``engine_manifest.json`` dict
"""
import json
import logging
import os
import re
import shutil
import subprocess
import sys

logger = logging.getLogger("pandamap")

_TIME_RE = re.compile(r"^\s*([0-9]+(?:\.[0-9]+)?)\s*(fs|ps|ns|us|µs|μs|ms)?\s*$", re.IGNORECASE)
_TIME_TO_PS = {
    "fs": 1e-3,
    "ps": 1.0,
    "ns": 1e3,
    "us": 1e6,
    "µs": 1e6,
    "μs": 1e6,
    "ms": 1e9,
}


def parse_time_ps(value):
    """
    Parse a human-readable MD time value into picoseconds.

    Examples
    --------
    >>> parse_time_ps("10ns")
    10000.0
    >>> parse_time_ps("1ns")
    1000.0
    >>> parse_time_ps("500ps")
    500.0
    >>> parse_time_ps(2500)
    2500.0
    """
    if value is None:
        return None
    if isinstance(value, (int, float)):
        val = float(value)
        if val < 0:
            raise ValueError(f"Time value must be non-negative, got {value!r}")
        return val
    m = _TIME_RE.match(str(value))
    if not m:
        raise ValueError(
            f"Invalid time specification {value!r}. Use e.g. '1ns', '10ns', '500ps', or a number in ps."
        )
    num = float(m.group(1))
    unit = (m.group(2) or "ps").lower()
    return num * _TIME_TO_PS[unit]


def _repo_root():
    """Return repository root when running from or inside a source checkout, else None."""
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.abspath(os.path.join(here, "..", "..")),
        os.path.abspath(os.getcwd()),
    ]
    # Also walk up from cwd (covers running pytest/CLI from a subdirectory of a checkout)
    cur = os.path.abspath(os.getcwd())
    for _ in range(4):
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        candidates.append(parent)
        cur = parent

    for candidate in candidates:
        if os.path.isfile(os.path.join(candidate, "rust", "pandamap-engine", "Cargo.toml")):
            return candidate
    return None


def find_rust_engine(auto_build=True):
    """
    Locate the compiled ``pandamap-engine`` binary.

    Search order:
    1. ``PANDAMAP_RUST_ENGINE`` environment variable
    2. Packaged binary under ``pandamap/bin/pandamap-engine``
    3. Source checkout release build ``rust/pandamap-engine/target/release/pandamap-engine``
    4. ``pandamap-engine`` on ``PATH``
    5. Optional automatic ``cargo build --release`` when in a source checkout
    """
    exe_name = "pandamap-engine.exe" if sys.platform == "win32" else "pandamap-engine"

    env_path = os.environ.get("PANDAMAP_RUST_ENGINE")
    if env_path and os.path.isfile(env_path) and os.access(env_path, os.X_OK):
        return env_path

    pkg_bin = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bin", exe_name)
    if os.path.isfile(pkg_bin) and os.access(pkg_bin, os.X_OK):
        return pkg_bin

    root = _repo_root()
    if root:
        rel_bin = os.path.join(root, "rust", "pandamap-engine", "target", "release", exe_name)
        if os.path.isfile(rel_bin) and os.access(rel_bin, os.X_OK):
            return rel_bin

    on_path = shutil.which(exe_name)
    if on_path:
        return on_path

    if auto_build and root and shutil.which("cargo"):
        return build_rust_engine(root)

    return None


def build_rust_engine(repo_root=None):
    """Compile ``rust/pandamap-engine`` in release mode and return the binary path."""
    root = repo_root or _repo_root()
    if not root:
        raise RuntimeError("Cannot build pandamap-engine: rust/pandamap-engine source not found.")
    manifest = os.path.join(root, "rust", "pandamap-engine", "Cargo.toml")
    cargo = shutil.which("cargo")
    if not cargo:
        raise RuntimeError(
            "Rust toolchain ('cargo') not found on PATH. Install Rust via https://rustup.rs "
            "or install a prebuilt pandamap wheel with pandamap-engine."
        )
    logger.info("Building PandaMap Rust engine (one-time release compile)...")
    proc = subprocess.run(
        [cargo, "build", "--release", "--manifest-path", manifest],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"Failed to build pandamap-engine:\n{proc.stderr}")
    exe_name = "pandamap-engine.exe" if sys.platform == "win32" else "pandamap-engine"
    bin_path = os.path.join(root, "rust", "pandamap-engine", "target", "release", exe_name)
    if not os.path.isfile(bin_path):
        raise RuntimeError(f"Compiled binary not found at expected path: {bin_path}")
    return bin_path


def run_rust_trajectory(
    topology_file,
    traj_file=None,
    output_dir="pandamap_traj",
    ligand_resname=None,
    ligand_chain=None,
    ligand_resnum=None,
    every=None,
    stride=None,
    start=None,
    end=None,
    dt_ps=0.0,
    water_shell_a=10.0,
    threads=0,
    strict_halogens=False,
    strict_hbond_geometry=False,
    water_angle_check=True,
    write_fingerprints=True,
    export_snapshots=True,
):
    """
    Execute the compiled Rust engine on a topology + XTC trajectory (or
    multi-MODEL PDB / multi-pose PDBQT) and return the parsed manifest dict.
    """
    exe = find_rust_engine(auto_build=True)
    if not exe:
        raise RuntimeError(
            "PandaMap Rust engine binary 'pandamap-engine' was not found and 'cargo' is unavailable. "
            "Build it with: cargo build --release --manifest-path rust/pandamap-engine/Cargo.toml"
        )

    cmd = [exe, "analyze", "--topology", str(topology_file), "--outdir", str(output_dir)]
    if traj_file:
        cmd += ["--traj", str(traj_file)]
    if ligand_resname:
        cmd += ["--ligand", str(ligand_resname)]
    if ligand_chain:
        cmd += ["--chain", str(ligand_chain)]
    if ligand_resnum is not None:
        cmd += ["--resnum", str(int(ligand_resnum))]

    every_ps = parse_time_ps(every)
    start_ps = parse_time_ps(start)
    end_ps = parse_time_ps(end)
    if every_ps is not None:
        cmd += ["--every-ps", str(every_ps)]
    if stride is not None:
        cmd += ["--stride", str(int(stride))]
    if start_ps is not None:
        cmd += ["--start-ps", str(start_ps)]
    if end_ps is not None:
        cmd += ["--end-ps", str(end_ps)]
    if dt_ps:
        cmd += ["--dt-ps", str(float(dt_ps))]
    if water_shell_a is not None:
        cmd += ["--water-shell", str(float(water_shell_a))]
    if threads:
        cmd += ["--threads", str(int(threads))]
    if strict_halogens:
        cmd.append("--strict-halogens")
    if strict_hbond_geometry:
        cmd.append("--strict-hbond-geometry")
    if not water_angle_check:
        cmd.append("--water-distance-only")
    if not write_fingerprints:
        cmd.append("--no-fingerprints")
    if not export_snapshots:
        cmd.append("--no-snapshots")

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        msg = (proc.stderr or proc.stdout).strip()
        raise RuntimeError(f"Rust engine failed: {msg}")
    if proc.stdout.strip():
        logger.info(proc.stdout.strip())

    manifest_path = os.path.join(output_dir, "engine_manifest.json")
    with open(manifest_path, encoding="utf-8") as fh:
        manifest = json.load(fh)
    manifest["manifest_json"] = manifest_path
    return manifest

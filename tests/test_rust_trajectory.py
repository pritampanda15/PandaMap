"""
Tests for PandaMap's Rust spatial-hash + parallel XTC/PDBQT/GRO trajectory
engine and its CLI integration (--traj, --every 1ns/10ns, --start/--end,
snapshots folder, ML fingerprints, and timeline figure).
"""
import csv
import json
import os
import subprocess
import sys

import pytest

from pandamap import HybridProtLigMapper, analyze_md_trajectory
from pandamap.rust_engine import find_rust_engine, parse_time_ps


EXAMPLE_1M17 = os.path.join(os.path.dirname(__file__), "..", "example", "1m17.pdb")
needs_1m17 = pytest.mark.skipif(
    not os.path.exists(EXAMPLE_1M17), reason="example/1m17.pdb not available"
)


def test_parse_time_ps():
    assert parse_time_ps("1ns") == pytest.approx(1000.0)
    assert parse_time_ps("10ns") == pytest.approx(10000.0)
    assert parse_time_ps("500ps") == pytest.approx(500.0)
    assert parse_time_ps("0.5 ns") == pytest.approx(500.0)
    assert parse_time_ps("1us") == pytest.approx(1_000_000.0)
    assert parse_time_ps(2500) == pytest.approx(2500.0)
    assert parse_time_ps(None) is None
    with pytest.raises(ValueError):
        parse_time_ps("10lightyears")
    with pytest.raises(ValueError):
        parse_time_ps(-5)


@needs_1m17
def test_rust_engine_matches_python_mapper_on_1m17(tmp_path):
    """Rust spatial-hash detector must match Python HybridProtLigMapper counts
    and score on EGFR + erlotinib (1M17, including water bridges)."""
    py_mapper = HybridProtLigMapper(EXAMPLE_1M17, ligand_resname="AQ4")
    py_mapper.detect_interactions()
    py_counts = {k: len(v) for k, v in py_mapper.interactions.items() if v}
    py_score = py_mapper.interaction_score()["score"]

    res = analyze_md_trajectory(
        topology_file=EXAMPLE_1M17,
        output_dir=str(tmp_path / "parity"),
        ligand_resname="AQ4",
        diagrams=(),
        write_timeline=False,
    )
    assert res["n_frames"] == 1
    with open(res["summary_csv"], encoding="utf-8") as fh:
        row = next(csv.DictReader(fh))

    for itype, count in py_counts.items():
        assert int(row[itype]) == count, f"mismatch for {itype}: rust={row[itype]} vs py={count}"
    assert float(row["score"]) == pytest.approx(py_score, abs=0.05)


@needs_1m17
def test_cli_xtc_every_10ns_creates_snapshot_folders_and_outputs(tmp_path, monkeypatch):
    """End-to-end CLI test: 20 ns compressed XTC (201 frames @ 100 ps) with
    --every 10ns renders snapshots at 0, 10, 20 ns plus summary/occupancy/ML
    fingerprints and timeline PNG."""
    exe = find_rust_engine(auto_build=False)
    assert exe is not None, "pandamap-engine binary must be built"

    xtc_path = tmp_path / "prod_20ns.xtc"
    subprocess.run(
        [
            exe, "make-xtc",
            "--topology", EXAMPLE_1M17,
            "--output", str(xtc_path),
            "--frames", "201",
            "--dt-ps", "100.0",
            "--jitter", "0.20",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert xtc_path.stat().st_size > 1000

    outdir = tmp_path / "md_10ns"
    from pandamap.cli import main

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "pandamap",
            EXAMPLE_1M17,
            "--traj", str(xtc_path),
            "--ligand", "AQ4",
            "--every", "10ns",
            "--outdir", str(outdir),
            "--diagrams", "modern",
            "--report",
            "--json",
            "--dpi", "60",
            "--quiet",
        ],
    )
    assert main() == 0

    # Verify full-trajectory tables and timeline
    with open(outdir / "trajectory_summary.csv", encoding="utf-8") as fh:
        summary_rows = list(csv.DictReader(fh))
    assert len(summary_rows) == 201
    assert float(summary_rows[0]["time_ns"]) == pytest.approx(0.0)
    assert float(summary_rows[-1]["time_ns"]) == pytest.approx(20.0)

    with open(outdir / "residue_occupancy.csv", encoding="utf-8") as fh:
        occ_rows = list(csv.DictReader(fh))
    assert occ_rows, "residue_occupancy.csv is empty"
    labels = {r["label"] for r in occ_rows}
    assert any("MET769" in lbl for lbl in labels)

    with open(outdir / "trajectory_fingerprints.csv", encoding="utf-8") as fh:
        fp_rows = list(csv.DictReader(fh))
    assert len(fp_rows) == 201
    assert len(fp_rows[0]) > 4

    assert (outdir / "trajectory_timeline.png").stat().st_size > 5000

    # Verify per-interval snapshot directories: 0 ns, 10 ns, 20 ns
    snap_dirs = sorted(p.name for p in (outdir / "snapshots").iterdir() if p.is_dir())
    assert snap_dirs == ["t_00000ns", "t_00010ns", "t_00020ns"]
    for name in snap_dirs:
        d = outdir / "snapshots" / name
        assert (d / f"{name}.pdb").stat().st_size > 1000
        assert (d / f"{name}_interactions.png").stat().st_size > 1000
        assert (d / f"{name}_interactions_report.txt").stat().st_size > 200
        data = json.loads((d / f"{name}_interactions.json").read_text(encoding="utf-8"))
        assert data["ligand"]["name"] == "AQ4"


@needs_1m17
def test_xtc_every_1ns_with_time_window(tmp_path):
    """--every 1ns with --start 2ns --end 5ns selects exactly 2, 3, 4, 5 ns."""
    exe = find_rust_engine(auto_build=False)
    xtc_path = tmp_path / "prod_10ns.xtc"
    subprocess.run(
        [
            exe, "make-xtc",
            "--topology", EXAMPLE_1M17,
            "--output", str(xtc_path),
            "--frames", "101",
            "--dt-ps", "100.0",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    res = analyze_md_trajectory(
        topology_file=EXAMPLE_1M17,
        traj_file=str(xtc_path),
        output_dir=str(tmp_path / "win"),
        ligand_resname="AQ4",
        every="1ns",
        start="2ns",
        end="5ns",
        diagrams=(),
        json_export=True,
        write_timeline=False,
    )
    assert res["n_frames"] == 101
    labels = [s["label"] for s in res["snapshots"]]
    assert labels == ["t_00002ns", "t_00003ns", "t_00004ns", "t_00005ns"]


def test_rust_engine_multi_pose_pdbqt_and_gro(tmp_path):
    """Rust engine handles both multi-pose PDBQT docking outputs and GRO topologies."""
    pdbqt = """MODEL 1
ATOM      1  OG1 THR A 113      12.000  12.000  10.000  1.00 20.00    -0.400 OA
HETATM    2  BR1 LIG B   1      14.900  12.400  10.000  1.00 20.00     0.050 Br
HETATM    3  C1  LIG B   1      14.000  13.500  10.000  1.00 20.00     0.000 A
ENDMDL
MODEL 2
ATOM      1  OG1 THR A 113      12.000  12.000  10.000  1.00 20.00    -0.400 OA
HETATM    2  BR1 LIG B   1      15.100  12.400  10.000  1.00 20.00     0.050 Br
HETATM    3  C1  LIG B   1      14.200  13.500  10.000  1.00 20.00     0.000 A
ENDMDL
"""
    pq_path = tmp_path / "docked.pdbqt"
    pq_path.write_text(pdbqt)
    res = analyze_md_trajectory(
        topology_file=str(pq_path),
        output_dir=str(tmp_path / "pq_out"),
        ligand_resname="LIG",
        stride=1,
        diagrams=(),
        write_timeline=False,
    )
    assert res["n_frames"] == 2
    with open(res["summary_csv"], encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2
    assert int(rows[0]["halogen_bonds"]) == 1
    assert int(rows[1]["halogen_bonds"]) == 1

    # GRO topology (coordinates in nm: THR OG1 at 1.200,1.200,1.000 nm; LIG O1 at 1.490,1.200,1.000 nm -> 2.9 A H-bond)
    gro = (
        "PandaMap GRO test\n"
        "    3\n"
        "  113THR    OG1    1   1.200   1.200   1.000\n"
        "    1LIG     O1    2   1.490   1.200   1.000\n"
        "    1LIG     C1    3   1.590   1.300   1.000\n"
        "   5.00000   5.00000   5.00000\n"
    )
    gro_path = tmp_path / "sys.gro"
    gro_path.write_text(gro)
    res_gro = analyze_md_trajectory(
        topology_file=str(gro_path),
        output_dir=str(tmp_path / "gro_out"),
        ligand_resname="LIG",
        diagrams=(),
        write_timeline=False,
    )
    assert res_gro["n_frames"] == 1
    with open(res_gro["summary_csv"], encoding="utf-8") as fh:
        row = next(csv.DictReader(fh))
    assert int(row["hydrogen_bonds"]) == 1

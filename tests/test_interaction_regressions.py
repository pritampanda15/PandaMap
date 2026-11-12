"""
Regression tests for interaction-detection bugs fixed in this release.

Each test pins one previously-broken behaviour.  The fixtures are minimal
synthetic PDB files written at run time, so the suite needs no external
structure downloads and runs in well under a second.
"""
import os
import textwrap

import pytest

from pandamap.core import (
    HybridProtLigMapper,
    HALOGEN_ELEMENTS,
    STRICT_HALOGEN_ELEMENTS,
    COVALENT_LIGAND_ELEMENTS,
    HBOND_DON_ANGLE_MIN,
)


def _write(tmp_path, name, body):
    path = tmp_path / name
    path.write_text(textwrap.dedent(body).lstrip("\n"))
    return str(path)


# --------------------------------------------------------------------------
# Halogen element casing
# --------------------------------------------------------------------------

BROMINE_COMPLEX = """
ATOM      1  N   THR A 113      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  THR A 113      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  CB  THR A 113      11.500  11.000  10.000  1.00 20.00           C
ATOM      4  OG1 THR A 113      12.000  12.000  10.000  1.00 20.00           O
HETATM    5  C1  LIG B 001      14.000  13.500  10.000  1.00 20.00           C
HETATM    6  BR1 LIG B 001      14.900  12.400  10.000  1.00 20.00          BR
END
"""


def test_bromine_is_detected_as_halogen(tmp_path):
    """BioPython reports element 'BR'; a case-sensitive membership test
    silently dropped every Br/Cl/I, so only fluorine was ever detected."""
    mapper = HybridProtLigMapper(
        _write(tmp_path, "br.pdb", BROMINE_COMPLEX), ligand_resname="LIG"
    )
    elements = {a.element.strip().capitalize() for a in mapper.halogen_atoms}
    assert "Br" in elements, "bromine was not collected as a halogen donor"


def test_halogen_element_sets_are_consistent():
    assert STRICT_HALOGEN_ELEMENTS < HALOGEN_ELEMENTS
    assert "F" in HALOGEN_ELEMENTS and "F" not in STRICT_HALOGEN_ELEMENTS


def test_strict_halogens_excludes_fluorine(tmp_path):
    mapper = HybridProtLigMapper(
        _write(tmp_path, "br2.pdb", BROMINE_COMPLEX),
        ligand_resname="LIG",
        strict_halogens=True,
    )
    mapper.detect_interactions()
    for contact in mapper.interactions["halogen_bonds"]:
        assert contact["halogen_element"] != "F"


def test_halogen_contacts_record_their_element(tmp_path):
    mapper = HybridProtLigMapper(
        _write(tmp_path, "br3.pdb", BROMINE_COMPLEX), ligand_resname="LIG"
    )
    mapper.detect_interactions()
    for contact in mapper.interactions["halogen_bonds"]:
        assert contact.get("halogen_element") in HALOGEN_ELEMENTS


# --------------------------------------------------------------------------
# Double counting
# --------------------------------------------------------------------------

# Phosphonate rather than carboxylate: this group is perceived as anionic by
# both the old and the new charge rules, so the double-counting tests below
# exercise the real code path instead of passing vacuously because no ionic
# contact was detected in the first place.
SALT_BRIDGE_COMPLEX = """
ATOM      1  N   LYS A  10      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  LYS A  10      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  CB  LYS A  10      11.500  11.000  10.000  1.00 20.00           C
ATOM      4  NZ  LYS A  10      12.000  12.000  10.000  1.00 20.00           N
HETATM    5  P1  LIG B 001      15.200  13.500  10.000  1.00 20.00           P
HETATM    6  O1  LIG B 001      14.000  12.900  10.000  1.00 20.00           O
HETATM    7  O2  LIG B 001      15.900  14.700  10.000  1.00 20.00           O
HETATM    8  O3  LIG B 001      16.100  12.500  10.000  1.00 20.00           O
END
"""

# Carboxylate fixture, used for the formal-charge perception tests.
CARBOXYLATE_COMPLEX = """
ATOM      1  N   LYS A  10      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  LYS A  10      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  CB  LYS A  10      11.500  11.000  10.000  1.00 20.00           C
ATOM      4  NZ  LYS A  10      12.000  12.000  10.000  1.00 20.00           N
HETATM    5  C1  LIG B 001      14.500  13.500  10.000  1.00 20.00           C
HETATM    6  O1  LIG B 001      15.000  12.600  10.000  1.00 20.00           O
HETATM    7  O2  LIG B 001      15.100  14.500  10.000  1.00 20.00           O
END
"""


def test_ionic_contacts_are_not_also_salt_bridges(tmp_path):
    """Every ionic contact used to be appended to both keys and then scored
    twice (-1.80 and -1.50)."""
    mapper = HybridProtLigMapper(
        _write(tmp_path, "salt.pdb", SALT_BRIDGE_COMPLEX), ligand_resname="LIG"
    )
    mapper.detect_interactions()
    ionic = mapper.interactions["ionic"]
    bridges = mapper.interactions["salt_bridge"]
    assert not (ionic and bridges), (
        "ionic and salt_bridge are populated simultaneously; "
        "charge-charge contacts would be double-weighted"
    )


def test_salt_bridge_has_no_independent_weight(tmp_path):
    """'salt_bridge' must not carry its own dG weight while it aliases 'ionic'."""
    mapper = HybridProtLigMapper(
        _write(tmp_path, "salt2.pdb", SALT_BRIDGE_COMPLEX), ligand_resname="LIG"
    )
    mapper.detect_interactions()
    breakdown = mapper.estimate_binding_affinity()["breakdown"]
    assert not ("ionic" in breakdown and "salt_bridge" in breakdown)


def test_cation_pi_and_pi_cation_are_disjoint(tmp_path):
    """The two pi-system classes matched overlapping geometry, so a single
    protein-cation/ligand-ring contact scored -0.65 twice."""
    body = """
    ATOM      1  N   ARG A  20      10.000  10.000  10.000  1.00 20.00           N
    ATOM      2  CA  ARG A  20      11.000  10.000  10.000  1.00 20.00           C
    ATOM      3  NH1 ARG A  20      12.000  11.000  10.000  1.00 20.00           N
    HETATM    4  C1  LIG B 001      14.000  12.000  10.000  1.00 20.00           C
    HETATM    5  C2  LIG B 001      15.200  12.400  10.000  1.00 20.00           C
    HETATM    6  C3  LIG B 001      16.000  11.400  10.000  1.00 20.00           C
    HETATM    7  C4  LIG B 001      15.600  10.100  10.000  1.00 20.00           C
    HETATM    8  C5  LIG B 001      14.400   9.700  10.000  1.00 20.00           C
    HETATM    9  C6  LIG B 001      13.600  10.700  10.000  1.00 20.00           C
    END
    """
    mapper = HybridProtLigMapper(
        _write(tmp_path, "catpi.pdb", body), ligand_resname="LIG"
    )
    mapper.detect_interactions()

    def keys(contacts):
        return {
            (c["protein_atom"].get_serial_number(), c["ligand_atom"].get_serial_number())
            for c in contacts
        }

    overlap = keys(mapper.interactions["cation_pi"]) & keys(
        mapper.interactions["pi_cation"]
    )
    assert not overlap, f"{len(overlap)} contact(s) counted under both pi classes"


# --------------------------------------------------------------------------
# Ligand formal charge perception
# --------------------------------------------------------------------------


def test_carboxylate_oxygen_is_recognised_as_anionic(tmp_path):
    """The old rule ("oxygen with no carbon neighbour") excluded every
    carboxylate -- the commonest anionic group in drug-like ligands."""
    mapper = HybridProtLigMapper(
        _write(tmp_path, "carb.pdb", CARBOXYLATE_COMPLEX), ligand_resname="LIG"
    )
    carboxylate_os = [a for a in mapper.ligand_atoms if a.get_name() in ("O1", "O2")]
    assert carboxylate_os
    assert any(mapper._is_ligand_anion(a) for a in carboxylate_os)


def test_hydroxyl_oxygen_is_not_anionic(tmp_path):
    body = """
    HETATM    1  C1  LIG B 001      14.000  13.500  10.000  1.00 20.00           C
    HETATM    2  O1  LIG B 001      15.000  12.600  10.000  1.00 20.00           O
    END
    """
    mapper = HybridProtLigMapper(
        _write(tmp_path, "oh.pdb", body), ligand_resname="LIG"
    )
    hydroxyl = [a for a in mapper.ligand_atoms if a.get_name() == "O1"][0]
    assert not mapper._is_ligand_anion(hydroxyl)


# --------------------------------------------------------------------------
# Hydrogen-bond geometry
# --------------------------------------------------------------------------

# Donor N-H points directly at the acceptor: D-H...A near 180 deg.
GOOD_HBOND = """
ATOM      1  N   SER A  39      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  H   SER A  39      10.980  10.000  10.000  1.00 20.00           H
ATOM      3  CA  SER A  39       9.500   9.000  10.000  1.00 20.00           C
HETATM    4  O1  LIG B 001      12.900  10.000  10.000  1.00 20.00           O
HETATM    5  C1  LIG B 001      13.900  10.000  10.000  1.00 20.00           C
END
"""

# Same D...A distance, but the hydrogen points away: D-H...A well below 100 deg.
BAD_HBOND = """
ATOM      1  N   SER A  39      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  H   SER A  39       9.020  10.000  10.000  1.00 20.00           H
ATOM      3  CA  SER A  39       9.500   9.000  10.000  1.00 20.00           C
HETATM    4  O1  LIG B 001      12.900  10.000  10.000  1.00 20.00           O
HETATM    5  C1  LIG B 001      13.900  10.000  10.000  1.00 20.00           C
END
"""


def test_linear_hbond_is_kept(tmp_path):
    mapper = HybridProtLigMapper(
        _write(tmp_path, "hb_good.pdb", GOOD_HBOND), ligand_resname="LIG"
    )
    mapper.detect_interactions()
    bonds = mapper.interactions["hydrogen_bonds"]
    assert len(bonds) == 1
    assert bonds[0]["angle_source"] == "explicit_H"
    assert bonds[0]["angle"] >= HBOND_DON_ANGLE_MIN


def test_bent_hbond_is_rejected(tmp_path):
    """Distance alone would accept this contact; the D-H...A angle must not."""
    mapper = HybridProtLigMapper(
        _write(tmp_path, "hb_bad.pdb", BAD_HBOND), ligand_resname="LIG"
    )
    mapper.detect_interactions()
    assert mapper.interactions["hydrogen_bonds"] == []


def test_distance_only_when_no_hydrogens(tmp_path):
    """Structures without explicit H must keep the published distance-only
    behaviour by default -- the heavy-atom proxy is opt-in."""
    body = GOOD_HBOND.replace(
        "ATOM      2  H   SER A  39      10.980  10.000  10.000  1.00 20.00           H\n",
        "",
    )
    mapper = HybridProtLigMapper(
        _write(tmp_path, "hb_noh.pdb", body), ligand_resname="LIG"
    )
    mapper.detect_interactions()
    bonds = mapper.interactions["hydrogen_bonds"]
    assert len(bonds) == 1
    assert bonds[0]["angle"] is None
    assert bonds[0]["angle_source"] == "heavy_atom"


# --------------------------------------------------------------------------
# Covalent bonds must not swallow metal coordination
# --------------------------------------------------------------------------


def test_metal_contact_is_not_scored_as_covalent(tmp_path):
    """A Zn-His contact at ~2.0 A is coordination; filing it as covalent cost
    -8.0 kcal/mol in the scoring function."""
    body = """
    ATOM      1  N   HIS A  94      10.000  10.000  10.000  1.00 20.00           N
    ATOM      2  CA  HIS A  94      11.000  10.000  10.000  1.00 20.00           C
    ATOM      3  CG  HIS A  94      11.500  11.000  10.000  1.00 20.00           C
    ATOM      4  ND1 HIS A  94      12.000  12.000  10.000  1.00 20.00           N
    HETATM    5 ZN    ZN A 200      13.400  13.400  10.000  1.00 20.00          ZN
    END
    """
    mapper = HybridProtLigMapper(
        _write(tmp_path, "zn.pdb", body), ligand_resname="ZN"
    )
    mapper.detect_interactions()
    assert mapper.interactions["covalent"] == []


def test_covalent_ligand_elements_exclude_metals():
    for metal in ("Zn", "Fe", "Mg", "Mn", "Ca", "Cu"):
        assert metal not in COVALENT_LIGAND_ELEMENTS


# --------------------------------------------------------------------------
# Outputs
# --------------------------------------------------------------------------


def test_csv_export_round_trips(tmp_path):
    import csv

    mapper = HybridProtLigMapper(
        _write(tmp_path, "csv.pdb", SALT_BRIDGE_COMPLEX), ligand_resname="LIG"
    )
    mapper.detect_interactions()
    out = mapper.export_interactions_csv(str(tmp_path / "out.csv"))
    with open(out) as fh:
        rows = list(csv.DictReader(fh))
    assert rows, "no interactions exported"
    expected = {
        "interaction_type", "protein_residue", "protein_resnum", "protein_chain",
        "protein_atom", "ligand_atom", "ligand_element", "distance_A",
        "solvent_accessible", "halogen_element",
    }
    assert expected <= set(rows[0])
    assert all(float(r["distance_A"]) > 0 for r in rows)


def test_python_api_report_names_residues(tmp_path):
    """Passing mapper.interactions straight to generate_report() -- exactly
    what the README's Python API section shows -- used to render every line
    as 'UNK?? -- 2.79A -- LIG', because the formatter reads flat
    restype/resnr/reschain keys that only the CLI path populated."""
    from pandamap.improved_interaction_detection import ImprovedInteractionDetection

    mapper = HybridProtLigMapper(
        _write(tmp_path, "api.pdb", GOOD_HBOND), ligand_resname="LIG"
    )
    mapper.detect_interactions()

    out = tmp_path / "report.txt"
    ImprovedInteractionDetection().generate_report(
        ligand_info={
            "hetid": mapper.ligand_residue.resname,
            "chain": mapper.ligand_residue.parent.id,
            "position": mapper.ligand_residue.id[1],
            "longname": mapper.ligand_residue.resname,
            "type": "LIGAND",
        },
        interactions=mapper.interactions,
        output_file=str(out),
    )
    text = out.read_text()
    assert "UNK??" not in text, "residue identity lost in the Python API report"
    assert "SER39" in text.replace(" ", "")


def test_normalize_interactions_is_idempotent(tmp_path):
    """Records that already carry flat keys must survive untouched, so the
    helper is safe on both the CLI and the Python API path."""
    from pandamap.improved_interaction_detection import normalize_interactions

    already_flat = {
        "hydrogen_bonds": [
            {"restype": "LEU", "resnr": 83, "reschain": "A", "distance": 2.79}
        ]
    }
    once = normalize_interactions(already_flat)
    twice = normalize_interactions(once)
    assert once == twice == already_flat


def test_solvent_accessibility_tiers_do_not_recurse(tmp_path):
    """estimate_solvent_accessibility() delegates to the Shrake-Rupley tier,
    whose own fallback must not call back into it."""
    mapper = HybridProtLigMapper(
        _write(tmp_path, "sasa.pdb", SALT_BRIDGE_COMPLEX), ligand_resname="LIG"
    )
    mapper.detect_interactions()
    assert hasattr(mapper, "estimate_solvent_accessibility")
    mapper.estimate_solvent_accessibility()      # must terminate


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))


# --------------------------------------------------------------------------
# Ligand selection, input handling and encoding (user-reported problems)
# --------------------------------------------------------------------------

MIXED_HETATM_COMPLEX = """
ATOM      1  N   THR A 113      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  THR A 113      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  CB  THR A 113      11.500  11.000  10.000  1.00 20.00           C
ATOM      4  OG1 THR A 113      12.000  12.000  10.000  1.00 20.00           O
HETATM    5  C1  LIG B   1      14.000  13.500  10.000  1.00 20.00           C
HETATM    6  O1  LIG B   1      14.900  12.400  10.000  1.00 20.00           O
HETATM    7  C2  LIG B   1      15.500  14.000  10.000  1.00 20.00           C
HETATM    8  C1  LIG C   2      34.000  13.500  10.000  1.00 20.00           C
HETATM    9  O1  LIG C   2      34.900  12.400  10.000  1.00 20.00           O
HETATM   10  C2  LIG C   2      35.500  14.000  10.000  1.00 20.00           C
HETATM   11  S   SO4 D   3      20.000  20.000  20.000  1.00 20.00           S
HETATM   12  O1  SO4 D   3      21.000  20.000  20.000  1.00 20.00           O
HETATM   13  O2  SO4 D   3      19.000  20.000  20.000  1.00 20.00           O
HETATM   14  O3  SO4 D   3      20.000  21.000  20.000  1.00 20.00           O
HETATM   15  O4  SO4 D   3      20.000  19.000  20.000  1.00 20.00           O
HETATM   16 ZN    ZN E   4      25.000  25.000  25.000  1.00 20.00          ZN
HETATM   17  O   HOH W   5      13.000  12.500  11.000  1.00 20.00           O
HETATM   18  O   HOH W   6      12.500  13.000   9.000  1.00 20.00           O
END
"""


def _mixed(tmp_path):
    return _write(tmp_path, "mixed.pdb", MIXED_HETATM_COMPLEX)


def test_auto_detect_ignores_water_ions_and_additives(tmp_path):
    """With no ligand name every HETATM group (waters included) used to be
    merged into one 'ligand', inflating hydrogen-bond counts and the score."""
    mapper = HybridProtLigMapper(_mixed(tmp_path))
    assert mapper.ligand_residue.resname == "LIG"
    assert {a.get_parent().resname for a in mapper.ligand_atoms} == {"LIG"}
    assert len(mapper.ligand_atoms) == 3          # one copy, not both
    assert [a.element for a in mapper.metal_atoms] == ["ZN"]


def test_ligand_name_is_case_insensitive_and_copy_selectable(tmp_path):
    path = _mixed(tmp_path)
    assert HybridProtLigMapper(path, ligand_resname="lig").ligand_residue.get_parent().id == "B"
    second = HybridProtLigMapper(path, ligand_resname="LIG", ligand_chain="C")
    assert second.ligand_residue.id[1] == 2
    assert HybridProtLigMapper(path, ligand_resname="SO4").ligand_residue.resname == "SO4"


def test_unknown_ligand_error_lists_available_ligands(tmp_path):
    with pytest.raises(ValueError, match=r"'XXX' not found.*LIG \(chain B"):
        HybridProtLigMapper(_mixed(tmp_path), ligand_resname="XXX")


def test_list_ligands(tmp_path):
    from pandamap import list_ligands
    names = [lig["resname"] for lig in list_ligands(_mixed(tmp_path))]
    assert names[0] == "SO4" and "LIG" in names and "HOH" not in names


def test_unparseable_file_gives_clear_error(tmp_path):
    with pytest.raises(ValueError, match="No atoms could be read"):
        HybridProtLigMapper(_write(tmp_path, "junk.pdb", "not a structure\n"))


def test_unknown_extension_is_parsed_by_content(tmp_path):
    mapper = HybridProtLigMapper(_write(tmp_path, "c.ent", BROMINE_COMPLEX))
    assert mapper.ligand_residue.resname == "LIG"


def test_report_is_written_as_utf8_regardless_of_locale(tmp_path, monkeypatch):
    """On Windows the default cp1252 codec cannot encode the report's
    non-ASCII characters ('π', 'Å'), so it was silently not written
    (GitHub issue #5)."""
    import builtins
    real_open = builtins.open

    def ascii_default_open(file, mode="r", *args, **kwargs):
        if "b" not in mode and "encoding" not in kwargs and len(args) < 2:
            kwargs["encoding"] = "ascii"
        return real_open(file, mode, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", ascii_default_open)
    mapper = HybridProtLigMapper(_mixed(tmp_path), ligand_resname="LIG")
    mapper.detect_interactions()
    report = tmp_path / "report.txt"
    mapper.generate_interaction_report(output_file=str(report))
    assert "Å" in report.read_bytes().decode("utf-8")


def test_cli_version_needs_no_structure_file(capsys):
    from pandamap.cli import main
    import sys
    argv, sys.argv = sys.argv, ["pandamap", "--version"]
    try:
        with pytest.raises(SystemExit) as exc:
            main()
    finally:
        sys.argv = argv
    assert exc.value.code == 0
    assert "PandaMap version" in capsys.readouterr().out


# --------------------------------------------------------------------------
# Modern 2D diagram
# --------------------------------------------------------------------------

def test_spread_angles_keeps_order_and_minimum_gap():
    import math
    from pandamap.diagram import _spread_angles
    preferred = [0.0, 0.01, 0.02, 0.03, 2.0, -3.0]
    gap = math.radians(18)
    order, angles = _spread_angles(preferred, gap)
    assert sorted(order) == list(range(len(preferred)))
    diffs = [angles[i + 1] - angles[i] for i in range(len(angles) - 1)]
    diffs.append(angles[0] + 2 * math.pi - angles[-1])
    assert min(diffs) >= gap - 1e-6
    # the isolated residues stay where they bind
    assert abs(angles[order.index(4)] - 2.0) < 1e-6


def test_smallest_ring_finds_ring_and_ignores_chains():
    from pandamap.diagram import _smallest_ring
    # six-membered ring 0-5 with a two-atom tail 5-6-7
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (5, 6), (6, 7)]
    neighbours = {i: [] for i in range(8)}
    for a, b in edges:
        neighbours[a].append(b)
        neighbours[b].append(a)
    assert sorted(_smallest_ring(2, neighbours)) == [0, 1, 2, 3, 4, 5]
    assert _smallest_ring(7, neighbours) == []


@pytest.mark.parametrize("name", ["diagram.png", "diagram.svg"])
def test_modern_diagram_is_written(tmp_path, name):
    mapper = HybridProtLigMapper(_mixed(tmp_path), ligand_resname="LIG")
    mapper.detect_interactions()
    out = tmp_path / name
    assert mapper.visualize(output_file=str(out), dpi=60,
                            show_nearby_residues=True) == str(out)
    assert out.stat().st_size > 1000


def test_classic_style_still_available_and_unknown_style_rejected(tmp_path):
    mapper = HybridProtLigMapper(_mixed(tmp_path), ligand_resname="LIG")
    mapper.detect_interactions()
    out = tmp_path / "classic.png"
    mapper.visualize(output_file=str(out), dpi=60, style="classic")
    assert out.stat().st_size > 1000
    with pytest.raises(ValueError, match="Unknown diagram style"):
        mapper.visualize(output_file=str(out), style="fancy")


# --------------------------------------------------------------------------
# 3D contact render, distant metals, deterministic solvent accessibility
# --------------------------------------------------------------------------

DISTANT_ZINC_COMPLEX = MIXED_HETATM_COMPLEX.replace(
    "HETATM   16 ZN    ZN E   4      25.000  25.000  25.000",
    "ATOM     20  N   CYS A 200      40.000  40.000  41.500  1.00 20.00           N\n"
    "ATOM     21  CA  CYS A 200      40.000  41.000  40.500  1.00 20.00           C\n"
    "ATOM     22  CB  CYS A 200      40.000  40.500  39.000  1.00 20.00           C\n"
    "ATOM     23  SG  CYS A 200      40.000  38.700  38.800  1.00 20.00           S\n"
    "HETATM   16 ZN    ZN E   4      40.000  38.700  36.500")


def test_metal_far_from_ligand_is_not_a_ligand_interaction(tmp_path):
    """A structural zinc 30 A from the ligand was reported (and drawn) as a
    metal-coordination contact of the ligand."""
    mapper = HybridProtLigMapper(_write(tmp_path, "zn.pdb", DISTANT_ZINC_COMPLEX),
                                 ligand_resname="LIG")
    mapper.detect_interactions()
    assert [a.element for a in mapper.metal_atoms] == ["ZN"]
    assert mapper.interactions["metal_coordination"] == []
    assert ("CYS", 200) not in mapper.interacting_residues


@pytest.mark.parametrize("name", ["contacts.png", "contacts.svg"])
def test_contacts_3d_render_is_written(tmp_path, name):
    mapper = HybridProtLigMapper(_mixed(tmp_path), ligand_resname="LIG")
    mapper.detect_interactions()
    assert sum(len(v) for v in mapper.interactions.values()) > 0
    out = tmp_path / name
    assert mapper.render_contacts_3d(str(out), dpi=60) == str(out)
    assert out.stat().st_size > 1000


def test_solvent_accessibility_is_the_same_on_every_run():
    """Score ties were broken by set iteration order (string hashing is
    randomised per process), so the residues marked solvent accessible
    changed from one run to the next on the same file."""
    import os
    import subprocess
    import sys
    example = os.path.join(os.path.dirname(__file__), "..", "example", "1m17.pdb")
    if not os.path.exists(example):
        pytest.skip("example structure not available")
    code = (
        "import io, contextlib\n"
        "from pandamap.core import HybridProtLigMapper\n"
        "with contextlib.redirect_stdout(io.StringIO()):\n"
        f"    m = HybridProtLigMapper({example!r}, ligand_resname='AQ4')\n"
        "    m.detect_interactions()\n"
        "    m.calculate_realistic_solvent_accessibility()\n"
        "print(sorted(m.solvent_accessible))\n"
    )
    results = set()
    for seed in ("1", "2", "3"):
        env = dict(os.environ, PYTHONHASHSEED=seed, PANDAMAP_NO_UPDATE_CHECK="1")
        out = subprocess.run([sys.executable, "-c", code], env=env, check=True,
                             capture_output=True, text=True).stdout
        results.add(out.strip().splitlines()[-1])
    assert len(results) == 1


def test_run_analysis_writes_all_three_diagrams_by_default(tmp_path):
    mapper = HybridProtLigMapper(_mixed(tmp_path), ligand_resname="LIG")
    out = tmp_path / "site.png"
    primary = mapper.run_analysis(output_file=str(out), dpi=50, use_dssp=False)
    assert primary == str(out)
    assert mapper.diagram_files == {
        "modern": str(out),
        "classic": str(tmp_path / "site_classic.png"),
        "3d": str(tmp_path / "site_3d.png"),
    }
    for path in mapper.diagram_files.values():
        assert os.path.getsize(path) > 1000


def test_run_analysis_diagram_selection_never_overwrites(tmp_path):
    mapper = HybridProtLigMapper(_mixed(tmp_path), ligand_resname="LIG")
    out = tmp_path / "site.png"
    mapper.run_analysis(output_file=str(out), dpi=50, use_dssp=False,
                        diagrams=["3d", "modern"])
    assert mapper.diagram_files == {"3d": str(out),
                                    "modern": str(tmp_path / "site_modern.png")}
    assert not (tmp_path / "site_classic.png").exists()


def test_cli_diagrams_option(tmp_path, monkeypatch, capsys):
    import sys
    from pandamap.cli import main
    out = tmp_path / "only.png"
    monkeypatch.setattr(sys, "argv", ["pandamap", _mixed(tmp_path), "--ligand", "LIG",
                                      "--diagrams", "classic", "--dpi", "50",
                                      "-o", str(out)])
    assert main() == 0
    assert out.exists()
    assert sorted(p.name for p in tmp_path.glob("only*")) == ["only.png"]
    assert "classic" in capsys.readouterr().out


def test_visualize_3d_style(tmp_path):
    mapper = HybridProtLigMapper(_mixed(tmp_path), ligand_resname="LIG")
    mapper.detect_interactions()
    out = tmp_path / "v3d.png"
    assert mapper.visualize(output_file=str(out), dpi=50, style="3d") == str(out)
    assert out.stat().st_size > 1000

"""
Regression tests for interaction-detection bugs fixed in this release.

Each test pins one previously-broken behaviour.  The fixtures are minimal
synthetic PDB files written at run time, so the suite needs no external
structure downloads and runs in well under a second.
"""
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

#!/usr/bin/env python
"""
Core functionality for PandaMap: A Python package for visualizing 
protein-ligand interactions with 2D ligand structure representation.
"""

import logging
import os
import sys
import math
from collections import defaultdict
import tempfile
import subprocess
from Bio.PDB import PDBParser, MMCIFParser
import numpy as np
import subprocess
from Bio.PDB.DSSP import DSSP
from Bio.PDB import PDBIO
# Set matplotlib backend to 'Agg' BEFORE importing pyplot
import matplotlib
matplotlib.use('Agg')  # This must come before any pyplot import
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.colors import to_rgb
from matplotlib.patches import Circle, FancyArrowPatch, Wedge, Rectangle, Polygon
from matplotlib.lines import Line2D
import matplotlib.patheffects as path_effects
from .improved_interaction_detection import add_interaction_detection_methods


# BioPython imports
from Bio.PDB import PDBParser, NeighborSearch

logger = logging.getLogger("pandamap")

# Exact ring-atom names for each aromatic amino acid.
# Only these atoms are considered part of the aromatic pi-system.
# Halogen-bond donors.  Note that organic fluorine has a negligible sigma-hole
# and is not a halogen-bond donor in the strict sense of Auffinger et al. (2004),
# which defines C-X...O for X = Cl, Br, I.  F is collected by default (PLIP also
# reports C-F...O contacts) but every detected halogen bond carries a
# 'halogen_element' field so F contacts can be separated, and
# strict_halogens=True drops them entirely.
HALOGEN_ELEMENTS = frozenset(['F', 'Cl', 'Br', 'I'])
STRICT_HALOGEN_ELEMENTS = frozenset(['Cl', 'Br', 'I'])

# Ligand elements that can form a covalent bond to a protein nucleophile.
# Deliberately excludes metals, whose short contacts are coordination bonds.
COVALENT_LIGAND_ELEMENTS = frozenset(['C', 'N', 'O', 'S', 'P', 'Se', 'B'])

# Minimum D-H...A angle for a hydrogen bond, in degrees.  Matches PLIP's
# HBOND_DON_ANGLE_MIN default.  Applied exactly when explicit hydrogens are
# present; otherwise a heavy-atom direction test stands in.
HBOND_DON_ANGLE_MIN = 100.0

# The 2D/3D diagrams PandaMap can draw; run_analysis() writes all of them
# unless told otherwise.
DIAGRAM_STYLES = ('modern', 'classic', '3d')
_DIAGRAM_SUFFIX = {'modern': '_modern', 'classic': '_classic', '3d': '_3d'}

# Water bridges: both hydrogen bonds within this distance window (Angstrom) and
# the ligand-water-protein angle within this range (degrees), following PLIP.
WATER_BRIDGE_MIN, WATER_BRIDGE_MAX = 2.5, 3.5
WATER_BRIDGE_ANGLE = (71.0, 140.0)

SCORE_INTERPRETATION = ('Empirical interaction score on a kcal/mol scale. It is not '
                        'calibrated against measured affinities and is not a '
                        'predicted binding free energy.')
SCORE_NOTE = ('Use it to compare poses or close analogues of the same ligand, not to '
              'judge absolute binding strength. For affinities use FEP, MM-GBSA or '
              'experiment.')

# A metal ion counts as bridging protein and ligand only within this distance
# (in Angstrom) of a ligand atom.
METAL_LIGAND_CUTOFF = 3.0

AROMATIC_RING_ATOMS = {
    'PHE': frozenset(['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']),
    'TYR': frozenset(['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']),
    'TRP': frozenset(['CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'NE1']),
    'HIS': frozenset(['CG', 'ND1', 'CD2', 'CE1', 'NE2']),
}

# Define three_to_one conversion manually if import isn't available
try:
    from Bio.PDB.Polypeptide import three_to_one 
except ImportError:
    # Define the conversion dictionary manually
    _aa_index = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    
    def three_to_one(residue):
        """Convert amino acid three letter code to one letter code."""
        if residue in _aa_index:
            return _aa_index[residue]
        else:
            return "X"  # Unknown amino acid


# Kept importable from pandamap.core for existing code
from .structure_io import (parse_pdbqt, MultiFormatParser, _looks_like_mmcif,   # noqa: E402,F401
                           WATER_RESNAMES, NON_LIGAND_RESNAMES, _describe_ligand,
                           _ligand_candidates, list_ligands, _select_ligand)
from .classic_ligand import SimpleLigandStructure                                # noqa: E402,F401


class HybridProtLigMapper:
    """
    Class for analyzing protein-ligand interactions and creating 
    visualizations with a simplified ligand structure.
    """
    
    def __init__(self, structure_file, ligand_resname=None, strict_halogens=False,
                 strict_hbond_geometry=False, ligand_chain=None, ligand_resnum=None,
                 water_angle_check=True):
        """
        Initialize with a structure file containing a protein-ligand complex.

        Parameters:
        -----------
        structure_file : str
            Path to the structure file (PDB, mmCIF/CIF, or PDBQT format)
        ligand_resname : str, optional
            Residue name of the ligand (case-insensitive).  If omitted, the
            largest HETATM group that is not water, an ion or a common
            crystallisation additive is selected automatically.
        ligand_chain, ligand_resnum : str / int, optional
            Chain ID and residue number, to pick one copy when the structure
            contains the ligand more than once.
        strict_halogens : bool, optional
            If True, restrict halogen bonds to Cl/Br/I (Auffinger et al. 2004).
            Default False also reports C-F...O contacts, matching PLIP.
        strict_hbond_geometry : bool, optional
            If True, apply the heavy-atom angular proxy to structures that lack
            explicit hydrogens.  The exact D-H...A criterion is always applied
            when hydrogens are present, independent of this flag.
        water_angle_check : bool, optional
            Default True requires the angle at each bridging water to lie in
            71-140 degrees (PLIP).  False accepts any water bridge whose
            hydrogen bonds are all 2.5-3.5 A long (distance only, as in
            WatFinder-style analyses of placed or 3D-RISM waters).
        """
        self.strict_halogens = strict_halogens
        self.water_angle_check = water_angle_check
        self.strict_hbond_geometry = strict_hbond_geometry
        self.structure_file = structure_file
        self.ligand_resname = ligand_resname
        
        # Parse the structure file using the multi-format parser
        parser = MultiFormatParser()
        try:
            self.structure = parser.parse_structure(structure_file)
            models = list(self.structure.get_models())
        except (OSError, ValueError):
            raise
        except Exception as e:
            raise ValueError(f"Could not parse {structure_file}: {e}") from e
        if not models:
            raise ValueError(f"No atoms could be read from {structure_file}. "
                             "Is it a valid PDB, mmCIF or PDBQT file?")
        self.model = models[0]
        
        # Separate ligand from protein
        self.protein_atoms = []
        self.ligand_atoms = []
        self.protein_residues = {}
        self.ligand_residue = None
        
        # Store atom information including metals
        self.metal_atoms = []
        self.halogen_atoms = []
        self.water_atoms = []
        
        # Exactly one HETATM residue is the ligand.  (Previously, with no
        # ligand name given, every HETATM group - waters and ions included -
        # was merged into one "ligand".)
        selected_ligand = _select_ligand(self.model, structure_file,
                                         resname=ligand_resname,
                                         chain=ligand_chain, resnum=ligand_resnum)

        for residue in self.model.get_residues():
            # Store ligand atoms (HETATM records)
            if residue.id[0] != ' ':  # Non-standard residue (HETATM)
                if residue is selected_ligand:
                    for atom in residue:
                        self.ligand_atoms.append(atom)
                        # Identify halogen atoms in ligand.
                        # BioPython returns PDB element symbols upper-cased
                        # ('BR', 'CL'), so normalise before comparing or Br/Cl/I
                        # are silently dropped and only F is ever detected.
                        if atom.element.strip().capitalize() in HALOGEN_ELEMENTS:
                            self.halogen_atoms.append(atom)
                    if self.ligand_residue is None:
                        self.ligand_residue = residue
                elif residue.resname.strip().upper() in WATER_RESNAMES:
                    # Waters can bridge ligand and protein (see water bridges)
                    self.water_atoms.extend(
                        a for a in residue if (a.element or '').strip().upper() == 'O')
                else:
                    # Check for metal ions (typically single-atom HETATM residues)
                    if len(list(residue.get_atoms())) == 1:
                        atom = next(residue.get_atoms())
                        if atom.element in ['MG', 'ZN', 'CA', 'FE', 'MN', 'CU', 'NA', 'K', 'LI', 'CO', 'NI']:
                            self.metal_atoms.append(atom)
                    
            else:  # Standard residues (protein)
                res_id = (residue.resname, residue.id[1])
                self.protein_residues[res_id] = residue
                for atom in residue:
                    self.protein_atoms.append(atom)
        
        # Check if we found a ligand
        if not self.ligand_atoms:
            raise ValueError(f"No ligand (HETATM) found in the file: {structure_file}")
        
        # Storage for the interaction data - UPDATED with new interaction types
        self.interactions = {
            'hydrogen_bonds': [],
            'carbon_pi': [],
            'pi_pi_stacking': [],
            'donor_pi': [],
            'amide_pi': [],
            'hydrophobic': [],
            'ionic': [],           # Ionic interactions
            'halogen_bonds': [],   # Halogen bonds
            'cation_pi': [],       # Cation-pi interactions
            'metal_coordination': [], # Metal coordination
            'salt_bridge': [],     # Salt bridges
            'covalent': [],        # Covalent bonds
            'alkyl_pi': [],        # NEW: Alkyl-Pi interactions
            'attractive_charge': [], # NEW: Attractive charge interactions
            'pi_cation': [],         # NEW: Pi-cation interactions
            'repulsion': [],        # NEW: Repulsion interactions
            'water_bridges': []     # ligand ... water ... protein hydrogen bonds
        }
        
        # Will store residues that interact with the ligand
        self.interacting_residues = set()
        
        # Store interaction directionality (protein->ligand, ligand->protein, or both)
        self.interaction_direction = {}
        
        # For solvent accessibility information (simplified)
        self.solvent_accessible = set()
        
        # Create the simple ligand structure
        self.ligand_structure = SimpleLigandStructure(self.ligand_atoms)

        # Pre-compute ligand ring atoms for accurate pi-system detection
        self.ligand_ring_atoms = self.ligand_structure.detect_rings()

    def _attached_hydrogens(self, heavy_atom, max_dist=1.3):
        """Explicit hydrogens covalently attached to a heavy atom, if any."""
        parent = heavy_atom.get_parent()
        if parent is None:
            return []
        origin = heavy_atom.get_coord()
        return [a for a in parent.get_atoms()
                if a.element.strip().upper() == 'H'
                and np.linalg.norm(a.get_coord() - origin) < max_dist]

    def _best_hbond_angle(self, prot_atom, lig_atom):
        """
        Measure the true D-H···A angle when explicit hydrogens are present.

        Either partner may be the donor, so both directions are tested and the
        most favourable (largest) angle is returned.  Returns
        ``(angle_degrees, source)`` where source is 'explicit_H', or
        ``(None, 'heavy_atom')`` when neither candidate donor carries a
        hydrogen, signalling the caller to use the heavy-atom approximation.
        """
        best = None
        for donor, acceptor in ((prot_atom, lig_atom), (lig_atom, prot_atom)):
            for h in self._attached_hydrogens(donor):
                hv = donor.get_coord() - h.get_coord()
                av = acceptor.get_coord() - h.get_coord()
                nh, na = np.linalg.norm(hv), np.linalg.norm(av)
                if nh == 0 or na == 0:
                    continue
                cos_a = np.clip(np.dot(hv, av) / (nh * na), -1.0, 1.0)
                angle = math.degrees(math.acos(cos_a))
                if best is None or angle > best:
                    best = angle
        if best is None:
            return None, 'heavy_atom'
        return best, 'explicit_H'

    def _check_hbond_geometry(self, donor_atom, acceptor_atom, donor_is_protein=True):
        """
        Validate H-bond geometry without explicit hydrogen positions.

        Approximates the D-H direction as the vector *opposite* to the
        centroid of the donor's covalent neighbours (heavy atoms within
        1.8 Å in the same residue).  The interaction is accepted when the
        angle between this estimated H-direction and the D···A vector is
        < 60°, which corresponds to a D-H···A angle > 120° — consistent
        with crystallographically observed hydrogen-bond geometries.

        Parameters
        ----------
        donor_atom : Bio.PDB.Atom
            The heavy-atom donor (N or O).
        acceptor_atom : Bio.PDB.Atom
            The heavy-atom acceptor (O or N).
        donor_is_protein : bool
            True when the donor belongs to the protein; False for ligand donors.

        Returns
        -------
        bool
            True  → geometry is consistent with a valid H-bond.
            False → geometry is unlikely to support an H-bond.
        """
        donor_coord = donor_atom.get_coord()
        acceptor_coord = acceptor_atom.get_coord()

        # Collect covalently bonded heavy-atom neighbours of the donor
        donor_res = donor_atom.get_parent()
        cov_neighbor_coords = []

        if donor_is_protein:
            for atom in donor_res.get_atoms():
                if atom.get_id() == donor_atom.get_id():
                    continue
                if atom.element == 'H':
                    continue
                d = np.linalg.norm(atom.get_coord() - donor_coord)
                if 0.9 < d < 1.8:
                    cov_neighbor_coords.append(atom.get_coord())
        else:
            # Ligand donor — search within ligand atoms
            for atom in self.ligand_atoms:
                if atom.get_id() == donor_atom.get_id():
                    continue
                if atom.element == 'H':
                    continue
                d = np.linalg.norm(atom.get_coord() - donor_coord)
                if 0.9 < d < 1.8:
                    cov_neighbor_coords.append(atom.get_coord())

        if not cov_neighbor_coords:
            return True  # No geometry constraint possible; accept by default

        avg_neighbor = np.mean(cov_neighbor_coords, axis=0)

        # Estimated H-direction: opposite to the mean covalent-neighbour direction
        donor_to_avg_neighbor = avg_neighbor - donor_coord
        estimated_h_dir = -donor_to_avg_neighbor

        # Direction from donor to acceptor
        donor_to_acceptor = acceptor_coord - donor_coord

        n1 = np.linalg.norm(estimated_h_dir)
        n2 = np.linalg.norm(donor_to_acceptor)
        if n1 < 1e-6 or n2 < 1e-6:
            return True

        cos_angle = np.dot(estimated_h_dir, donor_to_acceptor) / (n1 * n2)
        cos_angle = float(np.clip(cos_angle, -1.0, 1.0))
        angle = math.degrees(math.acos(cos_angle))

        # Require estimated D-H···A angle > 110° (relaxed from ideal 120° to
        # tolerate small errors from the heavy-atom approximation)
        return angle > 110.0




    def detect_interactions(self,
                  h_bond_cutoff=3.5,
                  pi_stack_cutoff=5.5,
                  hydrophobic_cutoff=4.0,
                  ionic_cutoff=5.5,
                  halogen_bond_cutoff=3.5,
                  metal_coord_cutoff=2.8,
                  covalent_cutoff=2.1):
        """
        Detect all interactions between protein and ligand with PLIP-like thresholds.
        """
        # Use neighbor search for efficiency
        ns = NeighborSearch(self.protein_atoms)
        max_cutoff = max(h_bond_cutoff, pi_stack_cutoff, hydrophobic_cutoff, 
                    ionic_cutoff, halogen_bond_cutoff, metal_coord_cutoff,
                    covalent_cutoff)
        
        # Define amino acid categories
        aromatic_residues = {'PHE', 'TYR', 'TRP', 'HIS'}
        h_bond_donors = {'ARG', 'LYS', 'HIS', 'ASN', 'GLN', 'SER', 'THR', 'TYR', 'TRP'}
        h_bond_acceptors = {'ASP', 'GLU', 'ASN', 'GLN', 'HIS', 'SER', 'THR', 'TYR'}
        neg_charged = {'ASP', 'GLU'}
        pos_charged = {'ARG', 'LYS', 'HIS'}
        amide_residues = {'ASN', 'GLN'}
        hydrophobic_residues = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'TYR'}
        alkyl_residues = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO'}
        
        # Storage for the interaction data - make sure this includes the new types!
        self.interactions = {
            'hydrogen_bonds': [],
            'carbon_pi': [],
            'pi_pi_stacking': [],
            'donor_pi': [],
            'amide_pi': [],
            'hydrophobic': [],
            'ionic': [],
            'halogen_bonds': [],
            'cation_pi': [],
            'metal_coordination': [],
            'salt_bridge': [],
            'covalent': [],
            'alkyl_pi': [],
            'attractive_charge': [],
            'pi_cation': [],
            'repulsion': [],
            'water_bridges': []
        }
        
        # Helper functions
        def is_halogen_acceptor(atom):
            return atom.element in ['O', 'N', 'S'] or (atom.element == 'C' and atom.name in ['CE1', 'CD2'])
        
        def is_metal_coordinator(atom):
            return atom.element in ['O', 'N', 'S'] or atom.name in ['SD', 'OD1', 'OD2', 'OE1', 'OE2', 'NE2', 'ND1']
        
        # Initialize temporary storage for all interactions
        all_interactions = {key: [] for key in self.interactions.keys()}
        
        # Check each ligand atom for interactions
        for lig_atom in self.ligand_atoms:
            # Find protein atoms within cutoff distance
            nearby_atoms = ns.search(lig_atom.get_coord(), max_cutoff)
            
            # Track ionizable ligand atoms for ionizable interactions
            is_lig_pos_charged = self._is_ligand_cation(lig_atom)
            is_lig_neg_charged = self._is_ligand_anion(lig_atom)
            
            # Ligand aromaticity: use pre-computed ring membership (exact cycle detection)
            is_lig_aromatic = (lig_atom.element == 'C' and
                               lig_atom.get_id() in self.ligand_ring_atoms)
            
            for prot_atom in nearby_atoms:
                prot_res = prot_atom.get_parent()
                distance = lig_atom - prot_atom
                
                # Skip if distance is too small (likely a clash or error)
                if distance < 1.5:
                    continue
                    
                # Store interacting residue for later visualization
                res_id = (prot_res.resname, prot_res.id[1])
                self.interacting_residues.add(res_id)
                
                # 1. Hydrogen bonds — D···A distance window (2.5–3.5 Å) plus an
                # angular criterion.  When the structure carries explicit
                # hydrogens the true D-H···A angle is measured and must exceed
                # HBOND_DON_ANGLE_MIN (100°, PLIP's default).  Otherwise the
                # angle is approximated from the donor's heavy-atom geometry.
                if distance <= h_bond_cutoff and distance >= 2.5:
                    if lig_atom.element in ['N', 'O'] and prot_atom.element in ['N', 'O']:
                        res_id = (prot_res.resname, prot_res.id[1])

                        hb_angle, angle_source = self._best_hbond_angle(
                            prot_atom, lig_atom)
                        if hb_angle is not None:
                            # Explicit hydrogens available: enforce the true
                            # D-H···A criterion, exactly as PLIP does.
                            if hb_angle < HBOND_DON_ANGLE_MIN:
                                continue
                        elif self.strict_hbond_geometry:
                            # No explicit H.  The heavy-atom direction test is
                            # only a proxy and cannot distinguish donor from
                            # acceptor, so it is opt-in: applied as a hard
                            # filter it rejects many genuine bonds.
                            if not (self._check_hbond_geometry(prot_atom, lig_atom, True) or
                                    self._check_hbond_geometry(lig_atom, prot_atom, False)):
                                continue

                        interaction_info = {
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance,
                            'angle': hb_angle,
                            'angle_source': angle_source
                        }
                        all_interactions['hydrogen_bonds'].append(interaction_info)

                        # Determine directionality
                        interaction_key = (res_id, 'hydrogen_bonds')
                        is_donor_prot = prot_res.resname in h_bond_donors and prot_atom.element == 'N'
                        is_acceptor_prot = prot_res.resname in h_bond_acceptors and prot_atom.element in ['O', 'N']
                        is_donor_lig = lig_atom.element == 'N'
                        is_acceptor_lig = lig_atom.element in ['O', 'N']

                        if (is_donor_prot and is_acceptor_lig) and (is_donor_lig and is_acceptor_prot):
                            self.interaction_direction[interaction_key] = 'bidirectional'
                        elif is_donor_prot and is_acceptor_lig:
                            self.interaction_direction[interaction_key] = 'protein_to_ligand'
                        elif is_donor_lig and is_acceptor_prot:
                            self.interaction_direction[interaction_key] = 'ligand_to_protein'
                        else:
                            self.interaction_direction[interaction_key] = 'bidirectional'
                
                # 2. Pi-stacking - only between aromatic residues and aromatic ligand parts
                
                if distance <= pi_stack_cutoff:
                    is_prot_pos_charged = prot_res.resname in pos_charged and prot_atom.element in ['N']
                    is_lig_aromatic_ring = is_lig_aromatic and lig_atom.element == 'C'
                    
                    if is_prot_pos_charged and is_lig_aromatic_ring:
                        interaction_info = {
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance,
                            'type': 'protein_cation'
                        }
                        all_interactions['pi_cation'].append(interaction_info)
    
                # Protein atom is a genuine aromatic ring atom (atom-name based)
                # Defined here so it is available for all subsequent checks.
                is_prot_ring_atom = (
                    prot_res.resname in AROMATIC_RING_ATOMS and
                    prot_atom.get_name().strip() in AROMATIC_RING_ATOMS[prot_res.resname]
                )

                # 2. Alkyl-Pi interactions - between alkyl groups and aromatic systems
                if distance <= pi_stack_cutoff:
                    # Check for alkyl groups in protein interacting with aromatic ligand
                    is_prot_alkyl = prot_res.resname in alkyl_residues and prot_atom.element == 'C'
                    
                    if is_prot_alkyl and is_lig_aromatic:
                        interaction_info = {
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        }
                        all_interactions['alkyl_pi'].append(interaction_info)
            
                    # Check for aromatic groups in protein interacting with alkyl in ligand
                    is_prot_aromatic = is_prot_ring_atom and prot_atom.element == 'C'
                    is_lig_alkyl = lig_atom.element == 'C' and not is_lig_aromatic

                    # NB: this test was previously dedented out of the enclosing
                    # distance check, so it read is_prot_aromatic/is_lig_alkyl
                    # left over from an earlier loop iteration (and would raise
                    # NameError on the first pass if pi_stack_cutoff were ever
                    # lowered below the neighbour-search radius).
                    if is_prot_aromatic and is_lig_alkyl:
                        interaction_info = {
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        }
                        all_interactions['alkyl_pi'].append(interaction_info)
    
                # 3. Oppositely charged groups are recorded once, as 'ionic'
                # (section 7).  An identical test used to file every such
                # contact a second time under 'attractive_charge', which
                # inflated the contact count and the score.  The key is kept
                # (always empty) so existing code that reads it keeps working.

                # 4. Repulsion interactions - between likely same-charged groups
                # Fixed 4.0 Å cutoff: same-sign charges within this distance
                # produce genuine electrostatic repulsion (independent of ionic_cutoff).
                if distance <= 4.0:
                    # Define charges for protein residue explicitly
                    is_prot_pos_charged = prot_res.resname in pos_charged and prot_atom.element in ['N']
                    is_prot_neg_charged = prot_res.resname in neg_charged and prot_atom.element in ['O']
                    
                    if (is_prot_pos_charged and is_lig_pos_charged) or (is_prot_neg_charged and is_lig_neg_charged):
                        interaction_info = {
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        }
                        all_interactions['repulsion'].append(interaction_info)
                
                if distance <= pi_stack_cutoff:
                    if is_prot_ring_atom and is_lig_aromatic and \
                    lig_atom.element == 'C' and prot_atom.element == 'C':
                        interaction_info = {
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        }
                        all_interactions['pi_pi_stacking'].append(interaction_info)

                # 3. Carbon-Pi - between carbon atoms and aromatic ring atoms
                if distance <= pi_stack_cutoff:
                    if is_prot_ring_atom and lig_atom.element == 'C':
                        interaction_info = {
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        }
                        all_interactions['carbon_pi'].append(interaction_info)
                
                # 4. Donor-Pi - negatively charged residues with pi systems
                if distance <= pi_stack_cutoff:
                    if prot_res.resname in neg_charged and is_lig_aromatic:
                        interaction_info = {
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        }
                        all_interactions['donor_pi'].append(interaction_info)
                
                # 5. Amide-Pi - amide groups with pi systems
                if distance <= pi_stack_cutoff:
                    if prot_res.resname in amide_residues and is_lig_aromatic:
                        interaction_info = {
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        }
                        all_interactions['amide_pi'].append(interaction_info)
                
                # 6. Hydrophobic interactions - carbon atoms in hydrophobic residues
                if distance <= hydrophobic_cutoff:
                    if prot_res.resname in hydrophobic_residues and \
                    lig_atom.element == 'C' and prot_atom.element == 'C':
                        interaction_info = {
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        }
                        all_interactions['hydrophobic'].append(interaction_info)
                
                # 7. Ionic/salt bridge interactions - charged residues and groups
                if distance <= ionic_cutoff:
                    is_prot_pos = prot_res.resname in pos_charged and prot_atom.element == 'N'
                    is_prot_neg = prot_res.resname in neg_charged and prot_atom.element == 'O'
                    
                    if (is_prot_pos and is_lig_neg_charged) or (is_prot_neg and is_lig_pos_charged):
                        interaction_info = {
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        }
                        # A salt bridge IS the charge-charge contact detected
                        # here, not an additional one.  Appending to both keys
                        # made every ionic contact score twice (-1.80 and
                        # -1.50).  Record once, under 'ionic'; the report
                        # labels this class "Ionic / salt bridge".
                        all_interactions['ionic'].append(interaction_info)
                
                # 8. Cation-Pi - positive charged residues with aromatic systems
                # Directional and disjoint from 'pi_cation' above, which
                # handles the protein-cation -> ligand-ring direction.  This
                # branch is ligand-cation -> protein-ring only; previously it
                # also matched the protein-cation direction, so those contacts
                # were scored under both keys.
                if distance <= pi_stack_cutoff:
                    if is_prot_ring_atom and is_lig_pos_charged:
                        interaction_info = {
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        }
                        all_interactions['cation_pi'].append(interaction_info)
                
                # 9. Covalent bonds - very close interactions with specific atom types
                # Lower bound 1.2 A excludes physically impossible contacts.
                # The ligand atom must be a plausible organic electrophile: a
                # metal at 2.0-2.1 A from HIS ND1/NE2 or CYS SG is coordination,
                # not a covalent bond, and misfiling it costs -8.0 kcal/mol in
                # the scoring function.
                if 1.2 <= distance <= covalent_cutoff and \
                   lig_atom.element.strip().capitalize() in COVALENT_LIGAND_ELEMENTS:
                    # Only include likely covalent bonds involving specific residues and atom types
                    if ((prot_res.resname == 'CYS' and prot_atom.name == 'SG') or
                        (prot_res.resname == 'SER' and prot_atom.name == 'OG') or
                        (prot_res.resname == 'LYS' and prot_atom.name == 'NZ') or
                        (prot_res.resname == 'HIS' and prot_atom.name in ['ND1', 'NE2'])):
                        interaction_info = {
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        }
                        all_interactions['covalent'].append(interaction_info)
        
        # Handle halogen bonds separately
        for halogen_atom in self.halogen_atoms:
            hal_element = halogen_atom.element.strip().capitalize()
            if self.strict_halogens and hal_element not in STRICT_HALOGEN_ELEMENTS:
                continue
            nearby_atoms = ns.search(halogen_atom.get_coord(), halogen_bond_cutoff)
            for prot_atom in nearby_atoms:
                prot_res = prot_atom.get_parent()
                distance = halogen_atom - prot_atom

                if 2.5 < distance <= halogen_bond_cutoff and is_halogen_acceptor(prot_atom):
                    res_id = (prot_res.resname, prot_res.id[1])
                    self.interacting_residues.add(res_id)

                    interaction_info = {
                        'ligand_atom': halogen_atom,
                        'protein_atom': prot_atom,
                        'protein_residue': prot_res,
                        'distance': distance,
                        'halogen_element': hal_element
                    }
                    all_interactions['halogen_bonds'].append(interaction_info)
        
        # Handle metal coordination
        for metal_atom in self.metal_atoms:
            # Only a metal that is itself bound to the ligand mediates a
            # protein-ligand contact; a structural ion elsewhere in the
            # protein (e.g. a zinc finger 30 A away) is not an interaction.
            if min(metal_atom - lig_atom for lig_atom in self.ligand_atoms) > METAL_LIGAND_CUTOFF:
                continue
            for prot_atom in self.protein_atoms:
                prot_res = prot_atom.get_parent()
                distance_squared = sum((a-b)**2 for a, b in zip(metal_atom.get_coord(), prot_atom.get_coord()))
                
                if distance_squared <= metal_coord_cutoff**2 and is_metal_coordinator(prot_atom):
                    distance = math.sqrt(distance_squared)
                    res_id = (prot_res.resname, prot_res.id[1])
                    self.interacting_residues.add(res_id)
                    
                    interaction_info = {
                        'ligand_atom': metal_atom,
                        'protein_atom': prot_atom,
                        'protein_residue': prot_res,
                        'distance': distance
                    }
                    all_interactions['metal_coordination'].append(interaction_info)
        
        # ---- water bridges: ligand N/O ... water [... water] ... protein N/O
        def _angle_ok(a, centre, b):
            u, v = a.get_coord() - centre.get_coord(), b.get_coord() - centre.get_coord()
            cos = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
            omega = math.degrees(math.acos(max(-1.0, min(1.0, cos))))
            if not self.water_angle_check:
                return True, omega
            return WATER_BRIDGE_ANGLE[0] <= omega <= WATER_BRIDGE_ANGLE[1], omega

        def _hbond_range(a, b):
            return WATER_BRIDGE_MIN <= a - b <= WATER_BRIDGE_MAX

        def _protein_partners(water, previous):
            """Protein N/O hydrogen-bonded to `water` at an allowed angle."""
            for prot_atom in ns.search(water.get_coord(), WATER_BRIDGE_MAX):
                if prot_atom.element in ('N', 'O') and _hbond_range(prot_atom, water):
                    ok, omega = _angle_ok(previous, water, prot_atom)
                    if ok:
                        yield prot_atom, omega

        def _bridge(lig_atom, prot_atom, waters, legs, omega):
            prot_res = prot_atom.get_parent()
            self.interacting_residues.add((prot_res.resname, prot_res.id[1]))
            names = " + ".join(f"{w.get_parent().resname.strip()} {w.get_parent().id[1]}"
                               for w in waters)
            all_interactions['water_bridges'].append({
                'ligand_atom': lig_atom,
                'protein_atom': prot_atom,
                'protein_residue': prot_res,
                # the longest hydrogen bond limits the bridge
                'distance': max(legs),
                'distance_protein_water': legs[-1],
                'distance_ligand_water': legs[0],
                'angle': omega,
                'water_residue': waters[0].get_parent(),
                'water_residues': [w.get_parent() for w in waters],
                'n_waters': len(waters),
                'type': f"via {names} (" + " + ".join(f"{d:.2f}" for d in legs) + " Å)",
            })

        water_ns = NeighborSearch(self.water_atoms) if self.water_atoms else None
        for water in self.water_atoms:
            lig_side = [(a, water - a) for a in self.ligand_atoms
                        if a.element in ('N', 'O') and _hbond_range(a, water)]
            for lig_atom, d_lig in lig_side:
                # one water: ligand ... W ... protein
                for prot_atom, omega in _protein_partners(water, lig_atom):
                    _bridge(lig_atom, prot_atom, [water], [d_lig, prot_atom - water], omega)
                # two waters: ligand ... W1 ... W2 ... protein (each water
                # obeys the same distance and angle rules)
                for water2 in water_ns.search(water.get_coord(), WATER_BRIDGE_MAX):
                    if water2 is water or not _hbond_range(water2, water) \
                            or not _angle_ok(lig_atom, water, water2)[0]:
                        continue
                    for prot_atom, omega in _protein_partners(water2, water):
                        _bridge(lig_atom, prot_atom, [water, water2],
                                [d_lig, water2 - water, prot_atom - water2], omega)

        self.interactions = self._select_contacts(all_interactions)

    # Types where each distinct atom pair is a separate, directional contact:
    # a residue can make two hydrogen bonds, or coordinate a metal twice.
    _ATOM_PAIR_TYPES = ('hydrogen_bonds', 'halogen_bonds', 'metal_coordination', 'covalent')
    # Carbon-carbon contact classes, most specific first.  One atom pair is
    # reported under the first class that matches and not again further down.
    _CARBON_CONTACT_ORDER = ('pi_pi_stacking', 'alkyl_pi', 'carbon_pi', 'hydrophobic')

    @staticmethod
    def _residue_key(contact):
        res = contact['protein_residue']
        parent = res.get_parent()
        return (parent.id if parent is not None else '', res.resname, res.id)

    def _select_contacts(self, all_interactions):
        """
        Reduce the raw atom-atom hits to the contacts that are reported.

        * Directional, atom-specific types (hydrogen bonds, halogen bonds,
          metal coordination, covalent bonds): every distinct pair is kept,
          closest first, unless both of its atoms already take part in a
          shorter contact of the same type with the same residue.  A bidentate
          arginine therefore reports two hydrogen bonds, not one and not four.
        * All other types describe a group-level contact: the closest atom
          pair per residue is kept.
        * The same atom pair is never reported under two carbon-contact
          classes (for example alkyl-pi and hydrophobic).
        """
        selected = {}
        for itype, hits in all_interactions.items():
            by_residue = {}
            for hit in hits:
                by_residue.setdefault(self._residue_key(hit), []).append(hit)
            kept = []
            for group in by_residue.values():
                group.sort(key=lambda h: h['distance'])
                if itype in self._ATOM_PAIR_TYPES:
                    used_prot, used_lig = set(), set()
                    for hit in group:
                        p, l = id(hit['protein_atom']), id(hit['ligand_atom'])
                        if p in used_prot and l in used_lig:
                            continue
                        used_prot.add(p)
                        used_lig.add(l)
                        kept.append(hit)
                else:
                    kept.append(group[0])
            selected[itype] = kept

        seen_pairs = set()
        for itype in self._CARBON_CONTACT_ORDER:
            unique = []
            for hit in selected.get(itype, []):
                pair = (id(hit['protein_atom']), id(hit['ligand_atom']))
                if pair in seen_pairs:
                    continue
                seen_pairs.add(pair)
                unique.append(hit)
            selected[itype] = unique
        return selected

    def calculate_realistic_solvent_accessibility(self, probe_radius=1.4, exposure_threshold=0.25, max_percent=0.5):
        """
        Realistic solvent accessibility calculation with proper constraints on number of accessible residues.
        
        Parameters:
        -----------
        probe_radius : float
            Radius of solvent probe in Angstroms (default: 1.4)
        exposure_threshold : float
            Threshold ratio for considering a residue solvent accessible (default: 0.25)
        max_percent : float
            Maximum percentage of interacting residues that can be solvent accessible (default: 0.5)
        """
        self.solvent_accessible = set()

        # Define which residues are typically surface-exposed
        likely_exposed = {'ARG', 'LYS', 'ASP', 'GLU', 'ASN', 'GLN', 'HIS', 'SER', 'THR', 'TYR'}
        likely_buried = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'CYS', 'PRO'}

        # First get all protein atoms (including non-interacting ones)
        all_protein_atoms = []
        for residue in self.model.get_residues():
            if residue.id[0] == ' ':  # Standard amino acid
                for atom in residue:
                    all_protein_atoms.append(atom)
        
        # Calculate protein center
        protein_center = np.zeros(3)
        for atom in all_protein_atoms:
            protein_center += atom.get_coord()
        protein_center /= len(all_protein_atoms) if all_protein_atoms else 1
        
        # Store exposure scores for all residues to sort later
        exposure_scores = {}
        
        # For each interacting residue, estimate accessibility
        for res_id in self.interacting_residues:
            residue = self.protein_residues.get(res_id)
            if residue is None:
                continue
            
            # Calculate residue center
            residue_center = np.zeros(3)
            residue_atoms = list(residue.get_atoms())
            for atom in residue_atoms:
                residue_center += atom.get_coord()
            residue_center /= len(residue_atoms) if residue_atoms else 1
            
            # Calculate vector from protein center to residue center
            direction = residue_center - protein_center
            direction_length = np.linalg.norm(direction)
            if direction_length > 0:
                direction = direction / direction_length
            
            # Count exposed atoms
            exposed_atoms = 0
            total_atoms = len(residue_atoms)
            
            for atom in residue_atoms:
                atom_coord = atom.get_coord()
                
                # Check if atom is on the protein surface
                is_exposed = True
                nearby_atom_count = 0
                
                for other_atom in all_protein_atoms:
                    if other_atom.get_parent() == residue:
                        continue  # Skip atoms in same residue
                    
                    distance = np.linalg.norm(atom_coord - other_atom.get_coord())
                    
                    # Simple distance threshold for exposure
                    if distance < 3.0:
                        nearby_atom_count += 1
                    
                    # If more than 8 atoms are nearby (instead of 12), consider it buried
                    if nearby_atom_count > 8:
                        is_exposed = False
                        break
                
                if is_exposed:
                    exposed_atoms += 1
            
            # Calculate exposure ratio
            exposure_ratio = exposed_atoms / total_atoms if total_atoms > 0 else 0
            
            # Calculate more strict surface bias
            if residue.resname in likely_exposed:
                surface_bias = 1.5  # More bias for typically exposed residues
            elif residue.resname in likely_buried:
                surface_bias = 0.6  # Much lower bias for typically buried residues
            else:
                surface_bias = 1.0
                
            # Calculate distance from surface bias
            # Residues further from center are more likely exposed
            distance_from_center_ratio = min(direction_length / 15.0, 1.0)
            
            # Combined score for exposure - more weight on actual exposure ratio
            exposure_score = (exposure_ratio * 2.0 + surface_bias + distance_from_center_ratio) / 4.0
            
            # Store score for later ranking
            exposure_scores[res_id] = exposure_score
            
            if exposure_score > exposure_threshold:
                self.solvent_accessible.add(res_id)

        # Enforce min/max constraints
        min_expected = max(1, int(len(self.interacting_residues) * 0.1))
        max_expected = min(int(len(self.interacting_residues) * max_percent),
                           len(self.interacting_residues) - 1)

        if len(self.solvent_accessible) < min_expected:
            remaining = sorted(
                [(r, exposure_scores.get(r, 0.0)) for r in self.interacting_residues
                 if r not in self.solvent_accessible],
                key=lambda x: (x[1], x[0]), reverse=True
            )
            for res_id, _ in remaining:
                if len(self.solvent_accessible) >= min_expected:
                    break
                self.solvent_accessible.add(res_id)

        if len(self.solvent_accessible) > max_expected:
            to_evaluate = sorted(
                [(r, exposure_scores.get(r, 0.0)) for r in self.solvent_accessible],
                key=lambda x: (x[1], x[0])
            )
            for res_id, _ in to_evaluate[:len(self.solvent_accessible) - max_expected]:
                self.solvent_accessible.discard(res_id)

        return self.solvent_accessible

    def calculate_enhanced_solvent_accessibility(self, probe_radius=1.4, exposure_threshold=0.15):
        """
        Enhanced solvent accessibility calculation with better debugging and less stringent criteria.
        
        Parameters:
        -----------
        probe_radius : float
            Radius of solvent probe in Angstroms (default: 1.4)
        exposure_threshold : float
            Threshold ratio for considering a residue solvent accessible (default: 0.15)
        """
        self.solvent_accessible = set()

        # Define which residues are typically surface-exposed
        likely_exposed = {'ARG', 'LYS', 'ASP', 'GLU', 'ASN', 'GLN', 'HIS', 'SER', 'THR', 'TYR'}

        # First get all protein atoms (including non-interacting ones)
        all_protein_atoms = []
        for residue in self.model.get_residues():
            if residue.id[0] == ' ':  # Standard amino acid
                for atom in residue:
                    all_protein_atoms.append(atom)
        
        # Calculate protein center
        protein_center = np.zeros(3)
        for atom in all_protein_atoms:
            protein_center += atom.get_coord()
        protein_center /= len(all_protein_atoms) if all_protein_atoms else 1
        
        # For each interacting residue, estimate accessibility
        for res_id in self.interacting_residues:
            residue = self.protein_residues.get(res_id)
            if residue is None:
                continue
            
            # Calculate residue center
            residue_center = np.zeros(3)
            residue_atoms = list(residue.get_atoms())
            for atom in residue_atoms:
                residue_center += atom.get_coord()
            residue_center /= len(residue_atoms) if residue_atoms else 1
            
            # Calculate vector from protein center to residue center
            direction = residue_center - protein_center
            direction_length = np.linalg.norm(direction)
            if direction_length > 0:
                direction = direction / direction_length
            
            # Count exposed atoms
            exposed_atoms = 0
            total_atoms = len(residue_atoms)
            
            for atom in residue_atoms:
                atom_coord = atom.get_coord()
                
                # Check if atom is on the protein surface
                is_exposed = True
                nearby_atom_count = 0
                
                for other_atom in all_protein_atoms:
                    if other_atom.get_parent() == residue:
                        continue  # Skip atoms in same residue
                    
                    distance = np.linalg.norm(atom_coord - other_atom.get_coord())
                    
                    # Simple distance threshold for exposure
                    if distance < 3.0:
                        nearby_atom_count += 1
                    
                    # If more than 12 atoms are nearby, consider it buried
                    if nearby_atom_count > 12:
                        is_exposed = False
                        break
                
                if is_exposed:
                    exposed_atoms += 1
            
            # Calculate exposure ratio and make decision
            exposure_ratio = exposed_atoms / total_atoms if total_atoms > 0 else 0
            
            # Three criteria for solvent accessibility:
            # 1. Reasonable exposure ratio
            # 2. Residue type typically found on protein surface
            # 3. Residue closer to protein surface than center
            
            surface_bias = 1.2 if residue.resname in likely_exposed else 0.8
            distance_from_center_ratio = direction_length / 10.0  # Normalize to ~0-1 scale
            
            # Combined score for exposure
            exposure_score = exposure_ratio * surface_bias * min(1.5, distance_from_center_ratio)
            
            if exposure_score > exposure_threshold:
                self.solvent_accessible.add(res_id)

        min_expected = max(2, int(len(self.interacting_residues) * 0.15))
        if len(self.solvent_accessible) < min_expected:
            remaining = sorted(
                [r for r in self.interacting_residues if r not in self.solvent_accessible],
                key=lambda r: 2 if r[0] in likely_exposed else 1,
                reverse=True
            )
            for res_id in remaining:
                if len(self.solvent_accessible) >= min_expected:
                    break
                self.solvent_accessible.add(res_id)

        return self.solvent_accessible
    
    def calculate_dssp_solvent_accessibility(self, dssp_executable='dssp'):
        """
        Calculate solvent accessibility using DSSP.
        Requires DSSP executable to be installed and in PATH.
        
        Parameters:
        -----------
        dssp_executable : str
            Path to DSSP executable (default: 'dssp')
            
        Returns:
        --------
        dict
            Dictionary mapping (resname, resnum) to relative solvent accessibility (0-1)
        """
        self.solvent_accessible = set()
        
        try:
            import warnings as _warnings
            # Create a temporary PDB file for DSSP input
            with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp_pdb:
                pdb_io = PDBIO()
                pdb_io.set_structure(self.structure)
                pdb_io.save(tmp_pdb.name)

                # Newer BioPython DSSP emits a harmless mmCIF parse warning
                # when given a PDB-format temp file — suppress it cleanly.
                with _warnings.catch_warnings():
                    _warnings.simplefilter('ignore')
                    dssp = DSSP(self.model, tmp_pdb.name, dssp=dssp_executable)
                
                # Clean up temporary file
                try:
                    os.unlink(tmp_pdb.name)
                except:
                    pass
                
                # Process DSSP results
                for (chain_id, res_id), dssp_data in dssp.property_dict.items():
                    resname = dssp_data[0]
                    resnum = res_id[1]
                    res_key = (resname, resnum)
                    
                    # Get relative solvent accessibility (0-1)
                    rsa = dssp_data[3]  # Relative accessibility
                    
                    # Consider residues with >15% accessibility as solvent accessible
                    if rsa > 0.25 and res_key in self.interacting_residues:
                        self.solvent_accessible.add(res_key)
                        
        except Exception as e:
            logger.warning(f"Warning: DSSP calculation failed. Falling back to geometric estimation. Error: {str(e)}")
            self.estimate_solvent_accessibility()
        
        return self.solvent_accessible
    
    def _ligand_neighbours(self, atom, cutoff=1.75):
        """Heavy-atom covalent neighbours of a ligand atom."""
        out = []
        for other in self.ligand_atoms:
            if other is atom or other.element.strip().upper() == 'H':
                continue
            if np.linalg.norm(other.get_coord() - atom.get_coord()) < cutoff:
                out.append(other)
        return out

    def _is_ligand_anion(self, atom):
        """
        Formally anionic ligand oxygen.

        The previous rule ("oxygen with no carbon neighbour") excluded every
        carboxylate, which is the single most common anionic group in
        drug-like ligands, so carboxylate salt bridges were never detected.
        Here an oxygen counts as anionic if it is terminal on a carboxylate,
        phosphate/phosphonate, or sulfate/sulfonate group.
        """
        if atom.element.strip().upper() != 'O':
            return False
        neighbours = self._ligand_neighbours(atom)
        if len(neighbours) != 1:
            return False                      # bridging or hydroxyl-like
        centre = neighbours[0]
        centre_el = centre.element.strip().upper()
        if centre_el not in ('C', 'P', 'S', 'N'):
            return False
        # Count terminal oxygens on the central atom: >=2 implies a
        # delocalised, formally charged group (COO-, PO3-, SO3-, NO2-).
        terminal_o = sum(
            1 for nb in self._ligand_neighbours(centre)
            if nb.element.strip().upper() == 'O' and len(self._ligand_neighbours(nb)) == 1
        )
        return terminal_o >= 2

    def _is_ligand_cation(self, atom):
        """
        Formally cationic ligand nitrogen.

        Counts quaternary/protonated amines (three or more heavy-atom
        neighbours) and guanidinium/amidinium nitrogens (terminal N on a
        carbon bearing two or more nitrogens).
        """
        if atom.element.strip().upper() != 'N':
            return False
        neighbours = self._ligand_neighbours(atom)
        if len(neighbours) >= 3:
            return True                        # quaternary / protonated amine
        for nb in neighbours:
            if nb.element.strip().upper() != 'C':
                continue
            n_count = sum(1 for x in self._ligand_neighbours(nb)
                          if x.element.strip().upper() == 'N')
            if n_count >= 2:
                return True                    # guanidinium / amidinium
        return False

    def estimate_solvent_accessibility(self, probe_radius=1.4):
        """
        Tier-2 solvent-accessibility estimate (dependency-free).

        Previously every call site referenced this name but no such method
        existed, so the DSSP and pure-Python tiers raised AttributeError and
        silently degraded to the geometric estimator.  Routing it to the
        Shrake-Rupley implementation makes the documented three-tier fallback
        actually reachable.
        """
        return self.calculate_python_solvent_accessibility(probe_radius=probe_radius)

    def calculate_python_solvent_accessibility(self, probe_radius=1.4):
        """
        Enhanced solvent accessibility calculation in pure Python.
        Based on Shrake-Rupley algorithm but with practical approximations.
        
        Parameters:
        -----------
        probe_radius : float
            Radius of solvent probe in Angstroms (default: 1.4)
        """
        self.solvent_accessible = set()
        
        # First get all protein atoms (including non-interacting ones)
        all_protein_atoms = []
        for residue in self.model.get_residues():
            if residue.id[0] == ' ':  # Standard amino acid
                for atom in residue:
                    all_protein_atoms.append(atom)
        
        # For each interacting residue, estimate accessibility
        for res_id in self.interacting_residues:
            residue = self.protein_residues.get(res_id)
            if residue is None:
                continue
                
            # Count how many atoms are exposed
            exposed_atoms = 0
            total_atoms = 0
            residue_atoms = list(residue.get_atoms())
            
            for atom in residue_atoms:
                total_atoms += 1
                atom_coord = atom.get_coord()
                
                # Get atom radius (approximated)
                atom_radius = 1.8  # Default carbon radius
                if atom.element == 'O':
                    atom_radius = 1.4
                elif atom.element == 'N':
                    atom_radius = 1.5
                elif atom.element == 'S':
                    atom_radius = 1.8
                
                # Check if atom is buried by other protein atoms
                is_exposed = True
                nearby_atom_count = 0
                
                for other_atom in all_protein_atoms:
                    if other_atom.get_parent() == residue:
                        continue  # Skip atoms in same residue
                    
                    distance = np.linalg.norm(atom_coord - other_atom.get_coord())
                    other_radius = 1.8  # Default
                    
                    if distance < (atom_radius + other_radius + probe_radius):
                        nearby_atom_count += 1
                        if nearby_atom_count > 15:  # If too many atoms are nearby, consider it buried
                            is_exposed = False
                            break
                
                if is_exposed:
                    exposed_atoms += 1
            
            # Consider residue accessible if >20% of its atoms are exposed
            exposure_ratio = exposed_atoms / total_atoms if total_atoms > 0 else 0
            if exposure_ratio > 0.2:
                self.solvent_accessible.add(res_id)
        
        # Check if reasonable number are marked accessible
        if len(self.solvent_accessible) < max(2, len(self.interacting_residues) // 4):
            # If too few, drop to the tier-3 geometric estimator.  (Calling
            # estimate_solvent_accessibility() here would recurse, since that
            # method now delegates to this one.)
            return self.calculate_realistic_solvent_accessibility()

        return self.solvent_accessible
    # Future enhancements for solvent accessibility calculation

# 1. Add SASA (Solvent Accessible Surface Area) calculation using Shrake-Rupley algorithm
    def calculate_sasa(self, probe_radius=1.4, n_sphere_points=100):
        """
        Calculate solvent accessible surface area using Shrake-Rupley algorithm.
        
        Parameters:
        -----------
        probe_radius : float
            Radius of solvent probe in Angstroms (default: 1.4)
        n_sphere_points : int
            Number of points on the sphere for the algorithm (default: 100)
            
        Returns:
        --------
        dict
            Dictionary mapping residue IDs to their SASA values
        """
        logger.info("Calculating Solvent Accessible Surface Areas...")
        
        # Step 1: Generate points on a unit sphere (Fibonacci sphere method)
        points = []
        phi = (1 + math.sqrt(5)) / 2  # Golden ratio
        
        for i in range(n_sphere_points):
            y = 1 - (i / float(n_sphere_points - 1)) * 2
            radius = math.sqrt(1 - y*y)
            
            theta = 2 * math.pi * i / phi
            x = math.cos(theta) * radius
            z = math.sin(theta) * radius
            
            points.append(np.array([x, y, z]))
        
        # Step 2: Calculate SASA for each residue
        residue_sasa = {}
        
        for res_id in self.interacting_residues:
            residue = self.protein_residues.get(res_id)
            if residue is None:
                continue
            
            total_area = 0
            for atom in residue.get_atoms():
                # Get atom radius based on element
                element = atom.element
                if element == 'H':
                    radius = 1.2
                elif element == 'C':
                    radius = 1.7
                elif element == 'N':
                    radius = 1.55
                elif element == 'O':
                    radius = 1.52
                elif element == 'S':
                    radius = 1.8
                else:
                    radius = 1.7  # Default
                
                # Add probe radius to atom radius
                expanded_radius = radius + probe_radius
                atom_coord = atom.get_coord()
                
                # Place points on the expanded sphere around the atom
                sphere_points = [atom_coord + p * expanded_radius for p in points]
                
                # Count points not within any other atom
                accessible_points = 0
                for point in sphere_points:
                    is_accessible = True
                    
                    for other_atom in self.protein_atoms:
                        # Skip the current atom and atoms in the same residue
                        if other_atom == atom or other_atom.get_parent() == residue:
                            continue
                        
                        other_coord = other_atom.get_coord()
                        distance = np.linalg.norm(point - other_coord)
                        
                        # Get radius of other atom
                        other_element = other_atom.element
                        if other_element == 'H':
                            other_radius = 1.2
                        elif other_element == 'C':
                            other_radius = 1.7
                        elif other_element == 'N':
                            other_radius = 1.55
                        elif other_element == 'O':
                            other_radius = 1.52
                        elif other_element == 'S':
                            other_radius = 1.8
                        else:
                            other_radius = 1.7
                        
                        # Check if point is within other atom
                        if distance < other_radius + probe_radius:
                            is_accessible = False
                            break
                    
                    if is_accessible:
                        accessible_points += 1
                
                # Calculate area of this atom
                point_ratio = accessible_points / len(sphere_points)
                atom_area = 4 * math.pi * expanded_radius**2 * point_ratio
                total_area += atom_area
            
            # Store total SASA for this residue
            residue_sasa[res_id] = total_area
        
        # Step 3: Determine solvent accessibility based on SASA thresholds
        # (Typically: >20Å² for small residues, >40Å² for large residues)
        self.solvent_accessible = set()
        for res_id, sasa in residue_sasa.items():
            resname = res_id[0]
            
            # Define threshold based on residue size
            if resname in {'GLY', 'ALA', 'SER', 'CYS', 'THR'}:
                threshold = 20  # Small residues
            elif resname in {'ARG', 'LYS', 'TRP', 'TYR', 'PHE'}:
                threshold = 40  # Large residues
            else:
                threshold = 30  # Medium residues
            
            # Mark as solvent accessible if above threshold
            if sasa >= threshold:
                self.solvent_accessible.add(res_id)
        
        # Return SASA values for reference
        return residue_sasa

    # 2. Integration with FreeSASA or other external tools
    def calculate_external_sasa(self, executable_path=None):
                """
                Calculate solvent accessibility using external tools like FreeSASA.
                Requires FreeSASA or similar SASA calculation tool.
                
                Parameters:
                -----------
                executable_path : str, optional
                    Path to FreeSASA executable (if None, assumes it's in PATH)
                
                Returns:
                --------
                set
                    Set of solvent accessible residue IDs
                """
                if executable_path is None:
                    executable_path = "freesasa"  # assume in PATH
                
                try:
                    import tempfile
                    import subprocess
                    from Bio.PDB.PDBIO import PDBIO
                    
                    # Create temp file for PDB
                    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp_pdb:
                        pdb_path = tmp_pdb.name
                        # Write structure to PDB
                        io = PDBIO()
                        io.set_structure(self.structure)
                        io.save(pdb_path)
                    
                    # Run FreeSASA
                    cmd = [executable_path, "--format=rsa", pdb_path]
                    result = subprocess.run(cmd, capture_output=True, text=True)
                    
                    # Parse the RSA output
                    self.solvent_accessible = set()
                    
                    for line in result.stdout.split('\n'):
                        if line.startswith("RES"):
                            parts = line.split()
                            if len(parts) >= 5:
                                resname = parts[1]
                                resnum = int(parts[2])
                                rel_sasa = float(parts[4])
                                
                                # If relative SASA > 20%, consider as accessible
                                if rel_sasa > 20.0:
                                    res_id = (resname, resnum)
                                    if res_id in self.interacting_residues:
                                        self.solvent_accessible.add(res_id)
                    
                    # Cleanup
                    try:
                        os.unlink(pdb_path)
                    except:
                        pass
                        
                except Exception as e:
                    logger.warning(f"Error using external SASA tool: {str(e)}")
                    logger.info("Falling back to geometric method")
                    self.calculate_realistic_solvent_accessibility()
                
                return self.solvent_accessible    

    def _force_directed_layout(self, positions, ring_radius,
                               n_iter=120, k_spring=0.15,
                               k_repel=6000.0, k_radial=0.25):
        """
        Refine residue positions from an initial circular layout using a
        spring-embedder (Fruchterman-Reingold-style) force-directed algorithm.

        Residues that share an interaction type exert attractive spring forces
        on each other; all pairs exert repulsive forces; and a radial restoring
        force keeps every node close to *ring_radius* from the origin so the
        diagram does not collapse or explode.

        Parameters
        ----------
        positions : dict
            Mapping res_id → [x, y]  (modified in-place and returned).
        ring_radius : float
            Target distance from the origin for each node.
        n_iter : int
            Number of simulation steps.
        k_spring : float
            Spring constant for same-interaction-type attraction.
        k_repel : float
            Repulsion constant (applied to all pairs).
        k_radial : float
            Strength of the radial restoring force.

        Returns
        -------
        dict
            Updated positions mapping res_id → [x, y].
        """
        if len(positions) < 3:
            return positions  # Not worth running for 1-2 residues

        res_ids = list(positions.keys())
        pos = {r: np.array(positions[r], dtype=float) for r in res_ids}

        # Build adjacency: residues connected by a shared interaction type
        connected = set()
        for itype, ints in self.interactions.items():
            # Group residues that appear in the same interaction type
            in_type = [
                (i['protein_residue'].resname, i['protein_residue'].id[1])
                for i in ints
            ]
            for a in in_type:
                for b in in_type:
                    if a != b and a in pos and b in pos:
                        key = tuple(sorted([a, b]))
                        connected.add(key)

        temperature = ring_radius * 0.15  # initial step cap
        cooling = temperature / n_iter

        for iteration in range(n_iter):
            disp = {r: np.zeros(2) for r in res_ids}

            # Repulsive forces (all pairs)
            for i, r1 in enumerate(res_ids):
                for r2 in res_ids[i + 1:]:
                    delta = pos[r1] - pos[r2]
                    dist = max(np.linalg.norm(delta), 1.0)
                    force = k_repel / (dist * dist)
                    unit = delta / dist
                    disp[r1] += force * unit
                    disp[r2] -= force * unit

            # Attractive spring forces (connected pairs only)
            for r1, r2 in connected:
                if r1 not in pos or r2 not in pos:
                    continue
                delta = pos[r2] - pos[r1]
                dist = max(np.linalg.norm(delta), 1.0)
                force = k_spring * dist
                unit = delta / dist
                disp[r1] += force * unit
                disp[r2] -= force * unit

            # Radial restoring force — keeps nodes near ring_radius
            for r in res_ids:
                dist_from_origin = max(np.linalg.norm(pos[r]), 1.0)
                radial_error = dist_from_origin - ring_radius
                unit_out = pos[r] / dist_from_origin
                disp[r] -= k_radial * radial_error * unit_out

            # Apply displacements with temperature cap
            for r in res_ids:
                d = np.linalg.norm(disp[r])
                if d > 0:
                    pos[r] += disp[r] / d * min(d, temperature)

            temperature -= cooling

        return {r: pos[r].tolist() for r in res_ids}

    def visualize(self, output_file='protein_ligand_interactions.png', figsize=None,
                  dpi=300, title=None, show_3d_cues=True, style='modern',
                  show_distances=True, show_nearby_residues=False):
        """
        Generate a complete 2D visualization of protein-ligand interactions.

        Parameters:
        -----------
        output_file : str
            Image path; the extension selects the format (.png, .svg, .pdf).
        style : str
            'modern' (default): residues are placed next to the ligand atoms
            they contact, boxes are coloured by residue type and the legend
            sits below the plot.  'classic': the original circular PandaMap
            layout.  '3d': the 3D contact view (transparent residues around
            the ligand, contacting atoms highlighted).
        show_distances, show_nearby_residues : bool
            Modern style only: print hydrogen-bond distances; also show
            binding-site residues without a specific interaction.
        show_3d_cues : bool
            Classic style only: depth cues on the ligand.
        """
        if style not in DIAGRAM_STYLES:
            raise ValueError(f"Unknown diagram style '{style}'. "
                             f"Use one of: {', '.join(DIAGRAM_STYLES)}.")
        if style == '3d':
            return self.render_contacts_3d(output_file=output_file, dpi=dpi, title=title,
                                           show_distances=show_distances)
        if style == 'modern':
            from .diagram import draw_interaction_diagram
            return draw_interaction_diagram(
                self, output_file=output_file, figsize=figsize or (11, 11.6),
                dpi=dpi, title=title, show_distances=show_distances,
                show_nearby_residues=show_nearby_residues)
        return self._visualize_classic(output_file=output_file,
                                       figsize=figsize or (12, 12), dpi=dpi,
                                       title=title, show_3d_cues=show_3d_cues)

    def _visualize_classic(self, output_file='protein_ligand_interactions.png',
                           figsize=(12, 12), dpi=300, title=None, show_3d_cues=True):
        """The original circular-layout diagram (``style='classic'``)."""
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Validate solvent accessibility counts
        
        # Force reasonable solvent accessibility if all residues are marked or none are marked
        if len(self.solvent_accessible) == len(self.interacting_residues):
            # All residues marked accessible is implausible; filter to most-likely exposed
            
            # If all residues are marked, keep only about 40% 
            # Focus on residues that are typically exposed
            likely_exposed = {'ARG', 'LYS', 'ASP', 'GLU', 'ASN', 'GLN', 'HIS', 'SER', 'THR', 'TYR'}
            keep_residues = []
            
            for res_id in self.interacting_residues:
                if res_id[0] in likely_exposed:
                    keep_residues.append(res_id)
                    
            # If we don't have enough exposed residues, add some hydrophobic ones near the surface
            max_to_keep = max(2, int(len(self.interacting_residues) * 0.4))
            if len(keep_residues) < max_to_keep:
                for res_id in self.interacting_residues:
                    if res_id not in keep_residues and len(keep_residues) < max_to_keep:
                        keep_residues.append(res_id)
            
            # Replace the solvent_accessible set with our filtered version
            self.solvent_accessible = set(keep_residues[:max_to_keep])
            
        
        # Add light blue background for ligand
        ligand_radius = 90
        ligand_pos = (0, 0)
        ligand_circle = Circle(ligand_pos, ligand_radius, facecolor='#ADD8E6', 
                            edgecolor='none', alpha=0.4, zorder=1)
        ax.add_patch(ligand_circle)
        
        # Draw the simplified ligand structure
        atom_positions = self.ligand_structure.draw_on_axes(ax, center=ligand_pos, 
                                                        radius=ligand_radius*0.8,
                                                        show_3d_cues=show_3d_cues)
        
        # Place interacting residues in a circle around the ligand
        n_residues = len(self.interacting_residues)
        if n_residues == 0:
            logger.warning("Warning: No interacting residues detected.")
            n_residues = 1
            
        # Calculate positions for residues
        radius = 250  # Distance from center to residues
        residue_positions = {}
        rect_width, rect_height = 60, 30  # Residue box dimensions
        
        
        # Arrange residues in a circle
        for i, res_id in enumerate(sorted(self.interacting_residues)):
            angle = 2 * math.pi * i / n_residues
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)
            residue_positions[res_id] = (x, y)
            
            if res_id in self.solvent_accessible:
                solvent_circle = Circle((x, y), 40, facecolor='#ADD8E6',
                                        edgecolor='#87CEEB', alpha=0.5, zorder=1)
                ax.add_patch(solvent_circle)
            
            # Draw residue node as rectangle
            residue_box = Rectangle((x-rect_width/2, y-rect_height/2), rect_width, rect_height,
                                facecolor='white', edgecolor='black', linewidth=1.5,
                                zorder=2, alpha=1.0)
            ax.add_patch(residue_box)
            
            # Add residue label
            resname, resnum = res_id
            label = f"{resname} {resnum}"
            text = ax.text(x, y, label, ha='center', va='center',
                        fontsize=11, fontweight='bold', zorder=3)
            text.set_path_effects([path_effects.withStroke(linewidth=2, foreground='white')])


        # Define interaction styles
        interaction_styles = {
            'hydrogen_bonds': {'color': 'green', 'linestyle': '-', 'linewidth': 1.5, 
                            'marker_text': 'H', 'marker_bg': '#E0FFE0', 'name': 'Hydrogen Bond'},
            'carbon_pi': {'color': '#666666', 'linestyle': '--', 'linewidth': 1.5,
                        'marker_text': 'C-π', 'marker_bg': 'white', 'name': 'Carbon-Pi'},
            'pi_pi_stacking': {'color': '#9370DB', 'linestyle': '--', 'linewidth': 1.5,
                            'marker_text': 'π-π', 'marker_bg': 'white', 'name': 'Pi-Pi'},
            'donor_pi': {'color': '#FF69B4', 'linestyle': '--', 'linewidth': 1.5,
                        'marker_text': 'D', 'marker_bg': 'white', 'name': 'Donor-Pi'},
            'amide_pi': {'color': '#A52A2A', 'linestyle': '--', 'linewidth': 1.5,
                        'marker_text': 'A', 'marker_bg': 'white', 'name': 'Amide-Pi'},
            'hydrophobic': {'color': '#808080', 'linestyle': ':', 'linewidth': 1.0,
                        'marker_text': 'h', 'marker_bg': 'white', 'name': 'Hydrophobic'},
            'ionic': {'color': '#FF4500', 'linestyle': '-', 'linewidth': 1.5,
                    'marker_text': 'I', 'marker_bg': '#FFE4E1', 'name': 'Ionic'},
            'halogen_bonds': {'color': '#00CED1', 'linestyle': '-', 'linewidth': 1.5,
                        'marker_text': 'X', 'marker_bg': '#E0FFFF', 'name': 'Halogen Bond'},
            'cation_pi': {'color': '#FF00FF', 'linestyle': '--', 'linewidth': 1.5,
                    'marker_text': 'C+π', 'marker_bg': 'white', 'name': 'Cation-Pi'},
            'metal_coordination': {'color': '#FFD700', 'linestyle': '-', 'linewidth': 1.5,
                            'marker_text': 'M', 'marker_bg': '#FFFACD', 'name': 'Metal Coordination'},
            'salt_bridge': {'color': '#FF6347', 'linestyle': '-', 'linewidth': 1.5,
                        'marker_text': 'S', 'marker_bg': '#FFEFD5', 'name': 'Salt Bridge'},
            'covalent': {'color': '#000000', 'linestyle': '-', 'linewidth': 2.0,
                    'marker_text': 'COV', 'marker_bg': '#FFFFFF', 'name': 'Covalent Bond'},
            # NEW interaction styles
            'alkyl_pi': {'color': '#4682B4', 'linestyle': '--', 'linewidth': 1.5,
                    'marker_text': 'A-π', 'marker_bg': 'white', 'name': 'Alkyl-Pi'},
            'attractive_charge': {'color': '#1E90FF', 'linestyle': '-', 'linewidth': 1.5,
                            'marker_text': 'A+', 'marker_bg': '#E6E6FA', 'name': 'Attractive Charge'},
            'pi_cation': {'color': '#FF00FF', 'linestyle': '--', 'linewidth': 1.5,
                    'marker_text': 'π-C+', 'marker_bg': 'white', 'name': 'Pi-Cation'},
            'repulsion': {'color': '#DC143C', 'linestyle': '-', 'linewidth': 1.5,
                    'marker_text': 'R', 'marker_bg': '#FFC0CB', 'name': 'Repulsion'}
        }
    

        # Function to find box edge intersection
        def find_box_edge(box_center, target_point, width, height):
            """Find where a line from box center to target point intersects the box edge"""
            dx = target_point[0] - box_center[0]
            dy = target_point[1] - box_center[1]
            angle = math.atan2(dy, dx)
            
            half_width = width/2
            half_height = height/2
            
            if abs(dx) > abs(dy):
                x_intersect = box_center[0] + (half_width if dx > 0 else -half_width)
                y_intersect = box_center[1] + (x_intersect - box_center[0]) * dy/dx
                if abs(y_intersect - box_center[1]) > half_height:
                    y_intersect = box_center[1] + (half_height if dy > 0 else -half_height)
                    x_intersect = box_center[0] + (y_intersect - box_center[1]) * dx/dy
            else:
                y_intersect = box_center[1] + (half_height if dy > 0 else -half_height)
                x_intersect = box_center[0] + (y_intersect - box_center[1]) * dx/dy
                if abs(x_intersect - box_center[0]) > half_width:
                    x_intersect = box_center[0] + (half_width if dx > 0 else -half_width)
                    y_intersect = box_center[1] + (x_intersect - box_center[0]) * dy/dx
                    
            return (x_intersect, y_intersect)

        # Store interaction lines for marker placement
        interaction_lines = []

        # Track how many lines have already been drawn to each residue so we can
        # stagger curvature and prevent marker overlap when multiple interaction
        # types share the same target residue.
        residue_line_count = {}

        # Draw interaction lines with arrows
        for interaction_type, interactions in self.interactions.items():
            if interaction_type not in interaction_styles:
                continue

            style = interaction_styles[interaction_type]

            for interaction in interactions:
                res = interaction['protein_residue']
                res_id = (res.resname, res.id[1])
                lig_atom = interaction['ligand_atom']

                if res_id not in residue_positions:
                    continue

                res_pos = residue_positions[res_id]

                # Get ligand atom position
                if lig_atom.get_id() in atom_positions:
                    lig_pos = atom_positions[lig_atom.get_id()]
                else:
                    dx = res_pos[0] - ligand_pos[0]
                    dy = res_pos[1] - ligand_pos[1]
                    angle = math.atan2(dy, dx)
                    lig_pos = (ligand_pos[0] + ligand_radius * math.cos(angle),
                               ligand_pos[1] + ligand_radius * math.sin(angle))

                # Find box edge intersection
                box_edge_pos = find_box_edge(res_pos, lig_pos, rect_width, rect_height)

                # Base curvature from distance; stagger each additional line to
                # the same residue by ±0.15 so arcs fan out and markers don't pile up.
                # Formula: idx=0 → base, idx=1 → base-0.15, idx=2 → base+0.15, ...
                dx = res_pos[0] - lig_pos[0]
                dy = res_pos[1] - lig_pos[1]
                distance = math.sqrt(dx*dx + dy*dy)
                base_curvature = 0.08 * (200 / max(distance, 100))
                line_idx = residue_line_count.get(res_id, 0)
                residue_line_count[res_id] = line_idx + 1
                sign = 1 if line_idx % 2 == 0 else -1
                offset = ((line_idx + 1) // 2) * 0.15
                curvature = base_curvature + sign * offset
                
                # Store line parameters
                line_params = {
                    'start_pos': box_edge_pos,
                    'end_pos': lig_pos,
                    'curvature': curvature,
                    'style': style,
                    'interaction_type': interaction_type,
                    'res_id': res_id,
                    'key': f"{interaction_type}_{res_id[0]}_{res_id[1]}",
                    'distance': distance
                }
                interaction_lines.append(line_params)
                
                # Get directionality info
                interaction_key = (res_id, interaction_type)
                direction = self.interaction_direction.get(interaction_key, 'bidirectional')
                
                # Draw appropriate arrow based on direction
                if direction == 'protein_to_ligand':
                    # Arrow from protein to ligand
                    arrow = FancyArrowPatch(
                        box_edge_pos, lig_pos,
                        connectionstyle=f"arc3,rad={curvature}",
                        color=style['color'],
                        linestyle=style['linestyle'],
                        linewidth=style['linewidth'],
                        arrowstyle='-|>',
                        mutation_scale=10,
                        alpha=0.7,
                        zorder=4
                    )
                    ax.add_patch(arrow)
                    
                elif direction == 'ligand_to_protein':
                    # Arrow from ligand to protein
                    arrow = FancyArrowPatch(
                        lig_pos, box_edge_pos,
                        connectionstyle=f"arc3,rad={-curvature}",
                        color=style['color'],
                        linestyle=style['linestyle'],
                        linewidth=style['linewidth'],
                        arrowstyle='-|>',
                        mutation_scale=10,
                        alpha=0.7,
                        zorder=4
                    )
                    ax.add_patch(arrow)
                    
                else:  # bidirectional - single line with arrows on both ends
                    # Draw one bidirectional arrow
                    arrow = FancyArrowPatch(
                        box_edge_pos, lig_pos,
                        connectionstyle=f"arc3,rad={curvature}",
                        color=style['color'],
                        linestyle=style['linestyle'],
                        linewidth=style['linewidth'],
                        arrowstyle='<|-|>',  # Arrows on both ends
                        mutation_scale=10,
                        alpha=0.7,
                        zorder=4
                    )
                    ax.add_patch(arrow)

        # Place markers along interaction lines
        marker_positions = {}
        type_order = {'hydrogen_bonds': 0, 'ionic': 1, 'salt_bridge': 2, 'halogen_bonds': 3,
                     'metal_coordination': 4, 'pi_pi_stacking': 5, 'cation_pi': 6, 
                     'carbon_pi': 7, 'donor_pi': 8, 'amide_pi': 9, 'hydrophobic': 10}
        
        sorted_lines = sorted(interaction_lines,
                            key=lambda x: (type_order.get(x['interaction_type'], 999), x['distance']))
        
        for line_params in sorted_lines:
            start_pos = line_params['start_pos']
            end_pos = line_params['end_pos']
            curvature = line_params['curvature']
            style = line_params['style']
            key = line_params['key']
            res_id = line_params['res_id']
            interaction_type = line_params['interaction_type']
            
            # Get directionality for adjusting marker position
            direction = self.interaction_direction.get((res_id, interaction_type), 'bidirectional')
            
            # Calculate points along the curved path
            path_points = []
            steps = 20
            for i in range(steps + 1):
                t = i / steps
                # For bidirectional, use the protein->ligand curve for marker
                # For directional, offset slightly based on direction
                curve_adjust = curvature
                if direction == 'ligand_to_protein':
                    curve_adjust = -curvature
                
                control_x = (start_pos[0] + end_pos[0])/2 + curve_adjust * (end_pos[1] - start_pos[1]) * 2
                control_y = (start_pos[1] + end_pos[1])/2 - curve_adjust * (end_pos[0] - start_pos[0]) * 2
                x = (1-t)*(1-t)*start_pos[0] + 2*(1-t)*t*control_x + t*t*end_pos[0]
                y = (1-t)*(1-t)*start_pos[1] + 2*(1-t)*t*control_y + t*t*end_pos[1]
                path_points.append((x, y))
            
            # Find best marker position
            best_position = None
            best_score = float('-inf')
            
            for t in [0.5, 0.45, 0.55, 0.4, 0.6, 0.35, 0.65, 0.3, 0.7, 0.25, 0.75]:
                idx = int(t * steps)
                pos = path_points[idx]
                
                # Calculate distance to existing markers
                if marker_positions:  # Only if there are existing markers
                    min_dist = min(math.sqrt((pos[0]-p[0])**2 + (pos[1]-p[1])**2) 
                                for p in marker_positions.values())
                else:
                    min_dist = float('inf')
                
                text_len = len(style['marker_text'])
                min_req_dist = 25 + text_len * 2
                score = min(min_dist / min_req_dist, 2.0) + (1.0 - abs(t - 0.5))
                
                if score > best_score:
                    best_score = score
                    best_position = pos
            
            if best_position is None:
                best_position = path_points[len(path_points)//2]
            
            marker_positions[key] = best_position
            x, y = best_position
            
            # Draw marker shape
            marker_radius = 9 + (len(style['marker_text']) - 1) * 1.5
            if 'pi' in line_params['interaction_type']:
                angles = np.linspace(0, 2*np.pi, 7)[:-1]
                vertices = [(x + marker_radius * math.cos(a), y + marker_radius * math.sin(a)) 
                        for a in angles]
                marker = Polygon(vertices, closed=True, facecolor=style['marker_bg'],
                            edgecolor=style['color'], linewidth=1.5, zorder=5)
            else:
                marker = Circle((x, y), marker_radius, facecolor=style['marker_bg'],
                            edgecolor=style['color'], linewidth=1.5, zorder=5)
            ax.add_patch(marker)
            
            # Add marker text
            text = ax.text(x, y, style['marker_text'], ha='center', va='center',
                        fontsize=max(7, 9 - (len(style['marker_text']) - 1) * 0.8),
                        color=style['color'], fontweight='bold', zorder=6)
            text.set_path_effects([path_effects.withStroke(linewidth=1, foreground='white')])

        # Create legend
        legend_elements = [
            Rectangle((0, 0), 1, 1, facecolor='white', edgecolor='black',
                    label='Interacting structural groups')
        ]
        
        # Add interaction type markers to legend
        for int_type, style in interaction_styles.items():
            if self.interactions[int_type]:
                marker = 'h' if 'pi' in int_type else 'o'
                legend_elements.append(
                    Line2D([0], [0], color=style['color'], linestyle=style['linestyle'],
                        linewidth=style['linewidth'], marker=marker,
                        markerfacecolor=style['marker_bg'], markeredgecolor=style['color'],
                        markersize=8, label=style['name'])
                )
        
        # Add directionality to legend
        legend_elements.append(
            Line2D([0], [0], color='black', linestyle='-', marker='>',
                markerfacecolor='black', markersize=8, 
                label='Unidirectional interaction')
        )
        # For bidirectional arrow, use a simpler approach with a diamond symbol
        bidirectional = Line2D([0], [0], color='black', linestyle='-',
                            marker='d', markerfacecolor='black',
                            markersize=8, label='Bidirectional interaction')
        legend_elements.append(bidirectional)
        
        # Add solvent accessibility indicator
        if self.solvent_accessible:
            legend_elements.append(
                Rectangle((0, 0), 1, 1, facecolor='#ADD8E6', alpha=0.5,
                        edgecolor='#87CEEB', label='Solvent accessible')
            )
        
        # Draw legend
        legend = ax.legend(
            handles=legend_elements,
            title="Interacting structural groups",
            loc='upper right',
            frameon=True,
            framealpha=0.7,
            fontsize=9,
            title_fontsize=10
        )
        
        # Set plot limits and appearance
        max_coord = radius + 100
        ax.set_xlim(-max_coord, max_coord)
        ax.set_ylim(-max_coord, max_coord)
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Add title
        if title:
            plt.title(title, fontsize=16)
        else:
            base_name = os.path.splitext(os.path.basename(self.structure_file))[0]
            plt.title(f"Protein-Ligand Interactions: {base_name}", fontsize=16)
        
        # Save figure
        plt.tight_layout()
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Interaction diagram saved to {output_file}")
        return output_file
    
    
    def _count_ligand_rotatable_bonds(self):
        """
        Count freely rotatable bonds in the ligand for the flexibility penalty.

        Definition (standard, matching RDKit/Böhm): single bond between two
        non-ring heavy atoms where both atoms have degree ≥ 2 (not terminal).
        """
        bonds = self.ligand_structure.find_bonds(distance_threshold=2.0)
        # Build heavy-atom adjacency
        adj = {}
        for atom_id, info in self.ligand_structure.atom_coords.items():
            if info.get('element', '') not in ('H', ''):
                adj[atom_id] = []
        for a, b in bonds:
            ea = self.ligand_structure.atom_coords.get(a, {}).get('element', '')
            eb = self.ligand_structure.atom_coords.get(b, {}).get('element', '')
            if ea in ('H', '') or eb in ('H', ''):
                continue
            if a in adj and b in adj:
                adj[a].append(b)
                adj[b].append(a)
        ring_atoms = self.ligand_ring_atoms
        rotatable = set()
        for a, b in bonds:
            ea = self.ligand_structure.atom_coords.get(a, {}).get('element', '')
            eb = self.ligand_structure.atom_coords.get(b, {}).get('element', '')
            if ea in ('H', '') or eb in ('H', ''):
                continue
            if a not in adj or b not in adj:
                continue
            if len(adj[a]) < 2 or len(adj[b]) < 2:
                continue  # terminal atom — not rotatable
            if a in ring_atoms or b in ring_atoms:
                continue  # ring bond — not rotatable
            rotatable.add(tuple(sorted([a, b])))
        return len(rotatable)

    def estimate_binding_affinity(self):
        """
        Compute an empirical interaction score on a kcal/mol scale.

        Despite the historical method name this is NOT a predicted binding
        free energy: see ``interaction_score`` and the 'calibrated' flag.

        Three scientifically grounded improvements over a simple count-based sum:

        1. **Per-residue deduplication** — multiple atom-atom contacts to the
           same residue for the same interaction type count as one interaction
           (using the closest contact distance). This prevents PHE contributing
           7× for 7 alkyl atoms within range.

        2. **Distance decay** — each residue contribution is scaled by a linear
           decay from 1.0 at the ideal distance to 0.0 at the cutoff distance,
           following the Böhm (1994) functional form. Closer contacts score
           higher than barely-detected ones.

        3. **Rotatable bond penalty** — +0.5 kcal/mol per freely rotatable
           non-ring ligand bond, reflecting the conformational entropy cost of
           restricting the ligand on binding (Böhm 1994; AutoDock empirical).

        Limitations: weights are not regression-fitted to a PDB benchmark;
        the score does not track measured affinities and carries no Kd label.
        Not a substitute for FEP or MM-GBSA.

        Returns
        -------
        dict with keys:
            'score'          : float — the interaction score
            'dG_estimated'   : float — the same number (older name)
            'breakdown'      : dict  — per-type contribution after deduplication
            'calibrated'     : False
            'interpretation' : str   — what the number is and is not
            'note'           : str   — how to use it
        """
        # Empirical weights (kcal/mol per unique interacting residue)
        # Based on Böhm 1994 J Comput Aided Mol Des and X-Score Wang 2002.
        WEIGHTS = {
            'hydrogen_bonds':      -1.20,
            'covalent':            -8.00,
            'metal_coordination':  -2.50,
            # 'salt_bridge' is an alias of 'ionic' (same detection); it carries
            # no separate weight so charge-charge contacts are counted once.
            'ionic':               -1.80,
            'halogen_bonds':       -0.80,
            'pi_pi_stacking':      -0.70,
            'cation_pi':           -0.65,
            'pi_cation':           -0.65,
            'hydrophobic':         -0.30,
            'carbon_pi':           -0.25,
            'donor_pi':            -0.20,
            'amide_pi':            -0.20,
            'alkyl_pi':            -0.20,
            'attractive_charge':   -0.50,
            'repulsion':           +0.60,
        }

        # Ideal contact distance (full weight) and cutoff (zero weight) per type
        IDEAL_DIST = {
            'hydrogen_bonds': 2.9,  'covalent': 1.8,
            'metal_coordination': 2.2, 'ionic': 3.5, 'salt_bridge': 3.5,
            'halogen_bonds': 3.0,   'pi_pi_stacking': 3.8,
            'cation_pi': 4.0,       'pi_cation': 4.0,
            'hydrophobic': 3.5,     'carbon_pi': 4.0,
            'donor_pi': 3.5,        'amide_pi': 3.5,
            'alkyl_pi': 4.0,        'attractive_charge': 3.5,
            'repulsion': 3.0,
        }
        MAX_DIST = {
            'hydrogen_bonds': 3.5,  'covalent': 2.1,
            'metal_coordination': 2.8, 'ionic': 5.5, 'salt_bridge': 5.5,
            'halogen_bonds': 3.5,   'pi_pi_stacking': 5.5,
            'cation_pi': 5.5,       'pi_cation': 5.5,
            'hydrophobic': 4.0,     'carbon_pi': 5.5,
            'donor_pi': 5.5,        'amide_pi': 5.5,
            'alkyl_pi': 5.5,        'attractive_charge': 5.5,
            'repulsion': 4.0,
        }

        # Translational/rotational entropy + base desolvation constant (Böhm 1994)
        DELTA_G_CONST = +2.9  # kcal/mol

        # Ligand flexibility penalty: +0.5 kcal/mol per rotatable bond
        n_rotors = self._count_ligand_rotatable_bonds()
        rotor_penalty = n_rotors * 0.5

        breakdown = {}
        dG = DELTA_G_CONST + rotor_penalty

        for itype, weight in WEIGHTS.items():
            contacts = self.interactions.get(itype, [])
            if not contacts:
                continue

            # Deduplicate: keep closest contact per unique residue
            best_per_residue = {}
            for c in contacts:
                res = c['protein_residue']
                res_key = (res.resname, res.id[1])
                d = c['distance']
                if res_key not in best_per_residue or d < best_per_residue[res_key]:
                    best_per_residue[res_key] = d

            # Sum distance-decayed contributions across unique residues
            d_ideal = IDEAL_DIST.get(itype, 3.5)
            d_max = MAX_DIST.get(itype, 5.5)
            total = 0.0
            for d in best_per_residue.values():
                if d <= d_ideal:
                    decay = 1.0
                else:
                    decay = max(0.0, (d_max - d) / max(d_max - d_ideal, 1e-6))
                total += weight * decay

            dG += total
            breakdown[itype] = {
                'unique_residues': len(best_per_residue),
                'weight': weight,
                'contribution_kcal_mol': round(total, 2),
            }

        if rotor_penalty > 0:
            breakdown['rotatable_bond_penalty'] = {
                'unique_residues': n_rotors,
                'weight': +0.5,
                'contribution_kcal_mol': round(rotor_penalty, 2),
            }

        # No binding-strength label is attached.  The weights are literature
        # values, not fitted to measured affinities, and the sum does not
        # track them: erlotinib-EGFR (PDB 1M17, a nanomolar complex) scores
        # close to zero.  Earlier releases mapped this number to Kd ranges
        # ("very weak / no binding expected"), which was misleading.
        score = round(dG, 2)
        return {
            'score': score,
            'dG_estimated': score,      # same number; kept for older scripts
            'breakdown': breakdown,
            'calibrated': False,
            'interpretation': SCORE_INTERPRETATION,
            'note': SCORE_NOTE,
        }

    def interaction_score(self):
        """
        Empirical interaction score (kcal/mol scale).  Preferred name for
        ``estimate_binding_affinity``; both return the same dictionary.
        """
        return self.estimate_binding_affinity()

    def filter_interactions_directly(self):
        """Emergency filter to remove chemically implausible interactions."""
        # Define which residues can participate in which interaction types
        aromatic_residues = {'PHE', 'TYR', 'TRP', 'HIS'}
        charged_residues = {'ASP', 'GLU', 'LYS', 'ARG', 'HIS'}
        neg_charged = {'ASP', 'GLU'}
        pos_charged = {'LYS', 'ARG', 'HIS'}
        
        # Filter pi-stacking - require at least one aromatic residue
        self.interactions['pi_pi_stacking'] = [
            i for i in self.interactions['pi_pi_stacking']
            if i['protein_residue'].resname in aromatic_residues and i['distance'] < 5.5
        ]
        
        # Filter ionic/salt bridge - require charged residues
        self.interactions['ionic'] = [
            i for i in self.interactions['ionic']
            if i['protein_residue'].resname in charged_residues and i['distance'] < 4.0
        ]
        self.interactions['salt_bridge'] = [
            i for i in self.interactions['salt_bridge']
            if i['protein_residue'].resname in charged_residues and i['distance'] < 4.0
        ]
        
        # Filter covalent bonds - require very close distance
        self.interactions['covalent'] = [
            i for i in self.interactions['covalent']
            if i['distance'] < 2.1
        ]
        
        # Filter amide-pi - require amide residues
        self.interactions['amide_pi'] = [
            i for i in self.interactions['amide_pi']
            if i['protein_residue'].resname in {'ASN', 'GLN'} and i['distance'] < 5.5
        ]
        
        # Filter donor-pi - require negatively charged residues
        self.interactions['donor_pi'] = [
            i for i in self.interactions['donor_pi']
            if i['protein_residue'].resname in neg_charged and i['distance'] < 5.5
        ]
    def interactions_as_dict(self):
        """
        The analysis as plain Python data (JSON-serialisable): structure and
        ligand, a summary, every interaction with the atoms involved, and a
        per-residue roll-up.  ``export_interactions_json`` writes exactly this.
        """
        lig = self.ligand_residue
        contacts = []
        for itype, items in self.interactions.items():
            for c in items:
                res, pa, la = c.get('protein_residue'), c.get('protein_atom'), c.get('ligand_atom')
                if res is None or pa is None or la is None:
                    continue
                parent = res.get_parent()
                entry = {
                    'type': itype,
                    'residue': {'name': res.resname.strip(), 'number': int(res.id[1]),
                                'chain': parent.id if parent is not None else ''},
                    'protein_atom': pa.get_name().strip(),
                    'ligand_atom': la.get_name().strip(),
                    'ligand_element': (la.element or '').strip().capitalize(),
                    'distance': round(float(c['distance']), 3),
                    'solvent_accessible': (res.resname, res.id[1]) in self.solvent_accessible,
                }
                if c.get('angle') is not None:
                    entry['angle'] = round(float(c['angle']), 1)
                if c.get('type'):
                    entry['note'] = str(c['type']).replace('_', ' ')
                if c.get('halogen_element'):
                    entry['halogen_element'] = c['halogen_element']
                contacts.append(entry)
        contacts.sort(key=lambda e: (e['type'], e['distance']))

        by_type, residues = {}, {}
        for e in contacts:
            by_type[e['type']] = by_type.get(e['type'], 0) + 1
            r = e['residue']
            row = residues.setdefault((r['chain'], r['number'], r['name']), {
                **r, 'solvent_accessible': e['solvent_accessible'],
                'n_contacts': 0, 'closest_distance': e['distance'], 'types': []})
            row['n_contacts'] += 1
            row['closest_distance'] = min(row['closest_distance'], e['distance'])
            if e['type'] not in row['types']:
                row['types'].append(e['type'])

        try:
            from importlib.metadata import version as _version
            pandamap_version = _version('pandamap')
        except Exception:
            pandamap_version = None
        return {
            'pandamap_version': pandamap_version,
            'structure': os.path.basename(str(self.structure_file)),
            'ligand': {
                'name': lig.resname.strip(),
                'chain': lig.get_parent().id,
                'number': int(lig.id[1]),
                'n_heavy_atoms': sum(1 for a in self.ligand_atoms
                                     if (a.element or '').strip().upper() != 'H'),
            },
            'units': {'distance': 'angstrom', 'angle': 'degree'},
            'summary': {'n_interactions': len(contacts), 'n_residues': len(residues),
                        'by_type': by_type},
            'interactions': contacts,
            'residues': sorted(residues.values(),
                               key=lambda r: (-r['n_contacts'], r['closest_distance'])),
        }

    def export_interactions_json(self, json_file='interactions.json', include_score=False,
                                 include_affinity=False):
        """
        Write the analysis as JSON (see ``interactions_as_dict`` for the
        layout).  With ``include_score`` the empirical interaction score is
        added under ``interaction_score`` (``include_affinity`` is the older
        name for the same switch).  Returns the path written.
        """
        import json as _json

        data = self.interactions_as_dict()
        if include_score or include_affinity:
            try:
                data['interaction_score'] = _json.loads(
                    _json.dumps(self.estimate_binding_affinity(), default=float))
            except Exception as exc:
                data['interaction_score'] = {'error': str(exc)}
        with open(json_file, 'w', encoding='utf-8') as fh:
            _json.dump(data, fh, indent=2, ensure_ascii=False)
            fh.write('\n')
        logger.info(f"Interaction JSON written to {json_file} "
              f"({data['summary']['n_interactions']} interactions)")
        return json_file

    def export_interactions_csv(self, csv_file='interactions.csv'):
        """
        Write every detected interaction as one row of a tidy CSV.

        The text report is human-readable but not parseable; this gives a
        machine-readable table for downstream filtering, plotting, or
        comparison across structures.  Columns:

            interaction_type, protein_residue, protein_resnum, protein_chain,
            protein_atom, ligand_atom, ligand_element, distance_A,
            solvent_accessible, halogen_element

        Returns the path written.
        """
        import csv as _csv

        rows = []
        for itype, contacts in self.interactions.items():
            for c in contacts:
                res = c.get('protein_residue')
                pa = c.get('protein_atom')
                la = c.get('ligand_atom')
                if res is None or pa is None or la is None:
                    continue
                chain = res.get_parent().id if res.get_parent() is not None else ''
                rows.append({
                    'interaction_type': itype,
                    'protein_residue': res.resname,
                    'protein_resnum': res.id[1],
                    'protein_chain': chain,
                    'protein_atom': pa.get_name(),
                    'ligand_atom': la.get_name(),
                    'ligand_element': la.element.strip().capitalize(),
                    'distance_A': round(float(c['distance']), 2),
                    'solvent_accessible': (res.resname, res.id[1]) in self.solvent_accessible,
                    'halogen_element': c.get('halogen_element', ''),
                })

        rows.sort(key=lambda r: (r['interaction_type'], r['distance_A']))
        fields = ['interaction_type', 'protein_residue', 'protein_resnum',
                  'protein_chain', 'protein_atom', 'ligand_atom',
                  'ligand_element', 'distance_A', 'solvent_accessible',
                  'halogen_element']
        with open(csv_file, 'w', newline='', encoding='utf-8') as fh:
            writer = _csv.DictWriter(fh, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        logger.info(f"Interaction CSV written to {csv_file} ({len(rows)} rows)")
        return csv_file

    def render_contacts_3d(self, output_file='contacts_3d.png', **kwargs):
        """
        Static 3D render of the binding site: solid ligand, transparent
        contacting residues, and each contacting atom highlighted, labelled
        and joined to the ligand atom it touches.  See
        :func:`pandamap.render_3d.render_contacts_3d` for the options.
        """
        from .render_3d import render_contacts_3d
        return render_contacts_3d(self, output_file=output_file, **kwargs)

    def generate_interaction_plots(self, output_file='interaction_report.png', dpi=300):
        """
        Draw the graphical interaction report: one publication-style figure
        with five panels that explain the binding site at a glance.

          a  the 2D interaction map
          b  contacts by interaction type
          c  every contact distance against the typical distance for its type
          d  residue x interaction-type matrix, with solvent exposure
          e  per-type contribution to the empirical dG estimate

        The extension of ``output_file`` selects the format (.png, .pdf,
        .svg).  Returns the path written, or None when nothing was detected.
        """
        from .report_figure import draw_report_figure
        return draw_report_figure(self, output_file=output_file, dpi=dpi)

    def run_analysis(self, output_file=None, use_dssp=True, generate_report=False,
                     report_file=None, dpi=300, title=None, show_3d_cues=True,
                     diagrams=DIAGRAM_STYLES, show_nearby_residues=False):
        """
        Run the complete analysis pipeline.
        
        Parameters:
        -----------
        output_file : str, optional
            Path where the output image will be saved.
        use_dssp : bool
            Whether to use DSSP for solvent accessibility (default: True)
        diagrams : sequence of str
            Which diagrams to draw: any of 'modern', 'classic', '3d'
            (default: all three); an empty tuple draws nothing.  The first is
            written to ``output_file``,
            the others next to it with a ``_classic`` / ``_3d`` suffix.  All
            paths are stored in ``self.diagram_files``.
        generate_report : bool
            Whether to generate a PLIP-like text report (default: False)
        report_file : str, optional
            Path where the report will be saved. If None but generate_report is True,
            a default path will be used.
        """
        if output_file is None:
            base_name = os.path.splitext(os.path.basename(self.structure_file))[0]
            output_file = f"{base_name}_interactions.png"
        
        # Detect protein-ligand interactions
        logger.info("Detecting interactions...")
        self.detect_interactions() 
        
        logger.info(f"Before filtering: {sum(len(ints) for ints in self.interactions.values())} total interactions")
        
        # Calculate solvent accessibility
        logger.info("Calculating solvent accessibility...")
        try:
            import shutil as _shutil
            dssp_exe = _shutil.which('mkdssp') or _shutil.which('dssp')
            if use_dssp and dssp_exe is None:
                # DSSP is optional; the built-in estimate needs no extra tools
                self.calculate_realistic_solvent_accessibility()
            elif use_dssp:
                # Try DSSP first
                logger.info("Trying DSSP method first...")
                self.calculate_dssp_solvent_accessibility(dssp_executable=dssp_exe)
                
                # Check if DSSP found too many solvent accessible residues
                if len(self.solvent_accessible) > len(self.interacting_residues) * 0.5:
                    logger.info("DSSP found too many solvent accessible residues, using realistic method")
                    self.calculate_realistic_solvent_accessibility()
                # Check if DSSP found too few or no solvent accessible residues
                elif len(self.solvent_accessible) < 2:
                    logger.info("DSSP didn't find enough solvent accessible residues, using realistic method")
                    self.calculate_realistic_solvent_accessibility()
            else:
                # Use our realistic method
                self.calculate_realistic_solvent_accessibility()
        except Exception as e:
            logger.warning(f"Error calculating solvent accessibility: {str(e)}")
            logger.info("Falling back to realistic method")
            self.calculate_realistic_solvent_accessibility(exposure_threshold=0.2)
        
        # Make sure we don't have all residues marked as solvent accessible
        if len(self.solvent_accessible) == len(self.interacting_residues):
            logger.warning("WARNING: All residues marked as solvent accessible - this is likely incorrect")
            logger.info("Applying stricter filtering")
            
            # Keep only 40% max of interacting residues as solvent accessible
            self.calculate_realistic_solvent_accessibility(exposure_threshold=0.3, max_percent=0.4)
        
        # Generate visualization
        logger.info("Generating visualization...")
        # --dpi / --title / --no-3d-cues were accepted by the CLI but never
        # reached visualize(), so they silently did nothing.
        # All requested diagrams.  The first one is written to ``output_file``;
        # the others get a suffix (_classic, _3d) before the extension.
        if isinstance(diagrams, str):
            diagrams = (diagrams,)
        # an empty selection (diagrams=()) runs the analysis without drawing
        diagrams = list(dict.fromkeys(diagrams or ()))
        stem, ext = os.path.splitext(output_file)
        self.diagram_files = {}
        for style in diagrams:
            path = output_file if style == diagrams[0] else \
                f"{stem}{_DIAGRAM_SUFFIX.get(style, '_' + style)}{ext}"
            self.diagram_files[style] = self.visualize(
                output_file=path, dpi=dpi, title=title, show_3d_cues=show_3d_cues,
                style=style, show_nearby_residues=show_nearby_residues)
        viz_file = self.diagram_files[diagrams[0]] if diagrams else None
        
        # Generate text report if requested
        if generate_report:
            if report_file is None:
                base_name = os.path.splitext(os.path.basename(self.structure_file))[0]
                report_file = f"{base_name}_interactions_report.txt"
            
            logger.info("Generating interaction report...")
            try:
                if hasattr(self, 'generate_interaction_report'):
                    self.generate_interaction_report(output_file=report_file)
                else:
                    logger.warning("Warning: Report generation not available, skipping.")
            except Exception as e:
                logger.warning(f"Error generating report: {str(e)}")
        
        return viz_file


# Attach the report methods to the mapper class
HybridProtLigMapper = add_interaction_detection_methods(HybridProtLigMapper)

from .trajectory import analyze_trajectory                                       # noqa: E402,F401

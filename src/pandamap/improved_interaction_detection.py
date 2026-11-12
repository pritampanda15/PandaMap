#!/usr/bin/env python
"""
Improved interaction detection for PandaMap.
Includes PLIP-like thresholds and filtering plus text report generation.
Add this to your core.py module or create a new file and import it.
"""

import logging
import os
import math
import numpy as np
from collections import defaultdict, namedtuple

logger = logging.getLogger("pandamap")

# Constants for interaction detection (PLIP-like thresholds)
DISTANCE_CUTOFFS = {
    'hydrogen_bond': 3.5,       # Max distance for hydrogen bonds
    'hydrophobic': 4.0,         # Max distance for hydrophobic interactions
    'pi_stacking': 5.5,         # Max distance for pi-stacking
    'pi_cation': 6.0,           # Max distance for pi-cation interactions
    'salt_bridge': 5.5,         # Max distance for salt bridges
    'halogen': 3.5,             # Max distance for halogen bonds
    'water_bridge': 3.5,        # Max distance for water bridges
    'min_dist': 1.5             # Minimum distance for any interaction
}

ANGLE_CUTOFFS = {
    'pi_stack_parallel': 30.0,      # Max deviation from parallel planes (degrees)
    'pi_stack_perpendicular': 30.0, # Max deviation from perpendicular planes (degrees)
    'pi_cation': 30.0,              # Max deviation for pi-cation (degrees)
    'hbond_donor': 120.0,           # Min angle for hydrogen bond donor (degrees)
    'halogen_donor': 140.0,         # Optimal angle for halogen donor (C-X...Y)
    'halogen_acceptor': 120.0,      # Optimal angle for halogen acceptor (X...Y-C)
    'halogen_angle_dev': 30.0       # Max deviation from optimal halogen bond angles
}

OFFSET_CUTOFFS = {
    'pi_stack': 2.0,    # Max offset for pi-stacking
    'pi_cation': 2.0    # Max offset for pi-cation
}

def filter_interactions(interactions, thresholds):
    """Apply stricter filtering to detected interactions."""
    filtered = {}
    
    # Filter hydrogen bonds
    if 'hydrogen_bonds' in interactions:
        filtered['hydrogen_bonds'] = [
            hb for hb in interactions['hydrogen_bonds'] 
            #if (thresholds['min_dist'] < hb['distance'] < thresholds['hydrogen_bond'] and
                #('angle' not in hb or hb['angle'] > ANGLE_CUTOFFS['hbond_donor']))
        ]
        
    # Filter hydrophobic interactions
    if 'hydrophobic' in interactions:
        filtered['hydrophobic'] = [
            h for h in interactions['hydrophobic']
            if thresholds['min_dist'] < h['distance'] < thresholds['hydrophobic']
        ]
        
    # Filter pi-stacking interactions
    if 'pi_pi_stacking' in interactions:
        filtered['pi_pi_stacking'] = [
            pi for pi in interactions['pi_pi_stacking']
            if (thresholds['min_dist'] < pi['distance'] < thresholds['pi_stacking'] and
                (
                    # Parallel stacking (angle between ring normals < 30°)
                    (pi['angle'] < ANGLE_CUTOFFS['pi_stack_parallel'] and 
                     pi['offset'] < OFFSET_CUTOFFS['pi_stack']) or
                    # T-shaped stacking (angle between ring normals ~ 90°)
                    (abs(pi['angle'] - 90) < ANGLE_CUTOFFS['pi_stack_perpendicular'] and 
                     pi['offset'] < OFFSET_CUTOFFS['pi_stack'])
                )
            )
        ]
        
    # Filter pi-cation interactions
    if 'carbon_pi' in interactions:
        filtered['carbon_pi'] = [
            pic for pic in interactions['carbon_pi']
            if (thresholds['min_dist'] < pic['distance'] < thresholds['pi_cation'] and
                pic['offset'] < OFFSET_CUTOFFS['pi_cation'])
        ]
        
    # Filter ionic interactions (salt bridges)
    if 'donor_pi' in interactions:
        filtered['donor_pi'] = [
            sb for sb in interactions['donor_pi']
            if thresholds['min_dist'] < sb['distance'] < thresholds['salt_bridge']
        ]

    # Filter amide-pi interactions
    if 'amide_pi' in interactions:
        filtered['amide_pi'] = [
            ap for ap in interactions['amide_pi']
            if thresholds['min_dist'] < ap['distance'] < thresholds['pi_cation']
        ]
        
    # Filter halogen bonds
    if 'halogen_bonds' in interactions:
        filtered['halogen_bonds'] = [
            hal for hal in interactions['halogen_bonds']
            if (thresholds['min_dist'] < hal['distance'] < thresholds['halogen'] and
                abs(hal['don_angle'] - ANGLE_CUTOFFS['halogen_donor']) < ANGLE_CUTOFFS['halogen_angle_dev'] and
                abs(hal['acc_angle'] - ANGLE_CUTOFFS['halogen_acceptor']) < ANGLE_CUTOFFS['halogen_angle_dev'])
        ]
        
    # Filter water bridges
    if 'water_bridges' in interactions:
        filtered['water_bridges'] = [
            wb for wb in interactions['water_bridges']
            if (thresholds['min_dist'] < wb['distance_aw'] < thresholds['water_bridge'] and
                thresholds['min_dist'] < wb['distance_dw'] < thresholds['water_bridge'])
        ]
    
    # Filter alkyl-pi interactions
    if 'alkyl_pi' in interactions:
        filtered['alkyl_pi'] = [
            ap for ap in interactions['alkyl_pi']
            if thresholds['min_dist'] < ap['distance'] < thresholds['pi_stacking']
        ]
    
    # Filter pi-cation interactions
    if 'pi_cation' in interactions:
        filtered['pi_cation'] = [
            pc for pc in interactions['pi_cation']
            if thresholds['min_dist'] < pc['distance'] < thresholds['pi_cation'] and
            pc['offset'] < OFFSET_CUTOFFS['pi_cation']
        ]
    
    # Filter attractive charge interactions
    if 'attractive_charge' in interactions:
        filtered['attractive_charge'] = [
            ac for ac in interactions['attractive_charge']
            if thresholds['min_dist'] < ac['distance'] < thresholds['salt_bridge']
        ]
    
    # Filter repulsion interactions
    if 'repulsion' in interactions:
        filtered['repulsion'] = [
            r for r in interactions['repulsion']
            if thresholds['min_dist'] < r['distance'] < thresholds['salt_bridge'] * 1.5
        ]
    
    return filtered

def refine_hydrophobic_interactions(interactions):
    """Improve hydrophobic interaction filtering to avoid excessive detections.
    Similar to PLIP's approach of clustering and keeping only the strongest interactions.
    """
    if not interactions:
        return []
    
    # Group by protein residue
    by_residue = defaultdict(list)
    for interaction in interactions:
        res_id = (interaction['restype'], interaction['resnr'], interaction['reschain'])
        by_residue[res_id].append(interaction)
    
    # Keep only the closest interaction for each residue
    refined = []
    for res_interactions in by_residue.values():
        closest = min(res_interactions, key=lambda x: x['distance'])
        refined.append(closest)
    
    # Group by ligand atom to find clusters
    by_ligand_atom = defaultdict(list)
    for interaction in refined:
        lig_atom_id = interaction['ligatom_orig_idx']
        by_ligand_atom[lig_atom_id].append(interaction)
    
    # For each ligand atom, keep only the closest interaction
    final_refined = []
    for lig_interactions in by_ligand_atom.values():
        closest = min(lig_interactions, key=lambda x: x['distance'])
        final_refined.append(closest)
    
    return final_refined

def normalize_interactions(interactions):
    """
    Flatten interaction records into the keys the report formatter expects.

    ``HybridProtLigMapper.interactions`` stores the protein residue as a
    BioPython ``Residue`` object under ``'protein_residue'``.  The report
    formatter reads the flat keys ``restype`` / ``resnr`` / ``reschain``.
    When a caller passed ``mapper.interactions`` straight to
    :func:`generate_plip_report` — exactly what the documented Python API
    shows — none of those keys were present and every line rendered as
    ``UNK?? -- 2.79A -- LIG``.  The CLI looked correct only because it
    reached the same formatter through a different code path that did this
    conversion first.

    Records that already carry the flat keys are passed through unchanged,
    so this is safe to apply on either path.

    Parameters
    ----------
    interactions : dict
        Mapping of interaction type -> list of interaction records.

    Returns
    -------
    dict
        Same structure, with ``restype``, ``resnr`` and ``reschain``
        populated wherever they could be derived.
    """
    if not isinstance(interactions, dict):
        return interactions

    normalized = {}
    for itype, records in interactions.items():
        out = []
        for record in records or []:
            if not isinstance(record, dict):
                out.append(record)
                continue

            flat = dict(record)
            residue = record.get('protein_residue')

            if flat.get('restype') in (None, 'UNK') and hasattr(residue, 'resname'):
                flat['restype'] = residue.resname

            if flat.get('resnr') in (None, 0, '?'):
                res_id = getattr(residue, 'id', None)
                if isinstance(res_id, tuple) and len(res_id) > 1:
                    flat['resnr'] = res_id[1]

            if flat.get('reschain') in (None, 'X', '?'):
                parent = None
                if hasattr(residue, 'get_parent') and callable(residue.get_parent):
                    parent = residue.get_parent()
                elif hasattr(residue, 'parent'):
                    parent = residue.parent
                if parent is not None and getattr(parent, 'id', None) is not None:
                    flat['reschain'] = parent.id

            # Atom names, so the report can say which atoms are in contact
            for key, name_key in (('protein_atom', 'protein_atom_name'),
                                  ('ligand_atom', 'ligand_atom_name')):
                atom = record.get(key)
                if not flat.get(name_key) and hasattr(atom, 'get_name'):
                    flat[name_key] = atom.get_name().strip()

            out.append(flat)
        normalized[itype] = out
    return normalized


# interaction type -> (section title, short name used in the residue table)
REPORT_SECTIONS = {
    'hydrogen_bonds': ("Hydrogen bonds", "hydrogen bond"),
    'salt_bridge': ("Salt bridges", "salt bridge"),
    'ionic': ("Ionic interactions", "ionic"),
    'attractive_charge': ("Attractive charges", "attractive charge"),
    'halogen_bonds': ("Halogen bonds", "halogen bond"),
    'metal_coordination': ("Metal coordination", "metal"),
    'covalent': ("Covalent bonds", "covalent"),
    'pi_pi_stacking': ("π-π stacking", "π-π stacking"),
    'pi_cation': ("π-cation interactions", "π-cation"),
    'cation_pi': ("Cation-π interactions", "cation-π"),
    'donor_pi': ("Donor-π interactions", "donor-π"),
    'amide_pi': ("Amide-π interactions", "amide-π"),
    'carbon_pi': ("Carbon-π interactions", "carbon-π"),
    'alkyl_pi': ("Alkyl-π interactions", "alkyl-π"),
    'hydrophobic': ("Hydrophobic contacts", "hydrophobic"),
    'water_bridges': ("Water bridges", "water bridge"),
    'repulsion': ("Repulsion", "repulsion"),
}
_REPORT_WIDTH = 78


def generate_plip_report(ligand_info, interactions, output_file=None):
    """
    Write the plain-text interaction report.

    Layout: a header with the structure and ligand, a summary with one bar
    per interaction type, one aligned table per interaction type (residue,
    residue atom, ligand atom, distance - closest first) and a per-residue
    table.  ``ligand_info`` may carry the optional keys ``structure``,
    ``heavy_atoms`` and ``solvent_accessible`` (a set of (resname, resnum)).
    """
    # Accept raw HybridProtLigMapper.interactions as well as pre-flattened
    # records, so the Python API and the CLI produce identical reports.
    interactions = normalize_interactions(interactions)
    if not isinstance(interactions, dict):
        interactions = {}

    heavy_rule, light_rule = "=" * _REPORT_WIDTH, "-" * _REPORT_WIDTH
    try:
        from importlib.metadata import version as _version
        stamp = f"PandaMap {_version('pandamap')}"
    except Exception:
        stamp = "PandaMap"

    def records(itype):
        return [r for r in interactions.get(itype) or [] if isinstance(r, dict)]

    def residue_of(rec):
        return (f"{str(rec.get('restype', 'UNK')).strip()} {rec.get('resnr', '?')} "
                f"{rec.get('reschain', '?')}")

    order = [t for t in REPORT_SECTIONS if records(t)]
    order += [t for t in interactions if t not in REPORT_SECTIONS and records(t)]
    title_of = lambda t: REPORT_SECTIONS.get(t, (t.replace('_', ' ').capitalize(),))[0]
    short_of = lambda t: REPORT_SECTIONS.get(t, (None, t.replace('_', ' ')))[1]

    by_residue, chains = {}, set()
    for itype in order:
        for rec in records(itype):
            by_residue.setdefault(residue_of(rec), []).append((itype, rec))
            if rec.get('reschain') not in (None, '?', 'X'):
                chains.add(str(rec['reschain']))
    chains = sorted(chains or set(ligand_info.get('interacting_chains') or []))
    total = sum(len(records(t)) for t in order)
    solvent = {(str(n).strip(), num) for n, num in ligand_info.get('solvent_accessible') or []}

    hetid = str(ligand_info.get('hetid', 'UNK')).strip()
    ligand = f"{hetid}  (chain {ligand_info.get('chain', '?')}, residue {ligand_info.get('position', '?')}"
    if ligand_info.get('heavy_atoms'):
        ligand += f", {ligand_info['heavy_atoms']} heavy atoms"
    ligand += ")"

    report = [heavy_rule,
              f"  {'PandaMap Interaction Report':<{_REPORT_WIDTH - 4 - len(stamp)}}{stamp}",
              heavy_rule, ""]
    if ligand_info.get('structure'):
        report.append(f"  {'Structure':<15}{ligand_info['structure']}")
    report.append(f"  {'Ligand':<15}{ligand}")
    for field, label in [('molweight', 'Mol. weight'), ('num_rings', 'Rings'),
                         ('num_rot_bonds', 'Rot. bonds'), ('num_hba', 'H-bond acc.'),
                         ('num_hbd', 'H-bond don.'), ('logp', 'logP')]:
        if field in ligand_info:
            report.append(f"  {label:<15}{ligand_info[field]}")
    site = f"{len(by_residue)} residue{'s' if len(by_residue) != 1 else ''}"
    if chains:
        site += f" in chain{'s' if len(chains) > 1 else ''} {', '.join(chains)}"
    if solvent:
        exposed = sum(1 for key in by_residue
                      if (key.split()[0], _as_int(key.split()[1])) in solvent)
        site += f"; {exposed} solvent accessible"
    report.append(f"  {'Binding site':<15}{site}")
    report.append(f"  {'Interactions':<15}{total} contact{'s' if total != 1 else ''} "
                  f"of {len(order)} type{'s' if len(order) != 1 else ''}")

    # ---- summary -------------------------------------------------------------
    report += ["", "SUMMARY", light_rule]
    if not order:
        report.append("  No interactions detected.")
    widest = max((len(title_of(t)) for t in order), default=0)
    most = max((len(records(t)) for t in order), default=1)
    for itype in order:
        n = len(records(itype))
        bar = "█" * max(1, round(n / most * 30))
        report.append(f"  {title_of(itype):<{widest}}  {n:>3}  {bar}")

    # ---- one table per interaction type --------------------------------------
    for itype in order:
        rows = sorted(records(itype), key=lambda r: r.get('distance', 0.0))
        report += ["", f"{title_of(itype).upper().replace('Π', 'π')} ({len(rows)})", light_rule,
                   f"  {'#':>2}  {'Residue':<11}{'Atom':<7}{'Ligand atom':<14}{'Distance':>8}  Note"]
        for n, rec in enumerate(rows, 1):
            try:
                note = ""
                if rec.get('type'):
                    note = str(rec['type']).replace('_', ' ')
                elif itype == 'pi_pi_stacking' and rec.get('angle') is not None:
                    note = "parallel" if float(rec['angle']) < 30 else "T-shaped"
                if rec.get('halogen_element'):
                    note = (note + " " if note else "") + str(rec['halogen_element'])
                line = (f"  {n:>2}  {residue_of(rec):<11}"
                        f"{rec.get('protein_atom_name') or '-':<7}"
                        f"{hetid + ':' + rec['ligand_atom_name'] if rec.get('ligand_atom_name') else hetid:<14}"
                        f"{float(rec.get('distance', 0.0)):>6.2f} Å  {note}")
                report.append(line.rstrip())
            except Exception as e:
                report.append(f"  {n:>2}  (could not format this contact: {e})")

    # ---- per residue ----------------------------------------------------------
    if by_residue:
        report += ["", f"RESIDUES ({len(by_residue)})", light_rule,
                   f"  {'Residue':<11}{'Contacts':>8}  {'Closest':>8}  {'Exposed':<8} Types"]
        ranked = sorted(by_residue.items(),
                        key=lambda kv: (-len(kv[1]), min(r.get('distance', 0.0) for _, r in kv[1])))
        for key, items in ranked:
            kinds = ", ".join(dict.fromkeys(short_of(t) for t, _ in items))
            closest = min(float(r.get('distance', 0.0)) for _, r in items)
            exposed = "yes" if (key.split()[0], _as_int(key.split()[1])) in solvent else "-"
            report.append(f"  {key:<11}{len(items):>8}  {closest:>6.2f} Å  {exposed:<8} {kinds}")

    report += ["", heavy_rule,
               "  Distances are between the named atoms, in Å. Residues are written as",
               "  NAME NUMBER CHAIN. Exposed = solvent accessible. For a table you can",
               "  load into other tools, use --csv.",
               heavy_rule, ""]
    text = "\n".join(report)

    # Write to file if specified
    if output_file:
        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(text)
            logger.info(f"Report saved to {output_file}")
        except Exception as e:
            logger.warning(f"Error writing report to file: {str(e)}")

    return text


def _as_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return value


class ImprovedInteractionDetection:
    """Class to improve the interaction detection in PandaMap.
    Use this to process the interactions after they've been detected by the original methods.
    """
    
    def __init__(self, use_plip_thresholds=True):
        self.use_plip_thresholds = use_plip_thresholds
        self.thresholds = DISTANCE_CUTOFFS if use_plip_thresholds else {
            'hydrogen_bond': 3.5,
            'hydrophobic': 4.0,
            'pi_stacking': 7.0,  # Default PandaMap value
            'pi_cation': 7.0,    # Default PandaMap value
            'salt_bridge': 7.0,  # Default PandaMap value
            'halogen': 4.0,      # Default PandaMap value
            'water_bridge': 4.0, # Default PandaMap value
            'min_dist': 1.5
        }
    
    def refine_interactions(self, interactions):
        """Apply PLIP-like filtering to interactions."""
        # Convert interactions to dicts for easier processing
        dict_interactions = {}
        for itype, ilist in interactions.items():
            dict_interactions[itype] = [self._interaction_to_dict(interaction) for interaction in ilist]
        
        # Apply refined filtering
        filtered_interactions = filter_interactions(dict_interactions, self.thresholds)
        
        # Additional refinement for hydrophobic interactions
        if 'hydrophobic' in filtered_interactions:
            filtered_interactions['hydrophobic'] = refine_hydrophobic_interactions(filtered_interactions['hydrophobic'])
        
        return filtered_interactions
    
    def _interaction_to_dict(self, interaction):
        """Convert a namedtuple interaction to a dictionary."""
        # Convert any namedtuple or object to a dictionary
        if hasattr(interaction, '_asdict'):
            # It's a namedtuple
            return interaction._asdict()
        elif hasattr(interaction, '__dict__'):
            # It's an object with attributes
            return vars(interaction)
        else:
            # Direct dictionary conversion might be needed for some types
            return dict(interaction)
    
    def generate_report(self, ligand_info, interactions, output_file=None):
        """Generate a PLIP-like report for the interactions."""
        # Convert ligand info to a dictionary if it's not already
        if not isinstance(ligand_info, dict):
            ligand_info = vars(ligand_info) if hasattr(ligand_info, '__dict__') else {'hetid': 'UNK', 'chain': 'X', 'position': 0}
        
        return generate_plip_report(ligand_info, interactions, output_file)

# Integration with HybridProtLigMapper
def add_interaction_detection_methods(cls):
    """Add the improved interaction detection methods to the HybridProtLigMapper class."""
    
    def detect_interactions_improved(self, h_bond_cutoff=3.5, pi_stack_cutoff=5.5, hydrophobic_cutoff=4.0):
        """Enhanced version of detect_interactions with better filtering."""
        # First, use the original method to get all potential interactions
        self.detect_interactions(h_bond_cutoff, pi_stack_cutoff, hydrophobic_cutoff)
        
        # Convert namedtuples to dictionaries for easier processing
        interactions = {}
        for itype, ilist in self.interactions.items():
            interactions[itype] = [i._asdict() if hasattr(i, '_asdict') else vars(i) for i in ilist]
        
        # Apply improved filtering
        detector = ImprovedInteractionDetection()
        filtered = detector.refine_interactions(interactions)
        
        # Update the interactions
        for itype, ilist in filtered.items():
            if itype in self.interactions:
                # Convert back to namedtuples if needed
                if ilist and hasattr(self.interactions[itype][0], '_fields'):
                    fields = self.interactions[itype][0]._fields
                    namedtuple_type = type(self.interactions[itype][0])
                    self.interactions[itype] = [namedtuple_type(**item) for item in ilist]
                else:
                    self.interactions[itype] = ilist
        
        return self.interactions
    
    def generate_interaction_report(self, output_file=None):
        """Generate a PLIP-like text report for the interactions."""
        aromatic = {'PHE', 'TYR', 'TRP', 'HIS'}
        self.interactions['pi_pi_stacking'] = [i for i in self.interactions['pi_pi_stacking'] 
                                         if i['protein_residue'].resname in aromatic]
    
        # Remove salt bridges for non-charged residues
        charged = {'ASP', 'GLU', 'LYS', 'ARG', 'HIS'}
        self.interactions['salt_bridge'] = [i for i in self.interactions['salt_bridge'] 
                                      if i['protein_residue'].resname in charged]
        self.interactions['ionic'] = [i for i in self.interactions['ionic'] 
                                if i['protein_residue'].resname in charged]
    
        # Remove covalent bonds with distances > 2.1Å
        self.interactions['covalent'] = [i for i in self.interactions['covalent'] 
                                   if i['distance'] < 2.1]
        
        # Prepare ligand info
        if hasattr(self, 'ligand_residue') and self.ligand_residue:
            ligand_info = {
                'hetid': getattr(self.ligand_residue, 'resname', 'UNK'),
                'chain': 'X',
                'position': 0,
                'longname': getattr(self.ligand_residue, 'resname', 'Unknown'),
                'type': 'LIGAND',
                'interacting_chains': [],
                'interacting_res': []
            }
            
            # Try to get chain and position
            try:
                if hasattr(self.ligand_residue, 'parent') and self.ligand_residue.parent:
                    ligand_info['chain'] = self.ligand_residue.parent.id
                elif hasattr(self.ligand_residue, 'get_parent') and callable(self.ligand_residue.get_parent):
                    parent = self.ligand_residue.get_parent()
                    if parent and hasattr(parent, 'id'):
                        ligand_info['chain'] = parent.id
            except Exception as e:
                logger.warning(f"Error getting ligand chain: {str(e)}")
                
            try:
                if hasattr(self.ligand_residue, 'id') and isinstance(self.ligand_residue.id, tuple):
                    ligand_info['position'] = self.ligand_residue.id[1]
            except Exception as e:
                logger.warning(f"Error getting ligand position: {str(e)}")
        else:
            ligand_info = {'hetid': 'UNK', 'chain': 'X', 'position': 0, 
                        'longname': 'Unknown Ligand', 'type': 'LIGAND',
                        'interacting_chains': [], 'interacting_res': []}
        
        # Process interactions - tailored for dictionary-style interactions
        processed_interactions = {}
        interacting_chains = set()
        interacting_res = set()
        
        for itype, interactions_list in self.interactions.items():
            processed_interactions[itype] = []
            
            for interaction in interactions_list:
                try:
                    # Extract basic interaction data
                    interaction_data = {
                        'distance': 0.0,
                        'restype': 'UNK',
                        'resnr': 0,
                        'reschain': 'X'
                    }
                    
                    # Handle dictionary-style interactions
                    if isinstance(interaction, dict):
                        # Get the distance
                        if 'distance' in interaction:
                            interaction_data['distance'] = interaction['distance']
                        
                        # Extract protein residue info
                        if 'protein_residue' in interaction:
                            res = interaction['protein_residue']
                            if hasattr(res, 'resname'):
                                interaction_data['restype'] = res.resname
                            if hasattr(res, 'id') and isinstance(res.id, tuple) and len(res.id) > 1:
                                interaction_data['resnr'] = res.id[1]
                            if hasattr(res, 'get_parent') and callable(res.get_parent):
                                parent = res.get_parent()
                                if parent and hasattr(parent, 'id'):
                                    interaction_data['reschain'] = parent.id
                                    interacting_chains.add(parent.id)
                        
                        # Direct dictionary access
                        if interaction_data['restype'] == 'UNK' and 'restype' in interaction:
                            interaction_data['restype'] = interaction['restype']
                        if interaction_data['resnr'] == 0 and 'resnr' in interaction:
                            interaction_data['resnr'] = interaction['resnr']
                        if interaction_data['reschain'] == 'X' and 'reschain' in interaction:
                            interaction_data['reschain'] = interaction['reschain']
                            interacting_chains.add(interaction['reschain'])
                    
                    # Record interacting residue
                    res_id = f"{interaction_data['resnr']}{interaction_data['reschain']}"
                    if res_id != "0X":  # Only add if we have valid info
                        interacting_res.add(res_id)
                    
                    # Add specific attributes for certain interaction types
                    if itype == 'pi_pi_stacking':
                        if isinstance(interaction, dict) and 'angle' in interaction:
                            interaction_data['angle'] = interaction['angle']
                    if isinstance(interaction, dict) and 'type' in interaction:
                        interaction_data['type'] = interaction['type']

                    # Atom names, so the report says which atoms are in contact
                    if isinstance(interaction, dict):
                        for key, name_key in (('protein_atom', 'protein_atom_name'),
                                              ('ligand_atom', 'ligand_atom_name')):
                            atom = interaction.get(key)
                            if hasattr(atom, 'get_name'):
                                interaction_data[name_key] = atom.get_name().strip()

                    processed_interactions[itype].append(interaction_data)
                    
                except Exception as e:
                    logger.warning(f"Warning: Error processing interaction: {str(e)}")
                    continue
        
        # Context for the report header and the residue table
        try:
            import os as _os
            ligand_info['structure'] = _os.path.basename(str(self.structure_file))
            ligand_info['heavy_atoms'] = sum(
                1 for a in self.ligand_atoms if (a.element or '').strip().upper() != 'H')
            ligand_info['solvent_accessible'] = set(getattr(self, 'solvent_accessible', ()) or ())
        except Exception:
            pass

        # Update ligand info with collected data
        ligand_info['interacting_chains'] = list(interacting_chains) if interacting_chains else []
        ligand_info['interacting_res'] = list(interacting_res) if interacting_res else []
        
        # Get interacting residues from the interacting_residues set
        try:
            if hasattr(self, 'interacting_residues') and self.interacting_residues:
                for res_id in self.interacting_residues:
                    if isinstance(res_id, tuple) and len(res_id) == 2:
                        restype, resnr = res_id
                        if res_id in self.protein_residues:
                            res = self.protein_residues[res_id]
                            if hasattr(res, 'get_parent') and callable(res.get_parent):
                                parent = res.get_parent()
                                if parent and hasattr(parent, 'id'):
                                    chain = parent.id
                                    # Create a sample interaction entry for each residue
                                    #for itype in processed_interactions:
                                        #processed_interactions[itype].append({
                                            #'restype': restype,
                                            #'resnr': resnr,
                                            #'reschain': chain,
                                            #'distance': 3.5  # Use a reasonable default distance
                                        #})
                                    #interacting_res.add(f"{resnr}{chain}")
                                    #interacting_chains.add(chain)
        except Exception as e:
            logger.warning(f"Error processing interacting_residues set: {str(e)}")
        
        # Update ligand info with final data
        ligand_info['interacting_chains'] = list(interacting_chains) if interacting_chains else []
        ligand_info['interacting_res'] = list(interacting_res) if interacting_res else []
        
        # Generate report
        try:
            report = generate_plip_report(ligand_info, processed_interactions, output_file)
            return report
        except Exception as e:
            logger.warning(f"Error in report generation: {str(e)}")
            # Create simple fallback report
            simple_report = [
                "=============================================================================",
                "PandaMap Interaction Report",
                "=============================================================================",
                "",
                f"Ligand: {ligand_info['hetid']}:{ligand_info['chain']}:{ligand_info['position']}",
                "",
                "Interactions summary:",
            ]
            
            for itype, ilist in processed_interactions.items():
                if ilist:
                    simple_report.append(f"- {itype}: {len(ilist)}")
            
            simple_report.append("=============================================================================")
            
            if output_file:
                with open(output_file, 'w', encoding='utf-8') as f:
                    f.write('\n'.join(simple_report))
                    
            return '\n'.join(simple_report)
    
    
    def generate_enhanced_interaction_report(self, output_file=None):
        """Generate a comprehensive, well-organized interaction report with categorized sections."""
        
        if not output_file:
            base_name = os.path.splitext(os.path.basename(self.structure_file))[0]
            output_file = f"{base_name}_detailed_interactions.txt"
        
        # Prepare ligand info
        ligand_info = {
            'hetid': getattr(self.ligand_residue, 'resname', 'UNK'),
            'chain': 'X',
            'position': 0,
            'longname': getattr(self.ligand_residue, 'resname', 'Unknown'),
            'type': 'LIGAND',
            'interacting_chains': [],
            'interacting_res': []
        }
        
        # Try to get chain and position
        try:
            if hasattr(self.ligand_residue, 'parent') and self.ligand_residue.parent:
                ligand_info['chain'] = self.ligand_residue.parent.id
        except:
            pass
        
        try:
            if hasattr(self.ligand_residue, 'id') and isinstance(self.ligand_residue.id, tuple):
                ligand_info['position'] = self.ligand_residue.id[1]
        except:
            pass
        
        # Collect unique chains and residues
        chains = set()
        residues = set()
        
        for itype, interactions in self.interactions.items():
            for interaction in interactions:
                if 'protein_residue' in interaction:
                    res = interaction['protein_residue']
                    if hasattr(res, 'get_parent') and callable(res.get_parent):
                        parent = res.get_parent()
                        if parent and hasattr(parent, 'id'):
                            chains.add(parent.id)
                    if hasattr(res, 'id') and isinstance(res.id, tuple) and len(res.id) > 1:
                        residues.add(f"{res.resname}{res.id[1]}")
        
        ligand_info['interacting_chains'] = list(chains)
        ligand_info['interacting_res'] = list(residues)
        
        # Start building the report
        report = []
        report.append("=============================================================================")
        report.append(f"PandaMap Interaction Report")
        report.append("=============================================================================")
        report.append("")
        report.append(f"Ligand: {ligand_info['hetid']}:{ligand_info['chain']}:{ligand_info['position']}")
        report.append(f"File: {os.path.basename(self.structure_file)}")
        report.append("")
        
        # Summary statistics
        total_interactions = sum(len(ints) for ints in self.interactions.values())
        report.append(f"Total interactions detected: {total_interactions}")
        report.append(f"Interacting residues: {len(residues)}")
        report.append(f"Solvent accessible residues: {len(self.solvent_accessible)}")
        report.append("")
        
        # Group interactions by category
        categories = {
            "Electrostatic Interactions": ['hydrogen_bonds', 'ionic', 'salt_bridge', 'attractive_charge', 'repulsion'],
            "π-System Interactions": ['pi_pi_stacking', 'pi_cation', 'donor_pi', 'carbon_pi', 'amide_pi', 'alkyl_pi'],
            "Hydrophobic Interactions": ['hydrophobic'],
            "Metal Interactions": ['metal_coordination'],
            "Other Interactions": ['halogen_bonds', 'water_bridges', 'covalent']
        }
        
        # Print interaction summary by category
        report.append("INTERACTION SUMMARY BY CATEGORY")
        report.append("------------------------------")
        
        for category, itypes in categories.items():
            category_count = sum(len(self.interactions.get(it, [])) for it in itypes)
            if category_count > 0:
                report.append(f"{category}: {category_count}")
                # List specific types and counts
                for it in itypes:
                    if it in self.interactions and len(self.interactions[it]) > 0:
                        type_name = it.replace('_', ' ').title()
                        report.append(f"  - {type_name}: {len(self.interactions[it])}")
        
        report.append("")
        report.append("DETAILED INTERACTIONS")
        report.append("------------------------------")
        
        # Show detailed interactions for each category
        for category, itypes in categories.items():
            category_count = sum(len(self.interactions.get(it, [])) for it in itypes)
            if category_count > 0:
                report.append(f"\n{category}:")
                
                for it in itypes:
                    interactions = self.interactions.get(it, [])
                    if interactions:
                        type_name = it.replace('_', ' ').title()
                        report.append(f"\n  {type_name}:")
                        
                        for i, interaction in enumerate(interactions, 1):
                            res_info = "Unknown"
                            distance = "?.??"
                            
                            if 'protein_residue' in interaction:
                                res = interaction['protein_residue']
                                if hasattr(res, 'resname') and hasattr(res, 'id'):
                                    res_info = f"{res.resname}{res.id[1]}"
                                    if hasattr(res, 'get_parent') and callable(res.get_parent):
                                        parent = res.get_parent()
                                        if parent and hasattr(parent, 'id'):
                                            res_info = f"{res.resname}{res.id[1]}{parent.id}"
                            
                            if 'distance' in interaction:
                                distance = f"{interaction['distance']:.2f}Å"
                            
                            # Add solvent accessibility info
                            accessibility = ""
                            if 'protein_residue' in interaction:
                                res = interaction['protein_residue']
                                if hasattr(res, 'resname') and hasattr(res, 'id'):
                                    res_id = (res.resname, res.id[1])
                                    if res_id in self.solvent_accessible:
                                        accessibility = " (solvent accessible)"
                            
                            report.append(f"    {i}. {res_info} -- {distance} -- {ligand_info['hetid']}{accessibility}")
        
        # Add solvent accessibility section
        if self.solvent_accessible:
            report.append("\nSOLVENT ACCESSIBLE RESIDUES")
            report.append("------------------------------")
            for i, res_id in enumerate(sorted(self.solvent_accessible), 1):
                resname, resnum = res_id
                report.append(f"  {i}. {resname}{resnum}")
        
        # Write report to file
        report_text = "\n".join(report)
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(report_text)
        
        logger.info(f"Enhanced interaction report saved to: {output_file}")
        return report_text
    # Add methods to the class
    setattr(cls, 'detect_interactions_improved', detect_interactions_improved)
    setattr(cls, 'generate_interaction_report', generate_interaction_report)
    return cls

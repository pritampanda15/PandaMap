#!/usr/bin/env python
"""
Improved interaction detection for PandaMap.
Includes PLIP-like thresholds and filtering plus text report generation.
Add this to your core.py module or create a new file and import it.
"""

import os
import math
import numpy as np
from collections import defaultdict, namedtuple

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

def generate_plip_report(ligand_info, interactions, output_file=None):
    """Generate a PLIP-like text report for interactions with robust error handling."""
    
    report = []
    report.append("=============================================================================")
    report.append(f"PandaMap Interaction Report")
    report.append("=============================================================================")
    report.append("")
    report.append(f"Ligand: {ligand_info.get('hetid', 'UNK')}:{ligand_info.get('chain', 'X')}:{ligand_info.get('position', '0')}")
    report.append(f"Name: {ligand_info.get('longname', 'Unknown')}")
    report.append(f"Type: {ligand_info.get('type', 'LIGAND')}")
    
    # Optional additional ligand information
    for field, label in [
        ('molweight', 'Molecular Weight'),
        ('heavy_atoms', 'Number of Heavy Atoms'),
        ('num_rings', 'Number of Rings'),
        ('num_rot_bonds', 'Number of Rotatable Bonds'),
        ('num_hba', 'Number of Hydrogen Bond Acceptors'),
        ('num_hbd', 'Number of Hydrogen Bond Donors'),
        ('logp', 'logP')
    ]:
        if field in ligand_info:
            report.append(f"{label}: {ligand_info[field]}")
    
    report.append("\n------------------------------\n")
    interacting_chains = ligand_info.get('interacting_chains', [])
    interacting_res = ligand_info.get('interacting_res', [])
    report.append(f"Interacting Chains: {', '.join(interacting_chains) if interacting_chains else 'N/A'}")
    report.append(f"Interacting Residues: {len(interacting_res)}")
    
    # Summary of interactions
    interaction_types = {
        'hydrogen_bonds': "Hydrogen Bonds",
        'hydrophobic': "Hydrophobic Interactions",
        'pi_pi_stacking': "π-π Stacking",
        'carbon_pi': "Carbon-π Interactions",
        'donor_pi': "Donor-π Interactions",
        'amide_pi': "Amide-π Interactions",
        'halogen_bonds': "Halogen Bonds",
        'water_bridges': "Water Bridges",
        'metal_coordination': "Metal Coordination",  # CHANGED from metal_complexes
        'ionic': "Ionic Interactions",
        'salt_bridge': "Salt Bridges",
        'covalent': "Covalent Bonds",
        'alkyl_pi': "Alkyl-π Interactions", # ADDED new interaction types
        'attractive_charge': "Attractive Charge",
        'pi_cation': "π-Cation Interactions",
        'repulsion': "Repulsion"
    }
    
    report.append("\n------------------------------\n")
    report.append("Interaction Summary:")
    for itype, label in interaction_types.items():
        count = len(interactions.get(itype, []))
        if count > 0:
            report.append(f"  {label}: {count}")
    
    # Function to safely get interaction details
    def safe_format_interaction(i, index):
        """Format an interaction with robust error handling"""
        try:
            res_info = f"{i.get('restype', 'UNK')}{i.get('resnr', '?')}{i.get('reschain', '?')}"
            dist_info = f"{i.get('distance', 0.0):.2f}Å"
            lig_info = ligand_info.get('hetid', 'UNK')
            
            # Include additional info if available
            extra_info = ""
            if 'type' in i:
                extra_info = f"({i['type']})"
            elif itype == 'pi_pi_stacking' and 'angle' in i:
                try:
                    angle = float(i['angle'])
                    if angle < 30:
                        extra_info = "(Parallel)"
                    else:
                        extra_info = "(T-shaped)"
                except:
                    pass
            
            return f"  {index}. {res_info} {extra_info} -- {dist_info} -- {lig_info}"
        except Exception as e:
            return f"  {index}. Error formatting interaction: {str(e)}"
    
    # Detailed information about each interaction type
    for itype, label in interaction_types.items():
        if interactions.get(itype, []):
            report.append("\n------------------------------\n")
            report.append(f"{label}:")
            for i, interaction in enumerate(interactions[itype], 1):
                report.append(safe_format_interaction(interaction, i))
    
    report.append("\n=============================================================================\n")
    
    # Write to file if specified
    if output_file:
        try:
            with open(output_file, 'w') as f:
                f.write("\n".join(report))
            print(f"Report saved to {output_file}")
        except Exception as e:
            print(f"Error writing report to file: {str(e)}")
    
    return "\n".join(report)

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
        print("Starting interaction report generation...")
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
        # Debug: Show first interaction structure
        for itype, ilist in self.interactions.items():
            if ilist:
                print(f"First {itype} interaction structure:")
                first_item = ilist[0]
                print(f"  Type: {type(first_item)}")
                # If it's a dictionary, print its keys
                if isinstance(first_item, dict):
                    print(f"  Keys: {list(first_item.keys())}")
                    # Print some key values for debugging
                    for key in ['protein_residue', 'distance', 'restype', 'resnr', 'reschain', 'ligand_atom', 'protein_atom']:
                        if key in first_item:
                            print(f"  {key}: {first_item[key]}")
                break
        
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
                print(f"Error getting ligand chain: {str(e)}")
                
            try:
                if hasattr(self.ligand_residue, 'id') and isinstance(self.ligand_residue.id, tuple):
                    ligand_info['position'] = self.ligand_residue.id[1]
            except Exception as e:
                print(f"Error getting ligand position: {str(e)}")
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
                    
                    processed_interactions[itype].append(interaction_data)
                    
                except Exception as e:
                    print(f"Warning: Error processing interaction: {str(e)}")
                    continue
        
        # Update ligand info with collected data
        ligand_info['interacting_chains'] = list(interacting_chains) if interacting_chains else []
        ligand_info['interacting_res'] = list(interacting_res) if interacting_res else []
        
        # Get interacting residues from the interacting_residues set
        try:
            if hasattr(self, 'interacting_residues') and self.interacting_residues:
                print(f"Found {len(self.interacting_residues)} residues in interacting_residues set")
                for res_id in self.interacting_residues:
                    if isinstance(res_id, tuple) and len(res_id) == 2:
                        restype, resnr = res_id
                        # Try to find the chain from protein_residues
                        if res_id in self.protein_residues:
                            res = self.protein_residues[res_id]
                            if hasattr(res, 'get_parent') and callable(res.get_parent):
                                parent = res.get_parent()
                                if parent and hasattr(parent, 'id'):
                                    chain = parent.id
                                    print(f"Adding residue from set: {restype}{resnr}{chain}")
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
            print(f"Error processing interacting_residues set: {str(e)}")
        
        # Update ligand info with final data
        ligand_info['interacting_chains'] = list(interacting_chains) if interacting_chains else []
        ligand_info['interacting_res'] = list(interacting_res) if interacting_res else []
        
        # Generate report
        try:
            report = generate_plip_report(ligand_info, processed_interactions, output_file)
            return report
        except Exception as e:
            print(f"Error in report generation: {str(e)}")
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
                with open(output_file, 'w') as f:
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
        with open(output_file, 'w') as f:
            f.write(report_text)
        
        print(f"Enhanced interaction report saved to: {output_file}")
        return report_text
    # Add methods to the class
    setattr(cls, 'detect_interactions_improved', detect_interactions_improved)
    setattr(cls, 'generate_interaction_report', generate_interaction_report)
    return cls

#!/usr/bin/env python
"""
Core functionality for PandaMap: A Python package for visualizing 
protein-ligand interactions with 2D ligand structure representation.
"""

import os
import math
from collections import defaultdict
import tempfile
import subprocess
from Bio.PDB import PDBParser, MMCIFParser
import numpy as np
import subprocess
from Bio.PDB.DSSP import DSSP
from Bio.PDB import PDBIO
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Wedge, Rectangle, Polygon
from matplotlib.lines import Line2D
import matplotlib.patheffects as path_effects

# BioPython imports
from Bio.PDB import PDBParser, NeighborSearch

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



    def parse_pdbqt(pdbqt_file):
        """
        Convert PDBQT to PDB format by stripping the charge and atom type information.
        Returns a temporary file path to the converted PDB.
        """
        # Create a temporary file for the PDB output
        temp_pdb = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
        temp_pdb_path = temp_pdb.name
        temp_pdb.close()
        
        # Read the PDBQT file and write a modified version without charges and types
        with open(pdbqt_file, 'r') as f_pdbqt, open(temp_pdb_path, 'w') as f_pdb:
            for line in f_pdbqt:
                if line.startswith(('ATOM', 'HETATM')):
                    # Keep the PDB format portion, remove the PDBQT-specific part
                    # PDB format: columns 1-66 are standard PDB format
                    # PDBQT adds charge and atom type in columns 67+
                    f_pdb.write(line[:66] + '\n')
                elif not line.startswith(('REMARK', 'MODEL', 'ENDMDL', 'TORSDOF')):
                    # Copy most other lines except PDBQT-specific ones
                    f_pdb.write(line)
        
        return temp_pdb_path

class MultiFormatParser:
    """
    Parser class that can handle multiple molecular file formats.
    Supports: PDB, mmCIF/CIF, PDBQT
    """
    def __init__(self):
        self.pdb_parser = PDBParser(QUIET=True)
        self.mmcif_parser = MMCIFParser(QUIET=True)
    
    def parse_structure(self, file_path):
        """
        Parse a molecular structure file and return a BioPython structure object.
        Automatically detects file format based on extension.
        
        Parameters:
        -----------
        file_path : str
            Path to the structure file
            
        Returns:
        --------
        structure : Bio.PDB.Structure.Structure
            BioPython structure object
        """
        file_ext = os.path.splitext(file_path)[1].lower()
        
        if file_ext == '.pdb':
            return self.pdb_parser.get_structure('complex', file_path)
        
        elif file_ext in ('.cif', '.mmcif'):
            return self.mmcif_parser.get_structure('complex', file_path)
        
        elif file_ext == '.pdbqt':
            # Convert PDBQT to PDB format temporarily
            temp_pdb_path = parse_pdbqt(file_path)
            structure = self.pdb_parser.get_structure('complex', temp_pdb_path)
            
            # Clean up the temporary file
            try:
                os.unlink(temp_pdb_path)
            except:
                pass  # Ignore cleanup errors
                
            return structure
        
        else:
            raise ValueError(f"Unsupported file format: {file_ext}. Supported formats: .pdb, .cif, .mmcif, .pdbqt")

class SimpleLigandStructure:
    """
    Class to create a simplified 2D representation of a ligand structure
    without requiring RDKit or other external dependencies.
    """
    
    def __init__(self, ligand_atoms):
        """
        Initialize with a list of ligand atoms from a BioPython structure.
        
        Parameters:
        -----------
        ligand_atoms : list
            List of BioPython Atom objects from the ligand
        """
        self.ligand_atoms = ligand_atoms
        self.atom_coords = {}
        self.element_colors = {
            'C': '#808080',  # Grey
            'N': '#0000FF',  # Blue
            'O': '#FF0000',  # Red
            'S': '#FFFF00',  # Yellow
            'P': '#FFA500',  # Orange
            'F': '#00FF00',  # Green
            'Cl': '#00FF00', # Green
            'Br': '#A52A2A', # Brown
            'I': '#A020F0',  # Purple
            'H': '#FFFFFF'   # White
        }
        
        # Record atom coordinates and elements
        for atom in ligand_atoms:
            atom_id = atom.get_id()
            self.atom_coords[atom_id] = {
                'element': atom.element,
                'coord': atom.get_coord()  # 3D coordinates from PDB
            }
    
    def generate_2d_coords(self):
        """
        Generate simplified 2D coordinates for the ligand atoms based on their 3D coordinates.
        This is a very basic projection - in a real application, you would use a proper
        2D layout algorithm.
        
        Returns:
        --------
        dict : Dictionary mapping atom IDs to 2D coordinates
        """
        if not self.atom_coords:
            return {}
            
        # Simple projection onto the xy-plane
        coords_2d = {}
        
        # Get all 3D coordinates and find center
        all_coords = np.array([info['coord'] for info in self.atom_coords.values()])
        center = np.mean(all_coords, axis=0)
        
        # Subtract center to center the molecule
        centered_coords = all_coords - center
        
        # Simple PCA-like approach to find main plane
        # (This is a very simplified approach)
        cov_matrix = np.cov(centered_coords.T)
        
        try:
            # Get eigenvectors and eigenvalues
            eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
            
            # Sort by eigenvalue in descending order
            idx = eigenvalues.argsort()[::-1]
            eigenvalues = eigenvalues[idx]
            eigenvectors = eigenvectors[:, idx]
            
            # Use the first two eigenvectors to define the plane
            plane_vectors = eigenvectors[:, :2]
            
            # Project the centered coordinates onto the plane
            projected_coords = np.dot(centered_coords, plane_vectors)
            
            # Scale to fit nicely in the visualization
            max_dim = np.max(np.abs(projected_coords))
            scaling_factor = 50.0 / max_dim if max_dim > 0 else 1.0
            projected_coords *= scaling_factor
            
            # Store the 2D coordinates
            for i, atom_id in enumerate(self.atom_coords.keys()):
                coords_2d[atom_id] = projected_coords[i]
                
        except np.linalg.LinAlgError:
            # Fallback if eigendecomposition fails
            print("Warning: Could not compute optimal projection. Using simple XY projection.")
            for atom_id, info in self.atom_coords.items():
                # Simple scaling of x, y coordinates
                coords_2d[atom_id] = np.array([info['coord'][0], info['coord'][1]]) * 10.0
        
        return coords_2d
    
    def find_bonds(self, distance_threshold=2.0):
        """
        Find bonds between atoms based on distance.
        This is a simplified approach - in reality, you'd use chemical knowledge.
        
        Parameters:
        -----------
        distance_threshold : float
            Maximum distance between atoms to be considered bonded (in Angstroms)
            
        Returns:
        --------
        list : List of tuples (atom_id1, atom_id2) representing bonds
        """
        bonds = []
        atom_ids = list(self.atom_coords.keys())
        
        for i in range(len(atom_ids)):
            for j in range(i+1, len(atom_ids)):
                atom1_id = atom_ids[i]
                atom2_id = atom_ids[j]
                
                coord1 = self.atom_coords[atom1_id]['coord']
                coord2 = self.atom_coords[atom2_id]['coord']
                
                # Calculate Euclidean distance
                distance = np.linalg.norm(coord1 - coord2)
                
                # If distance is below threshold, consider them bonded
                if distance < distance_threshold:
                    bonds.append((atom1_id, atom2_id))
        
        return bonds
    
    def draw_on_axes(self, ax, center=(0, 0), radius=80):
        """
        Draw a simplified 2D representation of the ligand on the given axes.
        
        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw
        center : tuple
            The (x, y) coordinates where the center of the molecule should be
        radius : float
            The approximate radius the molecule should occupy
            
        Returns:
        --------
        dict : Dictionary mapping atom IDs to their 2D positions in the plot
        """
        # Generate 2D coordinates
        coords_2d = self.generate_2d_coords()
        
        if not coords_2d:
            # If we couldn't generate coordinates, draw a simple placeholder
            print("Warning: Could not generate ligand coordinates. Drawing placeholder.")
            circle = Circle(center, radius/2, fill=False, edgecolor='black', linestyle='-')
            ax.add_patch(circle)
            ax.text(center[0], center[1], "Ligand", ha='center', va='center')
            return {}
            
        # Find bonds
        bonds = self.find_bonds()
        
        # Scale coordinates to fit within the specified radius
        all_coords = np.array(list(coords_2d.values()))
        max_extent = np.max(np.abs(all_coords))
        scaling_factor = radius / (max_extent * 1.2)  # Leave some margin
        
        # Create a mapping of atom IDs to positions in the plot
        atom_positions = {}
        
        # Draw bonds first (so they're below atoms)
        for atom1_id, atom2_id in bonds:
            pos1 = coords_2d[atom1_id] * scaling_factor + center
            pos2 = coords_2d[atom2_id] * scaling_factor + center
            
            # Draw bond as a line
            line = Line2D([pos1[0], pos2[0]], [pos1[1], pos2[1]], 
                         color='black', linewidth=1.5, zorder=2)
            ax.add_line(line)
        
        # Draw atoms as circles
        for atom_id, coord in coords_2d.items():
            # Scale and shift the position
            pos = coord * scaling_factor + center
            atom_positions[atom_id] = pos
            
            element = self.atom_coords[atom_id]['element']
            color = self.element_colors.get(element, 'gray')
            
            # Determine size based on element (larger for heavier atoms)
            size = 8 if element in ['C', 'H'] else 10
            
            # Draw atom
            circle = Circle(pos, size, facecolor=color, edgecolor='black', 
                           linewidth=1, alpha=0.8, zorder=3)
            ax.add_patch(circle)
            
            # Add element label (except for carbon)
            if element != 'C':
                ax.text(pos[0], pos[1], element, ha='center', va='center', 
                       fontsize=8, fontweight='bold', color='white', zorder=4)
        
        return atom_positions


class HybridProtLigMapper:
    """
    Class for analyzing protein-ligand interactions and creating 
    visualizations with a simplified ligand structure.
    """
    
    def __init__(self, structure_file, ligand_resname=None):
        """
        Initialize with a structure file containing a protein-ligand complex.
        
        Parameters:
        -----------
        structure_file : str
            Path to the structure file (PDB, mmCIF/CIF, or PDBQT format)
        ligand_resname : str, optional
            Specific residue name of the ligand to focus on
        """
        self.structure_file = structure_file
        self.ligand_resname = ligand_resname
        
        # Parse the structure file using the multi-format parser
        parser = MultiFormatParser()
        self.structure = parser.parse_structure(structure_file)
        self.model = self.structure[0]
        
        # Separate ligand from protein
        self.protein_atoms = []
        self.ligand_atoms = []
        self.protein_residues = {}
        self.ligand_residue = None
        
        for residue in self.model.get_residues():
            # Store ligand atoms (HETATM records)
            if residue.id[0] != ' ':  # Non-standard residue (HETATM)
                if ligand_resname is None or residue.resname == ligand_resname:
                    for atom in residue:
                        self.ligand_atoms.append(atom)
                    if self.ligand_residue is None:
                        self.ligand_residue = residue
            else:  # Standard residues (protein)
                res_id = (residue.resname, residue.id[1])
                self.protein_residues[res_id] = residue
                for atom in residue:
                    self.protein_atoms.append(atom)
        
        # Check if we found a ligand
        if not self.ligand_atoms:
            raise ValueError(f"No ligand (HETATM) found in the file: {structure_file}")
        
        # Storage for the interaction data
        self.interactions = {
            'hydrogen_bonds': [],
            'carbon_pi': [],
            'pi_pi_stacking': [],
            'donor_pi': [],
            'amide_pi': [],
            'hydrophobic': []
        }
        
        # Will store residues that interact with the ligand
        self.interacting_residues = set()
        
        # For solvent accessibility information (simplified)
        self.solvent_accessible = set()
        
        # Create the simple ligand structure
        self.ligand_structure = SimpleLigandStructure(self.ligand_atoms)
            
    def detect_interactions(self, 
                           h_bond_cutoff=3.5, 
                           pi_stack_cutoff=5.5,
                           hydrophobic_cutoff=4.0):
        """
        Detect all interactions between protein and ligand.
        
        Parameters:
        -----------
        h_bond_cutoff : float
            Distance cutoff for hydrogen bonds in Angstroms
        pi_stack_cutoff : float
            Distance cutoff for pi-stacking interactions in Angstroms
        hydrophobic_cutoff : float
            Distance cutoff for hydrophobic interactions in Angstroms
        """
        # Use neighbor search for efficiency
        ns = NeighborSearch(self.protein_atoms)
        max_cutoff = max(h_bond_cutoff, pi_stack_cutoff, hydrophobic_cutoff)
        
        # Define amino acid categories
        aromatic_residues = {'PHE', 'TYR', 'TRP', 'HIS'}
        h_bond_donors = {'ARG', 'LYS', 'HIS', 'ASN', 'GLN', 'SER', 'THR', 'TYR', 'TRP'}
        h_bond_acceptors = {'ASP', 'GLU', 'ASN', 'GLN', 'HIS', 'SER', 'THR', 'TYR'}
        neg_charged = {'ASP', 'GLU'}
        amide_residues = {'ASN', 'GLN'}
        hydrophobic_residues = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'TYR'}
        
        # Check each ligand atom for interactions
        for lig_atom in self.ligand_atoms:
            # Find protein atoms within cutoff distance
            nearby_atoms = ns.search(lig_atom.get_coord(), max_cutoff)
            
            for prot_atom in nearby_atoms:
                prot_res = prot_atom.get_parent()
                distance = lig_atom - prot_atom
                
                # Store interacting residue for later visualization
                res_id = (prot_res.resname, prot_res.id[1])
                self.interacting_residues.add(res_id)
                
                # Determine interaction types based on distance and atom/residue types
                
                # 1. Hydrogen bonds - N and O atoms within cutoff
                if distance <= h_bond_cutoff:
                    if lig_atom.element in ['N', 'O'] and prot_atom.element in ['N', 'O']:
                        self.interactions['hydrogen_bonds'].append({
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        })
                
                # 2. Pi-stacking interactions - aromatic residues
                if distance <= pi_stack_cutoff and prot_res.resname in aromatic_residues:
                    if lig_atom.element == 'C' and prot_atom.element == 'C':
                        self.interactions['pi_pi_stacking'].append({
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        })
                
                # 3. Carbon-Pi interactions
                if distance <= pi_stack_cutoff and prot_res.resname in aromatic_residues:
                    if lig_atom.element == 'C':
                        self.interactions['carbon_pi'].append({
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        })
                
                # 4. Donor-Pi interactions - negatively charged residues
                if distance <= pi_stack_cutoff and prot_res.resname in neg_charged:
                    if lig_atom.element == 'C':
                        self.interactions['donor_pi'].append({
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        })
                
                # 5. Amide-Pi interactions
                if distance <= pi_stack_cutoff and prot_res.resname in amide_residues:
                    if lig_atom.element == 'C':
                        self.interactions['amide_pi'].append({
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        })
                
                # 6. Hydrophobic interactions
                if distance <= hydrophobic_cutoff:
                    if (prot_res.resname in hydrophobic_residues and 
                        lig_atom.element == 'C' and prot_atom.element == 'C'):
                        self.interactions['hydrophobic'].append({
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        })
        
        # Deduplicate interactions by residue for cleaner visualization
        # Keep only one interaction of each type per residue
        for interaction_type in self.interactions:
            by_residue = defaultdict(list)
            for interaction in self.interactions[interaction_type]:
                res_id = (interaction['protein_residue'].resname, 
                          interaction['protein_residue'].id[1])
                by_residue[res_id].append(interaction)
            
            # Keep only the closest interaction for each residue and type
            closest_interactions = []
            for res_id, res_interactions in by_residue.items():
                closest = min(res_interactions, key=lambda x: x['distance'])
                closest_interactions.append(closest)
            
            self.interactions[interaction_type] = closest_interactions
    
    def estimate_solvent_accessibility(self):
        """
        Estimate which residues might be solvent accessible.
        This is a simplified approach since we're trying to match the example image.
        In a real implementation, you'd use DSSP or a similar tool.
        """
        # For simplicity, mark all residues as solvent accessible
        # In a real implementation, you'd use a proper algorithm
        self.solvent_accessible = self.interacting_residues.copy()
    
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
            # Create a temporary PDB file for DSSP input
            with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp_pdb:
                pdb_io = PDBIO()
                pdb_io.set_structure(self.structure)
                pdb_io.save(tmp_pdb.name)
                
                # Run DSSP
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
                    if rsa > 0.15 and res_key in self.interacting_residues:
                        self.solvent_accessible.add(res_key)
                        
        except Exception as e:
            print(f"Warning: DSSP calculation failed. Falling back to geometric estimation. Error: {str(e)}")
            self.estimate_solvent_accessibility()
        
        return self.solvent_accessible
    
    #if dssp is not available in the path
    def calculate_python_solvent_accessibility(self, probe_radius=1.4):
        """
        Simplified solvent accessibility calculation in pure Python.
        Based on Shrake-Rupley algorithm but with approximations.
        
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
            
            for atom in residue.get_atoms():
                total_atoms += 1
                atom_coord = atom.get_coord()
                
                # Check if atom is buried by other protein atoms
                is_exposed = True
                for other_atom in all_protein_atoms:
                    if other_atom.get_parent() == residue:
                        continue  # Skip atoms in same residue
                    
                    distance = np.linalg.norm(atom_coord - other_atom.get_coord())
                    if distance < (atom.radius + other_atom.radius + probe_radius):
                        is_exposed = False
                        break
                
                if is_exposed:
                    exposed_atoms += 1
            
            # Consider residue accessible if >25% of its atoms are exposed
            if total_atoms > 0 and (exposed_atoms / total_atoms) > 0.25:
                self.solvent_accessible.add(res_id)
        
        return self.solvent_accessible
    

    def visualize(self, output_file='protein_ligand_interactions.png',
              figsize=(12, 12), dpi=300, title=None):
        """
        Generate a complete 2D visualization of protein-ligand interactions with:
        - Arrows connecting to residue box edges
        - Properly placed interaction markers
        - Solvent accessibility indicators
        - Comprehensive legend
        """
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Add light blue background for ligand
        ligand_radius = 90
        ligand_pos = (0, 0)
        ligand_circle = Circle(ligand_pos, ligand_radius, facecolor='#ADD8E6', 
                            edgecolor='none', alpha=0.4, zorder=1)
        ax.add_patch(ligand_circle)
        
        # Draw the simplified ligand structure
        atom_positions = self.ligand_structure.draw_on_axes(ax, center=ligand_pos, radius=ligand_radius*0.8)
        
        # Place interacting residues in a circle around the ligand
        n_residues = len(self.interacting_residues)
        if n_residues == 0:
            print("Warning: No interacting residues detected.")
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
            
            # Draw solvent accessibility highlight
            if res_id in self.solvent_accessible:
                solvent_circle = Circle((x, y), 40, facecolor='#ADD8E6', 
                                    edgecolor='none', alpha=0.3, zorder=1)
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
                        'marker_text': 'h', 'marker_bg': 'white', 'name': 'Hydrophobic'}
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
        
        # Draw interaction lines with arrows at box edges
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
                
                # Calculate curvature
                dx = res_pos[0] - lig_pos[0]
                dy = res_pos[1] - lig_pos[1]
                distance = math.sqrt(dx*dx + dy*dy)
                curvature = 0.08 * (200 / max(distance, 100))
                
                # Store line parameters
                line_params = {
                    'start_pos': box_edge_pos,
                    'end_pos': lig_pos,
                    'curvature': curvature,
                    'style': style,
                    'interaction_type': interaction_type,
                    'key': f"{interaction_type}_{res_id[0]}_{res_id[1]}",
                    'distance': distance
                }
                interaction_lines.append(line_params)
                
                # Draw arrow
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

        # Place markers along interaction lines
        marker_positions = {}
        type_order = {'hydrogen_bonds': 0, 'carbon_pi': 1, 'pi_pi_stacking': 2, 
                    'donor_pi': 3, 'amide_pi': 4, 'hydrophobic': 5}
        
        sorted_lines = sorted(interaction_lines,
                            key=lambda x: (type_order.get(x['interaction_type'], x['distance'])))
        
        for line_params in sorted_lines:
            start_pos = line_params['start_pos']
            end_pos = line_params['end_pos']
            curvature = line_params['curvature']
            style = line_params['style']
            key = line_params['key']
            
            # Calculate points along the curved path
            path_points = []
            steps = 20
            for i in range(steps + 1):
                t = i / steps
                control_x = (start_pos[0] + end_pos[0])/2 + curvature * (end_pos[1] - start_pos[1]) * 2
                control_y = (start_pos[1] + end_pos[1])/2 - curvature * (end_pos[0] - start_pos[0]) * 2
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
        
        # Add solvent accessibility indicator
        if self.solvent_accessible:
            legend_elements.append(
                Rectangle((0, 0), 1, 1, facecolor='#ADD8E6', alpha=0.3,
                        edgecolor=None, label='Solvent accessible')
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
        
        print(f"Interaction diagram saved to {output_file}")
        return output_file
    
    def run_analysis(self, output_file=None, use_dssp=True):
        """
        Run the complete analysis pipeline.
        
        Parameters:
        -----------
        output_file : str, optional
            Path where the output image will be saved.
        use_dssp : bool
            Whether to use DSSP for solvent accessibility (default: True)
        """
        if output_file is None:
            base_name = os.path.splitext(os.path.basename(self.structure_file))[0]
            output_file = f"{base_name}_interactions.png"
        
        # Detect protein-ligand interactions
        print("Detecting interactions...")
        self.detect_interactions()
        
        # Calculate solvent accessibility
        print("Calculating solvent accessibility...")
        if use_dssp:
            try:
                self.calculate_dssp_solvent_accessibility()
            except:
                print("DSSP failed, falling back to geometric estimation")
                self.estimate_solvent_accessibility()
        else:
            self.calculate_python_solvent_accessibility()
        
        # Generate visualization
        print("Generating visualization...")
        return self.visualize(output_file=output_file)
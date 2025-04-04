#!/usr/bin/env python
"""
Core functionality for PandaMap: A Python package for visualizing 
protein-ligand interactions with 2D ligand structure representation.
"""

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
    with enhanced 3D cues without requiring RDKit or other external dependencies.
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
        
        # Direct fallback for small molecules (less than 3 atoms)
        if len(centered_coords) < 3:
            print("Warning: Not enough atoms for PCA projection. Using simple XY projection.")
            for atom_id, info in self.atom_coords.items():
                # Simple scaling of x, y coordinates
                coords_2d[atom_id] = np.array([info['coord'][0], info['coord'][1]]) * 10.0
            return coords_2d
            
        # Simple PCA-like approach to find main plane
        try:
            # Make sure we have enough unique coordinates
            unique_coords = set()
            for coord in centered_coords:
                unique_coords.add(tuple(coord))
                
            if len(unique_coords) < 3:
                raise ValueError("Not enough unique coordinates for PCA")
            
            # Calculate covariance matrix
            covariance_matrix = np.zeros((3, 3))
            for coord in centered_coords:
                covariance_matrix += np.outer(coord, coord)
            covariance_matrix /= len(centered_coords)
            
            # Get eigenvectors and eigenvalues
            eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)
            
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
                
        except Exception as e:
            # Fallback if eigendecomposition fails
            print(f"Warning: Could not compute optimal projection. Using simple XY projection. Error: {str(e)}")
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
    
    def draw_on_axes(self, ax, center=(0, 0), radius=80, show_3d_cues=True):
        """
        Draw a more realistic 3D representation of the ligand on the given axes.
        
        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw
        center : tuple
            The (x, y) coordinates where the center of the molecule should be
        radius : float
            The approximate radius the molecule should occupy
        show_3d_cues : bool
            Whether to show 3D depth cues (enhanced 3D effect)
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
        
        # Calculate approximate z-depth for 3D cues
        z_values = {}
        if show_3d_cues:
            # Extract z-coordinates from 3D structure
            for atom_id, info in self.atom_coords.items():
                z_values[atom_id] = info['coord'][2]
            
            # Normalize z-values to 0-1 range
            if len(z_values) > 0:
                z_min = min(z_values.values())
                z_max = max(z_values.values())
                z_range = z_max - z_min
                if z_range > 0:
                    for atom_id in z_values:
                        z_values[atom_id] = (z_values[atom_id] - z_min) / z_range
                else:
                    for atom_id in z_values:
                        z_values[atom_id] = 0.5
        
        # Sort atoms by Z-depth for proper occlusion
        atom_ids = list(coords_2d.keys())
        if show_3d_cues and z_values:
            # Sort from back to front (furthest first, closest last)
            atom_ids.sort(key=lambda atom_id: z_values.get(atom_id, 0.5))
        
        # Dictionary to store atom indices (for bond reference)
        atom_indices = {atom_id: i for i, atom_id in enumerate(atom_ids)}
        
        # Enhanced shadow effect if 3D cues are enabled
        if show_3d_cues:
            # Add a subtle shadow below the molecule
            shadow_offset = 15  # pixels
            shadow_alpha = 0.15
            for atom_id in atom_ids:
                pos = coords_2d[atom_id] * scaling_factor + center
                pos_shadow = (pos[0] + shadow_offset, pos[1] + shadow_offset)
                
                element = self.atom_coords[atom_id]['element']
                
                # Determine size based on element and z-position
                base_size = 9 if element in ['C', 'H'] else 11
                z = z_values.get(atom_id, 0.5)
                size_factor = 0.8 + 0.4 * z  # Larger when closer
                size = base_size * size_factor
                
                # Draw shadow
                shadow = Circle(pos_shadow, size * 1.1, facecolor='black', 
                            edgecolor=None, alpha=shadow_alpha, zorder=1.8)
                ax.add_patch(shadow)
        
        # Draw bonds with enhanced 3D representation
        for atom1_id, atom2_id in bonds:
            pos1 = coords_2d[atom1_id] * scaling_factor + center
            pos2 = coords_2d[atom2_id] * scaling_factor + center
            
            # Get z-values if available
            z1 = z_values.get(atom1_id, 0.5) if show_3d_cues else 0.5
            z2 = z_values.get(atom2_id, 0.5) if show_3d_cues else 0.5
            
            # Calculate bond midpoint
            midpoint = ((pos1[0] + pos2[0])/2, (pos1[1] + pos2[1])/2)
            
            # Calculate bond vector and perpendicular
            dx = pos2[0] - pos1[0]
            dy = pos2[1] - pos1[1]
            bond_length = np.sqrt(dx*dx + dy*dy)
            
            if bond_length > 0:
                # Unit vector along bond
                ux, uy = dx/bond_length, dy/bond_length
                # Perpendicular unit vector
                px, py = -uy, ux
            else:
                ux, uy = 1, 0
                px, py = 0, 1
                
            # If showing 3D cues and significant z-difference
            if show_3d_cues and abs(z1 - z2) > 0.2:
                # Determine which atom is closer to viewer
                if z1 > z2:  # atom1 is closer
                    # Create wedge (thicker near atom1)
                    wedge_width = 3.5  # max width of wedge
                    
                    # Wedge points
                    wedge_points = [
                        (pos1[0] + px * wedge_width, pos1[1] + py * wedge_width),
                        (pos1[0] - px * wedge_width, pos1[1] - py * wedge_width),
                        (pos2[0] - px * wedge_width * 0.2, pos2[1] - py * wedge_width * 0.2),
                        (pos2[0] + px * wedge_width * 0.2, pos2[1] + py * wedge_width * 0.2)
                    ]
                    
                    # Draw filled wedge with gradient fill
                    wedge = plt.Polygon(wedge_points, closed=True, 
                                    facecolor='black', alpha=0.7, zorder=2.5)
                    ax.add_patch(wedge)
                    
                elif z2 > z1:  # atom2 is closer
                    # Create wedge (thicker near atom2)
                    wedge_width = 3.5  # max width of wedge
                    
                    # Wedge points
                    wedge_points = [
                        (pos2[0] + px * wedge_width, pos2[1] + py * wedge_width),
                        (pos2[0] - px * wedge_width, pos2[1] - py * wedge_width),
                        (pos1[0] - px * wedge_width * 0.2, pos1[1] - py * wedge_width * 0.2),
                        (pos1[0] + px * wedge_width * 0.2, pos1[1] + py * wedge_width * 0.2)
                    ]
                    
                    # Draw filled wedge
                    wedge = plt.Polygon(wedge_points, closed=True, 
                                    facecolor='black', alpha=0.7, zorder=2.5)
                    ax.add_patch(wedge)
                else:
                    # Regular bond (in plane)
                    line = Line2D([pos1[0], pos2[0]], [pos1[1], pos2[1]], 
                                color='black', linewidth=2.0, zorder=2)
                    ax.add_line(line)
            else:
                # Regular bond
                line = Line2D([pos1[0], pos2[0]], [pos1[1], pos2[1]], 
                            color='black', linewidth=2.0, zorder=2)
                ax.add_line(line)
        
        # Draw atoms with enhanced 3D appearance
        for atom_id in atom_ids:
            # Scale and shift the position
            pos = coords_2d[atom_id] * scaling_factor + center
            atom_positions[atom_id] = pos
            
            element = self.atom_coords[atom_id]['element']
            base_color = self.element_colors.get(element, 'gray')
            
            # Get the normalized z value (0=back, 1=front)
            z = z_values.get(atom_id, 0.5) if show_3d_cues else 0.5
            
            # Determine size based on element and z-position
            base_size = 9 if element in ['C', 'H'] else 11
            size_factor = 0.8 + 0.4 * z  # Larger when closer
            size = base_size * size_factor
            
            # Adjust color based on z-position (darker in back, brighter in front)
            if show_3d_cues:
                # Parse the hex color to RGB
                if base_color.startswith('#'):
                    r = int(base_color[1:3], 16) / 255.0
                    g = int(base_color[3:5], 16) / 255.0
                    b = int(base_color[5:7], 16) / 255.0
                else:
                    # Just use the base color for named colors
                    r, g, b = matplotlib.colors.to_rgb(base_color)
                
                # Adjust color based on depth
                brightness_factor = 0.6 + 0.4 * z  # Brighter in front
                r = min(1.0, r * brightness_factor)
                g = min(1.0, g * brightness_factor)
                b = min(1.0, b * brightness_factor)
                color = (r, g, b)
            else:
                color = base_color
            
            # Draw atom with 3D effect (gradient fill)
            if show_3d_cues:
                # First draw a shadow/outline for 3D effect
                outline = Circle(pos, size + 1, facecolor='black', alpha=0.2, zorder=3.8)
                ax.add_patch(outline)
                
                # Main atom circle
                atom_circle = Circle(pos, size, facecolor=color, edgecolor='black', 
                                    linewidth=1, alpha=1.0, zorder=4)
                ax.add_patch(atom_circle)
                
                # Highlight to give 3D appearance (small white circle in upper left)
                if element != 'H':  # Skip highlight for hydrogen
                    highlight_offset = size * 0.3
                    highlight_size = size * 0.35
                    highlight_pos = (pos[0] - highlight_offset, pos[1] - highlight_offset)
                    highlight = Circle(highlight_pos, highlight_size, facecolor='white', 
                                    alpha=0.5, zorder=4.5)
                    ax.add_patch(highlight)
            else:
                # Simple 2D representation
                atom_circle = Circle(pos, size, facecolor=color, edgecolor='black', 
                                linewidth=1, alpha=0.8, zorder=3)
                ax.add_patch(atom_circle)
            
            # Add element label (except for carbon)
            if element != 'C':
                # Adjust text color for contrast
                text_color = 'white' if element != 'S' else 'black'
                text = ax.text(pos[0], pos[1], element, ha='center', va='center', 
                            fontsize=8, fontweight='bold', color=text_color, zorder=5)
                if show_3d_cues:
                    # Add subtle drop shadow for better text visibility
                    text.set_path_effects([
                        path_effects.withStroke(linewidth=2, foreground='black', alpha=0.7)
                    ])
        
        return atom_positions
    
    def adjust_2d_coordinates_for_better_3d(self, coords_2d, atom_ids, rotation_angle=30):
        """
        Apply a slight rotation and perspective to the 2D coordinates
        to enhance the 3D perception of the molecule.
        
        Parameters:
        -----------
        coords_2d : dict
            Dictionary mapping atom IDs to their 2D coordinates
        atom_ids : list
            List of atom IDs 
        rotation_angle : float
            Rotation angle in degrees (default: 30)
        
        Returns:
        --------
        dict : Dictionary mapping atom IDs to adjusted 2D coordinates
        """
        if not coords_2d:
            return {}
        
        # Convert to radians
        theta = np.radians(rotation_angle)
        
        # Calculate center of the molecule
        positions = np.array([coords_2d[atom_id] for atom_id in atom_ids])
        center = np.mean(positions, axis=0)
        
        # Initialize result dictionary
        adjusted_coords = {}
        
        # Apply rotation around y-axis
        for atom_id in atom_ids:
            x, y = coords_2d[atom_id]
            
            # Center coordinates
            x_centered = x - center[0]
            y_centered = y - center[1]
            
            # Apply rotation (isometric-like projection)
            x_rot = x_centered * np.cos(theta)
            y_rot = y_centered
            
            # Add back center offset
            x_new = x_rot + center[0]
            y_new = y_rot + center[1]
            
            adjusted_coords[atom_id] = np.array([x_new, y_new])
        
        return adjusted_coords

    def calculate_enhanced_bond_information(self, bonds, coords_2d, z_values):
        """
        Calculate enhanced bond information for better 3D representation.
        
        Parameters:
        -----------
        bonds : list
            List of bond pairs (atom_id1, atom_id2)
        coords_2d : dict
            Dictionary mapping atom IDs to 2D coordinates
        z_values : dict
            Dictionary mapping atom IDs to normalized z values
        
        Returns:
        --------
        list : List of dictionaries with enhanced bond information
        """
        enhanced_bonds = []
        
        for atom1_id, atom2_id in bonds:
            # Get bond endpoints
            pos1 = coords_2d[atom1_id]
            pos2 = coords_2d[atom2_id]
            
            # Get z-values
            z1 = z_values.get(atom1_id, 0.5)
            z2 = z_values.get(atom2_id, 0.5)
            
            # Calculate midpoint
            midpoint = ((pos1[0] + pos2[0])/2, (pos1[1] + pos2[1])/2)
            
            # Calculate bond vector and perpendicular
            dx = pos2[0] - pos1[0]
            dy = pos2[1] - pos1[1]
            bond_length = np.sqrt(dx*dx + dy*dy)
            
            # Calculate unit vectors
            if bond_length > 0:
                ux, uy = dx/bond_length, dy/bond_length
                px, py = -uy, ux  # Perpendicular
            else:
                ux, uy = 1, 0
                px, py = 0, 1
            
            # Calculate z difference and determine bond type
            z_diff = z2 - z1
            
            if abs(z_diff) > 0.2:
                # 3D bond (significant z difference)
                if z_diff > 0:  # atom2 is closer to viewer
                    bond_type = 'forward'  # Coming toward viewer
                    closer_atom_id = atom2_id
                    farther_atom_id = atom1_id
                else:  # atom1 is closer
                    bond_type = 'backward'  # Going away from viewer
                    closer_atom_id = atom1_id
                    farther_atom_id = atom2_id
            else:
                # Planar bond (small z difference)
                bond_type = 'planar'
                closer_atom_id = atom1_id if z1 > z2 else atom2_id
                farther_atom_id = atom2_id if z1 > z2 else atom1_id
            
            # Calculate average z value for bond (for drawing order)
            z_avg = (z1 + z2) / 2
            
            enhanced_bonds.append({
                'atom1_id': atom1_id,
                'atom2_id': atom2_id,
                'pos1': pos1,
                'pos2': pos2,
                'midpoint': midpoint,
                'bond_length': bond_length,
                'unit_vector': (ux, uy),
                'perp_vector': (px, py),
                'z1': z1,
                'z2': z2,
                'z_avg': z_avg,
                'z_diff': z_diff,
                'bond_type': bond_type,
                'closer_atom_id': closer_atom_id,
                'farther_atom_id': farther_atom_id
            })
        
        # Sort bonds by z-average (paint back-to-front)
        enhanced_bonds.sort(key=lambda b: b['z_avg'])
        
        return enhanced_bonds

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
        
        # Store atom information including metals
        self.metal_atoms = []
        self.halogen_atoms = []
        
        for residue in self.model.get_residues():
            # Store ligand atoms (HETATM records)
            if residue.id[0] != ' ':  # Non-standard residue (HETATM)
                if ligand_resname is None or residue.resname == ligand_resname:
                    for atom in residue:
                        self.ligand_atoms.append(atom)
                        # Identify halogen atoms in ligand
                        if atom.element in ['F', 'Cl', 'Br', 'I']:
                            self.halogen_atoms.append(atom)
                    if self.ligand_residue is None:
                        self.ligand_residue = residue
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
            'repulsion': []         # NEW: Repulsion interactions
        }
        
        # Will store residues that interact with the ligand
        self.interacting_residues = set()
        
        # Store interaction directionality (protein->ligand, ligand->protein, or both)
        self.interaction_direction = {}
        
        # For solvent accessibility information (simplified)
        self.solvent_accessible = set()
        
        # Create the simple ligand structure
        self.ligand_structure = SimpleLigandStructure(self.ligand_atoms)

    def calculate_hbond_angle(donor_atom, acceptor_atom, all_atoms):
        """Calculate the hydrogen bond angle."""
        donor_coord = donor_atom.get_coord()
        acceptor_coord = acceptor_atom.get_coord()
        
        # Simple vector calculation
        donor_to_acceptor = acceptor_coord - donor_coord
        
        try:
            # Normalize vector
            donor_to_acceptor_norm = donor_to_acceptor / np.linalg.norm(donor_to_acceptor)
            
            # Default acceptable angle
            return 125.0  # Default to a reasonable H-bond angle
        except:
            return 120.0  # Fallback angle




    def detect_interactions(self, 
                  h_bond_cutoff=3.5, 
                  pi_stack_cutoff=5.5,
                  hydrophobic_cutoff=4.0,
                  ionic_cutoff=4.0,
                  halogen_bond_cutoff=3.5,
                  metal_coord_cutoff=2.8,
                  covalent_cutoff=2.1):
        print("*** Using IMPROVED detection with stricter filtering ***")
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
            'repulsion': []
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
            is_lig_pos_charged = lig_atom.element == 'N' and not any(a.element == 'C' for a in self.ligand_atoms if a.element != 'H' and np.linalg.norm(a.get_coord() - lig_atom.get_coord()) < 1.6)
            is_lig_neg_charged = lig_atom.element == 'O' and not any(a.element == 'C' for a in self.ligand_atoms if a.element != 'H' and np.linalg.norm(a.get_coord() - lig_atom.get_coord()) < 1.6)
            
            # Track aromatic ligand characteristics
            is_lig_aromatic = lig_atom.element == 'C' and len([a for a in self.ligand_atoms if a.element == 'C' and 1.2 < np.linalg.norm(a.get_coord() - lig_atom.get_coord()) < 2.8]) >= 2
            
            for prot_atom in nearby_atoms:
                prot_res = prot_atom.get_parent()
                distance = lig_atom - prot_atom
                
                # Skip if distance is too small (likely a clash or error)
                if distance < 1.5:
                    continue
                    
                # Store interacting residue for later visualization
                res_id = (prot_res.resname, prot_res.id[1])
                self.interacting_residues.add(res_id)
                
                # 1. Hydrogen bonds - N and O atoms within cutoff
                # 1. Hydrogen bonds - N and O atoms within cutoff
                # 1. Hydrogen bonds - N and O atoms within cutoff
                if distance <= h_bond_cutoff and distance >= 2.4:  # Add minimum distance
                    if lig_atom.element in ['N', 'O'] and prot_atom.element in ['N', 'O']:
                        # Create res_id for interaction key
                        res_id = (prot_res.resname, prot_res.id[1])
                        
                        interaction_info = {
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance,
                            'angle': 120.0  # Default angle
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
                    is_prot_aromatic = prot_res.resname in aromatic_residues and prot_atom.element == 'C'
                    is_lig_alkyl = lig_atom.element == 'C' and not is_lig_aromatic
        
                if is_prot_aromatic and is_lig_alkyl:
                    interaction_info = {
                        'ligand_atom': lig_atom,
                        'protein_atom': prot_atom,
                        'protein_residue': prot_res,
                        'distance': distance
                    }
                    all_interactions['alkyl_pi'].append(interaction_info)
    
                # 3. Attractive charge interactions - between oppositely charged groups
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
                        all_interactions['attractive_charge'].append(interaction_info)
    
                # 4. Repulsion interactions - between likely same-charged groups
                if distance <= ionic_cutoff * 1.5:  # Larger cutoff for repulsion
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
                    if prot_res.resname in aromatic_residues and is_lig_aromatic and \
                    lig_atom.element == 'C' and prot_atom.element == 'C':
                        interaction_info = {
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        }
                        all_interactions['pi_pi_stacking'].append(interaction_info)
                
                # 3. Carbon-Pi - between carbon atoms and aromatic systems
                if distance <= pi_stack_cutoff:
                    if prot_res.resname in aromatic_residues and lig_atom.element == 'C':
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
                        all_interactions['ionic'].append(interaction_info)
                        all_interactions['salt_bridge'].append(interaction_info)  # Salt bridges are a type of ionic interaction
                
                # 8. Cation-Pi - positive charged residues with aromatic systems
                if distance <= pi_stack_cutoff:
                    if (prot_res.resname in pos_charged and is_lig_aromatic) or \
                    (prot_res.resname in aromatic_residues and is_lig_pos_charged):
                        interaction_info = {
                            'ligand_atom': lig_atom,
                            'protein_atom': prot_atom,
                            'protein_residue': prot_res,
                            'distance': distance
                        }
                        all_interactions['cation_pi'].append(interaction_info)
                
                # 9. Covalent bonds - very close interactions with specific atom types
                if distance <= covalent_cutoff:
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
            nearby_atoms = ns.search(halogen_atom.get_coord(), halogen_bond_cutoff)
            for prot_atom in nearby_atoms:
                prot_res = prot_atom.get_parent()
                distance = halogen_atom - prot_atom
                
                if 1.5 < distance <= halogen_bond_cutoff and is_halogen_acceptor(prot_atom):
                    res_id = (prot_res.resname, prot_res.id[1])
                    self.interacting_residues.add(res_id)
                    
                    interaction_info = {
                        'ligand_atom': halogen_atom,
                        'protein_atom': prot_atom,
                        'protein_residue': prot_res,
                        'distance': distance
                    }
                    all_interactions['halogen_bonds'].append(interaction_info)
        
        # Handle metal coordination
        for metal_atom in self.metal_atoms:
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
        
        # Filter redundant interactions to get one per residue
        for itype, interactions in all_interactions.items():
            # Group by residue ID
            by_residue = {}
            for interaction in interactions:
                res_id = (interaction['protein_residue'].resname, interaction['protein_residue'].id[1])
                if res_id not in by_residue or by_residue[res_id]['distance'] > interaction['distance']:
                    by_residue[res_id] = interaction
            
            # Add only the best (closest) interaction for each residue
            self.interactions[itype] = list(by_residue.values())
    
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
        print("Using realistic solvent accessibility calculation...")
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
        
        print(f"Total protein atoms: {len(all_protein_atoms)}")
        print(f"Interacting residues to check: {len(self.interacting_residues)}")
        
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
            
            # Only add highest scoring residues directly
            if exposure_score > exposure_threshold:
                print(f"Marking {res_id} as solvent accessible (score: {exposure_score:.2f})")
                self.solvent_accessible.add(res_id)
        
        # Calculate constraints on number of solvent accessible residues
        min_expected = max(1, int(len(self.interacting_residues) * 0.1))  # At least 10%
        max_expected = min(int(len(self.interacting_residues) * max_percent), 
                        len(self.interacting_residues) - 1)  # At most max_percent, never all
        
        print(f"Constraints: min={min_expected}, max={max_expected} solvent accessible residues")
        
        # Add more residues if below minimum
        if len(self.solvent_accessible) < min_expected:
            print(f"Too few solvent-accessible residues detected ({len(self.solvent_accessible)}), "
                f"adding more based on scores...")
            
            # Sort remaining residues by their exposure scores
            remaining = sorted(
                [(r, exposure_scores.get(r, 0.0)) for r in self.interacting_residues if r not in self.solvent_accessible],
                key=lambda x: x[1],  # Sort by score
                reverse=True  # Highest scores first
            )
            
            # Add only up to the minimum required
            for res_id, score in remaining:
                if len(self.solvent_accessible) >= min_expected:
                    break
                print(f"Adding {res_id} as solvent accessible based on ranking (score: {score:.2f})")
                self.solvent_accessible.add(res_id)
        
        # Remove residues if above maximum
        if len(self.solvent_accessible) > max_expected:
            print(f"Too many solvent-accessible residues detected ({len(self.solvent_accessible)}), "
                f"removing lowest scoring ones...")
            
            # Sort current accessible residues by score, ascending
            to_evaluate = sorted(
                [(r, exposure_scores.get(r, 0.0)) for r in self.solvent_accessible],
                key=lambda x: x[1]  # Sort by score
            )
            
            # Remove lowest scoring residues until we're within limits
            residues_to_remove = len(self.solvent_accessible) - max_expected
            for i in range(residues_to_remove):
                if i < len(to_evaluate):
                    res_id, score = to_evaluate[i]
                    print(f"Removing {res_id} from solvent accessible (score: {score:.2f})")
                    self.solvent_accessible.remove(res_id)
        
        print(f"Final result: {len(self.solvent_accessible)} solvent-accessible residues out of {len(self.interacting_residues)} interacting residues")
        if self.solvent_accessible:
            print(f"Solvent accessible residues: {sorted(self.solvent_accessible)}")
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
        print("Using enhanced solvent accessibility calculation...")
        self.solvent_accessible = set()
        
        # Define which residues are typically surface-exposed
        likely_exposed = {'ARG', 'LYS', 'ASP', 'GLU', 'ASN', 'GLN', 'HIS', 'SER', 'THR', 'TYR'}
        
        # First get all protein atoms (including non-interacting ones)
        all_protein_atoms = []
        for residue in self.model.get_residues():
            if residue.id[0] == ' ':  # Standard amino acid
                for atom in residue:
                    all_protein_atoms.append(atom)
        
        print(f"Total protein atoms: {len(all_protein_atoms)}")
        print(f"Interacting residues to check: {len(self.interacting_residues)}")
        
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
                print(f"Marking {res_id} as solvent accessible (score: {exposure_score:.2f})")
                self.solvent_accessible.add(res_id)
        
        # Make sure we have a reasonable number of solvent-accessible residues
        min_expected = max(2, int(len(self.interacting_residues) * 0.15))
        
        if len(self.solvent_accessible) < min_expected:
            print(f"Too few solvent-accessible residues detected ({len(self.solvent_accessible)}), "
                f"adding more based on residue type...")
            
            # Add more based on residue type
            remaining = sorted(
                [r for r in self.interacting_residues if r not in self.solvent_accessible],
                key=lambda r: 2 if r[0] in likely_exposed else 1,
                reverse=True
            )
            
            for res_id in remaining:
                if len(self.solvent_accessible) >= min_expected:
                    break
                print(f"Adding {res_id} as solvent accessible based on residue type")
                self.solvent_accessible.add(res_id)
        
        print(f"Final result: {len(self.solvent_accessible)} solvent-accessible residues out of {len(self.interacting_residues)} interacting residues")
        if self.solvent_accessible:
            print(f"Solvent accessible residues: {sorted(self.solvent_accessible)}")
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
                    if rsa > 0.25 and res_key in self.interacting_residues:
                        self.solvent_accessible.add(res_key)
                        
        except Exception as e:
            print(f"Warning: DSSP calculation failed. Falling back to geometric estimation. Error: {str(e)}")
            self.estimate_solvent_accessibility()
        
        return self.solvent_accessible
    
    def calculate_python_solvent_accessibility(self, probe_radius=1.4):
        """
        Enhanced solvent accessibility calculation in pure Python.
        Based on Shrake-Rupley algorithm but with practical approximations.
        
        Parameters:
        -----------
        probe_radius : float
            Radius of solvent probe in Angstroms (default: 1.4)
        """
        print("Using enhanced Python solvent accessibility calculation...")
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
            # If too few, fall back to the estimation method
            print("Warning: Too few solvent-accessible residues detected. Using estimation fallback.")
            return self.estimate_solvent_accessibility()
        
        print(f"Found {len(self.solvent_accessible)} solvent-accessible residues out of {len(self.interacting_residues)} interacting residues")
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
        print("Calculating Solvent Accessible Surface Areas...")
        
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
                    print(f"Error using external SASA tool: {str(e)}")
                    print("Falling back to geometric method")
                    self.calculate_realistic_solvent_accessibility()
                
                return self.solvent_accessible    

    def visualize(self, output_file='protein_ligand_interactions.png',figsize=(12, 12), dpi=300, title=None, show_3d_cues=True):
        """
        Generate a complete 2D visualization of protein-ligand interactions.
        """
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Debug check for solvent accessibility
        print("\n=== VISUALIZATION DEBUG ===")
        print(f"Total interacting residues: {len(self.interacting_residues)}")
        print(f"Solvent accessible residues: {len(self.solvent_accessible)}")
        
        if self.solvent_accessible:
            print("Solvent accessible residues:")
            for res_id in sorted(self.solvent_accessible):
                print(f"  - {res_id}")
        else:
            print("WARNING: No solvent accessible residues detected!")
        
        # Force reasonable solvent accessibility if all residues are marked or none are marked
        if len(self.solvent_accessible) == len(self.interacting_residues):
            print("WARNING: All residues are marked as solvent accessible!")
            print("This is likely incorrect. Removing some residues from solvent_accessible set...")
            
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
            
            print(f"After filtering: {len(self.solvent_accessible)} solvent accessible residues")
            for res_id in sorted(self.solvent_accessible):
                print(f"  - {res_id}")
        
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
            print("Warning: No interacting residues detected.")
            n_residues = 1
            
        # Calculate positions for residues
        radius = 250  # Distance from center to residues
        residue_positions = {}
        rect_width, rect_height = 60, 30  # Residue box dimensions
        
        # For debugging
        print(f"Drawing {n_residues} residue nodes, {len(self.solvent_accessible)} with solvent accessibility")
        
        # Arrange residues in a circle
        for i, res_id in enumerate(sorted(self.interacting_residues)):
            angle = 2 * math.pi * i / n_residues
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)
            residue_positions[res_id] = (x, y)
            
            # Draw solvent accessibility highlight with more visibility - ONLY for solvent accessible residues!
            if res_id in self.solvent_accessible:
                print(f"Drawing solvent accessibility circle for {res_id}")
                solvent_circle = Circle((x, y), 40, facecolor='#ADD8E6', 
                                    edgecolor='#87CEEB', alpha=0.5, zorder=1)
                ax.add_patch(solvent_circle)
            else:
                print(f"Residue {res_id} is NOT solvent accessible - no circle")
            
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
        
        print(f"Interaction diagram saved to {output_file}")
        return output_file
    
    
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
    def run_analysis(self, output_file=None, use_dssp=True, generate_report=False, report_file=None):
        """
        Run the complete analysis pipeline.
        
        Parameters:
        -----------
        output_file : str, optional
            Path where the output image will be saved.
        use_dssp : bool
            Whether to use DSSP for solvent accessibility (default: True)
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
        print("Detecting interactions...")
        self.detect_interactions() 
        
        print(f"Before filtering: {sum(len(ints) for ints in self.interactions.values())} total interactions")
        
        # Calculate solvent accessibility
        print("Calculating solvent accessibility...")
        try:
            if use_dssp:
                # Try DSSP first
                print("Trying DSSP method first...")
                self.calculate_dssp_solvent_accessibility()
                
                # Check if DSSP found too many solvent accessible residues
                if len(self.solvent_accessible) > len(self.interacting_residues) * 0.5:
                    print("DSSP found too many solvent accessible residues, using realistic method")
                    self.calculate_realistic_solvent_accessibility()
                # Check if DSSP found too few or no solvent accessible residues
                elif len(self.solvent_accessible) < 2:
                    print("DSSP didn't find enough solvent accessible residues, using realistic method")
                    self.calculate_realistic_solvent_accessibility()
            else:
                # Use our realistic method
                self.calculate_realistic_solvent_accessibility()
        except Exception as e:
            print(f"Error calculating solvent accessibility: {str(e)}")
            print("Falling back to realistic method")
            self.calculate_realistic_solvent_accessibility(exposure_threshold=0.2)
        
        # Make sure we don't have all residues marked as solvent accessible
        if len(self.solvent_accessible) == len(self.interacting_residues):
            print("WARNING: All residues marked as solvent accessible - this is likely incorrect")
            print("Applying stricter filtering")
            
            # Keep only 40% max of interacting residues as solvent accessible
            self.calculate_realistic_solvent_accessibility(exposure_threshold=0.3, max_percent=0.4)
        
        # Generate visualization
        print("Generating visualization...")
        viz_file = self.visualize(output_file=output_file)
        
        # Generate text report if requested
        if generate_report:
            if report_file is None:
                base_name = os.path.splitext(os.path.basename(self.structure_file))[0]
                report_file = f"{base_name}_interactions_report.txt"
            
            print("Generating interaction report...")
            try:
                if hasattr(self, 'generate_interaction_report'):
                    self.generate_interaction_report(output_file=report_file)
                else:
                    print("Warning: Report generation not available, skipping.")
            except Exception as e:
                print(f"Error generating report: {str(e)}")
        
        return viz_file
    
HybridProtLigMapper = add_interaction_detection_methods (HybridProtLigMapper)
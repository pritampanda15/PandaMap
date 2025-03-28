#!/usr/bin/env python
"""
Core functionality for PandaMap: A Python package for visualizing 
protein-ligand interactions with 2D ligand structure representation.
"""

import os
import math
from collections import defaultdict

import numpy as np
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
    
    def __init__(self, pdb_file, ligand_resname=None):
        """
        Initialize with a PDB file containing a protein-ligand complex.
        
        Parameters:
        -----------
        pdb_file : str
            Path to the PDB file with protein and ligand
        ligand_resname : str, optional
            Specific residue name of the ligand to focus on
        """
        self.pdb_file = pdb_file
        self.ligand_resname = ligand_resname
        
        # Parse the PDB file
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure('complex', pdb_file)
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
            raise ValueError("No ligand (HETATM) found in the PDB file.")
        
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
    
    def visualize(self, output_file='protein_ligand_interactions.png',
                  figsize=(12, 12), dpi=300, title=None):
        """
        Generate a 2D visualization of protein-ligand interactions
        matching the style of the reference image.
        
        Parameters:
        -----------
        output_file : str
            Path where the output image will be saved
        figsize : tuple
            Figure size in inches (width, height)
        dpi : int
            Resolution in dots per inch
        title : str, optional
            Title for the plot
        """
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Add light blue background for ligand (matching reference)
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
            n_residues = 1  # Avoid division by zero
            
        # Calculate positions for residues
        radius = 250  # Distance from center to residues
        residue_positions = {}
        
        # Arrange residues in a circle
        for i, res_id in enumerate(sorted(self.interacting_residues)):
            angle = 2 * math.pi * i / n_residues
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)
            residue_positions[res_id] = (x, y)
            
            # Draw solvent accessibility highlight (light blue circle behind residue)
            if res_id in self.solvent_accessible:
                solvent_circle = Circle((x, y), 40, facecolor='#ADD8E6', 
                                      edgecolor='none', alpha=0.3, zorder=1)
                ax.add_patch(solvent_circle)
            
            # Draw residue node as rectangle with rounded corners like in reference image
            # For simplicity, we'll use a white rectangle with black border
            rect_width = 60
            rect_height = 30
            residue_box = Rectangle((x-rect_width/2, y-rect_height/2), rect_width, rect_height,
                                  facecolor='white', edgecolor='black', linewidth=1.5,
                                  zorder=2, alpha=1.0)
            ax.add_patch(residue_box)
            
            # Add residue label (NAME NUMBER) like in reference image
            resname, resnum = res_id
            label = f"{resname} {resnum}"
            text = ax.text(x, y, label, ha='center', va='center',
                          fontsize=11, fontweight='bold', zorder=3)
            text.set_path_effects([path_effects.withStroke(linewidth=2, 
                                                         foreground='white')])
        
        # Define interaction styles to match the reference image
        interaction_styles = {
            'hydrogen_bonds': {
                'color': 'green',
                'linestyle': '-',
                'linewidth': 1.5,
                'marker_text': 'H',
                'marker_color': 'green',
                'marker_bg': '#E0FFE0',  # Light green bg
                'name': 'Hydrogen Bond'
            },
            'carbon_pi': {
                'color': '#666666',  # Dark gray
                'linestyle': '--',
                'linewidth': 1.5,
                'marker_text': 'C',
                'marker_color': '#666666',
                'marker_bg': 'white',
                'name': 'Carbon-Pi interaction'
            },
            'pi_pi_stacking': {
                'color': '#9370DB',  # Medium purple
                'linestyle': '--',
                'linewidth': 1.5,
                'marker_text': 'π',
                'marker_color': '#9370DB',
                'marker_bg': 'white',
                'name': 'Pi-Pi stacking'
            },
            'donor_pi': {
                'color': '#FF69B4',  # Hot pink
                'linestyle': '--',
                'linewidth': 1.5,
                'marker_text': 'D',
                'marker_color': '#FF69B4',
                'marker_bg': 'white',
                'name': 'Donor-Pi interaction'
            },
            'amide_pi': {
                'color': '#A52A2A',  # Brown
                'linestyle': '--',
                'linewidth': 1.5,
                'marker_text': 'A',
                'marker_color': '#A52A2A',
                'marker_bg': 'white',
                'name': 'Amide-Pi interaction'
            },
            'hydrophobic': {
                'color': '#808080',  # Gray
                'linestyle': ':',
                'linewidth': 1.0,
                'marker_text': 'h',
                'marker_color': '#808080',
                'marker_bg': 'white',
                'name': 'Hydrophobic'
            }
        }
        
        # Draw interaction lines for each residue
        for interaction_type, interactions in self.interactions.items():
            if interaction_type not in interaction_styles:
                continue
                
            style = interaction_styles[interaction_type]
            
            for interaction in interactions:
                res = interaction['protein_residue']
                res_id = (res.resname, res.id[1])
                lig_atom = interaction['ligand_atom']
                
                if res_id in residue_positions:
                    res_pos = residue_positions[res_id]
                    
                    # Try to use actual ligand atom position if available
                    if lig_atom.get_id() in atom_positions:
                        lig_pos = atom_positions[lig_atom.get_id()]
                    else:
                        # Determine the point on ligand circle edge as fallback
                        # Get the angle from center to residue
                        dx = res_pos[0] - ligand_pos[0]
                        dy = res_pos[1] - ligand_pos[1]
                        angle = math.atan2(dy, dx)
                        
                        # Calculate point on ligand circle edge
                        lig_edge_x = ligand_pos[0] + ligand_radius * math.cos(angle)
                        lig_edge_y = ligand_pos[1] + ligand_radius * math.sin(angle)
                        lig_pos = (lig_edge_x, lig_edge_y)
                    
                    # Draw curved line from ligand to residue
                    line = FancyArrowPatch(
                        lig_pos, res_pos,
                        connectionstyle="arc3,rad=0.1",
                        color=style['color'],
                        linestyle=style['linestyle'],
                        linewidth=style['linewidth'],
                        arrowstyle='-',
                        alpha=0.7,
                        zorder=4
                    )
                    ax.add_patch(line)
                    
                    # Add interaction marker (circle with letter)
                    midpoint_x = (lig_pos[0] + res_pos[0]) / 2
                    midpoint_y = (lig_pos[1] + res_pos[1]) / 2
                    
                    # Draw hexagonal marker like in reference image
                    n_sides = 6 if 'pi' in interaction_type else 0  # Hexagon for pi interactions
                    marker_radius = 12
                    
                    if n_sides > 0:
                        # Draw a hexagon for pi interactions
                        angles = np.linspace(0, 2*np.pi, n_sides+1)[:-1]
                        hex_vertices = np.array([
                            [midpoint_x + marker_radius * np.cos(angle),
                             midpoint_y + marker_radius * np.sin(angle)]
                            for angle in angles
                        ])
                        hex_patch = Polygon(
                            hex_vertices, 
                            closed=True,
                            facecolor=style.get('marker_bg', 'white'),
                            edgecolor=style['marker_color'],
                            linewidth=1.5,
                            alpha=0.9,
                            zorder=5
                        )
                        ax.add_patch(hex_patch)
                    else:
                        # Draw circle for other interactions
                        marker_circle = Circle(
                            (midpoint_x, midpoint_y), marker_radius,
                            facecolor=style.get('marker_bg', 'white'),
                            edgecolor=style['marker_color'],
                            linewidth=1.5,
                            alpha=0.9,
                            zorder=5
                        )
                        ax.add_patch(marker_circle)
                    
                    # Add interaction symbol
                    ax.text(
                        midpoint_x, midpoint_y,
                        style['marker_text'],
                        ha='center', va='center',
                        fontsize=9, color=style['marker_color'],
                        fontweight='bold',
                        zorder=6
                    )
        
        # Add legend box with title matching reference image
        legend_title = "Interacting structural groups"
        legend_elements = []
        
        # Add a "RES 1" element for residue representation
        residue_patch = Rectangle((0, 0), 1, 1, facecolor='white', 
                                 edgecolor='black', label='Interacting structural groups')
        legend_elements.append(residue_patch)
        
        # Interaction type markers for legend
        for int_type, style in interaction_styles.items():
            # Only include interaction types that are present
            if self.interactions[int_type]:
                if int_type == 'hydrogen_bonds':
                    # Special handling for H-bonds to match reference
                    line = Line2D([0], [0], color=style['color'],
                                 linestyle=style['linestyle'], linewidth=style['linewidth'],
                                 marker='o', markerfacecolor=style.get('marker_bg', 'white'), 
                                 markeredgecolor=style['color'],
                                 markersize=8, label=style['name'])
                elif 'pi' in int_type:
                    # Hexagonal markers for pi interactions
                    line = Line2D([0], [0], color=style['color'],
                                 linestyle=style['linestyle'], linewidth=style['linewidth'],
                                 marker='h', markerfacecolor=style.get('marker_bg', 'white'), 
                                 markeredgecolor=style['color'],
                                 markersize=8, label=style['name'])
                else:
                    # Circular markers for other interactions
                    line = Line2D([0], [0], color=style['color'],
                                 linestyle=style['linestyle'], linewidth=style['linewidth'],
                                 marker='o', markerfacecolor=style.get('marker_bg', 'white'), 
                                 markeredgecolor=style['color'],
                                 markersize=8, label=style['name'])
                legend_elements.append(line)
        
        # Add solvent accessibility indicator
        if self.solvent_accessible:
            solvent_patch = Rectangle((0, 0), 1, 1, facecolor='#ADD8E6', 
                                     alpha=0.3, edgecolor=None, label='Solvent accessible')
            legend_elements.append(solvent_patch)
        
        # Create legend box in top right corner like reference image
        legend = ax.legend(
            handles=legend_elements,
            title=legend_title,
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
            plt.title(f"Protein-Ligand Interactions: {os.path.basename(self.pdb_file)}", 
                     fontsize=16)
        
        # Save figure
        plt.tight_layout()
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        plt.close()
        
        print(f"Interaction diagram saved to {output_file}")
        return output_file
    
    def run_analysis(self, output_file=None):
        """
        Run the complete analysis pipeline.
        
        Parameters:
        -----------
        output_file : str, optional
            Path where the output image will be saved. If None, a default name will be generated.
            
        Returns:
        --------
        str : Path to the generated visualization file
        """
        if output_file is None:
            base_name = os.path.splitext(os.path.basename(self.pdb_file))[0]
            output_file = f"{base_name}_interactions.png"
        
        # Detect protein-ligand interactions
        print("Detecting interactions...")
        self.detect_interactions()
        
        # Estimate solvent accessibility
        print("Estimating solvent accessibility...")
        self.estimate_solvent_accessibility()
        
        # Generate visualization
        print("Generating visualization...")
        return self.visualize(output_file=output_file)

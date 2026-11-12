"""
Ligand drawing used by the classic (circular) PandaMap diagram.

The modern diagram and the report figure use ``pandamap.ligand2d`` instead.
"""
import logging
import math

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.colors import to_rgb
from matplotlib.patches import Circle, FancyArrowPatch, Wedge, Rectangle, Polygon
from matplotlib.lines import Line2D
import matplotlib.patheffects as path_effects

logger = logging.getLogger("pandamap")


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

        # ── Attempt RDKit coordinate generation first ──────────────────────
        rdkit_coords = self._try_rdkit_2d_coords()
        if rdkit_coords:
            return rdkit_coords
        # ── Fallback: PCA-based projection ─────────────────────────────────

        # Simple projection onto the xy-plane
        coords_2d = {}
        
        # Get all 3D coordinates and find center
        all_coords = np.array([info['coord'] for info in self.atom_coords.values()])
        center = np.mean(all_coords, axis=0)
        
        # Subtract center to center the molecule
        centered_coords = all_coords - center
        
        # Direct fallback for small molecules (less than 3 atoms)
        if len(centered_coords) < 3:
            logger.warning("Warning: Not enough atoms for PCA projection. Using simple XY projection.")
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
            logger.warning(f"Warning: Could not compute optimal projection. Using simple XY projection. Error: {str(e)}")
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
    
    def detect_rings(self, bond_threshold=2.0):
        """
        Identify ligand atoms that belong to ring systems using iterative
        leaf-node pruning on the bond graph (works for any ring size).

        Returns
        -------
        set
            Atom IDs of ring atoms (atoms that are part of at least one cycle).
        """
        atom_ids = list(self.atom_coords.keys())
        bonds = self.find_bonds(distance_threshold=bond_threshold)

        # Build adjacency list
        adj = {a: [] for a in atom_ids}
        for a, b in bonds:
            adj[a].append(b)
            adj[b].append(a)

        # Degree map (ignore H for ring detection)
        degree = {a: sum(1 for nb in adj[a]
                         if self.atom_coords[nb]['element'] not in ('H', ''))
                  for a in atom_ids
                  if self.atom_coords[a]['element'] not in ('H', '')}

        # Iteratively prune leaf nodes (degree == 1)
        changed = True
        while changed:
            changed = False
            for a in list(degree.keys()):
                if degree[a] == 1:
                    degree[a] = 0
                    for nb in adj[a]:
                        if nb in degree and degree[nb] > 0:
                            degree[nb] -= 1
                    changed = True

        # Atoms still with degree >= 2 are in rings
        return {a for a, d in degree.items() if d >= 2}

    def _try_rdkit_2d_coords(self):
        """
        Attempt to generate chemically accurate 2D coordinates using RDKit.
        Builds an RDKit molecule from bond topology (single bonds assumed) and
        calls AllChem.Compute2DCoords().  Returns None if RDKit is unavailable
        or if coordinate generation fails.

        Returns
        -------
        dict or None
            Mapping atom_id → np.ndarray([x, y]) scaled to match PCA output,
            or None on failure.
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem

            atom_ids = list(self.atom_coords.keys())
            bonds = self.find_bonds()

            mol = Chem.RWMol()
            idx_map = {}  # atom_id → RDKit atom index

            for atom_id in atom_ids:
                element = (self.atom_coords[atom_id]['element'] or 'C').capitalize()
                try:
                    rdatom = Chem.Atom(element)
                except Exception:
                    rdatom = Chem.Atom('C')
                idx = mol.AddAtom(rdatom)
                idx_map[atom_id] = idx

            for a1, a2 in bonds:
                if a1 in idx_map and a2 in idx_map:
                    try:
                        mol.AddBond(idx_map[a1], idx_map[a2], Chem.BondType.SINGLE)
                    except Exception:
                        pass

            mol = mol.GetMol()
            AllChem.Compute2DCoords(mol)
            conf = mol.GetConformer()

            coords_2d = {}
            for atom_id, idx in idx_map.items():
                pos = conf.GetAtomPosition(idx)
                coords_2d[atom_id] = np.array([pos.x, pos.y]) * 30.0  # scale to ~PCA units
            return coords_2d

        except Exception:
            return None

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
            logger.warning("Warning: Could not generate ligand coordinates. Drawing placeholder.")
            circle = Circle(center, radius/2, fill=False, edgecolor='black', linestyle='-')
            ax.add_patch(circle)
            ax.text(center[0], center[1], "Ligand", ha='center', va='center')
            return {}
            
        # Find bonds
        bonds = self.find_bonds()
        
        # Scale coordinates to fit within the specified radius.
        #
        # Previously this used max(|coord|) over both axes of *un-centred*
        # coordinates, so an off-centre layout (RDKit emits these routinely)
        # was scaled by its distance from the origin rather than its true
        # size, shrinking the depiction.  Centre on the bounding box first,
        # then scale by the larger half-range.
        all_coords = np.array(list(coords_2d.values()), dtype=float)
        bbox_centre = (all_coords.max(axis=0) + all_coords.min(axis=0)) / 2.0
        all_coords = all_coords - bbox_centre
        coords_2d = {aid: np.asarray(c, dtype=float) - bbox_centre
                     for aid, c in coords_2d.items()}
        max_extent = np.max(np.abs(all_coords))
        if max_extent <= 0:
            max_extent = 1.0
        scaling_factor = radius / (max_extent * 1.2)  # Leave some margin

        # Atom marker radius must track the drawn bond length, otherwise a
        # ligand with many atoms is scaled down until neighbouring circles
        # overlap and the depiction collapses into an unreadable blob.
        bond_lengths = [
            np.linalg.norm((coords_2d[a] - coords_2d[b]) * scaling_factor)
            for a, b in bonds if a in coords_2d and b in coords_2d
        ]
        if bond_lengths:
            median_bond = float(np.median(bond_lengths))
        else:
            median_bond = radius / 4.0
        # 0.30 x bond length keeps a clear gap between bonded atoms.
        atom_unit = float(np.clip(median_bond * 0.30, 2.0, 11.0))
        
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
                base_size = atom_unit if element in ['C', 'H'] else atom_unit * 1.22
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
            base_size = atom_unit if element in ['C', 'H'] else atom_unit * 1.22
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

"""
2D depiction of the ligand for the interaction diagrams.

``layout_ligand`` turns the ligand's 3D atoms into what a skeletal formula
needs: 2D coordinates without overlapping atoms, the bonds, which bonds are
double, and which rings are aromatic.

Bond orders and aromaticity are *perceived from the 3D geometry* (bond
lengths, ring planarity, number of neighbours), because structure files do
not store them.  The rules are deliberately conservative: a bond is drawn
double only when it is clearly shorter than a single bond between the same
elements.  Treat the result as a readable depiction, not as an authoritative
assignment of bond orders.
"""
from collections import defaultdict, deque

import numpy as np

_COVALENT_RADII = {'H': 0.31, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57,
                   'SI': 1.11, 'P': 1.07, 'S': 1.05, 'CL': 1.02, 'SE': 1.20,
                   'BR': 1.20, 'I': 1.39}
_METALS = {'LI', 'NA', 'K', 'MG', 'CA', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'MO', 'W',
           'PT', 'PD', 'RU', 'RH', 'IR', 'CD', 'HG', 'V', 'CR', 'AL', 'GA'}
_AROMATIC_ELEMENTS = {'C', 'N', 'O', 'S'}

# A bond shorter than this (Angstrom) between the two elements is drawn double.
# Single-bond lengths for comparison: C-C 1.54, C-N 1.47, C-O 1.43, amide C-N 1.33.
_DOUBLE_BOND_MAX = {
    frozenset(('C', 'C')): 1.38,
    frozenset(('C', 'N')): 1.31,
    frozenset(('C', 'O')): 1.28,
    frozenset(('N', 'N')): 1.28,
    frozenset(('N', 'O')): 1.28,
    frozenset(('C', 'S')): 1.72,
    frozenset(('P', 'O')): 1.52,
    frozenset(('S', 'O')): 1.50,
}
# how many double bonds an atom of this element may carry in the drawing
_MAX_DOUBLE = {'S': 2, 'P': 1}


def element_of(atom):
    return (atom.element or atom.get_name()[:1]).strip().upper()


def find_bonds(atoms):
    """Index pairs of covalently bonded atoms, from 3D distances."""
    bonds = []
    for i, a in enumerate(atoms):
        ra = _COVALENT_RADII.get(element_of(a), 1.3)
        for j in range(i + 1, len(atoms)):
            b = atoms[j]
            limit = ra + _COVALENT_RADII.get(element_of(b), 1.3) + 0.45
            if np.linalg.norm(a.coord - b.coord) < limit:
                bonds.append((i, j))
    return bonds


def smallest_ring(start, neighbours, max_size=8):
    """Atom indices of the smallest ring through ``start`` (empty if none)."""
    best = []
    for first in neighbours[start]:
        # shortest path first -> start that does not use the direct bond
        parent = {first: None}
        queue = deque([first])
        while queue:
            node = queue.popleft()
            if node == start:
                break
            for nxt in neighbours[node]:
                if nxt in parent or (node == first and nxt == start):
                    continue
                parent[nxt] = node
                queue.append(nxt)
        if start not in parent:
            continue
        ring, node = [], start
        while node is not None:
            ring.append(node)
            node = parent[node]
        if 3 <= len(ring) <= max_size and (not best or len(ring) < len(best)):
            best = ring
    return best


def find_rings(n_atoms, neighbours, max_size=8):
    """The distinct smallest rings of the molecule, each as an ordered index list."""
    rings = {}
    for i in range(n_atoms):
        ring = smallest_ring(i, neighbours, max_size)
        if ring:
            rings.setdefault(frozenset(ring), ring)
    return list(rings.values())


def _is_aromatic(ring, atoms, xyz, neighbours):
    """Planar 5- or 6-membered ring of sp2-like C/N/O/S atoms with short bonds."""
    if len(ring) not in (5, 6):
        return False
    if any(element_of(atoms[i]) not in _AROMATIC_ELEMENTS for i in ring):
        return False
    if any(len(neighbours[i]) > 3 for i in ring):          # sp3 centre in the ring
        return False
    pts = xyz[ring] - xyz[ring].mean(axis=0)
    normal = np.linalg.svd(pts)[2][-1]
    if np.abs(pts @ normal).max() > 0.20:                  # not planar
        return False
    lengths = [np.linalg.norm(xyz[ring[k]] - xyz[ring[(k + 1) % len(ring)]])
               for k in range(len(ring))]
    return float(np.mean(lengths)) < 1.45


def perceive_bond_orders(atoms, bonds, neighbours, aromatic_rings):
    """Set of bonds (i, j) to draw as double, perceived from bond lengths."""
    xyz = np.array([a.coord for a in atoms], dtype=float)
    aromatic_bonds = set()
    for ring in aromatic_rings:
        for k in range(len(ring)):
            aromatic_bonds.add(frozenset((ring[k], ring[(k + 1) % len(ring)])))

    candidates = []
    for i, j in bonds:
        if frozenset((i, j)) in aromatic_bonds:
            continue
        ei, ej = element_of(atoms[i]), element_of(atoms[j])
        limit = _DOUBLE_BOND_MAX.get(frozenset((ei, ej)))
        if limit is None:
            continue
        # an atom with four neighbours has no room for a double bond (S and P excepted)
        if any(len(neighbours[k]) > 3 and element_of(atoms[k]) not in ('S', 'P') for k in (i, j)):
            continue
        length = float(np.linalg.norm(xyz[i] - xyz[j]))
        if length < limit:
            candidates.append((length - limit, i, j))      # most clearly short first

    double, count = set(), defaultdict(int)
    for _, i, j in sorted(candidates):
        if all(count[k] < _MAX_DOUBLE.get(element_of(atoms[k]), 1) for k in (i, j)):
            double.add((i, j))
            count[i] += 1
            count[j] += 1
    return double


def _project(xyz):
    """Orthographic projection onto the plane of largest spread."""
    centred = xyz - xyz.mean(axis=0)
    if len(xyz) < 3:
        return centred[:, :2].copy()
    _, vectors = np.linalg.eigh(np.cov(centred.T))
    return centred @ vectors[:, ::-1][:, :2]


def _rdkit_coords(atoms, bonds, double):
    """2D coordinates from RDKit, or None when it is unavailable or fails."""
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import AllChem
        RDLogger.DisableLog('rdApp.*')
        mol = Chem.RWMol()
        for atom in atoms:
            try:
                mol.AddAtom(Chem.Atom(element_of(atom).capitalize()))
            except Exception:
                mol.AddAtom(Chem.Atom('C'))
        for i, j in bonds:
            order = Chem.BondType.DOUBLE if (i, j) in double else Chem.BondType.SINGLE
            mol.AddBond(i, j, order)
        mol = mol.GetMol()
        mol.UpdatePropertyCache(strict=False)
        AllChem.Compute2DCoords(mol)
        conf = mol.GetConformer()
        return np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y]
                         for i in range(len(atoms))], dtype=float)
    except Exception:
        return None


def _layout_quality(pts, bonds):
    """(median bond length, shortest non-bonded distance / median bond,
    coefficient of variation of the bond lengths)."""
    if not bonds:
        return 1.0, 9.9, 0.0
    lengths = np.array([np.linalg.norm(pts[i] - pts[j]) for i, j in bonds])
    median = float(np.median(lengths)) or 1.0
    bonded = {frozenset(b) for b in bonds}
    closest = 9.9
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if frozenset((i, j)) not in bonded:
                closest = min(closest, float(np.linalg.norm(pts[i] - pts[j])) / median)
    return median, closest, float(lengths.std() / median)


def relax_layout(pts, bonds, steps=400):
    """
    Tidy a projected layout: pull bonds towards one common length and push
    apart atoms that the projection placed on top of each other.  Ring shapes
    and the overall arrangement of the projection are preserved.
    """
    pts = np.array(pts, dtype=float)
    n = len(pts)
    if n < 3 or not bonds:
        return pts
    median = _layout_quality(pts, bonds)[0]
    pts = pts / median                                     # unit bond length
    bonded = {frozenset(b) for b in bonds}
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)
             if frozenset((i, j)) not in bonded]
    rng = np.random.RandomState(0)
    for step in range(steps):
        move = np.zeros_like(pts)
        for i, j in bonds:
            d = pts[j] - pts[i]
            length = np.linalg.norm(d)
            if length < 1e-6:
                d, length = rng.rand(2) - 0.5, 1e-3
            shift = 0.25 * (length - 1.0) * d / length
            move[i] += shift
            move[j] -= shift
        for i, j in pairs:
            d = pts[j] - pts[i]
            length = np.linalg.norm(d)
            if length < 1.15:
                if length < 1e-6:
                    d, length = rng.rand(2) - 0.5, 1e-3
                shift = 0.25 * (1.15 - length) * d / length
                move[i] -= shift
                move[j] += shift
        size = np.linalg.norm(move, axis=1, keepdims=True)
        move = np.where(size > 0.12, move * 0.12 / np.maximum(size, 1e-9), move)
        pts += move
        if size.max() < 1e-3:
            break
    return pts - pts.mean(axis=0)


def layout_ligand(atoms):
    """
    Lay a ligand out in 2D.

    Parameters
    ----------
    atoms : list of Bio.PDB atoms (heavy atoms of the ligand)

    Returns
    -------
    dict with
        ``coords``    N x 2 array, centred, median bond length 1
        ``bonds``     list of (i, j)
        ``double``    set of (i, j) drawn as double bonds
        ``aromatic``  list of rings (index lists) drawn with an inner circle
        ``rings``     all rings found
        ``neighbours`` dict index -> list of bonded indices
        ``method``    'rdkit', 'projection' or 'projection+relaxed'
    """
    xyz = np.array([a.coord for a in atoms], dtype=float)
    bonds = find_bonds(atoms)
    neighbours = defaultdict(list)
    for i, j in bonds:
        neighbours[i].append(j)
        neighbours[j].append(i)
    rings = find_rings(len(atoms), neighbours)
    aromatic = [r for r in rings if _is_aromatic(r, atoms, xyz, neighbours)]
    double = perceive_bond_orders(atoms, bonds, neighbours, aromatic)

    pts, method = None, 'projection'
    has_metal = any(element_of(a) in _METALS for a in atoms)
    if not has_metal and len(atoms) >= 3:
        # RDKit lays out organic molecules well, but distorts metal complexes
        candidate = _rdkit_coords(atoms, bonds, double)
        if candidate is not None and _layout_quality(candidate, bonds)[1] > 0.45:
            pts, method = candidate, 'rdkit'
    if pts is None:
        pts = _project(xyz)
        _, closest, spread = _layout_quality(pts, bonds)
        if closest < 0.6 or spread > 0.22:
            pts, method = relax_layout(pts, bonds), 'projection+relaxed'

    pts = pts - pts.mean(axis=0)
    median = _layout_quality(pts, bonds)[0]
    return {'coords': pts / median, 'bonds': bonds, 'double': double,
            'aromatic': aromatic, 'rings': rings, 'neighbours': neighbours,
            'method': method}

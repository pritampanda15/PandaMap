"""
Modern 2D interaction diagram for PandaMap.

Layout rules
------------
* The ligand is drawn as a skeletal formula in the centre.
* Every residue is placed on the side of the ligand it actually contacts, so
  interaction lines are short and rarely cross.
* Residue boxes are coloured by residue type; a blue outline marks solvent
  accessible residues.
* Contacts from one residue to the same ligand atom/ring share one line that
  carries one badge per interaction type.
* The legend sits below the plot and lists only what is present.
"""
import logging
import math
import os
from collections import defaultdict, deque

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyBboxPatch

from .ligand2d import element_of, find_bonds, layout_ligand, smallest_ring

logger = logging.getLogger("pandamap")

# interaction type -> (legend label, colour, line style, badge text), in priority order
INTERACTION_STYLE = {
    'covalent': ('Covalent', '#212529', '-', 'cov'),
    'metal_coordination': ('Metal coordination', '#b8860b', '-', 'M'),
    'hydrogen_bonds': ('Hydrogen bond', '#1a8f3c', '-', 'H'),
    'water_bridges': ('Water bridge', '#0b7285', '-', 'W'),
    'salt_bridge': ('Salt bridge', '#d9480f', '-', 'S'),
    'ionic': ('Ionic', '#e8590c', '-', 'I'),
    'attractive_charge': ('Attractive charge', '#1c7ed6', '-', '+'),
    'halogen_bonds': ('Halogen bond', '#0c8599', '-', 'X'),
    'pi_pi_stacking': ('π–π stacking', '#7048e8', '--', 'π'),
    'pi_cation': ('π–cation', '#c2255c', '--', 'π+'),
    'cation_pi': ('Cation–π', '#c2255c', '--', 'π+'),
    'donor_pi': ('Donor–π', '#e64980', '--', 'D'),
    'amide_pi': ('Amide–π', '#a0522d', '--', 'A'),
    'carbon_pi': ('Carbon–π', '#868e96', '--', 'C'),
    'alkyl_pi': ('Alkyl–π', '#4c7fb0', '--', 'Aπ'),
    'hydrophobic': ('Hydrophobic', '#868e96', ':', 'h'),
    'repulsion': ('Repulsion', '#e03131', ':', '!'),
}
_PRIORITY = {t: i for i, t in enumerate(INTERACTION_STYLE)}
# interactions that involve a ligand ring: their line points at the ring centre
_RING_TYPES = {'pi_pi_stacking', 'pi_cation', 'donor_pi', 'amide_pi', 'carbon_pi', 'alkyl_pi'}

RESIDUE_CLASSES = (
    ('Hydrophobic residue', '#e6f4d7', 'ALA VAL LEU ILE MET PHE TRP PRO'),
    ('Polar residue', '#e3eefb', 'SER THR ASN GLN TYR CYS GLY'),
    ('Positive residue', '#dcdcff', 'LYS ARG HIS'),
    ('Negative residue', '#ffdede', 'ASP GLU'),
)
_RESIDUE_CLASS = {name: (label, colour)
                  for label, colour, names in RESIDUE_CLASSES for name in names.split()}
_OTHER_RESIDUE = ('Other residue', '#f1f3f5')

ATOM_COLOURS = {'N': '#2845c7', 'O': '#d62828', 'S': '#b8930a', 'P': '#e8590c',
                'F': '#2b8a3e', 'CL': '#2b8a3e', 'BR': '#8b3a3a', 'I': '#7b2cbf'}
_SOLVENT_EDGE = '#3b9ad9'
_BOX_W, _BOX_H = 2.25, 0.86          # residue box, data units
_R_LIGAND = 3.4                      # half-extent of the ligand drawing
_R_RESIDUES = 8.6                    # minimum radius of the residue ring
_MERGE_DISTANCE = 2.3                # targets closer than this share one line
_MERGE_ANGLE = 13                    # ... or lines closer than this many degrees

# shared with the 3D contact view
_element, _ligand_bonds, _smallest_ring = element_of, find_bonds, smallest_ring


def _spread_angles(preferred, min_gap):
    """Move angles apart until neighbours are at least ``min_gap`` apart,
    staying as close as possible to the preferred directions.
    Returns (order, angles) with ``order`` indexing into ``preferred``."""
    n = len(preferred)
    order = sorted(range(n), key=lambda i: preferred[i])
    angles = np.array([preferred[i] for i in order], dtype=float)
    if n < 2:
        return order, angles
    min_gap = min(min_gap, 2 * math.pi / n * 0.97)
    for _ in range(1000):
        moved = False
        for i in range(n):
            j = (i + 1) % n
            gap = angles[j] - angles[i] + (2 * math.pi if j == 0 else 0.0)
            if gap < min_gap - 1e-9:
                shift = (min_gap - gap) / 2
                angles[i] -= shift
                angles[j] += shift
                moved = True
        if not moved:
            break
    return order, angles


def _box_edge(centre, target, pad=0.08):
    """Point where the segment centre->target leaves the residue box."""
    d = np.asarray(target, dtype=float) - centre
    norm = np.linalg.norm(d)
    if norm == 0:
        return np.array(centre, dtype=float)
    d /= norm
    tx = (_BOX_W / 2 + pad) / abs(d[0]) if abs(d[0]) > 1e-9 else float('inf')
    ty = (_BOX_H / 2 + pad) / abs(d[1]) if abs(d[1]) > 1e-9 else float('inf')
    return centre + d * min(tx, ty)


def draw_interaction_diagram(mapper, output_file='protein_ligand_interactions.png',
                             figsize=(11, 11.6), dpi=300, title=None,
                             show_distances=True, show_nearby_residues=False,
                             ax=None, theme=None):
    """
    Draw the 2D protein-ligand interaction diagram for an analysed mapper.

    Parameters
    ----------
    mapper : HybridProtLigMapper
        A mapper on which ``detect_interactions()`` has been run.
    output_file : str
        Image path; the extension selects the format (.png, .svg, .pdf, ...).
    figsize, dpi, title :
        Figure size in inches, resolution, and an optional custom title.
    show_distances : bool
        Print the donor-acceptor distance on hydrogen bonds.
    show_nearby_residues : bool
        Also show binding-site residues that have no specific interaction
        (drawn with a dashed outline and no line).

    ax, theme :
        For embedding the diagram in a larger figure: draw on ``ax`` instead
        of creating and saving a figure (no title or legend is added), with
        optional ``theme`` overrides - ``scale`` and ``min_font`` for text,
        ``colours`` (interaction type -> colour), ``ink``, ``box_fill``,
        ``box_edge``, ``solvent_edge``, ``weight`` and ``atom_colours``.
        A dict describing what was drawn is returned instead of a path.

    Returns the path written.
    """
    # ---- ligand layout ---------------------------------------------------
    all_atoms = list(mapper.ligand_atoms)
    heavy = [a for a in all_atoms if _element(a) != 'H'] or all_atoms
    layout = layout_ligand(heavy)
    pts = np.array(layout['coords'], dtype=float)
    extent = np.abs(pts).max()
    unit = _R_LIGAND / extent if extent > 0 else 1.0     # drawing units per bond length
    unit = min(unit, 1.15)                               # small ligands are not blown up
    pts *= unit
    index = {id(a): i for i, a in enumerate(heavy)}
    heavy_xyz = np.array([a.coord for a in heavy])

    def atom_index(atom):
        """Index of the drawn atom for a ligand atom (hydrogens map to the
        heavy atom they are attached to)."""
        i = index.get(id(atom))
        if i is None:
            i = int(np.argmin(np.linalg.norm(heavy_xyz - atom.coord, axis=1)))
        return i

    bonds, neighbours = layout['bonds'], layout['neighbours']

    ligand_ids = {id(a) for a in all_atoms}
    partners = {}                     # non-ligand partner atoms (metal ions)

    def target_point(itype, contact):
        """(key, xy) of the point on the ligand a contact line should reach."""
        atom = contact['ligand_atom']
        i = atom_index(atom)
        if id(atom) not in ligand_ids:
            # e.g. a metal ion bridging protein and ligand: draw it as its
            # own node just outside the ligand atom it sits next to
            direction = pts[i] / (np.linalg.norm(pts[i]) or 1.0)
            if not direction.any():
                direction = np.array([0.0, 1.0])
            xy = pts[i] + direction * 1.5
            partners[id(atom)] = (_element(atom).capitalize(), xy)
            return ('partner', id(atom)), xy
        if itype in _RING_TYPES:
            ring = _smallest_ring(i, neighbours)
            if ring:
                return ('ring', tuple(sorted(ring))), pts[ring].mean(axis=0)
        return ('atom', i), pts[i]

    # ---- collect contacts per residue ------------------------------------
    contacts = defaultdict(lambda: defaultdict(list))   # residue -> target -> [(type, contact)]
    targets = {}
    for itype, items in mapper.interactions.items():
        for contact in items:
            res = contact.get('protein_residue')
            if res is None or contact.get('ligand_atom') is None:
                continue
            key = (res.resname.strip(), res.id[1], res.get_parent().id)
            tkey, xy = target_point(itype, contact)
            targets[tkey] = xy
            contacts[key][tkey].append((itype, contact))

    residues = list(contacts)
    preferred = []
    for key in residues:
        v = np.mean([targets[t] for t in contacts[key]], axis=0)
        preferred.append(math.atan2(v[1], v[0]))

    nearby = []
    if show_nearby_residues:
        shown = {(k[0], k[1]) for k in residues}
        for res_id in sorted(mapper.interacting_residues):
            residue = mapper.protein_residues.get(res_id)
            if res_id in shown or residue is None:
                continue
            xyz = np.array([a.coord for a in residue])
            dist = np.linalg.norm(xyz[:, None, :] - heavy_xyz[None, :, :], axis=2)
            v = pts[int(np.unravel_index(dist.argmin(), dist.shape)[1])]
            nearby.append((res_id[0].strip(), res_id[1], residue.get_parent().id))
            preferred.append(math.atan2(v[1], v[0]))
    everyone = residues + nearby

    # ---- place residues on a ring, near the atoms they touch --------------
    n = len(everyone)
    radius = max(_R_RESIDUES, n * (_BOX_W + 0.45) / (2 * math.pi))
    limit = radius + 1.7
    font = (_R_RESIDUES + 1.7) / limit            # keep text in step with the boxes
    theme = dict(theme or {})
    font *= theme.get('scale', 1.0)
    fs = lambda size: max(theme.get('min_font', 0), size * font)
    weight = theme.get('weight', 'bold')
    ink = theme.get('ink', '#343a40')
    atom_colours = theme.get('atom_colours', ATOM_COLOURS)

    def colour_of(itype):
        default = INTERACTION_STYLE.get(itype, (itype, '#495057', '-', '?'))[1]
        return (theme.get('colours') or {}).get(itype, default)
    order, angles = _spread_angles(preferred, math.radians(18))
    centre = {everyone[i]: np.array([math.cos(a), math.sin(a)]) * radius
              for i, a in zip(order, angles)}

    # Contacts of one residue whose lines would run on top of each other -
    # targets close together on the ligand, or seen from the residue under a
    # small angle - share a single line to the highest-priority target.
    def _overlapping(c, a, b):
        if np.linalg.norm(targets[a] - targets[b]) < _MERGE_DISTANCE:
            return True
        u, v = targets[a] - c, targets[b] - c
        cos = np.dot(u, v) / ((np.linalg.norm(u) * np.linalg.norm(v)) or 1.0)
        return cos > math.cos(math.radians(_MERGE_ANGLE))

    for key, by_target in contacts.items():
        ranked = sorted(by_target, key=lambda t: min(_PRIORITY.get(i, 99)
                                                     for i, _ in by_target[t]))
        merged = {}
        for tkey in ranked:
            home = next((m for m in merged
                         if m[0] != 'partner' and tkey[0] != 'partner'
                         and _overlapping(centre[key], m, tkey)), None)
            merged.setdefault(home or tkey, []).extend(by_target[tkey])
        contacts[key] = merged

    own_figure = ax is None
    if own_figure:
        fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect('equal')
    ax.axis('off')

    # ---- ligand ------------------------------------------------------------
    ring_centres = [pts[ring].mean(axis=0) for ring in layout['rings']]
    for i, j in bonds:
        ax.plot([pts[i][0], pts[j][0]], [pts[i][1], pts[j][1]], color=ink,
                lw=2.6 * font, solid_capstyle='round', zorder=3)
        if (i, j) in layout['double'] or (j, i) in layout['double']:
            # second line of a double bond: shorter, on the side of the ring
            # the bond belongs to (or of the rest of the molecule)
            d = pts[j] - pts[i]
            normal = np.array([-d[1], d[0]]) / (np.linalg.norm(d) or 1.0)
            mid = (pts[i] + pts[j]) / 2
            inside = [c for ring, c in zip(layout['rings'], ring_centres)
                      if i in ring and j in ring]
            towards = (inside[0] if inside else pts.mean(axis=0)) - mid
            if np.dot(normal, towards) < 0:
                normal = -normal
            p, q = pts[i] + d * 0.14 + normal * 0.17 * unit, pts[j] - d * 0.14 + normal * 0.17 * unit
            ax.plot([p[0], q[0]], [p[1], q[1]], color=ink, lw=1.9 * font,
                    solid_capstyle='round', zorder=3)
    for ring in layout['aromatic']:
        centre_xy = pts[ring].mean(axis=0)
        radius = np.linalg.norm(pts[ring] - centre_xy, axis=1).mean()
        ax.add_patch(Circle(centre_xy, radius * 0.58, fill=False, ec=ink,
                            lw=1.5 * font, zorder=3))
    for i, atom in enumerate(heavy):
        element = _element(atom)
        if element != 'C' or not neighbours[i]:
            ax.text(pts[i][0], pts[i][1], element.capitalize(), ha='center', va='center',
                    fontsize=fs(13), fontweight='bold', zorder=5,
                    color=atom_colours.get(element, ink),
                    bbox=dict(boxstyle='circle,pad=0.12', fc='white', ec='none'))

    for symbol, xy in partners.values():
        ax.text(xy[0], xy[1], symbol, ha='center', va='center', fontsize=fs(11),
                fontweight='bold', color='#7a5c00', zorder=9,   # above lines and labels
                bbox=dict(boxstyle='circle,pad=0.25', fc='#fff3bf', ec='#b8860b',
                          lw=1.5 * font))

    # ---- residues and interaction lines ------------------------------------
    multi_chain = len({k[2] for k in everyone}) > 1
    solvent = set(getattr(mapper, 'solvent_accessible', ()) or ())
    used_types, used_classes, any_solvent = set(), [], False

    for key in everyone:
        resname, resnum, chain = key
        c = centre[key]
        label, fill = _RESIDUE_CLASS.get(resname, _OTHER_RESIDUE)
        if (label, fill) not in used_classes:
            used_classes.append((label, fill))
        accessible = (resname, resnum) in solvent
        any_solvent = any_solvent or accessible
        ax.add_patch(FancyBboxPatch(
            (c[0] - _BOX_W / 2, c[1] - _BOX_H / 2), _BOX_W, _BOX_H,
            boxstyle='round,pad=0.02,rounding_size=0.22', zorder=6,
            fc=theme.get('box_fill', fill),
            ec=theme.get('solvent_edge', _SOLVENT_EDGE) if accessible
            else theme.get('box_edge', '#495057'),
            lw=(3.2 if accessible else 1.3) * font,
            ls='-' if key in contacts else (0, (3, 2))))
        text = f"{resname} {resnum}" + (f" {chain}" if multi_chain else "")
        ax.text(c[0], c[1], text, ha='center', va='center',
                fontsize=fs(10.8 if multi_chain else 12.5),
                fontweight=weight, color=theme.get('ink', '#212529'), zorder=7)

        groups = sorted(contacts.get(key, {}).items(),
                        key=lambda kv: min(_PRIORITY.get(t, 99) for t, _ in kv[1]))
        # badges of one residue step inward one after another, across all of
        # its lines, so they never sit on top of each other
        n_badges = sum(len({t for t, _ in members}) for _, members in groups)
        shortest = min((np.linalg.norm(np.asarray(targets[t]) - c) for t, _ in groups),
                       default=0.0) - _BOX_W / 2
        step = min(theme.get('badge_step', 0.62), max(0.0, shortest - 1.2) / max(1, n_badges))
        first = theme.get('badge_start', 0.55)
        slot = 0
        for g, (tkey, members) in enumerate(groups):
            by_type = defaultdict(list)
            for itype, contact in members:
                by_type[itype].append(contact)
            types = sorted(by_type, key=lambda t: _PRIORITY.get(t, 99))
            used_types.update(types)
            _, _, linestyle, _ = INTERACTION_STYLE.get(
                types[0], (types[0], '#495057', '-', '?'))
            colour = colour_of(types[0])

            end = np.asarray(targets[tkey], dtype=float)
            start = _box_edge(c, end)
            direction = end - start
            length = np.linalg.norm(direction)
            if length < 1e-6:
                continue
            direction /= length
            # several lines from one residue leave the box side by side
            start = start + np.array([-direction[1], direction[0]]) * 0.22 * (
                g - (len(groups) - 1) / 2)
            end = end - direction * {'atom': 0.30, 'partner': 0.42}.get(tkey[0], 0.0)
            length = np.linalg.norm(end - start)
            ax.plot([start[0], end[0]], [start[1], end[1]], color=colour, ls=linestyle,
                    lw=(2.0 if linestyle == '-' else 1.6) * font, alpha=0.9, zorder=2)
            if tkey[0] == 'ring':
                ax.plot(end[0], end[1], 'o', ms=5 * font, color=colour, zorder=4)

            for itype in types:
                badge = INTERACTION_STYLE.get(itype, (itype, '#495057', '-', '?'))[3]
                bcolour = colour_of(itype)
                p = start + direction * (first + step * slot)
                slot += 1
                ax.text(p[0], p[1], badge, ha='center', va='center',
                        fontsize=fs(9.5), fontweight=weight, color=bcolour, zorder=8,
                        bbox=dict(boxstyle='circle,pad=0.28', fc='white', ec=bcolour,
                                  lw=1.5 * font))
            if show_distances and 'hydrogen_bonds' in by_type:
                best = min(c_['distance'] for c_ in by_type['hydrogen_bonds'])
                p = start + direction * max(first + step * slot + 0.35, length * 0.6)
                if np.linalg.norm(p - start) < length - 0.3:
                    ax.text(p[0], p[1], f"{best:.1f} Å", ha='center', va='center',
                            fontsize=fs(8.5), color=colour_of('hydrogen_bonds'),
                            zorder=8, bbox=dict(boxstyle='round,pad=0.15', fc='white',
                                                ec='none', alpha=0.9))

    if not residues:
        ax.text(0, -_R_LIGAND - 1.2, 'No interactions detected', ha='center',
                fontsize=fs(12), color='#868e96')

    if not own_figure:
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        return {'used_types': used_types, 'residue_classes': used_classes,
                'any_solvent': any_solvent, 'nearby': bool(nearby), 'limit': limit}

    # ---- legend: only what is on the plot ----------------------------------
    handles = [Line2D([0], [0], color=INTERACTION_STYLE[t][1], ls=INTERACTION_STYLE[t][2],
                      lw=2, label=INTERACTION_STYLE[t][0])
               for t in INTERACTION_STYLE if t in used_types]
    for label, fill in used_classes:
        handles.append(Line2D([0], [0], marker='s', ls='', ms=13, mfc=fill,
                              mec='#495057', label=label))
    if any_solvent:
        handles.append(Line2D([0], [0], marker='s', ls='', ms=13, mfc='white',
                              mec=_SOLVENT_EDGE, mew=3, label='Solvent accessible'))
    if nearby:
        handles.append(Line2D([0], [0], marker='s', ls='', ms=13, mfc='white',
                              mec='#495057', label='Nearby, no specific contact'))

    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    if not title:
        base = os.path.splitext(os.path.basename(str(mapper.structure_file)))[0]
        ligand = getattr(mapper.ligand_residue, 'resname', '').strip()
        title = f"Protein–ligand interactions: {base}" + (f" · {ligand}" if ligand else "")
    ax.set_title(title, fontsize=17, pad=6)
    rows = math.ceil(len(handles) / 4) if handles else 0
    if handles:
        fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False,
                   fontsize=11, bbox_to_anchor=(0.5, 0.005), handlelength=2.6,
                   columnspacing=1.6)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.03 + 0.024 * rows)
    fig.savefig(output_file, dpi=dpi, facecolor='white')
    plt.close(fig)

    logger.info(f"Interaction diagram saved to {output_file}")
    return output_file

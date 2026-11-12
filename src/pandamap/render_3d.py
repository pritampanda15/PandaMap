"""
3D contact view for PandaMap.

An "exploded" view of the binding site, laid out as cleanly as a 2D
interaction diagram:

* the ligand is drawn as solid 3D sticks in the centre, its best plane
  turned towards the viewer;
* every contacting residue is drawn as a transparent 3D stick model, moved
  straight outwards onto a ring around the ligand (and turned in the plane of
  the page so its contacting atoms face the ligand), so nothing is drawn on
  top of anything else;
* the residue atom that makes each contact is a solid sphere labelled with
  its atom name and joined to the ligand atom (or ring centre) it touches.

It answers "which atom of LEU 694 touches the ligand, and where?".  Because
the residues are moved apart, line lengths are not to scale; the distances
printed on the lines are the real ones.
"""
import logging
import math
import os
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.lines import Line2D

logger = logging.getLogger("pandamap")

from .diagram import (INTERACTION_STYLE, ATOM_COLOURS, _PRIORITY, _RING_TYPES,
                      _RESIDUE_CLASS, _OTHER_RESIDUE, _element, _ligand_bonds,
                      _smallest_ring, _spread_angles)

# stick colour of the transparent residues, by residue class label
_STICK_COLOURS = {
    'Hydrophobic residue': '#5c940d',
    'Polar residue': '#1971c2',
    'Positive residue': '#5f3dc4',
    'Negative residue': '#c92a2a',
    'Other residue': '#495057',
}
_CARBON = '#343a40'
_GAP = 3.2            # clear space between the ligand and the residue ring (A)
_SLOT = 6.6           # arc length reserved for one residue on the ring (A)


def _rotation(tilt_x, tilt_y):
    ax_, ay_ = math.radians(tilt_x), math.radians(tilt_y)
    rx = np.array([[1, 0, 0],
                   [0, math.cos(ax_), -math.sin(ax_)],
                   [0, math.sin(ax_), math.cos(ax_)]])
    ry = np.array([[math.cos(ay_), 0, math.sin(ay_)],
                   [0, 1, 0],
                   [-math.sin(ay_), 0, math.cos(ay_)]])
    return rx @ ry


def render_contacts_3d(mapper, output_file='contacts_3d.png', figsize=(11, 11.6),
                       dpi=300, title=None, tilt=(-10, 6), residue_alpha=0.38,
                       label_atoms=True, show_distances=True):
    """
    Draw the 3D contact view: solid ligand, transparent contacting residues
    arranged around it, and the contacting atoms highlighted and labelled.

    Parameters
    ----------
    mapper : HybridProtLigMapper
        A mapper on which ``detect_interactions()`` has been run.
    output_file : str
        Image path; the extension selects the format (.png, .svg, .pdf).
    tilt : (float, float)
        Rotation in degrees about the horizontal and vertical screen axes,
        applied after the ligand's best plane is turned to face the viewer.
    residue_alpha : float
        Opacity of the residue sticks (0 = invisible, 1 = solid).
    label_atoms : bool
        Label the contacting protein and ligand atoms with their atom names.
    show_distances : bool
        Print the real distance on hydrogen bonds.

    Returns the path written.
    """
    all_atoms = list(mapper.ligand_atoms)
    ligand_ids = {id(a) for a in all_atoms}
    heavy = [a for a in all_atoms if _element(a) != 'H'] or all_atoms
    xyz = np.array([a.coord for a in heavy], dtype=float)
    origin = xyz.mean(axis=0)

    # ---- view: ligand plane towards the viewer, then a small tilt ----------
    if len(heavy) >= 3:
        _, vectors = np.linalg.eigh(np.cov((xyz - origin).T))
        frame = vectors[:, ::-1].T                 # rows: long axis, short axis, normal
        if np.linalg.det(frame) < 0:
            frame[2] *= -1
    else:
        frame = np.eye(3)
    view = _rotation(*tilt) @ frame

    def project(coord):
        return view @ (np.asarray(coord, dtype=float) - origin)

    lig = np.array([project(a.coord) for a in heavy])
    index = {id(a): i for i, a in enumerate(heavy)}
    bonds = _ligand_bonds(heavy)
    neighbours = defaultdict(list)
    for i, j in bonds:
        neighbours[i].append(j)
        neighbours[j].append(i)

    # ---- one contact per residue and interaction type (the closest) --------
    best = {}
    for itype, items in mapper.interactions.items():
        for contact in items:
            res, patom, latom = (contact.get('protein_residue'),
                                 contact.get('protein_atom'), contact.get('ligand_atom'))
            if res is None or patom is None or latom is None:
                continue
            key = (res.get_parent().id, res.id, itype)
            if key not in best or contact['distance'] < best[key][1]['distance']:
                best[key] = (res, contact)

    residues = {}                                 # (chain, res.id) -> residue
    links = defaultdict(dict)                     # residue key -> {pair: link}
    # strongest interaction types first, so that an atom making several kinds
    # of contact gets one line, drawn for (and aimed like) the strongest
    for (chain, res_id, itype), (res, contact) in sorted(
            best.items(), key=lambda kv: _PRIORITY.get(kv[0][2], 99)):
        rkey = (chain, res_id)
        residues[rkey] = res
        patom, latom = contact['protein_atom'], contact['ligand_atom']
        end, lig_label, external = project(latom.coord), latom.get_name().strip(), None
        if id(latom) not in ligand_ids:
            external = latom                      # e.g. a bridging metal ion
        elif itype in _RING_TYPES and id(latom) in index:
            ring = _smallest_ring(index[id(latom)], neighbours)
            if ring:
                end, lig_label = lig[ring].mean(axis=0), None
        link = links[rkey].setdefault(patom.get_name(), {
            'types': [], 'atom_name': patom.get_name(), 'element': _element(patom),
            'end': end, 'lig_label': lig_label, 'external': external,
            'distance': contact['distance']})
        link['types'].append(itype)

    # ---- place the residues on a ring, next to the atoms they touch ---------
    keys = sorted(residues, key=lambda k: (k[0], k[1][1]))
    preferred = []
    for rkey in keys:
        v = np.mean([l['end'][:2] for l in links[rkey].values()], axis=0)
        preferred.append(math.atan2(v[1], v[0]))
    lig_radius = float(np.linalg.norm(lig[:, :2], axis=1).max())
    inner = max(lig_radius + _GAP, len(keys) * _SLOT / (2 * math.pi) - 2.5)
    order, angles = _spread_angles(preferred, _SLOT / (inner + 2.5))
    angle_of = {keys[i]: a for i, a in zip(order, angles)}

    placed = {}                                   # residue key -> drawing data
    outer = inner
    for rkey in keys:
        res = residues[rkey]
        atoms = [a for a in res if _element(a) != 'H']
        pos = np.array([project(a.coord) for a in atoms])
        pos[:, :2] -= pos[:, :2].mean(axis=0)
        names = [a.get_name() for a in atoms]
        theta = angle_of[rkey]
        out = np.array([math.cos(theta), math.sin(theta)])
        # turn the residue in the page so its contacting atoms face the ligand
        touching = [names.index(l['atom_name']) for l in links[rkey].values()
                    if l['atom_name'] in names]
        if touching:
            v = pos[touching, :2].mean(axis=0)
            if np.linalg.norm(v) > 0.3:
                turn = math.atan2(-out[1], -out[0]) - math.atan2(v[1], v[0])
                c, s = math.cos(turn), math.sin(turn)
                pos[:, :2] = pos[:, :2] @ np.array([[c, s], [-s, c]])
        along = pos[:, :2] @ out
        pos[:, :2] += out * (inner - along.min())
        far = inner + (along.max() - along.min())
        outer = max(outer, far)
        placed[rkey] = {'pos': pos, 'names': names, 'atoms': atoms, 'out': out, 'far': far}

    # ---- draw ---------------------------------------------------------------
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect('equal')
    ax.axis('off')
    limit = outer + 3.6
    scale = min(1.0, 21.0 / limit)                # keep text and sticks in proportion
    halo = [path_effects.withStroke(linewidth=2.6, foreground='white')]

    # residues: transparent sticks in the colour of their class, back to front
    used_classes = []
    multi_chain = len({chain for chain, _ in keys}) > 1
    for rkey in keys:
        res, data = residues[rkey], placed[rkey]
        pos = data['pos']
        label, fill = _RESIDUE_CLASS.get(res.resname.strip(), _OTHER_RESIDUE)
        if (label, fill) not in used_classes:
            used_classes.append((label, fill))
        colour = _STICK_COLOURS[label]
        sticks = sorted(_ligand_bonds(data['atoms']),
                        key=lambda b: pos[b[0]][2] + pos[b[1]][2])
        for i, j in sticks:
            ax.plot([pos[i][0], pos[j][0]], [pos[i][1], pos[j][1]], color=colour,
                    lw=5 * scale, alpha=residue_alpha, solid_capstyle='round', zorder=3)
        spot = data['out'] * (data['far'] + 2.0)
        text = f"{res.resname.strip()} {res.id[1]}" + (f" {rkey[0]}" if multi_chain else "")
        ax.text(spot[0], spot[1], text, ha='center', va='center', fontsize=11.5 * scale,
                fontweight='bold', color='#212529', zorder=9,
                bbox=dict(boxstyle='round,pad=0.35', fc=fill, ec=colour, lw=1.4))

    # ligand: solid sticks, each half coloured by its atom, back to front
    for i, j in sorted(bonds, key=lambda b: lig[b[0]][2] + lig[b[1]][2]):
        p, q = lig[i], lig[j]
        mid = (p + q) / 2
        ax.plot([p[0], q[0]], [p[1], q[1]], color='white', lw=7.5 * scale,
                solid_capstyle='round', zorder=5)
        ax.plot([p[0], mid[0]], [p[1], mid[1]], lw=5 * scale, solid_capstyle='round',
                color=ATOM_COLOURS.get(_element(heavy[i]), _CARBON), zorder=5.1)
        ax.plot([mid[0], q[0]], [mid[1], q[1]], lw=5 * scale, solid_capstyle='round',
                color=ATOM_COLOURS.get(_element(heavy[j]), _CARBON), zorder=5.1)
    for i, atom in enumerate(heavy):
        element = _element(atom)
        if element != 'C' or not neighbours[i]:
            ax.text(lig[i][0], lig[i][1], element.capitalize(), ha='center', va='center',
                    fontsize=11 * scale, fontweight='bold', zorder=7.5,
                    color=ATOM_COLOURS.get(element, _CARBON),
                    bbox=dict(boxstyle='circle,pad=0.1', fc='white', ec='none'))

    # contacts: solid contacting atom, line to the ligand atom / ring centre
    used_types = set()
    for rkey in keys:
        data = placed[rkey]
        for link in links[rkey].values():
            if link['atom_name'] not in data['names']:
                continue
            types = sorted(link['types'], key=lambda t: _PRIORITY.get(t, 99))
            used_types.update(types)
            _, colour, linestyle, _ = INTERACTION_STYLE.get(
                types[0], (types[0], '#495057', '-', '?'))
            p = data['pos'][data['names'].index(link['atom_name'])]
            q = link['end']
            ax.plot([p[0], q[0]], [p[1], q[1]], color='white', lw=4.2 * scale, zorder=6)
            ax.plot([p[0], q[0]], [p[1], q[1]], color=colour, ls=linestyle,
                    lw=2.2 * scale, zorder=6.1)
            ax.scatter(p[0], p[1], s=200 * scale ** 2, zorder=7, linewidths=2.2 * scale,
                       edgecolors=colour, color=ATOM_COLOURS.get(link['element'], '#868e96'))
            if label_atoms:
                side = np.array([-data['out'][1], data['out'][0]])
                ax.annotate(link['atom_name'].strip(), (p[0], p[1]),
                            xytext=tuple(side * 11 + data['out'] * 4),
                            textcoords='offset points', ha='center', va='center',
                            fontsize=9 * scale, fontweight='bold', color='#212529',
                            zorder=8, path_effects=halo)
            if link['external'] is not None:
                ax.scatter(q[0], q[1], s=330 * scale ** 2, color='#fff3bf',
                           edgecolors='#b8860b', linewidths=2, zorder=7)
                ax.text(q[0], q[1], _element(link['external']).capitalize(), ha='center',
                        va='center', fontsize=9 * scale, fontweight='bold',
                        color='#7a5c00', zorder=8)
            else:
                ax.scatter(q[0], q[1], s=55 * scale ** 2, color='white', edgecolors=colour,
                           linewidths=1.8 * scale, zorder=7)
                if label_atoms and link['lig_label'] and not link['lig_label'][:1] in 'NOSPF':
                    ax.annotate(link['lig_label'], (q[0], q[1]), xytext=(6, -10),
                                textcoords='offset points', fontsize=7.5 * scale,
                                color='#495057', zorder=8, path_effects=halo)
            if show_distances and 'hydrogen_bonds' in types:
                mid = p + (q - p) * 0.45
                ax.text(mid[0], mid[1], f"{link['distance']:.1f} Å", ha='center',
                        va='center', fontsize=8.5 * scale, color=colour, zorder=8,
                        bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none',
                                  alpha=0.92))

    if not keys:
        ax.text(0, lig[:, 1].min() - 2, 'No interactions detected', ha='center',
                fontsize=12, color='#868e96')
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)

    # ---- title, note and legend ---------------------------------------------
    if not title:
        base = os.path.splitext(os.path.basename(str(mapper.structure_file)))[0]
        ligand = getattr(mapper.ligand_residue, 'resname', '').strip()
        title = f"3D contacts: {base}" + (f" · {ligand}" if ligand else "")
    ax.set_title(title, fontsize=17, pad=6)
    handles = [Line2D([0], [0], color=INTERACTION_STYLE[t][1], ls=INTERACTION_STYLE[t][2],
                      lw=2.2, label=INTERACTION_STYLE[t][0])
               for t in INTERACTION_STYLE if t in used_types]
    for label, _ in used_classes:
        handles.append(Line2D([0], [0], color=_STICK_COLOURS[label], lw=5,
                              alpha=max(residue_alpha, 0.45), label=label))
    if used_types:
        handles.append(Line2D([0], [0], marker='o', ls='', ms=11, mfc='#868e96',
                              mec='#212529', mew=1.8, label='Contacting atom'))
    rows = math.ceil(len(handles) / 4) if handles else 0
    bottom = 0.03 + 0.024 * rows
    if handles:
        fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False,
                   fontsize=11, bbox_to_anchor=(0.5, 0.005), handlelength=2.6,
                   columnspacing=1.6)
    if keys:
        fig.text(0.5, bottom + 0.012,
                 'Residues are moved outwards for clarity; printed distances are the real ones.',
                 ha='center', fontsize=9.5, color='#868e96')
    fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=bottom + 0.03)
    fig.savefig(output_file, dpi=dpi, facecolor='white')
    plt.close(fig)

    logger.info(f"3D contact view saved to {output_file}")
    return output_file

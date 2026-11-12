"""
Graphical interaction report (``--plots``): one publication-style figure that
explains a binding site at a glance.

Panels
------
a  Interaction map      the 2D diagram: which residue touches which part of
                        the ligand, and how
b  Contacts by type     how many contacts of each kind
c  Contact distances    every contact against the typical distance for its kind
d  Residue profile      residue x interaction-type matrix, with solvent exposure
e  Estimated dG         contribution of each interaction type to the estimate

Styling follows common journal figure guidance: 183 mm wide, one sans-serif
family, 5-7 pt text, 8 pt bold lowercase panel letters, thin axes, no
decoration, and one colour per interaction family used identically in every
panel.
"""
import logging
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from .diagram import INTERACTION_STYLE, draw_interaction_diagram

logger = logging.getLogger("pandamap")

INK = '#424242'
GRID = '#E6E6E6'
EXPOSED = '#4DBBD5'

# interaction family -> (label, colour, member types); one colour per family
FAMILIES = (
    ('Hydrogen / halogen bond', '#00A087', ('hydrogen_bonds', 'halogen_bonds', 'water_bridges')),
    ('Electrostatic', '#E64B35', ('salt_bridge', 'ionic', 'attractive_charge', 'repulsion')),
    ('π interaction', '#3C5488', ('pi_pi_stacking', 'pi_cation', 'cation_pi', 'donor_pi',
                                  'amide_pi', 'carbon_pi', 'alkyl_pi')),
    ('Hydrophobic', '#F39B7F', ('hydrophobic',)),
    ('Metal / covalent', '#3D1B24', ('metal_coordination', 'covalent')),
)
FAMILY_COLOUR = {t: colour for _, colour, members in FAMILIES for t in members}
_TYPE_ORDER = [t for _, _, members in FAMILIES for t in members]

# typical contact distance per interaction type (A), marked in panel c
TYPICAL_DISTANCE = {
    'hydrogen_bonds': 2.9, 'covalent': 1.8, 'metal_coordination': 2.2, 'ionic': 3.5,
    'salt_bridge': 3.5, 'halogen_bonds': 3.0, 'pi_pi_stacking': 3.8, 'cation_pi': 4.0,
    'pi_cation': 4.0, 'hydrophobic': 3.5, 'carbon_pi': 4.0, 'donor_pi': 3.5,
    'amide_pi': 3.5, 'alkyl_pi': 4.0, 'attractive_charge': 3.5, 'repulsion': 3.0,
}

_RC = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'Nimbus Sans', 'DejaVu Sans'],
    'font.size': 6, 'axes.labelsize': 7, 'axes.titlesize': 7,
    'xtick.labelsize': 5.5, 'ytick.labelsize': 6, 'legend.fontsize': 5.5,
    'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK,
    'xtick.color': INK, 'ytick.color': INK,
    'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'xtick.major.size': 2, 'ytick.major.size': 2, 'xtick.major.pad': 2, 'ytick.major.pad': 2,
    'axes.spines.top': False, 'axes.spines.right': False,
    'svg.fonttype': 'none', 'pdf.fonttype': 42,
}


def _label(itype):
    return INTERACTION_STYLE[itype][0] if itype in INTERACTION_STYLE \
        else itype.replace('_', ' ').capitalize()


def _colour(itype):
    return FAMILY_COLOUR.get(itype, INK)


def _panel(fig, ax, letter, heading, dx=0.0):
    """Panel letter (8 pt bold) and a compact heading, at a fixed offset."""
    box = ax.get_position()
    x, y = box.x0 - dx, box.y1 + 0.012
    fig.text(x, y, letter, fontsize=8, fontweight='bold', ha='left', va='bottom')
    fig.text(x + 0.022, y + 0.0015, heading, fontsize=7, ha='left', va='bottom')


def draw_report_figure(mapper, output_file='interaction_report.png', dpi=300):
    """
    Draw the five-panel graphical report for an analysed mapper.

    ``output_file`` may end in .png, .pdf or .svg (text stays editable in the
    vector formats).  Returns the path written, or None when no interaction
    was detected.
    """
    present = {t: v for t, v in mapper.interactions.items() if v}
    if not present:
        logger.warning("No interactions detected; skipping plots.")
        return None
    order = [t for t in _TYPE_ORDER if t in present] + \
            [t for t in present if t not in _TYPE_ORDER]
    n_types = len(order)

    # ---- tidy records ------------------------------------------------------
    solvent = set(getattr(mapper, 'solvent_accessible', ()) or ())
    residues = {}                         # (chain, number, name) -> {type: closest distance}
    for itype in order:
        for c in present[itype]:
            res = c.get('protein_residue')
            if res is None:
                continue
            key = (res.get_parent().id, res.id[1], res.resname.strip())
            best = residues.setdefault(key, {})
            best[itype] = min(best.get(itype, 99.0), float(c['distance']))
    res_keys = sorted(residues)
    truncated = len(res_keys) > 32
    if truncated:                         # keep the residues with most contacts
        keep = sorted(res_keys, key=lambda k: (-len(residues[k]), min(residues[k].values())))[:32]
        res_keys = sorted(keep)
    multi_chain = len({k[0] for k in residues}) > 1
    res_label = lambda k: f"{k[2]} {k[1]}" + (f" {k[0]}" if multi_chain else "")
    n_contacts = sum(len(v) for v in present.values())

    try:
        affinity = mapper.estimate_binding_affinity()
    except Exception as exc:              # scoring is optional
        affinity, affinity_error = None, str(exc)

    with plt.rc_context(_RC):
        fig = plt.figure(figsize=(7.2, 7.7))          # 183 mm wide

        # ---- header: what this is, in one line -------------------------------
        lig = mapper.ligand_residue
        base = os.path.splitext(os.path.basename(str(mapper.structure_file)))[0]
        fig.text(0.02, 0.987, f"{base} · ligand {lig.resname.strip()} "
                 f"(chain {lig.get_parent().id}, residue {lig.id[1]})",
                 fontsize=7, fontweight='bold', va='top')
        facts = [f"{n_contacts} contact{'s' if n_contacts != 1 else ''} with "
                 f"{len(residues)} residue{'s' if len(residues) != 1 else ''}"]
        if 'hydrogen_bonds' in present:
            hb = min(present['hydrogen_bonds'], key=lambda c: c['distance'])
            r = hb['protein_residue']
            facts.append(f"{len(present['hydrogen_bonds'])} hydrogen bond"
                         f"{'s' if len(present['hydrogen_bonds']) != 1 else ''}, closest "
                         f"{r.resname.strip()} {r.id[1]} {hb['protein_atom'].get_name().strip()}"
                         f"–{hb['ligand_atom'].get_name().strip()} at {hb['distance']:.2f} Å")
        exposed = sum(1 for k in residues if (k[2], k[1]) in solvent)
        if solvent:
            facts.append(f"{exposed} of {len(residues)} contact residues solvent exposed")
        if affinity:
            facts.append(f"interaction score {affinity['score']:.1f}".replace('-', '−'))
        fig.text(0.02, 0.968, "  ·  ".join(facts), fontsize=6, va='top')

        # ---- a: the interaction map -------------------------------------------
        ax_a = fig.add_axes([0.015, 0.385, 0.585, 0.535])
        width_in = min(ax_a.get_position().width * fig.get_figwidth(),
                       ax_a.get_position().height * fig.get_figheight())
        info = draw_interaction_diagram(mapper, ax=ax_a, theme={
            'scale': width_in / 10.6, 'min_font': 5, 'badge_step': 0.9, 'badge_start': 0.75, 'weight': 'normal', 'ink': INK,
            'colours': {t: _colour(t) for t in INTERACTION_STYLE},
            'box_fill': '#F6F6F6', 'box_edge': '#9A9A9A', 'solvent_edge': EXPOSED,
            'atom_colours': {'N': '#3C5488', 'O': '#E64B35', 'S': '#F39B7F', 'P': '#F39B7F',
                             'F': '#00A087', 'CL': '#00A087', 'BR': '#00A087', 'I': '#00A087'},
        })
        _panel(fig, ax_a, 'a', 'Interaction map', dx=-0.005)

        # shared key under the map: families, badges, solvent exposure
        fams = [(lab, col) for lab, col, members in FAMILIES if any(t in present for t in members)]
        handles = [Line2D([0], [0], color=col, lw=1.6, label=lab) for lab, col in fams]
        if info.get('any_solvent'):
            handles.append(Line2D([0], [0], marker='s', ls='', ms=5, mfc='#F6F6F6', mec=EXPOSED,
                                  mew=1.2, label='Solvent-exposed residue'))
        fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.02, 0.378), ncol=3,
                   frameon=False, handlelength=1.6, columnspacing=1.4, handletextpad=0.5,
                   borderaxespad=0, labelspacing=0.35)
        badges = "   ".join(f"{INTERACTION_STYLE[t][3]} {_label(t).lower()}"
                            for t in order if t in INTERACTION_STYLE)
        fig.text(0.02, 0.338, "Badges:  " + badges, fontsize=5, va='top', color='#6b6b6b',
                 wrap=True)

        rows = np.arange(n_types)
        labels = [_label(t) for t in order]

        # ---- b: contacts by type -------------------------------------------------
        ax_b = fig.add_axes([0.735, 0.70, 0.245, 0.22])
        counts = [len(present[t]) for t in order]
        ax_b.barh(rows, counts, height=0.62, color=[_colour(t) for t in order], linewidth=0)
        for y, n in zip(rows, counts):
            ax_b.text(n + max(counts) * 0.03, y, str(n), va='center', fontsize=5.5)
        ax_b.set_yticks(rows)
        ax_b.set_yticklabels(labels)
        ax_b.invert_yaxis()
        ax_b.set_xlim(0, max(counts) * 1.15)
        ax_b.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
        ax_b.set_xlabel('Contacts (n)')
        ax_b.tick_params(axis='y', length=0)
        _panel(fig, ax_b, 'b', 'Contacts by type', dx=0.115)

        # ---- c: contact distances --------------------------------------------------
        ax_c = fig.add_axes([0.735, 0.39, 0.245, 0.22])
        for y, t in zip(rows, order):
            d = sorted(float(c['distance']) for c in present[t])
            offsets = (np.arange(len(d)) - (len(d) - 1) / 2) * min(0.14, 0.6 / max(len(d), 1))
            if t in TYPICAL_DISTANCE:
                ax_c.plot([TYPICAL_DISTANCE[t]] * 2, [y - 0.36, y + 0.36], color=INK, lw=0.6,
                          zorder=2)
            ax_c.scatter(d, y + offsets, s=9, color=_colour(t), edgecolors='white',
                         linewidths=0.3, zorder=3)
        ax_c.set_yticks(rows)
        ax_c.set_yticklabels(labels)
        ax_c.set_ylim(n_types - 0.4, -0.6)
        ax_c.set_xlabel('Distance (Å)')
        ax_c.tick_params(axis='y', length=0)
        ax_c.grid(axis='x', color=GRID, lw=0.4, zorder=0)
        ax_c.set_axisbelow(True)
        ax_c.text(1.0, 1.02, '│ typical distance', transform=ax_c.transAxes, ha='right',
                  va='bottom', fontsize=5, color='#6b6b6b')
        _panel(fig, ax_c, 'c', 'Contact distances', dx=0.115)

        # ---- d: residue profile ------------------------------------------------------
        ax_d = fig.add_axes([0.135, 0.08, 0.50, 0.185])
        cols = np.arange(len(res_keys))
        for y in range(n_types + 1):
            ax_d.axhline(y, color=GRID, lw=0.4, zorder=0)
        for x, key in zip(cols, res_keys):
            for y, t in enumerate(order):
                if t in residues[key]:
                    ax_d.scatter(x, y, s=20, color=_colour(t), linewidths=0, zorder=3)
            if (key[2], key[1]) in solvent:
                ax_d.scatter(x, n_types, s=16, facecolors='white', edgecolors=EXPOSED,
                             linewidths=0.9, zorder=3)
        ax_d.set_yticks(np.arange(n_types + 1))
        ax_d.set_yticklabels(labels + ['Solvent exposed'])
        ax_d.set_ylim(n_types + 0.6, -0.6)
        ax_d.set_xticks(cols)
        ax_d.set_xticklabels([res_label(k) for k in res_keys], rotation=90)
        ax_d.set_xlim(-0.7, len(res_keys) - 0.3)
        ax_d.tick_params(length=0)
        for side in ('left', 'bottom'):
            ax_d.spines[side].set_visible(False)
        heading = 'Residue profile' + (' (32 residues with most contacts)' if truncated else '')
        _panel(fig, ax_d, 'd', heading, dx=0.115)

        # ---- e: estimated free energy -------------------------------------------------
        ax_e = fig.add_axes([0.795, 0.08, 0.185, 0.185])
        if affinity:
            items = [(k, v['contribution_kcal_mol']) for k, v in affinity['breakdown'].items()
                     if isinstance(v, dict) and 'contribution_kcal_mol' in v]
            items.sort(key=lambda kv: (_TYPE_ORDER.index(kv[0]) if kv[0] in _TYPE_ORDER else 99))
            ys = np.arange(len(items))
            for y, (k, v) in zip(ys, items):
                penalty = k not in present
                ax_e.barh(y, v, height=0.62, linewidth=0.4 if penalty else 0,
                          color='white' if penalty else _colour(k),
                          edgecolor=INK, hatch='////' if penalty else None)
            ax_e.axvline(0, color=INK, lw=0.5)
            ax_e.set_yticks(ys)
            ax_e.set_yticklabels([_label(k) if k in present else
                                  k.replace('_', ' ').capitalize() for k, _ in items])
            ax_e.invert_yaxis()
            ax_e.tick_params(axis='y', length=0)
            ax_e.set_xlabel('Score contribution (kcal/mol scale)')
            total = f"{affinity['score']:.2f}".replace('-', '−')
            ax_e.text(1.0, 1.02, f"total {total}", transform=ax_e.transAxes, ha='right',
                      va='bottom', fontsize=5.5)
            fig.text(0.98, 0.012, 'Empirical, uncalibrated; not a binding free energy. Hatched = penalty',
                     fontsize=5, ha='right', va='bottom', color='#6b6b6b')
        else:
            ax_e.axis('off')
            ax_e.text(0.5, 0.5, 'Interaction score unavailable', ha='center', va='center',
                      transform=ax_e.transAxes, fontsize=6)
        _panel(fig, ax_e, 'e', 'Interaction score', dx=0.135)

        fig.savefig(output_file, dpi=dpi, facecolor='white')
        plt.close(fig)

    logger.info(f"Graphical interaction report saved to {output_file}")
    return output_file

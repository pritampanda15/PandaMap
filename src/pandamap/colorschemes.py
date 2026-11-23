"""
Colour schemes for the modern 2D interaction diagram (from PandaMap-Color).

A scheme is a ``theme`` dict for ``draw_interaction_diagram``:

    colours        interaction type -> line and badge colour
    atom_colours   element (upper case) -> ligand atom label colour
    residue_fills  residue class label -> box fill
    ink            bonds and text
    paper          background, and the fill behind labels and badges
    box_edge, solvent_edge   residue box outlines

Any key left out falls back to the default design.  Use a built-in name
(``default``, ``colorblind``, ``monochrome``, ``dark``, ``publication``) or a
JSON file with the keys above.  JSON files written for PandaMap-Color
(``interaction_styles``, ``element_colors``, ``residue_colors``,
``background_color``) are also accepted.
"""
import json
import os

_OKABE_ITO = {'orange': '#E69F00', 'sky': '#56B4E9', 'green': '#009E73',
              'yellow': '#F0E442', 'blue': '#0072B2', 'vermilion': '#D55E00',
              'purple': '#CC79A7', 'grey': '#666666'}

COLOR_SCHEMES = {
    'default': {},
    'colorblind': {
        # Okabe & Ito (2008); line style and badge letter carry the type too
        'colours': {
            'covalent': '#000000', 'metal_coordination': _OKABE_ITO['orange'],
            'hydrogen_bonds': _OKABE_ITO['green'], 'water_bridges': _OKABE_ITO['sky'],
            'salt_bridge': _OKABE_ITO['vermilion'], 'ionic': _OKABE_ITO['vermilion'],
            'attractive_charge': _OKABE_ITO['blue'], 'halogen_bonds': _OKABE_ITO['sky'],
            'pi_pi_stacking': _OKABE_ITO['blue'], 'pi_cation': _OKABE_ITO['purple'],
            'cation_pi': _OKABE_ITO['purple'], 'donor_pi': _OKABE_ITO['purple'],
            'amide_pi': _OKABE_ITO['orange'], 'carbon_pi': _OKABE_ITO['grey'],
            'alkyl_pi': _OKABE_ITO['blue'], 'hydrophobic': _OKABE_ITO['grey'],
            'repulsion': _OKABE_ITO['vermilion'],
        },
        'atom_colours': {'N': '#0072B2', 'O': '#D55E00', 'S': '#B8930A', 'P': '#CC79A7',
                         'F': '#009E73', 'CL': '#009E73', 'BR': '#882255', 'I': '#882255'},
        'residue_fills': {'Hydrophobic residue': '#FBF3C6', 'Polar residue': '#DDEFFA',
                          'Positive residue': '#D6E6F2', 'Negative residue': '#FBE3D2',
                          'Other residue': '#F2F2F2'},
        'solvent_edge': '#0072B2',
    },
    'monochrome': {
        'colours': {
            'covalent': '#000000', 'metal_coordination': '#333333',
            'hydrogen_bonds': '#111111', 'water_bridges': '#444444',
            'salt_bridge': '#222222', 'ionic': '#222222', 'attractive_charge': '#333333',
            'halogen_bonds': '#444444', 'pi_pi_stacking': '#333333', 'pi_cation': '#333333',
            'cation_pi': '#333333', 'donor_pi': '#555555', 'amide_pi': '#555555',
            'carbon_pi': '#777777', 'alkyl_pi': '#555555', 'hydrophobic': '#888888',
            'repulsion': '#000000'},
        'atom_colours': {e: '#000000' for e in ('N', 'O', 'S', 'P', 'F', 'CL', 'BR', 'I')},
        'residue_fills': {'Hydrophobic residue': '#EEEEEE', 'Polar residue': '#FFFFFF',
                          'Positive residue': '#DDDDDD', 'Negative residue': '#CCCCCC',
                          'Other residue': '#F7F7F7'},
        'ink': '#000000', 'box_edge': '#333333', 'solvent_edge': '#000000',
    },
    'dark': {
        'colours': {
            'covalent': '#F1F3F5', 'metal_coordination': '#FFD43B',
            'hydrogen_bonds': '#51CF66', 'water_bridges': '#3BC9DB',
            'salt_bridge': '#FF922B', 'ionic': '#FF922B', 'attractive_charge': '#4DABF7',
            'halogen_bonds': '#66D9E8', 'pi_pi_stacking': '#B197FC', 'pi_cation': '#F783AC',
            'cation_pi': '#F783AC', 'donor_pi': '#FAA2C1', 'amide_pi': '#E8A87C',
            'carbon_pi': '#ADB5BD', 'alkyl_pi': '#74C0FC', 'hydrophobic': '#ADB5BD',
            'repulsion': '#FF6B6B',
        },
        'atom_colours': {'N': '#74C0FC', 'O': '#FF8787', 'S': '#FFE066', 'P': '#FFA94D',
                         'F': '#8CE99A', 'CL': '#8CE99A', 'BR': '#E599A0', 'I': '#D0BFFF'},
        'residue_fills': {'Hydrophobic residue': '#2F3B24', 'Polar residue': '#22324A',
                          'Positive residue': '#2E2C52', 'Negative residue': '#4A2626',
                          'Other residue': '#2B2F33'},
        'ink': '#E9ECEF', 'paper': '#16181B', 'box_edge': '#868E96',
        'solvent_edge': '#4DABF7',
    },
    'publication': {
        # Nature Publishing Group palette, as in the summary figure
        'colours': {
            'covalent': '#3C3C3C', 'metal_coordination': '#B09C85',
            'hydrogen_bonds': '#00A087', 'water_bridges': '#4DBBD5',
            'salt_bridge': '#E64B35', 'ionic': '#E64B35', 'attractive_charge': '#3C5488',
            'halogen_bonds': '#4DBBD5', 'pi_pi_stacking': '#8491B4', 'pi_cation': '#DC0000',
            'cation_pi': '#DC0000', 'donor_pi': '#F39B7F', 'amide_pi': '#7E6148',
            'carbon_pi': '#999999', 'alkyl_pi': '#3C5488', 'hydrophobic': '#999999',
            'repulsion': '#DC0000',
        },
        'atom_colours': {'N': '#3C5488', 'O': '#E64B35', 'S': '#B09C85', 'P': '#F39B7F',
                         'F': '#00A087', 'CL': '#00A087', 'BR': '#7E6148', 'I': '#8491B4'},
        'residue_fills': {'Hydrophobic residue': '#F2F2F2', 'Polar residue': '#E6F0F5',
                          'Positive residue': '#E8EAF2', 'Negative residue': '#FBEAE7',
                          'Other residue': '#FFFFFF'},
        'ink': '#222222', 'box_edge': '#555555', 'solvent_edge': '#4DBBD5',
        'weight': 'normal',
    },
}

# PandaMap-Color residue_colors keys -> residue class labels used here
_LEGACY_RESIDUES = {'hydrophobic': 'Hydrophobic residue', 'polar': 'Polar residue',
                    'charged_pos': 'Positive residue', 'charged_neg': 'Negative residue',
                    'default': 'Other residue'}


def _from_legacy(data):
    """Translate a PandaMap-Color JSON scheme to theme keys."""
    theme = {k: v for k, v in data.items() if k in (
        'colours', 'atom_colours', 'residue_fills', 'ink', 'paper',
        'box_edge', 'solvent_edge', 'weight')}
    if 'interaction_styles' in data:
        theme.setdefault('colours', {}).update(
            {t: s['color'] for t, s in data['interaction_styles'].items() if 'color' in s})
    if 'element_colors' in data:
        theme.setdefault('atom_colours', {}).update(
            {e.upper(): c for e, c in data['element_colors'].items()})
    if 'residue_colors' in data:
        theme.setdefault('residue_fills', {}).update(
            {_LEGACY_RESIDUES[k]: c for k, c in data['residue_colors'].items()
             if k in _LEGACY_RESIDUES})
    if 'background_color' in data:
        theme.setdefault('paper', data['background_color'])
    return theme


def load_color_scheme(scheme=None):
    """
    Theme dict for a built-in scheme name or a JSON file path.
    ``None`` gives the default design.
    """
    if scheme is None or scheme == 'default':
        return {}
    if isinstance(scheme, dict):
        return _from_legacy(scheme)
    if scheme in COLOR_SCHEMES:
        return dict(COLOR_SCHEMES[scheme])
    if os.path.isfile(scheme):
        with open(scheme, encoding='utf-8') as fh:
            return _from_legacy(json.load(fh))
    raise ValueError(f"Unknown colour scheme '{scheme}': use one of "
                     f"{', '.join(COLOR_SCHEMES)} or a JSON file")

import os
import base64
import tempfile
import numpy as np
from Bio.PDB import PDBIO

def create_pandamap_3d_viz(mapper, output_file, width=800, height=600, show_surface=True):
    interaction_colors = {
        'hydrogen_bonds': 'green', 'carbon_pi': 'gray', 'pi_pi_stacking': 'purple',
        'donor_pi': 'pink', 'amide_pi': 'brown', 'hydrophobic': 'lightgray',
        'ionic': 'orange', 'halogen_bonds': 'cyan', 'cation_pi': 'magenta',
        'metal_coordination': 'gold', 'salt_bridge': 'tomato', 'covalent': 'black',
        'alkyl_pi': 'steelblue', 'attractive_charge': 'blue', 'pi_cation': 'magenta', 'repulsion': 'red'
    }

    temp_pdb = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
    io = PDBIO()
    io.set_structure(mapper.structure)
    io.save(temp_pdb.name)
    with open(temp_pdb.name, 'rb') as f:
        pdb_encoded = base64.b64encode(f.read()).decode()
    os.unlink(temp_pdb.name)

    ligand_resname = getattr(mapper.ligand_residue, 'resname', 'UNK')
    js_residues, js_solvent, js_interactions, interaction_data = [], [], [], []

    for res_id in mapper.interacting_residues:
        residue = mapper.protein_residues.get(res_id)
        if residue is None:
            continue
        res_num = str(residue.id[1])
        chain_id = getattr(residue.get_parent(), 'id', 'A')
        resname = residue.get_resname()

        js_residues.append(f"""
    viewer.setStyle({{chain: '{chain_id}', resi: {res_num}}},
                    {{stick: {{colorscheme: 'amino', radius: 0.2}}, cartoon: {{color: 'lightblue'}}}});""")

        js_residues.append(f"""
    viewer.addSphere({{
        center: viewer.getModel().selectedAtoms({{chain: '{chain_id}', resi: {res_num}}})[0],
        radius: 0.3,
        color: 'white',
        alpha: 0,
        clickable: true,
        callback: function() {{
            showTooltip("Residue: {resname} {res_num} (Chain {chain_id})",
                        viewer.getModel().selectedAtoms({{chain: '{chain_id}', resi: {res_num}}})[0]);
        }}
    }});""")

        if show_surface and res_id in mapper.solvent_accessible:
            js_solvent.append(f"""
    viewer.addSurface($3Dmol.VDW, {{opacity: 0.4, color: 'white'}}, {{chain: '{chain_id}', resi: {res_num}}});""")

    interaction_id = 0
    for itype, interactions in mapper.interactions.items():
        if itype not in interaction_colors or not interactions:
            continue
        color = interaction_colors[itype]
        is_dashed = itype in ['hydrophobic', 'carbon_pi', 'pi_pi_stacking']
        for interaction in interactions:
            pa = interaction.get('protein_atom')
            la = interaction.get('ligand_atom')
            if not pa or not la:
                continue

            prot_coords = pa.get_coord() if hasattr(pa, 'get_coord') else getattr(pa, 'coord', None)
            lig_coords = la.get_coord() if hasattr(la, 'get_coord') else getattr(la, 'coord', None)
            if prot_coords is None or lig_coords is None:
                continue

            prot_coords = list(prot_coords)
            lig_coords = list(lig_coords)
            length = round(np.linalg.norm(np.array(prot_coords) - np.array(lig_coords)), 2)
            if length > 5.0:
                continue
            label = itype.replace('_', ' ').title()
            mid_x, mid_y, mid_z = [(p + l) / 2 for p, l in zip(prot_coords, lig_coords)]

            js_interactions.append(f"""
    viewer.addCylinder({{
        start: {{x: {prot_coords[0]}, y: {prot_coords[1]}, z: {prot_coords[2]}}},
        end: {{x: {lig_coords[0]}, y: {lig_coords[1]}, z: {lig_coords[2]}}},
        radius: 0.1, color: '{color}', dashed: {str(is_dashed).lower()}
    }});
    viewer.addSphere({{
        center: {{x: {mid_x}, y: {mid_y}, z: {mid_z}}},
        radius: 0.15,
        color: '{color}',
        opacity: 1.0
    }});
    viewer.addSphere({{
        center: {{x: {lig_coords[0]}, y: {lig_coords[1]}, z: {lig_coords[2]}}},
        radius: 0.3,
        color: 'white',
        alpha: 0,
        clickable: true,
        callback: function() {{
            showTooltip("Ligand: {ligand_resname} {interaction['ligand_atom']}", {{x: {lig_coords[0]}, y: {lig_coords[1]}, z: {lig_coords[2]}}});
        }}
    }});
    viewer.addSphere({{
        center: {{x: {mid_x}, y: {mid_y}, z: {mid_z}}},
        radius: 0.3,
        color: 'white',
        alpha: 0,
        clickable: true,
        callback: function() {{ selectInteraction({interaction_id}); }}
    }});""")
            interaction_data.append(f"{{id: {interaction_id}, label: '{label} {length:.2f} Å', center: {{x: {mid_x}, y: {mid_y}, z: {mid_z}}}}}")
            interaction_id += 1

    js_legend = """
    var legend = document.createElement('div');
    legend.style.position = 'absolute';
    legend.style.top = '10px';
    legend.style.right = '10px';
    legend.style.background = 'rgba(255,255,255,0.8)';
    legend.style.padding = '10px';
    legend.style.borderRadius = '5px';
    legend.style.fontFamily = 'Arial, sans-serif';
    legend.style.zIndex = '10';
    legend.style.boxShadow = '0 0 5px rgba(0,0,0,0.2)';
    var title = document.createElement('h3');
    title.textContent = 'Interaction Types';
    legend.appendChild(title);
    """
    for itype in sorted(interaction_colors):
        color = interaction_colors[itype]
        name = itype.replace('_', ' ').title()
        js_legend += f"""
    var row = document.createElement('div');
    row.style.display = 'flex'; row.style.alignItems = 'center'; row.style.marginBottom = '5px';
    var colorBox = document.createElement('div');
    colorBox.style.width = '15px'; colorBox.style.height = '15px'; colorBox.style.backgroundColor = '{color}'; colorBox.style.marginRight = '5px';
    var label = document.createElement('span');
    label.textContent = '{name}';
    row.appendChild(colorBox); row.appendChild(label);
    legend.appendChild(row);
    """
    js_legend += "document.getElementById('container').appendChild(legend);"

    html = f"""<!DOCTYPE html>
<html><head>
  <title>PandaMap 3D Visualization</title>
  <script src='https://ajax.googleapis.com/ajax/libs/jquery/3.6.0/jquery.min.js'></script>
  <script src='https://3Dmol.org/build/3Dmol-min.js'></script>
  <style>
    html, body {{ margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden; }}
    #container {{ width: 100%; height: 100%; position: relative; }}
    #controls {{ position: absolute; bottom: 10px; left: 10px; background: white; padding: 10px; border-radius: 6px; box-shadow: 0 0 5px rgba(0,0,0,0.3); font-family: sans-serif; z-index: 20; }}
  </style>
</head>
<body>
  <div id='container'></div>
  <div id='controls'>
    <button onclick='viewer.zoomTo(); viewer.render();'>Reset View</button>
    <button onclick="exportImage();">Export PNG</button>
    <button onclick='clearSelection();'>Clear Selection</button>
    <button onclick='viewer.setStyle({{hetflag: true}}, {{stick: {{colorscheme: "elementColors", radius: 0.3}}}}); viewer.zoomTo({{hetflag: true}}); viewer.render();'>Show Ligand</button>


  </div>
  <script>
    var viewer;
    var pdbData = atob("{pdb_encoded}");
    var interactionList = [{','.join(interaction_data)}];
    var labelHandle = null;

    function showTooltip(text, position) {{
        if (labelHandle) viewer.removeLabel(labelHandle);
        labelHandle = viewer.addLabel(text, {{
            position: position,
            backgroundColor: 'white',
            fontColor: 'black',
            showBackground: true,
            inFront: true
        }});
        viewer.render();
    }}

    function selectInteraction(id) {{
      const data = interactionList.find(i => i.id === id);
      if (!data) return;
      showTooltip(data.label, data.center);
    }}

    function clearSelection() {{
        if (labelHandle) viewer.removeLabel(labelHandle);
        labelHandle = null;
        viewer.render();
    }}

    function selectLigand(resname) {{
        viewer.setStyle({{hetflag: true, resn: resname}}, {{stick: {{colorscheme: 'elementColors', radius: 0.3}}}});
        viewer.zoomTo({{hetflag: true, resn: resname}});
        viewer.render();
    }}

    function exportImage() {{
        viewer.getCanvas().toBlob(function(blob) {{
            var link = document.createElement('a');
            link.href = URL.createObjectURL(blob);
            link.download = 'pandamap_visualization.png';
            link.click();
        }});
        }}

    $(document).ready(function() {{
      viewer = $3Dmol.createViewer($('#container'), {{backgroundColor: 'white'}});
      viewer.addModel(pdbData, 'pdb');
      viewer.setStyle({{hetflag: false}}, {{cartoon: {{color: 'lightgray', opacity: 0.8}}}});
      viewer.setStyle({{hetflag: true, resn: '{ligand_resname}'}}, {{stick: {{colorscheme: 'elementColors', radius: 0.3}}}});
      {'viewer.addSurface($3Dmol.VDW, {opacity: 0.4, color: "white"}, {hetflag: false});' if show_surface else ''}
      {''.join(js_residues)}
      {''.join(js_solvent)}
      {''.join(js_interactions)}
      {js_legend}
      viewer.zoomTo();
      viewer.render();
    }});
  </script>
</body></html>"""

    with open(output_file, 'w') as f:
        f.write(html)
    print(f"3D visualization saved to {output_file}")
    return output_file
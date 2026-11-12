"""
Reading structure files and choosing the ligand.

Parsers for PDB, mmCIF and PDBQT, the list of HETATM groups that are never
the ligand of interest, and the rules for selecting one ligand residue.
"""
import logging
import os
import tempfile

from Bio.PDB import PDBParser, MMCIFParser

logger = logging.getLogger("pandamap")


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
    with open(pdbqt_file, 'r') as f_pdbqt, open(temp_pdb_path, 'w', encoding='utf-8') as f_pdb:
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
            # Unknown extension (e.g. .ent, .pdb1, no extension): decide from
            # the content rather than refusing or prompting the user.
            if _looks_like_mmcif(file_path):
                return self.mmcif_parser.get_structure('complex', file_path)
            return self.pdb_parser.get_structure('complex', file_path)


def _looks_like_mmcif(file_path):
    """True if the first meaningful line of the file is an mmCIF data block."""
    try:
        with open(file_path, 'r', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith('#'):
                    return line.startswith('data_')
    except OSError:
        pass
    return False


# HETATM residue names that are never the ligand of interest.
WATER_RESNAMES = frozenset({'HOH', 'WAT', 'H2O', 'DOD', 'TIP', 'TIP3', 'SOL'})
# Ions, buffer/cryo additives and modified amino acids: valid when requested
# explicitly, but skipped when the ligand is auto-detected.
NON_LIGAND_RESNAMES = frozenset({
    'NA', 'K', 'LI', 'CL', 'BR', 'IOD', 'MG', 'CA', 'ZN', 'FE', 'FE2', 'MN',
    'CU', 'CU1', 'CO', 'NI', 'CD', 'HG', 'SR', 'BA', 'CS', 'RB',
    'SO4', 'PO4', 'NO3', 'ACT', 'ACE', 'NH2', 'NH4', 'FMT', 'CO3', 'SCN', 'AZI',
    'GOL', 'EDO', 'PEG', 'PGE', 'PG4', '1PE', 'P6G', 'MPD', 'DMS', 'IPA',
    'EOH', 'MOH', 'BME', 'DTT', 'TRS', 'MES', 'EPE', 'CIT', 'IMD',
    'MSE', 'SEP', 'TPO', 'PTR', 'CSO', 'CME', 'KCX', 'MLY', 'HYP', 'PCA',
})


def _describe_ligand(residue):
    return (f"{residue.resname.strip()} (chain {residue.get_parent().id}, "
            f"residue {residue.id[1]}, {len(residue)} atoms)")


def _ligand_candidates(model):
    """All non-water HETATM residues of a model, largest first."""
    candidates = [r for r in model.get_residues()
                  if r.id[0] != ' ' and r.resname.strip().upper() not in WATER_RESNAMES]
    # sorted() is stable, so equal-sized copies keep file order
    return sorted(candidates, key=lambda r: -len(r))


def list_ligands(structure_file):
    """
    List the candidate ligands (non-water HETATM residues) in a structure file.

    Returns a list of dicts with keys ``resname``, ``chain``, ``resnum`` and
    ``n_atoms``, largest first.  Use it to find the name to pass as
    ``ligand_resname`` / ``--ligand``.
    """
    structure = MultiFormatParser().parse_structure(structure_file)
    models = list(structure.get_models())
    if not models:
        raise ValueError(f"No atoms could be read from {structure_file}. "
                         "Is it a valid PDB, mmCIF or PDBQT file?")
    return [{'resname': r.resname.strip(), 'chain': r.get_parent().id,
             'resnum': r.id[1], 'n_atoms': len(r)}
            for r in _ligand_candidates(models[0])]


def _select_ligand(model, structure_file, resname=None, chain=None, resnum=None):
    """Pick the single HETATM residue to analyse as the ligand."""
    candidates = _ligand_candidates(model)
    if not candidates:
        raise ValueError(f"No ligand (HETATM) found in the file: {structure_file}")
    available = ", ".join(_describe_ligand(r) for r in candidates[:15])

    matches = candidates
    if resname is not None:
        wanted = str(resname).strip().upper()
        matches = [r for r in matches if r.resname.strip().upper() == wanted]
    if chain is not None:
        matches = [r for r in matches if r.get_parent().id == str(chain)]
    if resnum is not None:
        matches = [r for r in matches if r.id[1] == int(resnum)]

    if not matches:
        wanted = " ".join(
            part for part in (
                f"'{resname}'" if resname is not None else "",
                f"chain {chain}" if chain is not None else "",
                f"residue {resnum}" if resnum is not None else "",
            ) if part)
        raise ValueError(f"Ligand {wanted} not found in {structure_file}. "
                         f"Available ligands: {available}")

    if resname is None:
        # Auto-detect: prefer real small molecules over ions and additives.
        preferred = [r for r in matches
                     if r.resname.strip().upper() not in NON_LIGAND_RESNAMES]
        matches = preferred or matches
        ligand = matches[0]
        others = [r for r in candidates if r is not ligand]
        if others and resnum is None:
            logger.info(f"Auto-selected ligand {_describe_ligand(ligand)}. "
                  f"Other HETATM groups: "
                  f"{', '.join(_describe_ligand(r) for r in others[:8])}. "
                  "Use --ligand / ligand_resname to choose another.")
        return ligand

    # Several copies of the named ligand: analyse one, in file order.
    match_ids = {id(r) for r in matches}
    matches = [r for r in model.get_residues() if id(r) in match_ids]
    ligand = matches[0]
    if len(matches) > 1:
        logger.info(f"Found {len(matches)} copies of {ligand.resname.strip()}; using "
              f"chain {ligand.get_parent().id} residue {ligand.id[1]}. "
              "Use --chain / --resnum (ligand_chain / ligand_resnum) to choose another.")
    return ligand

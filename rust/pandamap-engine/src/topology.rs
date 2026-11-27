//! Fast topology and multi-model reader (PDB, PDBQT, GRO) plus PandaMap
//! ligand selection, ring/aromaticity perception, formal charges, and snapshot
//! PDB writer.

use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::geometry::{is_planar, Vec3};

#[derive(Debug, Clone)]
pub struct AtomRecord {
    pub index: usize,
    pub serial: i32,
    pub name: String,
    pub resname: String,
    pub chain: String,
    pub resnum: i32,
    /// Monotonic residue block ID (disambiguates >99,999 wrapped GRO waters)
    pub res_uid: usize,
    pub element: String, // Uppercase normalized: "C", "N", "O", "CL", "BR", etc.
    pub is_hetatm: bool,
    pub coord: Vec3,
}

#[derive(Debug, Clone)]
pub struct Topology {
    pub atoms: Vec<AtomRecord>,
    /// Multi-model frames when the input is a multi-MODEL PDB or multi-pose PDBQT.
    pub models: Vec<Vec<[f32; 3]>>,
}

#[derive(Debug, Clone)]
pub struct PreparedSystem {
    pub atoms: Vec<AtomRecord>,
    pub protein_indices: Vec<usize>,
    pub ligand_indices: Vec<usize>,
    pub water_o_indices: Vec<usize>,
    pub metal_indices: Vec<usize>,
    pub halogen_indices: Vec<usize>,
    pub ligand_name: String,
    pub ligand_chain: String,
    pub ligand_resnum: i32,
    /// Per-atom flags indexed by global atom index (for O(1) hot-loop lookup)
    pub is_lig_aromatic: Vec<bool>,
    pub is_lig_cation: Vec<bool>,
    pub is_lig_anion: Vec<bool>,
    pub is_prot_ring_atom: Vec<bool>,
    /// Explicit hydrogen indices attached to each heavy atom (same residue, < 1.3 A)
    pub attached_h: Vec<Vec<usize>>,
    /// Heavy covalent neighbors in same residue / ligand (0.9..1.8 A) for strict H-bond proxy
    pub heavy_cov_neighbors: Vec<Vec<usize>>,
    /// Number of freely rotatable non-ring ligand bonds (for empirical score)
    pub n_rotatable_bonds: usize,
}

pub fn is_water_resname(resname: &str) -> bool {
    matches!(
        resname,
        "HOH" | "WAT" | "H2O" | "DOD" | "TIP" | "TIP3" | "SOL" | "SPC" | "TIP4"
    )
}

pub fn is_non_ligand_resname(resname: &str) -> bool {
    matches!(
        resname,
        "NA" | "K" | "LI" | "CL" | "BR" | "IOD" | "MG" | "CA" | "ZN" | "FE" | "FE2" | "MN"
            | "CU" | "CU1" | "CO" | "NI" | "CD" | "HG" | "SR" | "BA" | "CS" | "RB"
            | "SO4" | "PO4" | "NO3" | "ACT" | "ACE" | "NH2" | "NH4" | "FMT" | "CO3"
            | "SCN" | "AZI" | "GOL" | "EDO" | "PEG" | "PGE" | "PG4" | "1PE" | "P6G"
            | "MPD" | "DMS" | "IPA" | "EOH" | "MOH" | "BME" | "DTT" | "TRS" | "MES"
            | "EPE" | "CIT" | "IMD" | "MSE" | "SEP" | "TPO" | "PTR" | "CSO" | "CME"
            | "KCX" | "MLY" | "HYP" | "PCA" | "SOD" | "CLA" | "POT"
            // Common MD membrane lipids (CHARMM-GUI / Martini / Lipid21)
            | "POPC" | "POPE" | "POPS" | "POPG" | "DPPC" | "DMPC" | "DOPC" | "DOPE"
            | "DLPC" | "DSPC" | "SOPC" | "SDPC" | "CHL1" | "CHOL" | "SSM" | "PSM"
            | "POPI" | "SAPI" | "CARD"
    )
}

pub fn is_lipid_resname(resname: &str) -> bool {
    matches!(
        resname,
        "POPC" | "POPE" | "POPS" | "POPG" | "DPPC" | "DMPC" | "DOPC" | "DOPE"
            | "DLPC" | "DSPC" | "SOPC" | "SDPC" | "CHL1" | "CHOL" | "SSM" | "PSM"
            | "POPI" | "SAPI" | "CARD"
    )
}

pub fn is_standard_aa(resname: &str) -> bool {
    matches!(
        resname,
        "ALA" | "ARG" | "ASN" | "ASP" | "CYS" | "GLN" | "GLU" | "GLY" | "HIS"
            | "ILE" | "LEU" | "LYS" | "MET" | "PHE" | "PRO" | "SER" | "THR" | "TRP"
            | "TYR" | "VAL" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP"
            | "CYX" | "ASH" | "GLH" | "LYN"
    )
}

pub fn is_virtual_site_atom(atom_name: &str) -> bool {
    let n = atom_name.trim().to_uppercase();
    // CHARMM/CGenFF/AMBER/TIP4P lone-pair & dummy virtual sites (e.g. LP1, LP2, EP1, MW, DUM)
    if matches!(n.as_str(), "LP" | "EP" | "MW" | "DUM" | "DU" | "VS") {
        return true;
    }
    (n.starts_with("LP") || n.starts_with("EP") || n.starts_with("VS"))
        && n.len() >= 3
        && n.as_bytes()[2].is_ascii_digit()
}

pub fn is_metal_element(el: &str) -> bool {
    matches!(
        el,
        "MG" | "ZN" | "CA" | "FE" | "MN" | "CU" | "NA" | "K" | "LI" | "CO" | "NI"
    )
}

pub fn is_halogen_element(el: &str, strict: bool) -> bool {
    if strict {
        matches!(el, "CL" | "BR" | "I")
    } else {
        matches!(el, "F" | "CL" | "BR" | "I")
    }
}

pub fn is_covalent_ligand_element(el: &str) -> bool {
    matches!(el, "C" | "N" | "O" | "S" | "P" | "SE" | "B")
}

pub fn is_prot_aromatic_ring_atom(resname: &str, atom_name: &str) -> bool {
    match resname {
        "PHE" | "TYR" => matches!(atom_name, "CG" | "CD1" | "CD2" | "CE1" | "CE2" | "CZ"),
        "TRP" => matches!(
            atom_name,
            "CG" | "CD1" | "CD2" | "CE2" | "CE3" | "CZ2" | "CZ3" | "CH2" | "NE1"
        ),
        "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP" => {
            matches!(atom_name, "CG" | "ND1" | "CD2" | "CE1" | "NE2")
        }
        _ => false,
    }
}

fn map_pdbqt_element(ad_type: &str) -> Option<&'static str> {
    match ad_type {
        "A" | "C" => Some("C"),
        "N" | "NA" | "NS" => Some("N"),
        "O" | "OA" | "OS" => Some("O"),
        "S" | "SA" => Some("S"),
        "H" | "HD" | "HS" => Some("H"),
        "P" => Some("P"),
        "F" => Some("F"),
        "CL" => Some("CL"),
        "BR" => Some("BR"),
        "I" => Some("I"),
        "SI" => Some("SI"),
        "B" => Some("B"),
        "SE" => Some("SE"),
        "FE" => Some("FE"),
        "ZN" => Some("ZN"),
        "MG" => Some("MG"),
        "MN" => Some("MN"),
        "CA" => Some("CA"),
        "CU" => Some("CU"),
        "CO" => Some("CO"),
        "NI" => Some("NI"),
        "K" => Some("K"),
        "LI" => Some("LI"),
        _ => None,
    }
}

fn guess_element(atom_name: &str, raw_elem: &str, is_pdbqt: bool, is_hetatm: bool) -> String {
    let re = raw_elem.trim().to_uppercase();
    if is_pdbqt {
        if let Some(m) = map_pdbqt_element(&re) {
            return m.to_string();
        }
    }
    if !re.is_empty() {
        // Strip charge suffixes like "N1+" or "O1-"
        let alpha: String = re.chars().filter(|c| c.is_ascii_alphabetic()).collect();
        if !alpha.is_empty() {
            return alpha;
        }
    }
    // Guess from atom name
    let trimmed = atom_name.trim().to_uppercase();
    let letters: String = trimmed
        .chars()
        .skip_while(|c| c.is_ascii_digit())
        .take_while(|c| c.is_ascii_alphabetic())
        .collect();
    if is_hetatm && letters.len() >= 2 {
        let two = &letters[..2];
        if matches!(
            two,
            "CL" | "BR" | "ZN" | "MG" | "FE" | "CA" | "MN" | "CU" | "CO" | "NI" | "NA" | "SE" | "SI"
        ) {
            return two.to_string();
        }
    }
    letters.chars().next().map(|c| c.to_string()).unwrap_or_else(|| "C".to_string())
}

fn covalent_radius(el: &str) -> f64 {
    match el {
        "H" => 0.31,
        "B" => 0.84,
        "C" => 0.76,
        "N" => 0.71,
        "O" => 0.66,
        "F" => 0.57,
        "SI" => 1.11,
        "P" => 1.07,
        "S" => 1.05,
        "CL" => 1.02,
        "SE" => 1.20,
        "BR" => 1.20,
        "I" => 1.39,
        _ => 1.30,
    }
}

/// Load a PDB, multi-MODEL PDB, PDBQT (single or multi-pose), or GRO file.
pub fn load_topology<P: AsRef<Path>>(path: P) -> io::Result<Topology> {
    let p = path.as_ref();
    let ext = p
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    if ext == "gro" {
        return load_gro(p);
    }
    let is_pdbqt = ext == "pdbqt";
    load_pdb_or_pdbqt(p, is_pdbqt)
}

fn load_pdb_or_pdbqt(path: &Path, is_pdbqt: bool) -> io::Result<Topology> {
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(256 * 1024, file);

    let mut first_model_atoms: Vec<AtomRecord> = Vec::new();
    let mut models: Vec<Vec<[f32; 3]>> = Vec::new();
    let mut current_coords: Vec<[f32; 3]> = Vec::new();
    let mut in_first_model = true;
    let mut saw_model_record = false;

    let mut current_res_uid: usize = 0;
    let mut prev_res_sig: Option<(String, String, i32, bool)> = None;

    for line_res in reader.lines() {
        let line = line_res?;
        if line.starts_with("MODEL") {
            saw_model_record = true;
            current_coords.clear();
            continue;
        }
        if line.starts_with("ENDMDL") {
            if !current_coords.is_empty() {
                models.push(std::mem::take(&mut current_coords));
            }
            in_first_model = false;
            continue;
        }
        if !(line.starts_with("ATOM  ") || line.starts_with("HETATM")) {
            continue;
        }
        if line.len() < 54 {
            continue;
        }
        let x: f64 = line[30..38].trim().parse().unwrap_or(0.0);
        let y: f64 = line[38..46].trim().parse().unwrap_or(0.0);
        let z: f64 = line[46..54].trim().parse().unwrap_or(0.0);
        current_coords.push([x as f32, y as f32, z as f32]);

        if in_first_model {
            let is_hetatm = line.starts_with("HETATM");
            let serial: i32 = line[6..11].trim().parse().unwrap_or(first_model_atoms.len() as i32 + 1);
            let name = line[12..16].trim().to_string();
            let resname = line[17..20].trim().to_uppercase();
            let chain = line
                .get(21..22)
                .unwrap_or(" ")
                .trim()
                .to_string();
            let chain = if chain.is_empty() { " ".to_string() } else { chain };
            let resnum: i32 = line[22..26].trim().parse().unwrap_or(0);
            let sig = (chain.clone(), resname.clone(), resnum, is_hetatm);
            if prev_res_sig.as_ref() != Some(&sig) {
                current_res_uid += 1;
                prev_res_sig = Some(sig);
            }
            let raw_elem = if is_pdbqt {
                line.get(76..).unwrap_or("")
            } else {
                line.get(76..78).unwrap_or("")
            };
            let element = guess_element(&name, raw_elem, is_pdbqt, is_hetatm);
            let idx = first_model_atoms.len();
            first_model_atoms.push(AtomRecord {
                index: idx,
                serial,
                name,
                resname,
                chain,
                resnum,
                res_uid: current_res_uid,
                element,
                is_hetatm,
                coord: Vec3::new(x, y, z),
            });
        }
    }

    if !current_coords.is_empty() {
        models.push(current_coords);
    }
    if first_model_atoms.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("No atoms could be read from {}", path.display()),
        ));
    }
    if !saw_model_record && models.is_empty() {
        let coords = first_model_atoms
            .iter()
            .map(|a| [a.coord.x as f32, a.coord.y as f32, a.coord.z as f32])
            .collect();
        models.push(coords);
    }

    Ok(Topology {
        atoms: first_model_atoms,
        models,
    })
}

fn load_gro(path: &Path) -> io::Result<Topology> {
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(256 * 1024, file);
    let mut lines = reader.lines();
    let _title = lines.next().ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidData, "Empty GRO file")
    })??;
    let natoms_line = lines.next().ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidData, "Missing atom count in GRO file")
    })??;
    let natoms: usize = natoms_line.trim().parse().map_err(|e| {
        io::Error::new(io::ErrorKind::InvalidData, format!("Invalid GRO natoms: {}", e))
    })?;

    let mut atoms = Vec::with_capacity(natoms);
    let mut coords = Vec::with_capacity(natoms);
    let mut current_res_uid: usize = 0;
    let mut prev_res_sig: Option<(i32, String)> = None;

    for i in 0..natoms {
        let line = lines.next().ok_or_else(|| {
            io::Error::new(io::ErrorKind::UnexpectedEof, "Truncated GRO atom section")
        })??;
        if line.len() < 44 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Short GRO atom line {}", i + 3),
            ));
        }
        let resnum: i32 = line[0..5].trim().parse().unwrap_or(0);
        let resname = line[5..10].trim().to_uppercase();
        let sig = (resnum, resname.clone());
        if prev_res_sig.as_ref() != Some(&sig) {
            current_res_uid += 1;
            prev_res_sig = Some(sig);
        }
        let name = line[10..15].trim().to_string();
        let serial: i32 = line[15..20].trim().parse().unwrap_or(i as i32 + 1);
        // GRO coordinates are in nm -> convert to Angstrom
        let x = line[20..28].trim().parse::<f64>().unwrap_or(0.0) * 10.0;
        let y = line[28..36].trim().parse::<f64>().unwrap_or(0.0) * 10.0;
        let z = line[36..44].trim().parse::<f64>().unwrap_or(0.0) * 10.0;
        let is_hetatm = !is_standard_aa(&resname);
        let element = guess_element(&name, "", false, is_hetatm);
        coords.push([x as f32, y as f32, z as f32]);
        atoms.push(AtomRecord {
            index: i,
            serial,
            name,
            resname,
            chain: "A".to_string(),
            resnum,
            res_uid: current_res_uid,
            element,
            is_hetatm,
            coord: Vec3::new(x, y, z),
        });
    }

    Ok(Topology {
        atoms,
        models: vec![coords],
    })
}

/// Prepare a topology for fast per-frame interaction detection.
pub fn prepare_system(
    atoms: Vec<AtomRecord>,
    ligand_resname: Option<&str>,
    ligand_chain: Option<&str>,
    ligand_resnum: Option<i32>,
) -> io::Result<PreparedSystem> {
    // Group residues by unique contiguous block ID (handles >99,999 wrapped GRO waters)
    #[derive(Clone, Hash, PartialEq, Eq)]
    struct ResKey {
        res_uid: usize,
        chain: String,
        resname: String,
        resnum: i32,
        is_hetatm: bool,
    }

    let mut order: Vec<ResKey> = Vec::new();
    let mut groups: HashMap<ResKey, Vec<usize>> = HashMap::new();
    for (idx, a) in atoms.iter().enumerate() {
        let key = ResKey {
            res_uid: a.res_uid,
            chain: a.chain.clone(),
            resname: a.resname.clone(),
            resnum: a.resnum,
            is_hetatm: a.is_hetatm,
        };
        if !groups.contains_key(&key) {
            order.push(key.clone());
        }
        groups.entry(key).or_default().push(idx);
    }

    // Candidate ligands: non-water HETATM residues (or any matching explicitly named residue)
    let mut candidates: Vec<&ResKey> = order
        .iter()
        .filter(|k| k.is_hetatm && !is_water_resname(&k.resname))
        .collect();
    // If explicit ligand name was given and wasn't flagged HETATM (e.g. GRO), allow matching it
    if let Some(wanted) = ligand_resname {
        let w = wanted.trim().to_uppercase();
        if !candidates.iter().any(|k| k.resname == w) {
            candidates = order.iter().filter(|k| k.resname == w).collect();
        }
    }
    if candidates.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "No ligand (HETATM) found in topology",
        ));
    }

    // Sort candidates largest first (stable, matching PandaMap's _ligand_candidates)
    let mut by_size = candidates.clone();
    by_size.sort_by(|a, b| groups[*b].len().cmp(&groups[*a].len()));

    let mut matches: Vec<&ResKey> = by_size.clone();
    if let Some(wanted) = ligand_resname {
        let w = wanted.trim().to_uppercase();
        matches.retain(|k| k.resname == w);
    }
    if let Some(ch) = ligand_chain {
        let c = ch.trim();
        matches.retain(|k| k.chain.trim() == c);
    }
    if let Some(rn) = ligand_resnum {
        matches.retain(|k| k.resnum == rn);
    }

    if matches.is_empty() {
        let avail: Vec<String> = by_size
            .iter()
            .take(10)
            .map(|k| format!("{} (chain {}, res {}, {} atoms)", k.resname, k.chain.trim(), k.resnum, groups[*k].len()))
            .collect();
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Requested ligand not found. Available: {}", avail.join(", ")),
        ));
    }

    let selected: ResKey = if ligand_resname.is_none() {
        let preferred: Vec<&ResKey> = matches
            .iter()
            .copied()
            .filter(|k| !is_non_ligand_resname(&k.resname))
            .collect();
        if !preferred.is_empty() {
            preferred[0].clone()
        } else {
            matches[0].clone()
        }
    } else {
        // File order among matching copies (matching PandaMap's _select_ligand)
        let match_set: HashSet<&ResKey> = matches.into_iter().collect();
        order
            .iter()
            .find(|k| match_set.contains(k))
            .cloned()
            .unwrap()
    };

    let n = atoms.len();
    let mut protein_indices = Vec::new();
    let mut ligand_indices = Vec::new();
    let mut water_o_indices = Vec::new();
    let mut metal_indices = Vec::new();
    let mut halogen_indices = Vec::new();
    let mut is_prot_ring_atom = vec![false; n];

    for key in &order {
        let idxs = &groups[key];
        if *key == selected {
            for &idx in idxs {
                if is_virtual_site_atom(&atoms[idx].name) {
                    continue;
                }
                ligand_indices.push(idx);
                if is_halogen_element(&atoms[idx].element, false) {
                    halogen_indices.push(idx);
                }
            }
        } else if key.is_hetatm {
            if is_water_resname(&key.resname) {
                for &idx in idxs {
                    if atoms[idx].element == "O" {
                        water_o_indices.push(idx);
                    }
                }
            } else if idxs.len() == 1 && is_metal_element(&atoms[idxs[0]].element) {
                metal_indices.push(idxs[0]);
            }
        } else {
            for &idx in idxs {
                protein_indices.push(idx);
                if is_prot_aromatic_ring_atom(&atoms[idx].resname, &atoms[idx].name) {
                    is_prot_ring_atom[idx] = true;
                }
            }
        }
    }

    // Pre-compute ligand bonds (covalent radii + 0.45 for aromaticity; < 1.75 for charge; < 2.0 for rings/rotors)
    let mut is_lig_aromatic = vec![false; n];
    let mut is_lig_cation = vec![false; n];
    let mut is_lig_anion = vec![false; n];

    // Heavy ligand atoms
    let heavy_lig: Vec<usize> = ligand_indices
        .iter()
        .copied()
        .filter(|&i| atoms[i].element != "H")
        .collect();

    // Build covalent neighbor graph (ligand2d rules) for aromatic ring perception
    let m = heavy_lig.len();
    let mut lig_adj_cov: Vec<Vec<usize>> = vec![Vec::new(); m];
    for a in 0..m {
        let ra = covalent_radius(&atoms[heavy_lig[a]].element);
        for b in (a + 1)..m {
            let rb = covalent_radius(&atoms[heavy_lig[b]].element);
            if atoms[heavy_lig[a]].coord.dist(atoms[heavy_lig[b]].coord) < ra + rb + 0.45 {
                lig_adj_cov[a].push(b);
                lig_adj_cov[b].push(a);
            }
        }
    }

    // Find smallest rings (up to size 8) and mark planar aromatic 5/6 rings
    let mut seen_rings: HashSet<Vec<usize>> = HashSet::new();
    for start in 0..m {
        let ring = smallest_ring(start, &lig_adj_cov, 8);
        if ring.len() == 5 || ring.len() == 6 {
            let mut key = ring.clone();
            key.sort_unstable();
            if seen_rings.insert(key) {
                let all_arom_el = ring.iter().all(|&i| {
                    matches!(atoms[heavy_lig[i]].element.as_str(), "C" | "N" | "O" | "S")
                });
                let all_sp2_deg = ring.iter().all(|&i| lig_adj_cov[i].len() <= 3);
                if all_arom_el && all_sp2_deg {
                    let pts: Vec<Vec3> = ring.iter().map(|&i| atoms[heavy_lig[i]].coord).collect();
                    if is_planar(&pts, 0.20) {
                        let mean_len: f64 = (0..ring.len())
                            .map(|k| pts[k].dist(pts[(k + 1) % ring.len()]))
                            .sum::<f64>()
                            / (ring.len() as f64);
                        if mean_len < 1.45 {
                            for &i in &ring {
                                is_lig_aromatic[heavy_lig[i]] = true;
                            }
                        }
                    }
                }
            }
        }
    }

    // Build 1.75 A ligand neighbor lists for formal charge perception (matching core.py)
    let mut lig_nb_175: HashMap<usize, Vec<usize>> = HashMap::new();
    for &i in &heavy_lig {
        let mut nbs = Vec::new();
        for &j in &heavy_lig {
            if i != j && atoms[i].coord.dist(atoms[j].coord) < 1.75 {
                nbs.push(j);
            }
        }
        lig_nb_175.insert(i, nbs);
    }
    for &i in &heavy_lig {
        let el = atoms[i].element.as_str();
        let nbs = &lig_nb_175[&i];
        if el == "O" && nbs.len() == 1 {
            let centre = nbs[0];
            if matches!(atoms[centre].element.as_str(), "C" | "P" | "S" | "N") {
                let term_o = lig_nb_175[&centre]
                    .iter()
                    .filter(|&&nb| atoms[nb].element == "O" && lig_nb_175[&nb].len() == 1)
                    .count();
                if term_o >= 2 {
                    is_lig_anion[i] = true;
                }
            }
        }
        if el == "N" {
            if nbs.len() >= 3 {
                is_lig_cation[i] = true;
            } else {
                for &nb in nbs {
                    if atoms[nb].element == "C" {
                        let n_cnt = lig_nb_175[&nb]
                            .iter()
                            .filter(|&&x| atoms[x].element == "N")
                            .count();
                        if n_cnt >= 2 {
                            is_lig_cation[i] = true;
                            break;
                        }
                    }
                }
            }
        }
    }

    // Rotatable bonds (2.0 A bond graph + iterative leaf pruning for ring atoms, matching core.py)
    let mut adj20: Vec<Vec<usize>> = vec![Vec::new(); m];
    let mut bonds20: Vec<(usize, usize)> = Vec::new();
    for a in 0..m {
        for b in (a + 1)..m {
            if atoms[heavy_lig[a]].coord.dist(atoms[heavy_lig[b]].coord) < 2.0 {
                adj20[a].push(b);
                adj20[b].push(a);
                bonds20.push((a, b));
            }
        }
    }
    let mut deg: Vec<usize> = adj20.iter().map(|v| v.len()).collect();
    let mut changed = true;
    while changed {
        changed = false;
        for a in 0..m {
            if deg[a] == 1 {
                deg[a] = 0;
                for &nb in &adj20[a] {
                    if deg[nb] > 0 {
                        deg[nb] -= 1;
                    }
                }
                changed = true;
            }
        }
    }
    let n_rotatable_bonds = bonds20
        .iter()
        .filter(|&&(a, b)| adj20[a].len() >= 2 && adj20[b].len() >= 2 && deg[a] < 2 && deg[b] < 2)
        .count();

    // Pre-compute explicit attached H (< 1.3 A) and heavy covalent neighbors (0.9..1.8 A) per residue
    let mut attached_h: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut heavy_cov_neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
    for idxs in groups.values() {
        for &i in idxs {
            if atoms[i].element == "H" {
                continue;
            }
            for &j in idxs {
                if i == j {
                    continue;
                }
                let d = atoms[i].coord.dist(atoms[j].coord);
                if atoms[j].element == "H" {
                    if d < 1.3 {
                        attached_h[i].push(j);
                    }
                } else if d > 0.9 && d < 1.8 {
                    heavy_cov_neighbors[i].push(j);
                }
            }
        }
    }

    Ok(PreparedSystem {
        atoms,
        protein_indices,
        ligand_indices,
        water_o_indices,
        metal_indices,
        halogen_indices,
        ligand_name: selected.resname,
        ligand_chain: selected.chain,
        ligand_resnum: selected.resnum,
        is_lig_aromatic,
        is_lig_cation,
        is_lig_anion,
        is_prot_ring_atom,
        attached_h,
        heavy_cov_neighbors,
        n_rotatable_bonds,
    })
}

fn smallest_ring(start: usize, adj: &[Vec<usize>], max_size: usize) -> Vec<usize> {
    let mut best: Vec<usize> = Vec::new();
    for &first in &adj[start] {
        let mut parent: HashMap<usize, Option<usize>> = HashMap::new();
        parent.insert(first, None);
        let mut queue = VecDeque::new();
        queue.push_back(first);
        while let Some(node) = queue.pop_front() {
            if node == start {
                break;
            }
            for &nxt in &adj[node] {
                if parent.contains_key(&nxt) || (node == first && nxt == start) {
                    continue;
                }
                parent.insert(nxt, Some(node));
                queue.push_back(nxt);
            }
        }
        if !parent.contains_key(&start) {
            continue;
        }
        let mut ring = Vec::new();
        let mut cur = Some(start);
        while let Some(node) = cur {
            ring.push(node);
            cur = parent.get(&node).copied().flatten();
        }
        if ring.len() >= 3 && ring.len() <= max_size && (best.is_empty() || ring.len() < best.len()) {
            best = ring;
        }
    }
    best
}

/// Write a single snapshot PDB containing protein + ligand + metals + waters within `water_shell_a`
/// of the ligand, using the frame's coordinates (in Angstrom).
pub fn write_snapshot_pdb<P: AsRef<Path>>(
    sys: &PreparedSystem,
    coords: &[[f32; 3]],
    out_path: P,
    water_shell_a: f64,
    time_ps: f32,
) -> io::Result<()> {
    // Compute ligand AABB for fast solvent shell filtering
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut min_z = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    let mut max_z = f64::NEG_INFINITY;
    for &idx in &sys.ligand_indices {
        let c = coords[idx];
        min_x = min_x.min(c[0] as f64);
        min_y = min_y.min(c[1] as f64);
        min_z = min_z.min(c[2] as f64);
        max_x = max_x.max(c[0] as f64);
        max_y = max_y.max(c[1] as f64);
        max_z = max_z.max(c[2] as f64);
    }
    let pad = water_shell_a;

    // Collect water res_uids whose oxygen is within ligand AABB + pad
    let mut keep_waters: HashSet<usize> = HashSet::new();
    for &idx in &sys.water_o_indices {
        let c = coords[idx];
        let (x, y, z) = (c[0] as f64, c[1] as f64, c[2] as f64);
        if x >= min_x - pad
            && x <= max_x + pad
            && y >= min_y - pad
            && y <= max_y + pad
            && z >= min_z - pad
            && z <= max_z + pad
        {
            keep_waters.insert(sys.atoms[idx].res_uid);
        }
    }

    let lig_res_uids: HashSet<usize> = sys
        .ligand_indices
        .iter()
        .map(|&i| sys.atoms[i].res_uid)
        .collect();

    let file = File::create(out_path)?;
    let mut w = BufWriter::with_capacity(128 * 1024, file);
    writeln!(
        w,
        "REMARK   1 PandaMap MD snapshot  t = {:.3} ps ({:.4} ns)",
        time_ps,
        time_ps / 1000.0
    )?;

    let mut serial: i32 = 1;
    let mut water_resnum_remap: HashMap<usize, i32> = HashMap::new();
    let mut next_water_num: i32 = 1;
    for (idx, atom) in sys.atoms.iter().enumerate() {
        if is_virtual_site_atom(&atom.name) {
            continue;
        }
        let is_lig = lig_res_uids.contains(&atom.res_uid);
        if atom.is_hetatm && !is_lig {
            if is_water_resname(&atom.resname) {
                if !keep_waters.contains(&atom.res_uid) {
                    continue;
                }
            } else if is_lipid_resname(&atom.resname) {
                // Skip bulk membrane lipids in snapshot PDBs so BioPython parses quickly
                continue;
            } else {
                // Keep metals/ions only if within water_shell_a of ligand
                let c = coords[idx];
                let (x, y, z) = (c[0] as f64, c[1] as f64, c[2] as f64);
                if x < min_x - pad || x > max_x + pad || y < min_y - pad || y > max_y + pad || z < min_z - pad || z > max_z + pad {
                    continue;
                }
            }
        }
        let rec = if atom.is_hetatm { "HETATM" } else { "ATOM  " };
        let aname = if atom.name.len() >= 4 {
            format!("{:<4}", &atom.name[..4])
        } else {
            format!(" {:<3}", atom.name)
        };
        let (chain_ch, resnum_out, resname_out) = if atom.is_hetatm && !is_lig && is_water_resname(&atom.resname) {
            let wn = *water_resnum_remap.entry(atom.res_uid).or_insert_with(|| {
                let n = next_water_num;
                next_water_num += 1;
                n
            });
            ('W', wn % 10000, "HOH")
        } else {
            let ch = if is_lig {
                let c = atom.chain.chars().next().unwrap_or('L');
                if c == ' ' { 'L' } else { c }
            } else {
                let c = atom.chain.chars().next().unwrap_or('A');
                if c == ' ' { 'A' } else { c }
            };
            let rn = if atom.resname.len() > 3 { &atom.resname[..3] } else { &atom.resname };
            (ch, atom.resnum % 10000, rn)
        };
        let c = coords[idx];
        let serial_clamped = serial % 100000;
        let elem = if atom.element.len() > 2 {
            &atom.element[..2]
        } else {
            &atom.element
        };
        writeln!(
            w,
            "{}{:5} {} {:>3} {}{:4}    {:8.3}{:8.3}{:8.3}  1.00 20.00          {:>2}",
            rec,
            serial_clamped,
            aname,
            resname_out,
            chain_ch,
            resnum_out,
            c[0],
            c[1],
            c[2],
            elem
        )?;
        serial += 1;
    }
    writeln!(w, "END")?;
    w.flush()
}

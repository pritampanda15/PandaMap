//! High-performance PandaMap interaction detection, contact deduplication,
//! empirical scoring, and multi-core parallel trajectory / pose analysis.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::thread;

use crate::geometry::{angle_deg, Vec3};
use crate::spatial_hash::SpatialGrid;
use crate::topology::{
    is_covalent_ligand_element, is_halogen_element, write_snapshot_pdb, PreparedSystem,
};
use crate::xtc::RawXtcFrame;

pub const HBOND_MIN: f64 = 2.5;
pub const HBOND_MAX: f64 = 3.5;
pub const HBOND_DON_ANGLE_MIN: f64 = 100.0;
pub const PI_STACK_CUTOFF: f64 = 5.5;
pub const HYDROPHOBIC_CUTOFF: f64 = 4.0;
pub const IONIC_CUTOFF: f64 = 5.5;
pub const HALOGEN_CUTOFF: f64 = 3.5;
pub const METAL_COORD_CUTOFF: f64 = 2.8;
pub const METAL_LIGAND_CUTOFF: f64 = 3.0;
pub const COVALENT_MIN: f64 = 1.2;
pub const COVALENT_MAX: f64 = 2.1;
pub const REPULSION_CUTOFF: f64 = 4.0;
pub const WATER_BRIDGE_MIN: f64 = 2.5;
pub const WATER_BRIDGE_MAX: f64 = 3.5;
pub const WATER_BRIDGE_ANGLE_MIN: f64 = 71.0;
pub const WATER_BRIDGE_ANGLE_MAX: f64 = 140.0;

pub const INTERACTION_TYPES: [&str; 15] = [
    "hydrogen_bonds",
    "water_bridges",
    "ionic",
    "halogen_bonds",
    "metal_coordination",
    "covalent",
    "pi_pi_stacking",
    "pi_cation",
    "cation_pi",
    "donor_pi",
    "amide_pi",
    "carbon_pi",
    "alkyl_pi",
    "hydrophobic",
    "repulsion",
];

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(usize)]
pub enum InteractionKind {
    HydrogenBonds = 0,
    WaterBridges = 1,
    Ionic = 2,
    HalogenBonds = 3,
    MetalCoordination = 4,
    Covalent = 5,
    PiPiStacking = 6,
    PiCation = 7,
    CationPi = 8,
    DonorPi = 9,
    AmidePi = 10,
    CarbonPi = 11,
    AlkylPi = 12,
    Hydrophobic = 13,
    Repulsion = 14,
}

impl InteractionKind {
    pub const ALL: [InteractionKind; 15] = [
        Self::HydrogenBonds,
        Self::WaterBridges,
        Self::Ionic,
        Self::HalogenBonds,
        Self::MetalCoordination,
        Self::Covalent,
        Self::PiPiStacking,
        Self::PiCation,
        Self::CationPi,
        Self::DonorPi,
        Self::AmidePi,
        Self::CarbonPi,
        Self::AlkylPi,
        Self::Hydrophobic,
        Self::Repulsion,
    ];

    #[inline]
    pub fn as_str(self) -> &'static str {
        INTERACTION_TYPES[self as usize]
    }

    #[inline]
    pub fn is_atom_pair_type(self) -> bool {
        matches!(
            self,
            Self::HydrogenBonds
                | Self::HalogenBonds
                | Self::MetalCoordination
                | Self::Covalent
        )
    }

    #[inline]
    pub fn score_params(self) -> Option<(f64, f64, f64)> {
        // (weight, d_ideal, d_max) matching core.py estimate_binding_affinity
        match self {
            Self::HydrogenBonds => Some((-1.20, 2.9, 3.5)),
            Self::Covalent => Some((-8.00, 1.8, 2.1)),
            Self::MetalCoordination => Some((-2.50, 2.2, 2.8)),
            Self::Ionic => Some((-1.80, 3.5, 5.5)),
            Self::HalogenBonds => Some((-0.80, 3.0, 3.5)),
            Self::PiPiStacking => Some((-0.70, 3.8, 5.5)),
            Self::CationPi => Some((-0.65, 4.0, 5.5)),
            Self::PiCation => Some((-0.65, 4.0, 5.5)),
            Self::Hydrophobic => Some((-0.30, 3.5, 4.0)),
            Self::CarbonPi => Some((-0.25, 4.0, 5.5)),
            Self::DonorPi => Some((-0.20, 3.5, 5.5)),
            Self::AmidePi => Some((-0.20, 3.5, 5.5)),
            Self::AlkylPi => Some((-0.20, 4.0, 5.5)),
            Self::Repulsion => Some((0.60, 3.0, 4.0)),
            Self::WaterBridges => None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct DetectConfig {
    pub strict_halogens: bool,
    pub strict_hbond_geometry: bool,
    pub water_angle_check: bool,
    pub threads: usize,
}

impl Default for DetectConfig {
    fn default() -> Self {
        Self {
            strict_halogens: false,
            strict_hbond_geometry: false,
            water_angle_check: true,
            threads: 0, // 0 => auto-detect available parallelism
        }
    }
}

#[derive(Debug, Clone)]
pub struct RawHit {
    pub prot_atom_idx: usize,
    pub lig_atom_idx: usize,
    pub distance: f64,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ResidueId {
    pub chain: String,
    pub resname: String,
    pub resnum: i32,
}

impl ResidueId {
    pub fn label(&self) -> String {
        let ch = self.chain.trim();
        if ch.is_empty() {
            format!("{}{}", self.resname, self.resnum)
        } else {
            format!("{}{}:{}", self.resname, self.resnum, ch)
        }
    }
}

#[derive(Debug, Clone)]
pub struct FrameAnalysis {
    pub frame_idx: usize,
    pub step: i32,
    pub time_ps: f32,
    pub counts: [usize; 15],
    pub n_interactions: usize,
    pub n_residues: usize,
    pub score: f64,
    /// Per-(ResidueId, InteractionKind) minimum contact distance in this frame
    pub features: Vec<(ResidueId, InteractionKind, f64)>,
}

#[inline]
fn pos_of(coords: &[[f32; 3]], idx: usize) -> Vec3 {
    let c = coords[idx];
    Vec3::new(c[0] as f64, c[1] as f64, c[2] as f64)
}

#[inline]
fn is_neg_charged_res(r: &str) -> bool {
    matches!(r, "ASP" | "GLU")
}

#[inline]
fn is_pos_charged_res(r: &str) -> bool {
    matches!(r, "ARG" | "LYS" | "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP")
}

#[inline]
fn is_amide_res(r: &str) -> bool {
    matches!(r, "ASN" | "GLN")
}

#[inline]
fn is_hydrophobic_res(r: &str) -> bool {
    matches!(
        r,
        "ALA" | "VAL" | "LEU" | "ILE" | "MET" | "PHE" | "TRP" | "PRO" | "TYR"
    )
}

#[inline]
fn is_alkyl_res(r: &str) -> bool {
    matches!(r, "ALA" | "VAL" | "LEU" | "ILE" | "MET" | "PRO")
}

#[inline]
fn is_halogen_acceptor(el: &str, name: &str) -> bool {
    matches!(el, "O" | "N" | "S") || (el == "C" && matches!(name, "CE1" | "CD2"))
}

#[inline]
fn is_metal_coordinator(el: &str, name: &str) -> bool {
    matches!(el, "O" | "N" | "S")
        || matches!(name, "SD" | "OD1" | "OD2" | "OE1" | "OE2" | "NE2" | "ND1")
}

fn best_hbond_angle(
    sys: &PreparedSystem,
    coords: &[[f32; 3]],
    prot_idx: usize,
    lig_idx: usize,
) -> Option<f64> {
    let mut best: Option<f64> = None;
    let pairs = [(prot_idx, lig_idx), (lig_idx, prot_idx)];
    for (donor, acceptor) in pairs {
        let d_pos = pos_of(coords, donor);
        let a_pos = pos_of(coords, acceptor);
        for &h_idx in &sys.attached_h[donor] {
            let h_pos = pos_of(coords, h_idx);
            if let Some(ang) = angle_deg(d_pos, h_pos, a_pos) {
                best = Some(best.map_or(ang, |b| b.max(ang)));
            }
        }
    }
    best
}

fn check_heavy_hbond_geometry(
    sys: &PreparedSystem,
    coords: &[[f32; 3]],
    donor_idx: usize,
    acceptor_idx: usize,
) -> bool {
    let nbs = &sys.heavy_cov_neighbors[donor_idx];
    if nbs.is_empty() {
        return true;
    }
    let d_pos = pos_of(coords, donor_idx);
    let a_pos = pos_of(coords, acceptor_idx);
    let mut avg = Vec3::default();
    for &nb in nbs {
        avg = avg.add(pos_of(coords, nb));
    }
    avg = avg.scale(1.0 / nbs.len() as f64);
    let est_h_dir = d_pos.sub(avg);
    let d_to_a = a_pos.sub(d_pos);
    let n1 = est_h_dir.norm();
    let n2 = d_to_a.norm();
    if n1 < 1e-6 || n2 < 1e-6 {
        return true;
    }
    let cos_a = (est_h_dir.dot(d_to_a) / (n1 * n2)).clamp(-1.0, 1.0);
    cos_a.acos().to_degrees() > 110.0
}

/// Analyze a single frame using 3D spatial hashing + AABB pre-filtering.
pub fn analyze_frame(
    sys: &PreparedSystem,
    coords: &[[f32; 3]],
    frame_idx: usize,
    step: i32,
    time_ps: f32,
    cfg: &DetectConfig,
) -> FrameAnalysis {
    // 1. Ligand bounding box (enables O(1) rejection of distant protein & bulk water atoms)
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut min_z = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    let mut max_z = f64::NEG_INFINITY;
    for &idx in &sys.ligand_indices {
        let p = pos_of(coords, idx);
        min_x = min_x.min(p.x);
        min_y = min_y.min(p.y);
        min_z = min_z.min(p.z);
        max_x = max_x.max(p.x);
        max_y = max_y.max(p.y);
        max_z = max_z.max(p.z);
    }

    // Pad protein grid enough to cover direct contacts (5.5 A), 2-water bridges (10.5 A),
    // and bridging metals (3.0 + 2.8 = 5.8 A)
    let prot_pad = if sys.water_o_indices.is_empty() {
        6.0
    } else {
        11.0
    };
    let mut prot_grid = SpatialGrid::with_capacity(PI_STACK_CUTOFF, 256);
    for &pidx in &sys.protein_indices {
        let p = pos_of(coords, pidx);
        if p.x >= min_x - prot_pad
            && p.x <= max_x + prot_pad
            && p.y >= min_y - prot_pad
            && p.y <= max_y + prot_pad
            && p.z >= min_z - prot_pad
            && p.z <= max_z + prot_pad
        {
            prot_grid.insert(pidx, p);
        }
    }

    let mut raw: [Vec<RawHit>; 15] = Default::default();
    let max_cutoff_sq = PI_STACK_CUTOFF * PI_STACK_CUTOFF;

    // 2. Ligand atom vs protein spatial hash
    for &lidx in &sys.ligand_indices {
        let lpos = pos_of(coords, lidx);
        let latom = &sys.atoms[lidx];
        let lel = latom.element.as_str();
        let is_lig_pos = sys.is_lig_cation[lidx];
        let is_lig_neg = sys.is_lig_anion[lidx];
        let is_lig_arom = sys.is_lig_aromatic[lidx] && lel == "C";
        let is_lig_cov_el = is_covalent_ligand_element(lel);

        prot_grid.query_radius(lpos, PI_STACK_CUTOFF, |pidx| {
            let ppos = pos_of(coords, pidx);
            let d2 = lpos.dist_sq(ppos);
            if d2 > max_cutoff_sq {
                return;
            }
            let d = d2.sqrt();
            if d < 1.5 {
                // Still allow covalent bonds in [1.2, 1.5)
                if d >= COVALENT_MIN && d <= COVALENT_MAX && is_lig_cov_el {
                    let patom = &sys.atoms[pidx];
                    if (patom.resname == "CYS" && patom.name == "SG")
                        || (patom.resname == "SER" && patom.name == "OG")
                        || (patom.resname == "LYS" && patom.name == "NZ")
                        || (matches!(patom.resname.as_str(), "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP")
                            && matches!(patom.name.as_str(), "ND1" | "NE2"))
                    {
                        raw[InteractionKind::Covalent as usize].push(RawHit {
                            prot_atom_idx: pidx,
                            lig_atom_idx: lidx,
                            distance: d,
                        });
                    }
                }
                return;
            }

            let patom = &sys.atoms[pidx];
            let pel = patom.element.as_str();
            let pres = patom.resname.as_str();
            let is_prot_ring = sys.is_prot_ring_atom[pidx];

            // 1. Hydrogen bonds
            if d >= HBOND_MIN && d <= HBOND_MAX && matches!(lel, "N" | "O") && matches!(pel, "N" | "O") {
                let mut ok = true;
                if let Some(ang) = best_hbond_angle(sys, coords, pidx, lidx) {
                    if ang < HBOND_DON_ANGLE_MIN {
                        ok = false;
                    }
                } else if cfg.strict_hbond_geometry {
                    if !(check_heavy_hbond_geometry(sys, coords, pidx, lidx)
                        || check_heavy_hbond_geometry(sys, coords, lidx, pidx))
                    {
                        ok = false;
                    }
                }
                if ok {
                    raw[InteractionKind::HydrogenBonds as usize].push(RawHit {
                        prot_atom_idx: pidx,
                        lig_atom_idx: lidx,
                        distance: d,
                    });
                }
            }

            // 2. Pi-cation (protein cation -> ligand aromatic)
            if d <= PI_STACK_CUTOFF && is_pos_charged_res(pres) && pel == "N" && is_lig_arom && lel == "C" {
                raw[InteractionKind::PiCation as usize].push(RawHit {
                    prot_atom_idx: pidx,
                    lig_atom_idx: lidx,
                    distance: d,
                });
            }

            // 3. Alkyl-Pi
            if d <= PI_STACK_CUTOFF {
                if is_alkyl_res(pres) && pel == "C" && is_lig_arom {
                    raw[InteractionKind::AlkylPi as usize].push(RawHit {
                        prot_atom_idx: pidx,
                        lig_atom_idx: lidx,
                        distance: d,
                    });
                }
                if is_prot_ring && pel == "C" && lel == "C" && !is_lig_arom {
                    raw[InteractionKind::AlkylPi as usize].push(RawHit {
                        prot_atom_idx: pidx,
                        lig_atom_idx: lidx,
                        distance: d,
                    });
                }
            }

            // 4. Repulsion
            if d <= REPULSION_CUTOFF {
                let is_prot_pos = is_pos_charged_res(pres) && pel == "N";
                let is_prot_neg = is_neg_charged_res(pres) && pel == "O";
                if (is_prot_pos && is_lig_pos) || (is_prot_neg && is_lig_neg) {
                    raw[InteractionKind::Repulsion as usize].push(RawHit {
                        prot_atom_idx: pidx,
                        lig_atom_idx: lidx,
                        distance: d,
                    });
                }
            }

            // 5. Pi-Pi stacking
            if d <= PI_STACK_CUTOFF && is_prot_ring && is_lig_arom && lel == "C" && pel == "C" {
                raw[InteractionKind::PiPiStacking as usize].push(RawHit {
                    prot_atom_idx: pidx,
                    lig_atom_idx: lidx,
                    distance: d,
                });
            }

            // 6. Carbon-Pi
            if d <= PI_STACK_CUTOFF && is_prot_ring && lel == "C" {
                raw[InteractionKind::CarbonPi as usize].push(RawHit {
                    prot_atom_idx: pidx,
                    lig_atom_idx: lidx,
                    distance: d,
                });
            }

            // 7. Donor-Pi
            if d <= PI_STACK_CUTOFF && is_neg_charged_res(pres) && is_lig_arom {
                raw[InteractionKind::DonorPi as usize].push(RawHit {
                    prot_atom_idx: pidx,
                    lig_atom_idx: lidx,
                    distance: d,
                });
            }

            // 8. Amide-Pi
            if d <= PI_STACK_CUTOFF && is_amide_res(pres) && is_lig_arom {
                raw[InteractionKind::AmidePi as usize].push(RawHit {
                    prot_atom_idx: pidx,
                    lig_atom_idx: lidx,
                    distance: d,
                });
            }

            // 9. Hydrophobic
            if d <= HYDROPHOBIC_CUTOFF && is_hydrophobic_res(pres) && lel == "C" && pel == "C" {
                raw[InteractionKind::Hydrophobic as usize].push(RawHit {
                    prot_atom_idx: pidx,
                    lig_atom_idx: lidx,
                    distance: d,
                });
            }

            // 10. Ionic / salt bridge
            if d <= IONIC_CUTOFF {
                let is_prot_pos = is_pos_charged_res(pres) && pel == "N";
                let is_prot_neg = is_neg_charged_res(pres) && pel == "O";
                if (is_prot_pos && is_lig_neg) || (is_prot_neg && is_lig_pos) {
                    raw[InteractionKind::Ionic as usize].push(RawHit {
                        prot_atom_idx: pidx,
                        lig_atom_idx: lidx,
                        distance: d,
                    });
                }
            }

            // 11. Cation-Pi (ligand cation -> protein ring)
            if d <= PI_STACK_CUTOFF && is_prot_ring && is_lig_pos {
                raw[InteractionKind::CationPi as usize].push(RawHit {
                    prot_atom_idx: pidx,
                    lig_atom_idx: lidx,
                    distance: d,
                });
            }

            // 12. Covalent bonds
            if d >= COVALENT_MIN && d <= COVALENT_MAX && is_lig_cov_el {
                if (pres == "CYS" && patom.name == "SG")
                    || (pres == "SER" && patom.name == "OG")
                    || (pres == "LYS" && patom.name == "NZ")
                    || (matches!(pres, "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP")
                        && matches!(patom.name.as_str(), "ND1" | "NE2"))
                {
                    raw[InteractionKind::Covalent as usize].push(RawHit {
                        prot_atom_idx: pidx,
                        lig_atom_idx: lidx,
                        distance: d,
                    });
                }
            }
        });
    }

    // 3. Halogen bonds
    for &hidx in &sys.halogen_indices {
        let hel = sys.atoms[hidx].element.as_str();
        if !is_halogen_element(hel, cfg.strict_halogens) {
            continue;
        }
        let hpos = pos_of(coords, hidx);
        prot_grid.query_radius(hpos, HALOGEN_CUTOFF, |pidx| {
            let ppos = pos_of(coords, pidx);
            let d = hpos.dist(ppos);
            let patom = &sys.atoms[pidx];
            if d > 2.5 && d <= HALOGEN_CUTOFF && is_halogen_acceptor(&patom.element, &patom.name) {
                raw[InteractionKind::HalogenBonds as usize].push(RawHit {
                    prot_atom_idx: pidx,
                    lig_atom_idx: hidx,
                    distance: d,
                });
            }
        });
    }

    // 4. Metal coordination
    for &midx in &sys.metal_indices {
        let mpos = pos_of(coords, midx);
        let min_lig_d = sys
            .ligand_indices
            .iter()
            .map(|&lidx| mpos.dist(pos_of(coords, lidx)))
            .fold(f64::INFINITY, f64::min);
        if min_lig_d > METAL_LIGAND_CUTOFF {
            continue;
        }
        prot_grid.query_radius(mpos, METAL_COORD_CUTOFF, |pidx| {
            let ppos = pos_of(coords, pidx);
            let d = mpos.dist(ppos);
            let patom = &sys.atoms[pidx];
            if d <= METAL_COORD_CUTOFF && is_metal_coordinator(&patom.element, &patom.name) {
                raw[InteractionKind::MetalCoordination as usize].push(RawHit {
                    prot_atom_idx: pidx,
                    lig_atom_idx: midx,
                    distance: d,
                });
            }
        });
    }

    // 5. Water bridges (1-water and 2-water) using AABB-filtered water spatial hash
    if !sys.water_o_indices.is_empty() {
        let w_pad = 2.0 * WATER_BRIDGE_MAX;
        let mut active_waters: Vec<usize> = Vec::new();
        let mut water_grid = SpatialGrid::with_capacity(WATER_BRIDGE_MAX, 64);
        for &widx in &sys.water_o_indices {
            let wp = pos_of(coords, widx);
            if wp.x >= min_x - w_pad
                && wp.x <= max_x + w_pad
                && wp.y >= min_y - w_pad
                && wp.y <= max_y + w_pad
                && wp.z >= min_z - w_pad
                && wp.z <= max_z + w_pad
            {
                active_waters.push(widx);
                water_grid.insert(widx, wp);
            }
        }

        let lig_no: Vec<(usize, Vec3)> = sys
            .ligand_indices
            .iter()
            .copied()
            .filter(|&i| matches!(sys.atoms[i].element.as_str(), "N" | "O"))
            .map(|i| (i, pos_of(coords, i)))
            .collect();

        let angle_ok = |a: Vec3, c: Vec3, b: Vec3| -> bool {
            if !cfg.water_angle_check {
                return true;
            }
            match angle_deg(a, c, b) {
                Some(omega) => (WATER_BRIDGE_ANGLE_MIN..=WATER_BRIDGE_ANGLE_MAX).contains(&omega),
                None => false,
            }
        };

        for &w1 in &active_waters {
            let w1p = pos_of(coords, w1);
            let mut lig_side = Vec::new();
            for &(lidx, lpos) in &lig_no {
                let d = w1p.dist(lpos);
                if (WATER_BRIDGE_MIN..=WATER_BRIDGE_MAX).contains(&d) {
                    lig_side.push((lidx, lpos, d));
                }
            }
            if lig_side.is_empty() {
                continue;
            }

            // 1-water bridges: lig ... w1 ... prot
            for &(lidx, lpos, d_lig) in &lig_side {
                prot_grid.query_radius(w1p, WATER_BRIDGE_MAX, |pidx| {
                    let pel = sys.atoms[pidx].element.as_str();
                    if !matches!(pel, "N" | "O") {
                        return;
                    }
                    let ppos = pos_of(coords, pidx);
                    let d_prot = w1p.dist(ppos);
                    if (WATER_BRIDGE_MIN..=WATER_BRIDGE_MAX).contains(&d_prot)
                        && angle_ok(lpos, w1p, ppos)
                    {
                        raw[InteractionKind::WaterBridges as usize].push(RawHit {
                            prot_atom_idx: pidx,
                            lig_atom_idx: lidx,
                            distance: d_lig.max(d_prot),
                        });
                    }
                });

                // 2-water bridges: lig ... w1 ... w2 ... prot
                water_grid.query_radius(w1p, WATER_BRIDGE_MAX, |w2| {
                    if w2 == w1 {
                        return;
                    }
                    let w2p = pos_of(coords, w2);
                    let d_ww = w1p.dist(w2p);
                    if !(WATER_BRIDGE_MIN..=WATER_BRIDGE_MAX).contains(&d_ww)
                        || !angle_ok(lpos, w1p, w2p)
                    {
                        return;
                    }
                    prot_grid.query_radius(w2p, WATER_BRIDGE_MAX, |pidx| {
                        let pel = sys.atoms[pidx].element.as_str();
                        if !matches!(pel, "N" | "O") {
                            return;
                        }
                        let ppos = pos_of(coords, pidx);
                        let d_prot = w2p.dist(ppos);
                        if (WATER_BRIDGE_MIN..=WATER_BRIDGE_MAX).contains(&d_prot)
                            && angle_ok(w1p, w2p, ppos)
                        {
                            raw[InteractionKind::WaterBridges as usize].push(RawHit {
                                prot_atom_idx: pidx,
                                lig_atom_idx: lidx,
                                distance: d_lig.max(d_ww).max(d_prot),
                            });
                        }
                    });
                });
            }
        }
    }

    // 6. Exact _select_contacts deduplication matching core.py
    let mut selected: [Vec<RawHit>; 15] = Default::default();
    for kind in InteractionKind::ALL {
        let hits = &raw[kind as usize];
        if hits.is_empty() {
            continue;
        }
        let mut by_res: HashMap<(&str, &str, i32), Vec<RawHit>> = HashMap::new();
        for h in hits {
            let a = &sys.atoms[h.prot_atom_idx];
            by_res
                .entry((a.chain.as_str(), a.resname.as_str(), a.resnum))
                .or_default()
                .push(h.clone());
        }
        let out = &mut selected[kind as usize];
        for mut group in by_res.into_values() {
            group.sort_by(|a, b| a.distance.total_cmp(&b.distance));
            if kind.is_atom_pair_type() {
                let mut used_prot = HashSet::new();
                let mut used_lig = HashSet::new();
                for hit in group {
                    if used_prot.contains(&hit.prot_atom_idx) && used_lig.contains(&hit.lig_atom_idx) {
                        continue;
                    }
                    used_prot.insert(hit.prot_atom_idx);
                    used_lig.insert(hit.lig_atom_idx);
                    out.push(hit);
                }
            } else {
                out.push(group.swap_remove(0));
            }
        }
    }

    // Carbon-contact priority deduplication: pi_pi_stacking > alkyl_pi > carbon_pi > hydrophobic
    let carbon_order = [
        InteractionKind::PiPiStacking,
        InteractionKind::AlkylPi,
        InteractionKind::CarbonPi,
        InteractionKind::Hydrophobic,
    ];
    let mut seen_pairs: HashSet<(usize, usize)> = HashSet::new();
    for kind in carbon_order {
        let mut unique = Vec::new();
        for hit in std::mem::take(&mut selected[kind as usize]) {
            if seen_pairs.insert((hit.prot_atom_idx, hit.lig_atom_idx)) {
                unique.push(hit);
            }
        }
        selected[kind as usize] = unique;
    }

    // 7. Summarize counts, empirical score, and per-(residue, kind) feature map
    let mut counts = [0usize; 15];
    let mut n_interactions = 0usize;
    let mut unique_residues: HashSet<ResidueId> = HashSet::new();
    let mut features_map: BTreeMap<(ResidueId, InteractionKind), f64> = BTreeMap::new();

    let mut score = 2.9 + (sys.n_rotatable_bonds as f64) * 0.5;

    for kind in InteractionKind::ALL {
        let hits = &selected[kind as usize];
        counts[kind as usize] = hits.len();
        n_interactions += hits.len();
        if hits.is_empty() {
            continue;
        }
        let mut best_per_res: HashMap<ResidueId, f64> = HashMap::new();
        for hit in hits {
            let a = &sys.atoms[hit.prot_atom_idx];
            let rid = ResidueId {
                chain: a.chain.clone(),
                resname: a.resname.clone(),
                resnum: a.resnum,
            };
            unique_residues.insert(rid.clone());
            best_per_res
                .entry(rid.clone())
                .and_modify(|d| {
                    if hit.distance < *d {
                        *d = hit.distance;
                    }
                })
                .or_insert(hit.distance);
        }
        for (rid, best_d) in &best_per_res {
            features_map.insert((rid.clone(), kind), *best_d);
        }
        if let Some((weight, d_ideal, d_max)) = kind.score_params() {
            for &d in best_per_res.values() {
                let decay = if d <= d_ideal {
                    1.0
                } else {
                    ((d_max - d) / (d_max - d_ideal).max(1e-6)).max(0.0)
                };
                score += weight * decay;
            }
        }
    }

    let features: Vec<(ResidueId, InteractionKind, f64)> = features_map
        .into_iter()
        .map(|((rid, kind), d)| (rid, kind, d))
        .collect();

    FrameAnalysis {
        frame_idx,
        step,
        time_ps,
        counts,
        n_interactions,
        n_residues: unique_residues.len(),
        score: (score * 100.0).round() / 100.0,
        features,
    }
}

/// Parallel batch processing of raw XTC frames using `std::thread::scope`.
pub fn analyze_xtc_batch_parallel(
    sys: &PreparedSystem,
    batch: &[RawXtcFrame],
    cfg: &DetectConfig,
) -> io::Result<Vec<FrameAnalysis>> {
    if batch.is_empty() {
        return Ok(Vec::new());
    }
    let n_threads = if cfg.threads > 0 {
        cfg.threads
    } else {
        thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
    }
    .min(batch.len())
    .max(1);

    if n_threads == 1 {
        let mut out = Vec::with_capacity(batch.len());
        for raw in batch {
            let decoded = raw.decode()?;
            out.push(analyze_frame(
                sys,
                &decoded.coords_angstrom,
                decoded.frame_idx,
                decoded.step,
                decoded.time_ps,
                cfg,
            ));
        }
        return Ok(out);
    }

    let chunk_size = batch.len().div_ceil(n_threads);
    let results: Vec<io::Result<Vec<FrameAnalysis>>> = thread::scope(|scope| {
        let handles: Vec<_> = batch
            .chunks(chunk_size)
            .map(|chunk| {
                scope.spawn(move || {
                    let mut local = Vec::with_capacity(chunk.len());
                    for raw in chunk {
                        let decoded = raw.decode()?;
                        local.push(analyze_frame(
                            sys,
                            &decoded.coords_angstrom,
                            decoded.frame_idx,
                            decoded.step,
                            decoded.time_ps,
                            cfg,
                        ));
                    }
                    Ok(local)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });

    let mut out = Vec::with_capacity(batch.len());
    for part in results {
        out.extend(part?);
    }
    Ok(out)
}

/// Parallel analysis of in-memory coordinate frames (multi-MODEL PDB / multi-pose PDBQT).
pub fn analyze_models_parallel(
    sys: &PreparedSystem,
    models: &[Vec<[f32; 3]>],
    dt_ps: f32,
    cfg: &DetectConfig,
) -> Vec<FrameAnalysis> {
    if models.is_empty() {
        return Vec::new();
    }
    let n_threads = if cfg.threads > 0 {
        cfg.threads
    } else {
        thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
    }
    .min(models.len())
    .max(1);

    let chunk_size = models.len().div_ceil(n_threads);
    let parts: Vec<Vec<FrameAnalysis>> = thread::scope(|scope| {
        let handles: Vec<_> = models
            .chunks(chunk_size)
            .enumerate()
            .map(|(chunk_idx, chunk)| {
                let base = chunk_idx * chunk_size;
                scope.spawn(move || {
                    let mut local = Vec::with_capacity(chunk.len());
                    for (offset, coords) in chunk.iter().enumerate() {
                        let fidx = base + offset;
                        local.push(analyze_frame(
                            sys,
                            coords,
                            fidx,
                            fidx as i32,
                            fidx as f32 * dt_ps,
                            cfg,
                        ));
                    }
                    local
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });

    let mut out = Vec::with_capacity(models.len());
    for p in parts {
        out.extend(p);
    }
    out
}

#[derive(Debug, Clone)]
pub struct SnapshotSpec {
    pub every_ps: Option<f64>,
    pub stride: Option<usize>,
    pub start_ps: Option<f64>,
    pub end_ps: Option<f64>,
    pub water_shell_a: f64,
}

#[derive(Debug, Clone)]
pub struct ExportedSnapshot {
    pub frame_idx: usize,
    pub step: i32,
    pub time_ps: f32,
    pub label: String,
    pub pdb_path: PathBuf,
    pub score: f64,
    pub n_interactions: usize,
}

/// Choose which frame indices to export as snapshots given time interval or stride.
pub fn select_snapshot_indices(
    frames: &[FrameAnalysis],
    spec: &SnapshotSpec,
) -> Vec<usize> {
    if frames.is_empty() {
        return Vec::new();
    }
    let in_window: Vec<usize> = frames
        .iter()
        .enumerate()
        .filter(|(_, f)| {
            let t = f.time_ps as f64;
            spec.start_ps.map_or(true, |s| t >= s - 1e-3)
                && spec.end_ps.map_or(true, |e| t <= e + 1e-3)
        })
        .map(|(i, _)| i)
        .collect();
    if in_window.is_empty() {
        return Vec::new();
    }

    if let Some(every_ps) = spec.every_ps {
        if every_ps > 0.0 {
            let t0 = frames[in_window[0]].time_ps as f64;
            let t_end = frames[*in_window.last().unwrap()].time_ps as f64;
            // If timestamps are all identical (e.g. docking poses without time), fall back to stride
            if (t_end - t0).abs() < 1e-6 {
                return vec![in_window[0]];
            }
            let mut chosen = Vec::new();
            let mut seen = HashSet::new();
            let mut target = spec.start_ps.unwrap_or(t0);
            if target < t0 {
                target = t0;
            }
            let mut cursor = 0usize;
            while target <= t_end + 1e-3 {
                while cursor + 1 < in_window.len()
                    && (frames[in_window[cursor + 1]].time_ps as f64 - target).abs()
                        <= (frames[in_window[cursor]].time_ps as f64 - target).abs()
                {
                    cursor += 1;
                }
                let best_i = in_window[cursor];
                if (frames[best_i].time_ps as f64 - target).abs() <= every_ps * 0.55
                    && seen.insert(best_i)
                {
                    chosen.push(best_i);
                }
                target += every_ps;
            }
            return chosen;
        }
    }

    let stride = spec.stride.unwrap_or(1).max(1);
    in_window.into_iter().step_by(stride).collect()
}

pub fn format_time_label(time_ps: f32, frame_idx: usize, use_time: bool) -> String {
    if !use_time {
        return format!("frame_{:05}", frame_idx);
    }
    let ns = time_ps as f64 / 1000.0;
    if (ns - ns.round()).abs() < 1e-3 {
        format!("t_{:05}ns", ns.round() as i64)
    } else if ((ns * 10.0) - (ns * 10.0).round()).abs() < 1e-3 {
        format!("t_{:06.1}ns", ns)
    } else {
        format!("t_{:07.0}ps", time_ps.round())
    }
}

/// Write summary CSV, residue occupancy CSV, ML fingerprint CSV, and manifest JSON.
pub fn write_analysis_outputs(
    out_dir: &Path,
    sys: &PreparedSystem,
    frames: &[FrameAnalysis],
    snapshots: &[ExportedSnapshot],
    elapsed_sec: f64,
    write_fingerprints: bool,
) -> io::Result<()> {
    fs::create_dir_all(out_dir)?;

    // 1. trajectory_summary.csv
    let summary_path = out_dir.join("trajectory_summary.csv");
    {
        let mut w = BufWriter::new(File::create(&summary_path)?);
        write!(w, "frame_idx,step,time_ps,time_ns,score,n_interactions,n_residues")?;
        for t in INTERACTION_TYPES {
            write!(w, ",{}", t)?;
        }
        writeln!(w)?;
        for f in frames {
            write!(
                w,
                "{},{},{:.3},{:.6},{:.2},{},{}",
                f.frame_idx,
                f.step,
                f.time_ps,
                f.time_ps as f64 / 1000.0,
                f.score,
                f.n_interactions,
                f.n_residues
            )?;
            for c in f.counts {
                write!(w, ",{}", c)?;
            }
            writeln!(w)?;
        }
        w.flush()?;
    }

    // 2. Aggregate per-residue & per-(residue, type) occupancies
    let total_f = frames.len().max(1) as f64;
    let mut res_frame_count: BTreeMap<ResidueId, usize> = BTreeMap::new();
    let mut res_dist_sum: BTreeMap<ResidueId, f64> = BTreeMap::new();
    let mut res_type_counts: BTreeMap<ResidueId, [usize; 15]> = BTreeMap::new();
    let mut active_features: BTreeSet<(ResidueId, InteractionKind)> = BTreeSet::new();

    for f in frames {
        let mut frame_res_min: HashMap<&ResidueId, f64> = HashMap::new();
        for (rid, kind, dist) in &f.features {
            active_features.insert((rid.clone(), *kind));
            res_type_counts.entry(rid.clone()).or_insert([0; 15])[*kind as usize] += 1;
            frame_res_min
                .entry(rid)
                .and_modify(|d| {
                    if *dist < *d {
                        *d = *dist;
                    }
                })
                .or_insert(*dist);
        }
        for (rid, min_d) in frame_res_min {
            *res_frame_count.entry(rid.clone()).or_insert(0) += 1;
            *res_dist_sum.entry(rid.clone()).or_insert(0.0) += min_d;
        }
    }

    // Sort residues by descending occupancy
    let mut ranked: Vec<(ResidueId, usize)> = res_frame_count.into_iter().collect();
    ranked.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    let occ_path = out_dir.join("residue_occupancy.csv");
    {
        let mut w = BufWriter::new(File::create(&occ_path)?);
        write!(
            w,
            "chain,residue,resnum,label,n_frames,occupancy_pct,mean_min_dist_A,dominant_types"
        )?;
        for t in INTERACTION_TYPES {
            write!(w, ",{}_pct", t)?;
        }
        writeln!(w)?;
        for (rid, cnt) in &ranked {
            let pct = (*cnt as f64) * 100.0 / total_f;
            let mean_d = res_dist_sum[rid] / (*cnt as f64);
            let tcounts = res_type_counts[rid];
            let mut dom: Vec<String> = InteractionKind::ALL
                .iter()
                .filter(|&&k| tcounts[k as usize] > 0)
                .map(|&k| {
                    format!(
                        "{}({:.1}%)",
                        k.as_str(),
                        tcounts[k as usize] as f64 * 100.0 / total_f
                    )
                })
                .collect();
            dom.truncate(4);
            write!(
                w,
                "{},{},{},{},{},{:.2},{:.3},\"{}\"",
                rid.chain.trim(),
                rid.resname,
                rid.resnum,
                rid.label(),
                cnt,
                pct,
                mean_d,
                dom.join("; ")
            )?;
            for c in tcounts {
                write!(w, ",{:.2}", c as f64 * 100.0 / total_f)?;
            }
            writeln!(w)?;
        }
        w.flush()?;
    }

    // 3. ML-ready trajectory_fingerprints.csv
    if write_fingerprints {
        let fp_path = out_dir.join("trajectory_fingerprints.csv");
        let cols: Vec<(ResidueId, InteractionKind)> = active_features.into_iter().collect();
        let col_index: HashMap<&(ResidueId, InteractionKind), usize> =
            cols.iter().enumerate().map(|(i, c)| (c, i)).collect();
        let mut w = BufWriter::new(File::create(&fp_path)?);
        write!(w, "frame_idx,time_ns")?;
        for (rid, kind) in &cols {
            write!(w, ",{}:{}", rid.label(), kind.as_str())?;
        }
        writeln!(w)?;
        let mut row = vec![0u8; cols.len()];
        for f in frames {
            row.fill(0);
            for (rid, kind, _) in &f.features {
                if let Some(&ci) = col_index.get(&(rid.clone(), *kind)) {
                    row[ci] = 1;
                }
            }
            write!(w, "{},{:.6}", f.frame_idx, f.time_ps as f64 / 1000.0)?;
            for &bit in &row {
                write!(w, ",{}", bit)?;
            }
            writeln!(w)?;
        }
        w.flush()?;
    }

    // 4. engine_manifest.json (zero-dependency JSON writer for Python integration)
    let manifest_path = out_dir.join("engine_manifest.json");
    {
        let mut w = BufWriter::new(File::create(&manifest_path)?);
        writeln!(w, "{{")?;
        writeln!(w, "  \"engine\": \"pandamap-engine-rust\",")?;
        writeln!(w, "  \"n_frames\": {},", frames.len())?;
        writeln!(w, "  \"elapsed_sec\": {:.6},", elapsed_sec)?;
        writeln!(
            w,
            "  \"frames_per_sec\": {:.1},",
            if elapsed_sec > 0.0 {
                frames.len() as f64 / elapsed_sec
            } else {
                0.0
            }
        )?;
        writeln!(
            w,
            "  \"ligand\": {{\"name\": \"{}\", \"chain\": \"{}\", \"resnum\": {}, \"n_atoms\": {}}},",
            json_escape(&sys.ligand_name),
            json_escape(sys.ligand_chain.trim()),
            sys.ligand_resnum,
            sys.ligand_indices.len()
        )?;
        writeln!(
            w,
            "  \"summary_csv\": \"{}\",",
            json_escape(&summary_path.to_string_lossy())
        )?;
        writeln!(
            w,
            "  \"occupancy_csv\": \"{}\",",
            json_escape(&occ_path.to_string_lossy())
        )?;
        if write_fingerprints {
            writeln!(
                w,
                "  \"fingerprints_csv\": \"{}\",",
                json_escape(&out_dir.join("trajectory_fingerprints.csv").to_string_lossy())
            )?;
        } else {
            writeln!(w, "  \"fingerprints_csv\": null,")?;
        }
        writeln!(w, "  \"snapshots\": [")?;
        for (i, s) in snapshots.iter().enumerate() {
            let comma = if i + 1 < snapshots.len() { "," } else { "" };
            writeln!(
                w,
                "    {{\"frame_idx\": {}, \"step\": {}, \"time_ps\": {:.3}, \"time_ns\": {:.6}, \"label\": \"{}\", \"pdb_path\": \"{}\", \"score\": {:.2}, \"n_interactions\": {}}}{}",
                s.frame_idx,
                s.step,
                s.time_ps,
                s.time_ps as f64 / 1000.0,
                json_escape(&s.label),
                json_escape(&s.pdb_path.to_string_lossy()),
                s.score,
                s.n_interactions,
                comma
            )?;
        }
        writeln!(w, "  ]")?;
        writeln!(w, "}}")?;
        w.flush()?;
    }

    Ok(())
}

fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c => out.push(c),
        }
    }
    out
}

/// Export selected snapshot PDBs into `<out_dir>/snapshots/<label>/<label>.pdb`.
pub fn export_snapshot_pdbs(
    out_dir: &Path,
    sys: &PreparedSystem,
    selected: &[(usize, &FrameAnalysis, Vec<[f32; 3]>)],
    water_shell_a: f64,
    use_time_labels: bool,
) -> io::Result<Vec<ExportedSnapshot>> {
    let snap_root = out_dir.join("snapshots");
    if !selected.is_empty() {
        fs::create_dir_all(&snap_root)?;
    }
    let mut exported = Vec::with_capacity(selected.len());
    for &(_idx, fa, ref coords) in selected {
        let label = format_time_label(fa.time_ps, fa.frame_idx, use_time_labels);
        let dir = snap_root.join(&label);
        fs::create_dir_all(&dir)?;
        let pdb_path = dir.join(format!("{}.pdb", label));
        write_snapshot_pdb(sys, coords, &pdb_path, water_shell_a, fa.time_ps)?;
        exported.push(ExportedSnapshot {
            frame_idx: fa.frame_idx,
            step: fa.step,
            time_ps: fa.time_ps,
            label,
            pdb_path,
            score: fa.score,
            n_interactions: fa.n_interactions,
        });
    }
    Ok(exported)
}

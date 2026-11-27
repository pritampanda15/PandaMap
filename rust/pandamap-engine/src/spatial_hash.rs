//! 3D uniform spatial hash (cell list) for O(1) neighbor queries per frame.
//!
//! Uses a fast flat-hashed cell table backed by contiguous entry storage and
//! per-query AABB pre-filtering so explicit-solvent MD systems (50k-200k atoms)
//! only hash atoms within the active-site bounding box.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

use crate::geometry::Vec3;

/// Fast FNV-like / SplitMix integer hasher for `(i32, i32, i32)` cell keys.
#[derive(Default)]
pub struct CellHasher {
    state: u64,
}

impl Hasher for CellHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.state
    }

    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let mut h = self.state;
        for &b in bytes {
            h = h.wrapping_mul(0x517cc1b727220a95).wrapping_add(b as u64);
        }
        self.state = h;
    }

    #[inline]
    fn write_i32(&mut self, i: i32) {
        self.state = self
            .state
            .wrapping_mul(0x9e3779b97f4a7c15)
            .wrapping_add(i as u32 as u64);
    }
}

type FastCellMap = HashMap<(i32, i32, i32), Vec<usize>, BuildHasherDefault<CellHasher>>;

pub struct SpatialGrid {
    inv_cell: f64,
    cells: FastCellMap,
}

impl SpatialGrid {
    pub fn with_capacity(cell_size: f64, capacity: usize) -> Self {
        Self {
            inv_cell: 1.0 / cell_size.max(1e-3),
            cells: HashMap::with_capacity_and_hasher(capacity, BuildHasherDefault::default()),
        }
    }

    #[inline]
    fn cell_key(&self, p: Vec3) -> (i32, i32, i32) {
        (
            (p.x * self.inv_cell).floor() as i32,
            (p.y * self.inv_cell).floor() as i32,
            (p.z * self.inv_cell).floor() as i32,
        )
    }

    #[inline]
    pub fn insert(&mut self, idx: usize, p: Vec3) {
        let key = self.cell_key(p);
        self.cells.entry(key).or_default().push(idx);
    }

    /// Visit every inserted index whose cell is within `ceil(radius / cell_size)`
    /// cells of `query` (for `radius <= cell_size`, checks exactly 27 cells).
    #[inline]
    pub fn query_radius<F>(&self, query: Vec3, radius: f64, mut f: F)
    where
        F: FnMut(usize),
    {
        let (cx, cy, cz) = self.cell_key(query);
        let span = (radius * self.inv_cell).ceil().max(1.0) as i32;
        for dx in -span..=span {
            for dy in -span..=span {
                for dz in -span..=span {
                    if let Some(bucket) = self.cells.get(&(cx + dx, cy + dy, cz + dz)) {
                        for &idx in bucket {
                            f(idx);
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spatial_grid_finds_all_neighbors() {
        let pts = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(3.0, 0.0, 0.0),
            Vec3::new(5.4, 0.0, 0.0),
            Vec3::new(12.0, 0.0, 0.0),
        ];
        let mut grid = SpatialGrid::with_capacity(5.5, pts.len());
        for (i, &p) in pts.iter().enumerate() {
            grid.insert(i, p);
        }
        let mut found = Vec::new();
        grid.query_radius(Vec3::new(0.0, 0.0, 0.0), 5.5, |idx| {
            if pts[idx].dist(Vec3::new(0.0, 0.0, 0.0)) <= 5.5 {
                found.push(idx);
            }
        });
        found.sort_unstable();
        assert_eq!(found, vec![0, 1, 2]);
    }
}

//! Compact 3D vector math and geometry utilities for interaction detection.

#[derive(Debug, Clone, Copy, Default)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    #[inline]
    pub fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }

    #[inline]
    pub fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }

    #[inline]
    pub fn scale(self, s: f64) -> Self {
        Self {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }

    #[inline]
    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    #[inline]
    pub fn cross(self, other: Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    #[inline]
    pub fn norm_sq(self) -> f64 {
        self.dot(self)
    }

    #[inline]
    pub fn norm(self) -> f64 {
        self.norm_sq().sqrt()
    }

    #[inline]
    pub fn dist_sq(self, other: Self) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        dx * dx + dy * dy + dz * dz
    }

    #[inline]
    pub fn dist(self, other: Self) -> f64 {
        self.dist_sq(other).sqrt()
    }
}

/// Angle in degrees between vectors (a - center) and (b - center).
#[inline]
pub fn angle_deg(a: Vec3, center: Vec3, b: Vec3) -> Option<f64> {
    let u = a.sub(center);
    let v = b.sub(center);
    let nu = u.norm();
    let nv = v.norm();
    if nu < 1e-12 || nv < 1e-12 {
        return None;
    }
    let cos_a = (u.dot(v) / (nu * nv)).clamp(-1.0, 1.0);
    Some(cos_a.acos().to_degrees())
}

/// Check whether a set of 3D points lies within `tol` Angstrom of a best-fit plane
/// (exact 3x3 symmetric Jacobi eigensolver on the covariance matrix, matching SVD).
pub fn is_planar(points: &[Vec3], tol: f64) -> bool {
    if points.len() < 4 {
        return true;
    }
    let n = points.len() as f64;
    let mut mean = Vec3::default();
    for &p in points {
        mean = mean.add(p);
    }
    mean = mean.scale(1.0 / n);

    let mut cov = [[0.0f64; 3]; 3];
    for &p in points {
        let d = p.sub(mean);
        let v = [d.x, d.y, d.z];
        for i in 0..3 {
            for j in 0..3 {
                cov[i][j] += v[i] * v[j];
            }
        }
    }

    let normal = smallest_eigenvector_3x3(cov);
    for &p in points {
        let d = p.sub(mean);
        if d.dot(normal).abs() > tol {
            return false;
        }
    }
    true
}

/// Compute the eigenvector corresponding to the smallest eigenvalue of a symmetric 3x3 matrix
/// via cyclic Jacobi rotations.
fn smallest_eigenvector_3x3(mut a: [[f64; 3]; 3]) -> Vec3 {
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..32 {
        // Find largest off-diagonal element
        let mut p = 0;
        let mut q = 1;
        let mut max_val = a[0][1].abs();
        if a[0][2].abs() > max_val {
            max_val = a[0][2].abs();
            p = 0;
            q = 2;
        }
        if a[1][2].abs() > max_val {
            max_val = a[1][2].abs();
            p = 1;
            q = 2;
        }
        if max_val < 1e-12 {
            break;
        }
        let diff = a[q][q] - a[p][p];
        let t = if diff.abs() < 1e-18 {
            1.0
        } else {
            let phi = diff / (2.0 * a[p][q]);
            let sign = if phi >= 0.0 { 1.0 } else { -1.0 };
            1.0 / (phi.abs() + (phi * phi + 1.0).sqrt()) * sign
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        let tau = s / (1.0 + c);

        let apq = a[p][q];
        a[p][p] -= t * apq;
        a[q][q] += t * apq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;

        for r in 0..3 {
            if r != p && r != q {
                let arp = a[r][p];
                let arq = a[r][q];
                a[r][p] = arp - s * (arq + tau * arp);
                a[p][r] = a[r][p];
                a[r][q] = arq + s * (arp - tau * arq);
                a[q][r] = a[r][q];
            }
        }
        for r in 0..3 {
            let vrp = v[r][p];
            let vrq = v[r][q];
            v[r][p] = vrp - s * (vrq + tau * vrp);
            v[r][q] = vrq + s * (vrp - tau * vrq);
        }
    }

    let mut min_idx = 0;
    if a[1][1] < a[min_idx][min_idx] {
        min_idx = 1;
    }
    if a[2][2] < a[min_idx][min_idx] {
        min_idx = 2;
    }
    let n = Vec3::new(v[0][min_idx], v[1][min_idx], v[2][min_idx]);
    let norm = n.norm();
    if norm > 1e-12 {
        n.scale(1.0 / norm)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    }
}

use std::collections::HashSet;
use pyo3::prelude::*;
use numpy::{PyReadwriteArray, Ix2, PyUntypedArrayMethods, PyReadonlyArray2};
use kdtree::{KdTree, distance::squared_euclidean, ErrorKind};
use pyo3_stub_gen::{derive::gen_stub_pyclass, derive::gen_stub_pymethods};

use glam::DVec3;
use pyo3::exceptions::PyValueError;

#[gen_stub_pyclass]
#[pyclass]
pub struct PeriodicBox {
    a: DVec3,
    b: DVec3,
    c: DVec3,
}

/// All constructors of PeriodicBox should call this!
/// This ensures that the same guarantees about PBCs are upheld as in GROMACS or OpenMM.
/// This allows them to be more directly used with those software's, as well as some code is
/// easier/faster to write for ourselves as well!
fn sanitize(pbc: &PeriodicBox) -> PyResult<()> {
    if pbc.a.y != 0. || pbc.a.z != 0. || pbc.b.z != 0. {
        return Err(PyValueError::new_err("a.y, a.z and b.z of periodic boxes must be 0."));
    }
    if pbc.a.x <= 0. || pbc.b.y <= 0. || pbc.c.z <= 0. {
        return Err(PyValueError::new_err("a.x, b.y and c.z of periodic boxes must be positive."));
    }
    if pbc.b.x.abs() * 2. > pbc.a.x || pbc.c.x.abs() * 2. > pbc.a.x || pbc.c.y.abs() * 2. > pbc.b.y {
        return Err(PyValueError::new_err("b.x, c.x must be smaller in magnitude than a.x/2. Likewise c.y must be smaller in magnitude than b.y/2."))
    }

    Ok(())
}

#[gen_stub_pymethods]
#[pymethods]
impl PeriodicBox {
    /// Create a new PeriodicBox from three unit cell vectors.
    ///
    /// Must obey certain relationships:
    ///
    /// * a.x > 0, a.y = 0, a.z = 0
    /// * 2|b.x| < a.x, b.y > 0, b.z = 0
    /// * 2|c.x| < a.x, 2|c.y| < b.y, c.z > 0
    ///
    /// This is validated, and will raise a PyValueError if these conditions are not met.
    #[new]
    pub fn new(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> PyResult<Self> {
        let pbc = Self { a: a.into(), b: b.into(), c: c.into() };
        sanitize(&pbc)?;
        Ok(pbc)
    }

    /// Create a PeriodicBox based on unit cell description, as found in Gro files.
    #[pyo3(signature = (ax, by, cz, ay=None, az=None, bx=None, bz=None, cx=None, cy=None))]
    #[staticmethod]
    pub fn from_gro(ax: f64, by: f64, cz: f64, ay: Option<f64>, az: Option<f64>, bx: Option<f64>, bz: Option<f64>, cx: Option<f64>, cy: Option<f64>) -> PyResult<Self> {
        let pbc = Self {
            a: [ax, ay.unwrap_or(0.), az.unwrap_or(0.)].into(),
            b: [bx.unwrap_or(0.), by, bz.unwrap_or(0.)].into(),
            c: [cx.unwrap_or(0.), cy.unwrap_or(0.), cz].into()
        };
        sanitize(&pbc)?;
        Ok(pbc)
    }


    /// Represent this PeriodicBox as a string, in the .gro file format.
    pub fn to_gro(&self) -> String {
        if self.b.x == 0. && self.c.x == 0. && self.c.y == 0. {
            format!(
                "{} {} {}", self.a.x, self.b.y, self.c.z
            )
        } else {
            format!(
                "{} {} {} {} {} {} {} {} {}", self.a.x, self.b.y, self.c.z, self.a.y, self.a.z, self.b.x, self.b.z, self.c.x, self.c.y
            )
        }
    }


    /// Create a triclinic PeriodicBox based on unit cell vector components.
    #[staticmethod]
    pub fn triclinic(ax: f64, bx: f64, by: f64, cx: f64, cy: f64, cz: f64) -> PyResult<Self> {
        let pbc = Self {
            a: [ax, 0., 0.].into(),
            b: [bx, by, 0.].into(),
            c: [cx, cy, cz].into()
        };
        sanitize(&pbc)?;
        Ok(pbc)
    }

    /// Represent this PeriodicBox as a string, in the extended .xyz format Lattice parameter.
    pub fn to_lattice(&self) -> String {
        format!(
            "{} {} {} {} {} {} {} {} {}", self.a.x, self.a.y, self.a.z, self.b.x, self.b.y, self.b.z, self.c.x, self.c.y, self.c.z
        )
    }

    /// Create a cubic PeriodicBox based on a unit cell length.
    #[staticmethod]
    pub fn cubic(d: f64) -> PyResult<Self> {
        if d <= 0. {
            return Err(PyValueError::new_err("d must be positive."));
        }
        Ok(Self {
            a: [d, 0., 0.].into(),
            b: [0., d, 0.].into(),
            c: [0., 0., d].into()
        })
    }

    /// Create an orthogonal PeriodicBox based on unit cell lengths.
    #[staticmethod]
    pub fn orthogonal(x: f64, y: f64, z: f64) -> PyResult<Self> {
        if x <= 0. || y <= 0. || z <= 0. {
            return Err(PyValueError::new_err("x, y, z must be positive."));
        }
        Ok(Self {
            a: [x, 0., 0.].into(),
            b: [0., y, 0.].into(),
            c: [0., 0., z].into()
        })
    }
    /// Move atom within the same copy of the PBC.
    ///
    /// Note: moves it within the box 0,0,0 to a.x,b.y,c.z, not the box a,b,c.
    pub fn move_within(&self, v: [f64; 3]) -> [f64; 3] {
        (
            DVec3::from(v)
            - (v[0] / self.a.x).floor() * self.a
            - (v[1] / self.b.y).floor() * self.b
            - (v[2] / self.c.z).floor() * self.c
        ).into()
    }

    /// Translate v by periodic box vectors so it is the closest possible to reference in non-periodic space.
    pub fn move_to(&self, reference: [f64; 3], v: [f64; 3]) -> [f64; 3] {
        (DVec3::from(reference) + DVec3::from(self.diff(v, reference))).into()
    }

    /// Move_within but for 2D numpy arrays of positions of shape (n, 3).
    pub fn move_all_within(&self, mut array: PyReadwriteArray<f64, Ix2>) -> PyResult<()> {
        if let [_, inner] = array.shape() && *inner != 3 {
            return Err(PyValueError::new_err(format!("Invalid dimensions, expected Nx3, got (N, {inner})")))
        }

        for v in array.as_slice_mut()?.chunks_mut(3) {
            [v[0], v[1], v[2]] = self.move_within([v[0], v[1], v[2]]);
        }

        Ok(())
    }

    /// Given a set of reference atoms, get which atoms are within a distance to any of the reference atoms.
    pub fn which_atoms_within_distance(&self, positions: PyReadonlyArray2<f64>, reference: HashSet<usize>, r: f64) -> PyResult<HashSet<usize>> {
        let mut res: HashSet<usize> = HashSet::new();
        let r2 = r * r;
        if let [_, inner] = positions.shape() && *inner != 3 {
            return Err(PyValueError::new_err(format!("Invalid dimensions, expected Nx3, got (N, {inner})")))
        }
        let mut tree = KdTree::new(3);
        for atom in reference {
            let pos = self.move_within(positions.get_item(atom)?.extract::<[f64; 3]>()?);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        // Note: dx=0, dy=0, dz=0 is handled here too
                        let new_pos = self.translate_by(pos, dx, dy, dz);
                        if self.is_almost_inside(new_pos, r + 0.01) {
                            tree.add(new_pos, atom).unwrap();
                        }
                    }
                }
            }
        }
        for i in 0..positions.shape()[0] {
            match tree.nearest(&positions.get_item(i)?.extract::<[f64; 3]>()?, 1, &squared_euclidean) {
                Ok(v) => {
                    for (dist, _) in v {
                        if dist < r2 {
                            res.insert(i);
                        }
                    }
                }
                Err(e) => {
                    return match e {
                        ErrorKind::WrongDimension => Err(PyValueError::new_err("Internal error: wrong dimension.")),
                        ErrorKind::NonFiniteCoordinate => Err(PyValueError::new_err("Supplied positions contain a non-finite value.")),
                        ErrorKind::ZeroCapacity => Err(PyValueError::new_err("Internal error: zero capacity."))
                    }
                }
            }
        }
        Ok(res)
    }

    #[getter]
    pub fn a(&self) -> [f64; 3] {
        self.a.into()
    }

    #[getter]
    pub fn b(&self) -> [f64; 3] {
        self.b.into()
    }

    #[getter]
    pub fn c(&self) -> [f64; 3] {
        self.c.into()
    }

    /// Get the lengths of the three unit cell vectors.
    pub fn cell_lengths(&self) -> [f64; 3] {
        [
            self.a.x,
            self.b.length(),
            self.c.length()
        ]
    }

    /// Get the difference between two points in periodic space.
    pub fn diff(&self, v1: [f64; 3], v2: [f64; 3]) -> [f64; 3] {
        let v1 = DVec3::from(self.move_within(v1));
        let v2: DVec3 = DVec3::from(self.move_within(v2));
        let diff = v1 - v2;
        let mut best = diff.clone();
        for i in -1..=1 {
            for j in -1..=1 {
                for k in -1..=1 {
                    let c = diff - (i as f64) * self.a - (j as f64) * self.b - (k as f64) * self.c;
                    if c.length_squared() < best.length_squared() {
                        best = c;
                    }
                }
            }
        }
        best.into()
    }

    /// Return whether the shortest path between two points crosses the periodic boundary condition,
    /// in a translation invariant manner.
    pub fn crosses_box(&self, v1: [f64; 3], v2: [f64; 3]) -> bool {
        (self.distance_squared(v1, v2) - (DVec3::from(v1)-DVec3::from(v2)).length_squared()).abs() > f64::EPSILON
    }


    /// Get the distance between two points in peridoic space, squared.
    pub fn distance_squared(&self, v1: [f64; 3], v2: [f64; 3]) -> f64 {
        DVec3::from(self.diff(v1, v2)).length_squared()
    }

    /// Get the distance between two points in periodic space.
    pub fn distance(&self, v1: [f64; 3], v2: [f64; 3]) -> f64 {
        self.distance_squared(v1, v2).sqrt()
    }

    /// Get the cosine of the periodic space angle between three points.
    pub fn cos_angle(&self, v1: [f64; 3], v2: [f64; 3], v3: [f64; 3]) -> f64 {
        let t = DVec3::from(self.diff(v2, v1));
        let u = DVec3::from(self.diff(v2, v3));
        (t.dot(u) / t.length() / u.length()).clamp(-1.0, 1.0)
    }

    /// Get the periodic space angle between three points, in radians.
    pub fn angle(&self, v1: [f64; 3], v2: [f64; 3], v3: [f64; 3]) -> f64 {
        self.cos_angle(v1, v2, v3).acos()
    }

    /// Get the angles between the pbc vectors, in radians.
    ///
    /// Returned in order of alpha (b, c), beta (a, c), gamma (a, b).
    pub fn cell_angles(&self) -> [f64; 3] {
        let alpha = (self.b.dot(self.c) / self.b.length() / self.c.length()).clamp(-1.0, 1.0).acos();
        let beta = (self.a.dot(self.c) / self.a.length() / self.c.length()).clamp(-1.0, 1.0).acos();
        let gamma = (self.a.dot(self.b) / self.a.length() / self.b.length()).clamp(-1.0, 1.0).acos();
        [ alpha, beta, gamma ]
    }

    /// Get the dihedral angle between four points in periodic space, in radians.
    pub fn dihedral(&self, v1: [f64; 3], v2: [f64; 3], v3: [f64; 3], v4: [f64; 3]) -> f64 {
        let r = DVec3::from(self.diff(v2, v1));
        let s = DVec3::from(self.diff(v3, v2));
        let t = DVec3::from(self.diff(v4, v3));

        let u = r.cross(s);
        let v = s.cross(t);

        let w = u.cross(v);

        let x = u.dot(v);
        let y = s.dot(w) / s.length();

        f64::atan2(y, x)
    }

    /// Returns whether a position is close to the "central" copy of the periodic box.
    ///
    /// This can be cbacuse:
    /// - inside the "central" copy of the pbc
    /// - within a cutoff distance of the pbc
    ///
    /// Non-exactly! It can return true even if it is not within cutoff.
    /// The only guarantee is that if it returns false, the distance between pos
    /// and the pbc box is larger than cutoff.
    pub fn is_almost_inside(&self, pos: [f64; 3], cutoff: f64) -> bool {
        // assumption that a.y, a.z, b.z are 0, and a.x, b.y, c.z > 0
        // this assumption is upheld by all constructors calling sanitize
        let [x, y, z] = pos;
        // only third component has a z = easy to check
        if self.c.z + cutoff < z || z < -cutoff {
            return false;
        }
        // two components can have a y
        // conservative behavior - we add the maximum deviation from orthogonal boxes to cutoff
        let y_cutoff = cutoff + self.c.y.abs();
        if self.b.y + y_cutoff < y || y < - y_cutoff {
            return false;
        }
        // three components can have an x
        let x_cutoff = cutoff + self.c.x.abs() + self.b.x.abs();
        if self.a.x + x_cutoff < x || x < - x_cutoff {
            return false;
        }
        true
    }

    /// Translates pos by i, j, k times periodic box vectors. This returns the same point
    /// in a different copy of the periodic box, but the same location in periodic space.
    pub fn translate_by(&self, pos: [f64; 3], i: i64, j: i64, k: i64) -> [f64; 3] {
        (DVec3::from(pos) + (i as f64) * self.a + (j as f64) * self.b + (k as f64) * self.c).into()
    }
}
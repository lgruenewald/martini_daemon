use pyo3::prelude::*;
use numpy::{PyReadwriteArray, Ix2, PyUntypedArrayMethods};

use glam::DVec3;
use pyo3::exceptions::PyValueError;

#[pyclass]
pub struct PeriodicBox {
    a: DVec3,
    b: DVec3,
    c: DVec3,
}

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
#[pymethods]
impl PeriodicBox {
    #[new]
    pub fn new(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> PyResult<Self> {
        let pbc = Self { a: a.into(), b: b.into(), c: c.into() };
        sanitize(&pbc)?;
        Ok(pbc)
    }

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

    pub fn to_lattice(&self) -> String {
        format!(
            "{} {} {} {} {} {} {} {} {}", self.a.x, self.a.y, self.a.z, self.b.x, self.b.y, self.b.z, self.c.x, self.c.y, self.c.z
        )
    }

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

    pub fn move_within(&self, v: [f64; 3]) -> [f64; 3] {
        (
            DVec3::from(v)
            - (v[0] / self.a.x).floor() * self.a
            - (v[1] / self.b.y).floor() * self.b
            - (v[2] / self.c.z).floor() * self.c
        ).into()
    }

    pub fn move_to(&self, reference: [f64; 3], v: [f64; 3]) -> [f64; 3] {
        (DVec3::from(reference) + DVec3::from(self.diff(v, reference))).into()
    }

    pub fn move_all_within(&self, mut array: PyReadwriteArray<f64, Ix2>) -> PyResult<()> {
        if let [_, inner] = array.shape() && *inner != 3 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!("Invalid dimensions, expected Nx3, got (N, {inner})")))
        }

        for v in array.as_slice_mut()?.chunks_mut(3) {
            [v[0], v[1], v[2]] = self.move_within([v[0], v[1], v[2]]);
        }

        Ok(())
    }

    #[getter]
    fn a(&self) -> [f64; 3] {
        self.a.into()
    }

    #[getter]
    fn b(&self) -> [f64; 3] {
        self.b.into()
    }

    #[getter]
    fn c(&self) -> [f64; 3] {
        self.c.into()
    }

    pub fn diff(&self, v1: [f64; 3], v2: [f64; 3]) -> [f64; 3] {
        let v1 = DVec3::from(self.move_within(v1));
        let v2: DVec3 = DVec3::from(self.move_within(v2));
        let mut diff = v1 - v2;
        for i in -1..=1 {
            for j in -1..=1 {
                for k in -1..=1 {
                    let c = diff - (i as f64) * self.a - (j as f64) * self.b - (k as f64) * self.c;
                    if c.length_squared() < diff.length_squared() {
                        diff = c;
                    }
                }
            }
        }
        diff.into()
    }
    pub fn crosses_box(&self, v1: [f64; 3], v2: [f64; 3]) -> bool {
        (self.distance_squared(v1, v2) - (DVec3::from(v1)-DVec3::from(v2)).length_squared()).abs() > f64::EPSILON
    }

    pub fn distance_squared(&self, v1: [f64; 3], v2: [f64; 3]) -> f64 {
        DVec3::from(self.diff(v1, v2)).length_squared()
    }

    pub fn distance(&self, v1: [f64; 3], v2: [f64; 3]) -> f64 {
        self.distance_squared(v1, v2).sqrt()
    }

    pub fn cos_angle(&self, v1: [f64; 3], v2: [f64; 3], v3: [f64; 3]) -> f64 {
        let t = DVec3::from(self.diff(v2, v1));
        let u = DVec3::from(self.diff(v2, v3));
        (t.dot(u) / t.length() / u.length()).clamp(-1.0, 1.0)
    }

    pub fn angle(&self, v1: [f64; 3], v2: [f64; 3], v3: [f64; 3]) -> f64 {
        self.cos_angle(v1, v2, v3).acos()
    }

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
}

#[cfg(test)]
fn is_close(a: f64, b: f64) -> bool {
    (a-b).abs() <= f64::EPSILON
}
#[test]
fn test_pbc() {
    use std::f64::consts::{PI, SQRT_2};

    let pbc = PeriodicBox::orthogonal(5., 5., 5.).unwrap();
    let pbc2 = PeriodicBox::orthogonal(1., 2., 3.).unwrap();

    let origin = [0., 0., 0.];
    let v1 = [0., 0., 1.];
    let v2 = [0., 0., 4.];
    let v3 = [1., 0., 0.];
    let v4 = [2.5, 0., 0.];
    let v5 = [1., 0., 1.];

    assert!(is_close(pbc.distance(v1, v2), 2.));
    assert!(is_close(pbc.distance(v1, v3), SQRT_2));
    assert!(is_close(pbc2.distance(v3, v4), 0.5));

    assert!(is_close(pbc.cos_angle(v1, origin, v3), 0.));
    assert!(is_close(pbc.cos_angle(v1, origin, v2), -1.));
    assert!(is_close(pbc.cos_angle(v1, origin, v5), (PI / 4.).cos()));

    /*
     * TETRAHEDRON
     * corners of a tetrahedron constructed from the box centered at origin,
     * edge length 1, but a random multiple of 5 added to coordinates
     * t1 to t4 are essentially
     * 1/2 1/2 1/2
     * -1/2 -1/2 1/2
     * -1/2 1/2 -1/2
     * 1/2 -1/2 -1/2
     */

    let t1 = [5.5, 0.5, 0.5];
    let t2 = [-10.5, 9.5, -4.5];
    let t3 = [4.5, 5.5, -0.5];
    let t4 = [0.5, -0.5, -0.5];
    assert!(is_close(pbc.cos_angle(t1, t2, t3), (PI / 3.).cos()));
    assert!(is_close(pbc.cos_angle(t2, t3, t4), (PI / 3.).cos()));
    assert!(is_close(pbc.cos_angle(t1, t4, t3), (PI / 3.).cos()));
    assert!(is_close(pbc.cos_angle(t1, t2, t4), (PI / 3.).cos()));

    // single line
    let l1 = [1., 0., 0.];
    let l2 = [2., 0., 0.];
    let l3 = [0., 0., 0.];
    let l4 = [4., 0., 0.];

    assert!(is_close(pbc.dihedral(l1, l2, l3, l4), 0.));
    assert!(is_close(pbc.dihedral(l1, l2, l4, l3), 0.));

    // dihedrals
    // 1 0 0
    let d1 = [1., 5., 0.];
    // 2 0 1
    let d2 = [7., 0., 6.];
    // 3 0 1
    let d3 = [8., 5., 1.];
    // angle 0
    let d4_0 = [4., 0., 0.];
    // 45
    let d4_45 = [4., 1., 0.];
    // 90
    let d4_90 = [4., 1., 1.];
    // 135
    let d4_135 = [4., 1., 2.];
    // 180
    let d4_180 = [4., 0., 2.];
    // 225
    let d4_225 = [4., -1., 2.];
    // 270
    let d4_270 = [4., -1., 1.];

    assert!(is_close(pbc.dihedral(d1, d2, d3, d4_0), 0.));
    assert!(is_close(pbc.dihedral(d1, d2, d3, d4_45), PI / 4.));
    assert!(is_close(pbc.dihedral(d1, d2, d3, d4_90), PI / 2.));
    assert!(is_close(pbc.dihedral(d1, d2, d3, d4_135), 3. * PI / 4.));
    assert!(is_close(pbc.dihedral(d1, d2, d3, d4_180), PI));
    assert!(is_close(pbc.dihedral(d1, d2, d3, d4_225), - 3. * PI / 4.));
    assert!(is_close(pbc.dihedral(d1, d2, d3, d4_270), - PI / 2.));

    // test iupac convention on sign
    let e1 = [4., 0., 1.];
    let e2 = [3., 0., 0.];
    let e3 = [2., 0., 0.];
    let e4 = [1., 1., 1.];

    assert!(is_close(pbc.dihedral(e1, e2, e3, e4), PI / 4.));
    assert!(is_close(pbc.dihedral(e4, e3, e2, e1), PI / 4.));


    assert!(pbc.crosses_box(v1, v2));
    assert!(!pbc.crosses_box(v1, v3));
    // TODO: test all functions, test on non orthogonal boxes
}

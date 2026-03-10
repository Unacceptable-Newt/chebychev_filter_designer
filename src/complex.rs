// ─────────────────────────────────────────────────────────────────────────────
// Complex number helper (avoids pulling in a crate)
// ─────────────────────────────────────────────────────────────────────────────
#[derive(Debug, Clone, Copy)]
pub struct C(f64, f64); // C(re, im)

impl C {
    pub fn re(x: f64) -> Self { C(x, 0.0) }
    pub fn im(x: f64) -> Self { C(0.0, x) }
    pub fn norm_sq(self) -> f64 { self.0 * self.0 + self.1 * self.1 }
    pub fn abs(self) -> f64 { self.norm_sq().sqrt() }
    pub fn inv(self) -> Self {
        let d = self.norm_sq();
        C(self.0 / d, -self.1 / d)
    }
}

impl std::ops::Add for C {
    type Output = C;
    fn add(self, r: C) -> C { C(self.0 + r.0, self.1 + r.1) }
}
impl std::ops::Mul for C {
    type Output = C;
    fn mul(self, r: C) -> C {
        C(self.0 * r.0 - self.1 * r.1, self.0 * r.1 + self.1 * r.0)
    }
}
impl std::ops::Div for C {
    type Output = C;
    fn div(self, r: C) -> C { self * r.inv() }
}

// 2×2 complex matrix  [[a, b], [c, d]]
pub type Mat2 = [[C; 2]; 2];

pub fn mat_mul(a: Mat2, b: Mat2) -> Mat2 {
    [
        [a[0][0] * b[0][0] + a[0][1] * b[1][0],
         a[0][0] * b[0][1] + a[0][1] * b[1][1]],
        [a[1][0] * b[0][0] + a[1][1] * b[1][0],
         a[1][0] * b[0][1] + a[1][1] * b[1][1]],
    ]
}

pub fn identity() -> Mat2 {
    let o = C::re(0.0); let i = C::re(1.0);
    [[i, o], [o, i]]
}

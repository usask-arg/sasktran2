use std::ops::{Add, Div, Mul, Neg, Sub};

use super::{Geometry, LayerAdjoints, LayerInputs, RadianceBatch, SourceMode, View};

const LANES: usize = 4;
const FOUR_PI: f64 = 4.0 * std::f64::consts::PI;

#[cfg(feature = "simd")]
#[derive(Clone, Copy, Debug)]
struct Wide(std::simd::f64x4);

#[cfg(not(feature = "simd"))]
#[derive(Clone, Copy, Debug)]
struct Wide([f64; LANES]);

impl Wide {
    #[cfg(feature = "simd")]
    fn splat(value: f64) -> Self {
        Self(std::simd::f64x4::splat(value))
    }

    #[cfg(not(feature = "simd"))]
    fn splat(value: f64) -> Self {
        Self([value; LANES])
    }

    fn load(values: &[f64], base: usize, stride: usize) -> Self {
        let mut lanes = [0.0; LANES];
        for (lane, target) in lanes.iter_mut().enumerate() {
            let wave = base + lane;
            if wave < stride {
                *target = values[wave];
            }
        }
        Self::from_array(lanes)
    }

    fn load_row(values: &[f64], row: usize, base: usize, stride: usize) -> Self {
        Self::load(&values[row * stride..(row + 1) * stride], base, stride)
    }

    fn cotangent(values: &[f64], view: usize, base: usize, stride: usize) -> Self {
        Self::load_row(values, view, base, stride)
    }

    #[cfg(feature = "simd")]
    fn from_array(values: [f64; LANES]) -> Self {
        Self(std::simd::f64x4::from_array(values))
    }

    #[cfg(not(feature = "simd"))]
    fn from_array(values: [f64; LANES]) -> Self {
        Self(values)
    }

    #[cfg(feature = "simd")]
    fn to_array(self) -> [f64; LANES] {
        self.0.to_array()
    }

    #[cfg(not(feature = "simd"))]
    fn to_array(self) -> [f64; LANES] {
        self.0
    }

    #[cfg(feature = "simd")]
    fn exp(self) -> Self {
        use std::simd::StdFloat;
        Self(self.0.exp())
    }

    #[cfg(not(feature = "simd"))]
    fn exp(self) -> Self {
        Self(self.0.map(f64::exp))
    }

    #[cfg(feature = "simd")]
    fn sqrt(self) -> Self {
        use std::simd::StdFloat;
        Self(self.0.sqrt())
    }

    #[cfg(not(feature = "simd"))]
    fn sqrt(self) -> Self {
        Self(self.0.map(f64::sqrt))
    }
}

macro_rules! wide_binary {
    ($trait:ident, $method:ident, $op:tt) => {
        impl $trait for Wide {
            type Output = Self;
            fn $method(self, rhs: Self) -> Self {
                #[cfg(feature = "simd")]
                { Self(self.0 $op rhs.0) }
                #[cfg(not(feature = "simd"))]
                { Self(std::array::from_fn(|i| self.0[i] $op rhs.0[i])) }
            }
        }
    };
}

wide_binary!(Add, add, +);
wide_binary!(Sub, sub, -);
wide_binary!(Mul, mul, *);
wide_binary!(Div, div, /);

impl Neg for Wide {
    type Output = Self;
    fn neg(self) -> Self {
        #[cfg(feature = "simd")]
        {
            Self(-self.0)
        }
        #[cfg(not(feature = "simd"))]
        {
            Self(self.0.map(|value| -value))
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum Op {
    Leaf,
    Affine { input: usize, derivative: f64 },
    Add { left: usize, right: usize },
    Sub { left: usize, right: usize },
    Mul { left: usize, right: usize },
    Div { left: usize, right: usize },
    Exp { input: usize },
    Sqrt { input: usize },
    ConstDiv { denominator: usize, numerator: f64 },
}

#[derive(Clone, Debug, Default)]
struct Tape {
    values: Vec<Wide>,
    adjoints: Vec<Wide>,
    ops: Vec<Op>,
}

impl Tape {
    fn reset(&mut self, capacity: usize) {
        self.values.clear();
        self.adjoints.clear();
        self.ops.clear();
        self.values.reserve(capacity);
        self.adjoints.reserve(capacity);
        self.ops.reserve(capacity);
    }

    fn push(&mut self, value: Wide, op: Op) -> Var {
        let index = self.values.len();
        self.values.push(value);
        self.adjoints.push(Wide::splat(0.0));
        self.ops.push(op);
        Var { index, tape: self }
    }

    fn leaf(&mut self, value: Wide) -> Var {
        self.push(value, Op::Leaf)
    }

    fn constant(&mut self, value: f64) -> Var {
        self.leaf(Wide::splat(value))
    }

    fn reverse(&mut self, outputs: &[(usize, Wide)]) {
        for (index, seed) in outputs {
            self.adjoints[*index] = self.adjoints[*index] + *seed;
        }
        for index in (0..self.values.len()).rev() {
            let adjoint = self.adjoints[index];
            match self.ops[index] {
                Op::Leaf => {}
                Op::Affine { input, derivative } => {
                    self.adjoints[input] = self.adjoints[input] + adjoint * Wide::splat(derivative);
                }
                Op::Add { left, right } => {
                    self.adjoints[left] = self.adjoints[left] + adjoint;
                    self.adjoints[right] = self.adjoints[right] + adjoint;
                }
                Op::Sub { left, right } => {
                    self.adjoints[left] = self.adjoints[left] + adjoint;
                    self.adjoints[right] = self.adjoints[right] - adjoint;
                }
                Op::Mul { left, right } => {
                    self.adjoints[left] = self.adjoints[left] + adjoint * self.values[right];
                    self.adjoints[right] = self.adjoints[right] + adjoint * self.values[left];
                }
                Op::Div { left, right } => {
                    let denominator = self.values[right];
                    self.adjoints[left] = self.adjoints[left] + adjoint / denominator;
                    self.adjoints[right] = self.adjoints[right]
                        - adjoint * self.values[left] / (denominator * denominator);
                }
                Op::Exp { input } => {
                    self.adjoints[input] = self.adjoints[input] + adjoint * self.values[index];
                }
                Op::Sqrt { input } => {
                    self.adjoints[input] =
                        self.adjoints[input] + adjoint * Wide::splat(0.5) / self.values[index];
                }
                Op::ConstDiv {
                    denominator,
                    numerator,
                } => {
                    let value = self.values[denominator];
                    self.adjoints[denominator] = self.adjoints[denominator]
                        - adjoint * Wide::splat(numerator) / (value * value);
                }
            }
        }
    }

    fn adjoint(&self, index: usize) -> Wide {
        self.adjoints[index]
    }

    fn capacity_bytes(&self) -> usize {
        self.values.capacity() * std::mem::size_of::<Wide>()
            + self.adjoints.capacity() * std::mem::size_of::<Wide>()
            + self.ops.capacity() * std::mem::size_of::<Op>()
    }
}

#[derive(Clone, Debug, Default)]
pub(super) struct ReverseWorkspace {
    tape: Tape,
    optical_depth: Vec<Var>,
    single_scatter_albedo: Vec<Var>,
    first_legendre: Vec<Var>,
    transmission: Vec<Var>,
    average_secant: Vec<Var>,
    thermal_b0: Vec<Var>,
    thermal_b1: Vec<Var>,
    homogeneous: [Homogeneous; 2],
    particular: [Particular; 2],
    bvp: [Bvp; 2],
    source: Vec<Var>,
    outputs: Vec<Var>,
    seeds: Vec<(usize, Wide)>,
}

impl ReverseWorkspace {
    pub(super) fn capacity_bytes(&self) -> usize {
        let mut vars = self.optical_depth.capacity()
            + self.single_scatter_albedo.capacity()
            + self.first_legendre.capacity()
            + self.transmission.capacity()
            + self.average_secant.capacity()
            + self.thermal_b0.capacity()
            + self.thermal_b1.capacity()
            + self.source.capacity()
            + self.outputs.capacity();
        for h in &self.homogeneous {
            vars += h.k.capacity() + h.xp.capacity() + h.xm.capacity() + h.omega.capacity();
        }
        for p in &self.particular {
            vars += p.ap.capacity()
                + p.am.capacity()
                + p.at.capacity()
                + p.gpt.capacity()
                + p.gpb.capacity()
                + p.gmt.capacity()
                + p.gmb.capacity();
        }
        for q in &self.bvp {
            vars += q.e.capacity()
                + q.c.capacity()
                + q.d.capacity()
                + q.a.capacity()
                + q.b.capacity()
                + q.rhs.capacity()
                + q.alpha.capacity()
                + q.beta.capacity()
                + q.gamma.capacity()
                + q.mu.capacity()
                + q.z.capacity();
        }
        self.tape.capacity_bytes()
            + vars * std::mem::size_of::<Var>()
            + self.seeds.capacity() * std::mem::size_of::<(usize, Wide)>()
    }
}

#[derive(Clone, Copy, Debug)]
struct Var {
    index: usize,
    tape: *mut Tape,
}

impl Var {
    fn value(self) -> Wide {
        // SAFETY: Vars never escape `differentiate_chunk`; its Tape remains
        // alive and is not moved while the recorded computation is evaluated.
        unsafe { (&(*self.tape).values)[self.index] }
    }

    fn affine(self, value: Wide, derivative: f64) -> Self {
        unsafe {
            (*self.tape).push(
                value,
                Op::Affine {
                    input: self.index,
                    derivative,
                },
            )
        }
    }

    fn binary(self, rhs: Self, value: Wide, op: Op) -> Self {
        debug_assert_eq!(self.tape, rhs.tape);
        unsafe { (*self.tape).push(value, op) }
    }

    fn exp(self) -> Self {
        let value = self.value().exp();
        unsafe { (*self.tape).push(value, Op::Exp { input: self.index }) }
    }

    fn sqrt(self) -> Self {
        let value = self.value().sqrt();
        unsafe { (*self.tape).push(value, Op::Sqrt { input: self.index }) }
    }
}

impl Add for Var {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        self.binary(
            rhs,
            self.value() + rhs.value(),
            Op::Add {
                left: self.index,
                right: rhs.index,
            },
        )
    }
}
impl Sub for Var {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self.binary(
            rhs,
            self.value() - rhs.value(),
            Op::Sub {
                left: self.index,
                right: rhs.index,
            },
        )
    }
}
impl Mul for Var {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        self.binary(
            rhs,
            self.value() * rhs.value(),
            Op::Mul {
                left: self.index,
                right: rhs.index,
            },
        )
    }
}
impl Div for Var {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        self.binary(
            rhs,
            self.value() / rhs.value(),
            Op::Div {
                left: self.index,
                right: rhs.index,
            },
        )
    }
}
impl Neg for Var {
    type Output = Self;
    fn neg(self) -> Self {
        self.affine(-self.value(), -1.0)
    }
}

macro_rules! var_scalar {
    ($trait:ident, $method:ident, $value:expr, $derivative:expr) => {
        impl $trait<f64> for Var {
            type Output = Self;
            fn $method(self, rhs: f64) -> Self {
                self.affine($value(self.value(), rhs), $derivative(rhs))
            }
        }
    };
}
var_scalar!(Add, add, |a: Wide, b| a + Wide::splat(b), |_| 1.0);
var_scalar!(Sub, sub, |a: Wide, b| a - Wide::splat(b), |_| 1.0);
var_scalar!(Mul, mul, |a: Wide, b| a * Wide::splat(b), |b| b);
var_scalar!(Div, div, |a: Wide, b| a / Wide::splat(b), |b| 1.0 / b);

impl Add<Var> for f64 {
    type Output = Var;
    fn add(self, rhs: Var) -> Var {
        rhs + self
    }
}
impl Mul<Var> for f64 {
    type Output = Var;
    fn mul(self, rhs: Var) -> Var {
        rhs * self
    }
}
impl Sub<Var> for f64 {
    type Output = Var;
    fn sub(self, rhs: Var) -> Var {
        (-rhs) + self
    }
}
impl Div<Var> for f64 {
    type Output = Var;
    fn div(self, rhs: Var) -> Var {
        let value = Wide::splat(self) / rhs.value();
        unsafe {
            (*rhs.tape).push(
                value,
                Op::ConstDiv {
                    denominator: rhs.index,
                    numerator: self,
                },
            )
        }
    }
}

#[derive(Clone, Debug, Default)]
struct Homogeneous {
    k: Vec<Var>,
    xp: Vec<Var>,
    xm: Vec<Var>,
    omega: Vec<Var>,
}
#[derive(Clone, Debug, Default)]
struct Particular {
    ap: Vec<Var>,
    am: Vec<Var>,
    at: Vec<Var>,
    gpt: Vec<Var>,
    gpb: Vec<Var>,
    gmt: Vec<Var>,
    gmb: Vec<Var>,
}
#[derive(Clone, Debug, Default)]
struct Bvp {
    e: Vec<Var>,
    c: Vec<Var>,
    d: Vec<Var>,
    a: Vec<Var>,
    b: Vec<Var>,
    rhs: Vec<Var>,
    alpha: Vec<Var>,
    beta: Vec<Var>,
    gamma: Vec<Var>,
    mu: Vec<Var>,
    z: Vec<Var>,
}

fn resize_vars(values: &mut Vec<Var>, len: usize, zero: Var) {
    values.resize(len, zero);
    values.fill(zero);
}

pub(super) fn solve_vjp(
    geometry: &Geometry,
    mode: SourceMode,
    inputs: &LayerInputs,
    views: &[View],
    cotangent: &[f64],
    workspace: &mut ReverseWorkspace,
) -> (RadianceBatch, LayerAdjoints) {
    let n = inputs.num_layers;
    let nw = inputs.num_wavelengths;
    let mut result = LayerAdjoints {
        optical_depth: vec![0.0; n * nw],
        single_scatter_albedo: vec![0.0; n * nw],
        first_legendre: vec![0.0; n * nw],
        transmission: (mode == SourceMode::Solar).then(|| vec![0.0; (n + 1) * nw]),
        average_secant: (mode == SourceMode::Solar).then(|| vec![0.0; n * nw]),
        thermal_b0: (mode == SourceMode::Thermal).then(|| vec![0.0; n * nw]),
        thermal_b1: (mode == SourceMode::Thermal).then(|| vec![0.0; n * nw]),
        surface_albedo: vec![0.0; nw],
        surface_emission: (mode == SourceMode::Thermal).then(|| vec![0.0; nw]),
    };
    let mut radiance = RadianceBatch {
        num_views: views.len(),
        num_wavelengths: nw,
        values: vec![0.0; views.len() * nw],
    };
    for base in (0..nw).step_by(LANES) {
        differentiate_chunk(
            geometry,
            mode,
            inputs,
            views,
            cotangent,
            base,
            &mut radiance,
            &mut result,
            workspace,
        );
    }
    (radiance, result)
}

fn write_adjoint(target: &mut [f64], row: usize, base: usize, stride: usize, value: Wide) {
    for (lane, value) in value.to_array().into_iter().enumerate() {
        if base + lane < stride {
            target[row * stride + base + lane] = value;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn differentiate_chunk(
    geometry: &Geometry,
    mode: SourceMode,
    inputs: &LayerInputs,
    views: &[View],
    cotangent: &[f64],
    base: usize,
    radiance: &mut RadianceBatch,
    result: &mut LayerAdjoints,
    workspace: &mut ReverseWorkspace,
) {
    let n = inputs.num_layers;
    let nw = inputs.num_wavelengths;
    let ReverseWorkspace {
        tape,
        optical_depth: od,
        single_scatter_albedo: ssa,
        first_legendre: b1,
        transmission,
        average_secant: secant,
        thermal_b0: tb0,
        thermal_b1: tb1,
        homogeneous: hs,
        particular: ps,
        bvp: bvps,
        source,
        outputs,
        seeds,
    } = workspace;
    tape.reset(n * 220 + views.len() * n * 50);
    let zero = tape.constant(0.0);
    resize_vars(od, n, zero);
    resize_vars(ssa, n, zero);
    resize_vars(b1, n, zero);
    for layer in 0..n {
        od[layer] = tape.leaf(Wide::load_row(&inputs.optical_depth, layer, base, nw));
        ssa[layer] = tape.leaf(Wide::load_row(
            &inputs.single_scatter_albedo,
            layer,
            base,
            nw,
        ));
        b1[layer] = tape.leaf(Wide::load_row(&inputs.first_legendre, layer, base, nw));
    }
    let albedo = tape.leaf(Wide::load(&inputs.surface_albedo, base, nw));
    resize_vars(transmission, n + 1, zero);
    resize_vars(secant, n, zero);
    resize_vars(tb0, n, zero);
    resize_vars(tb1, n, zero);
    let thermal_surface;
    if mode == SourceMode::Solar {
        for (boundary, target) in transmission.iter_mut().enumerate() {
            let v = tape.leaf(Wide::load_row(
                inputs.transmission.as_ref().unwrap(),
                boundary,
                base,
                nw,
            ));
            *target = v;
        }
        for (layer, target) in secant.iter_mut().enumerate() {
            let v = tape.leaf(Wide::load_row(
                inputs.average_secant.as_ref().unwrap(),
                layer,
                base,
                nw,
            ));
            *target = v;
        }
        thermal_surface = zero;
    } else {
        for layer in 0..n {
            let v = tape.leaf(Wide::load_row(
                inputs.thermal_b0.as_ref().unwrap(),
                layer,
                base,
                nw,
            ));
            tb0[layer] = v;
            let v = tape.leaf(Wide::load_row(
                inputs.thermal_b1.as_ref().unwrap(),
                layer,
                base,
                nw,
            ));
            tb1[layer] = v;
        }
        thermal_surface = tape.leaf(Wide::load(
            inputs.surface_emission.as_ref().unwrap(),
            base,
            nw,
        ));
    }

    let naz = if mode == SourceMode::Solar { 2 } else { 1 };
    for az in 0..naz {
        let h = &mut hs[az];
        resize_vars(&mut h.k, n, zero);
        resize_vars(&mut h.xp, n, zero);
        resize_vars(&mut h.xm, n, zero);
        resize_vars(&mut h.omega, n, zero);
        let p = &mut ps[az];
        resize_vars(&mut p.ap, n, zero);
        resize_vars(&mut p.am, n, zero);
        resize_vars(&mut p.at, n, zero);
        resize_vars(&mut p.gpt, n, zero);
        resize_vars(&mut p.gpb, n, zero);
        resize_vars(&mut p.gmt, n, zero);
        resize_vars(&mut p.gmb, n, zero);
        for layer in 0..n {
            let mu = geometry.quadrature_cosine;
            let (d, s) = if az == 0 {
                (
                    ssa[layer] * b1[layer] * mu - 1.0 / mu,
                    (ssa[layer] - 1.0) / mu,
                )
            } else {
                (
                    -ssa[layer] * 0.0 - 1.0 / mu,
                    (ssa[layer] * b1[layer] * (1.0 - mu * mu) - 2.0) / (2.0 * mu),
                )
            };
            let k = (s * d).sqrt();
            let xp = 0.5 * (1.0 - s / k);
            let xm = 0.5 * (1.0 + s / k);
            let omega = (-k * od[layer]).exp();
            h.k[layer] = k;
            h.xp[layer] = xp;
            h.xm[layer] = xm;
            h.omega[layer] = omega;
            let norm = mu * (xp * xp - xm * xm);
            if mode == SourceMode::Solar {
                let (qp, qm) = if az == 0 {
                    (
                        ssa[layer] * (1.0 + b1[layer] * geometry.solar_cosine * mu) / FOUR_PI,
                        ssa[layer] * (1.0 - b1[layer] * geometry.solar_cosine * mu) / FOUR_PI,
                    )
                } else {
                    let q = ssa[layer]
                        * b1[layer]
                        * ((1.0 - mu * mu) * (1.0 - geometry.solar_cosine.powi(2))).sqrt()
                        / FOUR_PI;
                    (q, q)
                };
                let ap = (qp * xp + qm * xm) / norm;
                let am = (qm * xp + qp * xm) / norm;
                let ex = (-secant[layer] * od[layer]).exp();
                let cp = transmission[layer] * (omega - ex) / (secant[layer] - k);
                let cm = transmission[layer] * (1.0 - omega * ex) / (secant[layer] + k);
                p.ap[layer] = ap;
                p.am[layer] = am;
                p.gpt[layer] = am * cm * xm;
                p.gpb[layer] = ap * cp * xp;
                p.gmt[layer] = am * cm * xp;
                p.gmb[layer] = ap * cp * xm;
            } else {
                let at = (1.0 - ssa[layer]) * (xp + xm) / norm;
                let ex = (-tb1[layer] * od[layer]).exp();
                let cp = tb0[layer] * (omega - ex) / (tb1[layer] - k);
                let cm = tb0[layer] * (1.0 - omega * ex) / (tb1[layer] + k);
                p.at[layer] = at;
                p.gpt[layer] = at * cm * xm;
                p.gpb[layer] = at * cp * xp;
                p.gmt[layer] = at * cm * xp;
                p.gmb[layer] = at * cp * xm;
            }
        }
    }

    for az in 0..naz {
        let h = &hs[az];
        let p = &ps[az];
        let size = 2 * n;
        let q = &mut bvps[az];
        resize_vars(&mut q.e, size, zero);
        resize_vars(&mut q.c, size, zero);
        resize_vars(&mut q.d, size, zero);
        resize_vars(&mut q.a, size, zero);
        resize_vars(&mut q.b, size, zero);
        resize_vars(&mut q.rhs, size, zero);
        resize_vars(&mut q.alpha, size, zero);
        resize_vars(&mut q.beta, size, zero);
        resize_vars(&mut q.gamma, size, zero);
        resize_vars(&mut q.mu, size, zero);
        resize_vars(&mut q.z, size, zero);
        q.rhs[0] = -p.gpt[0];
        for l in 0..n - 1 {
            q.rhs[2 * l + 1] = p.gmt[l + 1] - p.gmb[l];
            q.rhs[2 * l + 2] = p.gpt[l + 1] - p.gpb[l];
        }
        let last = size - 1;
        let delta = if az == 0 { 1.0 } else { 0.0 };
        let direct_boundary_source = if mode == SourceMode::Solar {
            delta * geometry.solar_cosine * albedo / std::f64::consts::PI * transmission[n]
        } else {
            thermal_surface
        };
        q.rhs[last] = direct_boundary_source
            - (p.gmb[n - 1] - 2.0 * delta * geometry.quadrature_cosine * albedo * p.gpb[n - 1]);
        q.d[0] = h.xp[0];
        q.a[0] = h.xm[0] * h.omega[0];
        for l in 0..n - 1 {
            let r = 2 * l;
            q.c[r + 1] = h.xm[l] * h.omega[l];
            q.d[r + 1] = h.xp[l];
            q.a[r + 1] = -h.xm[l + 1];
            q.b[r + 1] = -h.xp[l + 1] * h.omega[l + 1];
            q.e[r + 2] = h.xp[l] * h.omega[l];
            q.c[r + 2] = h.xm[l];
            q.d[r + 2] = -h.xp[l + 1];
            q.a[r + 2] = -h.xm[l + 1] * h.omega[l + 1];
        }
        q.c[last] = (h.xm[n - 1] - 2.0 * geometry.quadrature_cosine * albedo * delta * h.xp[n - 1])
            * h.omega[n - 1];
        q.d[last] = h.xp[n - 1] - 2.0 * geometry.quadrature_cosine * albedo * delta * h.xm[n - 1];
        pentadiagonal(q);
    }

    outputs.clear();
    outputs.reserve(views.len());
    for view in views {
        resize_vars(source, n, zero);
        let mut attenuation = tape.constant(1.0);
        let mut integrated = zero;
        for l in 0..n {
            let beam = (-od[l] / view.cosine).exp();
            let exs = if mode == SourceMode::Solar {
                (-secant[l] * od[l]).exp()
            } else {
                zero
            };
            let ext = if mode == SourceMode::Thermal {
                (-tb1[l] * od[l]).exp()
            } else {
                zero
            };
            let es = if mode == SourceMode::Solar {
                (1.0 - exs * beam) / (1.0 + secant[l] * view.cosine)
            } else {
                zero
            };
            let et = if mode == SourceMode::Thermal {
                (1.0 - ext * beam) / (1.0 + tb1[l] * view.cosine)
            } else {
                zero
            };
            for az in 0..naz {
                let h = &hs[az];
                let p = &ps[az];
                let azi = (az as f64 * view.relative_azimuth).cos();
                let (lp, lm) = if az == 0 {
                    (
                        ssa[l] * (1.0 - b1[l] * view.cosine * geometry.quadrature_cosine) * 0.5,
                        ssa[l] * (1.0 + b1[l] * view.cosine * geometry.quadrature_cosine) * 0.5,
                    )
                } else if mode == SourceMode::Solar {
                    let x = ssa[l]
                        * b1[l]
                        * ((1.0 - view.cosine.powi(2))
                            * (1.0 - geometry.quadrature_cosine.powi(2)))
                        .sqrt()
                        * 0.25;
                    (x, x)
                } else {
                    (zero, zero)
                };
                let yp = lp * h.xp[l] + lm * h.xm[l];
                let ym = lp * h.xm[l] + lm * h.xp[l];
                let hm = (h.omega[l] - beam) / (1.0 - h.k[l] * view.cosine);
                let hp = (1.0 - h.omega[l] * beam) / (1.0 + h.k[l] * view.cosine);
                let v = if mode == SourceMode::Solar {
                    let dp = transmission[l] * (es - exs * hm) / (secant[l] + h.k[l]);
                    let dm = transmission[l] * (hp - es) / (secant[l] - h.k[l]);
                    p.ap[l] * yp * dm + p.am[l] * ym * dp
                } else {
                    let dp = tb0[l] * (et - ext * hm) / (tb1[l] + h.k[l]);
                    let dm = tb0[l] * (hp - et) / (tb1[l] - h.k[l]);
                    p.at[l] * (yp * dm + ym * dp)
                };
                source[l] = source[l]
                    + azi * (bvps[az].rhs[2 * l] * yp * hp + bvps[az].rhs[2 * l + 1] * ym * hm + v);
                if mode == SourceMode::Thermal {
                    source[l] = source[l] + tb0[l] * et * (1.0 - ssa[l]);
                }
            }
            integrated = integrated + source[l] * attenuation;
            attenuation = attenuation * beam;
        }
        let l = n - 1;
        let surface = (ps[0].gpb[l]
            + bvps[0].rhs[2 * l] * hs[0].xp[l] * hs[0].omega[l]
            + bvps[0].rhs[2 * l + 1] * hs[0].xm[l])
            * 2.0
            * geometry.quadrature_cosine
            * albedo;
        outputs.push(integrated + attenuation * (surface + thermal_surface));
    }
    for (view, output) in outputs.iter().enumerate() {
        write_adjoint(&mut radiance.values, view, base, nw, output.value());
    }
    seeds.clear();
    seeds.reserve(views.len());
    seeds.extend(
        outputs
            .iter()
            .enumerate()
            .map(|(view, out)| (out.index, Wide::cotangent(cotangent, view, base, nw))),
    );
    tape.reverse(seeds);
    for l in 0..n {
        write_adjoint(
            &mut result.optical_depth,
            l,
            base,
            nw,
            tape.adjoint(od[l].index),
        );
        write_adjoint(
            &mut result.single_scatter_albedo,
            l,
            base,
            nw,
            tape.adjoint(ssa[l].index),
        );
        write_adjoint(
            &mut result.first_legendre,
            l,
            base,
            nw,
            tape.adjoint(b1[l].index),
        );
    }
    write_adjoint(
        &mut result.surface_albedo,
        0,
        base,
        nw,
        tape.adjoint(albedo.index),
    );
    if mode == SourceMode::Solar {
        for (b, value) in transmission.iter().enumerate() {
            write_adjoint(
                result.transmission.as_mut().unwrap(),
                b,
                base,
                nw,
                tape.adjoint(value.index),
            );
        }
        for (l, value) in secant.iter().enumerate() {
            write_adjoint(
                result.average_secant.as_mut().unwrap(),
                l,
                base,
                nw,
                tape.adjoint(value.index),
            );
        }
    } else {
        for l in 0..n {
            write_adjoint(
                result.thermal_b0.as_mut().unwrap(),
                l,
                base,
                nw,
                tape.adjoint(tb0[l].index),
            );
            write_adjoint(
                result.thermal_b1.as_mut().unwrap(),
                l,
                base,
                nw,
                tape.adjoint(tb1[l].index),
            );
        }
        write_adjoint(
            result.surface_emission.as_mut().unwrap(),
            0,
            base,
            nw,
            tape.adjoint(thermal_surface.index),
        );
    }
}

fn pentadiagonal(q: &mut Bvp) {
    let n = q.d.len();
    q.mu[0] = q.d[0];
    q.alpha[0] = q.a[0] / q.mu[0];
    q.beta[0] = q.b[0] / q.mu[0];
    q.z[0] = q.rhs[0] / q.mu[0];
    if n > 1 {
        q.gamma[1] = q.c[1];
        q.mu[1] = q.d[1] - q.alpha[0] * q.gamma[1];
        q.alpha[1] = (q.a[1] - q.beta[0] * q.gamma[1]) / q.mu[1];
        q.beta[1] = q.b[1] / q.mu[1];
        q.z[1] = (q.rhs[1] - q.z[0] * q.gamma[1]) / q.mu[1];
    }
    for i in 2..n {
        q.gamma[i] = q.c[i] - q.alpha[i - 2] * q.e[i];
        q.mu[i] = q.d[i] - q.beta[i - 2] * q.e[i] - q.alpha[i - 1] * q.gamma[i];
        if i + 1 < n {
            q.alpha[i] = (q.a[i] - q.beta[i - 1] * q.gamma[i]) / q.mu[i];
        }
        if i + 2 < n {
            q.beta[i] = q.b[i] / q.mu[i];
        }
        q.z[i] = (q.rhs[i] - q.z[i - 2] * q.e[i] - q.z[i - 1] * q.gamma[i]) / q.mu[i];
    }
    q.rhs[n - 1] = q.z[n - 1];
    if n > 1 {
        q.rhs[n - 2] = q.z[n - 2] - q.alpha[n - 2] * q.rhs[n - 1];
    }
    for i in (0..n.saturating_sub(2)).rev() {
        q.rhs[i] = q.z[i] - q.alpha[i] * q.rhs[i + 1] - q.beta[i] * q.rhs[i + 2];
    }
}

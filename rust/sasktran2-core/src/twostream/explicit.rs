//! Hand-specialized reverse mode for the wavelength-batched two-stream solve.

use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

use super::{
    AtmosphereAdjoints, AtmosphereBatch, AtmosphereJacobians, ExecutionPolicy, Geometry,
    LayerAdjoints, LayerInputs, RadianceBatch, SourceMode, SphericalGeometry, View,
};

pub(super) const LANES: usize = 8;
const FOUR_PI: f64 = 4.0 * std::f64::consts::PI;

#[derive(Clone, Copy, Debug)]
struct Wide([f64; LANES]);

impl Wide {
    fn splat(value: f64) -> Self {
        Self([value; LANES])
    }

    fn load(values: &[f64], base: usize, stride: usize) -> Self {
        if base + LANES <= stride {
            return Self::from_array(values[base..base + LANES].try_into().unwrap());
        }
        let mut lanes = [0.0; LANES];
        for (lane, target) in lanes.iter_mut().enumerate() {
            if base + lane < stride {
                *target = values[base + lane];
            }
        }
        Self::from_array(lanes)
    }

    fn load_row(values: &[f64], row: usize, base: usize, stride: usize) -> Self {
        Self::load(&values[row * stride..(row + 1) * stride], base, stride)
    }

    fn from_array(values: [f64; LANES]) -> Self {
        Self(values)
    }

    fn to_array(self) -> [f64; LANES] {
        self.0
    }

    fn exp(self) -> Self {
        Self(self.0.map(f64::exp))
    }

    fn ln(self) -> Self {
        Self(self.0.map(f64::ln))
    }

    fn min(self, rhs: Self) -> Self {
        Self(std::array::from_fn(|i| self.0[i].min(rhs.0[i])))
    }

    fn sqrt(self) -> Self {
        Self(self.0.map(f64::sqrt))
    }

    fn transmission_ratio(top: Self, bottom: Self, exponent: Self) -> Self {
        Self(std::array::from_fn(|lane| {
            if top.0[lane].is_normal() {
                let ratio = bottom.0[lane] / top.0[lane];
                if ratio.is_finite() {
                    return ratio;
                }
            }
            exponent.0[lane].exp()
        }))
    }
}

macro_rules! wide_binary {
    ($trait:ident, $method:ident, $op:tt) => {
        impl $trait for Wide {
            type Output = Self;
            fn $method(self, rhs: Self) -> Self {
                Self(std::array::from_fn(|i| self.0[i] $op rhs.0[i]))
            }
        }
    };
}

wide_binary!(Add, add, +);
wide_binary!(Sub, sub, -);
wide_binary!(Mul, mul, *);
wide_binary!(Div, div, /);

impl AddAssign for Wide {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl SubAssign for Wide {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Neg for Wide {
    type Output = Self;
    fn neg(self) -> Self {
        Self(self.0.map(|value| -value))
    }
}

macro_rules! wide_scalar {
    ($trait:ident, $method:ident, $op:tt) => {
        impl $trait<f64> for Wide {
            type Output = Self;
            fn $method(self, rhs: f64) -> Self {
                self $op Self::splat(rhs)
            }
        }
    };
}

wide_scalar!(Add, add, +);
wide_scalar!(Sub, sub, -);
wide_scalar!(Mul, mul, *);
wide_scalar!(Div, div, /);

impl Add<Wide> for f64 {
    type Output = Wide;
    fn add(self, rhs: Wide) -> Wide {
        Wide::splat(self) + rhs
    }
}

impl Sub<Wide> for f64 {
    type Output = Wide;
    fn sub(self, rhs: Wide) -> Wide {
        Wide::splat(self) - rhs
    }
}

impl Mul<Wide> for f64 {
    type Output = Wide;
    fn mul(self, rhs: Wide) -> Wide {
        rhs * self
    }
}

impl Div<Wide> for f64 {
    type Output = Wide;
    fn div(self, rhs: Wide) -> Wide {
        Wide::splat(self) / rhs
    }
}

#[derive(Clone, Debug, Default)]
struct Homogeneous {
    d: Vec<Wide>,
    s: Vec<Wide>,
    k: Vec<Wide>,
    xp: Vec<Wide>,
    xm: Vec<Wide>,
    omega: Vec<Wide>,
    norm: Vec<Wide>,
}

#[derive(Clone, Debug, Default)]
struct HomogeneousAdjoint {
    k: Vec<Wide>,
    xp: Vec<Wide>,
    xm: Vec<Wide>,
    omega: Vec<Wide>,
}

#[derive(Clone, Debug, Default)]
struct Particular {
    qp: Vec<Wide>,
    qm: Vec<Wide>,
    ap: Vec<Wide>,
    am: Vec<Wide>,
    at: Vec<Wide>,
    exponential: Vec<Wide>,
    cp: Vec<Wide>,
    cm: Vec<Wide>,
    gpt: Vec<Wide>,
    gpb: Vec<Wide>,
    gmt: Vec<Wide>,
    gmb: Vec<Wide>,
}

#[derive(Clone, Debug, Default)]
struct ParticularAdjoint {
    ap: Vec<Wide>,
    am: Vec<Wide>,
    at: Vec<Wide>,
    gpt: Vec<Wide>,
    gpb: Vec<Wide>,
    gmt: Vec<Wide>,
    gmb: Vec<Wide>,
}

#[derive(Clone, Debug, Default)]
struct Bvp {
    e: Vec<Wide>,
    c: Vec<Wide>,
    d: Vec<Wide>,
    a: Vec<Wide>,
    b: Vec<Wide>,
    rhs: Vec<Wide>,
    inv_mu: Vec<Wide>,
    alpha: Vec<Wide>,
    beta: Vec<Wide>,
    gamma: Vec<Wide>,
    z: Vec<Wide>,
}

#[derive(Clone, Debug, Default)]
pub(super) struct ExplicitWorkspace {
    od: Vec<Wide>,
    ssa: Vec<Wide>,
    b1: Vec<Wide>,
    transmission: Vec<Wide>,
    secant: Vec<Wide>,
    thermal_b0: Vec<Wide>,
    thermal_b1: Vec<Wide>,
    d_od: Vec<Wide>,
    d_ssa: Vec<Wide>,
    d_b1: Vec<Wide>,
    d_transmission: Vec<Wide>,
    d_secant: Vec<Wide>,
    d_thermal_b0: Vec<Wide>,
    d_thermal_b1: Vec<Wide>,
    d_level_extinction: Vec<Wide>,
    d_level_ssa: Vec<Wide>,
    d_level_b1: Vec<Wide>,
    d_level_emission: Vec<Wide>,
    beam: Vec<Wide>,
    homogeneous: [Homogeneous; 2],
    d_homogeneous: [HomogeneousAdjoint; 2],
    particular: [Particular; 2],
    d_particular: [ParticularAdjoint; 2],
    bvp: [Bvp; 2],
    d_solution: [Vec<Wide>; 2],
    source: Vec<Wide>,
    attenuation: Vec<Wide>,
}

#[derive(Clone, Debug, Default)]
pub(super) struct SphericalExplicitWorkspace {
    columns: Vec<ExplicitWorkspace>,
    path_radiance: Vec<Wide>,
    path_transmission: Vec<Wide>,
    path_source: Vec<Wide>,
    d_level_extinction: Vec<Wide>,
    d_level_ssa: Vec<Wide>,
    d_level_b1: Vec<Wide>,
    d_level_emission: Vec<Wide>,
}

impl SphericalExplicitWorkspace {
    fn resize_columns<const SOLAR: bool>(&mut self, columns: usize, layers: usize) {
        self.columns
            .resize_with(columns, ExplicitWorkspace::default);
        for column in &mut self.columns {
            column.resize::<SOLAR>(layers);
        }
    }

    fn resize_path(&mut self, segments: usize) {
        resize_for_overwrite(&mut self.path_radiance, segments + 1);
        resize_for_overwrite(&mut self.path_transmission, segments);
        resize_for_overwrite(&mut self.path_source, segments);
    }

    fn zero_level_adjoints<const SOLAR: bool>(&mut self, levels: usize) {
        resize_and_zero(&mut self.d_level_extinction, levels);
        resize_and_zero(&mut self.d_level_ssa, levels);
        resize_and_zero(&mut self.d_level_b1, levels);
        if !SOLAR {
            resize_and_zero(&mut self.d_level_emission, levels);
        }
    }

    pub(super) fn capacity_bytes(&self) -> usize {
        let wide_capacity = self.path_radiance.capacity()
            + self.path_transmission.capacity()
            + self.path_source.capacity()
            + self.d_level_extinction.capacity()
            + self.d_level_ssa.capacity()
            + self.d_level_b1.capacity()
            + self.d_level_emission.capacity();
        wide_capacity * std::mem::size_of::<Wide>()
            + self
                .columns
                .iter()
                .map(ExplicitWorkspace::capacity_bytes)
                .sum::<usize>()
    }
}

fn resize_for_overwrite(values: &mut Vec<Wide>, len: usize) {
    values.resize(len, Wide::splat(0.0));
}

fn resize_and_zero(values: &mut Vec<Wide>, len: usize) {
    resize_for_overwrite(values, len);
    values.fill(Wide::splat(0.0));
}

impl ExplicitWorkspace {
    fn resize<const SOLAR: bool>(&mut self, n: usize) {
        for values in [
            &mut self.od,
            &mut self.ssa,
            &mut self.b1,
            &mut self.beam,
            &mut self.source,
        ] {
            resize_for_overwrite(values, n);
        }
        for values in [&mut self.d_od, &mut self.d_ssa, &mut self.d_b1] {
            resize_and_zero(values, n);
        }
        if SOLAR {
            resize_for_overwrite(&mut self.secant, n);
            resize_and_zero(&mut self.d_secant, n);
            resize_for_overwrite(&mut self.transmission, n + 1);
            resize_and_zero(&mut self.d_transmission, n + 1);
        } else {
            for values in [&mut self.thermal_b0, &mut self.thermal_b1] {
                resize_for_overwrite(values, n);
            }
            for values in [&mut self.d_thermal_b0, &mut self.d_thermal_b1] {
                resize_and_zero(values, n);
            }
        }
        resize_for_overwrite(&mut self.attenuation, n + 1);
        resize_and_zero(&mut self.d_level_extinction, n + 1);
        resize_and_zero(&mut self.d_level_ssa, n + 1);
        resize_and_zero(&mut self.d_level_b1, n + 1);
        if !SOLAR {
            resize_and_zero(&mut self.d_level_emission, n + 1);
        }
        let naz = if SOLAR { 2 } else { 1 };
        for az in 0..naz {
            let h = &mut self.homogeneous[az];
            for values in [
                &mut h.d,
                &mut h.s,
                &mut h.k,
                &mut h.xp,
                &mut h.xm,
                &mut h.omega,
                &mut h.norm,
            ] {
                resize_for_overwrite(values, n);
            }
            let dh = &mut self.d_homogeneous[az];
            for values in [&mut dh.k, &mut dh.xp, &mut dh.xm, &mut dh.omega] {
                resize_and_zero(values, n);
            }
            let p = &mut self.particular[az];
            if SOLAR {
                for values in [&mut p.qp, &mut p.qm, &mut p.ap, &mut p.am] {
                    resize_for_overwrite(values, n);
                }
            } else {
                resize_for_overwrite(&mut p.at, n);
            }
            for values in [
                &mut p.exponential,
                &mut p.cp,
                &mut p.cm,
                &mut p.gpt,
                &mut p.gpb,
                &mut p.gmt,
                &mut p.gmb,
            ] {
                resize_for_overwrite(values, n);
            }
            let dp = &mut self.d_particular[az];
            if SOLAR {
                resize_and_zero(&mut dp.ap, n);
                resize_and_zero(&mut dp.am, n);
            } else {
                resize_and_zero(&mut dp.at, n);
            }
            for values in [&mut dp.gpt, &mut dp.gpb, &mut dp.gmt, &mut dp.gmb] {
                resize_and_zero(values, n);
            }
            resize_bvp_for_overwrite(&mut self.bvp[az], 2 * n);
            resize_and_zero(&mut self.d_solution[az], 2 * n);
        }
    }

    pub(super) fn capacity_bytes(&self) -> usize {
        let mut capacity = self.od.capacity()
            + self.ssa.capacity()
            + self.b1.capacity()
            + self.transmission.capacity()
            + self.secant.capacity()
            + self.thermal_b0.capacity()
            + self.thermal_b1.capacity()
            + self.d_od.capacity()
            + self.d_ssa.capacity()
            + self.d_b1.capacity()
            + self.d_transmission.capacity()
            + self.d_secant.capacity()
            + self.d_thermal_b0.capacity()
            + self.d_thermal_b1.capacity()
            + self.d_level_extinction.capacity()
            + self.d_level_ssa.capacity()
            + self.d_level_b1.capacity()
            + self.d_level_emission.capacity()
            + self.beam.capacity()
            + self.source.capacity()
            + self.attenuation.capacity();
        for az in 0..2 {
            let h = &self.homogeneous[az];
            capacity += h.d.capacity()
                + h.s.capacity()
                + h.k.capacity()
                + h.xp.capacity()
                + h.xm.capacity()
                + h.omega.capacity()
                + h.norm.capacity();
            let dh = &self.d_homogeneous[az];
            capacity += dh.k.capacity() + dh.xp.capacity() + dh.xm.capacity() + dh.omega.capacity();
            let p = &self.particular[az];
            capacity += p.qp.capacity()
                + p.qm.capacity()
                + p.ap.capacity()
                + p.am.capacity()
                + p.at.capacity()
                + p.exponential.capacity()
                + p.cp.capacity()
                + p.cm.capacity()
                + p.gpt.capacity()
                + p.gpb.capacity()
                + p.gmt.capacity()
                + p.gmb.capacity();
            let dp = &self.d_particular[az];
            capacity += dp.ap.capacity()
                + dp.am.capacity()
                + dp.at.capacity()
                + dp.gpt.capacity()
                + dp.gpb.capacity()
                + dp.gmt.capacity()
                + dp.gmb.capacity();
            capacity += bvp_capacity(&self.bvp[az]) + self.d_solution[az].capacity();
        }
        capacity * std::mem::size_of::<Wide>()
    }
}

fn bvp_capacity(bvp: &Bvp) -> usize {
    bvp.e.capacity()
        + bvp.c.capacity()
        + bvp.d.capacity()
        + bvp.a.capacity()
        + bvp.b.capacity()
        + bvp.rhs.capacity()
        + bvp.inv_mu.capacity()
        + bvp.alpha.capacity()
        + bvp.beta.capacity()
        + bvp.gamma.capacity()
        + bvp.z.capacity()
}

fn resize_bvp_for_overwrite(bvp: &mut Bvp, n: usize) {
    for values in [
        &mut bvp.e,
        &mut bvp.c,
        &mut bvp.d,
        &mut bvp.a,
        &mut bvp.b,
        &mut bvp.rhs,
        &mut bvp.inv_mu,
        &mut bvp.alpha,
        &mut bvp.beta,
        &mut bvp.gamma,
        &mut bvp.z,
    ] {
        resize_for_overwrite(values, n);
    }
}

pub(super) fn solve_vjp(
    geometry: &Geometry,
    mode: SourceMode,
    inputs: &LayerInputs,
    views: &[View],
    cotangent: &[f64],
    execution: ExecutionPolicy,
    workspace: &mut ExplicitWorkspace,
) -> (RadianceBatch, LayerAdjoints) {
    let nw = inputs.num_wavelengths;
    if execution == ExecutionPolicy::Serial || nw < 2 * LANES {
        return solve_range(geometry, mode, inputs, views, cotangent, 0, nw, workspace);
    }

    use rayon::prelude::*;
    let threads = rayon::current_num_threads().min(nw.div_ceil(LANES));
    if threads == 1 {
        return solve_range(geometry, mode, inputs, views, cotangent, 0, nw, workspace);
    }
    let tile_size = nw.div_ceil(threads * LANES) * LANES;
    let ranges: Vec<_> = (0..nw)
        .step_by(tile_size)
        .map(|start| (start, (nw - start).min(tile_size)))
        .collect();
    let parts: Vec<_> = ranges
        .into_par_iter()
        .map_init(ExplicitWorkspace::default, |local, (start, len)| {
            (
                start,
                solve_range(geometry, mode, inputs, views, cotangent, start, len, local),
            )
        })
        .collect();
    let (mut radiance, mut adjoints) = empty_outputs(mode, inputs.num_layers, views.len(), nw);
    for (start, (part_radiance, part_adjoints)) in parts {
        merge_outputs(
            mode,
            inputs.num_layers,
            start,
            nw,
            &part_radiance,
            &part_adjoints,
            &mut radiance,
            &mut adjoints,
        );
    }
    (radiance, adjoints)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn solve_range(
    geometry: &Geometry,
    mode: SourceMode,
    inputs: &LayerInputs,
    views: &[View],
    cotangent: &[f64],
    input_start: usize,
    len: usize,
    workspace: &mut ExplicitWorkspace,
) -> (RadianceBatch, LayerAdjoints) {
    let (mut radiance, mut adjoints) = empty_outputs(mode, inputs.num_layers, views.len(), len);
    {
        let mut radiance_rows = OutputRows {
            rows: radiance.values.chunks_exact_mut(len).collect(),
        };
        for local_base in (0..len).step_by(LANES) {
            solve_chunk(
                geometry,
                mode,
                inputs,
                views,
                cotangent,
                input_start + local_base,
                local_base,
                len,
                &mut radiance_rows,
                &mut adjoints,
                workspace,
            );
        }
    }
    (radiance, adjoints)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn solve_atmosphere_range(
    geometry: &Geometry,
    mode: SourceMode,
    atmosphere: &AtmosphereBatch,
    inputs: &LayerInputs,
    views: &[View],
    cotangent: &[f64],
    input_start: usize,
    len: usize,
    workspace: &mut ExplicitWorkspace,
) -> (RadianceBatch, AtmosphereAdjoints) {
    let (mut radiance, mut adjoints) =
        empty_atmosphere_outputs(mode, inputs.num_layers, views.len(), len);
    let ranges = [(0, len)];
    let mut outputs = split_atmosphere_outputs(&ranges, &mut radiance, &mut adjoints);
    solve_atmosphere_tile(
        geometry,
        mode,
        atmosphere,
        inputs,
        views,
        cotangent,
        input_start,
        len,
        &mut outputs[0],
        workspace,
    );
    drop(outputs);
    (radiance, adjoints)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn solve_atmosphere_tile(
    geometry: &Geometry,
    mode: SourceMode,
    atmosphere: &AtmosphereBatch,
    inputs: &LayerInputs,
    views: &[View],
    cotangent: &[f64],
    input_start: usize,
    len: usize,
    outputs: &mut AtmosphereOutputTile<'_>,
    workspace: &mut ExplicitWorkspace,
) {
    for local_base in (0..len).step_by(LANES) {
        solve_atmosphere_chunk(
            geometry,
            mode,
            atmosphere,
            inputs,
            views,
            cotangent,
            input_start + local_base,
            local_base,
            outputs,
            workspace,
        );
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn solve_unprepared_atmosphere_tile(
    geometry: &Geometry,
    mode: SourceMode,
    atmosphere: &AtmosphereBatch,
    views: &[View],
    cotangent: &[f64],
    input_start: usize,
    len: usize,
    outputs: &mut AtmosphereOutputTile<'_>,
    workspace: &mut ExplicitWorkspace,
) {
    for local_base in (0..len).step_by(LANES) {
        solve_unprepared_atmosphere_chunk(
            geometry,
            mode,
            atmosphere,
            views,
            cotangent,
            input_start + local_base,
            local_base,
            outputs,
            workspace,
        );
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn solve_unprepared_atmosphere_forward_tile(
    geometry: &Geometry,
    mode: SourceMode,
    atmosphere: &AtmosphereBatch,
    views: &[View],
    input_start: usize,
    len: usize,
    radiance: &mut OutputRows<'_>,
    workspace: &mut ExplicitWorkspace,
) {
    for local_base in (0..len).step_by(LANES) {
        match mode {
            SourceMode::Solar => solve_unprepared_atmosphere_forward_chunk::<true>(
                geometry,
                atmosphere,
                views,
                input_start + local_base,
                local_base,
                radiance,
                workspace,
            ),
            SourceMode::Thermal => solve_unprepared_atmosphere_forward_chunk::<false>(
                geometry,
                atmosphere,
                views,
                input_start + local_base,
                local_base,
                radiance,
                workspace,
            ),
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn solve_unprepared_atmosphere_jacobian_tile(
    geometry: &Geometry,
    mode: SourceMode,
    atmosphere: &AtmosphereBatch,
    views: &[View],
    input_start: usize,
    len: usize,
    outputs: &mut AtmosphereOutputTile<'_>,
    workspace: &mut ExplicitWorkspace,
) {
    for local_base in (0..len).step_by(LANES) {
        match mode {
            SourceMode::Solar => solve_unprepared_atmosphere_jacobian_chunk::<true>(
                geometry,
                atmosphere,
                views,
                input_start + local_base,
                local_base,
                outputs,
                workspace,
            ),
            SourceMode::Thermal => solve_unprepared_atmosphere_jacobian_chunk::<false>(
                geometry,
                atmosphere,
                views,
                input_start + local_base,
                local_base,
                outputs,
                workspace,
            ),
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn solve_unprepared_spherical_forward_tile(
    spherical: &SphericalGeometry,
    mode: SourceMode,
    atmosphere: &AtmosphereBatch,
    input_start: usize,
    len: usize,
    radiance: &mut OutputRows<'_>,
    workspace: &mut SphericalExplicitWorkspace,
) {
    for local_base in (0..len).step_by(LANES) {
        match mode {
            SourceMode::Solar => solve_spherical_forward_chunk::<true>(
                spherical,
                atmosphere,
                input_start + local_base,
                local_base,
                radiance,
                workspace,
            ),
            SourceMode::Thermal => solve_spherical_forward_chunk::<false>(
                spherical,
                atmosphere,
                input_start + local_base,
                local_base,
                radiance,
                workspace,
            ),
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn solve_unprepared_spherical_jacobian_tile(
    spherical: &SphericalGeometry,
    mode: SourceMode,
    atmosphere: &AtmosphereBatch,
    input_start: usize,
    len: usize,
    outputs: &mut AtmosphereOutputTile<'_>,
    workspace: &mut SphericalExplicitWorkspace,
) {
    for local_base in (0..len).step_by(LANES) {
        match mode {
            SourceMode::Solar => solve_spherical_jacobian_chunk::<true>(
                spherical,
                atmosphere,
                input_start + local_base,
                local_base,
                outputs,
                workspace,
            ),
            SourceMode::Thermal => solve_spherical_jacobian_chunk::<false>(
                spherical,
                atmosphere,
                input_start + local_base,
                local_base,
                outputs,
                workspace,
            ),
        }
    }
}

fn empty_atmosphere_outputs(
    mode: SourceMode,
    n: usize,
    num_views: usize,
    nw: usize,
) -> (RadianceBatch, AtmosphereAdjoints) {
    let radiance = RadianceBatch {
        num_views,
        num_wavelengths: nw,
        values: vec![0.0; num_views * nw],
    };
    let levels = (n + 1) * nw;
    let adjoints = AtmosphereAdjoints {
        extinction: vec![0.0; levels],
        single_scatter_albedo: vec![0.0; levels],
        first_legendre: vec![0.0; levels],
        emission: (mode == SourceMode::Thermal).then(|| vec![0.0; levels]),
        surface_albedo: vec![0.0; nw],
        surface_emission: (mode == SourceMode::Thermal).then(|| vec![0.0; nw]),
    };
    (radiance, adjoints)
}

fn empty_outputs(
    mode: SourceMode,
    n: usize,
    num_views: usize,
    nw: usize,
) -> (RadianceBatch, LayerAdjoints) {
    let radiance = RadianceBatch {
        num_views,
        num_wavelengths: nw,
        values: vec![0.0; num_views * nw],
    };
    let adjoints = LayerAdjoints {
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
    (radiance, adjoints)
}

#[allow(clippy::too_many_arguments)]
fn merge_outputs(
    mode: SourceMode,
    n: usize,
    start: usize,
    nw: usize,
    part_radiance: &RadianceBatch,
    part: &LayerAdjoints,
    radiance: &mut RadianceBatch,
    output: &mut LayerAdjoints,
) {
    let len = part_radiance.num_wavelengths;
    copy_rows(&part_radiance.values, len, &mut radiance.values, nw, start);
    copy_rows(
        &part.optical_depth,
        len,
        &mut output.optical_depth,
        nw,
        start,
    );
    copy_rows(
        &part.single_scatter_albedo,
        len,
        &mut output.single_scatter_albedo,
        nw,
        start,
    );
    copy_rows(
        &part.first_legendre,
        len,
        &mut output.first_legendre,
        nw,
        start,
    );
    output.surface_albedo[start..start + len].copy_from_slice(&part.surface_albedo);
    if mode == SourceMode::Solar {
        copy_rows(
            part.transmission.as_ref().unwrap(),
            len,
            output.transmission.as_mut().unwrap(),
            nw,
            start,
        );
        copy_rows(
            part.average_secant.as_ref().unwrap(),
            len,
            output.average_secant.as_mut().unwrap(),
            nw,
            start,
        );
    } else {
        copy_rows(
            part.thermal_b0.as_ref().unwrap(),
            len,
            output.thermal_b0.as_mut().unwrap(),
            nw,
            start,
        );
        copy_rows(
            part.thermal_b1.as_ref().unwrap(),
            len,
            output.thermal_b1.as_mut().unwrap(),
            nw,
            start,
        );
        output.surface_emission.as_mut().unwrap()[start..start + len]
            .copy_from_slice(part.surface_emission.as_ref().unwrap());
    }
    debug_assert_eq!(output.optical_depth.len(), n * nw);
}

fn copy_rows(
    source: &[f64],
    source_stride: usize,
    target: &mut [f64],
    target_stride: usize,
    start: usize,
) {
    for (row, values) in source.chunks_exact(source_stride).enumerate() {
        target[row * target_stride + start..row * target_stride + start + source_stride]
            .copy_from_slice(values);
    }
}

pub(super) struct OutputRows<'a> {
    rows: Vec<&'a mut [f64]>,
}

impl OutputRows<'_> {
    fn write(&mut self, row: usize, base: usize, value: Wide) {
        let target = &mut self.rows[row];
        let values = value.to_array();
        if base + LANES <= target.len() {
            target[base..base + LANES].copy_from_slice(&values);
            return;
        }
        for (target, value) in target[base..].iter_mut().zip(values) {
            *target = value;
        }
    }
}

pub(super) struct AtmosphereOutputTile<'a> {
    radiance: OutputRows<'a>,
    extinction: OutputRows<'a>,
    single_scatter_albedo: OutputRows<'a>,
    first_legendre: OutputRows<'a>,
    emission: Option<OutputRows<'a>>,
    surface_albedo: OutputRows<'a>,
    surface_emission: Option<OutputRows<'a>>,
}

fn split_output_rows<'a>(
    values: &'a mut [f64],
    stride: usize,
    ranges: &[(usize, usize)],
) -> Vec<OutputRows<'a>> {
    let row_count = values.len() / stride;
    let mut tiles: Vec<Vec<&mut [f64]>> = (0..ranges.len())
        .map(|_| Vec::with_capacity(row_count))
        .collect();
    for row in values.chunks_exact_mut(stride) {
        let mut remainder = row;
        let mut cursor = 0;
        for (tile, &(start, len)) in tiles.iter_mut().zip(ranges) {
            debug_assert_eq!(start, cursor);
            let (values, next) = remainder.split_at_mut(len);
            tile.push(values);
            remainder = next;
            cursor += len;
        }
        debug_assert!(remainder.is_empty());
    }
    tiles.into_iter().map(|rows| OutputRows { rows }).collect()
}

pub(super) fn split_atmosphere_outputs<'a>(
    ranges: &[(usize, usize)],
    radiance: &'a mut RadianceBatch,
    adjoints: &'a mut AtmosphereAdjoints,
) -> Vec<AtmosphereOutputTile<'a>> {
    let nw = radiance.num_wavelengths;
    debug_assert_eq!(ranges.iter().map(|(_, len)| len).sum::<usize>(), nw);

    let mut radiance = split_output_rows(&mut radiance.values, nw, ranges).into_iter();
    let mut extinction = split_output_rows(&mut adjoints.extinction, nw, ranges).into_iter();
    let mut single_scatter_albedo =
        split_output_rows(&mut adjoints.single_scatter_albedo, nw, ranges).into_iter();
    let mut first_legendre =
        split_output_rows(&mut adjoints.first_legendre, nw, ranges).into_iter();
    let mut emission = adjoints
        .emission
        .as_mut()
        .map(|values| split_output_rows(values, nw, ranges).into_iter());
    let mut surface_albedo =
        split_output_rows(&mut adjoints.surface_albedo, nw, ranges).into_iter();
    let mut surface_emission = adjoints
        .surface_emission
        .as_mut()
        .map(|values| split_output_rows(values, nw, ranges).into_iter());

    (0..ranges.len())
        .map(|_| AtmosphereOutputTile {
            radiance: radiance.next().unwrap(),
            extinction: extinction.next().unwrap(),
            single_scatter_albedo: single_scatter_albedo.next().unwrap(),
            first_legendre: first_legendre.next().unwrap(),
            emission: emission.as_mut().map(|rows| rows.next().unwrap()),
            surface_albedo: surface_albedo.next().unwrap(),
            surface_emission: surface_emission.as_mut().map(|rows| rows.next().unwrap()),
        })
        .collect()
}

pub(super) fn split_radiance_outputs<'a>(
    ranges: &[(usize, usize)],
    radiance: &'a mut RadianceBatch,
) -> Vec<OutputRows<'a>> {
    split_output_rows(&mut radiance.values, radiance.num_wavelengths, ranges)
}

pub(super) fn split_atmosphere_jacobian_outputs<'a>(
    ranges: &[(usize, usize)],
    radiance: &'a mut RadianceBatch,
    jacobians: &'a mut AtmosphereJacobians,
) -> Vec<AtmosphereOutputTile<'a>> {
    let nw = radiance.num_wavelengths;
    debug_assert_eq!(jacobians.num_wavelengths, nw);
    debug_assert_eq!(jacobians.num_views, radiance.num_views);

    let mut radiance = split_output_rows(&mut radiance.values, nw, ranges).into_iter();
    let mut extinction = split_output_rows(&mut jacobians.extinction, nw, ranges).into_iter();
    let mut single_scatter_albedo =
        split_output_rows(&mut jacobians.single_scatter_albedo, nw, ranges).into_iter();
    let mut first_legendre =
        split_output_rows(&mut jacobians.first_legendre, nw, ranges).into_iter();
    let mut emission = jacobians
        .emission
        .as_mut()
        .map(|values| split_output_rows(values, nw, ranges).into_iter());
    let mut surface_albedo =
        split_output_rows(&mut jacobians.surface_albedo, nw, ranges).into_iter();
    let mut surface_emission = jacobians
        .surface_emission
        .as_mut()
        .map(|values| split_output_rows(values, nw, ranges).into_iter());

    (0..ranges.len())
        .map(|_| AtmosphereOutputTile {
            radiance: radiance.next().unwrap(),
            extinction: extinction.next().unwrap(),
            single_scatter_albedo: single_scatter_albedo.next().unwrap(),
            first_legendre: first_legendre.next().unwrap(),
            emission: emission.as_mut().map(|rows| rows.next().unwrap()),
            surface_albedo: surface_albedo.next().unwrap(),
            surface_emission: surface_emission.as_mut().map(|rows| rows.next().unwrap()),
        })
        .collect()
}

fn write_wide(target: &mut [f64], row: usize, base: usize, stride: usize, value: Wide) {
    let base = row * stride + base;
    let values = value.to_array();
    let remaining = stride - base % stride;
    if LANES <= remaining {
        target[base..base + LANES].copy_from_slice(&values);
        return;
    }
    for (target, value) in target[base..base + remaining].iter_mut().zip(values) {
        *target = value;
    }
}

#[allow(clippy::too_many_arguments)]
fn run_chunk(
    geometry: &Geometry,
    mode: SourceMode,
    inputs: &LayerInputs,
    views: &[View],
    cotangent: &[f64],
    input_base: usize,
    output_base: usize,
    radiance: &mut OutputRows<'_>,
    w: &mut ExplicitWorkspace,
) -> (Wide, Wide) {
    match mode {
        SourceMode::Solar => run_chunk_mode::<true>(
            geometry,
            inputs,
            views,
            cotangent,
            input_base,
            output_base,
            radiance,
            w,
        ),
        SourceMode::Thermal => run_chunk_mode::<false>(
            geometry,
            inputs,
            views,
            cotangent,
            input_base,
            output_base,
            radiance,
            w,
        ),
    }
}

#[allow(clippy::too_many_arguments)]
fn run_chunk_mode<const SOLAR: bool>(
    geometry: &Geometry,
    inputs: &LayerInputs,
    views: &[View],
    cotangent: &[f64],
    input_base: usize,
    output_base: usize,
    radiance: &mut OutputRows<'_>,
    w: &mut ExplicitWorkspace,
) -> (Wide, Wide) {
    let n = inputs.num_layers;
    let input_stride = inputs.num_wavelengths;
    w.resize::<SOLAR>(n);
    for layer in 0..n {
        w.od[layer] = Wide::load_row(&inputs.optical_depth, layer, input_base, input_stride);
        w.ssa[layer] = Wide::load_row(
            &inputs.single_scatter_albedo,
            layer,
            input_base,
            input_stride,
        );
        w.b1[layer] = Wide::load_row(&inputs.first_legendre, layer, input_base, input_stride);
        if SOLAR {
            w.secant[layer] = Wide::load_row(
                inputs.average_secant.as_ref().unwrap(),
                layer,
                input_base,
                input_stride,
            );
        } else {
            w.thermal_b0[layer] = Wide::load_row(
                inputs.thermal_b0.as_ref().unwrap(),
                layer,
                input_base,
                input_stride,
            );
            w.thermal_b1[layer] = Wide::load_row(
                inputs.thermal_b1.as_ref().unwrap(),
                layer,
                input_base,
                input_stride,
            );
        }
    }
    if SOLAR {
        for boundary in 0..=n {
            w.transmission[boundary] = Wide::load_row(
                inputs.transmission.as_ref().unwrap(),
                boundary,
                input_base,
                input_stride,
            );
        }
    }
    let albedo = Wide::load(&inputs.surface_albedo, input_base, input_stride);
    let thermal_surface = if SOLAR {
        Wide::splat(0.0)
    } else {
        Wide::load(
            inputs.surface_emission.as_ref().unwrap(),
            input_base,
            input_stride,
        )
    };

    run_loaded_chunk_mode::<SOLAR>(
        geometry,
        albedo,
        thermal_surface,
        views,
        cotangent,
        input_base,
        input_stride,
        output_base,
        radiance,
        w,
    )
}

#[allow(clippy::too_many_arguments)]
fn run_loaded_chunk_mode<const SOLAR: bool>(
    geometry: &Geometry,
    albedo: Wide,
    thermal_surface: Wide,
    views: &[View],
    cotangent: &[f64],
    input_base: usize,
    input_stride: usize,
    output_base: usize,
    radiance: &mut OutputRows<'_>,
    w: &mut ExplicitWorkspace,
) -> (Wide, Wide) {
    forward_layers::<SOLAR>(geometry, albedo, thermal_surface, w);
    let (mut d_albedo, mut d_thermal_surface) = reverse_views::<SOLAR>(
        geometry,
        albedo,
        thermal_surface,
        views,
        ViewSeed::Cotangent(cotangent),
        input_base,
        input_stride,
        output_base,
        0,
        radiance,
        w,
    );
    let (bvp_albedo, bvp_thermal_surface) =
        reverse_bvp::<SOLAR>(geometry, albedo, thermal_surface, w);
    d_albedo += bvp_albedo;
    d_thermal_surface += bvp_thermal_surface;
    reverse_layers::<SOLAR>(geometry, w);

    (d_albedo, d_thermal_surface)
}

#[derive(Clone, Copy)]
enum ViewSeed<'a> {
    Cotangent(&'a [f64]),
    Unit,
    None,
}

#[allow(clippy::too_many_arguments)]
fn run_atmosphere_chunk(
    geometry: &Geometry,
    mode: SourceMode,
    atmosphere: &AtmosphereBatch,
    views: &[View],
    cotangent: &[f64],
    input_base: usize,
    output_base: usize,
    radiance: &mut OutputRows<'_>,
    w: &mut ExplicitWorkspace,
) -> (Wide, Wide) {
    match mode {
        SourceMode::Solar => run_atmosphere_chunk_mode::<true>(
            geometry,
            atmosphere,
            views,
            cotangent,
            input_base,
            output_base,
            radiance,
            w,
        ),
        SourceMode::Thermal => run_atmosphere_chunk_mode::<false>(
            geometry,
            atmosphere,
            views,
            cotangent,
            input_base,
            output_base,
            radiance,
            w,
        ),
    }
}

#[allow(clippy::too_many_arguments)]
fn run_atmosphere_chunk_mode<const SOLAR: bool>(
    geometry: &Geometry,
    atmosphere: &AtmosphereBatch,
    views: &[View],
    cotangent: &[f64],
    input_base: usize,
    output_base: usize,
    radiance: &mut OutputRows<'_>,
    w: &mut ExplicitWorkspace,
) -> (Wide, Wide) {
    let input_stride = atmosphere.num_wavelengths;
    let (albedo, thermal_surface) =
        load_atmosphere_chunk_mode::<SOLAR>(geometry, atmosphere, input_base, w);
    run_loaded_chunk_mode::<SOLAR>(
        geometry,
        albedo,
        thermal_surface,
        views,
        cotangent,
        input_base,
        input_stride,
        output_base,
        radiance,
        w,
    )
}

fn load_atmosphere_chunk_mode<const SOLAR: bool>(
    geometry: &Geometry,
    atmosphere: &AtmosphereBatch,
    input_base: usize,
    w: &mut ExplicitWorkspace,
) -> (Wide, Wide) {
    let n = geometry.num_layers();
    let input_stride = atmosphere.num_wavelengths;
    w.resize::<SOLAR>(n);
    let mut top_extinction = Wide::load(&atmosphere.extinction, input_base, input_stride);
    let mut top_ssa = Wide::load(&atmosphere.single_scatter_albedo, input_base, input_stride);
    let mut top_b1 = Wide::load(&atmosphere.first_legendre, input_base, input_stride);
    let mut top_emission = if SOLAR {
        Wide::splat(0.0)
    } else {
        Wide::load(
            atmosphere.emission.as_ref().unwrap(),
            input_base,
            input_stride,
        )
    };
    for layer in 0..n {
        let bottom_extinction =
            Wide::load_row(&atmosphere.extinction, layer + 1, input_base, input_stride);
        let bottom_ssa = Wide::load_row(
            &atmosphere.single_scatter_albedo,
            layer + 1,
            input_base,
            input_stride,
        );
        let bottom_b1 = Wide::load_row(
            &atmosphere.first_legendre,
            layer + 1,
            input_base,
            input_stride,
        );
        let scattering_top = top_extinction * top_ssa;
        let scattering_bottom = bottom_extinction * bottom_ssa;
        let avg_extinction = 0.5 * (top_extinction + bottom_extinction);
        let avg_scattering = 0.5 * (scattering_top + scattering_bottom);
        w.od[layer] = avg_extinction * geometry.layer_thickness[layer];
        w.ssa[layer] = (avg_scattering / avg_extinction).min(Wide::splat(1.0 - 1.0e-9));
        w.b1[layer] =
            0.5 * (scattering_top * top_b1 + scattering_bottom * bottom_b1) / avg_scattering;

        if !SOLAR {
            let emission = atmosphere.emission.as_ref().unwrap();
            let bottom = Wide::load_row(emission, layer + 1, input_base, input_stride);
            w.thermal_b0[layer] = top_emission;
            w.thermal_b1[layer] = -(bottom / top_emission).ln() / geometry.layer_thickness[layer];
            top_emission = bottom;
        }
        top_extinction = bottom_extinction;
        top_ssa = bottom_ssa;
        top_b1 = bottom_b1;
    }

    if SOLAR {
        let irradiance = Wide::load(
            atmosphere.solar_irradiance.as_ref().unwrap(),
            input_base,
            input_stride,
        );
        w.attenuation[0] = Wide::splat(0.0);
        w.transmission[0] = irradiance;
        for boundary in 0..n {
            let mut slant = Wide::splat(0.0);
            for layer in 0..n {
                let factor = geometry.chapman_factors[boundary * n + layer];
                if factor != 0.0 {
                    slant += w.od[layer] * factor;
                }
            }
            w.attenuation[boundary + 1] = -slant;
            w.transmission[boundary + 1] = (-slant).exp() * irradiance;
        }
        for layer in 0..n {
            w.secant[layer] = (w.attenuation[layer] - w.attenuation[layer + 1]) / w.od[layer];
        }
    }

    let albedo = Wide::load(&atmosphere.surface_albedo, input_base, input_stride);
    let thermal_surface = if SOLAR {
        Wide::splat(0.0)
    } else {
        Wide::load(
            atmosphere.surface_emission.as_ref().unwrap(),
            input_base,
            input_stride,
        )
    };
    (albedo, thermal_surface)
}

#[allow(clippy::too_many_arguments)]
fn solve_unprepared_atmosphere_forward_chunk<const SOLAR: bool>(
    geometry: &Geometry,
    atmosphere: &AtmosphereBatch,
    views: &[View],
    input_base: usize,
    output_base: usize,
    radiance: &mut OutputRows<'_>,
    w: &mut ExplicitWorkspace,
) {
    let (albedo, thermal_surface) =
        load_atmosphere_chunk_mode::<SOLAR>(geometry, atmosphere, input_base, w);
    forward_layers::<SOLAR>(geometry, albedo, thermal_surface, w);
    reverse_views::<SOLAR>(
        geometry,
        albedo,
        thermal_surface,
        views,
        ViewSeed::None,
        input_base,
        atmosphere.num_wavelengths,
        output_base,
        0,
        radiance,
        w,
    );
}

#[allow(clippy::too_many_arguments)]
fn solve_unprepared_atmosphere_jacobian_chunk<const SOLAR: bool>(
    geometry: &Geometry,
    atmosphere: &AtmosphereBatch,
    views: &[View],
    input_base: usize,
    output_base: usize,
    outputs: &mut AtmosphereOutputTile<'_>,
    w: &mut ExplicitWorkspace,
) {
    let n = geometry.num_layers();
    let (albedo, thermal_surface) =
        load_atmosphere_chunk_mode::<SOLAR>(geometry, atmosphere, input_base, w);
    forward_layers::<SOLAR>(geometry, albedo, thermal_surface, w);

    for (view_index, view) in views.iter().enumerate() {
        // Preserve the shared forward/BVP state while clearing only adjoints
        // from the preceding view.
        w.resize::<SOLAR>(n);
        let (mut d_albedo, mut d_thermal_surface) = reverse_views::<SOLAR>(
            geometry,
            albedo,
            thermal_surface,
            std::slice::from_ref(view),
            ViewSeed::Unit,
            input_base,
            atmosphere.num_wavelengths,
            output_base,
            view_index,
            &mut outputs.radiance,
            w,
        );
        let (bvp_albedo, bvp_thermal_surface) =
            reverse_bvp::<SOLAR>(geometry, albedo, thermal_surface, w);
        d_albedo += bvp_albedo;
        d_thermal_surface += bvp_thermal_surface;
        reverse_layers::<SOLAR>(geometry, w);

        map_chunk_to_atmosphere(
            geometry,
            if SOLAR {
                SourceMode::Solar
            } else {
                SourceMode::Thermal
            },
            atmosphere,
            input_base,
            output_base,
            d_albedo,
            d_thermal_surface,
            view_index * (n + 1),
            view_index,
            outputs,
            w,
        );
    }
}

fn prepare_spherical_columns<const SOLAR: bool>(
    spherical: &SphericalGeometry,
    atmosphere: &AtmosphereBatch,
    input_base: usize,
    workspace: &mut SphericalExplicitWorkspace,
) -> (Wide, Wide) {
    let n = spherical.columns[0].num_layers();
    workspace.resize_columns::<SOLAR>(spherical.columns.len(), n);
    let (first, remaining) = workspace.columns.split_first_mut().unwrap();
    let boundary =
        load_atmosphere_chunk_mode::<SOLAR>(&spherical.columns[0], atmosphere, input_base, first);
    forward_layers::<SOLAR>(&spherical.columns[0], boundary.0, boundary.1, first);

    for (geometry, column) in spherical.columns.iter().skip(1).zip(remaining) {
        column.od.copy_from_slice(&first.od);
        column.ssa.copy_from_slice(&first.ssa);
        column.b1.copy_from_slice(&first.b1);
        if SOLAR {
            column.attenuation[0] = Wide::splat(0.0);
            column.transmission[0] = first.transmission[0];
            for boundary in 0..n {
                let mut slant = Wide::splat(0.0);
                for layer in 0..n {
                    let factor = geometry.chapman_factors[boundary * n + layer];
                    if factor != 0.0 {
                        slant += column.od[layer] * factor;
                    }
                }
                column.attenuation[boundary + 1] = -slant;
                column.transmission[boundary + 1] = (-slant).exp() * column.transmission[0];
            }
            for layer in 0..n {
                column.secant[layer] =
                    (column.attenuation[layer] - column.attenuation[layer + 1]) / column.od[layer];
            }
        } else {
            column.thermal_b0.copy_from_slice(&first.thermal_b0);
            column.thermal_b1.copy_from_slice(&first.thermal_b1);
        }
        forward_layers::<SOLAR>(geometry, boundary.0, boundary.1, column);
    }
    boundary
}

fn sza_interpolation(grid: &[f64], cosine: f64) -> (usize, usize, f64) {
    if grid.len() == 1 || cosine <= grid[0] {
        return (0, 0, 0.0);
    }
    let last = grid.len() - 1;
    if cosine >= grid[last] {
        return (last, last, 0.0);
    }
    let upper = grid.partition_point(|value| *value < cosine);
    let lower = upper - 1;
    let upper_weight = (cosine - grid[lower]) / (grid[upper] - grid[lower]);
    (lower, upper, upper_weight)
}

fn interpolated_local_source<const SOLAR: bool>(
    spherical: &SphericalGeometry,
    segment: usize,
    columns: &[ExplicitWorkspace],
) -> Wide {
    let rays = &spherical.rays;
    let (lower, upper, upper_weight) =
        sza_interpolation(&spherical.sza_grid, rays.segment_cos_sza[segment]);
    let evaluate = |column: usize| {
        local_source::<SOLAR>(
            &spherical.columns[column],
            rays.segment_layers[segment],
            rays.segment_fractions[segment],
            rays.segment_cosines[segment],
            rays.segment_relative_azimuths[segment],
            &columns[column],
        )
    };
    let lower_source = evaluate(lower);
    if lower == upper {
        lower_source
    } else {
        lower_source * (1.0 - upper_weight) + evaluate(upper) * upper_weight
    }
}

fn interpolated_surface_source<const SOLAR: bool>(
    spherical: &SphericalGeometry,
    view: usize,
    albedo: Wide,
    thermal_surface: Wide,
    columns: &[ExplicitWorkspace],
) -> Wide {
    let (lower, upper, upper_weight) =
        sza_interpolation(&spherical.sza_grid, spherical.rays.ground_cos_sza[view]);
    let lower_source = surface_source::<SOLAR>(
        &spherical.columns[lower],
        albedo,
        thermal_surface,
        &columns[lower],
    );
    if lower == upper {
        lower_source
    } else {
        lower_source * (1.0 - upper_weight)
            + surface_source::<SOLAR>(
                &spherical.columns[upper],
                albedo,
                thermal_surface,
                &columns[upper],
            ) * upper_weight
    }
}

fn segment_optical_depth(
    spherical: &SphericalGeometry,
    atmosphere: &AtmosphereBatch,
    segment: usize,
    input_base: usize,
) -> Wide {
    let mut optical_depth = Wide::splat(0.0);
    let start = spherical.rays.od_offsets[segment];
    let end = spherical.rays.od_offsets[segment + 1];
    for stencil in start..end {
        optical_depth += Wide::load_row(
            &atmosphere.extinction,
            spherical.rays.od_indices[stencil],
            input_base,
            atmosphere.num_wavelengths,
        ) * spherical.rays.od_weights[stencil];
    }
    optical_depth
}

#[allow(clippy::too_many_arguments)]
fn forward_spherical_ray<const SOLAR: bool>(
    spherical: &SphericalGeometry,
    atmosphere: &AtmosphereBatch,
    view: usize,
    input_base: usize,
    albedo: Wide,
    thermal_surface: Wide,
    workspace: &mut SphericalExplicitWorkspace,
) -> Wide {
    let start = spherical.rays.ray_offsets[view];
    let end = spherical.rays.ray_offsets[view + 1];
    let segments = end - start;
    workspace.resize_path(segments);
    workspace.path_radiance[0] = if spherical.rays.ground_hit[view] {
        interpolated_surface_source::<SOLAR>(
            spherical,
            view,
            albedo,
            thermal_surface,
            &workspace.columns,
        )
    } else {
        Wide::splat(0.0)
    };
    for (local, segment) in (start..end).enumerate() {
        let transmission =
            (-segment_optical_depth(spherical, atmosphere, segment, input_base)).exp();
        let source = interpolated_local_source::<SOLAR>(spherical, segment, &workspace.columns);
        workspace.path_transmission[local] = transmission;
        workspace.path_source[local] = source;
        workspace.path_radiance[local + 1] =
            workspace.path_radiance[local] * transmission + source * (1.0 - transmission);
    }
    workspace.path_radiance[segments]
}

fn solve_spherical_forward_chunk<const SOLAR: bool>(
    spherical: &SphericalGeometry,
    atmosphere: &AtmosphereBatch,
    input_base: usize,
    output_base: usize,
    radiance: &mut OutputRows<'_>,
    workspace: &mut SphericalExplicitWorkspace,
) {
    let (albedo, thermal_surface) =
        prepare_spherical_columns::<SOLAR>(spherical, atmosphere, input_base, workspace);
    for view in 0..spherical.rays.num_views() {
        let value = forward_spherical_ray::<SOLAR>(
            spherical,
            atmosphere,
            view,
            input_base,
            albedo,
            thermal_surface,
            workspace,
        );
        radiance.write(view, output_base, value);
    }
}

#[allow(clippy::too_many_arguments)]
fn reverse_interpolated_local_source<const SOLAR: bool>(
    spherical: &SphericalGeometry,
    segment: usize,
    seed: Wide,
    columns: &mut [ExplicitWorkspace],
) {
    let rays = &spherical.rays;
    let (lower, upper, upper_weight) =
        sza_interpolation(&spherical.sza_grid, rays.segment_cos_sza[segment]);
    reverse_local_source::<SOLAR>(
        &spherical.columns[lower],
        rays.segment_layers[segment],
        rays.segment_fractions[segment],
        rays.segment_cosines[segment],
        rays.segment_relative_azimuths[segment],
        seed * (1.0 - upper_weight),
        &mut columns[lower],
    );
    if lower != upper {
        reverse_local_source::<SOLAR>(
            &spherical.columns[upper],
            rays.segment_layers[segment],
            rays.segment_fractions[segment],
            rays.segment_cosines[segment],
            rays.segment_relative_azimuths[segment],
            seed * upper_weight,
            &mut columns[upper],
        );
    }
}

fn reverse_interpolated_surface<const SOLAR: bool>(
    spherical: &SphericalGeometry,
    view: usize,
    albedo: Wide,
    seed: Wide,
    columns: &mut [ExplicitWorkspace],
) -> (Wide, Wide) {
    let (lower, upper, upper_weight) =
        sza_interpolation(&spherical.sza_grid, spherical.rays.ground_cos_sza[view]);
    let mut result = reverse_surface_source::<SOLAR>(
        &spherical.columns[lower],
        albedo,
        seed * (1.0 - upper_weight),
        &mut columns[lower],
    );
    if lower != upper {
        let upper_result = reverse_surface_source::<SOLAR>(
            &spherical.columns[upper],
            albedo,
            seed * upper_weight,
            &mut columns[upper],
        );
        result.0 += upper_result.0;
        result.1 += upper_result.1;
    }
    result
}

#[allow(clippy::too_many_arguments)]
fn solve_spherical_jacobian_chunk<const SOLAR: bool>(
    spherical: &SphericalGeometry,
    atmosphere: &AtmosphereBatch,
    input_base: usize,
    output_base: usize,
    outputs: &mut AtmosphereOutputTile<'_>,
    workspace: &mut SphericalExplicitWorkspace,
) {
    let n = spherical.columns[0].num_layers();
    let levels = n + 1;
    let (albedo, thermal_surface) =
        prepare_spherical_columns::<SOLAR>(spherical, atmosphere, input_base, workspace);

    for view in 0..spherical.rays.num_views() {
        for column in &mut workspace.columns {
            column.resize::<SOLAR>(n);
        }
        workspace.zero_level_adjoints::<SOLAR>(levels);
        let value = forward_spherical_ray::<SOLAR>(
            spherical,
            atmosphere,
            view,
            input_base,
            albedo,
            thermal_surface,
            workspace,
        );
        outputs.radiance.write(view, output_base, value);

        let start = spherical.rays.ray_offsets[view];
        let end = spherical.rays.ray_offsets[view + 1];
        let mut d_radiance = Wide::splat(1.0);
        for local in (0..end - start).rev() {
            let segment = start + local;
            let transmission = workspace.path_transmission[local];
            let source = workspace.path_source[local];
            let incoming = workspace.path_radiance[local];
            let d_transmission = d_radiance * (incoming - source);
            let d_source = d_radiance * (1.0 - transmission);
            let d_optical_depth = -d_transmission * transmission;
            let od_start = spherical.rays.od_offsets[segment];
            let od_end = spherical.rays.od_offsets[segment + 1];
            for stencil in od_start..od_end {
                workspace.d_level_extinction[spherical.rays.od_indices[stencil]] +=
                    d_optical_depth * spherical.rays.od_weights[stencil];
            }
            reverse_interpolated_local_source::<SOLAR>(
                spherical,
                segment,
                d_source,
                &mut workspace.columns,
            );
            d_radiance = d_radiance * transmission;
        }

        let mut d_albedo = Wide::splat(0.0);
        let mut d_thermal_surface = Wide::splat(0.0);
        if spherical.rays.ground_hit[view] {
            let ground = reverse_interpolated_surface::<SOLAR>(
                spherical,
                view,
                albedo,
                d_radiance,
                &mut workspace.columns,
            );
            d_albedo += ground.0;
            d_thermal_surface += ground.1;
        }

        for (geometry, column) in spherical.columns.iter().zip(&mut workspace.columns) {
            let bvp = reverse_bvp::<SOLAR>(geometry, albedo, thermal_surface, column);
            d_albedo += bvp.0;
            d_thermal_surface += bvp.1;
            reverse_layers::<SOLAR>(geometry, column);
            prepare_spherical_column_atmosphere_adjoint::<SOLAR>(
                geometry, atmosphere, input_base, column,
            );
            for level in 0..levels {
                workspace.d_level_extinction[level] += column.d_level_extinction[level];
                workspace.d_level_ssa[level] += column.d_level_ssa[level];
                workspace.d_level_b1[level] += column.d_level_b1[level];
                if !SOLAR {
                    workspace.d_level_emission[level] += column.d_level_emission[level];
                }
            }
        }

        let level_row_offset = view * levels;
        for level in 0..levels {
            outputs.extinction.write(
                level_row_offset + level,
                output_base,
                workspace.d_level_extinction[level],
            );
            outputs.single_scatter_albedo.write(
                level_row_offset + level,
                output_base,
                workspace.d_level_ssa[level],
            );
            outputs.first_legendre.write(
                level_row_offset + level,
                output_base,
                workspace.d_level_b1[level],
            );
            if !SOLAR {
                outputs.emission.as_mut().unwrap().write(
                    level_row_offset + level,
                    output_base,
                    workspace.d_level_emission[level],
                );
            }
        }
        outputs.surface_albedo.write(view, output_base, d_albedo);
        if !SOLAR {
            outputs
                .surface_emission
                .as_mut()
                .unwrap()
                .write(view, output_base, d_thermal_surface);
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn solve_chunk(
    geometry: &Geometry,
    mode: SourceMode,
    inputs: &LayerInputs,
    views: &[View],
    cotangent: &[f64],
    input_base: usize,
    output_base: usize,
    output_stride: usize,
    radiance: &mut OutputRows<'_>,
    adjoints: &mut LayerAdjoints,
    w: &mut ExplicitWorkspace,
) {
    let (d_albedo, d_thermal_surface) = run_chunk(
        geometry,
        mode,
        inputs,
        views,
        cotangent,
        input_base,
        output_base,
        radiance,
        w,
    );
    let n = inputs.num_layers;
    for layer in 0..n {
        write_wide(
            &mut adjoints.optical_depth,
            layer,
            output_base,
            output_stride,
            w.d_od[layer],
        );
        write_wide(
            &mut adjoints.single_scatter_albedo,
            layer,
            output_base,
            output_stride,
            w.d_ssa[layer],
        );
        write_wide(
            &mut adjoints.first_legendre,
            layer,
            output_base,
            output_stride,
            w.d_b1[layer],
        );
    }
    write_wide(
        &mut adjoints.surface_albedo,
        0,
        output_base,
        output_stride,
        d_albedo,
    );
    if mode == SourceMode::Solar {
        for boundary in 0..=n {
            write_wide(
                adjoints.transmission.as_mut().unwrap(),
                boundary,
                output_base,
                output_stride,
                w.d_transmission[boundary],
            );
        }
        for layer in 0..n {
            write_wide(
                adjoints.average_secant.as_mut().unwrap(),
                layer,
                output_base,
                output_stride,
                w.d_secant[layer],
            );
        }
    } else {
        for layer in 0..n {
            write_wide(
                adjoints.thermal_b0.as_mut().unwrap(),
                layer,
                output_base,
                output_stride,
                w.d_thermal_b0[layer],
            );
            write_wide(
                adjoints.thermal_b1.as_mut().unwrap(),
                layer,
                output_base,
                output_stride,
                w.d_thermal_b1[layer],
            );
        }
        write_wide(
            adjoints.surface_emission.as_mut().unwrap(),
            0,
            output_base,
            output_stride,
            d_thermal_surface,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn solve_atmosphere_chunk(
    geometry: &Geometry,
    mode: SourceMode,
    atmosphere: &AtmosphereBatch,
    inputs: &LayerInputs,
    views: &[View],
    cotangent: &[f64],
    input_base: usize,
    output_base: usize,
    outputs: &mut AtmosphereOutputTile<'_>,
    w: &mut ExplicitWorkspace,
) {
    let (d_albedo, d_thermal_surface) = run_chunk(
        geometry,
        mode,
        inputs,
        views,
        cotangent,
        input_base,
        output_base,
        &mut outputs.radiance,
        w,
    );
    map_chunk_to_atmosphere(
        geometry,
        mode,
        atmosphere,
        input_base,
        output_base,
        d_albedo,
        d_thermal_surface,
        0,
        0,
        outputs,
        w,
    );
}

#[allow(clippy::too_many_arguments)]
fn solve_unprepared_atmosphere_chunk(
    geometry: &Geometry,
    mode: SourceMode,
    atmosphere: &AtmosphereBatch,
    views: &[View],
    cotangent: &[f64],
    input_base: usize,
    output_base: usize,
    outputs: &mut AtmosphereOutputTile<'_>,
    w: &mut ExplicitWorkspace,
) {
    let (d_albedo, d_thermal_surface) = run_atmosphere_chunk(
        geometry,
        mode,
        atmosphere,
        views,
        cotangent,
        input_base,
        output_base,
        &mut outputs.radiance,
        w,
    );
    map_chunk_to_atmosphere(
        geometry,
        mode,
        atmosphere,
        input_base,
        output_base,
        d_albedo,
        d_thermal_surface,
        0,
        0,
        outputs,
        w,
    );
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn level_optical_adjoint(
    extinction: Wide,
    ssa: Wide,
    b1: Wide,
    layer_ssa: Wide,
    layer_b1: Wide,
    common_od: Wide,
    common_ssa: Wide,
    common_b1: Wide,
) -> (Wide, Wide, Wide) {
    let delta_b1 = b1 - layer_b1;
    (
        common_od + common_ssa * (ssa - layer_ssa) + common_b1 * ssa * delta_b1,
        common_ssa * extinction + common_b1 * extinction * delta_b1,
        common_b1 * ssa * extinction,
    )
}

#[allow(clippy::too_many_arguments)]
fn map_chunk_to_atmosphere(
    geometry: &Geometry,
    mode: SourceMode,
    atmosphere: &AtmosphereBatch,
    input_base: usize,
    output_base: usize,
    d_albedo: Wide,
    d_thermal_surface: Wide,
    level_row_offset: usize,
    surface_row: usize,
    adjoints: &mut AtmosphereOutputTile<'_>,
    w: &mut ExplicitWorkspace,
) {
    let n = geometry.num_layers();
    let input_stride = atmosphere.num_wavelengths;

    if mode == SourceMode::Solar {
        w.attenuation.fill(Wide::splat(0.0));
        for boundary in 1..=n {
            w.attenuation[boundary] = w.d_transmission[boundary] * w.transmission[boundary];
        }
        for layer in 0..n {
            let inv_od = 1.0 / w.od[layer];
            let weight = w.d_secant[layer] * inv_od;
            w.attenuation[layer] += weight;
            w.attenuation[layer + 1] -= weight;
            w.d_od[layer] -= w.d_secant[layer] * w.secant[layer] * inv_od;
        }
        for layer in 0..n {
            let mut contribution = Wide::splat(0.0);
            for boundary in 1..=n {
                contribution +=
                    w.attenuation[boundary] * geometry.chapman_factors[(boundary - 1) * n + layer];
            }
            w.d_od[layer] -= contribution;
        }
    }

    let mut top_extinction = Wide::load(&atmosphere.extinction, input_base, input_stride);
    let mut top_ssa = Wide::load(&atmosphere.single_scatter_albedo, input_base, input_stride);
    let mut top_b1 = Wide::load(&atmosphere.first_legendre, input_base, input_stride);
    let mut top_emission = if mode == SourceMode::Thermal {
        Wide::load(
            atmosphere.emission.as_ref().unwrap(),
            input_base,
            input_stride,
        )
    } else {
        Wide::splat(0.0)
    };
    for layer in 0..n {
        let dz = geometry.layer_thickness[layer];
        let od_density = w.od[layer] / dz;
        let scattering_density = w.ssa[layer] * od_density;
        let inv_od_density = 1.0 / od_density;
        let inv_scattering_density = 1.0 / scattering_density;
        let bottom_extinction =
            Wide::load_row(&atmosphere.extinction, layer + 1, input_base, input_stride);
        let bottom_ssa = Wide::load_row(
            &atmosphere.single_scatter_albedo,
            layer + 1,
            input_base,
            input_stride,
        );
        let bottom_b1 = Wide::load_row(
            &atmosphere.first_legendre,
            layer + 1,
            input_base,
            input_stride,
        );
        let common_od = 0.5 * w.d_od[layer] * dz;
        let common_ssa = 0.5 * w.d_ssa[layer] * inv_od_density;
        let common_b1 = 0.5 * w.d_b1[layer] * inv_scattering_density;
        let (d_extinction, d_ssa, d_b1) = level_optical_adjoint(
            top_extinction,
            top_ssa,
            top_b1,
            w.ssa[layer],
            w.b1[layer],
            common_od,
            common_ssa,
            common_b1,
        );
        w.d_level_extinction[layer] += d_extinction;
        w.d_level_ssa[layer] += d_ssa;
        w.d_level_b1[layer] += d_b1;
        let (d_extinction, d_ssa, d_b1) = level_optical_adjoint(
            bottom_extinction,
            bottom_ssa,
            bottom_b1,
            w.ssa[layer],
            w.b1[layer],
            common_od,
            common_ssa,
            common_b1,
        );
        w.d_level_extinction[layer + 1] += d_extinction;
        w.d_level_ssa[layer + 1] += d_ssa;
        w.d_level_b1[layer + 1] += d_b1;

        if mode == SourceMode::Thermal {
            let emission = atmosphere.emission.as_ref().unwrap();
            let bottom = Wide::load_row(emission, layer + 1, input_base, input_stride);
            let d_b0 = w.d_thermal_b0[layer];
            let d_b1 = w.d_thermal_b1[layer];
            w.d_level_emission[layer] += d_b0 + d_b1 / (dz * top_emission);
            w.d_level_emission[layer + 1] -= d_b1 / (dz * bottom);
            // The physical layer thickness is geometry, not an atmospheric
            // variable.  Therefore thermal_b1 has no additional extinction
            // contribution when mapping back to level quantities.
            top_emission = bottom;
        }
        top_extinction = bottom_extinction;
        top_ssa = bottom_ssa;
        top_b1 = bottom_b1;
    }

    for level in 0..=n {
        adjoints.extinction.write(
            level_row_offset + level,
            output_base,
            w.d_level_extinction[level],
        );
        adjoints.single_scatter_albedo.write(
            level_row_offset + level,
            output_base,
            w.d_level_ssa[level],
        );
        adjoints
            .first_legendre
            .write(level_row_offset + level, output_base, w.d_level_b1[level]);
        if mode == SourceMode::Thermal {
            adjoints.emission.as_mut().unwrap().write(
                level_row_offset + level,
                output_base,
                w.d_level_emission[level],
            );
        }
    }
    adjoints
        .surface_albedo
        .write(surface_row, output_base, d_albedo);
    if mode == SourceMode::Thermal {
        adjoints.surface_emission.as_mut().unwrap().write(
            surface_row,
            output_base,
            d_thermal_surface,
        );
    }
}

/// Map one spherical column's layer adjoints to atmospheric levels without
/// writing them. The caller accumulates these level arrays across SZA columns
/// before emitting the per-view Jacobian.
fn prepare_spherical_column_atmosphere_adjoint<const SOLAR: bool>(
    geometry: &Geometry,
    atmosphere: &AtmosphereBatch,
    input_base: usize,
    w: &mut ExplicitWorkspace,
) {
    let n = geometry.num_layers();
    let input_stride = atmosphere.num_wavelengths;

    if SOLAR {
        w.attenuation.fill(Wide::splat(0.0));
        for boundary in 1..=n {
            w.attenuation[boundary] = w.d_transmission[boundary] * w.transmission[boundary];
        }
        for layer in 0..n {
            let inv_od = 1.0 / w.od[layer];
            let weight = w.d_secant[layer] * inv_od;
            w.attenuation[layer] += weight;
            w.attenuation[layer + 1] -= weight;
            w.d_od[layer] -= w.d_secant[layer] * w.secant[layer] * inv_od;
        }
        for layer in 0..n {
            let mut contribution = Wide::splat(0.0);
            for boundary in 1..=n {
                contribution +=
                    w.attenuation[boundary] * geometry.chapman_factors[(boundary - 1) * n + layer];
            }
            w.d_od[layer] -= contribution;
        }
    }

    let mut top_extinction = Wide::load(&atmosphere.extinction, input_base, input_stride);
    let mut top_ssa = Wide::load(&atmosphere.single_scatter_albedo, input_base, input_stride);
    let mut top_b1 = Wide::load(&atmosphere.first_legendre, input_base, input_stride);
    let mut top_emission = if SOLAR {
        Wide::splat(0.0)
    } else {
        Wide::load(
            atmosphere.emission.as_ref().unwrap(),
            input_base,
            input_stride,
        )
    };
    for layer in 0..n {
        let dz = geometry.layer_thickness[layer];
        let od_density = w.od[layer] / dz;
        let scattering_density = w.ssa[layer] * od_density;
        let bottom_extinction =
            Wide::load_row(&atmosphere.extinction, layer + 1, input_base, input_stride);
        let bottom_ssa = Wide::load_row(
            &atmosphere.single_scatter_albedo,
            layer + 1,
            input_base,
            input_stride,
        );
        let bottom_b1 = Wide::load_row(
            &atmosphere.first_legendre,
            layer + 1,
            input_base,
            input_stride,
        );
        let common_od = 0.5 * w.d_od[layer] * dz;
        let common_ssa = 0.5 * w.d_ssa[layer] / od_density;
        let common_b1 = 0.5 * w.d_b1[layer] / scattering_density;
        let (d_extinction, d_ssa, d_b1) = level_optical_adjoint(
            top_extinction,
            top_ssa,
            top_b1,
            w.ssa[layer],
            w.b1[layer],
            common_od,
            common_ssa,
            common_b1,
        );
        w.d_level_extinction[layer] += d_extinction;
        w.d_level_ssa[layer] += d_ssa;
        w.d_level_b1[layer] += d_b1;
        let (d_extinction, d_ssa, d_b1) = level_optical_adjoint(
            bottom_extinction,
            bottom_ssa,
            bottom_b1,
            w.ssa[layer],
            w.b1[layer],
            common_od,
            common_ssa,
            common_b1,
        );
        w.d_level_extinction[layer + 1] += d_extinction;
        w.d_level_ssa[layer + 1] += d_ssa;
        w.d_level_b1[layer + 1] += d_b1;

        if !SOLAR {
            let emission = atmosphere.emission.as_ref().unwrap();
            let bottom = Wide::load_row(emission, layer + 1, input_base, input_stride);
            let d_b0 = w.d_thermal_b0[layer];
            let d_b1 = w.d_thermal_b1[layer];
            w.d_level_emission[layer] += d_b0 + d_b1 / (dz * top_emission);
            w.d_level_emission[layer + 1] -= d_b1 / (dz * bottom);
            top_emission = bottom;
        }
        top_extinction = bottom_extinction;
        top_ssa = bottom_ssa;
        top_b1 = bottom_b1;
    }
}

fn forward_layers<const SOLAR: bool>(
    geometry: &Geometry,
    albedo: Wide,
    thermal_surface: Wide,
    w: &mut ExplicitWorkspace,
) {
    let n = w.od.len();
    let mu = geometry.quadrature_cosine;
    let naz = if SOLAR { 2 } else { 1 };
    let sin_product = ((1.0 - mu * mu) * (1.0 - geometry.solar_cosine.powi(2))).sqrt();
    for layer in 0..n {
        let exponential = if SOLAR {
            let exponent = -w.secant[layer] * w.od[layer];
            Wide::transmission_ratio(w.transmission[layer], w.transmission[layer + 1], exponent)
        } else {
            (-w.thermal_b1[layer] * w.od[layer]).exp()
        };
        for az in 0..naz {
            w.particular[az].exponential[layer] = exponential;
        }
    }
    for az in 0..naz {
        let h = &mut w.homogeneous[az];
        let p = &mut w.particular[az];
        for layer in 0..n {
            let (d, s) = if az == 0 {
                (
                    w.ssa[layer] * w.b1[layer] * mu - 1.0 / mu,
                    (w.ssa[layer] - 1.0) / mu,
                )
            } else {
                (
                    Wide::splat(-1.0 / mu),
                    (w.ssa[layer] * w.b1[layer] * (1.0 - mu * mu) - 2.0) / (2.0 * mu),
                )
            };
            let k = (s * d).sqrt();
            let s_over_k = s / k;
            let xp = 0.5 * (1.0 - s_over_k);
            let xm = 0.5 * (1.0 + s_over_k);
            let omega = (-k * w.od[layer]).exp();
            let norm = mu * (xp * xp - xm * xm);
            h.d[layer] = d;
            h.s[layer] = s;
            h.k[layer] = k;
            h.xp[layer] = xp;
            h.xm[layer] = xm;
            h.omega[layer] = omega;
            h.norm[layer] = norm;
            let inv_norm = 1.0 / norm;

            if SOLAR {
                let (qp, qm) = if az == 0 {
                    (
                        w.ssa[layer] * (1.0 + w.b1[layer] * geometry.solar_cosine * mu) / FOUR_PI,
                        w.ssa[layer] * (1.0 - w.b1[layer] * geometry.solar_cosine * mu) / FOUR_PI,
                    )
                } else {
                    let q = w.ssa[layer] * w.b1[layer] * sin_product / FOUR_PI;
                    (q, q)
                };
                let ap = (qp * xp + qm * xm) * inv_norm;
                let am = (qm * xp + qp * xm) * inv_norm;
                let exponential = p.exponential[layer];
                let cp = w.transmission[layer] * (omega - exponential) / (w.secant[layer] - k);
                let cm =
                    w.transmission[layer] * (1.0 - omega * exponential) / (w.secant[layer] + k);
                p.qp[layer] = qp;
                p.qm[layer] = qm;
                p.ap[layer] = ap;
                p.am[layer] = am;
                p.cp[layer] = cp;
                p.cm[layer] = cm;
                p.gpt[layer] = am * cm * xm;
                p.gpb[layer] = ap * cp * xp;
                p.gmt[layer] = am * cm * xp;
                p.gmb[layer] = ap * cp * xm;
            } else {
                let at = mu * (1.0 - w.ssa[layer]) * (xp + xm) * inv_norm;
                let exponential = p.exponential[layer];
                let cp = w.thermal_b0[layer] * (omega - exponential) / (w.thermal_b1[layer] - k);
                let cm =
                    w.thermal_b0[layer] * (1.0 - omega * exponential) / (w.thermal_b1[layer] + k);
                p.at[layer] = at;
                p.cp[layer] = cp;
                p.cm[layer] = cm;
                p.gpt[layer] = at * cm * xm;
                p.gpb[layer] = at * cp * xp;
                p.gmt[layer] = at * cm * xp;
                p.gmb[layer] = at * cp * xm;
            }
        }
    }
    for az in 0..naz {
        build_and_solve_bvp::<SOLAR>(geometry, az, albedo, thermal_surface, w);
    }
}

fn build_and_solve_bvp<const SOLAR: bool>(
    geometry: &Geometry,
    az: usize,
    albedo: Wide,
    thermal_surface: Wide,
    w: &mut ExplicitWorkspace,
) {
    let n = w.od.len();
    let size = 2 * n;
    let mu = geometry.quadrature_cosine;
    let h = &w.homogeneous[az];
    let p = &w.particular[az];
    let q = &mut w.bvp[az];
    q.e[0] = Wide::splat(0.0);
    q.c[0] = Wide::splat(0.0);
    q.b[0] = Wide::splat(0.0);
    q.rhs[0] = -p.gpt[0];
    for layer in 0..n - 1 {
        q.rhs[2 * layer + 1] = p.gmt[layer + 1] - p.gmb[layer];
        q.rhs[2 * layer + 2] = p.gpt[layer + 1] - p.gpb[layer];
    }
    let last = size - 1;
    let delta = if az == 0 { 1.0 } else { 0.0 };
    q.rhs[last] = if SOLAR {
        delta * geometry.solar_cosine * albedo / std::f64::consts::PI * w.transmission[n]
            - (p.gmb[n - 1] - 2.0 * delta * mu * albedo * p.gpb[n - 1])
    } else {
        thermal_surface
    };
    q.d[0] = h.xp[0];
    q.a[0] = h.xm[0] * h.omega[0];
    for layer in 0..n - 1 {
        let row = 2 * layer;
        q.e[row + 1] = Wide::splat(0.0);
        q.c[row + 1] = h.xm[layer] * h.omega[layer];
        q.d[row + 1] = h.xp[layer];
        q.a[row + 1] = -h.xm[layer + 1];
        q.b[row + 1] = -h.xp[layer + 1] * h.omega[layer + 1];
        q.e[row + 2] = h.xp[layer] * h.omega[layer];
        q.c[row + 2] = h.xm[layer];
        q.d[row + 2] = -h.xp[layer + 1];
        q.a[row + 2] = -h.xm[layer + 1] * h.omega[layer + 1];
        q.b[row + 2] = Wide::splat(0.0);
    }
    q.e[last] = Wide::splat(0.0);
    q.c[last] = (h.xm[n - 1] - 2.0 * mu * albedo * delta * h.xp[n - 1]) * h.omega[n - 1];
    q.d[last] = h.xp[n - 1] - 2.0 * mu * albedo * delta * h.xm[n - 1];
    q.a[last] = Wide::splat(0.0);
    q.b[last] = Wide::splat(0.0);
    pentadiagonal_solve(q);
}

fn pentadiagonal_solve(q: &mut Bvp) {
    let n = q.d.len();
    q.inv_mu[0] = 1.0 / q.d[0];
    q.alpha[0] = q.a[0] * q.inv_mu[0];
    q.beta[0] = q.b[0] * q.inv_mu[0];
    q.z[0] = q.rhs[0] * q.inv_mu[0];
    if n > 1 {
        q.gamma[1] = q.c[1];
        q.inv_mu[1] = 1.0 / (q.d[1] - q.alpha[0] * q.gamma[1]);
        q.alpha[1] = (q.a[1] - q.beta[0] * q.gamma[1]) * q.inv_mu[1];
        q.beta[1] = q.b[1] * q.inv_mu[1];
        q.z[1] = (q.rhs[1] - q.z[0] * q.gamma[1]) * q.inv_mu[1];
    }
    for i in 2..n {
        q.gamma[i] = q.c[i] - q.alpha[i - 2] * q.e[i];
        q.inv_mu[i] = 1.0 / (q.d[i] - q.beta[i - 2] * q.e[i] - q.alpha[i - 1] * q.gamma[i]);
        if i + 1 < n {
            q.alpha[i] = (q.a[i] - q.beta[i - 1] * q.gamma[i]) * q.inv_mu[i];
        }
        if i + 2 < n {
            q.beta[i] = q.b[i] * q.inv_mu[i];
        }
        q.z[i] = (q.rhs[i] - q.z[i - 2] * q.e[i] - q.z[i - 1] * q.gamma[i]) * q.inv_mu[i];
    }
    q.rhs[n - 1] = q.z[n - 1];
    if n > 1 {
        q.rhs[n - 2] = q.z[n - 2] - q.alpha[n - 2] * q.rhs[n - 1];
    }
    for i in (0..n.saturating_sub(2)).rev() {
        q.rhs[i] = q.z[i] - q.alpha[i] * q.rhs[i + 1] - q.beta[i] * q.rhs[i + 2];
    }
}

fn pentadiagonal_transpose_solve(q: &Bvp, rhs: &mut [Wide]) {
    let n = rhs.len();

    // A = L U in pentadiagonal_solve. First solve U^T z = rhs,
    // where U has unit diagonal and alpha/beta superdiagonals.
    if n > 1 {
        rhs[1] -= q.alpha[0] * rhs[0];
    }
    for i in 2..n {
        rhs[i] -= q.alpha[i - 1] * rhs[i - 1] + q.beta[i - 2] * rhs[i - 2];
    }

    // Then solve L^T x = z using the retained forward factors.
    rhs[n - 1] = rhs[n - 1] * q.inv_mu[n - 1];
    if n > 1 {
        rhs[n - 2] = (rhs[n - 2] - q.gamma[n - 1] * rhs[n - 1]) * q.inv_mu[n - 2];
    }
    for i in (0..n.saturating_sub(2)).rev() {
        rhs[i] = (rhs[i] - q.gamma[i + 1] * rhs[i + 1] - q.e[i + 2] * rhs[i + 2]) * q.inv_mu[i];
    }
}

fn local_source<const SOLAR: bool>(
    geometry: &Geometry,
    layer: usize,
    fraction_from_top: f64,
    propagation_cosine: f64,
    relative_azimuth: f64,
    w: &ExplicitWorkspace,
) -> Wide {
    let mu = geometry.quadrature_cosine;
    let naz = if SOLAR { 2 } else { 1 };
    let phase_mu = propagation_cosine * mu;
    let phase_sine = 0.25
        * ((1.0 - propagation_cosine.powi(2)) * (1.0 - mu * mu))
            .max(0.0)
            .sqrt();
    let azimuth_weight = [1.0, relative_azimuth.cos()];
    let fraction_from_bottom = 1.0 - fraction_from_top;
    let mut source = Wide::splat(0.0);

    for (az, &azimuth_weight) in azimuth_weight.iter().enumerate().take(naz) {
        let h = &w.homogeneous[az];
        let p = &w.particular[az];
        let (lp, lm) = lpsum::<SOLAR>(az, phase_mu, phase_sine, w.ssa[layer], w.b1[layer]);
        let yp = lp * h.xp[layer] + lm * h.xm[layer];
        let ym = lp * h.xm[layer] + lm * h.xp[layer];
        let top_exponential = (-h.k[layer] * w.od[layer] * fraction_from_top).exp();
        let bottom_exponential = (-h.k[layer] * w.od[layer] * fraction_from_bottom).exp();
        let mut value = w.bvp[az].rhs[2 * layer] * yp * top_exponential
            + w.bvp[az].rhs[2 * layer + 1] * ym * bottom_exponential;

        if SOLAR {
            let beam = (-w.secant[layer] * w.od[layer] * fraction_from_top).exp();
            let plus =
                w.transmission[layer] * (top_exponential - beam) / (w.secant[layer] - h.k[layer]);
            let minus = w.transmission[layer] * (beam - p.exponential[layer] * bottom_exponential)
                / (w.secant[layer] + h.k[layer]);
            value += p.ap[layer] * yp * plus + p.am[layer] * ym * minus;
        } else {
            let thermal = (-w.thermal_b1[layer] * w.od[layer] * fraction_from_top).exp();
            let plus = w.thermal_b0[layer] * (top_exponential - thermal)
                / (w.thermal_b1[layer] - h.k[layer]);
            let minus = w.thermal_b0[layer] * (thermal - p.exponential[layer] * bottom_exponential)
                / (w.thermal_b1[layer] + h.k[layer]);
            value += p.at[layer] * (yp * plus + ym * minus);
        }
        source += azimuth_weight * value;
    }

    if !SOLAR {
        let thermal = (-w.thermal_b1[layer] * w.od[layer] * fraction_from_top).exp();
        source += w.thermal_b0[layer] * thermal * (1.0 - w.ssa[layer]);
    }
    source
}

#[allow(clippy::too_many_arguments)]
fn reverse_local_source<const SOLAR: bool>(
    geometry: &Geometry,
    layer: usize,
    fraction_from_top: f64,
    propagation_cosine: f64,
    relative_azimuth: f64,
    seed: Wide,
    w: &mut ExplicitWorkspace,
) {
    let mu = geometry.quadrature_cosine;
    let naz = if SOLAR { 2 } else { 1 };
    let phase_mu = propagation_cosine * mu;
    let phase_sine = 0.25
        * ((1.0 - propagation_cosine.powi(2)) * (1.0 - mu * mu))
            .max(0.0)
            .sqrt();
    let azimuth_weight = [1.0, relative_azimuth.cos()];
    let fraction_from_bottom = 1.0 - fraction_from_top;

    for (az, &azimuth_weight) in azimuth_weight.iter().enumerate().take(naz) {
        let h = &w.homogeneous[az];
        let p = &w.particular[az];
        let (lp, lm) = lpsum::<SOLAR>(az, phase_mu, phase_sine, w.ssa[layer], w.b1[layer]);
        let yp = lp * h.xp[layer] + lm * h.xm[layer];
        let ym = lp * h.xm[layer] + lm * h.xp[layer];
        let top_exponential = (-h.k[layer] * w.od[layer] * fraction_from_top).exp();
        let bottom_exponential = (-h.k[layer] * w.od[layer] * fraction_from_bottom).exp();
        let source_scale = seed * azimuth_weight;
        let top_solution = w.bvp[az].rhs[2 * layer];
        let bottom_solution = w.bvp[az].rhs[2 * layer + 1];
        let mut d_yp = source_scale * top_solution * top_exponential;
        let mut d_ym = source_scale * bottom_solution * bottom_exponential;
        let mut d_top_exponential = source_scale * top_solution * yp;
        let mut d_bottom_exponential = source_scale * bottom_solution * ym;
        w.d_solution[az][2 * layer] += source_scale * yp * top_exponential;
        w.d_solution[az][2 * layer + 1] += source_scale * ym * bottom_exponential;
        let mut d_k = Wide::splat(0.0);

        if SOLAR {
            let beam = (-w.secant[layer] * w.od[layer] * fraction_from_top).exp();
            let plus_denominator = w.secant[layer] - h.k[layer];
            let minus_denominator = w.secant[layer] + h.k[layer];
            let inv_plus_denominator = 1.0 / plus_denominator;
            let inv_minus_denominator = 1.0 / minus_denominator;
            let plus_numerator = top_exponential - beam;
            let minus_numerator = beam - p.exponential[layer] * bottom_exponential;
            let plus = w.transmission[layer] * plus_numerator * inv_plus_denominator;
            let minus = w.transmission[layer] * minus_numerator * inv_minus_denominator;

            w.d_particular[az].ap[layer] += source_scale * yp * plus;
            d_yp += source_scale * p.ap[layer] * plus;
            let d_plus = source_scale * p.ap[layer] * yp;
            w.d_particular[az].am[layer] += source_scale * ym * minus;
            d_ym += source_scale * p.am[layer] * minus;
            let d_minus = source_scale * p.am[layer] * ym;

            w.d_transmission[layer] += d_plus * plus_numerator * inv_plus_denominator;
            let d_plus_numerator = d_plus * w.transmission[layer] * inv_plus_denominator;
            let d_plus_denominator = -d_plus * plus * inv_plus_denominator;
            d_top_exponential += d_plus_numerator;
            let mut d_beam = -d_plus_numerator;
            w.d_secant[layer] += d_plus_denominator;
            d_k -= d_plus_denominator;

            w.d_transmission[layer] += d_minus * minus_numerator * inv_minus_denominator;
            let d_minus_numerator = d_minus * w.transmission[layer] * inv_minus_denominator;
            let d_minus_denominator = -d_minus * minus * inv_minus_denominator;
            d_beam += d_minus_numerator;
            let d_layer_exponential = -d_minus_numerator * bottom_exponential;
            d_bottom_exponential -= d_minus_numerator * p.exponential[layer];
            w.d_secant[layer] += d_minus_denominator;
            d_k += d_minus_denominator;

            w.d_secant[layer] -= d_beam * beam * w.od[layer] * fraction_from_top;
            w.d_od[layer] -= d_beam * beam * w.secant[layer] * fraction_from_top;
            w.d_secant[layer] -= d_layer_exponential * p.exponential[layer] * w.od[layer];
            w.d_od[layer] -= d_layer_exponential * p.exponential[layer] * w.secant[layer];
        } else {
            let thermal = (-w.thermal_b1[layer] * w.od[layer] * fraction_from_top).exp();
            let plus_denominator = w.thermal_b1[layer] - h.k[layer];
            let minus_denominator = w.thermal_b1[layer] + h.k[layer];
            let inv_plus_denominator = 1.0 / plus_denominator;
            let inv_minus_denominator = 1.0 / minus_denominator;
            let plus_numerator = top_exponential - thermal;
            let minus_numerator = thermal - p.exponential[layer] * bottom_exponential;
            let plus = w.thermal_b0[layer] * plus_numerator * inv_plus_denominator;
            let minus = w.thermal_b0[layer] * minus_numerator * inv_minus_denominator;

            w.d_particular[az].at[layer] += source_scale * (yp * plus + ym * minus);
            d_yp += source_scale * p.at[layer] * plus;
            d_ym += source_scale * p.at[layer] * minus;
            let d_plus = source_scale * p.at[layer] * yp;
            let d_minus = source_scale * p.at[layer] * ym;

            w.d_thermal_b0[layer] += d_plus * plus_numerator * inv_plus_denominator;
            let d_plus_numerator = d_plus * w.thermal_b0[layer] * inv_plus_denominator;
            let d_plus_denominator = -d_plus * plus * inv_plus_denominator;
            d_top_exponential += d_plus_numerator;
            let mut d_thermal = -d_plus_numerator;
            w.d_thermal_b1[layer] += d_plus_denominator;
            d_k -= d_plus_denominator;

            w.d_thermal_b0[layer] += d_minus * minus_numerator * inv_minus_denominator;
            let d_minus_numerator = d_minus * w.thermal_b0[layer] * inv_minus_denominator;
            let d_minus_denominator = -d_minus * minus * inv_minus_denominator;
            d_thermal += d_minus_numerator;
            let d_layer_exponential = -d_minus_numerator * bottom_exponential;
            d_bottom_exponential -= d_minus_numerator * p.exponential[layer];
            w.d_thermal_b1[layer] += d_minus_denominator;
            d_k += d_minus_denominator;

            w.d_thermal_b1[layer] -= d_thermal * thermal * w.od[layer] * fraction_from_top;
            w.d_od[layer] -= d_thermal * thermal * w.thermal_b1[layer] * fraction_from_top;
            w.d_thermal_b1[layer] -= d_layer_exponential * p.exponential[layer] * w.od[layer];
            w.d_od[layer] -= d_layer_exponential * p.exponential[layer] * w.thermal_b1[layer];
        }

        let d_lp = d_yp * h.xp[layer] + d_ym * h.xm[layer];
        let d_lm = d_yp * h.xm[layer] + d_ym * h.xp[layer];
        w.d_homogeneous[az].xp[layer] += d_yp * lp + d_ym * lm;
        w.d_homogeneous[az].xm[layer] += d_yp * lm + d_ym * lp;
        if SOLAR {
            if az == 0 {
                w.d_ssa[layer] += d_lp * 0.5 * (1.0 - w.b1[layer] * phase_mu)
                    + d_lm * 0.5 * (1.0 + w.b1[layer] * phase_mu);
                w.d_b1[layer] += d_lp * (-0.5 * w.ssa[layer] * phase_mu)
                    + d_lm * (0.5 * w.ssa[layer] * phase_mu);
            } else {
                w.d_ssa[layer] += (d_lp + d_lm) * w.b1[layer] * phase_sine;
                w.d_b1[layer] += (d_lp + d_lm) * w.ssa[layer] * phase_sine;
            }
        }

        d_k -= d_top_exponential * top_exponential * w.od[layer] * fraction_from_top;
        w.d_od[layer] -= d_top_exponential * top_exponential * h.k[layer] * fraction_from_top;
        d_k -= d_bottom_exponential * bottom_exponential * w.od[layer] * fraction_from_bottom;
        w.d_od[layer] -=
            d_bottom_exponential * bottom_exponential * h.k[layer] * fraction_from_bottom;
        w.d_homogeneous[az].k[layer] += d_k;
    }

    if !SOLAR {
        let thermal = (-w.thermal_b1[layer] * w.od[layer] * fraction_from_top).exp();
        w.d_thermal_b0[layer] += seed * thermal * (1.0 - w.ssa[layer]);
        let d_thermal = seed * w.thermal_b0[layer] * (1.0 - w.ssa[layer]);
        w.d_ssa[layer] -= seed * w.thermal_b0[layer] * thermal;
        w.d_thermal_b1[layer] -= d_thermal * thermal * w.od[layer] * fraction_from_top;
        w.d_od[layer] -= d_thermal * thermal * w.thermal_b1[layer] * fraction_from_top;
    }
}

fn surface_source<const SOLAR: bool>(
    geometry: &Geometry,
    albedo: Wide,
    thermal_surface: Wide,
    w: &ExplicitWorkspace,
) -> Wide {
    let last = geometry.num_layers() - 1;
    let base = w.particular[0].gpb[last]
        + w.bvp[0].rhs[2 * last] * w.homogeneous[0].xp[last] * w.homogeneous[0].omega[last]
        + w.bvp[0].rhs[2 * last + 1] * w.homogeneous[0].xm[last];
    base * (2.0 * geometry.quadrature_cosine) * albedo
        + if SOLAR {
            Wide::splat(0.0)
        } else {
            thermal_surface
        }
}

fn reverse_surface_source<const SOLAR: bool>(
    geometry: &Geometry,
    albedo: Wide,
    seed: Wide,
    w: &mut ExplicitWorkspace,
) -> (Wide, Wide) {
    let last = geometry.num_layers() - 1;
    let base = w.particular[0].gpb[last]
        + w.bvp[0].rhs[2 * last] * w.homogeneous[0].xp[last] * w.homogeneous[0].omega[last]
        + w.bvp[0].rhs[2 * last + 1] * w.homogeneous[0].xm[last];
    let scale = 2.0 * geometry.quadrature_cosine;
    let d_base = seed * scale * albedo;
    w.d_particular[0].gpb[last] += d_base;
    w.d_solution[0][2 * last] += d_base * w.homogeneous[0].xp[last] * w.homogeneous[0].omega[last];
    w.d_homogeneous[0].xp[last] += d_base * w.bvp[0].rhs[2 * last] * w.homogeneous[0].omega[last];
    w.d_homogeneous[0].omega[last] += d_base * w.bvp[0].rhs[2 * last] * w.homogeneous[0].xp[last];
    w.d_solution[0][2 * last + 1] += d_base * w.homogeneous[0].xm[last];
    w.d_homogeneous[0].xm[last] += d_base * w.bvp[0].rhs[2 * last + 1];
    (
        seed * scale * base,
        if SOLAR { Wide::splat(0.0) } else { seed },
    )
}

#[allow(clippy::too_many_arguments)]
fn reverse_views<const SOLAR: bool>(
    geometry: &Geometry,
    albedo: Wide,
    thermal_surface: Wide,
    views: &[View],
    seed: ViewSeed<'_>,
    input_base: usize,
    input_stride: usize,
    output_base: usize,
    output_view_offset: usize,
    radiance: &mut OutputRows<'_>,
    w: &mut ExplicitWorkspace,
) -> (Wide, Wide) {
    let n = w.od.len();
    let mu = geometry.quadrature_cosine;
    let naz = if SOLAR { 2 } else { 1 };
    let mut d_albedo = Wide::splat(0.0);
    let mut d_thermal_surface = Wide::splat(0.0);
    for (view_index, view) in views.iter().enumerate() {
        let vmu = view.cosine;
        let inv_vmu = 1.0 / vmu;
        let phase_mu = vmu * mu;
        let phase_sine = 0.25 * ((1.0 - vmu * vmu) * (1.0 - mu * mu)).sqrt();
        let azimuth_weight = [1.0, view.relative_azimuth.cos()];
        w.source.fill(Wide::splat(0.0));
        w.attenuation[0] = Wide::splat(1.0);
        let mut integrated = Wide::splat(0.0);
        for layer in 0..n {
            let beam = (-w.od[layer] * inv_vmu).exp();
            w.beam[layer] = beam;
            let exponential = w.particular[0].exponential[layer];
            let denominator = if SOLAR {
                1.0 + w.secant[layer] * vmu
            } else {
                1.0 + w.thermal_b1[layer] * vmu
            };
            let inv_denominator = 1.0 / denominator;
            let source_integral = (1.0 - exponential * beam) * inv_denominator;
            for (az, &azimuth_weight) in azimuth_weight.iter().enumerate().take(naz) {
                let h = &w.homogeneous[az];
                let p = &w.particular[az];
                let (lp, lm) = lpsum::<SOLAR>(az, phase_mu, phase_sine, w.ssa[layer], w.b1[layer]);
                let yp = lp * h.xp[layer] + lm * h.xm[layer];
                let ym = lp * h.xm[layer] + lm * h.xp[layer];
                let hm = (h.omega[layer] - beam) / (1.0 - h.k[layer] * vmu);
                let hp = (1.0 - h.omega[layer] * beam) / (1.0 + h.k[layer] * vmu);
                let v = if SOLAR {
                    let dp = w.transmission[layer] * (source_integral - exponential * hm)
                        / (w.secant[layer] + h.k[layer]);
                    let dm = w.transmission[layer] * (hp - source_integral)
                        / (w.secant[layer] - h.k[layer]);
                    p.ap[layer] * yp * dm + p.am[layer] * ym * dp
                } else {
                    let dp = w.thermal_b0[layer] * (source_integral - exponential * hm)
                        / (w.thermal_b1[layer] + h.k[layer]);
                    let dm = w.thermal_b0[layer] * (hp - source_integral)
                        / (w.thermal_b1[layer] - h.k[layer]);
                    p.at[layer] * (yp * dm + ym * dp)
                };
                w.source[layer] += azimuth_weight
                    * (w.bvp[az].rhs[2 * layer] * yp * hp
                        + w.bvp[az].rhs[2 * layer + 1] * ym * hm
                        + v);
                if !SOLAR {
                    w.source[layer] += w.thermal_b0[layer] * source_integral * (1.0 - w.ssa[layer]);
                }
            }
            integrated += w.source[layer] * w.attenuation[layer];
            w.attenuation[layer + 1] = w.attenuation[layer] * beam;
        }

        let last = n - 1;
        let base_surface = w.particular[0].gpb[last]
            + w.bvp[0].rhs[2 * last] * w.homogeneous[0].xp[last] * w.homogeneous[0].omega[last]
            + w.bvp[0].rhs[2 * last + 1] * w.homogeneous[0].xm[last];
        let surface = base_surface * (2.0 * mu) * albedo;
        let output = integrated + w.attenuation[n] * (surface + thermal_surface);
        let output_view = output_view_offset + view_index;
        radiance.write(output_view, output_base, output);

        let seed = match seed {
            ViewSeed::Cotangent(cotangent) => {
                Wide::load_row(cotangent, output_view, input_base, input_stride)
            }
            ViewSeed::Unit => Wide::splat(1.0),
            ViewSeed::None => continue,
        };
        let d_surface = seed * w.attenuation[n];
        d_thermal_surface += seed * w.attenuation[n];
        let d_base_surface = d_surface * (2.0 * mu) * albedo;
        d_albedo += d_surface * (2.0 * mu) * base_surface;
        w.d_particular[0].gpb[last] += d_base_surface;
        w.d_solution[0][2 * last] +=
            d_base_surface * w.homogeneous[0].xp[last] * w.homogeneous[0].omega[last];
        w.d_homogeneous[0].xp[last] +=
            d_base_surface * w.bvp[0].rhs[2 * last] * w.homogeneous[0].omega[last];
        w.d_homogeneous[0].omega[last] +=
            d_base_surface * w.bvp[0].rhs[2 * last] * w.homogeneous[0].xp[last];
        w.d_solution[0][2 * last + 1] += d_base_surface * w.homogeneous[0].xm[last];
        w.d_homogeneous[0].xm[last] += d_base_surface * w.bvp[0].rhs[2 * last + 1];

        let mut d_attenuation = seed * (surface + thermal_surface);
        for layer in (0..n).rev() {
            let beam = w.beam[layer];
            let mut d_beam = d_attenuation * w.attenuation[layer];
            d_attenuation = d_attenuation * beam + seed * w.source[layer];
            let d_source = seed * w.attenuation[layer];
            let exponential = w.particular[0].exponential[layer];
            let denominator = if SOLAR {
                1.0 + w.secant[layer] * vmu
            } else {
                1.0 + w.thermal_b1[layer] * vmu
            };
            let inv_denominator = 1.0 / denominator;
            let source_integral = (1.0 - exponential * beam) * inv_denominator;
            let mut d_exponential = Wide::splat(0.0);
            let mut d_source_integral = Wide::splat(0.0);
            for (az, &azimuth_weight) in azimuth_weight.iter().enumerate().take(naz) {
                reverse_view_azimuth::<SOLAR>(
                    az,
                    *view,
                    azimuth_weight,
                    phase_mu,
                    phase_sine,
                    layer,
                    beam,
                    exponential,
                    source_integral,
                    d_source,
                    &mut d_beam,
                    &mut d_exponential,
                    &mut d_source_integral,
                    w,
                );
                if !SOLAR {
                    w.d_thermal_b0[layer] += d_source * source_integral * (1.0 - w.ssa[layer]);
                    d_source_integral += d_source * w.thermal_b0[layer] * (1.0 - w.ssa[layer]);
                    w.d_ssa[layer] -= d_source * w.thermal_b0[layer] * source_integral;
                }
            }
            let d_numerator = d_source_integral * inv_denominator;
            let d_denominator = -d_source_integral * source_integral * inv_denominator;
            d_exponential -= d_numerator * beam;
            d_beam -= d_numerator * exponential;
            if SOLAR {
                w.d_secant[layer] += d_denominator * vmu;
                w.d_secant[layer] -= d_exponential * exponential * w.od[layer];
                w.d_od[layer] -= d_exponential * exponential * w.secant[layer];
            } else {
                w.d_thermal_b1[layer] += d_denominator * vmu;
                w.d_thermal_b1[layer] -= d_exponential * exponential * w.od[layer];
                w.d_od[layer] -= d_exponential * exponential * w.thermal_b1[layer];
            }
            w.d_od[layer] -= d_beam * beam * inv_vmu;
        }
    }
    (d_albedo, d_thermal_surface)
}

fn lpsum<const SOLAR: bool>(
    az: usize,
    phase_mu: f64,
    phase_sine: f64,
    ssa: Wide,
    b1: Wide,
) -> (Wide, Wide) {
    if !SOLAR {
        return (Wide::splat(0.0), Wide::splat(0.0));
    }
    if az == 0 {
        (
            0.5 * ssa * (1.0 - b1 * phase_mu),
            0.5 * ssa * (1.0 + b1 * phase_mu),
        )
    } else {
        let value = ssa * b1 * phase_sine;
        (value, value)
    }
}

#[allow(clippy::too_many_arguments)]
fn reverse_view_azimuth<const SOLAR: bool>(
    az: usize,
    view: View,
    azimuth_weight: f64,
    phase_mu: f64,
    phase_sine: f64,
    layer: usize,
    beam: Wide,
    exponential: Wide,
    source_integral: Wide,
    d_source: Wide,
    d_beam: &mut Wide,
    d_exponential: &mut Wide,
    d_source_integral: &mut Wide,
    w: &mut ExplicitWorkspace,
) {
    let vmu = view.cosine;
    let h = &w.homogeneous[az];
    let p = &w.particular[az];
    let (lp, lm) = lpsum::<SOLAR>(az, phase_mu, phase_sine, w.ssa[layer], w.b1[layer]);
    let yp = lp * h.xp[layer] + lm * h.xm[layer];
    let ym = lp * h.xm[layer] + lm * h.xp[layer];
    let hm_denominator = 1.0 - h.k[layer] * vmu;
    let hp_denominator = 1.0 + h.k[layer] * vmu;
    let inv_hm_denominator = 1.0 / hm_denominator;
    let inv_hp_denominator = 1.0 / hp_denominator;
    let hm = (h.omega[layer] - beam) * inv_hm_denominator;
    let hp = (1.0 - h.omega[layer] * beam) * inv_hp_denominator;
    let source_scale = d_source * azimuth_weight;
    let mut d_yp = source_scale * w.bvp[az].rhs[2 * layer] * hp;
    let mut d_ym = source_scale * w.bvp[az].rhs[2 * layer + 1] * hm;
    let mut d_hp = source_scale * w.bvp[az].rhs[2 * layer] * yp;
    let mut d_hm = source_scale * w.bvp[az].rhs[2 * layer + 1] * ym;
    w.d_solution[az][2 * layer] += source_scale * yp * hp;
    w.d_solution[az][2 * layer + 1] += source_scale * ym * hm;

    if SOLAR {
        let dp_denominator = w.secant[layer] + h.k[layer];
        let dm_denominator = w.secant[layer] - h.k[layer];
        let inv_dp_denominator = 1.0 / dp_denominator;
        let inv_dm_denominator = 1.0 / dm_denominator;
        let dp = w.transmission[layer] * (source_integral - exponential * hm) * inv_dp_denominator;
        let dm = w.transmission[layer] * (hp - source_integral) * inv_dm_denominator;
        w.d_particular[az].ap[layer] += source_scale * yp * dm;
        d_yp += source_scale * p.ap[layer] * dm;
        let d_dm = source_scale * p.ap[layer] * yp;
        w.d_particular[az].am[layer] += source_scale * ym * dp;
        d_ym += source_scale * p.am[layer] * dp;
        let d_dp = source_scale * p.am[layer] * ym;

        let dp_numerator = source_integral - exponential * hm;
        w.d_transmission[layer] += d_dp * dp_numerator * inv_dp_denominator;
        let d_dp_numerator = d_dp * w.transmission[layer] * inv_dp_denominator;
        let d_dp_denominator = -d_dp * dp * inv_dp_denominator;
        *d_source_integral += d_dp_numerator;
        *d_exponential -= d_dp_numerator * hm;
        d_hm -= d_dp_numerator * exponential;
        w.d_secant[layer] += d_dp_denominator;
        w.d_homogeneous[az].k[layer] += d_dp_denominator;

        let dm_numerator = hp - source_integral;
        w.d_transmission[layer] += d_dm * dm_numerator * inv_dm_denominator;
        let d_dm_numerator = d_dm * w.transmission[layer] * inv_dm_denominator;
        let d_dm_denominator = -d_dm * dm * inv_dm_denominator;
        d_hp += d_dm_numerator;
        *d_source_integral -= d_dm_numerator;
        w.d_secant[layer] += d_dm_denominator;
        w.d_homogeneous[az].k[layer] -= d_dm_denominator;
    } else {
        let dp_denominator = w.thermal_b1[layer] + h.k[layer];
        let dm_denominator = w.thermal_b1[layer] - h.k[layer];
        let inv_dp_denominator = 1.0 / dp_denominator;
        let inv_dm_denominator = 1.0 / dm_denominator;
        let dp = w.thermal_b0[layer] * (source_integral - exponential * hm) * inv_dp_denominator;
        let dm = w.thermal_b0[layer] * (hp - source_integral) * inv_dm_denominator;
        w.d_particular[az].at[layer] += source_scale * (yp * dm + ym * dp);
        d_yp += source_scale * p.at[layer] * dm;
        d_ym += source_scale * p.at[layer] * dp;
        let d_dm = source_scale * p.at[layer] * yp;
        let d_dp = source_scale * p.at[layer] * ym;

        let dp_numerator = source_integral - exponential * hm;
        w.d_thermal_b0[layer] += d_dp * dp_numerator * inv_dp_denominator;
        let d_dp_numerator = d_dp * w.thermal_b0[layer] * inv_dp_denominator;
        let d_dp_denominator = -d_dp * dp * inv_dp_denominator;
        *d_source_integral += d_dp_numerator;
        *d_exponential -= d_dp_numerator * hm;
        d_hm -= d_dp_numerator * exponential;
        w.d_thermal_b1[layer] += d_dp_denominator;
        w.d_homogeneous[az].k[layer] += d_dp_denominator;

        let dm_numerator = hp - source_integral;
        w.d_thermal_b0[layer] += d_dm * dm_numerator * inv_dm_denominator;
        let d_dm_numerator = d_dm * w.thermal_b0[layer] * inv_dm_denominator;
        let d_dm_denominator = -d_dm * dm * inv_dm_denominator;
        d_hp += d_dm_numerator;
        *d_source_integral -= d_dm_numerator;
        w.d_thermal_b1[layer] += d_dm_denominator;
        w.d_homogeneous[az].k[layer] -= d_dm_denominator;
    }

    let d_hp_numerator = d_hp * inv_hp_denominator;
    let d_hp_denominator = -d_hp * hp * inv_hp_denominator;
    w.d_homogeneous[az].omega[layer] -= d_hp_numerator * beam;
    *d_beam -= d_hp_numerator * h.omega[layer];
    w.d_homogeneous[az].k[layer] += d_hp_denominator * vmu;

    let d_hm_numerator = d_hm * inv_hm_denominator;
    let d_hm_denominator = -d_hm * hm * inv_hm_denominator;
    w.d_homogeneous[az].omega[layer] += d_hm_numerator;
    *d_beam -= d_hm_numerator;
    w.d_homogeneous[az].k[layer] -= d_hm_denominator * vmu;

    let d_lp = d_yp * h.xp[layer] + d_ym * h.xm[layer];
    let d_lm = d_yp * h.xm[layer] + d_ym * h.xp[layer];
    w.d_homogeneous[az].xp[layer] += d_yp * lp + d_ym * lm;
    w.d_homogeneous[az].xm[layer] += d_yp * lm + d_ym * lp;
    if SOLAR {
        if az == 0 {
            w.d_ssa[layer] += d_lp * 0.5 * (1.0 - w.b1[layer] * phase_mu)
                + d_lm * 0.5 * (1.0 + w.b1[layer] * phase_mu);
            w.d_b1[layer] +=
                d_lp * (-0.5 * w.ssa[layer] * phase_mu) + d_lm * (0.5 * w.ssa[layer] * phase_mu);
        } else {
            w.d_ssa[layer] += (d_lp + d_lm) * w.b1[layer] * phase_sine;
            w.d_b1[layer] += (d_lp + d_lm) * w.ssa[layer] * phase_sine;
        }
    }
}

fn reverse_bvp<const SOLAR: bool>(
    geometry: &Geometry,
    albedo: Wide,
    _thermal_surface: Wide,
    w: &mut ExplicitWorkspace,
) -> (Wide, Wide) {
    let n = w.od.len();
    let size = 2 * n;
    let mu = geometry.quadrature_cosine;
    let naz = if SOLAR { 2 } else { 1 };
    let mut d_albedo = Wide::splat(0.0);
    let mut d_thermal_surface = Wide::splat(0.0);
    for az in 0..naz {
        pentadiagonal_transpose_solve(&w.bvp[az], &mut w.d_solution[az]);

        let lambda = &w.d_solution[az];
        let solution = &w.bvp[az].rhs;
        let delta = if az == 0 { 1.0 } else { 0.0 };

        // Reverse the BVP right-hand side.
        w.d_particular[az].gpt[0] -= lambda[0];
        for layer in 0..n - 1 {
            w.d_particular[az].gmt[layer + 1] += lambda[2 * layer + 1];
            w.d_particular[az].gmb[layer] -= lambda[2 * layer + 1];
            w.d_particular[az].gpt[layer + 1] += lambda[2 * layer + 2];
            w.d_particular[az].gpb[layer] -= lambda[2 * layer + 2];
        }
        let last_row = size - 1;
        if SOLAR {
            w.d_transmission[n] +=
                lambda[last_row] * delta * geometry.solar_cosine * albedo / std::f64::consts::PI;
            d_albedo += lambda[last_row]
                * (delta * geometry.solar_cosine / std::f64::consts::PI * w.transmission[n]
                    + 2.0 * delta * mu * w.particular[az].gpb[n - 1]);
            w.d_particular[az].gmb[n - 1] -= lambda[last_row];
            w.d_particular[az].gpb[n - 1] += lambda[last_row] * 2.0 * delta * mu * albedo;
        } else {
            d_thermal_surface += lambda[last_row];
        }

        // Reverse the pentadiagonal matrix using dA = -lambda * solution^T.
        let g_d0 = -lambda[0] * solution[0];
        let g_a0 = -lambda[0] * solution[1];
        w.d_homogeneous[az].xp[0] += g_d0;
        w.d_homogeneous[az].xm[0] += g_a0 * w.homogeneous[az].omega[0];
        w.d_homogeneous[az].omega[0] += g_a0 * w.homogeneous[az].xm[0];

        for layer in 0..n - 1 {
            let row = 2 * layer + 1;
            let g_c = -lambda[row] * solution[row - 1];
            let g_d = -lambda[row] * solution[row];
            let g_a = -lambda[row] * solution[row + 1];
            let g_b = -lambda[row] * solution[row + 2];
            w.d_homogeneous[az].xm[layer] += g_c * w.homogeneous[az].omega[layer];
            w.d_homogeneous[az].omega[layer] += g_c * w.homogeneous[az].xm[layer];
            w.d_homogeneous[az].xp[layer] += g_d;
            w.d_homogeneous[az].xm[layer + 1] -= g_a;
            w.d_homogeneous[az].xp[layer + 1] -= g_b * w.homogeneous[az].omega[layer + 1];
            w.d_homogeneous[az].omega[layer + 1] -= g_b * w.homogeneous[az].xp[layer + 1];

            let row = 2 * layer + 2;
            let g_e = -lambda[row] * solution[row - 2];
            let g_c = -lambda[row] * solution[row - 1];
            let g_d = -lambda[row] * solution[row];
            let g_a = if row + 1 < size {
                -lambda[row] * solution[row + 1]
            } else {
                Wide::splat(0.0)
            };
            w.d_homogeneous[az].xp[layer] += g_e * w.homogeneous[az].omega[layer];
            w.d_homogeneous[az].omega[layer] += g_e * w.homogeneous[az].xp[layer];
            w.d_homogeneous[az].xm[layer] += g_c;
            w.d_homogeneous[az].xp[layer + 1] -= g_d;
            w.d_homogeneous[az].xm[layer + 1] -= g_a * w.homogeneous[az].omega[layer + 1];
            w.d_homogeneous[az].omega[layer + 1] -= g_a * w.homogeneous[az].xm[layer + 1];
        }

        let layer = n - 1;
        let g_c = -lambda[last_row] * solution[last_row - 1];
        let g_d = -lambda[last_row] * solution[last_row];
        let c_base =
            w.homogeneous[az].xm[layer] - 2.0 * mu * albedo * delta * w.homogeneous[az].xp[layer];
        w.d_homogeneous[az].xm[layer] += g_c * w.homogeneous[az].omega[layer];
        w.d_homogeneous[az].xp[layer] -=
            g_c * 2.0 * mu * albedo * delta * w.homogeneous[az].omega[layer];
        w.d_homogeneous[az].omega[layer] += g_c * c_base;
        d_albedo -=
            g_c * 2.0 * mu * delta * w.homogeneous[az].xp[layer] * w.homogeneous[az].omega[layer];
        w.d_homogeneous[az].xp[layer] += g_d;
        w.d_homogeneous[az].xm[layer] -= g_d * 2.0 * mu * albedo * delta;
        d_albedo -= g_d * 2.0 * mu * delta * w.homogeneous[az].xm[layer];
    }
    (d_albedo, d_thermal_surface)
}

fn reverse_layers<const SOLAR: bool>(geometry: &Geometry, w: &mut ExplicitWorkspace) {
    let n = w.od.len();
    let mu = geometry.quadrature_cosine;
    let naz = if SOLAR { 2 } else { 1 };
    let solar_sine = ((1.0 - mu * mu) * (1.0 - geometry.solar_cosine.powi(2))).sqrt();
    for az in 0..naz {
        for layer in 0..n {
            let h = &w.homogeneous[az];
            let p = &w.particular[az];
            let mut d_xp = w.d_homogeneous[az].xp[layer];
            let mut d_xm = w.d_homogeneous[az].xm[layer];
            let mut d_omega = w.d_homogeneous[az].omega[layer];
            let mut d_k = w.d_homogeneous[az].k[layer];
            let mut d_norm = Wide::splat(0.0);
            let mut d_cp = Wide::splat(0.0);
            let mut d_cm = Wide::splat(0.0);

            if SOLAR {
                let mut d_ap = w.d_particular[az].ap[layer];
                let mut d_am = w.d_particular[az].am[layer];
                let d_gpt = w.d_particular[az].gpt[layer];
                let d_gpb = w.d_particular[az].gpb[layer];
                let d_gmt = w.d_particular[az].gmt[layer];
                let d_gmb = w.d_particular[az].gmb[layer];

                d_am += d_gpt * p.cm[layer] * h.xm[layer] + d_gmt * p.cm[layer] * h.xp[layer];
                d_cm += d_gpt * p.am[layer] * h.xm[layer] + d_gmt * p.am[layer] * h.xp[layer];
                d_xm += d_gpt * p.am[layer] * p.cm[layer];
                d_xp += d_gmt * p.am[layer] * p.cm[layer];
                d_ap += d_gpb * p.cp[layer] * h.xp[layer] + d_gmb * p.cp[layer] * h.xm[layer];
                d_cp += d_gpb * p.ap[layer] * h.xp[layer] + d_gmb * p.ap[layer] * h.xm[layer];
                d_xp += d_gpb * p.ap[layer] * p.cp[layer];
                d_xm += d_gmb * p.ap[layer] * p.cp[layer];

                let mut d_qp = Wide::splat(0.0);
                let mut d_qm = Wide::splat(0.0);
                let inv_norm = 1.0 / h.norm[layer];
                let d_ap_numerator = d_ap * inv_norm;
                d_norm -= d_ap * p.ap[layer] * inv_norm;
                d_qp += d_ap_numerator * h.xp[layer];
                d_xp += d_ap_numerator * p.qp[layer];
                d_qm += d_ap_numerator * h.xm[layer];
                d_xm += d_ap_numerator * p.qm[layer];
                let d_am_numerator = d_am * inv_norm;
                d_norm -= d_am * p.am[layer] * inv_norm;
                d_qm += d_am_numerator * h.xp[layer];
                d_xp += d_am_numerator * p.qm[layer];
                d_qp += d_am_numerator * h.xm[layer];
                d_xm += d_am_numerator * p.qp[layer];

                if az == 0 {
                    w.d_ssa[layer] += d_qp * (1.0 + w.b1[layer] * geometry.solar_cosine * mu)
                        / FOUR_PI
                        + d_qm * (1.0 - w.b1[layer] * geometry.solar_cosine * mu) / FOUR_PI;
                    w.d_b1[layer] += d_qp * w.ssa[layer] * geometry.solar_cosine * mu / FOUR_PI
                        - d_qm * w.ssa[layer] * geometry.solar_cosine * mu / FOUR_PI;
                } else {
                    w.d_ssa[layer] += (d_qp + d_qm) * w.b1[layer] * solar_sine / FOUR_PI;
                    w.d_b1[layer] += (d_qp + d_qm) * w.ssa[layer] * solar_sine / FOUR_PI;
                }
            } else {
                let mut d_at = w.d_particular[az].at[layer];
                let d_gpt = w.d_particular[az].gpt[layer];
                let d_gpb = w.d_particular[az].gpb[layer];
                let d_gmt = w.d_particular[az].gmt[layer];
                let d_gmb = w.d_particular[az].gmb[layer];
                d_at += d_gpt * p.cm[layer] * h.xm[layer]
                    + d_gpb * p.cp[layer] * h.xp[layer]
                    + d_gmt * p.cm[layer] * h.xp[layer]
                    + d_gmb * p.cp[layer] * h.xm[layer];
                d_cm += d_gpt * p.at[layer] * h.xm[layer] + d_gmt * p.at[layer] * h.xp[layer];
                d_cp += d_gpb * p.at[layer] * h.xp[layer] + d_gmb * p.at[layer] * h.xm[layer];
                d_xm += d_gpt * p.at[layer] * p.cm[layer] + d_gmb * p.at[layer] * p.cp[layer];
                d_xp += d_gpb * p.at[layer] * p.cp[layer] + d_gmt * p.at[layer] * p.cm[layer];

                let inv_norm = 1.0 / h.norm[layer];
                let d_at_numerator = d_at * inv_norm;
                d_norm -= d_at * p.at[layer] * inv_norm;
                w.d_ssa[layer] -= d_at_numerator * mu * (h.xp[layer] + h.xm[layer]);
                d_xp += d_at_numerator * mu * (1.0 - w.ssa[layer]);
                d_xm += d_at_numerator * mu * (1.0 - w.ssa[layer]);
            }

            let mut d_exponential = Wide::splat(0.0);
            if SOLAR {
                let cp_denominator = w.secant[layer] - h.k[layer];
                let inv_cp_denominator = 1.0 / cp_denominator;
                let cp_numerator = h.omega[layer] - p.exponential[layer];
                w.d_transmission[layer] += d_cp * cp_numerator * inv_cp_denominator;
                let d_cp_numerator = d_cp * w.transmission[layer] * inv_cp_denominator;
                let d_cp_denominator = -d_cp * p.cp[layer] * inv_cp_denominator;
                d_omega += d_cp_numerator;
                d_exponential -= d_cp_numerator;
                w.d_secant[layer] += d_cp_denominator;
                d_k -= d_cp_denominator;

                let cm_denominator = w.secant[layer] + h.k[layer];
                let inv_cm_denominator = 1.0 / cm_denominator;
                let cm_numerator = 1.0 - h.omega[layer] * p.exponential[layer];
                w.d_transmission[layer] += d_cm * cm_numerator * inv_cm_denominator;
                let d_cm_numerator = d_cm * w.transmission[layer] * inv_cm_denominator;
                let d_cm_denominator = -d_cm * p.cm[layer] * inv_cm_denominator;
                d_omega -= d_cm_numerator * p.exponential[layer];
                d_exponential -= d_cm_numerator * h.omega[layer];
                w.d_secant[layer] += d_cm_denominator;
                d_k += d_cm_denominator;
                w.d_secant[layer] -= d_exponential * p.exponential[layer] * w.od[layer];
                w.d_od[layer] -= d_exponential * p.exponential[layer] * w.secant[layer];
            } else {
                let cp_denominator = w.thermal_b1[layer] - h.k[layer];
                let inv_cp_denominator = 1.0 / cp_denominator;
                let cp_numerator = h.omega[layer] - p.exponential[layer];
                w.d_thermal_b0[layer] += d_cp * cp_numerator * inv_cp_denominator;
                let d_cp_numerator = d_cp * w.thermal_b0[layer] * inv_cp_denominator;
                let d_cp_denominator = -d_cp * p.cp[layer] * inv_cp_denominator;
                d_omega += d_cp_numerator;
                d_exponential -= d_cp_numerator;
                w.d_thermal_b1[layer] += d_cp_denominator;
                d_k -= d_cp_denominator;

                let cm_denominator = w.thermal_b1[layer] + h.k[layer];
                let inv_cm_denominator = 1.0 / cm_denominator;
                let cm_numerator = 1.0 - h.omega[layer] * p.exponential[layer];
                w.d_thermal_b0[layer] += d_cm * cm_numerator * inv_cm_denominator;
                let d_cm_numerator = d_cm * w.thermal_b0[layer] * inv_cm_denominator;
                let d_cm_denominator = -d_cm * p.cm[layer] * inv_cm_denominator;
                d_omega -= d_cm_numerator * p.exponential[layer];
                d_exponential -= d_cm_numerator * h.omega[layer];
                w.d_thermal_b1[layer] += d_cm_denominator;
                d_k += d_cm_denominator;
                w.d_thermal_b1[layer] -= d_exponential * p.exponential[layer] * w.od[layer];
                w.d_od[layer] -= d_exponential * p.exponential[layer] * w.thermal_b1[layer];
            }

            d_xp += d_norm * 2.0 * mu * h.xp[layer];
            d_xm -= d_norm * 2.0 * mu * h.xm[layer];
            w.d_od[layer] -= d_omega * h.omega[layer] * h.k[layer];
            d_k -= d_omega * h.omega[layer] * w.od[layer];
            let inv_k = 1.0 / h.k[layer];
            let mut d_s = (-d_xp + d_xm) * 0.5 * inv_k;
            d_k += (d_xp - d_xm) * 0.5 * h.s[layer] * inv_k * inv_k;
            d_s += d_k * 0.5 * h.d[layer] * inv_k;
            let d_d = d_k * 0.5 * h.s[layer] * inv_k;
            if az == 0 {
                w.d_ssa[layer] += d_d * w.b1[layer] * mu + d_s / mu;
                w.d_b1[layer] += d_d * w.ssa[layer] * mu;
            } else {
                w.d_ssa[layer] += d_s * w.b1[layer] * (1.0 - mu * mu) / (2.0 * mu);
                w.d_b1[layer] += d_s * w.ssa[layer] * (1.0 - mu * mu) / (2.0 * mu);
            }
        }
    }
}

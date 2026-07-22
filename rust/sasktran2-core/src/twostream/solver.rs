use super::types::{
    AtmosphereAdjoints, AtmosphereBatch, AtmosphereJacobians, Geometry, LayerAdjoints, LayerInputs,
    RadianceBatch, SourceMode, SphericalGeometry, TwoStreamError, View,
};

const FOUR_PI: f64 = 4.0 * std::f64::consts::PI;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum ExecutionPolicy {
    /// Do not create parallel work. Intended for embedding in an outer engine.
    Serial,
    /// Parallelize independent wavelength tiles with Rayon.
    #[default]
    Rayon,
}

/// Reusable allocation owned by the caller.
#[derive(Clone, Debug, Default)]
pub struct Workspace {
    scratch: ScalarWorkspace,
    explicit: super::explicit::ExplicitWorkspace,
    spherical: super::explicit::SphericalExplicitWorkspace,
    #[cfg(test)]
    pub(super) reverse: super::reverse::ReverseWorkspace,
}

impl Workspace {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn capacity_bytes(&self) -> usize {
        let bytes = self.scratch.capacity_bytes()
            + self.explicit.capacity_bytes()
            + self.spherical.capacity_bytes();
        #[cfg(test)]
        {
            bytes + self.reverse.capacity_bytes()
        }
        #[cfg(not(test))]
        {
            bytes
        }
    }
}

#[derive(Clone, Debug)]
pub struct TwoStreamSolver {
    geometry: Geometry,
    mode: SourceMode,
    execution: ExecutionPolicy,
}

impl TwoStreamSolver {
    pub fn new(geometry: Geometry, mode: SourceMode) -> Result<Self, TwoStreamError> {
        validate_geometry(&geometry)?;
        Ok(Self {
            geometry,
            mode,
            execution: ExecutionPolicy::Rayon,
        })
    }

    pub fn with_execution_policy(mut self, execution: ExecutionPolicy) -> Self {
        self.execution = execution;
        self
    }

    pub fn geometry(&self) -> &Geometry {
        &self.geometry
    }

    pub fn source_mode(&self) -> SourceMode {
        self.mode
    }

    /// Convert level quantities to the exact layer parameterization used by
    /// the C++ two-stream source.
    pub fn prepare(&self, atmosphere: &AtmosphereBatch) -> Result<LayerInputs, TwoStreamError> {
        validate_atmosphere(&self.geometry, self.mode, atmosphere)?;
        let n = self.geometry.num_layers();
        let w = atmosphere.num_wavelengths;
        let mut od = vec![0.0; n * w];
        let mut ssa = vec![0.0; n * w];
        let mut b1 = vec![0.0; n * w];
        let mut thermal_b0 = (self.mode == SourceMode::Thermal).then(|| vec![0.0; n * w]);
        let mut thermal_b1 = (self.mode == SourceMode::Thermal).then(|| vec![0.0; n * w]);

        #[cfg(feature = "simd")]
        super::simd::prepare_optical(
            self.execution == ExecutionPolicy::Rayon,
            n,
            w,
            &self.geometry.layer_thickness,
            &atmosphere.extinction,
            &atmosphere.single_scatter_albedo,
            &atmosphere.first_legendre,
            &mut od,
            &mut ssa,
            &mut b1,
        );

        #[cfg(not(feature = "simd"))]
        for layer in 0..n {
            let dz = self.geometry.layer_thickness[layer];
            for wave in 0..w {
                let top = layer * w + wave;
                let bottom = (layer + 1) * w + wave;
                let out = layer * w + wave;
                let k_top = atmosphere.extinction[top];
                let k_bottom = atmosphere.extinction[bottom];
                let ks_top = k_top * atmosphere.single_scatter_albedo[top];
                let ks_bottom = k_bottom * atmosphere.single_scatter_albedo[bottom];
                let avg_k = 0.5 * (k_top + k_bottom);
                let avg_ks = 0.5 * (ks_top + ks_bottom);
                let avg_b1 = 0.5
                    * (ks_top * atmosphere.first_legendre[top]
                        + ks_bottom * atmosphere.first_legendre[bottom]);
                od[out] = avg_k * dz;
                ssa[out] = if avg_k > 0.0 {
                    (avg_ks / avg_k).min(1.0 - 1.0e-9)
                } else {
                    0.0
                };
                b1[out] = if avg_ks > 0.0 { avg_b1 / avg_ks } else { 0.0 };
            }
        }

        if self.mode == SourceMode::Thermal {
            let emission = atmosphere.emission.as_ref().expect("validated");
            for layer in 0..n {
                for wave in 0..w {
                    let top = layer * w + wave;
                    let bottom = (layer + 1) * w + wave;
                    let out = layer * w + wave;
                    thermal_b0.as_mut().unwrap()[out] = emission[top];
                    thermal_b1.as_mut().unwrap()[out] =
                        super::thermal_profile_slope(emission[top], emission[bottom], od[out]);
                }
            }
        }

        let (transmission, average_secant) = if self.mode == SourceMode::Solar {
            let irradiance = atmosphere.solar_irradiance.as_ref().expect("validated");
            let mut transmission = vec![0.0; (n + 1) * w];
            let mut average_secant = vec![0.0; n * w];
            let mut slant = vec![0.0; (n + 1) * w];
            transmission[..w].copy_from_slice(irradiance);
            match self.execution {
                ExecutionPolicy::Serial => {
                    for boundary in 0..n {
                        prepare_solar_boundary(
                            boundary,
                            n,
                            w,
                            &self.geometry.chapman_factors,
                            &od,
                            irradiance,
                            &mut slant[(boundary + 1) * w..(boundary + 2) * w],
                            &mut transmission[(boundary + 1) * w..(boundary + 2) * w],
                        );
                    }
                }
                ExecutionPolicy::Rayon => {
                    use rayon::prelude::*;
                    slant[w..]
                        .par_chunks_mut(w)
                        .zip(transmission[w..].par_chunks_mut(w))
                        .enumerate()
                        .for_each(|(boundary, (slant_row, transmission_row))| {
                            prepare_solar_boundary(
                                boundary,
                                n,
                                w,
                                &self.geometry.chapman_factors,
                                &od,
                                irradiance,
                                slant_row,
                                transmission_row,
                            );
                        });
                }
            }
            match self.execution {
                ExecutionPolicy::Serial => {
                    for layer in 0..n {
                        prepare_average_secant_row(
                            layer,
                            w,
                            &slant,
                            &od,
                            &mut average_secant[layer * w..(layer + 1) * w],
                        );
                    }
                }
                ExecutionPolicy::Rayon => {
                    use rayon::prelude::*;
                    average_secant
                        .par_chunks_mut(w)
                        .enumerate()
                        .for_each(|(layer, output)| {
                            prepare_average_secant_row(layer, w, &slant, &od, output);
                        });
                }
            }
            (Some(transmission), Some(average_secant))
        } else {
            (None, None)
        };

        Ok(LayerInputs {
            num_layers: n,
            num_wavelengths: w,
            optical_depth: od,
            single_scatter_albedo: ssa,
            first_legendre: b1,
            transmission,
            average_secant,
            thermal_b0,
            thermal_b1,
            surface_albedo: atmosphere.surface_albedo.clone(),
            surface_emission: atmosphere.surface_emission.clone(),
        })
    }

    pub fn solve(
        &self,
        inputs: &LayerInputs,
        views: &[View],
        workspace: &mut Workspace,
    ) -> Result<RadianceBatch, TwoStreamError> {
        validate_layer_inputs(&self.geometry, self.mode, inputs, views)?;
        let mut values = vec![0.0; views.len() * inputs.num_wavelengths];
        let mut wavelength_major = vec![0.0; views.len() * inputs.num_wavelengths];

        // The serial path owns the reusable workspace. The parallel path uses
        // one workspace per Rayon worker and writes disjoint wavelengths.
        match self.execution {
            ExecutionPolicy::Serial => {
                for wave in 0..inputs.num_wavelengths {
                    solve_wavelength_column(
                        &self.geometry,
                        self.mode,
                        inputs,
                        views,
                        wave,
                        &mut workspace.scratch,
                        &mut wavelength_major[wave * views.len()..(wave + 1) * views.len()],
                    );
                }
            }
            ExecutionPolicy::Rayon => {
                use rayon::prelude::*;
                wavelength_major
                    .par_chunks_mut(views.len())
                    .enumerate()
                    .for_each_init(ScalarWorkspace::default, |scratch, (wave, column)| {
                        solve_wavelength_column(
                            &self.geometry,
                            self.mode,
                            inputs,
                            views,
                            wave,
                            scratch,
                            column,
                        );
                    });
            }
        }
        for wave in 0..inputs.num_wavelengths {
            for view in 0..views.len() {
                values[view * inputs.num_wavelengths + wave] =
                    wavelength_major[wave * views.len() + view];
            }
        }

        Ok(RadianceBatch {
            num_views: views.len(),
            num_wavelengths: inputs.num_wavelengths,
            values,
        })
    }

    /// Vector-Jacobian product for the layer API.
    ///
    /// The explicit adjoint advances four wavelength lanes together and
    /// reverses the complete solver, including a transposed pentadiagonal BVP
    /// solve and line-of-sight attenuation, in linear time.
    pub fn vjp_layers(
        &self,
        inputs: &LayerInputs,
        views: &[View],
        cotangent: &[f64],
        workspace: &mut Workspace,
    ) -> Result<LayerAdjoints, TwoStreamError> {
        validate_layer_inputs(&self.geometry, self.mode, inputs, views)?;
        if cotangent.len() != views.len() * inputs.num_wavelengths {
            return Err(TwoStreamError::invalid("cotangent has the wrong shape"));
        }
        Ok(super::explicit::solve_vjp(
            &self.geometry,
            self.mode,
            inputs,
            views,
            cotangent,
            self.execution,
            &mut workspace.explicit,
        )
        .1)
    }

    /// Fused radiance and layer-adjoint evaluation.
    ///
    /// This is the preferred derivative path: radiance values are retained
    /// from the SIMD forward sweep, avoiding a second scalar forward solve.
    pub fn solve_with_vjp(
        &self,
        inputs: &LayerInputs,
        views: &[View],
        cotangent: &[f64],
        workspace: &mut Workspace,
    ) -> Result<(RadianceBatch, LayerAdjoints), TwoStreamError> {
        validate_layer_inputs(&self.geometry, self.mode, inputs, views)?;
        if cotangent.len() != views.len() * inputs.num_wavelengths {
            return Err(TwoStreamError::invalid("cotangent has the wrong shape"));
        }
        Ok(super::explicit::solve_vjp(
            &self.geometry,
            self.mode,
            inputs,
            views,
            cotangent,
            self.execution,
            &mut workspace.explicit,
        ))
    }

    /// Fused radiance and atmospheric-level adjoint evaluation.
    ///
    /// In the parallel path, each worker maps its wavelength tile to
    /// atmospheric quantities before returning it. This avoids materializing
    /// and merging a batch-wide set of intermediate layer adjoints.
    pub fn solve_with_atmosphere_vjp(
        &self,
        atmosphere: &AtmosphereBatch,
        inputs: &LayerInputs,
        views: &[View],
        cotangent: &[f64],
        workspace: &mut Workspace,
    ) -> Result<(RadianceBatch, AtmosphereAdjoints), TwoStreamError> {
        validate_atmosphere(&self.geometry, self.mode, atmosphere)?;
        validate_layer_inputs(&self.geometry, self.mode, inputs, views)?;
        let nw = inputs.num_wavelengths;
        if atmosphere.num_wavelengths != nw {
            return Err(TwoStreamError::invalid(
                "atmosphere and layer inputs have different wavelength counts",
            ));
        }
        if cotangent.len() != views.len() * nw {
            return Err(TwoStreamError::invalid("cotangent has the wrong shape"));
        }

        if self.execution == ExecutionPolicy::Serial || nw < 2 * super::explicit::LANES {
            return Ok(solve_atmosphere_vjp_range(
                &self.geometry,
                self.mode,
                atmosphere,
                inputs,
                views,
                cotangent,
                0,
                nw,
                &mut workspace.explicit,
            ));
        }

        use rayon::prelude::*;
        let threads = rayon::current_num_threads().min(nw.div_ceil(super::explicit::LANES));
        if threads == 1 {
            return Ok(solve_atmosphere_vjp_range(
                &self.geometry,
                self.mode,
                atmosphere,
                inputs,
                views,
                cotangent,
                0,
                nw,
                &mut workspace.explicit,
            ));
        }
        // More tiles than workers let Rayon rebalance heterogeneous cores.
        // Keep cache-line alignment as well as SIMD alignment. Adjacent tiles
        // write disjoint columns of the same rows; splitting a cache line
        // would otherwise introduce false sharing between workers.
        let target_tiles = (4 * threads).min(nw.div_ceil(super::explicit::LANES));
        const CACHE_LINE_VALUES: usize = 64 / std::mem::size_of::<f64>();
        let tile_alignment = CACHE_LINE_VALUES.max(super::explicit::LANES);
        let tile_size = nw.div_ceil(target_tiles * tile_alignment) * tile_alignment;
        let ranges: Vec<_> = (0..nw)
            .step_by(tile_size)
            .map(|start| (start, (nw - start).min(tile_size)))
            .collect();

        let mut radiance = RadianceBatch {
            num_views: views.len(),
            num_wavelengths: nw,
            values: vec![0.0; views.len() * nw],
        };
        let mut atmosphere_gradient =
            zero_atmosphere_adjoint(self.mode, self.geometry.num_layers(), nw);
        let output_tiles = super::explicit::split_atmosphere_outputs(
            &ranges,
            &mut radiance,
            &mut atmosphere_gradient,
        );
        ranges
            .into_par_iter()
            .zip(output_tiles.into_par_iter())
            .for_each_init(
                super::explicit::ExplicitWorkspace::default,
                |local, ((start, len), mut output)| {
                    super::explicit::solve_atmosphere_tile(
                        &self.geometry,
                        self.mode,
                        atmosphere,
                        inputs,
                        views,
                        cotangent,
                        start,
                        len,
                        &mut output,
                        local,
                    );
                },
            );
        Ok((radiance, atmosphere_gradient))
    }

    /// Fused atmospheric preparation, radiance, and atmospheric VJP.
    ///
    /// This avoids materializing batch-wide [`LayerInputs`]. Each SIMD tile
    /// prepares its layer quantities locally immediately before the forward
    /// and reverse sweeps.
    pub fn solve_atmosphere_with_vjp(
        &self,
        atmosphere: &AtmosphereBatch,
        views: &[View],
        cotangent: &[f64],
        workspace: &mut Workspace,
    ) -> Result<(RadianceBatch, AtmosphereAdjoints), TwoStreamError> {
        validate_atmosphere(&self.geometry, self.mode, atmosphere)?;
        validate_views(views)?;
        let nw = atmosphere.num_wavelengths;
        if cotangent.len() != views.len() * nw {
            return Err(TwoStreamError::invalid("cotangent has the wrong shape"));
        }
        let mut radiance = RadianceBatch {
            num_views: views.len(),
            num_wavelengths: nw,
            values: vec![0.0; views.len() * nw],
        };
        let mut atmosphere_gradient =
            zero_atmosphere_adjoint(self.mode, self.geometry.num_layers(), nw);

        if self.execution == ExecutionPolicy::Serial || nw < 2 * super::explicit::LANES {
            let ranges = [(0, nw)];
            let mut outputs = super::explicit::split_atmosphere_outputs(
                &ranges,
                &mut radiance,
                &mut atmosphere_gradient,
            );
            super::explicit::solve_unprepared_atmosphere_tile(
                &self.geometry,
                self.mode,
                atmosphere,
                views,
                cotangent,
                0,
                nw,
                &mut outputs[0],
                &mut workspace.explicit,
            );
            drop(outputs);
            return Ok((radiance, atmosphere_gradient));
        }

        use rayon::prelude::*;
        let threads = rayon::current_num_threads().min(nw.div_ceil(super::explicit::LANES));
        if threads == 1 {
            let ranges = [(0, nw)];
            let mut outputs = super::explicit::split_atmosphere_outputs(
                &ranges,
                &mut radiance,
                &mut atmosphere_gradient,
            );
            super::explicit::solve_unprepared_atmosphere_tile(
                &self.geometry,
                self.mode,
                atmosphere,
                views,
                cotangent,
                0,
                nw,
                &mut outputs[0],
                &mut workspace.explicit,
            );
            drop(outputs);
            return Ok((radiance, atmosphere_gradient));
        }
        let target_tiles = (4 * threads).min(nw.div_ceil(super::explicit::LANES));
        const CACHE_LINE_VALUES: usize = 64 / std::mem::size_of::<f64>();
        let tile_alignment = CACHE_LINE_VALUES.max(super::explicit::LANES);
        let tile_size = nw.div_ceil(target_tiles * tile_alignment) * tile_alignment;
        let ranges: Vec<_> = (0..nw)
            .step_by(tile_size)
            .map(|start| (start, (nw - start).min(tile_size)))
            .collect();

        let output_tiles = super::explicit::split_atmosphere_outputs(
            &ranges,
            &mut radiance,
            &mut atmosphere_gradient,
        );
        ranges
            .into_par_iter()
            .zip(output_tiles.into_par_iter())
            .for_each_init(
                super::explicit::ExplicitWorkspace::default,
                |local, ((start, len), mut output)| {
                    super::explicit::solve_unprepared_atmosphere_tile(
                        &self.geometry,
                        self.mode,
                        atmosphere,
                        views,
                        cotangent,
                        start,
                        len,
                        &mut output,
                        local,
                    );
                },
            );
        Ok((radiance, atmosphere_gradient))
    }

    /// Fused atmospheric preparation and SIMD radiance evaluation.
    ///
    /// Unlike [`Self::prepare`] followed by [`Self::solve`], this path never
    /// materializes batch-wide layer inputs and is intended for embedding in
    /// the main engine.
    pub fn solve_atmosphere(
        &self,
        atmosphere: &AtmosphereBatch,
        views: &[View],
        workspace: &mut Workspace,
    ) -> Result<RadianceBatch, TwoStreamError> {
        validate_atmosphere(&self.geometry, self.mode, atmosphere)?;
        validate_views(views)?;
        let nw = atmosphere.num_wavelengths;
        let value_count = views
            .len()
            .checked_mul(nw)
            .ok_or_else(|| TwoStreamError::invalid("radiance allocation is too large"))?;
        let mut radiance = RadianceBatch {
            num_views: views.len(),
            num_wavelengths: nw,
            values: try_zeroed(value_count, "radiance")?,
        };
        let ranges = wavelength_tile_ranges(self.execution, nw);
        let output_tiles = super::explicit::split_radiance_outputs(&ranges, &mut radiance);

        if ranges.len() == 1 {
            let (start, len) = ranges[0];
            let mut outputs = output_tiles;
            super::explicit::solve_unprepared_atmosphere_forward_tile(
                &self.geometry,
                self.mode,
                atmosphere,
                views,
                start,
                len,
                &mut outputs[0],
                &mut workspace.explicit,
            );
            drop(outputs);
            return Ok(radiance);
        }

        use rayon::prelude::*;
        ranges
            .into_par_iter()
            .zip(output_tiles.into_par_iter())
            .for_each_init(
                super::explicit::ExplicitWorkspace::default,
                |local, ((start, len), mut output)| {
                    super::explicit::solve_unprepared_atmosphere_forward_tile(
                        &self.geometry,
                        self.mode,
                        atmosphere,
                        views,
                        start,
                        len,
                        &mut output,
                        local,
                    );
                },
            );
        Ok(radiance)
    }

    /// Fused atmospheric preparation, radiance, and one Jacobian per view.
    ///
    /// The atmospheric preparation, layer forward sweep, and BVP
    /// factorization are shared by every view in a wavelength tile.  Only the
    /// line-of-sight and explicit-adjoint sweeps are repeated per view.
    pub fn solve_atmosphere_with_jacobians(
        &self,
        atmosphere: &AtmosphereBatch,
        views: &[View],
        workspace: &mut Workspace,
    ) -> Result<(RadianceBatch, AtmosphereJacobians), TwoStreamError> {
        validate_atmosphere(&self.geometry, self.mode, atmosphere)?;
        validate_views(views)?;
        let nw = atmosphere.num_wavelengths;
        let num_levels = self.geometry.num_layers() + 1;
        let surface_count = views
            .len()
            .checked_mul(nw)
            .ok_or_else(|| TwoStreamError::invalid("surface Jacobian allocation is too large"))?;
        let level_count = surface_count
            .checked_mul(num_levels)
            .ok_or_else(|| TwoStreamError::invalid("level Jacobian allocation is too large"))?;
        let mut radiance = RadianceBatch {
            num_views: views.len(),
            num_wavelengths: nw,
            values: try_zeroed(surface_count, "radiance")?,
        };
        let mut jacobians = AtmosphereJacobians {
            num_views: views.len(),
            num_levels,
            num_wavelengths: nw,
            extinction: try_zeroed(level_count, "extinction Jacobian")?,
            single_scatter_albedo: try_zeroed(level_count, "SSA Jacobian")?,
            first_legendre: try_zeroed(level_count, "phase Jacobian")?,
            emission: if self.mode == SourceMode::Thermal {
                Some(try_zeroed(level_count, "emission Jacobian")?)
            } else {
                None
            },
            surface_albedo: try_zeroed(surface_count, "surface-albedo Jacobian")?,
            surface_emission: if self.mode == SourceMode::Thermal {
                Some(try_zeroed(surface_count, "surface-emission Jacobian")?)
            } else {
                None
            },
        };
        let ranges = wavelength_tile_ranges(self.execution, nw);
        let output_tiles = super::explicit::split_atmosphere_jacobian_outputs(
            &ranges,
            &mut radiance,
            &mut jacobians,
        );

        if ranges.len() == 1 {
            let (start, len) = ranges[0];
            let mut outputs = output_tiles;
            super::explicit::solve_unprepared_atmosphere_jacobian_tile(
                &self.geometry,
                self.mode,
                atmosphere,
                views,
                start,
                len,
                &mut outputs[0],
                &mut workspace.explicit,
            );
            drop(outputs);
            return Ok((radiance, jacobians));
        }

        use rayon::prelude::*;
        ranges
            .into_par_iter()
            .zip(output_tiles.into_par_iter())
            .for_each_init(
                super::explicit::ExplicitWorkspace::default,
                |local, ((start, len), mut output)| {
                    super::explicit::solve_unprepared_atmosphere_jacobian_tile(
                        &self.geometry,
                        self.mode,
                        atmosphere,
                        views,
                        start,
                        len,
                        &mut output,
                        local,
                    );
                },
            );
        Ok((radiance, jacobians))
    }

    /// Evaluate the local two-stream source along traced spherical rays.
    pub fn solve_spherical_atmosphere(
        &self,
        atmosphere: &AtmosphereBatch,
        spherical: &SphericalGeometry,
        workspace: &mut Workspace,
    ) -> Result<RadianceBatch, TwoStreamError> {
        validate_spherical_geometry(&self.geometry, self.mode, spherical)?;
        validate_atmosphere(&spherical.columns[0], self.mode, atmosphere)?;
        let nw = atmosphere.num_wavelengths;
        let num_views = spherical.rays.num_views();
        let value_count = num_views
            .checked_mul(nw)
            .ok_or_else(|| TwoStreamError::invalid("radiance allocation is too large"))?;
        let mut radiance = RadianceBatch {
            num_views,
            num_wavelengths: nw,
            values: try_zeroed(value_count, "radiance")?,
        };
        let ranges = wavelength_tile_ranges(self.execution, nw);
        let output_tiles = super::explicit::split_radiance_outputs(&ranges, &mut radiance);
        if ranges.len() == 1 {
            let (start, len) = ranges[0];
            let mut outputs = output_tiles;
            super::explicit::solve_unprepared_spherical_forward_tile(
                spherical,
                self.mode,
                atmosphere,
                start,
                len,
                &mut outputs[0],
                &mut workspace.spherical,
            );
            drop(outputs);
            return Ok(radiance);
        }

        use rayon::prelude::*;
        ranges
            .into_par_iter()
            .zip(output_tiles.into_par_iter())
            .for_each_init(
                super::explicit::SphericalExplicitWorkspace::default,
                |local, ((start, len), mut output)| {
                    super::explicit::solve_unprepared_spherical_forward_tile(
                        spherical,
                        self.mode,
                        atmosphere,
                        start,
                        len,
                        &mut output,
                        local,
                    );
                },
            );
        Ok(radiance)
    }

    /// Spherical radiance and one explicit atmospheric Jacobian per ray.
    pub fn solve_spherical_atmosphere_with_jacobians(
        &self,
        atmosphere: &AtmosphereBatch,
        spherical: &SphericalGeometry,
        workspace: &mut Workspace,
    ) -> Result<(RadianceBatch, AtmosphereJacobians), TwoStreamError> {
        validate_spherical_geometry(&self.geometry, self.mode, spherical)?;
        validate_atmosphere(&spherical.columns[0], self.mode, atmosphere)?;
        let nw = atmosphere.num_wavelengths;
        let num_views = spherical.rays.num_views();
        let num_levels = spherical.rays.num_levels;
        let surface_count = num_views
            .checked_mul(nw)
            .ok_or_else(|| TwoStreamError::invalid("surface Jacobian allocation is too large"))?;
        let level_count = surface_count
            .checked_mul(num_levels)
            .ok_or_else(|| TwoStreamError::invalid("level Jacobian allocation is too large"))?;
        let mut radiance = RadianceBatch {
            num_views,
            num_wavelengths: nw,
            values: try_zeroed(surface_count, "radiance")?,
        };
        let mut jacobians = AtmosphereJacobians {
            num_views,
            num_levels,
            num_wavelengths: nw,
            extinction: try_zeroed(level_count, "extinction Jacobian")?,
            single_scatter_albedo: try_zeroed(level_count, "SSA Jacobian")?,
            first_legendre: try_zeroed(level_count, "phase Jacobian")?,
            emission: if self.mode == SourceMode::Thermal {
                Some(try_zeroed(level_count, "emission Jacobian")?)
            } else {
                None
            },
            surface_albedo: try_zeroed(surface_count, "surface-albedo Jacobian")?,
            surface_emission: if self.mode == SourceMode::Thermal {
                Some(try_zeroed(surface_count, "surface-emission Jacobian")?)
            } else {
                None
            },
        };
        let ranges = wavelength_tile_ranges(self.execution, nw);
        let output_tiles = super::explicit::split_atmosphere_jacobian_outputs(
            &ranges,
            &mut radiance,
            &mut jacobians,
        );
        if ranges.len() == 1 {
            let (start, len) = ranges[0];
            let mut outputs = output_tiles;
            super::explicit::solve_unprepared_spherical_jacobian_tile(
                spherical,
                self.mode,
                atmosphere,
                start,
                len,
                &mut outputs[0],
                &mut workspace.spherical,
            );
            drop(outputs);
            return Ok((radiance, jacobians));
        }

        use rayon::prelude::*;
        ranges
            .into_par_iter()
            .zip(output_tiles.into_par_iter())
            .for_each_init(
                super::explicit::SphericalExplicitWorkspace::default,
                |local, ((start, len), mut output)| {
                    super::explicit::solve_unprepared_spherical_jacobian_tile(
                        spherical,
                        self.mode,
                        atmosphere,
                        start,
                        len,
                        &mut output,
                        local,
                    );
                },
            );
        Ok((radiance, jacobians))
    }

    /// Map layer-input adjoints back to atmospheric level quantities.
    pub fn map_adjoint_to_atmosphere(
        &self,
        atmosphere: &AtmosphereBatch,
        inputs: &LayerInputs,
        layer: &LayerAdjoints,
    ) -> Result<AtmosphereAdjoints, TwoStreamError> {
        validate_atmosphere(&self.geometry, self.mode, atmosphere)?;
        map_layer_adjoint(
            &self.geometry,
            self.mode,
            atmosphere,
            inputs,
            layer,
            self.execution,
        )
    }
}

fn try_zeroed(count: usize, name: &str) -> Result<Vec<f64>, TwoStreamError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| TwoStreamError::invalid(format!("unable to allocate {name}")))?;
    values.resize(count, 0.0);
    Ok(values)
}

fn wavelength_tile_ranges(execution: ExecutionPolicy, nw: usize) -> Vec<(usize, usize)> {
    if execution == ExecutionPolicy::Serial || nw < 2 * super::explicit::LANES {
        return vec![(0, nw)];
    }
    let threads = rayon::current_num_threads().min(nw.div_ceil(super::explicit::LANES));
    if threads == 1 {
        return vec![(0, nw)];
    }
    let target_tiles = (4 * threads).min(nw.div_ceil(super::explicit::LANES));
    const CACHE_LINE_VALUES: usize = 64 / std::mem::size_of::<f64>();
    let tile_alignment = CACHE_LINE_VALUES.max(super::explicit::LANES);
    let tile_size = nw.div_ceil(target_tiles * tile_alignment) * tile_alignment;
    (0..nw)
        .step_by(tile_size)
        .map(|start| (start, (nw - start).min(tile_size)))
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn prepare_solar_boundary(
    boundary: usize,
    n: usize,
    nw: usize,
    chapman: &[f64],
    optical_depth: &[f64],
    irradiance: &[f64],
    slant: &mut [f64],
    transmission: &mut [f64],
) {
    #[cfg(feature = "simd")]
    super::simd::prepare_solar_boundary(
        &chapman[boundary * n..(boundary + 1) * n],
        nw,
        optical_depth,
        irradiance,
        slant,
        transmission,
    );

    #[cfg(not(feature = "simd"))]
    for layer in 0..n {
        let factor = chapman[boundary * n + layer];
        if factor == 0.0 {
            continue;
        }
        for (tau, optical_depth) in slant
            .iter_mut()
            .zip(&optical_depth[layer * nw..(layer + 1) * nw])
        {
            *tau += factor * optical_depth;
        }
    }
    #[cfg(not(feature = "simd"))]
    for ((transmission, slant), irradiance) in transmission.iter_mut().zip(slant).zip(irradiance) {
        *slant = -*slant;
        *transmission = slant.exp() * irradiance;
    }
}

fn prepare_average_secant_row(
    layer: usize,
    nw: usize,
    slant: &[f64],
    optical_depth: &[f64],
    output: &mut [f64],
) {
    let top = &slant[layer * nw..(layer + 1) * nw];
    let bottom = &slant[(layer + 1) * nw..(layer + 2) * nw];
    let optical_depth = &optical_depth[layer * nw..(layer + 1) * nw];
    #[cfg(feature = "simd")]
    super::simd::prepare_average_secant(top, bottom, optical_depth, output);
    #[cfg(not(feature = "simd"))]
    for (((output, top), bottom), optical_depth) in
        output.iter_mut().zip(top).zip(bottom).zip(optical_depth)
    {
        *output = if *optical_depth > 0.0 {
            (top - bottom) / optical_depth
        } else {
            0.0
        };
    }
}

#[allow(clippy::too_many_arguments)]
fn solve_atmosphere_vjp_range(
    geometry: &Geometry,
    mode: SourceMode,
    atmosphere: &AtmosphereBatch,
    inputs: &LayerInputs,
    views: &[View],
    cotangent: &[f64],
    start: usize,
    len: usize,
    workspace: &mut super::explicit::ExplicitWorkspace,
) -> (RadianceBatch, AtmosphereAdjoints) {
    super::explicit::solve_atmosphere_range(
        geometry, mode, atmosphere, inputs, views, cotangent, start, len, workspace,
    )
}

#[derive(Clone, Debug, Default)]
struct Homogeneous {
    k: Vec<f64>,
    x_plus: Vec<f64>,
    x_minus: Vec<f64>,
    omega: Vec<f64>,
}

#[derive(Clone, Debug, Default)]
struct Particular {
    a_plus: Vec<f64>,
    a_minus: Vec<f64>,
    a_thermal: Vec<f64>,
    g_plus_top: Vec<f64>,
    g_plus_bottom: Vec<f64>,
    g_minus_top: Vec<f64>,
    g_minus_bottom: Vec<f64>,
}

#[derive(Clone, Debug, Default)]
struct Bvp {
    e: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
    a: Vec<f64>,
    b: Vec<f64>,
    rhs: Vec<f64>,
    mu: Vec<f64>,
    alpha: Vec<f64>,
    beta: Vec<f64>,
    gamma: Vec<f64>,
    z: Vec<f64>,
}

#[derive(Clone, Debug, Default)]
struct ScalarWorkspace {
    od: Vec<f64>,
    ssa: Vec<f64>,
    b1: Vec<f64>,
    transmission: Vec<f64>,
    secant: Vec<f64>,
    thermal_b0: Vec<f64>,
    thermal_b1: Vec<f64>,
    homogeneous: [Homogeneous; 2],
    particular: [Particular; 2],
    bvp: [Bvp; 2],
    source: Vec<f64>,
}

impl ScalarWorkspace {
    fn resize(&mut self, n: usize, mode: SourceMode) {
        resize_zero(&mut self.od, n);
        resize_zero(&mut self.ssa, n);
        resize_zero(&mut self.b1, n);
        resize_zero(&mut self.source, n);
        resize_zero(&mut self.transmission, n + 1);
        resize_zero(&mut self.secant, n);
        resize_zero(&mut self.thermal_b0, n);
        resize_zero(&mut self.thermal_b1, n);
        let naz = if mode == SourceMode::Solar { 2 } else { 1 };
        for az in 0..naz {
            let h = &mut self.homogeneous[az];
            resize_zero(&mut h.k, n);
            resize_zero(&mut h.x_plus, n);
            resize_zero(&mut h.x_minus, n);
            resize_zero(&mut h.omega, n);
            let p = &mut self.particular[az];
            resize_zero(&mut p.a_plus, n);
            resize_zero(&mut p.a_minus, n);
            resize_zero(&mut p.a_thermal, n);
            resize_zero(&mut p.g_plus_top, n);
            resize_zero(&mut p.g_plus_bottom, n);
            resize_zero(&mut p.g_minus_top, n);
            resize_zero(&mut p.g_minus_bottom, n);
            let bvp = &mut self.bvp[az];
            for v in [
                &mut bvp.e,
                &mut bvp.c,
                &mut bvp.d,
                &mut bvp.a,
                &mut bvp.b,
                &mut bvp.rhs,
                &mut bvp.mu,
                &mut bvp.alpha,
                &mut bvp.beta,
                &mut bvp.gamma,
                &mut bvp.z,
            ] {
                resize_zero(v, 2 * n);
            }
        }
    }

    fn capacity_bytes(&self) -> usize {
        let mut count = self.od.capacity()
            + self.ssa.capacity()
            + self.b1.capacity()
            + self.transmission.capacity()
            + self.secant.capacity()
            + self.thermal_b0.capacity()
            + self.thermal_b1.capacity()
            + self.source.capacity();
        for h in &self.homogeneous {
            count +=
                h.k.capacity() + h.x_plus.capacity() + h.x_minus.capacity() + h.omega.capacity();
        }
        for p in &self.particular {
            count += p.a_plus.capacity()
                + p.a_minus.capacity()
                + p.a_thermal.capacity()
                + p.g_plus_top.capacity()
                + p.g_plus_bottom.capacity()
                + p.g_minus_top.capacity()
                + p.g_minus_bottom.capacity();
        }
        for b in &self.bvp {
            count += b.e.capacity()
                + b.c.capacity()
                + b.d.capacity()
                + b.a.capacity()
                + b.b.capacity()
                + b.rhs.capacity()
                + b.mu.capacity()
                + b.alpha.capacity()
                + b.beta.capacity()
                + b.gamma.capacity()
                + b.z.capacity();
        }
        count * std::mem::size_of::<f64>()
    }
}

fn resize_zero(values: &mut Vec<f64>, len: usize) {
    values.resize(len, 0.0);
    values.fill(0.0);
}

fn solve_wavelength_column(
    geometry: &Geometry,
    mode: SourceMode,
    inputs: &LayerInputs,
    views: &[View],
    wave: usize,
    scratch: &mut ScalarWorkspace,
    output: &mut [f64],
) {
    let n = inputs.num_layers;
    let nw = inputs.num_wavelengths;
    scratch.resize(n, mode);
    for layer in 0..n {
        let idx = layer * nw + wave;
        scratch.od[layer] = inputs.optical_depth[idx];
        scratch.ssa[layer] = inputs.single_scatter_albedo[idx];
        scratch.b1[layer] = inputs.first_legendre[idx];
        if mode == SourceMode::Solar {
            scratch.secant[layer] = inputs.average_secant.as_ref().unwrap()[idx];
        } else {
            scratch.thermal_b0[layer] = inputs.thermal_b0.as_ref().unwrap()[idx];
            scratch.thermal_b1[layer] = inputs.thermal_b1.as_ref().unwrap()[idx];
        }
    }
    if mode == SourceMode::Solar {
        for boundary in 0..=n {
            scratch.transmission[boundary] =
                inputs.transmission.as_ref().unwrap()[boundary * nw + wave];
        }
    }

    solve_layers(geometry, mode, inputs.surface_albedo[wave], scratch);
    solve_bvp(
        geometry,
        mode,
        inputs.surface_albedo[wave],
        inputs
            .surface_emission
            .as_ref()
            .map_or(0.0, |values| values[wave]),
        scratch,
    );
    for (view, target) in views.iter().zip(output.iter_mut()) {
        *target = post_process(
            geometry,
            mode,
            inputs.surface_albedo[wave],
            inputs
                .surface_emission
                .as_ref()
                .map_or(0.0, |values| values[wave]),
            *view,
            scratch,
        );
    }
}

fn solve_layers(geometry: &Geometry, mode: SourceMode, _albedo: f64, w: &mut ScalarWorkspace) {
    let n = w.od.len();
    let mu = geometry.quadrature_cosine;
    let naz = if mode == SourceMode::Solar { 2 } else { 1 };
    let sin_product = ((1.0 - mu * mu) * (1.0 - geometry.solar_cosine.powi(2))).sqrt();
    for az in 0..naz {
        for layer in 0..n {
            let ssa = w.ssa[layer];
            let b1 = w.b1[layer];
            let (d, s) = if az == 0 {
                (ssa * b1 * mu - 1.0 / mu, (ssa - 1.0) / mu)
            } else {
                (-1.0 / mu, (ssa * b1 * (1.0 - mu * mu) - 2.0) / (2.0 * mu))
            };
            let k = (s * d).sqrt();
            let xp = 0.5 * (1.0 - s / k);
            let xm = 0.5 * (1.0 + s / k);
            let omega = (-k * w.od[layer]).exp();
            w.homogeneous[az].k[layer] = k;
            w.homogeneous[az].x_plus[layer] = xp;
            w.homogeneous[az].x_minus[layer] = xm;
            w.homogeneous[az].omega[layer] = omega;
            let norm = mu * (xp * xp - xm * xm);

            if mode == SourceMode::Solar {
                let (qp, qm) = if az == 0 {
                    (
                        ssa * (1.0 + b1 * geometry.solar_cosine * mu) / FOUR_PI,
                        ssa * (1.0 - b1 * geometry.solar_cosine * mu) / FOUR_PI,
                    )
                } else {
                    let q = ssa * b1 * sin_product / FOUR_PI;
                    (q, q)
                };
                let ap = (qp * xp + qm * xm) / norm;
                let am = (qm * xp + qp * xm) / norm;
                let expsec = (-w.secant[layer] * w.od[layer]).exp();
                let cp = w.transmission[layer]
                    * super::explicit::exp_difference_ratio_scalar(
                        omega,
                        expsec,
                        k,
                        w.secant[layer],
                        w.od[layer],
                    )
                    .0;
                let cm = w.transmission[layer]
                    * super::explicit::exp_difference_ratio_scalar(
                        1.0,
                        omega * expsec,
                        0.0,
                        w.secant[layer] + k,
                        w.od[layer],
                    )
                    .0;
                let p = &mut w.particular[az];
                p.a_plus[layer] = ap;
                p.a_minus[layer] = am;
                p.g_plus_top[layer] = am * cm * xm;
                p.g_plus_bottom[layer] = ap * cp * xp;
                p.g_minus_top[layer] = am * cm * xp;
                p.g_minus_bottom[layer] = ap * cp * xm;
            } else {
                let at = (1.0 - ssa) * (xp + xm) / norm;
                let exp_thermal = (-w.thermal_b1[layer] * w.od[layer]).exp();
                let cp = w.thermal_b0[layer]
                    * super::explicit::exp_difference_ratio_scalar(
                        omega,
                        exp_thermal,
                        k,
                        w.thermal_b1[layer],
                        w.od[layer],
                    )
                    .0;
                let cm = w.thermal_b0[layer]
                    * super::explicit::exp_difference_ratio_scalar(
                        1.0,
                        omega * exp_thermal,
                        0.0,
                        w.thermal_b1[layer] + k,
                        w.od[layer],
                    )
                    .0;
                let p = &mut w.particular[az];
                p.a_thermal[layer] = at;
                p.g_plus_top[layer] = at * cm * xm;
                p.g_plus_bottom[layer] = at * cp * xp;
                p.g_minus_top[layer] = at * cm * xp;
                p.g_minus_bottom[layer] = at * cp * xm;
            }
        }
    }
}

fn solve_bvp(
    geometry: &Geometry,
    mode: SourceMode,
    albedo: f64,
    thermal_surface: f64,
    w: &mut ScalarWorkspace,
) {
    let n = w.od.len();
    let mu = geometry.quadrature_cosine;
    let naz = if mode == SourceMode::Solar { 2 } else { 1 };
    for az in 0..naz {
        let h = &w.homogeneous[az];
        let p = &w.particular[az];
        let bvp = &mut w.bvp[az];
        bvp.rhs[0] = -p.g_plus_top[0];
        for layer in 0..n - 1 {
            bvp.rhs[2 * layer + 1] = p.g_minus_top[layer + 1] - p.g_minus_bottom[layer];
            bvp.rhs[2 * layer + 2] = p.g_plus_top[layer + 1] - p.g_plus_bottom[layer];
        }
        let last = 2 * n - 1;
        let delta = if az == 0 { 1.0 } else { 0.0 };
        let direct_boundary_source = if mode == SourceMode::Solar {
            delta * geometry.solar_cosine * albedo / std::f64::consts::PI * w.transmission[n]
        } else {
            thermal_surface
        };
        bvp.rhs[last] = direct_boundary_source
            - (p.g_minus_bottom[n - 1] - 2.0 * delta * mu * albedo * p.g_plus_bottom[n - 1]);

        bvp.d[0] = h.x_plus[0];
        bvp.a[0] = h.x_minus[0] * h.omega[0];
        for layer in 0..n - 1 {
            let row = 2 * layer;
            bvp.c[row + 1] = h.x_minus[layer] * h.omega[layer];
            bvp.d[row + 1] = h.x_plus[layer];
            bvp.a[row + 1] = -h.x_minus[layer + 1];
            bvp.b[row + 1] = -h.x_plus[layer + 1] * h.omega[layer + 1];
            bvp.e[row + 2] = h.x_plus[layer] * h.omega[layer];
            bvp.c[row + 2] = h.x_minus[layer];
            bvp.d[row + 2] = -h.x_plus[layer + 1];
            bvp.a[row + 2] = -h.x_minus[layer + 1] * h.omega[layer + 1];
        }
        bvp.c[last] =
            (h.x_minus[n - 1] - 2.0 * mu * albedo * delta * h.x_plus[n - 1]) * h.omega[n - 1];
        bvp.d[last] = h.x_plus[n - 1] - 2.0 * mu * albedo * delta * h.x_minus[n - 1];
        pentadiagonal_solve(bvp);
    }
}

fn pentadiagonal_solve(bvp: &mut Bvp) {
    let n = bvp.d.len();
    bvp.mu[0] = bvp.d[0];
    bvp.alpha[0] = bvp.a[0] / bvp.mu[0];
    bvp.beta[0] = bvp.b[0] / bvp.mu[0];
    bvp.z[0] = bvp.rhs[0] / bvp.mu[0];
    if n > 1 {
        bvp.gamma[1] = bvp.c[1];
        bvp.mu[1] = bvp.d[1] - bvp.alpha[0] * bvp.gamma[1];
        bvp.alpha[1] = (bvp.a[1] - bvp.beta[0] * bvp.gamma[1]) / bvp.mu[1];
        bvp.beta[1] = bvp.b[1] / bvp.mu[1];
        bvp.z[1] = (bvp.rhs[1] - bvp.z[0] * bvp.gamma[1]) / bvp.mu[1];
    }
    for i in 2..n {
        bvp.gamma[i] = bvp.c[i] - bvp.alpha[i - 2] * bvp.e[i];
        bvp.mu[i] = bvp.d[i] - bvp.beta[i - 2] * bvp.e[i] - bvp.alpha[i - 1] * bvp.gamma[i];
        if i + 1 < n {
            bvp.alpha[i] = (bvp.a[i] - bvp.beta[i - 1] * bvp.gamma[i]) / bvp.mu[i];
        }
        if i + 2 < n {
            bvp.beta[i] = bvp.b[i] / bvp.mu[i];
        }
        bvp.z[i] = (bvp.rhs[i] - bvp.z[i - 2] * bvp.e[i] - bvp.z[i - 1] * bvp.gamma[i]) / bvp.mu[i];
    }
    bvp.rhs[n - 1] = bvp.z[n - 1];
    if n > 1 {
        bvp.rhs[n - 2] = bvp.z[n - 2] - bvp.alpha[n - 2] * bvp.rhs[n - 1];
    }
    for i in (0..n.saturating_sub(2)).rev() {
        bvp.rhs[i] = bvp.z[i] - bvp.alpha[i] * bvp.rhs[i + 1] - bvp.beta[i] * bvp.rhs[i + 2];
    }
}

fn post_process(
    geometry: &Geometry,
    mode: SourceMode,
    albedo: f64,
    thermal_surface: f64,
    view: View,
    w: &mut ScalarWorkspace,
) -> f64 {
    let n = w.od.len();
    let mu = geometry.quadrature_cosine;
    let vmu = view.cosine;
    w.source.fill(0.0);
    let mut attenuation = 1.0;
    let mut integrated = 0.0;
    let naz = if mode == SourceMode::Solar { 2 } else { 1 };

    for layer in 0..n {
        let beam = (-w.od[layer] / vmu).exp();
        let expsec = if mode == SourceMode::Solar {
            (-w.secant[layer] * w.od[layer]).exp()
        } else {
            0.0
        };
        let expthermal = if mode == SourceMode::Thermal {
            (-w.thermal_b1[layer] * w.od[layer]).exp()
        } else {
            0.0
        };
        let rate = if mode == SourceMode::Solar {
            w.secant[layer]
        } else {
            w.thermal_b1[layer]
        };
        let exponential = if mode == SourceMode::Solar {
            expsec
        } else {
            expthermal
        };
        let inverse_view = 1.0 / vmu;
        let source_integral = inverse_view
            * super::explicit::exp_difference_ratio_scalar(
                1.0,
                exponential * beam,
                0.0,
                rate + inverse_view,
                w.od[layer],
            )
            .0;

        for az in 0..naz {
            let h = &w.homogeneous[az];
            let p = &w.particular[az];
            let azi = (az as f64 * view.relative_azimuth).cos();
            let (lp, lm) = if az == 0 {
                (
                    0.5 * w.ssa[layer] * (1.0 - w.b1[layer] * vmu * mu),
                    0.5 * w.ssa[layer] * (1.0 + w.b1[layer] * vmu * mu),
                )
            } else if mode == SourceMode::Solar {
                let value = 0.25
                    * w.ssa[layer]
                    * w.b1[layer]
                    * ((1.0 - vmu * vmu) * (1.0 - mu * mu)).sqrt();
                (value, value)
            } else {
                (0.0, 0.0)
            };
            let yp = lp * h.x_plus[layer] + lm * h.x_minus[layer];
            let ym = lp * h.x_minus[layer] + lm * h.x_plus[layer];
            let hm = inverse_view
                * super::explicit::exp_difference_ratio_scalar(
                    h.omega[layer],
                    beam,
                    h.k[layer],
                    inverse_view,
                    w.od[layer],
                )
                .0;
            let hp = inverse_view
                * super::explicit::exp_difference_ratio_scalar(
                    1.0,
                    h.omega[layer] * beam,
                    0.0,
                    h.k[layer] + inverse_view,
                    w.od[layer],
                )
                .0;
            let dp_ratio = inverse_view
                * super::explicit::integrated_exp_difference_ratio_scalar(
                    exponential * h.omega[layer],
                    exponential * beam,
                    rate + h.k[layer],
                    rate + inverse_view,
                    w.od[layer],
                )
                .0;
            let dm_ratio = inverse_view
                * super::explicit::integrated_exp_difference_ratio_scalar(
                    h.omega[layer] * beam,
                    exponential * beam,
                    h.k[layer] + inverse_view,
                    rate + inverse_view,
                    w.od[layer],
                )
                .0;
            let v = if mode == SourceMode::Solar {
                let dp = w.transmission[layer] * dp_ratio;
                let dm = w.transmission[layer] * dm_ratio;
                p.a_plus[layer] * yp * dm + p.a_minus[layer] * ym * dp
            } else {
                let dp = w.thermal_b0[layer] * dp_ratio;
                let dm = w.thermal_b0[layer] * dm_ratio;
                p.a_thermal[layer] * (yp * dm + ym * dp)
            };
            w.source[layer] += azi
                * (w.bvp[az].rhs[2 * layer] * yp * hp + w.bvp[az].rhs[2 * layer + 1] * ym * hm + v);
            if mode == SourceMode::Thermal {
                // This is inside the azimuth loop in the C++ implementation.
                w.source[layer] += w.thermal_b0[layer] * source_integral * (1.0 - w.ssa[layer]);
            }
        }
        integrated += w.source[layer] * attenuation;
        attenuation *= beam;
    }

    let last = n - 1;
    let surface_diffuse = (w.particular[0].g_plus_bottom[last]
        + w.bvp[0].rhs[2 * last] * w.homogeneous[0].x_plus[last] * w.homogeneous[0].omega[last]
        + w.bvp[0].rhs[2 * last + 1] * w.homogeneous[0].x_minus[last])
        * 2.0
        * mu
        * albedo;
    integrated + attenuation * (surface_diffuse + thermal_surface)
}

#[cfg(test)]
#[allow(dead_code)]
fn finite_difference_vjp(
    solver: &TwoStreamSolver,
    inputs: &LayerInputs,
    views: &[View],
    cotangent: &[f64],
    workspace: &mut Workspace,
) -> Result<LayerAdjoints, TwoStreamError> {
    let mut perturbed = inputs.clone();
    let layer_len = inputs.num_layers * inputs.num_wavelengths;
    let mut result = LayerAdjoints {
        optical_depth: vec![0.0; layer_len],
        single_scatter_albedo: vec![0.0; layer_len],
        first_legendre: vec![0.0; layer_len],
        transmission: (solver.mode == SourceMode::Solar)
            .then(|| vec![0.0; (inputs.num_layers + 1) * inputs.num_wavelengths]),
        average_secant: (solver.mode == SourceMode::Solar).then(|| vec![0.0; layer_len]),
        thermal_b0: (solver.mode == SourceMode::Thermal).then(|| vec![0.0; layer_len]),
        thermal_b1: (solver.mode == SourceMode::Thermal).then(|| vec![0.0; layer_len]),
        surface_albedo: vec![0.0; inputs.num_wavelengths],
        surface_emission: (solver.mode == SourceMode::Thermal)
            .then(|| vec![0.0; inputs.num_wavelengths]),
    };

    fn objective_wave(
        solver: &TwoStreamSolver,
        inputs: &LayerInputs,
        views: &[View],
        cotangent: &[f64],
        wave: usize,
        workspace: &mut Workspace,
    ) -> f64 {
        let mut column = vec![0.0; views.len()];
        solve_wavelength_column(
            &solver.geometry,
            solver.mode,
            inputs,
            views,
            wave,
            &mut workspace.scratch,
            &mut column,
        );
        column
            .iter()
            .enumerate()
            .map(|(view, value)| value * cotangent[view * inputs.num_wavelengths + wave])
            .sum()
    }

    fn step(value: f64) -> f64 {
        f64::EPSILON.cbrt() * value.abs().max(1.0)
    }

    macro_rules! differentiate_vec {
        ($field:ident, $out:expr) => {{
            for index in 0..perturbed.$field.len() {
                let original = perturbed.$field[index];
                let h = step(original);
                let wave = index % inputs.num_wavelengths;
                perturbed.$field[index] = original + h;
                let plus = objective_wave(solver, &perturbed, views, cotangent, wave, workspace);
                perturbed.$field[index] = original - h;
                let minus = objective_wave(solver, &perturbed, views, cotangent, wave, workspace);
                perturbed.$field[index] = original;
                ($out)[index] = (plus - minus) / (2.0 * h);
            }
        }};
    }

    differentiate_vec!(optical_depth, result.optical_depth);
    differentiate_vec!(single_scatter_albedo, result.single_scatter_albedo);
    differentiate_vec!(first_legendre, result.first_legendre);
    differentiate_vec!(surface_albedo, result.surface_albedo);

    if solver.mode == SourceMode::Solar {
        for index in 0..perturbed.transmission.as_ref().unwrap().len() {
            let original = perturbed.transmission.as_ref().unwrap()[index];
            let h = step(original);
            let wave = index % inputs.num_wavelengths;
            perturbed.transmission.as_mut().unwrap()[index] = original + h;
            let plus = objective_wave(solver, &perturbed, views, cotangent, wave, workspace);
            perturbed.transmission.as_mut().unwrap()[index] = original - h;
            let minus = objective_wave(solver, &perturbed, views, cotangent, wave, workspace);
            perturbed.transmission.as_mut().unwrap()[index] = original;
            result.transmission.as_mut().unwrap()[index] = (plus - minus) / (2.0 * h);
        }
        for index in 0..layer_len {
            let original = perturbed.average_secant.as_ref().unwrap()[index];
            let h = step(original);
            let wave = index % inputs.num_wavelengths;
            perturbed.average_secant.as_mut().unwrap()[index] = original + h;
            let plus = objective_wave(solver, &perturbed, views, cotangent, wave, workspace);
            perturbed.average_secant.as_mut().unwrap()[index] = original - h;
            let minus = objective_wave(solver, &perturbed, views, cotangent, wave, workspace);
            perturbed.average_secant.as_mut().unwrap()[index] = original;
            result.average_secant.as_mut().unwrap()[index] = (plus - minus) / (2.0 * h);
        }
    }

    if solver.mode == SourceMode::Thermal {
        for index in 0..layer_len {
            let original = perturbed.thermal_b0.as_ref().unwrap()[index];
            let h = step(original);
            let wave = index % inputs.num_wavelengths;
            perturbed.thermal_b0.as_mut().unwrap()[index] = original + h;
            let plus = objective_wave(solver, &perturbed, views, cotangent, wave, workspace);
            perturbed.thermal_b0.as_mut().unwrap()[index] = original - h;
            let minus = objective_wave(solver, &perturbed, views, cotangent, wave, workspace);
            perturbed.thermal_b0.as_mut().unwrap()[index] = original;
            result.thermal_b0.as_mut().unwrap()[index] = (plus - minus) / (2.0 * h);

            let original = perturbed.thermal_b1.as_ref().unwrap()[index];
            let h = step(original);
            perturbed.thermal_b1.as_mut().unwrap()[index] = original + h;
            let plus = objective_wave(solver, &perturbed, views, cotangent, wave, workspace);
            perturbed.thermal_b1.as_mut().unwrap()[index] = original - h;
            let minus = objective_wave(solver, &perturbed, views, cotangent, wave, workspace);
            perturbed.thermal_b1.as_mut().unwrap()[index] = original;
            result.thermal_b1.as_mut().unwrap()[index] = (plus - minus) / (2.0 * h);
        }
        for wave in 0..inputs.num_wavelengths {
            let original = perturbed.surface_emission.as_ref().unwrap()[wave];
            let h = step(original);
            perturbed.surface_emission.as_mut().unwrap()[wave] = original + h;
            let plus = objective_wave(solver, &perturbed, views, cotangent, wave, workspace);
            perturbed.surface_emission.as_mut().unwrap()[wave] = original - h;
            let minus = objective_wave(solver, &perturbed, views, cotangent, wave, workspace);
            perturbed.surface_emission.as_mut().unwrap()[wave] = original;
            result.surface_emission.as_mut().unwrap()[wave] = (plus - minus) / (2.0 * h);
        }
    }
    Ok(result)
}

fn map_layer_adjoint(
    geometry: &Geometry,
    mode: SourceMode,
    atmosphere: &AtmosphereBatch,
    inputs: &LayerInputs,
    layer: &LayerAdjoints,
    execution: ExecutionPolicy,
) -> Result<AtmosphereAdjoints, TwoStreamError> {
    let n = geometry.num_layers();
    let nw = atmosphere.num_wavelengths;
    if execution == ExecutionPolicy::Serial || nw < 2 {
        return Ok(map_layer_adjoint_range(
            geometry, mode, atmosphere, inputs, layer, 0, nw, 0, nw,
        ));
    }

    use rayon::prelude::*;
    let threads = rayon::current_num_threads().min(nw);
    if threads == 1 {
        return Ok(map_layer_adjoint_range(
            geometry, mode, atmosphere, inputs, layer, 0, nw, 0, nw,
        ));
    }
    let tile_size = nw.div_ceil(threads);
    let parts: Vec<_> = (0..nw)
        .step_by(tile_size)
        .map(|start| (start, (nw - start).min(tile_size)))
        .collect::<Vec<_>>()
        .into_par_iter()
        .map(|(start, len)| {
            (
                start,
                map_layer_adjoint_range(
                    geometry, mode, atmosphere, inputs, layer, start, len, start, nw,
                ),
            )
        })
        .collect();
    let mut out = zero_atmosphere_adjoint(mode, n, nw);
    for (start, part) in parts {
        merge_atmosphere_adjoint(n, start, nw, &part, &mut out);
    }
    Ok(out)
}

fn empty_atmosphere_adjoint(
    mode: SourceMode,
    n: usize,
    nw: usize,
    layer: &LayerAdjoints,
    layer_start: usize,
) -> AtmosphereAdjoints {
    let mut output = zero_atmosphere_adjoint(mode, n, nw);
    output
        .surface_albedo
        .copy_from_slice(&layer.surface_albedo[layer_start..layer_start + nw]);
    if let (Some(output_emission), Some(layer_emission)) =
        (&mut output.surface_emission, &layer.surface_emission)
    {
        output_emission.copy_from_slice(&layer_emission[layer_start..layer_start + nw]);
    }
    output
}

fn zero_atmosphere_adjoint(mode: SourceMode, n: usize, nw: usize) -> AtmosphereAdjoints {
    let levels = (n + 1) * nw;
    AtmosphereAdjoints {
        extinction: vec![0.0; levels],
        single_scatter_albedo: vec![0.0; levels],
        first_legendre: vec![0.0; levels],
        emission: (mode == SourceMode::Thermal).then(|| vec![0.0; levels]),
        surface_albedo: vec![0.0; nw],
        surface_emission: (mode == SourceMode::Thermal).then(|| vec![0.0; nw]),
    }
}

#[allow(clippy::too_many_arguments)]
fn map_layer_adjoint_range(
    geometry: &Geometry,
    mode: SourceMode,
    atmosphere: &AtmosphereBatch,
    inputs: &LayerInputs,
    layer: &LayerAdjoints,
    start: usize,
    len: usize,
    layer_start: usize,
    layer_nw: usize,
) -> AtmosphereAdjoints {
    let n = geometry.num_layers();
    let input_nw = atmosphere.num_wavelengths;
    let mut out = empty_atmosphere_adjoint(mode, n, len, layer, layer_start);
    let mut d_od = vec![0.0; n * len];
    for l in 0..n {
        d_od[l * len..(l + 1) * len].copy_from_slice(
            &layer.optical_depth[l * layer_nw + layer_start..l * layer_nw + layer_start + len],
        );
    }

    if mode == SourceMode::Solar {
        let d_transmission = layer.transmission.as_ref().unwrap();
        let d_secant = layer.average_secant.as_ref().unwrap();
        let transmission = inputs.transmission.as_ref().unwrap();
        let average_secant = inputs.average_secant.as_ref().unwrap();
        let mut d_slant = vec![0.0; n + 1];
        for output_wave in 0..len {
            d_slant.fill(0.0);
            let wave = start + output_wave;
            let layer_wave = layer_start + output_wave;
            for boundary in 1..=n {
                d_slant[boundary] += d_transmission[boundary * layer_nw + layer_wave]
                    * transmission[boundary * input_nw + wave];
            }
            for l in 0..n {
                let li = l * input_nw + wave;
                let layer_li = l * layer_nw + layer_wave;
                let output_li = l * len + output_wave;
                let inv_od = if inputs.optical_depth[li] > 0.0 {
                    1.0 / inputs.optical_depth[li]
                } else {
                    0.0
                };
                let weight = d_secant[layer_li] * inv_od;
                d_slant[l] += weight;
                d_slant[l + 1] -= weight;
                d_od[output_li] -= d_secant[layer_li] * average_secant[li] * inv_od;
            }
            for l in 0..n {
                for (boundary, d_slant) in d_slant.iter().enumerate().skip(1) {
                    d_od[l * len + output_wave] -=
                        *d_slant * geometry.chapman_factors[(boundary - 1) * n + l];
                }
            }
        }
    }
    for l in 0..n {
        let dz = geometry.layer_thickness[l];
        for output_wave in 0..len {
            let wave = start + output_wave;
            let li = l * input_nw + wave;
            let layer_li = l * layer_nw + layer_start + output_wave;
            let output_li = l * len + output_wave;

            if mode == SourceMode::Thermal {
                let emission = atmosphere.emission.as_ref().unwrap();
                let top = emission[l * input_nw + wave];
                let bottom = emission[(l + 1) * input_nw + wave];
                let d_b0 = layer.thermal_b0.as_ref().unwrap()[layer_li];
                let d_b1 = layer.thermal_b1.as_ref().unwrap()[layer_li];
                let (d_top, d_bottom, d_layer_od) = super::thermal_profile_slope_adjoint(
                    top,
                    bottom,
                    inputs.optical_depth[li],
                    inputs.thermal_b1.as_ref().unwrap()[li],
                    d_b1,
                );
                out.emission.as_mut().unwrap()[l * len + output_wave] += d_b0 + d_top;
                out.emission.as_mut().unwrap()[(l + 1) * len + output_wave] += d_bottom;
                d_od[output_li] += d_layer_od;
            }

            let od_density = inputs.optical_depth[li] / dz;
            let scattering_density = inputs.single_scatter_albedo[li] * od_density;
            let inv_od_density = if od_density > 0.0 {
                1.0 / od_density
            } else {
                0.0
            };
            let inv_scattering_density = if scattering_density > 0.0 {
                1.0 / scattering_density
            } else {
                0.0
            };
            for level in [l, l + 1] {
                let i = level * input_nw + wave;
                let output_i = level * len + output_wave;
                out.extinction[output_i] += 0.5 * d_od[output_li] * dz;
                out.single_scatter_albedo[output_i] += 0.5
                    * layer.single_scatter_albedo[layer_li]
                    * atmosphere.extinction[i]
                    * inv_od_density;
                out.extinction[output_i] += 0.5
                    * layer.single_scatter_albedo[layer_li]
                    * (atmosphere.single_scatter_albedo[i] - inputs.single_scatter_albedo[li])
                    * inv_od_density;

                out.first_legendre[output_i] += 0.5
                    * layer.first_legendre[layer_li]
                    * atmosphere.single_scatter_albedo[i]
                    * atmosphere.extinction[i]
                    * inv_scattering_density;
                out.extinction[output_i] += 0.5
                    * layer.first_legendre[layer_li]
                    * atmosphere.single_scatter_albedo[i]
                    * (atmosphere.first_legendre[i] - inputs.first_legendre[li])
                    * inv_scattering_density;
                out.single_scatter_albedo[output_i] += 0.5
                    * layer.first_legendre[layer_li]
                    * atmosphere.extinction[i]
                    * (atmosphere.first_legendre[i] - inputs.first_legendre[li])
                    * inv_scattering_density;
            }
        }
    }
    out
}

fn merge_atmosphere_adjoint(
    n: usize,
    start: usize,
    output_nw: usize,
    part: &AtmosphereAdjoints,
    output: &mut AtmosphereAdjoints,
) {
    let part_nw = part.surface_albedo.len();
    for level in 0..=n {
        let output_range = level * output_nw + start..level * output_nw + start + part_nw;
        let part_range = level * part_nw..(level + 1) * part_nw;
        output.extinction[output_range.clone()]
            .copy_from_slice(&part.extinction[part_range.clone()]);
        output.single_scatter_albedo[output_range.clone()]
            .copy_from_slice(&part.single_scatter_albedo[part_range.clone()]);
        output.first_legendre[output_range.clone()]
            .copy_from_slice(&part.first_legendre[part_range.clone()]);
        if let (Some(output_emission), Some(part_emission)) = (&mut output.emission, &part.emission)
        {
            output_emission[output_range].copy_from_slice(&part_emission[part_range]);
        }
    }
    output.surface_albedo[start..start + part_nw].copy_from_slice(&part.surface_albedo);
    if let (Some(output_emission), Some(part_emission)) =
        (&mut output.surface_emission, &part.surface_emission)
    {
        output_emission[start..start + part_nw].copy_from_slice(part_emission);
    }
}

fn validate_geometry(geometry: &Geometry) -> Result<(), TwoStreamError> {
    let n = geometry.num_layers();
    if n == 0 {
        return Err(TwoStreamError::invalid(
            "geometry must contain at least one layer",
        ));
    }
    if geometry.chapman_factors.len() != n * n {
        return Err(TwoStreamError::invalid(
            "Chapman matrix must have shape [layer, layer]",
        ));
    }
    if geometry
        .layer_thickness
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(TwoStreamError::invalid(
            "layer thicknesses must be finite and positive",
        ));
    }
    if !(-1.0..=1.0).contains(&geometry.solar_cosine)
        || !(0.0..=1.0).contains(&geometry.quadrature_cosine)
        || geometry.quadrature_cosine == 0.0
    {
        return Err(TwoStreamError::invalid("invalid angular cosine"));
    }
    Ok(())
}

fn validate_spherical_geometry(
    solver_geometry: &Geometry,
    mode: SourceMode,
    spherical: &SphericalGeometry,
) -> Result<(), TwoStreamError> {
    if spherical.columns.is_empty() || spherical.columns.len() != spherical.sza_grid.len() {
        return Err(TwoStreamError::invalid(
            "spherical SZA columns have inconsistent shapes",
        ));
    }
    for (index, column) in spherical.columns.iter().enumerate() {
        validate_geometry(column)?;
        if column.layer_thickness != solver_geometry.layer_thickness
            || column.quadrature_cosine != solver_geometry.quadrature_cosine
        {
            return Err(TwoStreamError::invalid(
                "spherical columns do not match the solver vertical grid",
            ));
        }
        if mode == SourceMode::Solar && column.solar_cosine != spherical.sza_grid[index] {
            return Err(TwoStreamError::invalid(
                "spherical column cosine does not match the SZA grid",
            ));
        }
    }
    if spherical.rays.num_levels != solver_geometry.num_layers() + 1 {
        return Err(TwoStreamError::invalid(
            "spherical rays do not match the solver vertical grid",
        ));
    }
    Ok(())
}

fn validate_atmosphere(
    geometry: &Geometry,
    mode: SourceMode,
    atmosphere: &AtmosphereBatch,
) -> Result<(), TwoStreamError> {
    let nw = atmosphere.num_wavelengths;
    let level_len = (geometry.num_layers() + 1) * nw;
    if nw == 0 {
        return Err(TwoStreamError::invalid(
            "atmosphere must contain wavelengths",
        ));
    }
    for (name, len) in [
        ("extinction", atmosphere.extinction.len()),
        (
            "single_scatter_albedo",
            atmosphere.single_scatter_albedo.len(),
        ),
        ("first_legendre", atmosphere.first_legendre.len()),
    ] {
        if len != level_len {
            return Err(TwoStreamError::invalid(format!(
                "{name} has the wrong shape"
            )));
        }
    }
    if atmosphere.surface_albedo.len() != nw {
        return Err(TwoStreamError::invalid(
            "surface_albedo has the wrong shape",
        ));
    }
    match mode {
        SourceMode::Solar if atmosphere.solar_irradiance.as_ref().map(Vec::len) != Some(nw) => {
            return Err(TwoStreamError::invalid("solar mode requires irradiance"));
        }
        SourceMode::Thermal
            if atmosphere.emission.as_ref().map(Vec::len) != Some(level_len)
                || atmosphere.surface_emission.as_ref().map(Vec::len) != Some(nw) =>
        {
            return Err(TwoStreamError::invalid(
                "thermal mode requires level and surface emission",
            ));
        }
        _ => {}
    }
    Ok(())
}

fn validate_layer_inputs(
    geometry: &Geometry,
    mode: SourceMode,
    inputs: &LayerInputs,
    views: &[View],
) -> Result<(), TwoStreamError> {
    if inputs.num_layers != geometry.num_layers() || inputs.num_wavelengths == 0 {
        return Err(TwoStreamError::invalid(
            "layer input dimensions do not match geometry",
        ));
    }
    let len = inputs.num_layers * inputs.num_wavelengths;
    for (name, actual) in [
        ("optical_depth", inputs.optical_depth.len()),
        ("single_scatter_albedo", inputs.single_scatter_albedo.len()),
        ("first_legendre", inputs.first_legendre.len()),
    ] {
        if actual != len {
            return Err(TwoStreamError::invalid(format!(
                "{name} has the wrong shape"
            )));
        }
    }
    if inputs.surface_albedo.len() != inputs.num_wavelengths {
        return Err(TwoStreamError::invalid(
            "surface_albedo has the wrong shape",
        ));
    }
    match mode {
        SourceMode::Solar
            if inputs.transmission.as_ref().map(Vec::len)
                != Some((inputs.num_layers + 1) * inputs.num_wavelengths)
                || inputs.average_secant.as_ref().map(Vec::len) != Some(len) =>
        {
            return Err(TwoStreamError::invalid("solar layer inputs are incomplete"));
        }
        SourceMode::Thermal
            if inputs.thermal_b0.as_ref().map(Vec::len) != Some(len)
                || inputs.thermal_b1.as_ref().map(Vec::len) != Some(len)
                || inputs.surface_emission.as_ref().map(Vec::len)
                    != Some(inputs.num_wavelengths) =>
        {
            return Err(TwoStreamError::invalid(
                "thermal layer inputs are incomplete",
            ));
        }
        _ => {}
    }
    validate_views(views)
}

fn validate_views(views: &[View]) -> Result<(), TwoStreamError> {
    if views.is_empty()
        || views
            .iter()
            .any(|view| !view.cosine.is_finite() || view.cosine <= 0.0 || view.cosine > 1.0)
    {
        return Err(TwoStreamError::invalid(
            "views must be finite and upwelling",
        ));
    }
    Ok(())
}

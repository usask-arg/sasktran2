use crate::atmosphere::*;
use crate::constituent::traits::*;
use crate::prelude::*;

pub struct Manual {
    pub extinction: Array2<f64>,
    pub ssa: Array2<f64>,
    pub legendre_moments: Option<Array3<f64>>,
}

impl Manual {
    pub fn new(
        extinction: Array2<f64>,
        ssa: Array2<f64>,
        legendre_moments: Option<Array3<f64>>,
    ) -> Self {
        Self {
            extinction,
            ssa,
            legendre_moments,
        }
    }
}

impl Constituent for Manual {
    fn add_to_atmosphere(&self, storage: &mut impl AtmosphereStorageAccess) -> Result<()> {
        let (_inputs, outputs) = storage.split_inputs_outputs();

        let mut outputs = outputs.mut_view();
        let (num_output_locations, num_wavel) = outputs.total_extinction.dim();
        let (num_input_locations, input_num_wavel) = self.extinction.dim();

        anyhow::ensure!(
            self.ssa.dim() == self.extinction.dim(),
            "manual extinction and single-scatter albedo shapes must match"
        );
        anyhow::ensure!(
            num_input_locations > 0,
            "manual constituent must contain at least one spatial location"
        );
        anyhow::ensure!(
            input_num_wavel == num_wavel,
            "manual wavelength dimension ({input_num_wavel}) does not match atmosphere ({num_wavel})"
        );
        anyhow::ensure!(
            num_output_locations % num_input_locations == 0,
            "manual spatial dimension ({num_input_locations}) does not match atmosphere ({num_output_locations})"
        );

        if let Some(legendre_moments) = &self.legendre_moments {
            anyhow::ensure!(
                legendre_moments.dim()
                    == (outputs.legendre.dim().0, num_input_locations, num_wavel,),
                "manual Legendre moment shape {:?} does not match atmosphere ({}, {}, {})",
                legendre_moments.dim(),
                outputs.legendre.dim().0,
                num_input_locations,
                num_wavel
            );
        }

        let thread_pool = crate::threading::thread_pool()?;

        thread_pool.install(|| {
            Zip::from(outputs.total_extinction.columns_mut())
                .and(outputs.ssa.columns_mut())
                .and(self.extinction.columns())
                .and(self.ssa.columns())
                .par_for_each(|k, w, k_input, w_input| {
                    Zip::indexed(k).and(w).for_each(|location, k, w| {
                        let input_location = location % num_input_locations;
                        let k_self = k_input[input_location];
                        let w_self = w_input[input_location];
                        *k += k_self;
                        *w += k_self * w_self;
                    })
                })
        });

        if let Some(legendre_moments) = &self.legendre_moments {
            thread_pool.install(|| {
                Zip::from(outputs.legendre.axis_iter_mut(Axis(2)))
                    .and(self.extinction.columns())
                    .and(self.ssa.columns())
                    .and(legendre_moments.axis_iter(Axis(2)))
                    .par_for_each(|mut bl_assign, k_input, w_input, bl_input| {
                        Zip::indexed(bl_assign.columns_mut()).for_each(|location, bl_assign| {
                            let input_location = location % num_input_locations;
                            let scat_ext = k_input[input_location] * w_input[input_location];
                            let bl = bl_input.column(input_location);

                            Zip::from(bl_assign).and(bl).for_each(|bl_assign, bl| {
                                *bl_assign += scat_ext * *bl;
                            });
                        })
                    })
            })
        }

        Ok(())
    }

    fn register_derivatives(
        &self,
        _storage: &mut impl AtmosphereStorageAccess,
        _constituent_name: &str,
    ) -> Result<()> {
        // No derivatives for this one?
        Ok(())
    }
}

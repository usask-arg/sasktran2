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
        // Ensure that the arrays have compatible shapes
        let (_inputs, outputs) = storage.split_inputs_outputs();

        let mut outputs = outputs.mut_view();

        let thread_pool = crate::threading::thread_pool()?;

        thread_pool.install(|| {
            Zip::from(outputs.total_extinction.columns_mut())
                .and(outputs.ssa.columns_mut())
                .and(self.extinction.columns())
                .and(self.ssa.columns())
                .par_for_each(|k, w, k_self, w_self| {
                    Zip::from(k)
                        .and(w)
                        .and(k_self)
                        .and(w_self)
                        .for_each(|k, w, k_self, w_self| {
                            *k += *k_self;
                            *w += *k_self * *w_self;
                        })
                })
        });

        if let Some(legendre_moments) = &self.legendre_moments {
            thread_pool.install(|| {
                Zip::from(outputs.legendre.axis_iter_mut(Axis(2)))
                    .and(self.extinction.columns())
                    .and(self.ssa.columns())
                    .and(legendre_moments.axis_iter(Axis(2)))
                    .par_for_each(|mut bl_assign, k, w, bl| {
                        Zip::from(bl_assign.columns_mut())
                            .and(k)
                            .and(w)
                            .and(bl.columns())
                            .for_each(|bl_assign, k, w, bl| {
                                let scat_ext = *k * *w;

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

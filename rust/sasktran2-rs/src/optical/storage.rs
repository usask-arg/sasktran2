use ndarray::{IntoDimension, ShapeBuilder};

use crate::prelude::*;

pub struct OpticalQuantities {
    pub cross_section: Array2<f64>,
    pub ssa: Array2<f64>,
    pub legendre: Option<Array3<f64>>,
    pub fortran_ordering: bool,
}

impl Default for OpticalQuantities {
    fn default() -> Self {
        Self {
            cross_section: Array2::zeros((0, 0)),
            ssa: Array2::zeros((0, 0)),
            legendre: None,
            fortran_ordering: false,
        }
    }
}

impl OpticalQuantities {
    pub fn new(num_geometry: usize, num_wavelengths: usize, fortran_ordering: bool) -> Self {
        let mut instance = Self {
            cross_section: Array2::zeros((0, 0)),
            ssa: Array2::zeros((0, 0)),
            legendre: None,
            fortran_ordering,
        };
        instance.resize(num_geometry, num_wavelengths);

        instance
    }

    pub fn resize(&mut self, num_geometry: usize, num_wavelengths: usize) -> &mut Self {
        let dims = match self.fortran_ordering {
            true => (num_geometry, num_wavelengths).f(),
            false => (num_geometry, num_wavelengths).into_dimension().into(),
        };

        if self.cross_section.dim() != (num_geometry, num_wavelengths) {
            self.cross_section = Array2::zeros(dims);
        }
        if self.ssa.dim() != (num_geometry, num_wavelengths) {
            self.ssa = Array2::zeros(dims);
        }

        self
    }

    // Here num_legendre is the STACKED dimension (i.e., num_orders * stokes_factor)
    pub fn with_scatterer(&mut self, num_legendre: usize, _num_stokes: usize) -> &mut Self {
        if self.legendre.is_none() {
            self.legendre = Some(Array3::zeros((
                self.cross_section.dim().0,
                self.cross_section.dim().1,
                num_legendre,
            )));
        }

        self
    }

    pub fn set_zero(&mut self) {
        self.cross_section.fill(0.0);
        self.ssa.fill(0.0);
        if let Some(legendre) = &mut self.legendre {
            legendre.fill(0.0);
        }
    }

    pub fn mut_split(&mut self) -> (&mut Array2<f64>, &mut Array2<f64>, Option<&mut Array3<f64>>) {
        (
            &mut self.cross_section,
            &mut self.ssa,
            self.legendre.as_mut(),
        )
    }
}

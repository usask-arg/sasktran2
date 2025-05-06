use crate::prelude::*;

pub struct OpticalQuantities {
    pub cross_section: Array2<f64>,
    pub ssa: Array2<f64>,
    pub legendre: Option<Array3<f64>>,
}

impl Default for OpticalQuantities {
    fn default() -> Self {
        Self {
            cross_section: Array2::zeros((0, 0)),
            ssa: Array2::zeros((0, 0)),
            legendre: None,
        }
    }
}

impl OpticalQuantities {
    pub fn new(num_geometry: usize, num_wavelengths: usize) -> Self {
        let mut default = Self::default();
        default.resize(num_geometry, num_wavelengths);

        default
    }

    pub fn resize(&mut self, num_geometry: usize, num_wavelengths: usize) -> &mut Self {
        if self.cross_section.dim() != (num_geometry, num_wavelengths) {
            self.cross_section = Array2::zeros((num_geometry, num_wavelengths));
        }
        if self.ssa.dim() != (num_geometry, num_wavelengths) {
            self.ssa = Array2::zeros((num_geometry, num_wavelengths));
        }

        self
    }

    pub fn with_scatterer(&mut self, num_legendre: usize, num_stokes: usize) -> &mut Self {
        let num_params = match num_stokes {
            1 => 1,
            3 => 4,
            4 => 6,
            _ => panic!("Invalid number of Stokes parameters"),
        };

        if self.legendre.is_none() {
            self.legendre = Some(Array3::zeros((
                num_params * num_legendre,
                self.cross_section.dim().0,
                self.cross_section.dim().1,
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

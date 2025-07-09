use crate::interpolation::linear::*;
use crate::interpolation::*;
use crate::prelude::*;

pub struct Grid1D {
    pub x: Array1<f64>,
    is_uniform: bool,
    start: Option<f64>,
    dx: Option<f64>,
}

pub struct Grid1DView<'a> {
    pub x: &'a [f64],
    is_uniform: bool,
    start: Option<f64>,
    dx: Option<f64>,
}

impl Grid1D {
    pub fn new(x: Array1<f64>) -> Self {
        // Check if the distance between successive elements is uniform
        let first_diff = x[1] - x[0];
        let is_uniform = x
            .windows(2)
            .into_iter()
            .all(|w| (w[1] - w[0] - first_diff).abs() < 1e-10);

        let start = x.first().copied();
        let dx = if is_uniform { Some(x[1] - x[0]) } else { None };
        Self {
            x,
            is_uniform,
            start,
            dx,
        }
    }

    pub fn view(&self) -> Grid1DView<'_> {
        Grid1DView {
            x: self.x.as_slice().unwrap(),
            is_uniform: self.is_uniform,
            start: self.start,
            dx: self.dx,
        }
    }
}

impl<'a> Grid1DView<'a> {
    pub fn new(x: &'a [f64]) -> Self {
        // Check if the distance between successive elements is uniform
        let first_diff = x[1] - x[0];
        let is_uniform = x
            .windows(2)
            .all(|w| (w[1] - w[0] - first_diff).abs() < 1e-10);

        let start = x.first().copied();
        let dx = if is_uniform { Some(x[1] - x[0]) } else { None };
        Self {
            x,
            is_uniform,
            start,
            dx,
        }
    }

    pub fn slice(&self, start: usize, end: usize) -> Self {
        let x = &self.x[start..end];
        let is_uniform = self.is_uniform;
        let start = x.first().copied();
        let dx = if is_uniform {
            Some(self.x[1] - self.x[0])
        } else {
            None
        };
        Self {
            x,
            is_uniform,
            start,
            dx,
        }
    }

    pub fn lower_bound(&self, to_x: f64) -> usize {
        if self.is_uniform {
            if to_x < self.x[0] {
                return 0;
            } else if to_x > self.x[self.x.len() - 1] {
                return self.x.len() - 1;
            }
            ((to_x - self.start.unwrap()) / self.dx.unwrap()).floor() as usize
        } else {
            self.x
                .binary_search_by(|x| x.partial_cmp(&to_x).unwrap())
                .unwrap_or_else(|i| i)
        }
    }
}

impl Interp1Weights for Grid1D {
    fn interp1_weights(
        &self,
        to_x: f64,
        out_of_bounds_mode: OutOfBoundsMode,
    ) -> [(usize, f64, f64); 2] {
        if self.is_uniform {
            if to_x < *self.x.first().unwrap() {
                match out_of_bounds_mode {
                    OutOfBoundsMode::Zero => return [(0, 0.0, 0.0), (0, 0.0, 0.0)],
                    OutOfBoundsMode::Extend => return [(0, 1.0, 0.0), (0, 0.0, 0.0)],
                }
            } else if to_x > *self.x.last().unwrap() {
                match out_of_bounds_mode {
                    OutOfBoundsMode::Zero => return [(0, 0.0, 0.0), (0, 0.0, 0.0)],
                    OutOfBoundsMode::Extend => {
                        return [(self.x.len() - 1, 1.0, 0.0), (self.x.len() - 1, 0.0, 0.0)];
                    }
                }
            }
            let idx = ((to_x - self.start.unwrap()) / self.dx.unwrap()).floor() as usize;
            let w =
                (to_x - self.start.unwrap() - (idx as f64) * self.dx.unwrap()) / self.dx.unwrap();
            let d_w = 1.0 / self.dx.unwrap();
            [(idx, 1.0 - w, -d_w), (idx + 1, w, d_w)]
        } else {
            self.x.interp1_weights(to_x, out_of_bounds_mode)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interp1_weights_uniform_grid() {
        let grid = Grid1D::new(array![0.0, 1.0, 2.0, 3.0]);

        let weights = grid.interp1_weights(1.5, OutOfBoundsMode::Zero);
        assert_eq!(weights, [(1, 0.5, -1.0), (2, 0.5, 1.0)]);
    }

    #[test]
    fn test_interp1_weights_non_uniform_grid() {
        let grid = Grid1D::new(array![0.0, 1.0, 2.5, 4.0]);

        let weights = grid.interp1_weights(2.0, OutOfBoundsMode::Zero);
        let expected_weights = [
            (1, 0.3333333333333333, -0.6666666666666666),
            (2, 0.6666666666666666, 0.6666666666666666),
        ];
        for (actual, expected) in weights.iter().zip(expected_weights.iter()) {
            assert_eq!(actual.0, expected.0);
            assert!((actual.1 - expected.1).abs() < 1e-10);
            assert!((actual.2 - expected.2).abs() < 1e-10);
        }
    }

    #[test]
    fn test_interp1_weights_out_of_bounds() {
        let grid = Grid1D::new(array![0.0, 1.0, 2.0]);

        let weights = grid.interp1_weights(-1.0, OutOfBoundsMode::Extend);
        assert_eq!(weights, [(0, 1.0, 0.0), (0, 0.0, 0.0)]);

        let weights = grid.interp1_weights(3.0, OutOfBoundsMode::Extend);
        assert_eq!(weights, [(2, 1.0, 0.0), (2, 0.0, 0.0)]);
    }
}

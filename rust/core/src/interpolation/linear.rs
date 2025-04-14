use super::OutOfBoundsMode;
use ndarray::*;

pub fn linear_interpolating_matrix<S1, S2>(
    from_grid: &ArrayBase<S1, Ix1>,
    to_grid: &ArrayBase<S2, Ix1>,
    out_of_bounds_mode: OutOfBoundsMode,
) -> Array2<f64>
where
    S1: Data<Elem = f64>,
    S2: Data<Elem = f64>,
{
    let mut result: Array2<f64> = Array2::zeros((to_grid.len(), from_grid.len()));

    Zip::from(result.rows_mut())
        .and(to_grid)
        .for_each(|mut row, to_val| {
            if to_val < from_grid.first().unwrap() {
                match out_of_bounds_mode {
                    OutOfBoundsMode::Zero => {} // Do nothing
                    OutOfBoundsMode::Extend => row.first_mut().unwrap().assign_elem(1.0),
                }
            } else if to_val > from_grid.last().unwrap() {
                match out_of_bounds_mode {
                    OutOfBoundsMode::Zero => {} // Do nothing
                    OutOfBoundsMode::Extend => row.last_mut().unwrap().assign_elem(1.0),
                }
            } else {
                let idx = from_grid.iter().position(|&x| x >= *to_val).unwrap();

                let w = (*to_val - from_grid[idx - 1]) / (from_grid[idx] - from_grid[idx - 1]);
                row[idx - 1] = 1.0 - w;
                row[idx] = w;
            }
        });
    result
}

pub struct Grid {
    pub x: Array1<f64>,
    is_uniform: bool,
    start: Option<f64>,
    dx: Option<f64>,
}

impl Grid {
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
}

pub trait Interp1<S: Data<Elem = f64>> {
    fn interp1(
        &self,
        x_grid: &ArrayBase<S, Ix1>,
        to_x: f64,
        out_of_bounds_mode: OutOfBoundsMode,
    ) -> f64;
}
pub trait Interp1Weights {
    fn interp1_weights(
        &self,
        to_x: f64,
        out_of_bounds_mode: OutOfBoundsMode,
    ) -> [(usize, f64, f64); 2];
}

impl<S> Interp1<S> for ArrayBase<S, Ix1>
where
    S: Data<Elem = f64>,
{
    fn interp1(
        &self,
        x_grid: &ArrayBase<S, Ix1>,
        to_x: f64,
        out_of_bounds_mode: OutOfBoundsMode,
    ) -> f64 {
        if to_x < *x_grid.first().unwrap() {
            match out_of_bounds_mode {
                OutOfBoundsMode::Zero => return 0.0,
                OutOfBoundsMode::Extend => return *self.first().unwrap(),
            }
        } else if to_x > *x_grid.last().unwrap() {
            match out_of_bounds_mode {
                OutOfBoundsMode::Zero => return 0.0,
                OutOfBoundsMode::Extend => return *self.last().unwrap(),
            }
        }
        // Binary search for the first element >= to_x
        let idx = match x_grid
            .as_slice()
            .unwrap()
            .binary_search_by(|x| x.partial_cmp(&to_x).unwrap())
        {
            Ok(i) => i,  // Exact match found
            Err(i) => i, // Not found, i is the index where `to_x` would be inserted
        };

        if idx == 0 {
            return *self.first().unwrap();
        } else if idx == x_grid.len() {
            return *self.last().unwrap();
        }

        let w = (to_x - x_grid[idx - 1]) / (x_grid[idx] - x_grid[idx - 1]);
        let y1 = *self.get(idx - 1).unwrap();
        let y2 = *self.get(idx).unwrap();
        y1 + (y2 - y1) * w
    }
}
impl Interp1Weights for Array1<f64> {
    fn interp1_weights(
        &self,
        to_x: f64,
        out_of_bounds_mode: OutOfBoundsMode,
    ) -> [(usize, f64, f64); 2] {
        if to_x < *self.first().unwrap() {
            match out_of_bounds_mode {
                OutOfBoundsMode::Zero => return [(0, 0.0, 0.0), (0, 0.0, 0.0)],
                OutOfBoundsMode::Extend => return [(0, 1.0, 0.0), (0, 0.0, 0.0)],
            }
        } else if to_x > *self.last().unwrap() {
            match out_of_bounds_mode {
                OutOfBoundsMode::Zero => return [(0, 0.0, 0.0), (0, 0.0, 0.0)],
                OutOfBoundsMode::Extend => {
                    return [(self.len() - 1, 1.0, 0.0), (self.len() - 1, 0.0, 0.0)];
                }
            }
        }
        // Binary search for the first element >= to_x
        let idx = match self
            .as_slice()
            .unwrap()
            .binary_search_by(|x| x.partial_cmp(&to_x).unwrap())
        {
            Ok(i) => i,  // Exact match found
            Err(i) => i, // Not found, i is the index where `to_x` would be inserted
        };

        let w = (to_x - self[idx - 1]) / (self[idx] - self[idx - 1]);

        // Also consider dw/dx = 1/(x[idx] - x[idx-1])
        let d_w = 1.0 / (self[idx] - self[idx - 1]);

        [(idx - 1, 1.0 - w, -d_w), (idx, w, d_w)]
    }
}

impl Interp1Weights for Grid {
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
    fn test_linear_interpolating_matrix() {
        let from_grid = array![0.0, 1.0, 2.0];
        let to_grid = array![0.5, 1.5];

        let result = linear_interpolating_matrix(&from_grid, &to_grid, OutOfBoundsMode::Zero);

        println!("{:?}", result);
    }

    #[test]
    fn test_linear_interpolating_matrix_with_extend() {
        let from_grid = array![0.0, 1.0, 2.0];
        let to_grid = array![-0.5, 0.5, 1.5, 2.5];

        let result = linear_interpolating_matrix(&from_grid, &to_grid, OutOfBoundsMode::Extend);

        assert_eq!(
            result,
            array![
                [1.0, 0.0, 0.0],
                [0.5, 0.5, 0.0],
                [0.0, 0.5, 0.5],
                [0.0, 0.0, 1.0]
            ]
        );
    }

    #[test]
    fn test_interp1_with_zero_out_of_bounds() {
        let x_grid = array![0.0, 1.0, 2.0];
        let y_values = array![10.0, 20.0, 30.0];

        let result = y_values.interp1(&x_grid, -1.0, OutOfBoundsMode::Zero);
        assert_eq!(result, 0.0);

        let result = y_values.interp1(&x_grid, 3.0, OutOfBoundsMode::Zero);
        assert_eq!(result, 0.0);
    }

    #[test]
    fn test_interp1_with_extend_out_of_bounds() {
        let x_grid = array![0.0, 1.0, 2.0];
        let y_values = array![10.0, 20.0, 30.0];

        let result = y_values.interp1(&x_grid, -1.0, OutOfBoundsMode::Extend);
        assert_eq!(result, 10.0);

        let result = y_values.interp1(&x_grid, 3.0, OutOfBoundsMode::Extend);
        assert_eq!(result, 30.0);
    }

    #[test]
    fn test_interp1_weights_uniform_grid() {
        let grid = Grid::new(array![0.0, 1.0, 2.0, 3.0]);

        let weights = grid.interp1_weights(1.5, OutOfBoundsMode::Zero);
        assert_eq!(weights, [(1, 0.5, -1.0), (2, 0.5, 1.0)]);
    }

    #[test]
    fn test_interp1_weights_non_uniform_grid() {
        let grid = Grid::new(array![0.0, 1.0, 2.5, 4.0]);

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
        let grid = Grid::new(array![0.0, 1.0, 2.0]);

        let weights = grid.interp1_weights(-1.0, OutOfBoundsMode::Extend);
        assert_eq!(weights, [(0, 1.0, 0.0), (0, 0.0, 0.0)]);

        let weights = grid.interp1_weights(3.0, OutOfBoundsMode::Extend);
        assert_eq!(weights, [(2, 1.0, 0.0), (2, 0.0, 0.0)]);
    }
}

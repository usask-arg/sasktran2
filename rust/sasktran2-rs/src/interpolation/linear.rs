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
            let weights = from_grid.interp1_weights(*to_val, out_of_bounds_mode);
            row[weights[0].0] += weights[0].1;
            row[weights[1].0] += weights[1].1;
        });
    result
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
impl<S> Interp1Weights for ArrayBase<S, Ix1>
where
    S: Data<Elem = f64>,
{
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

        if idx == 0 {
            return [(0, 1.0, 0.0), (0, 0.0, 0.0)];
        } else if idx == self.len() {
            return [(self.len() - 1, 1.0, 0.0), (self.len() - 1, 0.0, 0.0)];
        }

        let w = (to_x - self[idx - 1]) / (self[idx] - self[idx - 1]);

        // Also consider dw/dx = 1/(x[idx] - x[idx-1])
        let d_w = 1.0 / (self[idx] - self[idx - 1]);

        [(idx - 1, 1.0 - w, -d_w), (idx, w, d_w)]
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

        println!("{result:?}");
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
    fn test_linear_interpolating_matrix_oob() {
        let from_grid = array![0.0, 0.000004];
        let to_grid = array![-0.5, 0.0, 1.5, 2.5];

        let _result = linear_interpolating_matrix(&from_grid, &to_grid, OutOfBoundsMode::Extend);
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
}

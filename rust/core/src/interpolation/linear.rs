use ndarray::*;
use super::OutOfBoundsMode;

pub fn linear_interpolating_matrix<S>(
    from_grid: &ArrayBase<S, Ix1>,
    to_grid: &ArrayBase<S, Ix1>,
    out_of_bounds_mode: OutOfBoundsMode,
) -> Array2<f64>
where
    S: Data<Elem = f64>,
{
    let mut result: Array2<f64> = Array2::zeros((to_grid.len(), from_grid.len()));

    Zip::from(result.rows_mut())
        .and(to_grid)
        .for_each(|mut row, to_val| {
            if to_val < from_grid.first().unwrap() {
                match out_of_bounds_mode {
                    OutOfBoundsMode::Zero =>  {}, // Do nothing
                    OutOfBoundsMode::Extend => row.first_mut().unwrap().assign_elem(1.0),
                }
            } else if to_val > from_grid.last().unwrap() {
                match out_of_bounds_mode {
                    OutOfBoundsMode::Zero => {}, // Do nothing
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

pub trait Interp1<S: Data<Elem = f64>> {
    fn interp1(
        &self,
        x_grid: &ArrayBase<S, Ix1>,
        to_x: f64,
        out_of_bounds_mode: OutOfBoundsMode,
    ) -> f64;
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
        let idx = x_grid.iter().position(|&x| x >= to_x).unwrap();
        let w = (to_x - x_grid[idx - 1]) / (x_grid[idx] - x_grid[idx - 1]);
        let y1 = *self.get(idx - 1).unwrap();
        let y2 = *self.get(idx).unwrap();
        y1 + (y2 - y1) * w
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
}

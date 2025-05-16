use ndarray::*;

pub struct LegendreAccess<'a, D: Dimension, const NSTOKES: usize> {
    pub a1: ArrayView<'a, f64, D>,
    pub a2: Option<ArrayView<'a, f64, D>>,
    pub a3: Option<ArrayView<'a, f64, D>>,
    pub a4: Option<ArrayView<'a, f64, D>>,
    pub b1: Option<ArrayView<'a, f64, D>>,
    pub b2: Option<ArrayView<'a, f64, D>>,
}

impl<'a, D: Dimension> LegendreAccess<'a, D, 1> {
    pub fn new(legendre: ArrayView<'a, f64, D>) -> Self {
        let a1 = legendre.clone();

        Self {
            a1,
            a3: None,
            a2: None,
            a4: None,
            b1: None,
            b2: None,
        }
    }
}

impl<'a, D: Dimension> LegendreAccess<'a, D, 3> {
    pub fn new(legendre: ArrayView<'a, f64, D>) -> Self {
        Self {
            a1: legendre
                .clone()
                .slice_axis_move(Axis(0), Slice::new(0, None, 4)),
            a2: Some(
                legendre
                    .clone()
                    .slice_axis_move(Axis(0), Slice::new(1, None, 4)),
            ),
            a3: Some(
                legendre
                    .clone()
                    .slice_axis_move(Axis(0), Slice::new(2, None, 4)),
            ),
            a4: None,
            b1: Some(
                legendre
                    .clone()
                    .slice_axis_move(Axis(0), Slice::new(3, None, 4)),
            ),
            b2: None,
        }
    }
}

impl<'a, D: Dimension> LegendreAccess<'a, D, 4> {
    pub fn new(legendre: ArrayView<'a, f64, D>) -> Self {
        Self {
            a1: legendre
                .clone()
                .slice_axis_move(Axis(D::NDIM.unwrap() - 1), Slice::new(0, None, 6)),
            a2: Some(
                legendre
                    .clone()
                    .slice_axis_move(Axis(D::NDIM.unwrap() - 1), Slice::new(1, None, 6)),
            ),
            a3: Some(
                legendre
                    .clone()
                    .slice_axis_move(Axis(D::NDIM.unwrap() - 1), Slice::new(2, None, 6)),
            ),
            a4: Some(
                legendre
                    .clone()
                    .slice_axis_move(Axis(D::NDIM.unwrap() - 1), Slice::new(3, None, 6)),
            ),
            b1: Some(
                legendre
                    .clone()
                    .slice_axis_move(Axis(D::NDIM.unwrap() - 1), Slice::new(4, None, 6)),
            ),
            b2: Some(
                legendre
                    .clone()
                    .slice_axis_move(Axis(D::NDIM.unwrap() - 1), Slice::new(5, None, 6)),
            ),
        }
    }
}

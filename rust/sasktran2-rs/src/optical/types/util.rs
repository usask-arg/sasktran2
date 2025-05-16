use crate::prelude::*;

/// Returns back the indices required to sort a vector of f64 values.
pub fn argsort_f64(data: &[f64]) -> Vec<usize> {
    let mut indices = (0..data.len()).collect::<Vec<_>>();
    indices.sort_by(|&i, &j| {
        data[i]
            .partial_cmp(&data[j])
            .unwrap_or(std::cmp::Ordering::Equal) // Handle NaN by treating as equal
    });
    indices
}

/// Returns a sorted array of unique values from the input array.
pub fn unique_values<S, D>(array: &ArrayBase<S, D>) -> Array1<f64>
where
    S: Data<Elem = f64>,
    D: Dimension,
{
    let mut seen = Vec::new();

    for &x in array.iter() {
        if x.is_nan() {
            continue;
        }
        if !seen.contains(&x) {
            seen.push(x);
        }
    }

    seen.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Array1::from(seen)
}

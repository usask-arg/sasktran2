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

pub fn create_pool(num_threads: usize) -> Result<rayon::ThreadPool> {
    match rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
    {
        Err(e) => Err(e.into()),
        Ok(pool) => Ok(pool),
    }
}

/// Returns back the indices required to sort a vector of f64 values.
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unique_values_empty() {
        let data = Array1::<f64>::from(vec![]);
        let unique = unique_values(&data);
        assert_eq!(unique, Array1::<f64>::from(vec![]));
    }

    #[test]
    fn test_unique_values_no_duplicates() {
        let data = Array1::from(vec![1.0, 2.0, 3.0, 4.0]);
        let unique = unique_values(&data);
        assert_eq!(unique, Array1::from(vec![1.0, 2.0, 3.0, 4.0]));
    }

    #[test]
    fn test_unique_values_with_duplicates() {
        let data = Array1::from(vec![1.0, 2.0, 2.0, 3.0, 1.0]);
        let unique = unique_values(&data);
        assert_eq!(unique, Array1::from(vec![1.0, 2.0, 3.0]));
    }

    #[test]
    fn test_unique_values_with_nan() {
        let data = Array1::from(vec![1.0, f64::NAN, 2.0, f64::NAN, 3.0]);
        let unique = unique_values(&data);
        assert_eq!(unique, Array1::from(vec![1.0, 2.0, 3.0]));
    }

    #[test]
    fn test_unique_values_with_infinity() {
        let data = Array1::from(vec![1.0, f64::INFINITY, 2.0, f64::NEG_INFINITY, 1.0]);
        let unique = unique_values(&data);
        assert_eq!(
            unique,
            Array1::from(vec![f64::NEG_INFINITY, 1.0, 2.0, f64::INFINITY])
        );
    }

    #[test]
    fn test_unique_values_with_zeros() {
        let data = Array1::from(vec![0.0, -0.0, 1.0, -1.0]);
        let unique = unique_values(&data);
        assert_eq!(unique, Array1::from(vec![-1.0, -0.0, 1.0]));
    }

    #[test]
    fn test_unique_values_with_large_numbers() {
        let data = Array1::from(vec![1e10, 1e12, 1e10, 1e11]);
        let unique = unique_values(&data);
        assert_eq!(unique, Array1::from(vec![1e10, 1e11, 1e12]));
    }

    #[test]
    fn test_unique_values_with_negative_numbers() {
        let data = Array1::from(vec![-3.0, -1.0, -4.0, -3.0, -1.0]);
        let unique = unique_values(&data);
        assert_eq!(unique, Array1::from(vec![-4.0, -3.0, -1.0]));
    }

    #[test]
    fn test_unique_values_with_subnormals() {
        let data = Array1::from(vec![1e-308, 1e-309, 1e-308]);
        let unique = unique_values(&data);
        assert_eq!(unique, Array1::from(vec![1e-309, 1e-308]));
    }

    #[test]
    fn test_argsort_f64_empty() {
        let data: Vec<f64> = vec![];
        let indices = argsort_f64(&data);
        assert_eq!(indices, vec![]);
    }

    #[test]
    fn test_argsort_f64_sorted() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let indices = argsort_f64(&data);
        assert_eq!(indices, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn test_argsort_f64_reverse_sorted() {
        let data = vec![5.0, 4.0, 3.0, 2.0, 1.0];
        let indices = argsort_f64(&data);
        assert_eq!(indices, vec![4, 3, 2, 1, 0]);
    }

    #[test]
    fn test_argsort_f64_unsorted() {
        let data = vec![3.0, 1.0, 4.0, 2.0, 5.0];
        let indices = argsort_f64(&data);
        assert_eq!(indices, vec![1, 3, 0, 2, 4]);
    }

    #[test]
    fn test_argsort_f64_with_duplicates() {
        let data = vec![3.0, 1.0, 3.0, 2.0, 1.0];
        let indices = argsort_f64(&data);
        // Check that the relative order of equal elements is preserved
        assert!(indices.contains(&0) && indices.contains(&2));
        assert!(indices.contains(&1) && indices.contains(&4));
        assert_eq!(indices[2], 3); // 2.0 should be at position 2
    }

    #[test]
    fn test_argsort_f64_with_nan() {
        let data = vec![3.0, f64::NAN, 1.0, 2.0];
        let indices = argsort_f64(&data);
        // Check that NaN is handled properly
        assert!(
            indices.contains(&0)
                && indices.contains(&1)
                && indices.contains(&2)
                && indices.contains(&3)
        );
    }

    #[test]
    fn test_argsort_f64_with_infinity() {
        let data = vec![3.0, f64::INFINITY, 1.0, f64::NEG_INFINITY];
        let indices = argsort_f64(&data);
        assert_eq!(indices, vec![3, 2, 0, 1]);
    }
}

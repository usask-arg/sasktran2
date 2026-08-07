//! Four-lane wavelength values used by the successive-orders batch kernels.
//!
//! The storage layout is array-of-structures: four consecutive `f64` values
//! belong to one spatial/angular element. That keeps every sparse transport
//! gather scalar while the wavelength arithmetic is contiguous and packed.

use std::ops::{Add, Div, Mul, Neg, Sub};

pub(crate) const LANES: usize = 4;

#[cfg(feature = "simd")]
#[derive(Clone, Copy, Debug)]
pub(crate) struct Batch4(std::simd::f64x4);

#[cfg(not(feature = "simd"))]
#[derive(Clone, Copy, Debug)]
pub(crate) struct Batch4([f64; LANES]);

impl Batch4 {
    #[inline(always)]
    pub(crate) fn splat(value: f64) -> Self {
        #[cfg(feature = "simd")]
        {
            Self(std::simd::f64x4::splat(value))
        }
        #[cfg(not(feature = "simd"))]
        {
            Self([value; LANES])
        }
    }

    #[inline(always)]
    pub(crate) fn from_array(values: [f64; LANES]) -> Self {
        #[cfg(feature = "simd")]
        {
            Self(std::simd::f64x4::from_array(values))
        }
        #[cfg(not(feature = "simd"))]
        {
            Self(values)
        }
    }

    #[inline(always)]
    pub(crate) fn to_array(self) -> [f64; LANES] {
        #[cfg(feature = "simd")]
        {
            self.0.to_array()
        }
        #[cfg(not(feature = "simd"))]
        {
            self.0
        }
    }

    #[inline(always)]
    pub(crate) fn load(values: &[f64], element: usize) -> Self {
        let start = element * LANES;
        debug_assert!(start + LANES <= values.len());
        #[cfg(feature = "simd")]
        {
            Self(std::simd::f64x4::from_slice(&values[start..start + LANES]))
        }
        #[cfg(not(feature = "simd"))]
        {
            Self(values[start..start + LANES].try_into().unwrap())
        }
    }

    #[inline(always)]
    pub(crate) fn gather(values: &[f64], element: usize, lane_stride: usize) -> Self {
        debug_assert!(element < lane_stride);
        debug_assert!(values.len() >= lane_stride * LANES);
        Self::from_array(std::array::from_fn(|lane| {
            values[lane * lane_stride + element]
        }))
    }

    #[inline(always)]
    pub(crate) fn store(self, values: &mut [f64], element: usize) {
        let start = element * LANES;
        debug_assert!(start + LANES <= values.len());
        #[cfg(feature = "simd")]
        self.0.copy_to_slice(&mut values[start..start + LANES]);
        #[cfg(not(feature = "simd"))]
        values[start..start + LANES].copy_from_slice(&self.0);
    }

    #[inline(always)]
    pub(crate) fn select(mask: [bool; LANES], when_true: Self, when_false: Self) -> Self {
        #[cfg(feature = "simd")]
        {
            Self(std::simd::Mask::<i64, LANES>::from_array(mask).select(when_true.0, when_false.0))
        }
        #[cfg(not(feature = "simd"))]
        {
            let when_true = when_true.to_array();
            let when_false = when_false.to_array();
            Self::from_array(std::array::from_fn(|lane| {
                if mask[lane] {
                    when_true[lane]
                } else {
                    when_false[lane]
                }
            }))
        }
    }

    #[inline(always)]
    pub(crate) fn select_abs_lt(
        value: Self,
        limit: f64,
        when_true: Self,
        when_false: Self,
    ) -> Self {
        #[cfg(feature = "simd")]
        {
            use std::simd::{cmp::SimdPartialOrd, num::SimdFloat};
            Self(
                value
                    .0
                    .abs()
                    .simd_lt(std::simd::f64x4::splat(limit))
                    .select(when_true.0, when_false.0),
            )
        }
        #[cfg(not(feature = "simd"))]
        {
            let value = value.to_array();
            let when_true = when_true.to_array();
            let when_false = when_false.to_array();
            Self::from_array(std::array::from_fn(|lane| {
                if value[lane].abs() < limit {
                    when_true[lane]
                } else {
                    when_false[lane]
                }
            }))
        }
    }

    #[inline(always)]
    pub(crate) fn exp(self) -> Self {
        #[cfg(feature = "simd")]
        {
            use std::simd::StdFloat;
            Self(self.0.exp())
        }
        #[cfg(not(feature = "simd"))]
        {
            Self(self.0.map(f64::exp))
        }
    }
}

macro_rules! batch_binary {
    ($trait:ident, $method:ident, $operator:tt) => {
        impl $trait for Batch4 {
            type Output = Self;

            #[inline(always)]
            fn $method(self, rhs: Self) -> Self {
                #[cfg(feature = "simd")]
                {
                    Self(self.0 $operator rhs.0)
                }
                #[cfg(not(feature = "simd"))]
                {
                    Self(std::array::from_fn(|lane| {
                        self.0[lane] $operator rhs.0[lane]
                    }))
                }
            }
        }
    };
}

batch_binary!(Add, add, +);
batch_binary!(Sub, sub, -);
batch_binary!(Mul, mul, *);
batch_binary!(Div, div, /);

impl Neg for Batch4 {
    type Output = Self;

    #[inline(always)]
    fn neg(self) -> Self {
        #[cfg(feature = "simd")]
        {
            Self(-self.0)
        }
        #[cfg(not(feature = "simd"))]
        {
            Self(self.0.map(|value| -value))
        }
    }
}

#[inline]
pub(crate) fn interleave_wavelength_major(source: &[f64], lane_stride: usize, target: &mut [f64]) {
    debug_assert_eq!(source.len(), lane_stride * LANES);
    debug_assert_eq!(target.len(), source.len());
    for element in 0..lane_stride {
        Batch4::gather(source, element, lane_stride).store(target, element);
    }
}

#[inline]
#[cfg_attr(test, allow(dead_code))]
pub(crate) fn interleave_lane(source: &[f64], lane: usize, target: &mut [f64]) {
    debug_assert!(lane < LANES);
    debug_assert_eq!(target.len(), source.len() * LANES);
    for (element, &value) in source.iter().enumerate() {
        target[element * LANES + lane] = value;
    }
}

#[inline]
#[cfg_attr(test, allow(dead_code))]
pub(crate) fn extract_lane(source: &[f64], lane: usize, target: &mut [f64]) {
    debug_assert!(lane < LANES);
    debug_assert_eq!(source.len(), target.len() * LANES);
    for (element, value) in target.iter_mut().enumerate() {
        *value = source[element * LANES + lane];
    }
}

#[inline]
pub(crate) fn infinity_norm(values: &[f64]) -> [f64; LANES] {
    debug_assert!(values.len().is_multiple_of(LANES));
    let mut norm = [0.0_f64; LANES];
    for element in 0..values.len() / LANES {
        let value = Batch4::load(values, element).to_array();
        for lane in 0..LANES {
            norm[lane] = norm[lane].max(value[lane].abs());
        }
    }
    norm
}

#[inline]
pub(crate) fn dot_product(left: &[f64], right: &[f64]) -> Batch4 {
    debug_assert_eq!(left.len(), right.len());
    debug_assert!(left.len().is_multiple_of(LANES));
    let mut result = Batch4::splat(0.0);
    for element in 0..left.len() / LANES {
        result = result + Batch4::load(left, element) * Batch4::load(right, element);
    }
    result
}

use integrate::gauss_quadrature::legendre_rule;


use crate::basis::*;

pub fn integrate_rectangle_delta(rect: &Rectangle, delta: &Delta) -> f64 {
    // If delta function is in support of rectangle, return 1.0, otherwise 0.0

    let delta_center = delta.lower_limit();

    if delta_center > rect.lower_limit() && delta_center < rect.upper_limit() {
        1.0
    } else {
        0.0
    }
}

pub fn integrate_triangle_delta(triangle: &Triangle, delta: &Delta) -> f64 {
    // If delta function is in support of triangle, return value of triangle at delta center, otherwise 0.0

    let delta_center = delta.lower_limit();

    if delta_center > triangle.lower_limit() && delta_center < triangle.upper_limit() {
        triangle.evaluate(delta_center)
    } else {
        0.0
    }
}

pub fn integrate_gaussian_delta(gaussian: &Gaussian, delta: &Delta) -> f64 {
    // Return value of gaussian at delta center

    let delta_center = delta.lower_limit();

    if delta_center < gaussian.lower_limit() || delta_center > gaussian.upper_limit() {
        return 0.0;
    }

    gaussian.evaluate(delta_center)
}

pub fn integrate_numeric(basis1: &BasisType, basis2: &BasisType) -> f64 {

    let min1 = basis1.lower_limit();
    let min2 = basis2.lower_limit();

    let max1 = basis2.upper_limit();
    let max2 = basis2.upper_limit();

    // overlap region is
    let integration_min = match min1 < min2 {
        true => min2,
        false => min1
    };

    let integration_max = match max1 < max2 {
        true => max1,
        false => max2
    };

    if integration_max <= integration_min {
        return 0.0; // No overlap
    }

    legendre_rule(|x: f64| basis1.evaluate(x) * basis2.evaluate(x), integration_min, integration_max, 21_usize)
}
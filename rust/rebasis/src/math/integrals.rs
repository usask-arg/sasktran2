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

    let max1 = basis1.upper_limit();
    let max2 = basis2.upper_limit();

    // overlap region is
    let integration_min = match min1 < min2 {
        true => min2,
        false => min1,
    };

    let integration_max = match max1 < max2 {
        true => max1,
        false => max2,
    };

    if integration_max <= integration_min {
        return 0.0; // No overlap
    }

    legendre_rule(
        |x: f64| basis1.evaluate(x) * basis2.evaluate(x),
        integration_min,
        integration_max,
        21_usize,
    )
}

pub fn integrate_triangle_triangle(t1: &Triangle, t2: &Triangle) -> f64 {
    let a1 = t1.lower_limit();
    let c1 = t1.center();
    let b1 = t1.upper_limit();

    let a2 = t2.lower_limit();
    let c2 = t2.center();
    let b2 = t2.upper_limit();

    // Overlap interval
    let l = a1.max(a2);
    let r = b1.min(b2);

    if r <= l {
        return 0.0;
    }

    // Breakpoints: overlap bounds + centers
    let mut knots = [l, r, c1, c2];
    knots.sort_by(|x, y| x.partial_cmp(y).unwrap());

    let mut result = 0.0;

    for w in knots.windows(2) {
        let x0 = w[0];
        let x1 = w[1];

        if x1 <= x0 {
            continue;
        }

        let xm = 0.5 * (x0 + x1);

        // Linear params for triangle 1: T(x) = p1*x + q1
        let (p1, q1) = triangle_linear_params(t1, xm);
        let (p2, q2) = triangle_linear_params(t2, xm);

        // Skip zero regions
        if (p1 == 0.0 && q1 == 0.0) || (p2 == 0.0 && q2 == 0.0) {
            continue;
        }

        result += integrate_linear_product(p1, q1, p2, q2, x0, x1);
    }

    result
}

#[inline]
fn triangle_linear_params(t: &Triangle, x: f64) -> (f64, f64) {
    let a = t.lower_limit();
    let c = t.center();
    let b = t.upper_limit();

    if x <= a || x >= b {
        return (0.0, 0.0);
    }

    if x <= c {
        // rising side
        let p = 1.0 / (c - a);
        let q = -a / (c - a);
        (p, q)
    } else {
        // falling side
        let p = -1.0 / (b - c);
        let q = b / (b - c);
        (p, q)
    }
}

#[inline]
fn integrate_linear_product(p1: f64, q1: f64, p2: f64, q2: f64, x0: f64, x1: f64) -> f64 {
    // ∫ (p1 x + q1)(p2 x + q2) dx over [x0,x1]

    let a = p1 * p2;
    let b = p1 * q2 + p2 * q1;
    let c = q1 * q2;

    let x0_2 = x0 * x0;
    let x1_2 = x1 * x1;
    let x0_3 = x0_2 * x0;
    let x1_3 = x1_2 * x1;

    a * (x1_3 - x0_3) / 3.0 + b * (x1_2 - x0_2) / 2.0 + c * (x1 - x0)
}

use crate::math::integrals::*;

pub trait Basis {
    fn lower_limit(&self) -> f64;
    fn upper_limit(&self) -> f64;
    fn center(&self) -> f64;
    fn evaluate(&self, x: f64) -> f64;
}

#[derive(Debug, Clone)]
pub enum BasisType {
    Rectangle(Rectangle),
    Delta(Delta),
    Gaussian(Gaussian),
    Triangle(Triangle),
}

impl BasisType {
    pub fn overlap_integral(&self, other: &BasisType) -> f64 {
        match (self, other) {
            (BasisType::Rectangle(r), BasisType::Delta(d)) => integrate_rectangle_delta(r, d),
            (BasisType::Delta(d), BasisType::Rectangle(r)) => integrate_rectangle_delta(r, d),
            (BasisType::Triangle(t), BasisType::Delta(d)) => integrate_triangle_delta(t, d),
            (BasisType::Delta(d), BasisType::Triangle(t)) => integrate_triangle_delta(t, d),
            (BasisType::Gaussian(g), BasisType::Delta(d)) => integrate_gaussian_delta(g, d),
            (BasisType::Delta(d), BasisType::Gaussian(g)) => integrate_gaussian_delta(g, d),
            (BasisType::Triangle(t1), BasisType::Triangle(t2)) => {
                integrate_triangle_triangle(t1, t2)
            }
            (b1, b2) => integrate_numeric(b1, b2),
        }
    }

    pub fn lower_limit(&self) -> f64 {
        match self {
            BasisType::Rectangle(r) => r.lower_limit(),
            BasisType::Delta(d) => d.lower_limit(),
            BasisType::Gaussian(g) => g.lower_limit(),
            BasisType::Triangle(t) => t.lower_limit(),
        }
    }

    pub fn upper_limit(&self) -> f64 {
        match self {
            BasisType::Rectangle(r) => r.upper_limit(),
            BasisType::Delta(d) => d.upper_limit(),
            BasisType::Gaussian(g) => g.upper_limit(),
            BasisType::Triangle(t) => t.upper_limit(),
        }
    }

    pub fn evaluate(&self, x: f64) -> f64 {
        match self {
            BasisType::Rectangle(r) => r.evaluate(x),
            BasisType::Delta(d) => d.evaluate(x),
            BasisType::Gaussian(g) => g.evaluate(x),
            BasisType::Triangle(t) => t.evaluate(x),
        }
    }

    pub fn center(&self) -> f64 {
        match self {
            BasisType::Rectangle(r) => r.center(),
            BasisType::Delta(d) => d.center(),
            BasisType::Gaussian(g) => g.center(),
            BasisType::Triangle(t) => t.center(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Rectangle {
    left: f64,
    right: f64,
    norm_factor: f64,
}

impl Rectangle {
    pub fn new(left: f64, right: f64) -> Self {
        let norm_factor = 1.0 / (right - left);

        Rectangle {
            left,
            right,
            norm_factor,
        }
    }
}

impl Basis for Rectangle {
    fn lower_limit(&self) -> f64 {
        self.left
    }

    fn upper_limit(&self) -> f64 {
        self.right
    }

    fn evaluate(&self, x: f64) -> f64 {
        if x > self.left && x < self.right {
            1.0 * self.norm_factor
        } else {
            0.0
        }
    }

    fn center(&self) -> f64 {
        0.5 * (self.left + self.right)
    }
}

#[derive(Debug, Clone)]
pub struct Delta {
    center: f64,
}

impl Delta {
    pub fn new(center: f64) -> Self {
        Delta { center }
    }
}

impl Basis for Delta {
    fn lower_limit(&self) -> f64 {
        self.center
    }

    fn upper_limit(&self) -> f64 {
        self.center
    }

    fn evaluate(&self, _x: f64) -> f64 {
        0.0
    }

    fn center(&self) -> f64 {
        self.center
    }
}

#[derive(Debug, Clone)]
pub struct Gaussian {
    center: f64,
    stdev: f64,
    max_stdev: i32,
    norm_factor: f64,
}

impl Gaussian {
    pub fn new(center: f64, stdev: f64, max_stdev: i32) -> Self {
        let norm_factor = 1.0 / (stdev * (2.0 * std::f64::consts::PI).sqrt());

        Gaussian {
            center,
            stdev,
            max_stdev,
            norm_factor,
        }
    }
}

impl Basis for Gaussian {
    fn lower_limit(&self) -> f64 {
        self.center - self.max_stdev as f64 * self.stdev
    }

    fn upper_limit(&self) -> f64 {
        self.center + self.max_stdev as f64 * self.stdev
    }

    fn evaluate(&self, x: f64) -> f64 {
        let exponent = -((x - self.center).powi(2)) / (2.0 * self.stdev.powi(2));
        self.norm_factor * exponent.exp()
    }

    fn center(&self) -> f64 {
        self.center
    }
}

#[derive(Debug, Clone)]
pub struct Triangle {
    left: f64,
    right: f64,
    center: f64,
    norm_factor: f64,
}

impl Triangle {
    pub fn new(left: f64, right: f64, center: f64) -> Self {
        // integral will be
        // (center - left) * height / 2
        // + (right - center) * height / 2
        // = (right - left) * height / 2
        // then if we want this to equal 1, we need to set height = 2 / (right - left)

        Triangle {
            left,
            right,
            center,
            norm_factor: 2.0 / (right - left),
        }
    }
}

impl Basis for Triangle {
    fn lower_limit(&self) -> f64 {
        self.left
    }

    fn upper_limit(&self) -> f64 {
        self.right
    }

    fn evaluate(&self, x: f64) -> f64 {
        if x < self.left || x > self.right {
            return 0.0;
        }
        if x < self.center {
            (x - self.left) / (self.center - self.left) * self.norm_factor
        } else {
            (self.right - x) / (self.right - self.center) * self.norm_factor
        }
    }

    fn center(&self) -> f64 {
        self.center
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basis_type() {
        let _ = BasisType::Rectangle(Rectangle::new(0.0, 1.0));
    }

    #[test]
    fn test_dispatch() {
        let b = BasisType::Rectangle(Rectangle::new(0.0, 1.0));

        let rect2 = BasisType::Rectangle(Rectangle::new(0.0, 1.0));
        let delta = BasisType::Delta(Delta { center: 0.5 });

        let result = b.overlap_integral(&rect2);
        assert!((result - 1.0).abs() < 1e-8);

        let result = b.overlap_integral(&delta);
        assert!((result - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_rectangle() {
        let _ = Rectangle::new(0.0, 1.0);
    }
}

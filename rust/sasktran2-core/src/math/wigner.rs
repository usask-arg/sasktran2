use std::f64;

/// A struct for computing Wigner D functions (and Legendre polynomials).
#[derive(Debug)]
pub struct WignerDCalculator {
    m_m: i32,
    m_n: i32,
    m_lmin: i32,
    m_recurrence_start_factor: f64,
    m_zeta: i32,
}

impl WignerDCalculator {
    /// Creates a new WignerDCalculator for d^l_{m,n}.
    pub fn new(m: i32, n: i32) -> Self {
        // Calculate zeta from F.2.4
        let zeta = if n >= m {
            1
        } else {
            // If (m - n) is even => zeta = 1, else -1
            if (m - n) % 2 == 0 { 1 } else { -1 }
        };

        let lmin = std::cmp::max(m.abs(), n.abs());

        let mut calc = WignerDCalculator {
            m_m: m,
            m_n: n,
            m_lmin: lmin,
            m_recurrence_start_factor: 0.0,
            m_zeta: zeta,
        };

        calc.m_recurrence_start_factor = calc.recurrence_start_factor();
        calc
    }

    /// Computes d^l_{m,n}(theta).
    pub fn d(&self, theta: f64, l: i32) -> f64 {
        if l < self.m_lmin {
            0.0
        } else {
            let x = theta.cos();
            // Value at current l
            let mut val_l = self.recurrence_start(theta);
            // Value at current (l - 1)
            let mut val_lm1 = 0.0;

            // If n == 0, use the associated Legendre function recurrence.
            if self.m_n == 0 {
                for lidx in (self.m_lmin + 1)..=l {
                    let lidx_f = lidx as f64;
                    let mm = (lidx * lidx - self.m_m * self.m_m) as f64;
                    // multiplier = 1 / (sqrt(lidx*lidx - m_m*m_m) * lidx)
                    let multiplier = 1.0 / (mm.sqrt() * lidx_f);

                    // (2*lidx - 1)*(lidx * x)
                    let curfactor = (2 * lidx - 1) as f64 * (lidx_f * x);

                    // lidx * sqrt((lidx - 1)^2 - m_m^2)
                    let priorfactor =
                        lidx_f * (((lidx - 1) * (lidx - 1) - self.m_m * self.m_m) as f64).sqrt();

                    let temp = val_l;
                    val_l = multiplier * (curfactor * val_l - priorfactor * val_lm1);
                    val_lm1 = temp;
                }
            } else {
                // Use F.4.1 from Mischenko
                for lidx in (self.m_lmin + 1)..=l {
                    let lidx_f = lidx as f64;
                    let mm = (lidx * lidx - self.m_m * self.m_m) as f64;
                    let nn = (lidx * lidx - self.m_n * self.m_n) as f64;

                    // multiplier = 1 / ((lidx - 1) * sqrt(lidx^2 - m_m^2) * sqrt(lidx^2 - m_n^2))
                    let multiplier = 1.0 / (((lidx - 1) as f64) * mm.sqrt() * nn.sqrt());

                    // (2*lidx - 1)* (lidx * (lidx-1) * x - m_n*m_m)
                    let curfactor = (2 * lidx - 1) as f64
                        * (lidx_f * (lidx_f - 1.0) * x - (self.m_n * self.m_m) as f64);

                    // lidx * sqrt((lidx-1)^2 - m_m^2) * sqrt((lidx-1)^2 - m_n^2)
                    let priorfactor = lidx_f
                        * (((lidx - 1) * (lidx - 1) - self.m_m * self.m_m) as f64).sqrt()
                        * (((lidx - 1) * (lidx - 1) - self.m_n * self.m_n) as f64).sqrt();

                    let temp = val_l;
                    val_l = multiplier * (curfactor * val_l - priorfactor * val_lm1);
                    val_lm1 = temp;
                }
            }

            val_l
        }
    }

    /// Computes d^l_{m,n}(theta).
    pub fn vector_d(&self, theta: f64, output: &mut [f64]) {
        for (l, elem) in output.iter_mut().enumerate() {
            if (l as i32) < self.m_lmin {
                *elem = 0.0;
            } else {
                let x = theta.cos();
                // Value at current l
                let mut val_l = self.recurrence_start(theta);
                // Value at current (l - 1)
                let mut val_lm1 = 0.0;

                // If n == 0, use the associated Legendre function recurrence.
                if self.m_n == 0 {
                    for lidx in (self.m_lmin + 1)..=(l as i32) {
                        let lidx_f = lidx as f64;
                        let mm = (lidx * lidx - self.m_m * self.m_m) as f64;
                        // multiplier = 1 / (sqrt(lidx*lidx - m_m*m_m) * lidx)
                        let multiplier = 1.0 / (mm.sqrt() * lidx_f);

                        // (2*lidx - 1)*(lidx * x)
                        let curfactor = (2 * lidx - 1) as f64 * (lidx_f * x);

                        // lidx * sqrt((lidx - 1)^2 - m_m^2)
                        let priorfactor = lidx_f
                            * (((lidx - 1) * (lidx - 1) - self.m_m * self.m_m) as f64).sqrt();

                        let temp = val_l;
                        val_l = multiplier * (curfactor * val_l - priorfactor * val_lm1);
                        val_lm1 = temp;
                    }
                } else {
                    // Use F.4.1 from Mischenko
                    for lidx in (self.m_lmin + 1)..=(l as i32) {
                        let lidx_f = lidx as f64;
                        let mm = (lidx * lidx - self.m_m * self.m_m) as f64;
                        let nn = (lidx * lidx - self.m_n * self.m_n) as f64;

                        // multiplier = 1 / ((lidx - 1) * sqrt(lidx^2 - m_m^2) * sqrt(lidx^2 - m_n^2))
                        let multiplier = 1.0 / (((lidx - 1) as f64) * mm.sqrt() * nn.sqrt());

                        // (2*lidx - 1)* (lidx * (lidx-1) * x - m_n*m_m)
                        let curfactor = (2 * lidx - 1) as f64
                            * (lidx_f * (lidx_f - 1.0) * x - (self.m_n * self.m_m) as f64);

                        // lidx * sqrt((lidx-1)^2 - m_m^2) * sqrt((lidx-1)^2 - m_n^2)
                        let priorfactor = lidx_f
                            * (((lidx - 1) * (lidx - 1) - self.m_m * self.m_m) as f64).sqrt()
                            * (((lidx - 1) * (lidx - 1) - self.m_n * self.m_n) as f64).sqrt();

                        let temp = val_l;
                        val_l = multiplier * (curfactor * val_l - priorfactor * val_lm1);
                        val_lm1 = temp;
                    }
                }

                *elem = val_l;
            }
        }
    }

    /// Internal method to compute the start scaling factor of the recurrence.
    fn recurrence_start_factor(&self) -> f64 {
        // Start of the recurrence, i.e. d^lmin_{m,n}(theta)
        // Equation F.4.3: the factor that doesn't depend on x.

        // Numerator = (2*m_lmin)!
        // Denominator = |m_m - m_n|! * |m_m + m_n|!
        let num = 2 * self.m_lmin;
        let den1 = (self.m_m - self.m_n).abs();
        let den2 = (self.m_m + self.m_n).abs();

        // We'll replicate the factorial logic from the C++ code:
        //   factorial(num)! / [ factorial(den1)! * factorial(den2)! ]
        // done by iteratively multiplying and dividing.
        let mut factorial = 1.0;
        let mut d1 = den1;
        let mut d2 = den2;

        // Loop down from num to 2
        for i in (2..=num).rev() {
            factorial *= i as f64;
            if d1 >= i {
                factorial /= i as f64;
                d1 -= 1;
            }
            if d2 >= i {
                factorial /= i as f64;
                d2 -= 1;
            }
        }

        // otherfactor = 2^(-m_lmin)
        let otherfactor = 2.0_f64.powf(-(self.m_lmin as f64));

        // Return zeta * 2^(-lmin) * sqrt(factorial)
        (self.m_zeta as f64) * otherfactor * factorial.sqrt()
    }

    /// Internal method to get the initial recurrence value (theta-dependent).
    fn recurrence_start(&self, theta: f64) -> f64 {
        let x = theta.cos();
        // (1 - x)^{|m_m - m_n|/2} * (1 + x)^{|m_m + m_n|/2}
        let xfactor = (1.0 - x).powf(((self.m_m - self.m_n).abs() as f64) / 2.0)
            * (1.0 + x).powf(((self.m_m + self.m_n).abs() as f64) / 2.0);

        self.m_recurrence_start_factor * xfactor
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wigner_d_basic() {
        let calc = WignerDCalculator::new(1, 1);
        let val = calc.d(std::f64::consts::PI / 3.0, 2);
        // For a real test, you would compare with some known expected result
        println!("d^2_{{1,1}}(π/3) = {val}");
    }
}

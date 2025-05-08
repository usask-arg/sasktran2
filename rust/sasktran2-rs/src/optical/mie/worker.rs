#![allow(non_snake_case)]

use ndarray::{Array1, Array2, Zip};
use num::abs;
use num::complex::Complex64;

pub struct MieOutput {
    pub Qext: Array1<f64>,
    pub Qsca: Array1<f64>,
    pub S1: Array2<Complex64>,
    pub S2: Array2<Complex64>,
}

fn max_order(size_param: f64) -> usize {
    (size_param + 4.05 * size_param.powf(0.3333) + 2.0).round() as usize
}

fn Dn_lentz(z: Complex64, N: usize) -> Complex64 {
    let mut z_inv = 2.0 / z;
    let alpha_1 = (N as f64 + 0.5) * z_inv;

    let mut a_j = -(N as f64 + 1.5) * z_inv;
    let mut alpha_j1 = a_j + 1.0 / alpha_1;
    let mut alpha_j2 = a_j;

    let mut cur_ratio = alpha_j1 / alpha_j2;
    let mut overall_ratio = alpha_1 * cur_ratio;

    while abs(cur_ratio.norm() - 1.0) > 1e-12 {
        a_j = z_inv - a_j;
        alpha_j1 = a_j + 1.0 / alpha_j1;
        alpha_j2 = a_j + 1.0 / alpha_j2;

        z_inv *= -1.0;
        cur_ratio = alpha_j1 / alpha_j2;
        overall_ratio *= cur_ratio;
    }

    -(N as f64) / z + overall_ratio
}

fn tau_pi_matrices(cos_angles: &Array1<f64>, max_order: usize) -> (Array2<f64>, Array2<f64>) {
    let mut tau = Array2::<f64>::zeros((max_order, cos_angles.len()));
    let mut pi = Array2::<f64>::zeros((max_order, cos_angles.len()));

    let N = max_order;

    pi.slice_mut(ndarray::s![0, ..]).fill(1.0);
    tau.slice_mut(ndarray::s![0, ..]).assign(cos_angles);

    if N > 1 {
        for (r, &a) in pi
            .slice_mut(ndarray::s![1, ..])
            .iter_mut()
            .zip(cos_angles.iter())
        {
            *r = a * 3.0;
        }

        for ((r, &p), &a) in tau
            .slice_mut(ndarray::s![1, ..])
            .iter_mut()
            .zip(pi.slice(ndarray::s![1, ..]).iter())
            .zip(cos_angles.iter())
        {
            *r = 2.0 * a * p - 3.0;
        }
    }

    for n in 3..N + 1 {
        let nf = n as f64;

        // Copy out rows that we need for reads
        let pi_n2 = pi.slice(ndarray::s![n - 2, ..]).to_owned();
        let pi_n3 = pi.slice(ndarray::s![n - 3, ..]).to_owned();

        // Compute new values into a temporary buffer
        let mut new_pi_n1 = pi_n2.clone(); // same shape
        Zip::from(&mut new_pi_n1)
            .and(&pi_n2)
            .and(&pi_n3)
            .and(cos_angles)
            .for_each(|r, &p2, &p3, &a| {
                *r = (2.0 * nf - 1.0) * a * p2 - nf * p3 / (nf - 1.0);
            });

        // Now do a single mutable borrow to write back
        pi.slice_mut(ndarray::s![n - 1, ..]).assign(&new_pi_n1);

        // Repeat the pattern for tau
        let pi_n1_read = pi.slice(ndarray::s![n - 1, ..]).to_owned();
        let pi_n2_read = pi.slice(ndarray::s![n - 2, ..]).to_owned();
        let mut new_tau_n1 = pi_n1_read.clone();
        Zip::from(&mut new_tau_n1)
            .and(&pi_n1_read)
            .and(&pi_n2_read)
            .and(cos_angles)
            .for_each(|r, &p1, &p2, &a| {
                *r = nf * a * p1 - (nf + 1.0) * p2;
            });
        tau.slice_mut(ndarray::s![n - 1, ..]).assign(&new_tau_n1);
    }

    (tau, pi)
}

pub struct Mie {
    An: Vec<Complex64>,
    Bn: Vec<Complex64>,
    Dn: Vec<Complex64>,
    cos_angles: Array1<f64>,
    tau: Array2<f64>,
    pi: Array2<f64>,
    S1: Array1<Complex64>,
    S2: Array1<Complex64>,
    recalculate_tau: bool,
}

impl Default for Mie {
    fn default() -> Self {
        Mie::new()
    }
}

impl Mie {
    pub fn new() -> Self {
        Mie {
            An: Vec::new(),
            Bn: Vec::new(),
            Dn: Vec::new(),
            cos_angles: Array1::from_vec(vec![]),
            tau: Array2::zeros((0, 0)),
            pi: Array2::zeros((0, 0)),
            S1: Array1::from_vec(vec![]),
            S2: Array1::from_vec(vec![]),
            recalculate_tau: false,
        }
    }

    pub fn with_cos_angles(mut self, cos_angles: Vec<f64>) -> Self {
        self.cos_angles = Array1::from_vec(cos_angles);
        self.S1 = Array1::from_vec(vec![Complex64::new(0.0, 0.0); self.cos_angles.len()]);
        self.S2 = Array1::from_vec(vec![Complex64::new(0.0, 0.0); self.cos_angles.len()]);
        self.recalculate_tau = true;
        self
    }

    fn allocate(&mut self, N: usize) {
        self.An.resize(N, Complex64::new(0.0, 0.0));
        self.Bn.resize(N, Complex64::new(0.0, 0.0));
        self.Dn.resize(N, Complex64::new(0.0, 0.0));
    }

    pub fn calculate(&mut self, size_param: f64, refractive_index: Complex64) -> (f64, f64) {
        let N = max_order(size_param);

        if N > self.current_N() {
            self.allocate(N);
            self.recalculate_tau = true;
        }

        if self.recalculate_tau {
            (self.tau, self.pi) = tau_pi_matrices(&self.cos_angles, N);
            self.recalculate_tau = false;
        }

        if refractive_index.norm() * size_param < 0.1 {
            self.small_Q_S(refractive_index, size_param)
        } else {
            self.Dn(refractive_index, size_param);
            self.An_Bn(refractive_index, size_param);
            self.regular_Q_S(size_param)
        }
    }

    fn current_N(&self) -> usize {
        self.An.len()
    }

    fn Dn_downwards(&mut self, z: Complex64) {
        let N = self.current_N();

        self.Dn[N - 1] = Dn_lentz(z, N);
        for n in (2..N + 1).rev() {
            self.Dn[n - 2] = (n as f64) / z - 1.0 / ((n as f64) / z + self.Dn[n - 1]);
        }
    }

    fn Dn_upwards(&mut self, z: Complex64) {
        let N = self.current_N();

        let j: Complex64 = Complex64::new(0.0, 1.0);
        let exp_2 = (-2.0 * j * z).exp();

        self.Dn[0] = -1.0 / z + (1.0 - exp_2) / ((1.0 - exp_2) / z - j * (1.0 + exp_2));

        for n in 2..N + 1 {
            self.Dn[n - 1] = -(n as f64) / z + 1.0 / ((n as f64) / z - self.Dn[n - 2]);
        }
    }

    fn Dn(&mut self, refractive_index: Complex64, size_param: f64) {
        let z = refractive_index * size_param;

        if refractive_index.re < 1.0
            || refractive_index.re > 10.0
            || abs(refractive_index.im) > 10.0
        {
            self.Dn_downwards(z);
        } else {
            let temp = 3.9 - 10.8 * refractive_index.re + 13.78 * refractive_index.re.powf(2.0);

            if abs(refractive_index.im) * size_param >= temp {
                self.Dn_downwards(z);
            } else {
                self.Dn_upwards(z);
            }
        }
    }

    fn An_Bn(&mut self, refractive_index: Complex64, size_param: f64) {
        let N = self.current_N();

        // Convenience
        let An = &mut self.An;
        let Bn = &mut self.Bn;
        let Dn = &self.Dn;

        let j = Complex64::new(0.0, 1.0);

        let mut psi_n_1: Complex64 = (size_param).sin().into();
        let mut psi_n: Complex64 = (size_param.sin() / size_param - size_param.cos()).into();

        let mut xi_n_1 = size_param.sin() + j * size_param.cos();
        let mut xi_n = size_param.sin() / size_param - size_param.cos()
            + j * (size_param.cos() / size_param + size_param.sin());

        let mut temp: Complex64 = xi_n_1;

        An[0] = ((Dn[0] / refractive_index + 1.0 / size_param) * psi_n - psi_n_1)
            / ((Dn[0] / refractive_index + 1.0 / size_param) * xi_n - xi_n_1);

        Bn[0] = ((Dn[0] * refractive_index + 1.0 / size_param) * psi_n - psi_n_1)
            / ((Dn[0] * refractive_index + 1.0 / size_param) * xi_n - xi_n_1);

        xi_n_1 = xi_n;
        xi_n = (2.0 + 1.0) / size_param * xi_n_1 - temp;

        psi_n_1 = psi_n;
        psi_n = xi_n.re.into();

        for n in 2..N + 1 {
            let nf = n as f64;
            An[n - 1] = ((Dn[n - 1] / refractive_index + nf / size_param) * psi_n - psi_n_1)
                / ((Dn[n - 1] / refractive_index + nf / size_param) * xi_n - xi_n_1);

            Bn[n - 1] = ((Dn[n - 1] * refractive_index + nf / size_param) * psi_n - psi_n_1)
                / ((Dn[n - 1] * refractive_index + nf / size_param) * xi_n - xi_n_1);

            temp = xi_n_1;
            xi_n_1 = xi_n;
            xi_n = (2.0 * nf + 1.0) / size_param * xi_n_1 - temp;

            psi_n_1 = psi_n;
            psi_n = xi_n.re.into();
        }
    }

    fn regular_Q_S(&mut self, size_param: f64) -> (f64, f64) {
        let Qext: f64 = self
            .An
            .iter()
            .zip(self.Bn.iter())
            .enumerate()
            .map(|(n, (An, Bn))| {
                let n = n as f64;

                (2.0 * (n + 1.0) + 1.0) * (An.re + Bn.re)
            })
            .sum();

        let Qsca: f64 = self
            .An
            .iter()
            .zip(self.Bn.iter())
            .enumerate()
            .map(|(n, (An, Bn))| {
                let n = n as f64;

                (2.0 * (n + 1.0) + 1.0) * (An.norm_sqr() + Bn.norm_sqr())
            })
            .sum();

        self.S1.fill(Complex64::new(0.0, 0.0));
        self.S2.fill(Complex64::new(0.0, 0.0));

        for i in 0..self.current_N() {
            let n = i as f64;
            let n_factor = (2.0 * (n + 1.0) + 1.0) / (n + 1.0) / (n + 2.0);

            Zip::from(&mut self.S1)
                .and(&mut self.S2)
                .and(&self.tau.row(i))
                .and(&self.pi.row(i))
                .for_each(|s1, s2, &tau, &pi| {
                    *s1 += n_factor * (self.An[i] * pi + self.Bn[i] * tau);
                    *s2 += n_factor * (self.An[i] * tau + self.Bn[i] * pi);
                });
        }

        (
            Qext * 2.0 / (size_param * size_param),
            Qsca * 2.0 / (size_param * size_param),
        )
    }

    fn small_Q_S(&mut self, refractive_index: Complex64, size_param: f64) -> (f64, f64) {
        let m_2 = refractive_index * refractive_index;
        let x_2 = size_param * size_param;
        let j = Complex64::new(0.0, 1.0);

        let N_1 = 1.0 - 0.1 * x_2 + (4.0 * m_2 + 5.0) * x_2 * x_2 / 1400.0;

        let mut D_1 = m_2 + 2.0 + (1.0 - 0.7 * m_2) * x_2
            - (8.0 * m_2 * m_2 - 385.0 * m_2 + 350.0) * x_2 * x_2 / 1400.0;

        D_1 += 2.0 * j * (m_2 - 1.0) * size_param.powf(3.0) * (1.0 - 0.1 * x_2) / 3.0;

        let a_hat1 = 2.0 * j * (m_2 - 1.0) / 3.0 * N_1 / D_1;

        let b_hat1 = (j * x_2 * (m_2 - 1.0) / 45.0 * (1.0 + (2.0 * m_2 - 5.0) / 70.0 * x_2))
            / (1.0 - x_2 * (2.0 * m_2 - 5.0) / 30.0);

        let a_hat2 = (j * x_2 * (m_2 - 1.0) / 15.0 * (1.0 - x_2 / 14.0))
            / (2.0 * m_2 + 3.0 - x_2 * (2.0 * m_2 - 7.0) / 14.0);

        let T = a_hat1.norm_sqr() + b_hat1.norm_sqr() + 5.0 / 3.0 * a_hat2.norm_sqr();

        let Qsca = 6.0 * x_2 * x_2 * T;

        // Angle calculation
        Zip::from(&mut self.S1)
            .and(&mut self.S2)
            .and(&self.cos_angles)
            .for_each(|s1, s2, &cos_theta| {
                *s1 = 3.0 / 2.0
                    * size_param.powf(3.0)
                    * (a_hat1 + cos_theta * (b_hat1 + 5.0 / 3.0 * a_hat2));
                *s2 = 3.0 / 2.0
                    * size_param.powf(3.0)
                    * (b_hat1
                        + a_hat1 * cos_theta
                        + 5.0 / 3.0 * a_hat2 * (2.0 * cos_theta * cos_theta - 1.0));
            });

        if refractive_index.im == 0.0 {
            (Qsca, Qsca)
        } else {
            (
                6.0 * size_param * ((a_hat1 + b_hat1).re + 5.0 / 3.0 * a_hat2.re),
                Qsca,
            )
        }
    }

    // Accessors
    pub fn S1(&self) -> &Array1<Complex64> {
        &self.S1
    }
    pub fn S2(&self) -> &Array1<Complex64> {
        &self.S2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::complex::Complex64;

    fn default_cos_angles() -> Vec<f64> {
        let mut cos_angles = Vec::new();
        for i in 0..7 {
            let angle = (i as f64) * std::f64::consts::PI / 6.0;
            cos_angles.push(angle.cos());
        }
        cos_angles
    }

    #[test]
    fn test_basic_execution_mie() {
        let mut mie = Mie::new();
        let size_param = 1.0;
        let refractive_index = Complex64::new(1.5, 0.0);
        mie.calculate(size_param, refractive_index);
    }

    #[test]
    fn test_basic_non_absorbing() {
        let lambda = 0.6328;
        let radius = 0.525;

        let size_param = 2.0 * std::f64::consts::PI * radius / lambda;
        let refractive_index = Complex64::new(1.55, 0.0);

        let mut mie = Mie::new();
        let (Qext, Qsca) = mie.calculate(size_param, refractive_index);

        assert!((Qext - 3.10543).abs() < 1e-5);
        assert!((Qsca - 3.10543).abs() < 1e-5);
    }

    #[test]
    fn test_miev0_case_5() {
        let size_param = 0.099;
        let refractive_index = Complex64::new(0.75, 0.0);

        let mut mie = Mie::new().with_cos_angles(default_cos_angles());
        let (_Qext, Qsca) = mie.calculate(size_param, refractive_index);

        assert!((Qsca - 0.000007).abs() < 1e-5);

        assert!((mie.S1()[0].re - 1.817558e-8).abs() < 1e-8);
        assert!((mie.S1()[0].im + 1.654225e-4).abs() < 1e-8);
        assert!((mie.S2()[0].re - 1.817558e-8).abs() < 1e-8);
        assert!((mie.S2()[0].im + 1.654225e-4).abs() < 1e-8);

        assert!((mie.S1()[1].re - 1.817558e-8).abs() < 1e-8);
        assert!((mie.S1()[1].im + 1.653815e-4).abs() < 1e-6);
        assert!((mie.S2()[1].re - 1.574051e-8).abs() < 1e-8);
        assert!((mie.S2()[1].im + 1.432172e-4).abs() < 1e-8);

        assert!((mie.S1()[2].re - 1.817558e-8).abs() < 1e-8);
        assert!((mie.S1()[2].im + 1.652694e-4).abs() < 1e-8);
        assert!((mie.S2()[2].re - 9.087788e-9).abs() < 1e-8);
        assert!((mie.S2()[2].im + 8.261265e-5).abs() < 1e-8);

        assert!((mie.S1()[3].re - 1.817558e-8).abs() < 1e-8);
        assert!((mie.S1()[3].im + 1.651163e-4).abs() < 1e-8);
        assert!((mie.S2()[3].re - 9.797186e-23).abs() < 1e-8);
        assert!((mie.S2()[3].im - 2.938374e-8).abs() < 1e-8);

        assert!((mie.S1()[4].re - 1.817558e-8).abs() < 1e-8);
        assert!((mie.S1()[4].im + 1.649631e-4).abs() < 1e-8);
        assert!((mie.S2()[4].re + 9.087788e-9).abs() < 1e-8);
        assert!((mie.S2()[4].im - 8.250360e-5).abs() < 1e-8);

        assert!((mie.S1()[5].re - 1.817558e-8).abs() < 1e-8);
        assert!((mie.S1()[5].im + 1.648510e-4).abs() < 1e-8);
        assert!((mie.S2()[5].re + 1.574051e-8).abs() < 1e-8);
        assert!((mie.S2()[5].im - 1.427725e-4).abs() < 1e-8);

        assert!((mie.S1()[6].re - 1.817558e-8).abs() < 1e-8);
        assert!((mie.S1()[6].im + 1.648100e-4).abs() < 1e-8);
        assert!((mie.S2()[6].re + 1.817558e-8).abs() < 1e-8);
        assert!((mie.S2()[6].im - 1.648100e-4).abs() < 1e-8);
    }

    #[test]
    fn test_miev0_case_6() {
        let size_param = 0.101;
        let refractive_index = Complex64::new(0.75, 0.0);

        let mut mie = Mie::new();
        let (Qext, Qsca) = mie.calculate(size_param, refractive_index);
        assert!((Qsca - 0.000008).abs() < 1e-6);
        assert!((Qext - 0.000008).abs() < 1e-6);
    }

    #[test]
    fn test_miev0_case_7() {
        let size_param = 10.0;
        let refractive_index = Complex64::new(0.75, 0.0);

        let mut mie = Mie::new();
        let (_Qext, Qsca) = mie.calculate(size_param, refractive_index);

        assert!((Qsca - 2.232265).abs() < 1e-6);
    }

    #[test]
    fn test_miev0_case_8() {
        let size_param = 1000.0;
        let refractive_index = Complex64::new(0.75, 0.0);

        let mut mie = Mie::new();
        let (_Qext, Qsca) = mie.calculate(size_param, refractive_index);

        assert!((Qsca - 1.997908).abs() < 1e-6);
    }

    #[test]
    fn test_miev0_case_1() {
        let size_param = 10.0;
        let refractive_index = Complex64::new(1.5, 0.0);

        let mut mie = Mie::new();
        let (_Qext, Qsca) = mie.calculate(size_param, refractive_index);
        assert!((Qsca - 2.8820).abs() < 1e-4);
    }

    #[test]
    fn test_miev0_case_2() {
        let size_param = 100.0;
        let refractive_index = Complex64::new(1.5, 0.0);

        let mut mie = Mie::new();
        let (_Qext, Qsca) = mie.calculate(size_param, refractive_index);
        assert!((Qsca - 2.0944).abs() < 1e-4);
    }

    #[test]
    fn test_miev0_case_3() {
        let size_param = 1000.0;
        let refractive_index = Complex64::new(1.5, 0.0);

        let mut mie = Mie::new();
        let (_Qext, Qsca) = mie.calculate(size_param, refractive_index);
        assert!((Qsca - 2.0139).abs() < 1e-4);
    }

    #[test]
    fn test_miev0_case_4() {
        let size_param = 5000.0;
        let refractive_index = Complex64::new(1.5, 0.0);

        let mut mie = Mie::new();
        let (_Qext, Qsca) = mie.calculate(size_param, refractive_index);
        assert!((Qsca - 2.0086).abs() < 1e-4);
    }

    #[test]
    fn test_non_dialetric() {
        let lambda = 0.6328;
        let radius = 0.525;

        let size_param = 2.0 * std::f64::consts::PI * radius / lambda;
        let refractive_index = Complex64::new(1.55, -0.1);

        let mut mie = Mie::new();
        let (Qext, Qsca) = mie.calculate(size_param, refractive_index);

        assert!((Qext - 2.86165188243).abs() < 1e-7);
        assert!((Qsca - 1.66424911991).abs() < 1e-7);
    }
}

#![allow(non_snake_case)]

use ndarray::{Array1, Array2, ArrayView1, ArrayView2, ArrayViewMut1, Zip};
use num::abs;
use num::complex::Complex64;

pub struct MieOutput {
    pub Qext: Array1<f64>,
    pub Qsca: Array1<f64>,
    pub S1: Array2<Complex64>,
    pub S2: Array2<Complex64>,
    pub size_param: Array1<f64>,
    pub refractive_index: Complex64,
    pub cos_angles: Array1<f64>,
}

impl MieOutput {
    pub fn new(
        size_param: ArrayView1<f64>,
        refractive_index: Complex64,
        cos_angles: ArrayView1<f64>,
    ) -> Self {
        MieOutput {
            Qext: Array1::zeros(size_param.len()),
            Qsca: Array1::zeros(size_param.len()),
            S1: Array2::zeros((size_param.len(), cos_angles.len())),
            S2: Array2::zeros((size_param.len(), cos_angles.len())),
            size_param: size_param.to_owned(),
            refractive_index,
            cos_angles: cos_angles.to_owned(),
        }
    }
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

fn Dn_downwards(z: Complex64, N: usize, mut Dn: ArrayViewMut1<Complex64>) {
    Dn[N - 1] = Dn_lentz(z, N);
    for n in (2..N + 1).rev() {
        Dn[n - 2] = (n as f64) / z - 1.0 / ((n as f64) / z + Dn[n - 1]);
    }
}

fn Dn_upwards(z: Complex64, N: usize, mut Dn: ArrayViewMut1<Complex64>) {
    let j: Complex64 = Complex64::new(0.0, 1.0);
    let exp_2 = (-2.0 * j * z).exp();

    Dn[0] = -1.0 / z + (1.0 - exp_2) / ((1.0 - exp_2) / z - j * (1.0 + exp_2));

    for n in 2..N + 1 {
        Dn[n - 1] = -(n as f64) / z + 1.0 / ((n as f64) / z - Dn[n - 2]);
    }
}

fn calc_Dn(
    refractive_index: Complex64,
    size_param: f64,
    N: usize,
    Dn_array: ArrayViewMut1<Complex64>,
) {
    let z = refractive_index * size_param;

    if refractive_index.re < 1.0 || refractive_index.re > 10.0 || abs(refractive_index.im) > 10.0 {
        Dn_downwards(z, N, Dn_array);
    } else {
        let temp = 3.9 - 10.8 * refractive_index.re + 13.78 * refractive_index.re.powf(2.0);

        if abs(refractive_index.im) * size_param >= temp {
            Dn_downwards(z, N, Dn_array);
        } else {
            Dn_upwards(z, N, Dn_array);
        }
    }
}

fn An_Bn(
    refractive_index: Complex64,
    size_param: f64,
    N: usize,
    mut An: ArrayViewMut1<Complex64>,
    mut Bn: ArrayViewMut1<Complex64>,
    Dn: ArrayView1<Complex64>,
) {
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

fn tau_pi_matrices(cos_angles: ArrayView1<f64>, max_order: usize) -> (Array2<f64>, Array2<f64>) {
    let mut tau = Array2::<f64>::zeros((max_order, cos_angles.len()));
    let mut pi = Array2::<f64>::zeros((max_order, cos_angles.len()));

    let N = max_order;

    pi.slice_mut(ndarray::s![0, ..]).fill(1.0);
    tau.slice_mut(ndarray::s![0, ..]).assign(&cos_angles);

    if N > 1 {
        Zip::from(pi.slice_mut(ndarray::s![1, ..]))
            .and(tau.slice_mut(ndarray::s![1, ..]))
            .and(cos_angles)
            .for_each(|pi, tau, cos_angle| {
                *pi = 3.0 * cos_angle;
                *tau = 2.0 * cos_angle * *pi - 3.0;
            });
    }

    for n in 3..N + 1 {
        let nf = n as f64;

        // Copy out rows that we need for reads
        let pi_n2 = pi.slice(ndarray::s![n - 2, ..]).to_owned();
        let pi_n3 = pi.slice(ndarray::s![n - 3, ..]).to_owned();

        // Compute new values into a temporary buffer
        Zip::from(pi.slice_mut(ndarray::s![n - 1, ..]))
            .and(&pi_n2)
            .and(&pi_n3)
            .and(cos_angles)
            .for_each(|r, &p2, &p3, &a| {
                *r = ((2.0 * nf - 1.0) * a * p2 - nf * p3) / (nf - 1.0);
            });

        // Repeat the pattern for tau
        let pi_n1_read = pi.slice(ndarray::s![n - 1, ..]).to_owned();
        let pi_n2_read = pi.slice(ndarray::s![n - 2, ..]).to_owned();
        Zip::from(tau.slice_mut(ndarray::s![n - 1, ..]))
            .and(&pi_n1_read)
            .and(&pi_n2_read)
            .and(cos_angles)
            .for_each(|r, &p1, &p2, &a| {
                *r = nf * a * p1 - (nf + 1.0) * p2;
            });
    }

    (tau, pi)
}

#[allow(clippy::too_many_arguments)]
fn regular_Q_S(
    size_param: f64,
    N: usize,
    An: ArrayView1<Complex64>,
    Bn: ArrayView1<Complex64>,
    tau: ArrayView2<f64>,
    pi: ArrayView2<f64>,
    mut S1: ArrayViewMut1<Complex64>,
    mut S2: ArrayViewMut1<Complex64>,
) -> (f64, f64) {
    let Qext: f64 = An
        .iter()
        .zip(Bn.iter())
        .enumerate()
        .map(|(n, (An, Bn))| {
            let n = n as f64;

            (2.0 * (n + 1.0) + 1.0) * (An.re + Bn.re)
        })
        .sum();

    let Qsca: f64 = An
        .iter()
        .zip(Bn.iter())
        .enumerate()
        .map(|(n, (An, Bn))| {
            let n = n as f64;

            (2.0 * (n + 1.0) + 1.0) * (An.norm_sqr() + Bn.norm_sqr())
        })
        .sum();

    S1.fill(Complex64::new(0.0, 0.0));
    S2.fill(Complex64::new(0.0, 0.0));

    for i in 0..N {
        let n = i as f64;
        let n_factor = (2.0 * (n + 1.0) + 1.0) / (n + 1.0) / (n + 2.0);

        let An = An[i];
        let Bn = Bn[i];

        Zip::from(&mut S1)
            .and(&mut S2)
            .and(&tau.row(i))
            .and(&pi.row(i))
            .for_each(|s1, s2, &tau, &pi| {
                *s1 += n_factor * (An * pi + Bn * tau);
                *s2 += n_factor * (An * tau + Bn * pi);
            });
    }

    (
        Qext * 2.0 / (size_param * size_param),
        Qsca * 2.0 / (size_param * size_param),
    )
}

fn small_Q_S(
    refractive_index: Complex64,
    size_param: f64,
    cos_angles: ArrayView1<f64>,
    mut S1: ArrayViewMut1<Complex64>,
    mut S2: ArrayViewMut1<Complex64>,
) -> (f64, f64) {
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
    Zip::from(&mut S1)
        .and(&mut S2)
        .and(&cos_angles)
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

pub fn mie(
    size_param: ArrayView1<f64>,
    refractive_index: Complex64,
    cos_angles: ArrayView1<f64>,
) -> MieOutput {
    let mut output = MieOutput::new(size_param, refractive_index, cos_angles);

    let max_x = output
        .size_param
        .clone()
        .into_iter()
        .reduce(f64::max)
        .unwrap_or(0.0);

    let N = max_order(max_x);

    // Local memory
    let mut Dn = Array1::<Complex64>::zeros(N);
    let mut An = Array1::<Complex64>::zeros(N);
    let mut Bn = Array1::<Complex64>::zeros(N);

    let (tau, pi) = tau_pi_matrices(cos_angles, N);

    Zip::indexed(size_param)
        .and(&mut output.Qext)
        .and(&mut output.Qsca)
        .for_each(|i, &x, Qext, Qsca| {
            let S1_slice = output.S1.slice_mut(ndarray::s![i, ..]);
            let S2_slice = output.S2.slice_mut(ndarray::s![i, ..]);
            if refractive_index.norm() * x < 0.1 {
                (*Qext, *Qsca) = small_Q_S(refractive_index, x, cos_angles, S1_slice, S2_slice);
            } else {
                let current_N = max_order(x);

                let An_slice = An.slice_mut(ndarray::s![..current_N]);
                let Bn_slice = Bn.slice_mut(ndarray::s![..current_N]);
                let Dn_slice = Dn.slice_mut(ndarray::s![..current_N]);

                calc_Dn(refractive_index, x, current_N, Dn_slice);
                let Dn_slice = Dn.slice(ndarray::s![..current_N]);

                An_Bn(refractive_index, x, current_N, An_slice, Bn_slice, Dn_slice);
                let An_slice = An.slice(ndarray::s![..current_N]);
                let Bn_slice = Bn.slice(ndarray::s![..current_N]);

                (*Qext, *Qsca) = regular_Q_S(
                    x,
                    current_N,
                    An_slice,
                    Bn_slice,
                    tau.view(),
                    pi.view(),
                    S1_slice,
                    S2_slice,
                );
            }
        });

    output
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_cos_angles() -> Array1<f64> {
        let mut cos_angles = Vec::new();
        for i in 0..7 {
            let angle = (i as f64) * std::f64::consts::PI / 6.0;
            cos_angles.push(angle.cos());
        }
        Array1::from_vec(cos_angles)
    }

    #[test]
    fn test_miev0_case_8() {
        let size_param = Array1::from_elem(1, 1000.0);
        let refractive_index = Complex64::new(0.75, 0.0);
        let cos_angles = default_cos_angles();

        let mie = mie(size_param.view(), refractive_index, cos_angles.view());

        // Qsca
        assert!((mie.Qsca[[0]] - 1.997908).abs() < 1e-6);

        // S1(0) and S2(0)
        assert!((mie.S1[[0, 0]].re - 4.994770e5).abs() < 1e-1);
        assert!((mie.S1[[0, 0]].im + 1.336502e4).abs() < 1e-2);
        assert!((mie.S2[[0, 0]].re - 4.994770e5).abs() < 1e-1);
        assert!((mie.S2[[0, 0]].im + 1.336502e4).abs() < 1e-2);

        // S1(1) and S2(1)
        assert!((mie.S1[[0, 1]].re + 3.999296e2).abs() < 1e-4);
        assert!((mie.S1[[0, 1]].im + 3.316361e2).abs() < 1e-4);
        assert!((mie.S2[[0, 1]].re + 3.946018e2).abs() < 1e-4);
        assert!((mie.S2[[0, 1]].im + 1.147791e2).abs() < 1e-4);

        // S1(2) and S2(2)
        assert!((mie.S1[[0, 2]].re + 5.209852e2).abs() < 1e-4);
        assert!((mie.S1[[0, 2]].im + 5.776614e2).abs() < 1e-4);
        assert!((mie.S2[[0, 2]].re + 1.970767e2).abs() < 1e-4);
        assert!((mie.S2[[0, 2]].im + 6.937470e2).abs() < 1e-4);

        // S1(3) and S2(3)
        assert!((mie.S1[[0, 3]].re + 1.600887e2).abs() < 1e-4);
        assert!((mie.S1[[0, 3]].im - 1.348013e2).abs() < 1e-4);
        assert!((mie.S2[[0, 3]].re + 4.152365e1).abs() < 1e-5);
        assert!((mie.S2[[0, 3]].im - 1.143000e2).abs() < 1e-4);

        // S1(4) and S2(4)
        assert!((mie.S1[[0, 4]].re - 8.431720e1).abs() < 1e-5);
        assert!((mie.S1[[0, 4]].im + 1.209493e2).abs() < 1e-4);
        assert!((mie.S2[[0, 4]].re + 4.261732e1).abs() < 1e-5);
        assert!((mie.S2[[0, 4]].im - 5.535055e1).abs() < 1e-5);

        // S1(5) and S2(5)
        assert!((mie.S1[[0, 5]].re + 7.556092e1).abs() < 1e-5);
        assert!((mie.S1[[0, 5]].im + 8.134810e1).abs() < 1e-5);
        assert!((mie.S2[[0, 5]].re - 4.218303e1).abs() < 1e-5);
        assert!((mie.S2[[0, 5]].im - 9.100831e1).abs() < 1e-5);

        // S1(6) and S2(6)
        assert!((mie.S1[[0, 6]].re - 1.705778e1).abs() < 1e-5);
        assert!((mie.S1[[0, 6]].im - 4.842510e2).abs() < 1e-4);
        assert!((mie.S2[[0, 6]].re + 1.705778e1).abs() < 1e-5);
        assert!((mie.S2[[0, 6]].im + 4.842510e2).abs() < 1e-4);
    }

    #[test]
    fn test_miev0_case_5() {
        let size_param = 0.099;
        let refractive_index = Complex64::new(0.75, 0.0);

        let size_param = Array1::from_elem(1, size_param);
        let cos_angles = default_cos_angles();

        let mie = mie(size_param.view(), refractive_index, cos_angles.view());

        assert!((mie.Qsca[[0]] - 0.000007).abs() < 1e-5);

        assert!((mie.S1[[0, 0]].re - 1.817558e-8).abs() < 1e-8);
        assert!((mie.S1[[0, 0]].im + 1.654225e-4).abs() < 1e-8);
        /*
        assert!((mie.S2[0].re - 1.817558e-8).abs() < 1e-8);
        assert!((mie.S2[0].im + 1.654225e-4).abs() < 1e-8);

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
        */
    }
}

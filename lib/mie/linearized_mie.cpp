#include <sasktran2/mie/mie.h>
#include <sasktran2/mie/linearized_mie.h>

namespace sasktran2::mie {

    LinearizedMie::LinearizedMie() {}

    void LinearizedMie::internal_calculate(
        const Eigen::VectorXd& size_param,
        const std::complex<double>& refractive_index,
        const Eigen::VectorXd& cos_angles, bool calculate_derivative,
        MieOutput& output) {
        // this is where my code is going to go: first start with calculating
        // the parameters Q ext, Q sca, S1, S2 size param is x refractive index
        // is m, complex number Angles needed S1 and S2

        // output is a pointer to the structure MieOutput, so we need to fill
        // those values
        // TODO: actually type in the equations here
        MieData data_values;
        Eigen::VectorXd Qext;
        Eigen::VectorXd Qsca;
        Eigen::MatrixXcd S1; // Scattering amplitude1 [nsize, nangles]
        Eigen::MatrixXcd S2; // Scattering amplitude2 [nsize, nangles]

        auto res = (abs(refractive_index) * size_param.array() < 0.1);
        // wherever res is true, should use small mie approx, wherever res if
        // False, should do regular
        auto iterator = std::find(res.begin(), res.end(), false);
        int index;
        if (iterator == res.end() || refractive_index.real() <= 0.0) {
            Small_Q_S(refractive_index, size_param, cos_angles, Qext, Qsca, S1,
                      S2);
            // small mie for all
        } else if (iterator == res.begin()) {
            Regular_Q_S(refractive_index, size_param, cos_angles, Qext, Qsca,
                        S1, S2);
            // regular method
        } else {
            // we need to split it
            index = iterator - res.begin();
            Eigen::VectorXd small_array = size_param(Eigen::seq(0, index - 1));
            Eigen::VectorXd reg_array =
                size_param(Eigen::seq(index, size_param.size() - 1));

            // regular
            Eigen::VectorXd Qext_reg;
            Eigen::VectorXd Qsca_reg;
            Eigen::MatrixXcd S1_reg;
            Eigen::MatrixXcd S2_reg;
            Regular_Q_S(refractive_index, reg_array, cos_angles, Qext_reg,
                        Qsca_reg, S1_reg, S2_reg);
            // small
            Eigen::VectorXd Qext_small;
            Eigen::VectorXd Qsca_small;
            Eigen::MatrixXcd S1_small;
            Eigen::MatrixXcd S2_small;
            Small_Q_S(refractive_index, small_array, cos_angles, Qext_small,
                      Qsca_small, S1_small, S2_small);
            // rejoin down here

            Qext.resize(size_param.size());
            ;
            Qext << Qext_small, Qext_reg;

            Qsca.resize(size_param.size());
            ;
            Qsca << Qsca_small, Qsca_reg;

            S1.resize(size_param.size(), cos_angles.size());
            S1 << S1_small, S1_reg;
            S2.resize(size_param.size(), cos_angles.size());
            S2 << S2_small, S2_reg;
        }

        data_values.Qext = Qext;
        data_values.Qsca = Qsca;
        data_values.S1 = S1;
        data_values.S2 = S2;

        // S1 S2 calculation

        // first need 3D tensor that represents the multiplication
        // the n terms: (2n+1)/(n(n+1))

        // need to transform an and bn from 2D matrix to 3D tensor - just need
        // to replicate it however many times there are angles

        // then need the tau and pi functions - these need to be calculated.
        // These are angle dependent, n dependent, but not x dependent

        // then sum along the N dimension to represent summing over n=1 to N

        output.values = data_values;
    }

    void
    LinearizedMie::Regular_Q_S(const std::complex<double>& refractive_index,
                               const Eigen::VectorXd& size_param,
                               const Eigen::VectorXd& cos_angles,
                               Eigen::VectorXd& Qext, Eigen::VectorXd& Qsca,
                               Eigen::MatrixXcd& S1, Eigen::MatrixXcd& S2) {

        // An Bn calculation - need for both
        double x = size_param.maxCoeff(); // largest size parameter, use to do
                                          // the calculations for highest N
        int N = int(x + 4.05 * pow(x, 0.33333) + 2.0) + 2;
        Eigen::MatrixXcd An_matrix;
        Eigen::MatrixXcd Bn_matrix;
        An_Bn(refractive_index, size_param, N, An_matrix, Bn_matrix);

        Eigen::VectorXd n2_1 =
            Eigen::VectorXd::LinSpaced(N, 1, N).array() * 2.0 + 1.0;
        Eigen::MatrixXd n2_1_matrix;
        n2_1_matrix.resize(N, size_param.size());

        for (int i = 0; i < size_param.size(); i++) {
            n2_1_matrix(Eigen::all, i) = n2_1.transpose();
        }

        // Qext Qsca calculation
        Eigen::MatrixXd qext_temp_matrix =
            n2_1_matrix.array() *
            (An_matrix.array().real() + Bn_matrix.array().real());
        Qext = qext_temp_matrix.colwise().sum();
        Qext.array() *= 2.0 / size_param.array().square();

        Eigen::MatrixXd qsca_temp_matrix =
            n2_1_matrix.array() *
            (abs(An_matrix.array()).square() + abs(Bn_matrix.array()).square());
        Qsca = qsca_temp_matrix.colwise().sum();
        Qsca.array() *= 2.0 / size_param.array().square();

        // S1 S2 calculation
        Eigen::Tensor<std::complex<double>, 3> S1_temp_tensor;
        S1_temp_tensor.resize(N, size_param.size(), cos_angles.size());
        Eigen::Tensor<std::complex<double>, 3> S2_temp_tensor;
        S2_temp_tensor.resize(N, size_param.size(), cos_angles.size());

        Eigen::Tensor<double, 3> n_tensor;
        n_tensor.resize(N, size_param.size(), cos_angles.size());
        Eigen::Tensor<std::complex<double>, 3> an_tensor;
        an_tensor.resize(N, size_param.size(), cos_angles.size());
        Eigen::Tensor<std::complex<double>, 3> bn_tensor;
        bn_tensor.resize(N, size_param.size(), cos_angles.size());

        Eigen::VectorXd n = Eigen::VectorXd::LinSpaced(N, 1, N);
        n.array() = (2.0 * n.array() + 1.0) / n.array() / (n.array() + 1.0);
        Eigen::MatrixXd temp = n.rowwise().replicate(size_param.size());

        Eigen::TensorMap<Eigen::Tensor<double, 2>> n_matrix(temp.data(), N,
                                                            size_param.size());
        Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 2>> An_matrix_map(
            An_matrix.data(), N, size_param.size());
        Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 2>> Bn_matrix_map(
            Bn_matrix.data(), N, size_param.size());
        for (int i = 0; i < cos_angles.size(); i++) {
            n_tensor.chip(i, 2) = n_matrix;
            an_tensor.chip(i, 2) = An_matrix_map;
            bn_tensor.chip(i, 2) = Bn_matrix_map;
        }

        Eigen::Tensor<double, 3> tau_tensor;
        tau_tensor.resize(N, size_param.size(), cos_angles.size());
        Eigen::Tensor<double, 3> pi_tensor;
        pi_tensor.resize(N, size_param.size(), cos_angles.size());
        Tau_Pi(size_param, cos_angles, N, tau_tensor, pi_tensor);

        S1_temp_tensor =
            (an_tensor * pi_tensor + bn_tensor * tau_tensor); //*n_tensor;
        S2_temp_tensor =
            (an_tensor * tau_tensor + bn_tensor * pi_tensor); //*n_tensor;
        S1_temp_tensor = S1_temp_tensor * n_tensor;
        S2_temp_tensor = S2_temp_tensor * n_tensor;

        // adding up along the N dimension
        Eigen::Tensor<std::complex<double>, 2> S1_matrix =
            S1_temp_tensor.sum(Eigen::array<int, 1>({0}));
        Eigen::Tensor<std::complex<double>, 2> S2_matrix =
            S2_temp_tensor.sum(Eigen::array<int, 1>({0}));
        S1 = Eigen::Map<Eigen::MatrixXcd>(S1_matrix.data(), size_param.size(),
                                          cos_angles.size());
        S2 = Eigen::Map<Eigen::MatrixXcd>(S2_matrix.data(), size_param.size(),
                                          cos_angles.size());
    }

    void LinearizedMie::Small_Q_S(const std::complex<double>& refractive_index,
                                  const Eigen::VectorXd& size_param,
                                  const Eigen::VectorXd& cos_angles,
                                  Eigen::VectorXd& Qext, Eigen::VectorXd& Qsca,
                                  Eigen::MatrixXcd& S1, Eigen::MatrixXcd& S2) {
        // Mie Scattering Calculations: Advances in Technique and Fast,
        // Vector-Speed Computer Codes, Section 4 for small spheres

        std::complex<double> m_2 = pow(refractive_index, 2);
        Eigen::VectorXd x_2 = size_param.array().square();
        std::complex<double> j(0.0, 1.0);

        Eigen::VectorXcd N_1 =
            (1.0 - 0.1 * x_2.array() +
             (4.0 * m_2 + 5.0) * x_2.array().square() / 1400.0);

        Eigen::VectorXcd D_1 =
            m_2 + 2.0 + (1.0 - 0.7 * m_2) * x_2.array() -
            (8.0 * pow(refractive_index, 4) - 385.0 * m_2 + 350.0) *
                x_2.array().square() / 1400.0;
        D_1 = D_1.array() + 2.0 * j * (m_2 - 1.0) * pow(size_param.array(), 3) *
                                (1.0 - 0.1 * x_2.array()) / 3.0;

        Eigen::VectorXcd a_hat1 =
            2.0 * j * (m_2 - 1.0) / 3.0 * N_1.array() / D_1.array();

        Eigen::VectorXcd b_hat1 =
            (j * x_2.array() * (m_2 - 1.0) / 45.0 *
             (1.0 + (2.0 * m_2 - 5.0) / 70.0 * x_2.array())) /
            (1.0 - x_2.array() * (2.0 * m_2 - 5.0) / 30.0);

        Eigen::VectorXcd a_hat2 =
            (j * x_2.array() * (m_2 - 1.0) / 15.0 *
             (1.0 - x_2.array() / 14.0)) /
            (2.0 * m_2 + 3.0 - x_2.array() * (2.0 * m_2 - 7.0) / 14.0);

        Eigen::VectorXd T = pow(abs(a_hat1.array()), 2) +
                            pow(abs(b_hat1.array()), 2) +
                            5.0 / 3.0 * pow(abs(a_hat2.array()), 2);
        Qsca.resize(size_param.size());
        Qsca = 6.0 * x_2.array().square() * T.array();

        if (refractive_index.imag() == 0) {
            Qext = Qsca;
        } else {
            Qext = 6.0 * size_param.array() *
                   ((a_hat1.array() + b_hat1.array()).real() +
                    5.0 / 3.0 * a_hat2.array().real());
        }

        // need to deal with matrices here
        // S1 and S2 are of size size_param.size x cos_angles size
        S1.resize(size_param.size(), cos_angles.size());
        S2.resize(size_param.size(), cos_angles.size());

        Eigen::MatrixXd size_param_matrix =
            size_param.rowwise().replicate(cos_angles.size());
        Eigen::MatrixXcd a_hat1_matrix =
            a_hat1.rowwise().replicate(cos_angles.size());
        Eigen::MatrixXcd a_hat2_matrix =
            a_hat2.rowwise().replicate(cos_angles.size());
        Eigen::MatrixXcd b_hat1_matrix =
            b_hat1.rowwise().replicate(cos_angles.size());
        Eigen::MatrixXd mu_matrix =
            cos_angles.rowwise().replicate(size_param.size());
        mu_matrix.array() = mu_matrix.transpose();

        S1 = 3.0 / 2.0 * pow(size_param_matrix.array(), 3) *
             (a_hat1_matrix.array() +
              mu_matrix.array() *
                  (b_hat1_matrix.array() + 5.0 / 3.0 * a_hat2_matrix.array()));
        S2 =
            3.0 / 2.0 * pow(size_param_matrix.array(), 3) *
            (b_hat1_matrix.array() + a_hat1_matrix.array() * mu_matrix.array() +
             5.0 / 3.0 * a_hat2_matrix.array() *
                 (2.0 * pow(mu_matrix.array(), 2) - 1.0));
    }

    void LinearizedMie::An_Bn(const std::complex<double>& refractive_index,
                              const Eigen::VectorXd& size_param, const int N,
                              Eigen::MatrixXcd& An_matrix,
                              Eigen::MatrixXcd& Bn_matrix) {
        // here we are calculating An and bn coefficients
        // start by finding Dns.
        An_matrix.resize(N, size_param.size());
        Bn_matrix.resize(N, size_param.size());

        std::complex<double> j(0.0, 1.0);
        Eigen::VectorXcd
            temp; // temporary variable, will need when propagating psis and xis
        Eigen::VectorXcd psi_n_1 = sin(size_param.array()); // this is n=0
        Eigen::VectorXcd psi_n = sin(size_param.array()) / size_param.array() -
                                 cos(size_param.array()); // this is n=1
        Eigen::VectorXcd xi_n_1 = sin(size_param.array()) +
                                  j * cos(size_param.array()); // this is n=0
        Eigen::VectorXcd xi_n =
            sin(size_param.array()) / size_param.array() -
            cos(size_param.array()) +
            j * (cos(size_param.array()) / size_param.array() +
                 sin(size_param.array())); // this is n=1

        Eigen::MatrixXcd Dn_matrix;
        Dn(Dn_matrix, refractive_index, size_param, N);

        // initial values for psi and xi functions

        // now calculate An and Bn. need from n =1 to n=N
        Eigen::VectorXcd An =
            ((Dn_matrix.row(0).transpose().array() / refractive_index +
              1.0 / size_param.array()) *
                 psi_n.array() -
             psi_n_1.array()) /
            ((Dn_matrix.row(0).transpose().array() / refractive_index +
              1.0 / size_param.array()) *
                 xi_n.array() -
             xi_n_1.array());
        An_matrix(0, Eigen::all) = An; // this is for n=1, stored in row 0
        Eigen::VectorXcd Bn =
            ((Dn_matrix.row(0).transpose().array() * refractive_index +
              1.0 / size_param.array()) *
                 psi_n.array() -
             psi_n_1.array()) /
            ((Dn_matrix.row(0).transpose().array() * refractive_index +
              1.0 / size_param.array()) *
                 xi_n.array() -
             xi_n_1.array());
        Bn_matrix(0, Eigen::all) = Bn; // this is for n=1, stored in row 0

        temp = xi_n_1; // this is n=0
        xi_n_1 = xi_n; // this is n=1
        xi_n = (2.0 + 1.0) / size_param.array() * xi_n_1.array() -
               temp.array(); // this is n=2

        psi_n_1 = psi_n;     // this is n=1
        psi_n = xi_n.real(); // this is n=2

        for (int n = 2; n < N + 1; n++) {

            An = ((Dn_matrix.row(n - 1).transpose().array() / refractive_index +
                   n / size_param.array()) *
                      psi_n.array() -
                  psi_n_1.array()) /
                 ((Dn_matrix.row(n - 1).transpose().array() / refractive_index +
                   n / size_param.array()) *
                      xi_n.array() -
                  xi_n_1.array());
            An_matrix(n - 1, Eigen::all) =
                An; // in matrix, the value for n is stored in row n-1
            Bn = ((Dn_matrix.row(n - 1).transpose().array() * refractive_index +
                   n / size_param.array()) *
                      psi_n.array() -
                  psi_n_1.array()) /
                 ((Dn_matrix.row(n - 1).transpose().array() * refractive_index +
                   n / size_param.array()) *
                      xi_n.array() -
                  xi_n_1.array());
            Bn_matrix(n - 1, Eigen::all) =
                Bn; // in matrix, the value for n is stored in row n-1

            // update psis and xis
            temp = xi_n_1;
            xi_n_1 = xi_n;
            xi_n = (2 * n + 1) / size_param.array() * xi_n_1.array() -
                   temp.array();

            psi_n_1 = psi_n;
            psi_n = xi_n.real(); // allowed since xi_n = psi_n + i*X_n
        }
    }

    void LinearizedMie::Tau_Pi(const Eigen::VectorXd& size_param,
                               const Eigen::VectorXd& cos_angles, const int N,
                               Eigen::Tensor<double, 3>& tau_tensor,
                               Eigen::Tensor<double, 3>& pi_tensor) {

        Eigen::MatrixXd Pi_matrix;
        Eigen::MatrixXd Tau_matrix;
        Pi_matrix.resize(N, cos_angles.size());
        Tau_matrix.resize(N, cos_angles.size());

        // have to start at 0, and then make our way up all the way to N
        Eigen::VectorXd temp;
        Eigen::VectorXd pi_n_1 =
            Eigen::VectorXd::Constant(cos_angles.size(), 0); // n=0
        Eigen::VectorXd pi_n =
            Eigen::VectorXd::Constant(cos_angles.size(), 1); // n=1
        Eigen::VectorXd tau_n = 1 * cos_angles.array() * pi_n.array() -
                                (1 + 1) * pi_n_1.array(); // n=1

        Pi_matrix(0, Eigen::all) = pi_n;   // this is for n=1, stored in row 0
        Tau_matrix(0, Eigen::all) = tau_n; // this is for n=1, stored in row 0

        for (int n = 2; n < N + 1; n++) {
            temp = pi_n;
            pi_n = ((2 * n - 1) * cos_angles.array() * pi_n.array() -
                    n * pi_n_1.array()) /
                   (n - 1);
            pi_n_1 = temp;
            tau_n = n * cos_angles.array() * pi_n.array() -
                    (n + 1) * pi_n_1.array();
            Pi_matrix(n - 1, Eigen::all) =
                pi_n; // this is for n, stored in row n-1
            Tau_matrix(n - 1, Eigen::all) =
                tau_n; // this is for n, stored in row n-1
        }

        Eigen::TensorMap<Eigen::Tensor<double, 2>> tau_map(Tau_matrix.data(), N,
                                                           cos_angles.size());
        Eigen::TensorMap<Eigen::Tensor<double, 2>> pi_map(Pi_matrix.data(), N,
                                                          cos_angles.size());

        for (int i = 0; i < size_param.size(); i++) {
            tau_tensor.chip(i, 1) = tau_map;
            pi_tensor.chip(i, 1) = pi_map;
        }
    }

    void LinearizedMie::Dn(Eigen::MatrixXcd& Dn_matrix,
                           const std::complex<double>& refractive_index,
                           const Eigen::VectorXd& size_param, const int N) {
        // check which way to go - upwards or downwards, then go for it,
        // concatenating arrays at the end if real part of m is less than 1 or
        // bigger than 10, abs of imaginary is bigger than 10, or x*
        // m_imag>= 3.9 - 10.8 * mreal + 13.78 * mreal**2
        double m_real = refractive_index.real();
        double m_imag = refractive_index.imag();

        Eigen::VectorXcd z_array = size_param.array() * refractive_index;

        if (m_real < 1 || m_real > 10 || abs(m_imag) > 10) {
            // downwards
            Dn_downwards(z_array, N, Dn_matrix);

        } else {
            // need to split and figure out what's up and what's down
            // assuming the size parameters are ordered smallest to largest
            double temp = 3.9 - 10.8 * m_real + 13.78 * pow(m_real, 2);
            auto res = size_param.array() * abs(m_imag) >=
                       3.9 - 10.8 * m_real + 13.78 * pow(m_real, 2);
            // wherever res is true, should do downwards, wherever res if False,
            // should do upwards
            auto iterator = std::find(res.begin(), res.end(), true);
            int index;
            if (iterator == res.end()) {
                Dn_upwards(z_array, N, Dn_matrix);
            } else if (iterator == res.begin()) {
                Dn_downwards(z_array, N, Dn_matrix);
            } else {
                // we need to split it
                index = iterator - res.begin();
                Eigen::VectorXcd up_array = z_array(Eigen::seq(0, index - 1));
                Eigen::VectorXcd down_array =
                    z_array(Eigen::seq(index, z_array.size() - 1));

                // upwards
                Eigen::MatrixXcd Dn_matrix_up;
                Dn_upwards(up_array, N, Dn_matrix_up);
                // downwards
                Eigen::MatrixXcd Dn_matrix_down;
                Dn_downwards(down_array, N, Dn_matrix_down);
                // TODO rejoin down here

                Dn_matrix.resize(N, z_array.size());
                ;
                Dn_matrix << Dn_matrix_up, Dn_matrix_down;
            }
        }
    }

    void LinearizedMie::Dn_upwards(const Eigen::VectorXcd& z, const int N,
                                   Eigen::MatrixXcd& Dn_matrix) {
        // here we are using upward recurrence to calculate Dn

        // first calculate D1
        std::complex<double> j(0.0, 1.0);
        Eigen::VectorXcd temp = -2.0 * j * z.array();
        Eigen::VectorXcd exp_2 = temp.array().exp();
        Eigen::VectorXcd prev_Dn =
            -1.0 / z.array() +
            (1.0 - exp_2.array()) /
                ((1.0 - exp_2.array()) / z.array() - j * (1.0 + exp_2.array()));

        // now calculate each successive term and store it in the Dn array
        Dn_matrix.resize(N, z.size());
        Dn_matrix(0, Eigen::all) = prev_Dn;

        for (int n = 2; n < N + 1; n++) {
            prev_Dn = -n / z.array() + 1.0 / (n / z.array() - prev_Dn.array());
            // put Dns in the matrix
            Dn_matrix(n - 1, Eigen::all) = prev_Dn;
        }
    }

    void LinearizedMie::Dn_downwards(const Eigen::VectorXcd& z, const int N,
                                     Eigen::MatrixXcd& Dn_matrix) {
        // here we are using downward recurrence to calculate Dn

        // first calculate the last Dn we will need
        Eigen::VectorXcd prev_Dn;
        Dn_Lentz(z, N, prev_Dn);
        Dn_matrix.resize(N, z.size());
        Dn_matrix(N - 1, Eigen::all) = prev_Dn;
        // now calculate each term before it, and store it in the big Dn thing
        for (int n = N; n > 1; n--) {
            prev_Dn = n / z.array() - 1 / (n / z.array() + prev_Dn.array());
            // TODO put Dns in the matrix
            Dn_matrix(n - 2, Eigen::all) = prev_Dn;
        }
    }

    void LinearizedMie::Dn_Lentz(const Eigen::VectorXcd& z, const int N,
                                 Eigen::VectorXcd& Dn_array) {
        // keeping track
        Eigen::VectorXcd z_inv = 2.0 / z.array();

        Eigen::VectorXcd alpha_1 = (N + 0.5) * z_inv.array(); // [a1], sam as a1
        Eigen::VectorXcd a_j =
            -(N + 1.5) *
            z_inv
                .array(); // aj term to be used for alpha j1 and j2 calculations
        Eigen::VectorXcd alpha_j1 =
            a_j.array() + 1.0 / alpha_1.array(); // [aj,... a1]
        Eigen::VectorXcd alpha_j2 = a_j.array(); // [aj,... a2]

        Eigen::VectorXcd cur_ratio = alpha_j1.array() / alpha_j2.array();
        Eigen::VectorXcd overall_ratio = alpha_1.array() * cur_ratio.array();

        while ((cur_ratio.array().abs() - 1).abs().maxCoeff() > 1e-12) {
            // iterating until we are within some number of decimal places

            a_j = z_inv - a_j;
            alpha_j1 = 1.0 / alpha_j1.array() + a_j.array();
            alpha_j2 = 1.0 / alpha_j2.array() + a_j.array();

            z_inv *= -1;
            cur_ratio = alpha_j1.array() / alpha_j2.array();
            overall_ratio = overall_ratio.array() * cur_ratio.array();
        }
        Dn_array = -N / z.array() + overall_ratio.array();
    }

} // namespace sasktran2::mie

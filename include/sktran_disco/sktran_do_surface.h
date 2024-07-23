#pragma once
#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_quadrature.h"
#include "sktran_disco/sktran_do_types.h"
#include <sasktran2/atmosphere/surface.h>

namespace sasktran_disco {
    template <int NSTOKES, int CNSTR = -1> struct SurfaceStorage {
        static constexpr int CSIZE = CNSTR == -1 ? Eigen::Dynamic : CNSTR / 2;

        struct SurfaceExpansion {
            Eigen::Matrix<double, -1, CSIZE> los_stream;
            Eigen::Matrix<double, CSIZE, CSIZE> stream_stream;
            Eigen::Matrix<double, -1, 1> los_solar;
            Eigen::Matrix<double, CSIZE, 1> stream_solar;
        };

        SurfaceExpansion brdf;
        std::vector<SurfaceExpansion> d_brdf;

        Eigen::VectorXd m_quadrature_phi;
        Eigen::VectorXd m_quadrature_weight;

        // Geometry factors
        double csz;                          /** Cosine solar zenith */
        const std::vector<double>* mu;       /** Stream angles */
        const std::vector<LineOfSight>* los; /** Los for the LOS angles*/

        void resize(int nstr, int nlos, int nderiv) {
            brdf.los_stream.resize(nlos, nstr / 2);
            brdf.stream_stream.resize(nstr / 2, nstr / 2);
            brdf.los_solar.resize(nlos, 1);
            brdf.stream_solar.resize(nstr / 2, 1);

            d_brdf.resize(nderiv);
            for (auto& d_b : d_brdf) {
                d_b.los_stream.resize(nlos, nstr / 2);
                d_b.stream_stream.resize(nstr / 2, nstr / 2);
                d_b.los_solar.resize(nlos, 1);
                d_b.stream_solar.resize(nstr / 2, 1);
            }

            m_quadrature_phi = Eigen::Map<const Eigen::VectorXd>(
                getQuadratureAbscissae(512), 512);
            m_quadrature_weight = Eigen::Map<const Eigen::VectorXd>(
                getQuadratureWeights(512), 512);
        }

        double compute_expansion(
            AEOrder m,
            const sasktran2::atmosphere::Surface<NSTOKES>& sk2_surface,
            int wavel_idx, double mu_out, double mu_in) const {
            if (sk2_surface.max_azimuthal_order() == 1) {
                // Special case
                return sk2_surface.brdf(wavel_idx, mu_in, mu_out, 0)(0, 0) *
                       EIGEN_PI;
            }
            if (m >= sk2_surface.max_azimuthal_order()) {
                return 0;
            }
            // else
            double result = 0;

            for (int i = 0; i < m_quadrature_phi.size() / 2; i++) {
                double a1 = 0.5 * m_quadrature_phi(i) + 0.5;
                double a2 = -0.5 * m_quadrature_phi(i) + 0.5;
                double a3 = 0.5 * m_quadrature_phi(i) - 0.5;
                double a4 = -0.5 * m_quadrature_phi(i) - 0.5;
                double w = 0.5 * m_quadrature_weight(i);

                result += w *
                          sk2_surface.brdf(wavel_idx, mu_in, mu_out,
                                           EIGEN_PI * a1)(0, 0) *
                          cos(m * EIGEN_PI * a1);
                result += w *
                          sk2_surface.brdf(wavel_idx, mu_in, mu_out,
                                           EIGEN_PI * a2)(0, 0) *
                          cos(m * EIGEN_PI * a2);
                result += w *
                          sk2_surface.brdf(wavel_idx, mu_in, mu_out,
                                           EIGEN_PI * a3)(0, 0) *
                          cos(m * EIGEN_PI * a3);
                result += w *
                          sk2_surface.brdf(wavel_idx, mu_in, mu_out,
                                           EIGEN_PI * a4)(0, 0) *
                          cos(m * EIGEN_PI * a4);
            }
            return result * 0.5 * EIGEN_PI * (2.0 - kronDelta(m, 0));
        }

        double d_compute_expansion(
            AEOrder m,
            const sasktran2::atmosphere::Surface<NSTOKES>& sk2_surface,
            int wavel_idx, double mu_out, double mu_in, int deriv_index) const {
            if (sk2_surface.max_azimuthal_order() == 1) {
                // Special case
                return sk2_surface.d_brdf(wavel_idx, mu_in, mu_out, 0,
                                          deriv_index)(0, 0) *
                       EIGEN_PI;
            }
            if (m >= sk2_surface.max_azimuthal_order()) {
                return 0;
            }
            // else
            double result = 0;

            for (int i = 0; i < m_quadrature_phi.size() / 2; i++) {
                double a1 = 0.5 * m_quadrature_phi(i) + 0.5;
                double a2 = -0.5 * m_quadrature_phi(i) + 0.5;
                double a3 = 0.5 * m_quadrature_phi(i) - 0.5;
                double a4 = -0.5 * m_quadrature_phi(i) - 0.5;
                double w = 0.5 * m_quadrature_weight(i);

                result += w *
                          sk2_surface.d_brdf(wavel_idx, mu_in, mu_out,
                                             EIGEN_PI * a1, deriv_index)(0, 0) *
                          cos(m * EIGEN_PI * a1);
                result += w *
                          sk2_surface.d_brdf(wavel_idx, mu_in, mu_out,
                                             EIGEN_PI * a2, deriv_index)(0, 0) *
                          cos(m * EIGEN_PI * a2);
                result += w *
                          sk2_surface.d_brdf(wavel_idx, mu_in, mu_out,
                                             EIGEN_PI * a3, deriv_index)(0, 0) *
                          cos(m * EIGEN_PI * a3);
                result += w *
                          sk2_surface.d_brdf(wavel_idx, mu_in, mu_out,
                                             EIGEN_PI * a4, deriv_index)(0, 0) *
                          cos(m * EIGEN_PI * a4);
            }
            return result * 0.5 * EIGEN_PI * (2.0 - kronDelta(m, 0));
        }
    };

    template <int NSTOKES, int CNSTR = -1> class Surface {
      private:
        SurfaceStorage<NSTOKES, CNSTR>& m_storage;

        const sasktran2::atmosphere::Surface<NSTOKES>& m_sk2_surface;
        int m_wavel_idx;

      public:
        Surface(SurfaceStorage<NSTOKES, CNSTR>& storage,
                const sasktran2::atmosphere::Surface<NSTOKES>& sk2_surface,
                int wavel_idx)
            : m_storage(storage), m_sk2_surface(sk2_surface),
              m_wavel_idx(wavel_idx) {}

        const sasktran2::atmosphere::Surface<NSTOKES>& sk2_surface() const {
            return m_sk2_surface;
        }

        void calculate(AEOrder m) {
            int numderiv = m_storage.d_brdf.size();
            // For each stream
            for (int i = 0; i < m_storage.mu->size() / 2; i++) {
                // LOS stream
                for (int j = 0; j < m_storage.los->size(); j++) {
                    m_storage.brdf.los_stream(
                        (*m_storage.los)[j].unsorted_index, i) =
                        m_storage.compute_expansion(
                            m, m_sk2_surface, m_wavel_idx,
                            (*m_storage.los)[j].coszenith, (*m_storage.mu)[i]);

                    for (int k = 0; k < numderiv; ++k) {
                        m_storage.d_brdf[k].los_stream(
                            (*m_storage.los)[j].unsorted_index, i) =
                            m_storage.d_compute_expansion(
                                m, m_sk2_surface, m_wavel_idx,
                                (*m_storage.los)[j].coszenith,
                                (*m_storage.mu)[i], k);
                    }
                }

                // Stream Stream
                for (int j = 0; j < m_storage.mu->size() / 2; j++) {
                    m_storage.brdf.stream_stream(j, i) =
                        m_storage.compute_expansion(
                            m, m_sk2_surface, m_wavel_idx, (*m_storage.mu)[j],
                            (*m_storage.mu)[i]);

                    for (int k = 0; k < numderiv; ++k) {
                        m_storage.d_brdf[k].stream_stream(j, i) =
                            m_storage.d_compute_expansion(
                                m, m_sk2_surface, m_wavel_idx,
                                (*m_storage.mu)[j], (*m_storage.mu)[i], k);
                    }
                }

                // Stream solar
                m_storage.brdf.stream_solar(i) = m_storage.compute_expansion(
                    m, m_sk2_surface, m_wavel_idx, (*m_storage.mu)[i],
                    m_storage.csz);

                for (int k = 0; k < numderiv; ++k) {
                    m_storage.d_brdf[k].stream_solar(i) =
                        m_storage.d_compute_expansion(
                            m, m_sk2_surface, m_wavel_idx, (*m_storage.mu)[i],
                            m_storage.csz, k);
                }
            }

            // Los solar
            for (int i = 0; i < m_storage.los->size(); i++) {
                m_storage.brdf.los_solar((*m_storage.los)[i].unsorted_index) =
                    m_storage.compute_expansion(m, m_sk2_surface, m_wavel_idx,
                                                (*m_storage.los)[i].coszenith,
                                                m_storage.csz);

                for (int k = 0; k < numderiv; ++k) {
                    m_storage.d_brdf[k].los_solar(
                        (*m_storage.los)[i].unsorted_index) =
                        m_storage.d_compute_expansion(
                            m, m_sk2_surface, m_wavel_idx,
                            (*m_storage.los)[i].coszenith, m_storage.csz, k);
                }
            }
        }
        const SurfaceStorage<NSTOKES, CNSTR>& storage() const {
            return m_storage;
        }
        SurfaceStorage<NSTOKES, CNSTR>& storage() { return m_storage; }
    };

} // namespace sasktran_disco
